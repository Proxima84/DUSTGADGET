#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */
/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
    int i, j;
#ifdef SAVETEMP_LTR
    double mass, *tmass;
    mass = 0;
    tmass = malloc(NTask * sizeof(double));
#endif
    All.Time = All.TimeBegin;
    read_ic(All.InitCondFile);
    All.Time = All.TimeBegin;
    All.Ti_Current = 0;
    All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
    set_softenings();
    All.DustSize /= All.UnitLength_in_cm;
    All.DustRho /= All.UnitDensity_in_cgs;
    All.NumCurrentTiStep = 0; /* setup some counters */
    All.SnapshotFileCount = 0;
    All.NumForcesSinceLastDomainDecomp = 0;
    for (i = 0; i < NumPart; i++) /*  start-up initialization */
    {
        for (j = 0; j < 3; j++)
            P[i].GravAccel[j] = 0;
        P[i].Accretion = 0;
        P[i].Ti_endstep = 0;
        P[i].Ti_begstep = 0;
        P[i].OldAcc = 0;
        P[i].Radius = 0;
#ifdef NOSTARSELFGRAV
        P[i].ItsStar = 0;
#endif
        if (P[i].Type == 1)
        {
            P[i].Radius = All.DustSize;
        }
        else if (P[i].Type == 2)
        {
            P[i].Radius = 21.5 * pow(All.G * P[i].Mass / 4 / M_PI / M_PI * All.MinSizeTimestep
                                         * All.MinSizeTimestep,
                                     1. / 3.);
            printf("Accretion Radius:%f,Mass:%f\n", P[i].Radius, P[i].Mass);
#ifdef NOSTARSELFGRAV
            P[i].ItsStar = (P[i].Mass >= (All.StarMassIndicator * SOLAR_MASS / All.UnitMass_in_g));
#endif
#ifdef SAVETEMP_LTR
            mass += P[i].Mass;
#endif
        }
    }

#ifdef SAVETEMP_LTR
    MPI_Allgather(&mass, 1, MPI_DOUBLE, tmass, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    All.CentralMass = 0;
    for (i = 0; i < NTask; i++)
        All.CentralMass += tmass[i];
    free(tmass);
#endif
    for (i = 0; i < N_gas; i++) /* initialize sph_properties */
    {

        for (j = 0; j < 3; j++)
        {
            SphP[i].VelPred[j] = P[i].Vel[j];
            SphP[i].HydroAccel[j] = 0;
        }

        SphP[i].DtEntropy = 0;

        if (RestartFlag == 0)
        {
            SphP[i].Hsml = All.HsmlConstant;
            SphP[i].Density = -1;
        }
    }
    ngb_treeallocate(MAX_NGB);

    force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

    All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

    Flag_FullStep = 1; /* to ensure that Peano-Hilber order is done */
    DomainDecompositionFlag = 0;

    domain_Decomposition(); /* do initial domain decomposition (gives equal numbers of particles) */

    ngb_treebuild(); /* will build tree */


    setup_smoothinglengths();

    TreeReconstructFlag = 1;

/* at this point, the entropy variable normally contains the
 * internal energy, read in from the initial conditions file, unless the file
 * explicitly signals that the initial conditions contain the entropy directly.
 * Once the density has been computed, we can convert thermal energy to entropy.
 */

#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
    if (header.flag_entropy_instead_u == 0)
        for (i = 0; i < N_gas; i++)
            SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density, GAMMA_MINUS1);
#endif
}

/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
    int i, no, p;

    if (RestartFlag == 0)
    {
        for (i = 0; i < N_gas; i++)
        {
            no = Father[i];

            while (10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
            {
                p = Nodes[no].u.d.father;

                if (p < 0)
                    break;

                no = p;
            }
#ifdef HSMLCONSTANT
            SphP[i].Hsml = All.HsmlConstant;
#else
            SphP[i].Hsml
                = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3)
                * Nodes[no].len;
#endif
        }
    }

    density();
#ifdef HSMLCONSTANT
    double Num, *Ntemp;
    Ntemp = malloc(NTask * sizeof(double));
    Num = 0;
    for (i = 0; i < N_gas; i++)
    {
        Num += SphP[i].NumNgb;
    }
    Num /= N_gas;
    MPI_Allgather(&Num, 1, MPI_DOUBLE, Ntemp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    if (ThisTask == 0)
    {
        for (i = 1; i < NTask; i++)
            Num += Ntemp[i];
        Num /= NTask;
        printf("\nAverage number of neighboring particles: %f\n", Num);
    }
    free(Ntemp);
    dustGas();
    Ntemp = malloc(NTask * sizeof(double));
    Num = 0;
    for (i = 0; i < N_gas; i++)
    {
        Num += SphP[i].DustGasNgb;
    }
    Num /= N_gas;
    MPI_Allgather(&Num, 1, MPI_DOUBLE, Ntemp, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    if (ThisTask == 0)
    {
        for (i = 1; i < NTask; i++)
            Num += Ntemp[i];
        Num /= NTask;
        printf("\nAverage number of dust gas interaction particles: %f\n", Num);
    }
// for (i = 0; i < N_gas; i++)
//  printf("%i %f %f\n",P[i].ID,SphP[i].NumNgb,SphP[i].DustGasNgb);
#endif
}
#ifdef HSMLCONSTANT
void dustGas()
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int i, j, k, ndone, maxfill, source;
    int level, ngrp, recvTask, place, nexport;
    MPI_Status status;

    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);
    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&N_gas, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
    for (i = 0, ntot = 0; i < NTask; i++)
        ntot += numlist[i];
    free(numlist);
    /* we will repeat the whole thing for those particles where we didn't
     * find enough neighbours
     */
    i = 0; /* beginn with this index */
    ntotleft = ntot; /* particles left for all tasks together */
    while (ntotleft > 0)
    {
        for (j = 0; j < NTask; j++)
            nsend_local[j] = 0;

        /* do local particles and prepare export list */
        for (nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHydro - NTask; i++)
        {
            ndone++;
            for (j = 0; j < NTask; j++)
                Exportflag[j] = 0;
            dustGas_evaluate(i, 0);
            for (j = 0; j < NTask; j++)
            {
                if (Exportflag[j])
                {
                    for (k = 0; k < 3; k++)
                        DragDataIn[nexport].Pos[k] = P[i].Pos[k];
                    DragDataIn[nexport].Type = P[i].Type;
                    DragDataIn[nexport].Index = i;
                    DragDataIn[nexport].Task = j;
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
        qsort(DragDataIn, nexport, sizeof(struct dragdata_in), hydro_compare_key);

        for (j = 1, noffset[0] = 0; j < NTask; j++)
            noffset[j] = noffset[j - 1] + nsend_local[j - 1];

        MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

        /* now do the particles that need to be exported */

        for (level = 1; level < (1 << PTask); level++)
        {

            for (j = 0; j < NTask; j++)
                nbuffer[j] = 0;
            for (ngrp = level; ngrp < (1 << PTask); ngrp++)
            {
                maxfill = 0;
                for (j = 0; j < NTask; j++)
                {
                    if ((j ^ ngrp) < NTask)
                        if (maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
                            maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
                }
                if (maxfill >= All.BunchSizeDensity)
                    break;
                recvTask = ThisTask ^ ngrp;

                if (recvTask < NTask)
                {
                    if (nsend[ThisTask * NTask + recvTask] > 0
                        || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(&DragDataIn[noffset[recvTask]],
                            nsend_local[recvTask] * sizeof(struct dragdata_in), MPI_BYTE, recvTask,
                            TAG_DRAG_A, &DragDataGet[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct dragdata_in),
                            MPI_BYTE, recvTask, TAG_DRAG_A, MPI_COMM_WORLD, &status);
                    }
                }

                for (j = 0; j < NTask; j++)
                    if ((j ^ ngrp) < NTask)
                        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
            }


            for (j = 0; j < nbuffer[ThisTask]; j++)
                dustGas_evaluate(j, 1);

            MPI_Barrier(MPI_COMM_WORLD);


            /* get the result */

            for (j = 0; j < NTask; j++)
                nbuffer[j] = 0;
            for (ngrp = level; ngrp < (1 << PTask); ngrp++)
            {
                maxfill = 0;
                for (j = 0; j < NTask; j++)
                {
                    if ((j ^ ngrp) < NTask)
                        if (maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
                            maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
                }
                if (maxfill >= All.BunchSizeDensity)
                    break;

                recvTask = ThisTask ^ ngrp;

                if (recvTask < NTask)
                {
                    if (nsend[ThisTask * NTask + recvTask] > 0
                        || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(&DragDataResult[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct dragdata_out),
                            MPI_BYTE, recvTask, TAG_DRAG_B,
                            &DragDataPartialResult[noffset[recvTask]],
                            nsend_local[recvTask] * sizeof(struct dragdata_out), MPI_BYTE, recvTask,
                            TAG_DRAG_B, MPI_COMM_WORLD, &status);

                        /* add the result to the particles */
                        for (j = 0; j < nsend_local[recvTask]; j++)
                        {
                            source = j + noffset[recvTask];
                            place = DragDataIn[source].Index;
                            SphP[place].DustGasNgb += DragDataPartialResult[source].DustGasNgb;
                        }
                    }
                }
                for (j = 0; j < NTask; j++)
                    if ((j ^ ngrp) < NTask)
                        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
            }

            level = ngrp - 1;
        }

        MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);

        for (j = 0; j < NTask; j++)
            ntotleft -= ndonelist[j];
    }
    free(ndonelist);
    free(nsend);
    free(nsend_local);
    free(nbuffer);
    free(noffset);
}

void dustGas_evaluate(int target, int mode)
{
    int j, startnode, numngb, numngb_inbox;
    FLOAT* pos;
    FLOAT X[3];
    int ix[3], type;
    if (mode == 0)
    {
        pos = P[target].Pos;
        type = P[target].Type;
    }
    else
    {
        pos = DragDataGet[target].Pos;
        type = P[target].Type;
    }
#ifdef DIM1
    ix[0] = floor((All.BarrierDistance + pos[0]) / All.DustGasMechStep);
    X[0] = -All.BarrierDistance + (ix[0] + 0.5) * All.DustGasMechStep;
    X[1] = X[2] = 0;
#else
    for (j = 0; j < 3; j++)
    {
        ix[j] = floor((All.BarrierDistance + pos[j]) / All.DustGasMechStep);
        X[j] = -All.BarrierDistance + (ix[j] + 0.5) * All.DustGasMechStep;
    }
#endif
    /* initialize variables before SPH loop is started */
    startnode = All.MaxPart;
    numngb = 0;
    do
    {
        numngb_inbox = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, !type);
        numngb += numngb_inbox;
    } while (startnode >= 0);

    if (mode == 0)
        SphP[target].DustGasNgb = numngb;
    else
        DragDataResult[target].DustGasNgb = numngb;
}
#endif
