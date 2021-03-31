#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

/*! \file interact.c
 *  \calculate the interaction of the gas and dust
 */
void interaction(void)
{
    interactionGas();
    interactionDust();
}

#if !defined(MKSKHEMA) && !defined(LPSKHEMA) 	
void interactionGas()
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int i, j, k, n, ndone, maxfill, source;
    int level, ngrp, recvTask, place, nexport;
    MPI_Status status;

    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);
    int NumSphUpdate;
    for (n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
        if ((P[n].Type == 0) && (P[n].Ti_endstep == All.Ti_Current))
            NumSphUpdate++;
    }

    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
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
            if (P[i].Type == 0)
            {
                if (P[i].Ti_endstep == All.Ti_Current)
                {
                    ndone++;
                    for (j = 0; j < NTask; j++)
                        Exportflag[j] = 0;
                    interactGas_evaluate(i, 0);
                    for (j = 0; j < NTask; j++)
                    {
                        if (Exportflag[j])
                        {
                            for (k = 0; k < 3; k++)
                                DragDataIn[nexport].Pos[k] = P[i].Pos[k];
                            DragDataIn[nexport].id = P[i].ID;
                            DragDataIn[nexport].Index = i;
                            DragDataIn[nexport].Task = j;
                            nexport++;
                            nsend_local[j]++;
                        }
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
                interactGas_evaluate(j, 1);

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
                            for (k = 0; k < 3; k++)
                            {
                                SphP[place].GasVelMid[k]
                                    += DragDataPartialResult[source].GasVelMid[k];
                                SphP[place].HydroAccelMid[k]
                                    += DragDataPartialResult[source].HydroAccelMid[k];
                                SphP[place].GasGravMid[k]
                                    += DragDataPartialResult[source].GasGravMid[k];
                                SphP[place].DustVelMid[k]
                                        += DragDataPartialResult[source].DustVelMid[k];
 				SphP[place].DustGravMid[k]
                                        += DragDataPartialResult[source].DustGravMid[k];

                            }
                            SphP[place].GasDensity += DragDataPartialResult[source].GasDensity;
                            SphP[place].SoundSpeedMid
                                += DragDataPartialResult[source].SoundSpeedMid;                           
                            SphP[place].GasNumber += DragDataPartialResult[source].GasNumber;                             
                            SphP[place].DustNumber += DragDataPartialResult[source].DustNumber;
                            SphP[place].DustDensity
                                    += DragDataPartialResult[source].DustDensity;    
                                               
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
    double t_st, Ks, eps, tau;
    double xn[3], yn[3], xn1[3], yn1[3], U_ast[3], V_n1[3];
    for (i = 0; i < N_gas; i++)
        if ((P[i].Ti_endstep == All.Ti_Current) && (P[i].Type == 0))
        {
#ifdef SPH_BND_PARTICLES
            if (P[i].ID == 0)
            {
                SphP[i].DtDragEntropy = 0;
                for (k = 0; k < 3; k++)
                    SphP[i].DragAccel[k] = 0;
            }
            else
            {
#endif
                if (SphP[i].GasNumber > 0)
                {
                    
                    SphP[i].GasDensity /= SphP[i].GasNumber;
                    SphP[i].SoundSpeedMid /= SphP[i].GasNumber;
                    tau = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
                    SphP[i].DtDragEntropy = 0;
                    for (j = 0; j < 3; j++)
                    {
                        SphP[i].DragAccel[j] = 0;
                        SphP[i].GasVelMid[j] /= SphP[i].GasNumber;
                        SphP[i].HydroAccelMid[j] /= SphP[i].GasNumber;
                        SphP[i].GasGravMid[j] /= SphP[i].GasNumber;                        
                    }
                    if (SphP[i].DustNumber > 0)
                    {                        
                        SphP[i].DustDensity /= SphP[i].DustNumber;
                        eps = All.MassTable[1] * SphP[i].DustNumber * 1. / SphP[i].GasNumber/ All.MassTable[0]; 
                        SphP[i].GasDensity=SphP[i].DustDensity/eps;
                        t_st = All.DustRho * All.DustSize/SphP[i].SoundSpeedMid/SphP[i].GasDensity;     
                                                 
                        Ks = SphP[i].DustDensity / t_st;                       
                        for (j = 0; j < 3; j++)
                        {
                            SphP[i].DustVelMid[j] /= SphP[i].DustNumber;
                            SphP[i].DustGravMid[j] /= SphP[i].DustNumber;
                            xn[j] = SphP[i].GasVelMid[j] - SphP[i].DustVelMid[j];
                            yn[j] = SphP[i].GasVelMid[j] + eps * SphP[i].DustVelMid[j];
                            xn1[j] = (xn[j] + (SphP[i].HydroAccelMid[j] +SphP[i].GasGravMid[j]-SphP[i].DustGravMid[j])*tau)
                                / (1. + (eps + 1.) * tau / t_st);
                            yn1[j] = yn[j] + tau * (SphP[i].HydroAccelMid[j] +SphP[i].GasGravMid[j]-eps*SphP[i].DustGravMid[j]);
                            U_ast[j] = (yn1[j] - xn1[j]) / (1. + eps);
                            V_n1[j] = (SphP[i].VelPred[j] + U_ast[j] * tau * eps / t_st
                                          + SphP[i].HydroAccel[j] * tau + P[i].GravAccel[j] * tau)
                                / (1. + eps * tau / t_st);


                            SphP[i].DragAccel[j] -= Ks/SphP[i].GasDensity * (V_n1[j] - U_ast[j]);
#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
                           SphP[i].DtDragEntropy += Ks/SphP[i].GasDensity  * (U_ast[j] - V_n1[j]) * (U_ast[j] - V_n1[j]);
#endif
                        }
                    }

#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
                    SphP[i].DtDragEntropy *= GAMMA_MINUS1 / pow(SphP[i].Density, GAMMA_MINUS1);
                    SphP[i].DtEntropy += SphP[i].DtDragEntropy;
#endif
                }
                else
                    for (j = 0; j < 3; j++)
                        SphP[i].DragAccel[j] = 0;
#ifdef SPH_BND_PARTICLES
            }
#endif
        }
}

void interactGas_evaluate(int target, int mode)
{
    int j, k, n, startnode, numngb;
    FLOAT* pos;
    double V_a[3], H_a[3], rho_j, rho_a, U_j[3], soundspeed_a, G_a[3], G_j[3];
    int Ni, Li, id;
    if (mode == 0)
    {
        pos = P[target].Pos;
        id = P[target].ID;
    }
    else
    {
        id = DragDataGet[target].id;
        pos = DragDataGet[target].Pos;
    }
#ifdef DIM1
    FLOAT X[3];
    int ix[3];
    ix[0] = floor((All.BarrierDistance + pos[0]) / All.DustGasMechStep);
    X[0] = -All.BarrierDistance + (ix[0] + 0.5) * All.DustGasMechStep;
    X[1] = X[2] = 0;
#else
#ifdef DECART
    FLOAT X[3];
    int ix[3];
    for (j = 0; j < 3; j++)
    {
        ix[j] = floor((All.BarrierDistance + pos[j]) / All.DustGasMechStep);
        X[j] = -All.BarrierDistance + (ix[j] + 0.5) * All.DustGasMechStep;
    }
#endif
#endif
    /* initialize variables before SPH loop is started */
    for (j = 0; j < 3; j++)
    {
        V_a[j] = 0;
        H_a[j] = 0;
        U_j[j] = 0;
        G_a[j] = 0;
        G_j[j] = 0;
    }
    soundspeed_a = 0.;
    rho_a = 0.;
    Ni = 0;
    Li = 0;
    rho_j = 0;
    startnode = All.MaxPart;
    do
    {
#ifdef DIM1
        numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 1);
#else
#ifdef DECART
        numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 1);
#else
        numngb = ngb_treefind_rad(&pos[0], &startnode, 1);
#endif
#endif

        for (n = 0; n < numngb; n++)
        {
            j = Ngblist[n];
            rho_j += SphP[j].Density;
            for (k = 0; k < 3; k++)
            {
                U_j[k] += SphP[j].VelPred[k];
                G_j[k] += P[j].GravAccel[k];
            }
        }
        Li = numngb;
    } while (startnode >= 0);
    /* Now start the actual SPH computation for this particle */
    if (Li > 0)
    {
        startnode = All.MaxPart;
        do
        {
#ifdef DIM1
            numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 0);
#else
#ifdef DECART
            numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 0);
#else
            numngb = ngb_treefind_rad(&pos[0], &startnode, 0);
#endif
#endif
            for (n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                rho_a += SphP[j].Density;
                soundspeed_a += sqrt(GAMMA * SphP[j].Pressure / SphP[j].Density);
                for (k = 0; k < 3; k++)
                {
                    V_a[k] += SphP[j].VelPred[k];
                    H_a[k] += SphP[j].HydroAccel[k];
                    G_a[k] += P[j].GravAccel[k];
                }
            }
            Ni = numngb;
        } while (startnode >= 0);
    }
    if (mode == 0)
    {
        for (k = 0; k < 3; k++)
        {
            SphP[target].GasVelMid[k] = V_a[k];
            SphP[target].HydroAccelMid[k] = H_a[k];
            SphP[target].GasGravMid[k] = G_a[k];
            SphP[target].DustVelMid[k] = U_j[k];
            SphP[target].DustGravMid[k] = G_j[k];
        }
        SphP[target].GasDensity = rho_a;
        SphP[target].GasNumber = Ni;
        SphP[target].DustDensity = rho_j;
        SphP[target].DustNumber = Li;
        SphP[target].SoundSpeedMid = soundspeed_a;
    }
    else
    {
        for (k = 0; k < 3; k++)
        {
            DragDataResult[target].GasVelMid[k] = V_a[k];
            DragDataResult[target].GasGravMid[k] = G_a[k];
            DragDataResult[target].HydroAccelMid[k] = H_a[k];
            DragDataResult[target].DustVelMid[k] = U_j[k];
            DragDataResult[target].DustGravMid[k] = G_j[k];
        }
        DragDataResult[target].GasDensity = rho_a;
        DragDataResult[target].GasNumber = Ni;
        DragDataResult[target].DustDensity = rho_j;
        DragDataResult[target].DustNumber = Li;
        DragDataResult[target].SoundSpeedMid = soundspeed_a;
    }
}
void interactionDust()
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int i, j, k, n, ndone, maxfill, source;
    int level, ngrp, recvTask, place, nexport;
    MPI_Status status;

    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);
    int NumSphUpdate;
    for (n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
        if ((P[n].Type == 1) && (P[n].Ti_endstep == All.Ti_Current))
            NumSphUpdate++;
    }

    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
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
            if (P[i].Type == 1)
            {
                if (P[i].Ti_endstep == All.Ti_Current)
                {
                    ndone++;
                    for (j = 0; j < NTask; j++)
                        Exportflag[j] = 0;
                    interactDust_evaluate(i, 0);
                    for (j = 0; j < NTask; j++)
                    {
                        if (Exportflag[j])
                        {
                            for (k = 0; k < 3; k++)
                                DragDataIn[nexport].Pos[k] = P[i].Pos[k];
                            DragDataIn[nexport].Index = i;
                            DragDataIn[nexport].Task = j;
                            nexport++;
                            nsend_local[j]++;
                        }
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
                interactDust_evaluate(j, 1);

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

                // sendTask = ThisTask;
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
                            for (k = 0; k < 3; k++)
                            {
                                SphP[place].GasVelMid[k]
                                    += DragDataPartialResult[source].GasVelMid[k];                                
                                SphP[place].DustVelMid[k]
                                    += DragDataPartialResult[source].DustVelMid[k];
                                SphP[place].HydroAccelMid[k]
                                    += DragDataPartialResult[source].HydroAccelMid[k];
				SphP[place].GasGravMid[k]
                                    += DragDataPartialResult[source].GasGravMid[k];                                
                                SphP[place].DustGravMid[k]
                                    += DragDataPartialResult[source].DustGravMid[k];
                            }
                            SphP[place].GasDensity += DragDataPartialResult[source].GasDensity;
                            SphP[place].GasNumber += DragDataPartialResult[source].GasNumber;
                            SphP[place].SoundSpeedMid
                                += DragDataPartialResult[source].SoundSpeedMid;
                            SphP[place].DustNumber += DragDataPartialResult[source].DustNumber;
                            SphP[place].DustDensity += DragDataPartialResult[source].DustDensity;			                              
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
    double t_st, Ks, eps, tau, xn[3], yn[3], xn1[3], yn1[3], V_ast[3], U_n1[3];
    for (i = 0; i < N_gas; i++)
        if ((P[i].Ti_endstep == All.Ti_Current) && (P[i].Type == 1))
        {
#ifdef SPH_BND_PARTICLES
            if (P[i].ID == 0)
                for (k = 0; k < 3; k++)
                    SphP[i].DragAccel[k] = 0;
            else
            {
#endif
                if ((SphP[i].GasNumber > 0) && (SphP[i].DustNumber > 0))
                {
                    SphP[i].SoundSpeedMid /= SphP[i].GasNumber;
                    SphP[i].GasDensity /= SphP[i].GasNumber;
                    SphP[i].DustDensity /= SphP[i].DustNumber;
                    t_st = All.DustRho * P[i].Radius/ SphP[i].SoundSpeedMid/SphP[i].GasDensity;
                    eps = All.MassTable[1] * SphP[i].DustNumber * 1. / SphP[i].GasNumber/ All.MassTable[0]; 
 		    tau = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;                          
                    Ks = SphP[i].DustDensity / t_st;                   
                    for (j = 0; j < 3; j++)
                    {
                        SphP[i].GasVelMid[j] /= SphP[i].GasNumber;
                        SphP[i].HydroAccelMid[j] /= SphP[i].GasNumber;
                        SphP[i].GasGravMid[j] /= SphP[i].GasNumber;
                        SphP[i].DustVelMid[j] /= SphP[i].DustNumber;
			SphP[i].DustGravMid[j] /= SphP[i].DustNumber;
                        xn[j] = SphP[i].GasVelMid[j] - SphP[i].DustVelMid[j];
                        yn[j] = SphP[i].GasVelMid[j] + eps * SphP[i].DustVelMid[j];
                         xn1[j] = (xn[j] + (SphP[i].HydroAccelMid[j] +SphP[i].GasGravMid[j]-SphP[i].DustGravMid[j])*tau)
                                / (1. + (eps + 1.) * tau / t_st);
                            yn1[j] = yn[j] + tau * (SphP[i].HydroAccelMid[j] +SphP[i].GasGravMid[j]-eps*SphP[i].DustGravMid[j]);
                        V_ast[j] = (yn1[j] + eps * xn1[j]) / (1. + eps);
                        U_n1[j] = (SphP[i].VelPred[j] + tau * V_ast[j] / t_st + P[i].GravAccel[j]*tau)
                            / (1. + tau / t_st);
                        SphP[i].DragAccel[j] = Ks/SphP[i].DustDensity * (V_ast[j] - U_n1[j]);
                    }
                }
                else
                    for (j = 0; j < 3; j++)
                        SphP[i].DragAccel[j] = 0;
#ifdef SPH_BND_PARTICLES
            }
#endif
        }
}

void interactDust_evaluate(int target, int mode)
{
    int j, k, n, startnode, numngb;
    FLOAT* pos;
    double rho_a, V_a[3], H_a[3], U_j[3], rho_j, soundspeed_a, G_a[3], G_j[3];
    ;
    int Ni, Li;
    if (mode == 0)
    {
        pos = P[target].Pos;
    }
    else
    {
        pos = DragDataGet[target].Pos;
    }
/* initialize variables before SPH loop is started */
#ifdef DIM1
    FLOAT X[3];
    int ix[3];
    ix[0] = floor((All.BarrierDistance + pos[0]) / All.DustGasMechStep);
    X[0] = -All.BarrierDistance + (ix[0] + 0.5) * All.DustGasMechStep;
    X[1] = X[2] = 0;
#else
#ifdef DECART
    FLOAT X[3];
    int ix[3];
    for (j = 0; j < 3; j++)
    {
        ix[j] = floor((All.BarrierDistance + pos[j]) / All.DustGasMechStep);
        X[j] = -All.BarrierDistance + (ix[j] + 0.5) * All.DustGasMechStep;
    }
#endif
#endif

    rho_a = 0;
    rho_j = 0;
    soundspeed_a = 0;
    Ni = Li = 0;
    for (j = 0; j < 3; j++)
    {
        V_a[j] = 0;
        H_a[j] = 0;
        U_j[j] = 0;
        G_a[j] = 0;
        G_j[j] = 0;
    }
    startnode = All.MaxPart;
    do
    {
#ifdef DIM1
        numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 0);
#else
#ifdef DECART
        numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 0);
#else
        numngb = ngb_treefind_rad(&pos[0], &startnode, 0);
#endif
#endif
        for (n = 0; n < numngb; n++)
        {
            j = Ngblist[n];
            soundspeed_a += sqrt(GAMMA * SphP[j].Pressure / SphP[j].Density);
            rho_a += SphP[j].Density;
            for (k = 0; k < 3; k++)
            {
                V_a[k] += SphP[j].VelPred[k];
                H_a[k] += SphP[j].HydroAccel[k];
                G_a[k] += P[j].GravAccel[k];
            }
        }
        Ni += n;
    } while (startnode >= 0);
    /* Now start the actual SPH computation for this particle */
    if (Ni > 0)
    {
        startnode = All.MaxPart;
        do
        {
#ifdef DIM1
            numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 1);
#else
#ifdef DECART
            numngb = ngb_treefind_cell(&X[0], 0.5 * All.DustGasMechStep, &startnode, 1);
#else
            numngb = ngb_treefind_rad(&pos[0], &startnode, 1);
#endif
#endif
            for (n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                rho_j += SphP[j].Density;
                for (k = 0; k < 3; k++)
                {
                    U_j[k] += SphP[j].VelPred[k];
                    G_j[k] += P[j].GravAccel[k];
                }
            }
            Li += n;
        } while (startnode >= 0);
    }
    /* Now collect the result at the right place */
    if (mode == 0)
    {
        for (k = 0; k < 3; k++)
        {
            SphP[target].GasVelMid[k] = V_a[k];
            SphP[target].GasGravMid[k] = G_a[k];
            SphP[target].DustVelMid[k] = U_j[k];
            SphP[target].DustGravMid[k] = G_j[k];
            SphP[target].HydroAccelMid[k] = H_a[k];
        }
        SphP[target].GasDensity = rho_a;
        SphP[target].DustDensity = rho_j;
        SphP[target].SoundSpeedMid = soundspeed_a;
        SphP[target].DustNumber = Li;
        SphP[target].GasNumber = Ni;
    }
    else
    {
        for (k = 0; k < 3; k++)
        {
            DragDataResult[target].GasVelMid[k] = V_a[k];
            DragDataResult[target].GasGravMid[k] = G_a[k];
            DragDataResult[target].DustVelMid[k] = U_j[k];
            DragDataResult[target].DustGravMid[k] = G_j[k];
            DragDataResult[target].HydroAccelMid[k] = H_a[k];
        }
        DragDataResult[target].GasDensity = rho_a;
        DragDataResult[target].DustDensity = rho_j;
        DragDataResult[target].SoundSpeedMid = soundspeed_a;
        DragDataResult[target].DustNumber = Li;
        DragDataResult[target].GasNumber = Ni;
    }
}
#endif
#ifdef LPSKHEMA	

void interactionGas(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, k, n, ndone, maxfill, source;
  int level, ngrp, recvTask, place, nexport;
  MPI_Status status;

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);
  int NumSphUpdate;
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {    
     if((P[n].Type==0)&&(P[n].Ti_endstep == All.Ti_Current))
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);
  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
    i=0;		/* beginn with this index */
    ntotleft = ntot;		/* particles left for all tasks together */
    while(ntotleft > 0) 
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

	  /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
      {
       if(P[i].Type==0)
       {   	  
	if(P[i].Ti_endstep == All.Ti_Current)
	 {
	  ndone++;          
	  for(j = 0; j < NTask; j++)
	    Exportflag[j] = 0;
       	  interactGas_evaluate(i, 0);

	  for(j = 0; j < NTask; j++)
	  {
	    if(Exportflag[j])
	      {
	 	for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;
		   // HydroDataIn[nexport].Mass = P[i].Mass;	   
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;	   
		    HydroDataIn[nexport].Density = SphP[i].Density;	          
		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;              
		    nexport++;
		    nsend_local[j]++;		  
	      }
	    }
          }
         }
        }

	qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

	for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];


	MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  /* now do the particles that need to be exported */
        
	for(level = 1; level < (1 << PTask); level++)
	  {
	  
	    for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	    for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		 if((j ^ ngrp) < NTask)
		 if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		    maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
		if(maxfill >= All.BunchSizeDensity)
		  break;
		recvTask = ThisTask ^ ngrp;

  	        if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		  }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	     }
                 
	
	     for(j = 0; j < nbuffer[ThisTask]; j++)
	       interactGas_evaluate(j, 1);
	
           
	     MPI_Barrier(MPI_COMM_WORLD);
	

	      /* get the result */

	     for(j = 0; j < NTask; j++)
	      nbuffer[j] = 0;
	     for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	      {
		maxfill = 0;
		for(j = 0; j < NTask; j++)
		 {
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		 }
		 if(maxfill >= All.BunchSizeDensity)
		   break;

		 recvTask = ThisTask ^ ngrp;

		 if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* send the results */
		     MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
		    for(j = 0; j < nsend_local[recvTask]; j++)
		      {
		       source = j + noffset[recvTask];
		       place = HydroDataIn[source].Index;                                                               
		       for(k = 0; k < 3; k++)
                       {     			
                        SphP[place].DragAccel[k] += HydroDataPartialResult[source].AccDrag[k];
                       }
#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
    	               SphP[place].DtDragEntropy+=HydroDataPartialResult[source].DtEntropy;    
#endif                         
		      } 
		    }
		  }
             

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	      }
           
	

	      level = ngrp - 1;
	   }
          
	   MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
         
	   for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}

     free(ndonelist);
     free(nsend);
     free(nsend_local);
     free(nbuffer);
     free(noffset);
#if defined(SPH_BND_PARTICLES) || !defined(SAVETEMP_LTR) || !defined(ISOTHERM_EQS)
    for(i = 0; i < N_gas; i++)
     if(P[i].Ti_endstep == All.Ti_Current)
      {     
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].DtDragEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].DragAccel[k] = 0;
	  }
#endif
#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
      SphP[i].DtDragEntropy *= GAMMA_MINUS1 / pow(SphP[i].Density, GAMMA_MINUS1);
      SphP[i].DtEntropy+=SphP[i].DtDragEntropy;
#endif
      }
#endif
}



void interactGas_evaluate(int target, int mode)
{
  int j, k, n, startnode, numngb;
  FLOAT *pos, *vel;
  FLOAT h_i;
  double acc[3];
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv3;
  double  rho,soundspeed_i,pressure,DtEntropy;
  double hfc, vdotr;
  double wi,r, r2, u,Ks;

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;    
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;     
      soundspeed_i = sqrt(GAMMA * pressure / rho);
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;       
      rho=HydroDataGet[target].Density;  
      pressure=HydroDataGet[target].Pressure;
      soundspeed_i = sqrt(GAMMA * pressure / rho);     
    }

  DtEntropy=0;
  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = 0; 
  h_i2 = h_i * h_i;
  hinv = 1.0 / h_i;
#ifdef DIM1
    hinv3 = hinv;
#else
    hinv3 = hinv * hinv * hinv;
#endif
  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;

  Ks=NUMDIMS*sqrt(8./M_PI/GAMMA)/(All.DustSize*All.DustRho)*soundspeed_i; 
  do
    {
      numngb = ngb_treefind_variable(&pos[0], h_i, &startnode,1);
      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  dx = P[j].Pos[0]-pos[0];
	  dy = P[j].Pos[1]-pos[1];
	  dz = P[j].Pos[2]-pos[2];
          
	  r2 = dx * dx + dy * dy + dz * dz;	 
	  if(r2 < h_i2)
	    {
	      r = sqrt(r2);
	      if(r > 0)
		{		  
		  dvx =  SphP[j].VelPred[0]-vel[0];
		  dvy =  SphP[j].VelPred[1]-vel[1] ;
		  dvz =  SphP[j].VelPred[2]-vel[2];
		  vdotr = dx * dvx + dy * dvy + dz * dvz;                                               
		  u = r * hinv;
	          if(u < 0.5)
		  {
		   wi = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u)*u*u*KERNEL_COEFF_7;
		  }
	          else
		  {
		   wi = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u)*u*u*KERNEL_COEFF_7;
                  }                                                   		 		                  
                  hfc = P[j].Mass*Ks*vdotr*wi/r/r;//SphP[j].Ni;             	  
		  acc[0] += hfc*dx;
		  acc[1] += hfc*dy;
		  acc[2] += hfc*dz;    
                  DtEntropy +=hfc*vdotr;                                        
                 }
                }
	   }	
     }
    while(startnode >= 0); 
  /* Now collect the result at the right place */
   if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].DragAccel[k] = acc[k];    
      SphP[target].DtDragEntropy = DtEntropy;               
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].AccDrag[k] = acc[k];  
      HydroDataResult[target].DtEntropy = DtEntropy;         
    }         		
}

void interactionDust(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, k, n, ndone, maxfill, source;
  int level, ngrp, recvTask, place, nexport;
  MPI_Status status;

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);
  int NumSphUpdate;
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {    
     if((P[n].Type==1)&&(P[n].Ti_endstep == All.Ti_Current))
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);



  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
    i=0;		/* beginn with this index */
    ntotleft = ntot;		/* particles left for all tasks together */
    while(ntotleft > 0) 
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

	  /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
      {
       if(P[i].Type==1)
       {   	  
	if(P[i].Ti_endstep == All.Ti_Current)
	 {
	  ndone++;
       
	  for(j = 0; j < NTask; j++)
	    Exportflag[j] = 0;
       	  interactDust_evaluate(i, 0);

	  for(j = 0; j < NTask; j++)
	  {
	    if(Exportflag[j])
	      {
	 	for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;
		    HydroDataIn[nexport].id = P[i].ID;		   		                     
		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
		    nexport++;
		    nsend_local[j]++;		  
	      }
	    }
          }
         }
        }
	qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

	for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  /* now do the particles that need to be exported */
        
	for(level = 1; level < (1 << PTask); level++)
	  {

	    for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	    for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		 if((j ^ ngrp) < NTask)
		 if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		    maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
		if(maxfill >= All.BunchSizeDensity)
		  break;
		recvTask = ThisTask ^ ngrp;

  	        if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		  }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	     }
                 
	
	     for(j = 0; j < nbuffer[ThisTask]; j++)
	       interactDust_evaluate(j, 1);
        
	     MPI_Barrier(MPI_COMM_WORLD);
	  
	     for(j = 0; j < NTask; j++)
	      nbuffer[j] = 0;
	     for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	      {
		maxfill = 0;
		for(j = 0; j < NTask; j++)
		 {
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		 }
		 if(maxfill >= All.BunchSizeDensity)
		   break;

		 recvTask = ThisTask ^ ngrp;

		 if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* send the results */
		     MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
		    for(j = 0; j < nsend_local[recvTask]; j++)
		      {
		       source = j + noffset[recvTask];
		       place = HydroDataIn[source].Index;                                                                                  
		       for(k = 0; k < 3; k++)
                       {     			
                        SphP[place].DragAccel[k] +=  HydroDataPartialResult[source].AccDrag[k];
                       }
    	                     
		      } 
		    }
		  }
             

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	      }
           
	      level = ngrp - 1;
   }
          
	   MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
         
	   for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}
//for(i=0;i<N_gas;i++)
//if(P[i].ID==4924)
	//   printf("%1.16f %1.16f %1.16f %1.16f\n",SphP[i].DragAccel[0],P[i].Pos[0],P[i].Vel[0],SphP[i].Density);
     free(ndonelist);
     free(nsend);
     free(nsend_local);
     free(nbuffer);
     free(noffset);    
}



void interactDust_evaluate(int target, int mode)
{
  int j, k, n, startnode, numngb;//,id;
  FLOAT *pos, *vel;
  FLOAT h_i;
  double acc[3];
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv3;
  double  soundspeed_j;
  double hfc, vdotr;
  double h_j,wi, r, r2, u,Ks;

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;      
     
      //id=P[target].ID;
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;          
   
     // id=HydroDataGet[target].id;
    }


  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = 0; 
  h_i2 = h_i * h_i;
  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_variable(&pos[0], h_i, &startnode,0);
      //printf("numngb=%i\n",numngb);  
      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  dx =  pos[0]-P[j].Pos[0];
	  dy =  pos[1]-P[j].Pos[1];
	  dz =  pos[2]-P[j].Pos[2];
          
	  r2 = dx * dx + dy * dy + dz * dz;
	  h_j = SphP[j].Hsml;
	  if(r2 < h_i2)
	    {
	      r = sqrt(r2);
	      if((r > 0))
		{                 		  
		  soundspeed_j = sqrt(GAMMA * SphP[j].Pressure / SphP[j].Density);
		  dvx =  vel[0]-SphP[j].VelPred[0];
		  dvy =  vel[1]-SphP[j].VelPred[1];
		  dvz =  vel[2]-SphP[j].VelPred[2];
		  vdotr = dx * dvx + dy * dvy + dz * dvz;

                 if(r < h_j)
		    {
		      hinv = 1.0 / h_j;
#ifdef DIM1
		     hinv3 = hinv;
#else
  		     hinv3 = hinv * hinv * hinv;
#endif
		      u = r * hinv;
		      if(u < 0.5)
		      {
		       wi = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u)*u*u*KERNEL_COEFF_7;
		      }
	              else
		      {
		       wi = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u)*u*u*KERNEL_COEFF_7;
                      }   
		    }
		  else
		    {
		      wi = 0;
		    }

                                 
                   Ks=NUMDIMS*sqrt(8./M_PI/GAMMA)/(All.DustSize*All.DustRho)*soundspeed_j;    
         
                  hfc = P[j].Mass*Ks*vdotr*wi/r/r;	  
		  acc[0] -= hfc*dx;
		  acc[1] -= hfc*dy;
		  acc[2] -= hfc*dz;                                                             
		}
	    }
	}
    }
  while(startnode >= 0);

  /* Now collect the result at the right place */
 if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].DragAccel[k] = acc[k];                   
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].AccDrag[k] = acc[k];           
    }      
}
#endif

#ifdef MKSKHEMA

void interactionGas(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, k, n, ndone, maxfill, source;
  int level, ngrp, recvTask, place, nexport; 
  MPI_Status status;

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);
  int NumSphUpdate;
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
  {    
   if((P[n].Type==0)&&(P[n].Ti_endstep == All.Ti_Current))
   NumSphUpdate++;
  }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
  i=0;		/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */
  while(ntotleft > 0) 
  {
   for(j = 0; j < NTask; j++)
    nsend_local[j] = 0;
	  /* do local particles and prepare export list */
      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
      {
       if(P[i].Type==0)
       {   	  
	if(P[i].Ti_endstep == All.Ti_Current)
	 {
	  ndone++;          
	  for(j = 0; j < NTask; j++)
	    Exportflag[j] = 0;
       	  interactGas_evaluate(i, 0);

	  for(j = 0; j < NTask; j++)
	  {
	    if(Exportflag[j])
	      {
	 	for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;		  
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;	   
		    HydroDataIn[nexport].Density = SphP[i].Density;	          		 		          	   
		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;              
		    nexport++;
		    nsend_local[j]++;		  
	      }
	    }
          }
         }
        }

	qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

	for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];


	MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  /* now do the particles that need to be exported */
        
	for(level = 1; level < (1 << PTask); level++)
	  {
	  
	    for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	    for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		 if((j ^ ngrp) < NTask)
		 if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		    maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
		if(maxfill >= All.BunchSizeDensity)
		  break;		
		recvTask = ThisTask ^ ngrp;

  	        if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		  }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	     }
                 
	
	     for(j = 0; j < nbuffer[ThisTask]; j++)
	       interactGas_evaluate(j, 1);
	
           
	     MPI_Barrier(MPI_COMM_WORLD);
	

	      /* get the result */

	     for(j = 0; j < NTask; j++)
	      nbuffer[j] = 0;
	     for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	      {
		maxfill = 0;
		for(j = 0; j < NTask; j++)
		 {
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		 }
		 if(maxfill >= All.BunchSizeDensity)
		   break;		
		 recvTask = ThisTask ^ ngrp;

		 if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* send the results */
		     MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
		    for(j = 0; j < nsend_local[recvTask]; j++)
		      {
		       source = j + noffset[recvTask];
		       place = HydroDataIn[source].Index;                                                               
		       for(k = 0; k < 3; k++)
                       {     			
                        SphP[place].DragAccel[k] += HydroDataPartialResult[source].Acc[k];
                       }
#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
    	               SphP[place].DtDragEntropy+=HydroDataPartialResult[source].DtEntropy;    
#endif                         
		      } 
		    }
		  }
             

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	      }          
	      level = ngrp - 1;
	   }
          
	   MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
         
	   for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}

     free(ndonelist);
     free(nsend);
     free(nsend_local);
     free(nbuffer);
     free(noffset);
#if defined(SPH_BND_PARTICLES) || !defined(SAVETEMP_LTR) || !defined(ISOTHERM_EQS)
    for(i = 0; i < N_gas; i++)
     if(P[i].Ti_endstep == All.Ti_Current)
      {     
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].DtDragEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].DragAccel[k] = 0;
	  }
#endif
#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
      SphP[i].DtDragEntropy *= GAMMA_MINUS1 / pow(SphP[i].Density, GAMMA_MINUS1);
      SphP[i].DtEntropy+=SphP[i].DtDragEntropy;
#endif
      }
#endif
}



void interactGas_evaluate(int target, int mode)
{
  int j, k, n, startnode, numngb;
  FLOAT *pos, *vel;
  FLOAT h_i,h_i2,hinv,hinv3;;
  double acc[3],dpos[3],dvel[3];
  double  rho,soundspeed_i,pressure;
  double  Ka[3],Ke,K,dtEntropy,r2,u,r,wk,vdotr;
  int Ni;
  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;      
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;     
      soundspeed_i = sqrt(GAMMA * pressure / rho);     
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;          
      rho=HydroDataGet[target].Density;  
      pressure=HydroDataGet[target].Pressure;
      soundspeed_i = sqrt(GAMMA * pressure / rho);     
    } 
  h_i2 = h_i* h_i;
  hinv = 1.0 / h_i;
#ifdef DIM1
    hinv3 = hinv;
#else
    hinv3 = hinv * hinv * hinv;
#endif
  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = 0; 
   dtEntropy=0; 
  Ni=0;
  for(j=0;j<3;j++)
   Ka[j]=0;
  Ke=0;
  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart; 
  do
  {
   numngb = ngb_treefind_variable(&pos[0], h_i, &startnode,1);
   for(n = 0; n < numngb; n++)
   {
    j = Ngblist[n];	
    r2=0;
    for(k=0;k<3;k++)
    {
     dpos[k] =  P[j].Pos[k]-pos[k];              
     r2+= dpos[k] * dpos[k];
    }
    if(r2 < h_i2)
    {
     Ni++;	
     r = sqrt(r2);
     u = r * hinv;
     if(u < 0.5)
     {
       wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
     }
     else
     {
       wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);		
     }      		  
     vdotr=0;
     for(k=0;k<3;k++) 
     {       
       dvel[k] =  SphP[j].VelPred[k]-vel[k];
       vdotr+=dvel[k]*dpos[k];
     }
     K=P[j].Mass*soundspeed_i/(All.DustRho * All.DustSize)*vdotr/(r2+0.001*h_i2)*wk;
     for(k=0;k<3;k++) 
      Ka[k]+=K*dpos[k];
     Ke+=K*vdotr;      				                                                         
    }                                                       
   }
  }
  while(startnode >= 0);
  if(Ni>0)
  {                  
   for(k = 0; k < 3; k++) 
    acc[k]=SIGMA*Ka[k]; 
   dtEntropy+=SIGMA*Ke;
  }    
   /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].DragAccel[k] = acc[k];    
      SphP[target].DtDragEntropy = dtEntropy;               
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];  
      HydroDataResult[target].DtEntropy = dtEntropy;         
    }         		
}

void interactionDust(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, k, n, ndone, maxfill, source;
  int level, ngrp, recvTask, place, nexport;

  MPI_Status status;

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);
  int NumSphUpdate;
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {    
     if((P[n].Type==1)&&(P[n].Ti_endstep == All.Ti_Current))
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
    i=0;		/* beginn with this index */
    ntotleft = ntot;		/* particles left for all tasks together */
    while(ntotleft > 0) 
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

	  /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
      {
       if(P[i].Type==1)
       {   	  
	if(P[i].Ti_endstep == All.Ti_Current)
	 {
	  ndone++;
       
	  for(j = 0; j < NTask; j++)
	    Exportflag[j] = 0;
       	  interactDust_evaluate(i, 0);

	  for(j = 0; j < NTask; j++)
	  {
	    if(Exportflag[j])
	      {
	 	for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;		   	 		  
		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
		    nexport++;
		    nsend_local[j]++;		  
	      }
	    }
          }
         }
        }
	qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

	for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  /* now do the particles that need to be exported */
        
	for(level = 1; level < (1 << PTask); level++)
	  {

	    for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	    for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		 if((j ^ ngrp) < NTask)
		 if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		    maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
		if(maxfill >= All.BunchSizeDensity)
		  break;
		recvTask = ThisTask ^ ngrp;

  	        if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		  }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	     }
                 
	
	     for(j = 0; j < nbuffer[ThisTask]; j++)
	       interactDust_evaluate(j, 1);
        
	     MPI_Barrier(MPI_COMM_WORLD);
	  
	     for(j = 0; j < NTask; j++)
	      nbuffer[j] = 0;
	     for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	      {
		maxfill = 0;
		for(j = 0; j < NTask; j++)
		 {
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		 }
		 if(maxfill >= All.BunchSizeDensity)
		   break;
		 recvTask = ThisTask ^ ngrp;

		 if(recvTask < NTask)
		  {
		   if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			  /* send the results */
		     MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
		    for(j = 0; j < nsend_local[recvTask]; j++)
		      {
		       source = j + noffset[recvTask];
		       place = HydroDataIn[source].Index;                                                                                
		       for(k = 0; k < 3; k++)
                       {     			
                        SphP[place].DragAccel[k] += HydroDataPartialResult[source].Acc[k];
                       }
    	                     
		      } 
		    }
		  }
             

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	      }
           
	      level = ngrp - 1;
	   }
          
	   MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
         
	   for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}
     free(ndonelist);
     free(nsend);
     free(nsend_local);
     free(nbuffer);
     free(noffset);
}



void interactDust_evaluate(int target, int mode)
{
  int j, k, n, startnode, numngb;
  FLOAT *pos, *vel;
  FLOAT h_i,h_i2,hinv,hinv3;
  double acc[3];
  double  soundspeed_j,dvel[3],dpos[3],K[3];
  double  Ni,r,r2,u,wk,vdotr;
  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;      
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;            
    }
 // h_i=All.DustGasMechStep;
  h_i2 = h_i* h_i;
  hinv = 1.0 / h_i;
#ifdef DIM1
    hinv3 = hinv;
#else
    hinv3 = hinv * hinv * hinv;
#endif
  /* initialize variables before SPH loop is started */  
  acc[0] = acc[1] = acc[2] = 0;   
  
  soundspeed_j=0;
  Ni=0;
  for(j=0;j<3;j++)
   K[j]=0;
  Ni=0;
  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_variable(&pos[0], h_i, &startnode,0);
      for(n = 0; n < numngb; n++)
      {       
       j = Ngblist[n];
       r2=0;
       for(k=0;k<3;k++)
       {
        dpos[k] =  pos[k]-P[j].Pos[k];              
        r2+= dpos[k] * dpos[k];
       }
       if(r2 < h_i2)
       {
          Ni++;	
          r = sqrt(r2);
	  u = r * hinv;
	  if(u < 0.5)
	  {
	   wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
	  }
	  else
	  {
	   wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);		
	  }      		  
	  soundspeed_j = sqrt(GAMMA * SphP[j].Pressure / SphP[j].Density);  
          vdotr=0;
          for(k=0;k<3;k++) 
          {       
	   dvel[k] = vel[k]-SphP[j].VelPred[k];
           vdotr+=dvel[k]*dpos[k];
	  }
          for(k=0;k<3;k++)
           K[k]-=P[j].Mass*soundspeed_j/(All.DustRho * All.DustSize)*vdotr/(r2+0.001*h_i2)*dpos[k]*wk;				                                                         
	}
    }
  }
  while(startnode >= 0);
  if(Ni>0)
  {                      
   for(k = 0; k < 3; k++) 
    acc[k]=SIGMA*K[k]; 
   }    

   /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].DragAccel[k] = acc[k];                   
    }
  els
