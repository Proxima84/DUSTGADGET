#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */


/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void density(void)
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int i, j, n, ndone, npleft, maxfill, source, iter = 0;
    int level, ngrp, recvTask, place, nexport;
    double dt_entr;
    MPI_Status status;
    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);

    for (n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
        SphP[n].Left = SphP[n].Right = 0;

        if (P[n].Ti_endstep == All.Ti_Current)
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
    do
    {
        i = 0; /* beginn with this index */
        ntotleft = ntot; /* particles left for all tasks together */

        while (ntotleft > 0)
        {
            for (j = 0; j < NTask; j++)
                nsend_local[j] = 0;

            /* do local particles and prepare export list */
            for (nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
                if (P[i].Ti_endstep == All.Ti_Current)
                {
                    ndone++;

                    for (j = 0; j < NTask; j++)
                        Exportflag[j] = 0;

                    density_evaluate(i, 0);

                    for (j = 0; j < NTask; j++)
                    {
                        if (Exportflag[j])
                        {
                            DensDataIn[nexport].Pos[0] = P[i].Pos[0];
                            DensDataIn[nexport].Pos[1] = P[i].Pos[1];
                            DensDataIn[nexport].Pos[2] = P[i].Pos[2];
                            DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
                            DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
                            DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
                            DensDataIn[nexport].Hsml = SphP[i].Hsml;
                            DensDataIn[nexport].Type = P[i].Type;  
                            DensDataIn[nexport].Index = i;
                            DensDataIn[nexport].Task = j;
                            nexport++;
                            nsend_local[j]++;
                        }
                    }
                }

            qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

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
                            MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
                                nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
                                recvTask, TAG_DENS_A, &DensDataGet[nbuffer[ThisTask]],
                                nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
                                MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
                        }
                    }

                    for (j = 0; j < NTask; j++)
                        if ((j ^ ngrp) < NTask)
                            nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
                }
                for (j = 0; j < nbuffer[ThisTask]; j++)
                    density_evaluate(j, 1);
                MPI_Barrier(MPI_COMM_WORLD);
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
                            MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
                                nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
                                MPI_BYTE, recvTask, TAG_DENS_B,
                                &DensDataPartialResult[noffset[recvTask]],
                                nsend_local[recvTask] * sizeof(struct densdata_out), MPI_BYTE,
                                recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

                            /* add the result to the particles */
                            for (j = 0; j < nsend_local[recvTask]; j++)
                            {
                                source = j + noffset[recvTask];
                                place = DensDataIn[source].Index;

                                SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
                                SphP[place].Density += DensDataPartialResult[source].Rho;
                                SphP[place].DivVel += DensDataPartialResult[source].Div;
#ifndef HSMLCONSTANT
                                SphP[place].DhsmlDensityFactor
                                    += DensDataPartialResult[source].DhsmlDensity;
#endif
                                SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
                                SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
                                SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
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

        /* do final operations on results */

        for (i = 0, npleft = 0; i < N_gas; i++)
        {
            if (P[i].Ti_endstep == All.Ti_Current)
            {
                if (P[i].Type == 0)
                {
                    SphP[i].CurlVel
                        = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] + SphP[i].Rot[1] * SphP[i].Rot[1]
                              + SphP[i].Rot[2] * SphP[i].Rot[2])
                        / SphP[i].Density;

                    SphP[i].DivVel /= SphP[i].Density;

                    dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2)
                        * All.Timebase_interval;

                    SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr)
                        * pow(SphP[i].Density, GAMMA);
                }

/* now check whether we had enough neighbours */
#ifndef HSMLCONSTANT
                SphP[i].DhsmlDensityFactor = 1
                    / (1 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));
                if (SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation)
                    || (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
                           && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
                {
                    /* need to redo this particle */
                    npleft++;

                    if (SphP[i].Left > 0 && SphP[i].Right > 0)
                        if ((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].Ti_endstep = -P[i].Ti_endstep - 1; /* Mark as inactive */
                            continue;
                        }

                    if (SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
                        SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
                    else
                    {
                        if (SphP[i].Right != 0)
                        {
                            if (SphP[i].Hsml < SphP[i].Right)
                                SphP[i].Right = SphP[i].Hsml;
                        }
                        else
                            SphP[i].Right = SphP[i].Hsml;
                    }

                    if (SphP[i].Right > 0 && SphP[i].Left > 0)
                        SphP[i].Hsml
                            = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 1.0 / 3);
                    else
                    {
                        if (SphP[i].Right == 0 && SphP[i].Left == 0)
                        {
                            endrun(8188); /* can't occur */
                        }

                        if (SphP[i].Right == 0 && SphP[i].Left > 0)
                        {
                            if (P[i].Type < 2
                                && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
                            {
                                SphP[i].Hsml *= 1
                                    - (SphP[i].NumNgb - All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb)
                                        * SphP[i].DhsmlDensityFactor;
                            }
                            else
                                SphP[i].Hsml *= 1.26;
                        }

                        if (SphP[i].Right > 0 && SphP[i].Left == 0)
                        {
                            if (P[i].Type < 2
                                && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
                            {
                                SphP[i].Hsml *= 1
                                    - (SphP[i].NumNgb - All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb)
                                        * SphP[i].DhsmlDensityFactor;
                            }
                            else
                                SphP[i].Hsml /= 1.26;
                        }
                    }

                    if (SphP[i].Hsml < All.MinGasHsml)
                        SphP[i].Hsml = All.MinGasHsml;
                }
                else
                    P[i].Ti_endstep = -P[i].Ti_endstep - 1; /* Mark as inactive */
#endif
            }
        }

        numlist = malloc(NTask * sizeof(int) * NTask);
        MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
        for (i = 0, ntot = 0; i < NTask; i++)
            ntot += numlist[i];
        free(numlist);

        if (ntot > 0)
        {
            iter++;
            if (iter > MAXITER)
            {
                printf("failed to converge in neighbour iteration in density()\n");
                fflush(stdout);
                endrun(1155);
            }
        }
    } while (ntot > 0);
    /* mark as active again */
    for (i = 0; i < NumPart; i++)
        if (P[i].Ti_endstep < 0)
            P[i].Ti_endstep = -P[i].Ti_endstep - 1;

    free(ndonelist);
    free(nsend);
    free(nsend_local);
    free(nbuffer);
    free(noffset);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
    int j, n, startnode, numngb, numngb_inbox, type;
    double h, h2, fac, hinv, hinv3, hinv4;
    double rho, divv, wk, dwk;
    double dx, dy, dz, r, r2, u, mass_j;
    double dvx, dvy, dvz, rotv[3];
    double weighted_numngb;
    FLOAT *pos, *vel;

#ifndef HSMLCONSTANT
    double dhsmlrho;
    dhsmlrho = 0;
#endif
    if (mode == 0)
    {
        pos = P[target].Pos;
        vel = SphP[target].VelPred;
        h = SphP[target].Hsml;
        type = P[target].Type;

    }
    else
    {
        pos = DensDataGet[target].Pos;
        vel = DensDataGet[target].Vel;
        h = DensDataGet[target].Hsml;     
        type = DensDataGet[target].Type;
	
    }

    h2 = h * h;
    hinv = 1.0 / h;
#ifdef DIM1
    hinv3 = hinv;
#else
    hinv3 = hinv * hinv * hinv;
#endif
    hinv4 = hinv3 * hinv;
    rho = divv = rotv[0] = rotv[1] = rotv[2] = 0.;
    weighted_numngb = 0;
    startnode = All.MaxPart;

    numngb = 0;
    do
    {
        numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, type);

        for (n = 0; n < numngb_inbox; n++)
        {
            j = Ngblist[n];

            dx = pos[0] - P[j].Pos[0];
            dy = pos[1] - P[j].Pos[1];
            dz = pos[2] - P[j].Pos[2];
            r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < h2)
            {
                numngb++;

                r = sqrt(r2);

                u = r * hinv;

                if (u < 0.5)
                {
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1.0) * u * u);
                    dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                }
                else
                {
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                }

                mass_j = P[j].Mass;

                rho += mass_j * wk;

                weighted_numngb += NORM_COEFF * wk / hinv3;
#ifndef HSMLCONSTANT
                dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);
#endif
                if (r > 0)
                {
                    fac = mass_j * dwk / r;

                    dvx = vel[0] - SphP[j].VelPred[0];
                    dvy = vel[1] - SphP[j].VelPred[1];
                    dvz = vel[2] - SphP[j].VelPred[2];

                    divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

                    rotv[0] += fac * (dz * dvy - dy * dvz);
                    rotv[1] += fac * (dx * dvz - dz * dvx);
                    rotv[2] += fac * (dy * dvx - dx * dvy);
                    

                }
            }
        }
    } while (startnode >= 0);

    if (mode == 0)
    {
       
        SphP[target].Density = rho;
        SphP[target].DivVel = divv;
#ifndef HSMLCONSTANT
        SphP[target].DhsmlDensityFactor = dhsmlrho;
        SphP[target].NumNgb = weighted_numngb;
#else
        SphP[target].NumNgb = numngb;
#endif
        SphP[target].Rot[0] = rotv[0];
        SphP[target].Rot[1] = rotv[1];
        SphP[target].Rot[2] = rotv[2];
    }
    else
    {
        DensDataResult[target].Rho = rho;
        DensDataResult[target].Div = divv;

#ifndef HSMLCONSTANT
        DensDataResult[target].DhsmlDensity = dhsmlrho;
        DensDataResult[target].Ngb = weighted_numngb;
#else
        DensDataResult[target].Ngb = numngb;
#endif
        DensDataResult[target].Rot[0] = rotv[0];
        DensDataResult[target].Rot[1] = rotv[1];
        DensDataResult[target].Rot[2] = rotv[2];
    }
}

/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void* a, const void* b)
{
    if (((struct densdata_in*)a)->Task < (((struct densdata_in*)b)->Task))
        return -1;

    if (((struct densdata_in*)a)->Task > (((struct densdata_in*)b)->Task))
        return +1;

    return 0;
}
