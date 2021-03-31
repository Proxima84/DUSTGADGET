#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

/*! \file hydra.c
 *  \brief Computation of SPH forces and rate of entropy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */

/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void hydro_force(void)
{
    long long ntot, ntotleft;
    int i, j, k, n, ngrp, maxfill, source, ndone;
    int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
    int level, recvTask, nexport, place;
    double soundspeed_i;
    MPI_Status status;
    /* `NumSphUpdate' gives the number of particles on this processor that want a force update */
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


    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);


    i = 0; /* first particle for this task */
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

                    hydro_evaluate(i, 0);

                    for (j = 0; j < NTask; j++)
                    {
                        if (Exportflag[j])
                        {
                            for (k = 0; k < 3; k++)
                            {
                                HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
                                HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
                            }
                            HydroDataIn[nexport].Hsml = SphP[i].Hsml;
                            HydroDataIn[nexport].Mass = P[i].Mass;
                  
#ifndef HSMLCONSTANT
                            HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
#endif
                            HydroDataIn[nexport].Density = SphP[i].Density;
                            HydroDataIn[nexport].Pressure = SphP[i].Pressure;
                            HydroDataIn[nexport].Timestep = P[i].Ti_endstep - P[i].Ti_begstep;
                            /* calculation of F1 */
                            soundspeed_i = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
                            HydroDataIn[nexport].F1 = fabs(SphP[i].DivVel)
                                / (fabs(SphP[i].DivVel) + SphP[i].CurlVel
                                                          + 0.0001 * soundspeed_i / SphP[i].Hsml);

                            HydroDataIn[nexport].Index = i;
                            HydroDataIn[nexport].Task = j;
                            nexport++;
                            nsend_local[j]++;
                        }
                    }
                }
            }
        qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

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
                if (maxfill >= All.BunchSizeHydro)
                    break;

                recvTask = ThisTask ^ ngrp;

                if (recvTask < NTask)
                {
                    if (nsend[ThisTask * NTask + recvTask] > 0
                        || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
                            nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE, recvTask,
                            TAG_HYDRO_A, &HydroDataGet[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in),
                            MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
                    }
                }

                for (j = 0; j < NTask; j++)
                    if ((j ^ ngrp) < NTask)
                        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
            }
            /* now do the imported particles */
            for (j = 0; j < nbuffer[ThisTask]; j++)
                hydro_evaluate(j, 1);
            /* do a block to measure imbalance */
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
                if (maxfill >= All.BunchSizeHydro)
                    break;

                recvTask = ThisTask ^ ngrp;

                if (recvTask < NTask)
                {
                    if (nsend[ThisTask * NTask + recvTask] > 0
                        || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
                            MPI_BYTE, recvTask, TAG_HYDRO_B,
                            &HydroDataPartialResult[noffset[recvTask]],
                            nsend_local[recvTask] * sizeof(struct hydrodata_out), MPI_BYTE,
                            recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

                        /* add the result to the particles */
                        for (j = 0; j < nsend_local[recvTask]; j++)
                        {
                            source = j + noffset[recvTask];
                            place = HydroDataIn[source].Index;

                            for (k = 0; k < 3; k++)
                                SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];

                            SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;

                            if (SphP[place].MaxSignalVel
                                < HydroDataPartialResult[source].MaxSignalVel)
                                SphP[place].MaxSignalVel
                                    = HydroDataPartialResult[source].MaxSignalVel;
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

    /* do final operations on results */


    for (i = 0; i < N_gas; i++)
        if (P[i].Ti_endstep == All.Ti_Current)
        {
            SphP[i].DtEntropy *= GAMMA_MINUS1 / pow(SphP[i].Density, GAMMA_MINUS1);
#ifdef SPH_BND_PARTICLES
            if (P[i].ID == 0)
            {
                SphP[i].DtEntropy = 0;
                for (k = 0; k < 3; k++)
                    SphP[i].HydroAccel[k] = 0;
            }
#endif
        }
}

/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void hydro_evaluate(int target, int mode)
{
    int j, k, n, startnode, numngb;//,id;
    FLOAT *pos, *vel;
    FLOAT h_i, rho, pressure, f1, f2;
    double acc[3], dtEntropy;
    double dx, dy, dz, dvx, dvy, dvz;
    double h_i2, hinv, hinv4;
    double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
    double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
    double h_j, dwk_j, r, r2, u, hfc_visc;
    double dt, mass, maxSignalVel, timestep;
    
#ifndef HSMLCONSTANT
    FLOAT dhsmlDensityFactor;
#endif
    if (mode == 0)
    {
        pos = P[target].Pos;
        vel = SphP[target].VelPred;
        h_i = SphP[target].Hsml;
        mass = P[target].Mass;
#ifndef HSMLCONSTANT
        dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
#endif
        rho = SphP[target].Density;
        pressure = SphP[target].Pressure;
        timestep = P[target].Ti_endstep - P[target].Ti_begstep;
        soundspeed_i = sqrt(GAMMA * pressure / rho);
        f1 = fabs(SphP[target].DivVel) / (fabs(SphP[target].DivVel) + SphP[target].CurlVel
                                             + 0.0001 * soundspeed_i / SphP[target].Hsml);      
    }
    else
    {
        pos = HydroDataGet[target].Pos;
        vel = HydroDataGet[target].Vel;
        h_i = HydroDataGet[target].Hsml;
        mass = HydroDataGet[target].Mass;
#ifndef HSMLCONSTANT
        dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
#endif
        rho = HydroDataGet[target].Density;
        pressure = HydroDataGet[target].Pressure;
        timestep = HydroDataGet[target].Timestep;
        soundspeed_i = sqrt(GAMMA * pressure / rho);
        f1 = HydroDataGet[target].F1;
    }
    /* initialize variables before SPH loop is started */
    acc[0] = acc[1] = acc[2] = dtEntropy = 0.;
    maxSignalVel = 0.;
#ifdef HSMLCONSTANT
    p_over_rho2_i = pressure / (rho * rho);
#else
    p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
#endif
    h_i2 = h_i * h_i;

    /* Now start the actual SPH computation for this particle */
    startnode = All.MaxPart;
    do
    {
        numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode, 0);

        for (n = 0; n < numngb; n++)
        {
            j = Ngblist[n];

            dx = pos[0] - P[j].Pos[0];
            dy = pos[1] - P[j].Pos[1];
            dz = pos[2] - P[j].Pos[2];
            r2 = dx * dx + dy * dy + dz * dz;
            h_j = SphP[j].Hsml;
            if (r2 < h_i2 || r2 < h_j * h_j)
            {
                r = sqrt(r2);
                if (r > 0)
                {
                    p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
                    soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].Density);
                    dvx = vel[0] - SphP[j].VelPred[0];
                    dvy = vel[1] - SphP[j].VelPred[1];
                    dvz = vel[2] - SphP[j].VelPred[2];
                    vdotr = dx * dvx + dy * dvy + dz * dvz;
                    vdotr2 = vdotr;
                    if (r2 < h_i2)
                    {
                        hinv = 1.0 / h_i;
#ifdef DIM1
                        hinv4 = hinv * hinv;
#else
                        hinv4 = hinv * hinv * hinv * hinv;
#endif

                        u = r * hinv;
                        if (u < 0.5)
                            dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                        else
                            dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                    }
                    else
                    {
                        dwk_i = 0;
                    }

                    if (r2 < h_j * h_j)
                    {
                        hinv = 1.0 / h_j;
#ifdef DIM1
                        hinv4 = hinv * hinv;
#else
                        hinv4 = hinv * hinv * hinv * hinv;
#endif
                        u = r * hinv;
                        if (u < 0.5)
                            dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                        else
                            dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                    }
                    else
                    {
                        dwk_j = 0;
                    }

                    if (soundspeed_i + soundspeed_j > maxSignalVel)
                        maxSignalVel = soundspeed_i + soundspeed_j;

                    if (vdotr2 < 0) /* ... artificial viscosity */
                    {
                        mu_ij = 0.5 * (h_i + h_j) * vdotr2 / r2; /* note: this is negative! */

                        vsig = 0.5 * (soundspeed_i + soundspeed_j); // - 3 * mu_ij;

                        rho_ij = 0.5 * (rho + SphP[j].Density);
                        f2 = fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel
                                                        + 0.0001 * soundspeed_j / SphP[j].Hsml);

                        visc = -0.5 * All.ArtBulkViscConst * vsig * mu_ij / rho_ij * (f1 + f2);

                        /* .... end artificial viscosity evaluation */
                        /* make sure that viscous acceleration is not too large */
                        dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep))
                            * All.Timebase_interval;
                        if (dt > 0 && (dwk_i + dwk_j) < 0)
                        {
                            visc = dmin(visc, 0.5 * vdotr2
                                    / (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
                        }
                    }
                    else
                        visc = 0;
#ifndef HSMLCONSTANT
                    p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
#endif
                    hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;

                    hfc = hfc_visc
                        + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
                   
                    acc[0] -= hfc * dx;
                    acc[1] -= hfc * dy;
                    acc[2] -= hfc * dz;
 
#if !defined(SAVETEMP_LTR) && !defined(ISOTHERM_EQS)
                    dtEntropy += 0.5 * hfc_visc * vdotr2;
#endif
                }
            }
        }
    } while (startnode >= 0);

    /* Now collect the result at the right place */
    if (mode == 0)
    {
        for (k = 0; k < 3; k++)
            SphP[target].HydroAccel[k] = acc[k];
        SphP[target].DtEntropy = dtEntropy;
        SphP[target].MaxSignalVel = maxSignalVel;
    }
    else
    {
        for (k = 0; k < 3; k++)
            HydroDataResult[target].Acc[k] = acc[k];
        HydroDataResult[target].DtEntropy = dtEntropy;
        HydroDataResult[target].MaxSignalVel = maxSignalVel;
    }
}

/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int hydro_compare_key(const void* a, const void* b)
{
    if (((struct hydrodata_in*)a)->Task < (((struct hydrodata_in*)b)->Task))
        return -1;
    if (((struct hydrodata_in*)a)->Task > (((struct hydrodata_in*)b)->Task))
        return +1;
    return 0;
}
