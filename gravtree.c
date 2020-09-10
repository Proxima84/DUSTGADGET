#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gravity_tree(void)
{
    long long ntot;
    int i, j;
    int* numlist;

#ifndef NOGRAVITY
    int nexport, *ndonelist;
    int nexportsum = 0;
    int iter = 0;
    int *noffset, *nbuffer, *nsend, *nsend_local;
    long long ntotleft;
    int ndone, maxfill, ngrp;
    int k, place;
    int level, recvTask;
    double ax, ay, az;
    MPI_Status status;
#endif
    /* contruct tree if needed */
    if (TreeReconstructFlag)
    {
        force_treebuild(NumPart);
        TreeReconstructFlag = 0;
    }


    /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
    for (i = 0, ntot = 0; i < NTask; i++)
        ntot += numlist[i];
    free(numlist);

#ifndef NOGRAVITY
    /*if(ThisTask == 0)
      printf("Begin tree force.\n");*/
    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);

    i = 0; /* beginn with this index */
    ntotleft = ntot; /* particles left for all tasks together */
    while (ntotleft > 0)
    {
        iter++;
#ifdef NOSTARSELFGRAV
        star_gravity();
#endif
        for (j = 0; j < NTask; j++)
            nsend_local[j] = 0;
        /* do local particles and prepare export list */
        for (nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
            if (P[i].Ti_endstep == All.Ti_Current)
            {
                ndone++;
#ifdef NOSTARSELFGRAV
                if (P[i].ItsStar)
                    continue;
#endif
                for (j = 0; j < NTask; j++)
                    Exportflag[j] = 0;

                force_treeevaluate(i, 0);
                for (j = 0; j < NTask; j++)
                {
                    if (Exportflag[j])
                    {
                        for (k = 0; k < 3; k++)
                            GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
                        GravDataGet[nexport].Type = P[i].Type;
                        if (P[i].Type < 2)
                            GravDataGet[nexport].Soft = SphP[i].Hsml;
                        GravDataGet[nexport].w.OldAcc = P[i].OldAcc;
                        GravDataIndexTable[nexport].Task = j;
                        GravDataIndexTable[nexport].Index = i;
                        GravDataIndexTable[nexport].SortIndex = nexport;
                        nexport++;
                        nexportsum++;
                        nsend_local[j]++;
                    }
                }
            }
        qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

        for (j = 0; j < nexport; j++)
            GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

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
                if (maxfill >= All.BunchSizeForce)
                    break;
                recvTask = ThisTask ^ ngrp;

                if (recvTask < NTask)
                {
                    if (nsend[ThisTask * NTask + recvTask] > 0
                        || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
                            nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE, recvTask,
                            TAG_GRAV_A, &GravDataGet[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
                            MPI_BYTE, recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
                    }
                }

                for (j = 0; j < NTask; j++)
                    if ((j ^ ngrp) < NTask)
                        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
            }

            for (j = 0; j < nbuffer[ThisTask]; j++)
            {
                force_treeevaluate(j, 1);
            }
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
                if (maxfill >= All.BunchSizeForce)
                    break;

                recvTask = ThisTask ^ ngrp;
                if (recvTask < NTask)
                {
                    if (nsend[ThisTask * NTask + recvTask] > 0
                        || nsend[recvTask * NTask + ThisTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
                            MPI_BYTE, recvTask, TAG_GRAV_B, &GravDataOut[noffset[recvTask]],
                            nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE, recvTask,
                            TAG_GRAV_B, MPI_COMM_WORLD, &status);

                        /* add the result to the particles */
                        for (j = 0; j < nsend_local[recvTask]; j++)
                        {
                            place = GravDataIndexTable[noffset[recvTask] + j].Index;

                            for (k = 0; k < 3; k++)
                                P[place].GravAccel[k]
                                    += GravDataOut[j + noffset[recvTask]].u.Acc[k];

                            P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
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

    for (i = 0; i < NumPart; i++)
        if (P[i].Ti_endstep == All.Ti_Current)
        {
            ax = P[i].GravAccel[0];
            ay = P[i].GravAccel[1];
            az = P[i].GravAccel[2];
            P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
        }
    if (All.TypeOfOpeningCriterion == 1)
        All.ErrTolTheta
            = 0; /* This will switch to the relative opening criterion for the following force
                    computations */

    /*  muliply by G */
    for (i = 0; i < NumPart; i++)
        if (P[i].Ti_endstep == All.Ti_Current)
            for (j = 0; j < 3; j++)
                P[i].GravAccel[j] *= All.G;
/*  if(ThisTask == 0)
    printf("tree is done.\n");*/

#else /* gravity is switched off */

    for (i = 0; i < NumPart; i++)
        if (P[i].Ti_endstep == All.Ti_Current)
            for (j = 0; j < 3; j++)
                P[i].GravAccel[j] = 0;

#endif
}
#ifdef NOSTARSELFGRAV
void star_gravity()
{
    MPI_Status status;
    int StarsNum, n, i, stot, j, k;
    int *noffset, *numlist;
    double fac;
    for (n = N_gas, StarsNum = 0; n < NumPart; n++)
    {
        for (i = 0; i < 3; i++)
            GravDataIn[StarsNum].u.Pos[i] = P[n].Pos[i];
        GravDataIn[StarsNum].id = P[n].ID;
        GravDataIn[StarsNum].Mass = P[n].Mass;
        StarsNum++;
    }
    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&StarsNum, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
    for (i = 0, stot = 0; i < NTask; i++)
        stot += numlist[i];

    if (StarsNum > 0)
    {
        for (i = 0; i < NTask; i++)
            if (i != ThisTask)
                MPI_Send(&GravDataIn[0], StarsNum * sizeof(struct gravdata_in), MPI_BYTE, i,
                    TAG_GRAV_A, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (j = 1, noffset[0] = 0; j < NTask; j++)
    {
        if ((j - 1) == ThisTask)
            noffset[j] = noffset[j - 1];
        else
            noffset[j] = noffset[j - 1] + numlist[j - 1];
    }
    for (i = 0; i < NTask; i++)
    {
        if ((i != ThisTask) && (numlist[i] > 0))
            MPI_Recv(&GravDataGet[noffset[i]], numlist[i] * sizeof(struct gravdata_in), MPI_BYTE, i,
                TAG_GRAV_A, MPI_COMM_WORLD, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(noffset);
    free(numlist);
    k = stot - StarsNum;
    for (j = 0, n = StarsNum; j < k; j++)
    {

        for (i = 0; i < 3; i++)
            GravDataIn[n].u.Pos[i] = GravDataGet[j].u.Pos[i];
        GravDataIn[n].id = GravDataGet[j].id;
        GravDataIn[n].Mass = GravDataGet[j].Mass;
        n++;
    }
    double r2, dx, dy, dz, mass, r;
    double acc_x, acc_y, acc_z;
    for (i = N_gas; i < NumPart; i++)
    {
        if (P[i].Ti_endstep == All.Ti_Current)
        {
            if (P[i].ItsStar)
            {
                acc_x = 0;
                acc_y = 0;
                acc_z = 0;
                for (j = 0; j < stot; j++)
                {
                    if (P[i].ID == GravDataIn[j].id)
                        continue;
                    dx = GravDataIn[j].u.Pos[0] - P[i].Pos[0];
                    dy = GravDataIn[j].u.Pos[1] - P[i].Pos[1];
                    dz = GravDataIn[j].u.Pos[2] - P[i].Pos[2];
                    mass = GravDataIn[j].Mass;
                    r2 = dx * dx + dy * dy + dz * dz;
                    r = sqrt(r2);
                    fac = mass / (r2 * r);
                    acc_x += dx * fac;
                    acc_y += dy * fac;
                    acc_z += dz * fac;
                }
                P[i].GravAccel[0] = acc_x;
                P[i].GravAccel[1] = acc_y;
                P[i].GravAccel[2] = acc_z;
            }
        }
    }
}
#endif

/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
    int i;
    for (i = 0; i < 2; i++)
        All.ForceSoftening[i] = All.HsmlConstant;   
}

/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int grav_tree_compare_key(const void* a, const void* b)
{
    if (((struct gravdata_index*)a)->Task < (((struct gravdata_index*)b)->Task))
        return -1;

    if (((struct gravdata_index*)a)->Task > (((struct gravdata_index*)b)->Task))
        return +1;

    return 0;
}
