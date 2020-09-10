#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

/*! file accretion.c
 * calculate accretion on the stars
 */

/*! This function finds the shp particle inside accretion radius of stars, and excludes
 them from calculation.
 */
#ifndef NOGRAVITY
void accretion(void)
{
    long long ntot, ntotleft;
    int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
    int i, j, n, ndone, maxfill, source;
    int level, ngrp, recvTask, place, nexport;
    MPI_Status status;

    noffset = malloc(sizeof(int) * NTask); /* offsets of bunches in common list */
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);
    int StarUpdate;
    for (n = N_gas, StarUpdate = 0; n < NumPart; n++)
    {
        StarUpdate++;
    }

    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&StarUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
    for (i = 0, ntot = 0; i < NTask; i++)
        ntot += numlist[i];
    free(numlist);
    i = N_gas; /* beginn with this index */
    ntotleft = ntot; /* particles left for all tasks together */

    while (ntotleft > 0)
    {
        for (j = 0; j < NTask; j++)
            nsend_local[j] = 0;

        /* do local particles and prepare export list */

        for (nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeDensity - NTask; i++)
        {
            ndone++;

            for (j = 0; j < NTask; j++)
                Exportflag[j] = 0;
            accretion_evaluate(i, 0);

            for (j = 0; j < NTask; j++)
            {
                if (Exportflag[j])
                {
                    DensDataIn[nexport].Pos[0] = P[i].Pos[0];
                    DensDataIn[nexport].Pos[1] = P[i].Pos[1];
                    DensDataIn[nexport].Pos[2] = P[i].Pos[2];
                    DensDataIn[nexport].Hsml = All.AccretionRadius * P[i].Radius;
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
                            nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE, recvTask,
                            TAG_DENS_A, &DensDataGet[nbuffer[ThisTask]],
                            nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
                            MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
                    }
                }

                for (j = 0; j < NTask; j++)
                    if ((j ^ ngrp) < NTask)
                        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
            }

            for (j = 0; j < nbuffer[ThisTask]; j++)
                accretion_evaluate(j, 1);
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
                            nsend_local[recvTask] * sizeof(struct densdata_out), MPI_BYTE, recvTask,
                            TAG_DENS_B, MPI_COMM_WORLD, &status);

                        /* add the result to the particles */
                        for (j = 0; j < nsend_local[recvTask]; j++)
                        {
                            source = j + noffset[recvTask];
                            place = DensDataIn[source].Index;

                            P[place].Accretion += DensDataPartialResult[source].Ngb;
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

void remove_particles(int curtime)
{
    int NumAcc, i, j, ntot, Acc_glob;
    double R2, z, dist;
    NumAcc = 0;
    dist = All.BarrierDistance * All.BarrierDistance;
    for (i = 0; i < N_gas; i++)
    {
        R2 = P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1];
        z = P[i].Pos[2] * P[i].Pos[2];
        if ((P[i].Accretion) || (R2 > dist) || (z > (dist / 25.)))
        {
            NumAcc++;
            if (P[i].Ti_endstep == curtime)
            {
                NumForceUpdate--;
            }
            P[i].ID = -1;
        }
    }

    MPI_Allreduce(&NumAcc, &Acc_glob, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    DomainDecompositionFlag = Acc_glob;
    if (Acc_glob)
    {
        for (j = 0; j < NumAcc; j++)
        {
            for (i = 0; i < NumPart - 1; i++)
            {
                if (P[i].ID == -1)
                {
                    P[i] = P[i + 1];
                    if (i < (N_gas - 1))
                        SphP[i] = SphP[i + 1];
                    P[i + 1].ID = -1;
                    P[i + 1].Mass = 0;
                }
            }
            N_gas--;
            NumPart--;
        }
        int* list_NumPart;
        list_NumPart = malloc(sizeof(int) * NTask);
        MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MPI_COMM_WORLD);
        for (i = 0, ntot = 0; i < NTask; i++)
        {
            ntot += list_NumPart[i];
        }
        free(list_NumPart);
        All.TotNumPart = ntot;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void accretion_evaluate(int target, int mode)
{
    int j, n, startnode, numngb, numngb_inbox;
    double h, h2;
    double dx, dy, dz, r2;
    FLOAT* pos;
    if (mode == 0)
    {
        pos = P[target].Pos;
        h = All.AccretionRadius * P[target].Radius;
    }
    else
    {
        pos = DensDataGet[target].Pos;
        h = DensDataGet[target].Hsml;
    }

    h2 = h * h;

    startnode = All.MaxPart;
    numngb = 0;
    do
    {
        numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, 0);

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
                P[j].Accretion = 1;
            }
        }
    } while (startnode >= 0);
    startnode = All.MaxPart;
    do
    {
        numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, 1);

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
                P[j].Accretion = 1;
            }
        }
    } while (startnode >= 0);

    if (mode == 0)
    {
        P[target].Accretion += numngb;
    }
    else
    {
        DensDataResult[target].Ngb = numngb;
    }
}
#endif

