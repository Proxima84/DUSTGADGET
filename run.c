#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over single
 *  timesteps. The loop terminates when the cpu-time limit is reached, when a
 *  `stop' file is found in the output directory, or when the simulation ends
 *  because we arrived at TimeMax.
 */
void run(void)
{
    FILE* fd;
    int stopflag = 0;
    char stopfname[200], contfname[200];
    sprintf(stopfname, "%sstop", All.OutputDir);
    sprintf(contfname, "%scont", All.OutputDir);
    unlink(contfname);

    do /* main loop */
    {
        find_next_sync_point_and_drift(); /* find next synchronization point and drift particles to
                                           * this time.
                                           * If needed, this function will also write an output file
                                           * at the desired time.
                                           */
        if (ThisTask == 0)
        {
            printf("\nBegin Step %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
                All.TimeStep);
        }
        domain_Decomposition(); /* do domain decomposition if needed */
        compute_accelerations(0); /* compute accelerations for
                                   * the particles that are to be advanced
                                   */
        advance_and_find_timesteps(); /* 'kick' active particles in
                                       * momentum space and compute new
                                       * timesteps for them
                                       */
        All.NumCurrentTiStep++;

        /* Check whether we need to interrupt the run */
        if (ThisTask == 0)
        {
            /* Is the stop-file present? If yes, interrupt the run. */
            if ((fd = fopen(stopfname, "r")))
            {
                fclose(fd);
                stopflag = 1;
                unlink(stopfname);
            }
        }

        MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (stopflag)
        {
            restart(0); /* write restart file */
            MPI_Barrier(MPI_COMM_WORLD);

            if (stopflag == 2 && ThisTask == 0)
            {
                if ((fd = fopen(contfname, "w")))
                    fclose(fd);
            }
            return;
        }

    } while (All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

    restart(0);
    savepositions(All.SnapshotFileCount++);
    /* write a last snapshot
                                                     * file at final time (will
                                                     * be overwritten if
                                                     * All.TimeMax is increased
                                                     * and the run is continued)
                                                     */
}


/*! This function finds the next synchronization point of the system (i.e. the
 *  earliest point of time any of the particles needs a force computation),
 *  and drifts the system to this point of time.  If the system drifts over
 *  the desired time of a snapshot file, the function will drift to this
 *  moment, generate an output, and then resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
    int n, min, min_glob, flag, *temp;
    double timeold;
    timeold = All.Time;

    for (n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
        if (min > P[n].Ti_endstep)
            min = P[n].Ti_endstep;

    MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    /* We check whether this is a full step where all particles are synchronized */
    flag = 1;
    for (n = 0; n < NumPart; n++)
        if (P[n].Ti_endstep > min_glob)
            flag = 0;

    MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to
     * be active */
    for (n = 0, NumForceUpdate = 0; n < NumPart; n++)
    {
        if (P[n].Ti_endstep == min_glob)
            NumForceUpdate++;
    }
    /* note: NumForcesSinceLastDomainDecomp has type "long long" */
    temp = malloc(NTask * sizeof(int));
    MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
    for (n = 0; n < NTask; n++)
        All.NumForcesSinceLastDomainDecomp += temp[n];
    free(temp);
    while (min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
        move_particles(All.Ti_Current, All.Ti_nextoutput);
        All.Ti_Current = All.Ti_nextoutput;
#ifndef NOGRAVITY
        accretion();
        saveaccretion();
        remove_particles(min_glob);
#endif

        All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;
        savepositions(All.SnapshotFileCount++); /* write snapshot file */
        All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }
    move_particles(All.Ti_Current, min_glob);
#ifndef NOGRAVITY
    accretion();
    remove_particles(min_glob);
#endif

    All.Ti_Current = min_glob;
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;
    All.TimeStep = All.Time - timeold;
}


/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
int find_next_outputtime(int ti_curr)
{
    int ti, ti_next, iter = 0;
    double time;

    ti_next = -1;
    time = All.TimeOfFirstSnapshot;
    iter = 0;
    while (time < All.TimeBegin)
    {
        time += All.TimeBetSnapshot;
        iter++;
        if (iter > 1000000)
        {
            printf("Can't determine next output time.\n");
            endrun(110);
        }
    }
    while (time <= All.TimeMax)
    {
        ti = (time - All.TimeBegin) / All.Timebase_interval;
        if (ti >= ti_curr)
        {
            ti_next = ti;
            break;
        }
        time += All.TimeBetSnapshot;
        iter++;
        if (iter > 1000000)
        {
            printf("Can't determine next output time.\n");
            endrun(111);
        }
    }
    if (ti_next == -1)
    {
        ti_next = 2 * TIMEBASE; /* this will prevent any further output */

        if (ThisTask == 0)
            printf("\nThere is no valid time for a further snapshot file.\n");
    }

    return ti_next;
}

