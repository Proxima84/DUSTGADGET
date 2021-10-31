#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file read_ic.c
 *  \brief Read initial conditions in one of Gadget's file formats
 */

/*! This function reads initial conditions, in one of the three possible file
 *  formats currently supported by Gadget.  Note: When a snapshot file is
 *  started from initial conditions (start-option 0), not all the information
 *  in the header is used, in particular, the STARTING TIME needs to be set in
 *  the parameterfile.  Also, for gas particles, only the internal energy is
 *  read, the density and mean molecular weight will be recomputed by the
 *  code.  When InitGasTemp>0 is given, the gas temperature will be initialzed
 *  to this value assuming a mean colecular weight either corresponding to
 *  complete neutrality, or full ionization.
 *
 */
void read_ic(char* fname)
{
    int i, masterTask, lastTask;
    double u_init;
    char buf[500];
    NumPart = 0;
    N_gas = 0;
    All.TotNumPart = 0;
    masterTask = 0;
    lastTask = NTask - 1;
    sprintf(buf, "%s", fname);
    read_file(buf, masterTask, lastTask);
    MPI_Barrier(MPI_COMM_WORLD);
    if(N_gas==0)
    {
      printf("\nno gas or dust matter in initial condition file\n");
      fflush(stdout);
      endrun(115);
    }
    /* this makes sure that masses are initialized in the case that the mass-block
       is completely empty */
    for (i = 0; i < NumPart; i++)
    {
        if (All.MassTable[P[i].Type] == 0)
            All.MassTable[P[i].Type] = P[i].Mass;
    }

    if (RestartFlag == 0)
    {
        if (All.InitGasTemp > 0)
        {
            u_init = (BOLTZMANN / PROTONMASS) * All.InitGasTemp;
            u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs; /* unit conversion */

#if defined(ISOTHERM_EQS) || defined(SAVETEMP_LTR)
            u_init *= 1.0;
#else
            double molecular_weight;
            u_init *= (1.0 / GAMMA_MINUS1);

            if (All.InitGasTemp > 1.0e4) /* assuming FULL ionization */
                molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
            else /* assuming NEUTRAL GAS */
                molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

            u_init /= molecular_weight;
#endif

            for (i = 0; i < N_gas; i++)
            {
                if (SphP[i].Entropy == 0)
                    SphP[i].Entropy = u_init;

                /* Note: the coversion to entropy will be done in the function init(),
                   after the densities have been computed */
            }
        }
    }

   /* for (i = 0; i < N_gas; i++)
        SphP[i].Entropy = dmax(All.MinEgySpec, SphP[i].Entropy);*/

    MPI_Barrier(MPI_COMM_WORLD);
    if (ThisTask == 0)
    {
        printf("Total number of particles :  %d%09d\n\n", (int)(All.TotNumPart / 1000000000),
            (int)(All.TotNumPart % 1000000000));
        fflush(stdout);
    }
}


/*! This function reads out the buffer that was filled with particle data, and
 *  stores it at the appropriate place in the particle structures.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
    int n, k;
    float* fp;
    int* ip;
    fp = CommBuffer;
    ip = CommBuffer;

    switch (blocknr)
    {
        case IO_POS: /* positions */
            for (n = 0; n < pc; n++)
                for (k = 0; k < 3; k++)
                    P[offset + n].Pos[k] = *fp++;

            for (n = 0; n < pc; n++)
                P[offset + n].Type = type; /* initialize type here as well */
            break;

        case IO_VEL: /* velocities */
            for (n = 0; n < pc; n++)
                for (k = 0; k < 3; k++)
                    P[offset + n].Vel[k] = *fp++;
            break;

        case IO_ID: /* particle ID */
            for (n = 0; n < pc; n++)
                P[offset + n].ID = *ip++;
            break;

        case IO_MASS: /* particle mass */
            for (n = 0; n < pc; n++)
                P[offset + n].Mass = *fp++;
            break;

        case IO_U: /* temperature */
            for (n = 0; n < pc; n++)
                SphP[offset + n].Entropy = *fp++;
            break;

        case IO_RHO: /* density */
            for (n = 0; n < pc; n++)
                SphP[offset + n].Density = *fp++;
            break;


        case IO_HSML: /* SPH smoothing length */
            for (n = 0; n < pc; n++)
                SphP[offset + n].Hsml = *fp++;
            break;

        /* the other input fields (if present) are not needed to define the
           initial conditions of the code */

        case IO_POT:
        case IO_ACCEL:
        case IO_DTENTR:
        case IO_TSTP:
            break;
    }
}

/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
void read_file(char* fname, int readTask, int lastTask)
{
    int blockmaxlen;
    int i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task;
    int blksize1, blksize2;
    MPI_Status status;
    FILE* fd = 0;
    int nall;
    int type;
    int nstart, bytes_per_blockelement, npart, typelist[3];
    enum iofields blocknr;
#define SKIP                                                                                       \
    {                                                                                              \
        my_fread(&blksize1, sizeof(int), 1, fd);                                                   \
    }
#define SKIP2                                                                                      \
    {                                                                                              \
        my_fread(&blksize2, sizeof(int), 1, fd);                                                   \
    }

    if (ThisTask == readTask)
    {
        if (!(fd = fopen(fname, "r")))
        {
            printf("can't open file `%s' for reading initial conditions.\n", fname);
            endrun(123);
        }
        SKIP;
        my_fread(&header, sizeof(header), 1, fd);
        SKIP2;
        if (blksize1 != 256 || blksize2 != 256)
        {
            printf("incorrect header format\n");
            fflush(stdout);
            endrun(890);
        }
        for (task = readTask + 1; task <= lastTask; task++)
            MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
    }
    else
        MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);


    if (All.TotNumPart == 0)
    {
        for (i = 0; i < 3; i++)
            header.npartTotal[i] = header.npart[i];
        All.TotN_gas = 0;
        for (i = 0; i < 2; i++)
            All.TotN_gas
                += header.npartTotal[i]; // + (((long long) header.npartTotalHighWord[0]) << 32);

        for (i = 0, All.TotNumPart = 0; i < 3; i++)
        {
            All.TotNumPart += header.npartTotal[i];
            // All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
        }

        All.MaxPart = All.PartAllocFactor
            * (All.TotNumPart / NTask); /* sets the maximum number of particles that may */
        All.MaxPartSph = All.PartAllocFactor
            * (All.TotN_gas / NTask); /* sets the maximum number of particles that may
                                         reside on a processor */
        allocate_memory();

    }

    if (ThisTask == readTask)
    {
        for (i = 0, n_in_file = 0; i < 3; i++)
            n_in_file += header.npart[i];

        printf("\nreading file `%s' on task=%d (contains %d particles.)\n"
               "distributing this file to tasks %d-%d\n"
               "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 1 (dust):  %8d  (tot=%6d%09d) masstab=%g size=%g\n"
               "Type 2 (stars):  %8d  (tot=%6d%09d)\n",
            fname, ThisTask, n_in_file, readTask, lastTask, header.npart[0],
            (int)(header.npartTotal[0] / 1000000000), (int)(header.npartTotal[0] % 1000000000),
            All.MassTable[0], header.npart[1], (int)(header.npartTotal[1] / 1000000000),
            (int)(header.npartTotal[1] % 1000000000), All.MassTable[1], All.DustSize,
            header.npart[2], (int)(header.npartTotal[2] / 1000000000),
            (int)(header.npartTotal[2] % 1000000000));
        fflush(stdout);
    }

    ntask = lastTask - readTask + 1;


    /* to collect the gas particles all at the beginning (in case several
       snapshot files are read on the current CPU) we move the collisionless
       particles such that a gap of the right size is created */
    for (type = 0, nall = 0; type < 3; type++)
    {
        n_in_file = header.npart[type];

        n_for_this_task = n_in_file / ntask;
        if ((ThisTask - readTask) < (n_in_file % ntask))
            n_for_this_task++;

        nall += n_for_this_task;
    }
    memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));
    nstart = N_gas;

    for (blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
        if (blockpresent(blocknr))
        {
            if (RestartFlag == 0 && blocknr > IO_U)
                continue; /* ignore all other blocks in initial conditions */

            bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

            blockmaxlen = ((int)(All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

            npart = get_particles_in_block(blocknr, &typelist[0]);

            if (npart > 0)
            {
                if (ThisTask == readTask)
                {
                    SKIP;
                }

                for (type = 0, offset = 0; type < 3; type++)
                {
                    n_in_file = header.npart[type];
                    if (typelist[type] == 0)
                    {
                        n_for_this_task = n_in_file / ntask;
                        if ((ThisTask - readTask) < (n_in_file % ntask))
                            n_for_this_task++;

                        offset += n_for_this_task;
                    }
                    else
                    {
                        for (task = readTask; task <= lastTask; task++)
                        {
                            n_for_this_task = n_in_file / ntask;
                            if ((task - readTask) < (n_in_file % ntask))
                                n_for_this_task++;

                            if (task == ThisTask)
                                if (NumPart + n_for_this_task > All.MaxPart)
                                {
                                    printf("too many particles\n");
                                    endrun(1313);
                                }
                            do
                            {
                                pc = n_for_this_task;

                                if (pc > blockmaxlen)
                                    pc = blockmaxlen;

                                if (ThisTask == readTask)
                                {
                                    my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
                                }

                                if (ThisTask == readTask && task != readTask)
                                    MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE,
                                        task, TAG_PDATA, MPI_COMM_WORLD);

                                if (ThisTask != readTask && task == ThisTask)
                                    MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE,
                                        readTask, TAG_PDATA, MPI_COMM_WORLD, &status);

                                if (ThisTask == task)
                                {
                                    empty_read_buffer(blocknr, nstart + offset, pc, type);

                                    offset += pc;
                                }

                                n_for_this_task -= pc;
                            } while (n_for_this_task > 0);
                        }
                    }
                }
                if (ThisTask == readTask)
                {
                    SKIP2;
                    if (blksize1 != blksize2)
                    {
                        printf("incorrect block-sizes detected!\n");
                        printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask,
                            blocknr, blksize1, blksize2);
                        fflush(stdout);
                        endrun(1889);
                    }
                }
            }
        }
    }
    for (type = 0; type < 3; type++)
    {
        n_in_file = header.npart[type];

        n_for_this_task = n_in_file / ntask;
        if ((ThisTask - readTask) < (n_in_file % ntask))
            n_for_this_task++;

        NumPart += n_for_this_task;

        if (type < 2)
            N_gas += n_for_this_task;
    }

    if (ThisTask == readTask)
    {
        fclose(fd);
    }
}

/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */
void distribute_file(
    int nfiles, int firstfile, int firsttask, int lasttask, int* filenr, int* master, int* last)
{
    int ntask, filesleft, filesright, tasksleft;

    if (nfiles > 1)
    {
        ntask = lasttask - firsttask + 1;

        filesleft = (((double)(ntask / 2)) / ntask) * nfiles;
        if (filesleft <= 0)
            filesleft = 1;
        if (filesleft >= nfiles)
            filesleft = nfiles - 1;

        filesright = nfiles - filesleft;

        tasksleft = ntask / 2;

        distribute_file(
            filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
        distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr,
            master, last);
    }
    else
    {
        if (ThisTask >= firsttask && ThisTask <= lasttask)
        {
            *filenr = firstfile;
            *master = firsttask;
            *last = lasttask;
        }
    }
}
