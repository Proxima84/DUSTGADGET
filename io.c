#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"


/*! \file io.c
 *  \brief Routines for producing a snapshot file on disk.
 */

static int n_type[3];
static long long ntot_type_all[3];

/*! This function writes a snapshot of the particle distribution to one or
 *  several files using the selected file format.  If NumFilesPerSnapshot>1,
 *  the snapshot is distributed onto several files, several of them can be
 *  written simultaneously (up to NumFilesWrittenInParallel). Each file
 *  contains data from a group of processors.
 */
void savepositions(int num)
{
    char buf[500];
    int i, j, *temp, n, masterTask, lastTask;   
    /* determine global and local particle numbers */
    for (n = 0; n < 3; n++)
        n_type[n] = 0;

    for (n = 0; n < NumPart; n++)
        n_type[P[n].Type]++;

    /* because ntot_type_all[] is of type `long long', we cannot do a simple
     * MPI_Allreduce() to sum the total particle numbers
     */
    temp = malloc(NTask * 3 * sizeof(int));
    MPI_Allgather(n_type, 3, MPI_INT, temp, 3, MPI_INT, MPI_COMM_WORLD);
    for (i = 0; i < 3; i++)
    {
        ntot_type_all[i] = 0;
        for (j = 0; j < NTask; j++)
            ntot_type_all[i] += temp[j * 3 + i];
    }
    free(temp);

    /* assign processors to output files */
    //  distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);
    masterTask = 0;
    lastTask = NTask - 1;
    sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);
    write_file(buf, masterTask, lastTask);
    MPI_Barrier(MPI_COMM_WORLD);
    if(ThisTask == 0)
      printf("\nwrited snapshot file %i \n",num);
}

void saveaccretion()
{
    int StarsNum, i;
    int *accStars, *IdStars;
    int ntot;
    StarsNum = NumPart - N_gas;
    accStars = malloc(StarsNum * sizeof(int));
    IdStars = malloc(StarsNum * sizeof(int));

    for (i = N_gas; i < NumPart; i++)
    {
        accStars[i - N_gas] = P[i].Accretion;
        IdStars[i - N_gas] = P[i].ID;

        P[i].Accretion = 0;
    }

    int *displs, *rcounts;
    rcounts = malloc(sizeof(int) * NTask);
    MPI_Gather(&StarsNum, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0, ntot = 0; i < NTask; i++)
        ntot += rcounts[i];

    displs = (int*)malloc(sizeof(int) * NTask);
    displs[0] = 0;

    for (i = 1; i < NTask; ++i)
        displs[i] = displs[i - 1] + rcounts[i - 1];
    int *AllAccStars, *AllIdStars;
    AllAccStars = malloc(sizeof(int) * ntot);
    AllIdStars = malloc(sizeof(int) * ntot);

    MPI_Gatherv(
        accStars, StarsNum, MPI_INT, AllAccStars, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(
        IdStars, StarsNum, MPI_INT, AllIdStars, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (ThisTask == 0)
    {   
        
        int* SortStarsAcc,minID;        
        SortStarsAcc = malloc(sizeof(int) * ntot);
        minID=AllIdStars[0];
        for (i = 1; i < ntot; i++)
 	 if(minID>AllIdStars[i])
  	  minID=AllIdStars[i];
        for (i = 0; i < ntot; i++)
        {           
            SortStarsAcc[AllIdStars[i] - minID] = AllAccStars[i];
        }

        char buf[200];
        sprintf(buf, "%s%s", All.OutputDir, "accretion.dat");

        FILE* FdAccretion;
        FdAccretion = fopen(buf, "a");

        fprintf(FdAccretion, "%f", All.Time);
        for (i = 0; i < ntot; i++)
        {

            fprintf(FdAccretion, " %i", SortStarsAcc[i]);
        }

        fprintf(FdAccretion, "\n");

        fclose(FdAccretion);

        free(SortStarsAcc);
    }

    free(displs);
    free(AllAccStars);

    free(AllIdStars);
    free(rcounts);

    free(accStars);
    free(IdStars);
}


/*! This function fills the write buffer with particle data. New output blocks
 *  can in principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int* startindex, int pc, int type)
{
    int n, k, pindex;
    float* fp;
    int* ip;
#ifdef PMGRID
    double dt_gravkick_pm = 0;
#endif
    double dt_gravkick, dt_hydrokick;

#ifdef PMGRID
    dt_gravkick_pm
        = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif
    fp = CommBuffer;
    ip = CommBuffer;

    pindex = *startindex;

    switch (blocknr)
    {
        case IO_POS: /* positions */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    for (k = 0; k < 3; k++)
                    {
                        fp[k] = P[pindex].Pos[k];
                    }
                    n++;
                    fp += 3;
                }
            break;

        case IO_VEL: /* velocities */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    dt_gravkick = dt_hydrokick
                        = (All.Ti_Current - (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2)
                        * All.Timebase_interval;

                    for (k = 0; k < 3; k++)
                    {
                        fp[k] = P[pindex].Vel[k] + P[pindex].GravAccel[k] * dt_gravkick;
                        if (P[pindex].Type < 2)
                        {
                            fp[k] += SphP[pindex].DragAccel[k] * dt_hydrokick;
                            if (P[pindex].Type == 0)
                                fp[k] += SphP[pindex].HydroAccel[k] * dt_hydrokick;
                        }
                    }
                    n++;
                    fp += 3;
                }
            break;

        case IO_ID: /* particle ID */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *ip++ = P[pindex].ID;
                    n++;
                }
            break;

        case IO_MASS: /* particle mass */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = P[pindex].Mass;
                    n++;
                }
            break;

        case IO_U: /* internal energy */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
#ifdef ISOTHERM_EQS
                    *fp++ = SphP[pindex].Entropy;
#else
                    *fp++ = dmax(All.MinEgySpec, SphP[pindex].Entropy / GAMMA_MINUS1
                            * pow(SphP[pindex].Density, GAMMA_MINUS1));
#endif
                    n++;
                }
            break;

        case IO_RHO: /* density */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Density;
                    n++;
                }
            break;

        case IO_HSML: /* SPH smoothing length */
            for (n = 0; n < pc; pindex++)
                if (P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Hsml;
                    n++;
                }
            break;


        case IO_POT: /* gravitational potential */
            break;

        case IO_ACCEL: /* acceleration */
            break;

        case IO_DTENTR: /* rate of change of entropy */
            break;

        case IO_TSTP: /* timestep  */
            break;
    }

    *startindex = pindex;
}

/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file. If one wants to add a new output-block, this
 *  function should be augmented accordingly.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
    int bytes_per_blockelement = 0;

    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
            bytes_per_blockelement = 3 * sizeof(float);
            break;

        case IO_ID:
            bytes_per_blockelement = sizeof(int);
            break;

        case IO_MASS:
        case IO_U:
        case IO_RHO:
        case IO_HSML:
        case IO_POT:
        case IO_DTENTR:
        case IO_TSTP:
            bytes_per_blockelement = sizeof(float);
            break;
    }

    return bytes_per_blockelement;
}

/*! This function returns the type of the data contained in a given block of
 *  the output file. If one wants to add a new output-block, this function
 *  should be augmented accordingly.
 */
int get_datatype_in_block(enum iofields blocknr)
{
    int typekey;

    switch (blocknr)
    {
        case IO_ID:
            typekey = 0; /* native int */
            break;

        default:
            typekey = 1; /* native float */
            break;
    }

    return typekey;
}

/*! This function informs about the number of elements stored per particle for
 *  the given block of the output file. If one wants to add a new
 *  output-block, this function should be augmented accordingly.
 */
int get_values_per_blockelement(enum iofields blocknr)
{
    int values = 0;

    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
            values = 3;
            break;

        case IO_ID:
        case IO_MASS:
        case IO_U:
        case IO_RHO:
        case IO_HSML:
        case IO_POT:
        case IO_DTENTR:
        case IO_TSTP:
            values = 1;
            break;
    }

    return values;
}


/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array. If one wants to
 *  add a new output-block, this function should be augmented accordingly.
 */
int get_particles_in_block(enum iofields blocknr, int* typelist)
{
    int i, nall, ntot_withmasses, ngas;

    nall = 0;
    ntot_withmasses = 0;

    for (i = 0; i < 3; i++)
    {
        typelist[i] = 0;

        if (header.npart[i] > 0)
        {
            nall += header.npart[i];
            typelist[i] = 1;
        }
        ntot_withmasses += header.npart[i];
    }

    ngas = header.npart[0];
    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_TSTP:
        case IO_ID:
        case IO_POT:
            return nall;
            break;

        case IO_RHO:
        case IO_MASS:
            for (i = 0; i < 3; i++)
            {
                typelist[i] = 0;
                // if(All.MassTable[i] == 0 && header.npart[i] > 0)
                if (header.npart[i] > 0)
                    typelist[i] = 1;
            }
            return ntot_withmasses;
            break;

        case IO_U:
        case IO_HSML:
        case IO_DTENTR:
            for (i = 1; i < 3; i++)
                typelist[i] = 0;
            return ngas;
            break;
    }

    endrun(212);
    return 0;
}

/*! This function tells whether or not a given block in the output file is
 *  present, depending on the type of simulation run and the compile-time
 *  options. If one wants to add a new output-block, this function should be
 *  augmented accordingly.
 */
int blockpresent(enum iofields blocknr)
{


    if (blocknr == IO_POT)
        return 0;
    if (blocknr == IO_ACCEL)
        return 0;
    if (blocknr == IO_DTENTR)
        return 0;
    if (blocknr == IO_TSTP)
        return 0;
    if (blocknr == IO_HSML)
        return 0;

    return 1; /* default: present */
}


/*! This function writes an actual snapshot file containing the data from
 *  processors 'writeTask' to 'lastTask'. 'writeTask' is the one that actually
 *  writes.  Each snapshot file contains a header first, then particle
 *  positions, velocities and ID's.  Particle masses are written only for
 *  those particle types with zero entry in MassTable.  After that, first the
 *  internal energies u, and then the density is written for the SPH
 *  particles.  If cooling is enabled, mean molecular weight and neutral
 *  hydrogen abundance are written for the gas particles. This is followed by
 *  the SPH smoothing length and further blocks of information, depending on
 *  included physics and compile-time flags.  If HDF5 is used, the header is
 *  stored in a group called "/Header", and the particle data is stored
 *  separately for each particle type in groups calles "/PartType0",
 *  "/PartType1", etc. The sequence of the blocks is unimportant in this case.
 */
void write_file(char* fname, int writeTask, int lastTask)
{
    int type, bytes_per_blockelement, npart, typelist[3];
    int n_for_this_task, n, p, pc, offset = 0, task;
    int blockmaxlen, ntot_type[3], nn[3];
    enum iofields blocknr;
    int blksize;
    MPI_Status status;
    FILE* fd = 0;
#define SKIP                                                                                       \
    {                                                                                              \
        my_fwrite(&blksize, sizeof(int), 1, fd);                                                   \
    }
    /* determine particle numbers of each type in file */

    if (ThisTask == writeTask)
    {
        for (n = 0; n < 3; n++)
            ntot_type[n] = n_type[n];

        for (task = writeTask + 1; task <= lastTask; task++)
        {
            MPI_Recv(&nn[0], 3, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
            for (n = 0; n < 3; n++)
                ntot_type[n] += nn[n];
        }

        for (task = writeTask + 1; task <= lastTask; task++)
            MPI_Send(&ntot_type[0], 3, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&n_type[0], 3, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
        MPI_Recv(&ntot_type[0], 3, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }
    /* fill file header */

    for (n = 0; n < 3; n++)
    {
        header.npart[n] = ntot_type[n];
        header.npartTotal[n] = (unsigned int)ntot_type_all[n];
    }
    header.time = All.Time;
    /* open file and write header */

    if (ThisTask == writeTask)
    {
        if (!(fd = fopen(fname, "w")))
        {
            printf("can't open file `%s' for writing snapshot.\n", fname);
            endrun(123);
        }
        blksize = sizeof(header);
        SKIP;
        my_fwrite(&header, sizeof(header), 1, fd);
        SKIP;
    }

    for (blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
        if (blockpresent(blocknr))
        {
            bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

            blockmaxlen = ((int)(All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

            npart = get_particles_in_block(blocknr, &typelist[0]);

            if (npart > 0)
            {
                if (ThisTask == writeTask)
                {
                    blksize = npart * bytes_per_blockelement;
                    SKIP;
                }
                for (type = 0; type < 3; type++)
                {
                    if (typelist[type])
                    {
                        for (task = writeTask, offset = 0; task <= lastTask; task++)
                        {
                            if (task == ThisTask)
                            {
                                n_for_this_task = n_type[type];

                                for (p = writeTask; p <= lastTask; p++)
                                    if (p != ThisTask)
                                        MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK,
                                            MPI_COMM_WORLD);
                            }
                            else
                                MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK,
                                    MPI_COMM_WORLD, &status);

                            while (n_for_this_task > 0)
                            {
                                pc = n_for_this_task;

                                if (pc > blockmaxlen)
                                    pc = blockmaxlen;

                                if (ThisTask == task)
                                    fill_write_buffer(blocknr, &offset, pc, type);

                                if (ThisTask == writeTask && task != writeTask)
                                    MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE,
                                        task, TAG_PDATA, MPI_COMM_WORLD, &status);

                                if (ThisTask != writeTask && task == ThisTask)
                                    MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE,
                                        writeTask, TAG_PDATA, MPI_COMM_WORLD);

                                if (ThisTask == writeTask)
                                {
                                    my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                                }

                                n_for_this_task -= pc;
                            }
                        }
                    }
                }

                if (ThisTask == writeTask)
                {
                    SKIP;
                }
            }
        }
    }

    if (ThisTask == writeTask)
    {
        fclose(fd);
    }
}

/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void* ptr, size_t size, size_t nmemb, FILE* stream)
{
    size_t nwritten;

    if ((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
        printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
        fflush(stdout);
        endrun(777);
    }
    return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void* ptr, size_t size, size_t nmemb, FILE* stream)
{
    size_t nread;

    if ((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
        printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
        fflush(stdout);
        endrun(778);
    }
    return nread;
}
