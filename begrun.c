#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include "allvars.h"
#include "proto.h"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are
 *  initialized to their proper values.
 */


/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{

    struct global_data_all_processes all;
    read_parameter_file(ParameterFile); /* ... read in parameters for this run */
    allocate_commbuffers(); /* ... allocate buffer-memory for particle
                               exchange during force computation */
    set_units();
#ifndef NOGRAVITY
    open_outputfiles();
#endif

    if (RestartFlag == 0)
    {
        init(); /* ... read in initial model */
    }
    else
    {
        int i;
        all = All; /* save global variables. (will be read from restart file) */
        restart(
            RestartFlag); /* ... read restart file. Note: This also resets
                                     all variables in the struct `All'.
                                     However, during the run, some variables in the parameter
                                     file are allowed to be changed, if desired. These need to
                                     copied in the way below.
                                     Note:  All.PartAllocFactor is treated in restart() separately.
                                   */
        All.MinSizeTimestep = all.MinSizeTimestep;
        All.MaxSizeTimestep = all.MaxSizeTimestep;
        All.BufferSize = all.BufferSize;
        All.BunchSizeForce = all.BunchSizeForce;
        All.BunchSizeDensity = all.BunchSizeDensity;
        All.BunchSizeHydro = all.BunchSizeHydro;
        All.BunchSizeDomain = all.BunchSizeDomain;
        All.BunchSizeDrag = all.BunchSizeDrag;
        All.TimeBetSnapshot = all.TimeBetSnapshot;
        All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
        All.ErrTolForceAcc = all.ErrTolForceAcc;
        All.ErrTolTheta = all.ErrTolTheta;
        All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
        All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;
        All.DesNumNgb = all.DesNumNgb;
        All.MaxNumNgbDeviation = all.MaxNumNgbDeviation;
        All.ArtBulkViscConst = all.ArtBulkViscConst;
        All.CourantFac = all.CourantFac;
        All.InitGasTemp=all.InitGasTemp;
        All.MinGasTemp=all.MinGasTemp;
        for(i=0;i<3;i++)
           All.ForceSoftening[i] = all.ForceSoftening[i];
        All.DustGasMechStep = all.DustGasMechStep;
        All.DustGasMechAngle = all.DustGasMechAngle;
        All.DustSize = all.DustSize/All.UnitLength_in_cm;
        All.DustRho = all.DustRho/All.UnitDensity_in_cgs;
        All.HsmlConstant = all.HsmlConstant;
        All.MinGasHsml = all.MinGasHsml;
        All.BarrierDistance = all.BarrierDistance;
        All.AccretionRadius = all.AccretionRadius;
        All.StarMassIndicator = all.StarMassIndicator;
        All.CentralTemperature = all.CentralTemperature;
        All.CentralRadius = all.CentralRadius;
        strcpy(All.OutputDir, all.OutputDir);
        strcpy(All.RestartFile, all.RestartFile);
        strcpy(All.SnapshotFileBase, all.SnapshotFileBase);
        if (All.TimeMax != all.TimeMax)
            readjust_timebase(All.TimeMax, all.TimeMax);
    }
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

}

/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    if (All.GravityConstantInternal == 0)
        All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g
            * pow(All.UnitTime_in_s, 2);
    else
        All.G = All.GravityConstantInternal;
    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
    All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
    All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
    All.UnitEnergy_in_cgs
        = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);
    /* convert some physical input parameters to internal units */
    if (ThisTask == 0)
    {
        printf("\nG (internal units) = %g\n", All.G);
        printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
        printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
        printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
        printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
        printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
        printf("\n");
    }

#ifdef ISOTHERM_EQS
    All.MinEgySpec = 0;
#else
    double meanweight;
    meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC); /* note: we assume neutral gas here */
    All.MinEgySpec
        = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
    All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif
}

/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
    char mode[2], buf[200];

    if (ThisTask != 0) /* only the root processor writes to the log files */
        return;

    if (RestartFlag == 0)
        strcpy(mode, "w");
    else
        strcpy(mode, "a");
    FILE* FdAccretion;
    sprintf(buf, "%s%s", All.OutputDir, "accretion.dat");
    if (!(FdAccretion = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }
    fclose(FdAccretion);
}

/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char* fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300
    FILE *fd, *fdout;
    char buf[200], buf1[200], buf2[200], buf3[400];
    int i, j, nt;
    int id[MAXTAGS];
    void* addr[MAXTAGS];
    char tag[MAXTAGS][50];
    int errorFlag = 0;

    if (sizeof(long long) != 8)
    {
        if (ThisTask == 0)
            printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
        endrun(0);
    }

    if (sizeof(int) != 4)
    {
        if (ThisTask == 0)
            printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
        endrun(0);
    }

    if (sizeof(float) != 4)
    {
        if (ThisTask == 0)
            printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
        endrun(0);
    }

    if (sizeof(double) != 8)
    {
        if (ThisTask == 0)
            printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
        endrun(0);
    }


    if (ThisTask == 0) /* read parameter file on process 0 */
    {
        nt = 0;

        strcpy(tag[nt], "InitCondFile");
        addr[nt] = All.InitCondFile;
        id[nt++] = STRING;

        strcpy(tag[nt], "OutputDir");
        addr[nt] = All.OutputDir;
        id[nt++] = STRING;

        strcpy(tag[nt], "SnapshotFileBase");
        addr[nt] = All.SnapshotFileBase;
        id[nt++] = STRING;

        strcpy(tag[nt], "RestartFile");
        addr[nt] = All.RestartFile;
        id[nt++] = STRING;

        strcpy(tag[nt], "TimeOfFirstSnapshot");
        addr[nt] = &All.TimeOfFirstSnapshot;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "TimeBegin");
        addr[nt] = &All.TimeBegin;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "TimeMax");
        addr[nt] = &All.TimeMax;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "TimeBetSnapshot");
        addr[nt] = &All.TimeBetSnapshot;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
        addr[nt] = &All.UnitVelocity_in_cm_per_s;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "UnitLength_in_cm");
        addr[nt] = &All.UnitLength_in_cm;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "UnitMass_in_g");
        addr[nt] = &All.UnitMass_in_g;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "TreeDomainUpdateFrequency");
        addr[nt] = &All.TreeDomainUpdateFrequency;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "ErrTolIntAccuracy");
        addr[nt] = &All.ErrTolIntAccuracy;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "ErrTolTheta");
        addr[nt] = &All.ErrTolTheta;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "ErrTolForceAcc");
        addr[nt] = &All.ErrTolForceAcc;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "MinGasHsml");
        addr[nt] = &All.MinGasHsml;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "HsmlConstant");
        addr[nt] = &All.HsmlConstant;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "MaxSizeTimestep");
        addr[nt] = &All.MaxSizeTimestep;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "MinSizeTimestep");
        addr[nt] = &All.MinSizeTimestep;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "ArtBulkViscConst");
        addr[nt] = &All.ArtBulkViscConst;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "CourantFac");
        addr[nt] = &All.CourantFac;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "DesNumNgb");
        addr[nt] = &All.DesNumNgb;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "MaxNumNgbDeviation");
        addr[nt] = &All.MaxNumNgbDeviation;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "TypeOfOpeningCriterion");
        addr[nt] = &All.TypeOfOpeningCriterion;
        id[nt++] = INT;

        strcpy(tag[nt], "SofteningStars");
        addr[nt] = &All.ForceSoftening[2];
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "DustGasMechStep");
        addr[nt] = &All.DustGasMechStep;
        id[nt++] = DOUBLE;
 
        strcpy(tag[nt], "DustGasMechAngle");
        addr[nt] = &All.DustGasMechAngle;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "DustSize");
        addr[nt] = &All.DustSize;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "DustRho");
        addr[nt] = &All.DustRho;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "BufferSize");
        addr[nt] = &All.BufferSize;
        id[nt++] = INT;

        strcpy(tag[nt], "PartAllocFactor");
        addr[nt] = &All.PartAllocFactor;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "TreeAllocFactor");
        addr[nt] = &All.TreeAllocFactor;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "GravityConstantInternal");
        addr[nt] = &All.GravityConstantInternal;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "InitGasTemp");
        addr[nt] = &All.InitGasTemp;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "MinGasTemp");
        addr[nt] = &All.MinGasTemp;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "AccretionRadius");
        addr[nt] = &All.AccretionRadius;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "BarrierDistance");
        addr[nt] = &All.BarrierDistance;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "StarMassIndicator");
        addr[nt] = &All.StarMassIndicator;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "CentralTemperature");
        addr[nt] = &All.CentralTemperature;
        id[nt++] = DOUBLE;

        strcpy(tag[nt], "CentralRadius");
        addr[nt] = &All.CentralRadius;
        id[nt++] = DOUBLE;

        if ((fd = fopen(fname, "r")))
        {
            sprintf(buf, "%s%s", fname, "-usedvalues");
            if (!(fdout = fopen(buf, "w")))
            {
                printf("error opening file '%s' \n", buf);
                errorFlag = 1;
            }
            else
            {
                while (!feof(fd))
                {
                    *buf = 0;
                    if (!fgets(buf, 200, fd))
                    {
                        printf("fgets problem ...");
                        continue;
                    }
                    if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                        continue;

                    if (buf1[0] == '%')
                        continue;

                    for (i = 0, j = -1; i < nt; i++)
                        if (strcmp(buf1, tag[i]) == 0)
                        {
                            j = i;
                            tag[i][0] = 0;
                            break;
                        }

                    if (j >= 0)
                    {
                        switch (id[j])
                        {
                            case DOUBLE:
                                *((double*)addr[j]) = atof(buf2);
                                fprintf(fdout, "%-35s%g\n", buf1, *((double*)addr[j]));
                                break;
                            case STRING:
                                strcpy(addr[j], buf2);
                                fprintf(fdout, "%-35s%s\n", buf1, buf2);
                                break;
                            case INT:
                                *((int*)addr[j]) = atoi(buf2);
                                fprintf(fdout, "%-35s%d\n", buf1, *((int*)addr[j]));
                                break;
                        }
                    }
                    else
                    {
                        fprintf(stdout,
                            "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
                            fname, buf1);
                        errorFlag = 1;
                    }
                }
                fclose(fd);
                fclose(fdout);

                i = strlen(All.OutputDir);
                if (i > 0)
                    if (All.OutputDir[i - 1] != '/')
                        strcat(All.OutputDir, "/");

                sprintf(buf1, "%s%s", fname, "-usedvalues");
                sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
                sprintf(buf3, "cp %s %s", buf1, buf2);
                if (system(buf3))
                {
                    printf("system problem ...");
                    errorFlag = 2;
                }
            }
        }
        else
        {
            printf("\nParameter file %s not found.\n\n", fname);
        }

        if (errorFlag != 2)
            for (i = 0; i < nt; i++)
            {
                if (*tag[i])
                {
                    printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i],
                        fname);
                    errorFlag = 1;
                }
            }
    }

    MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (errorFlag)
    {
        MPI_Finalize();
        exit(0);
    }

    /* now communicate the relevant parameters to the other processes */
    MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}

/*! If a restart from restart-files is carried out where the TimeMax
 *  variable is increased, then the integer timeline needs to be
 *  adjusted. The approach taken here is to reduce the resolution of the
 *  integer timeline by factors of 2 until the new final time can be
 *  reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
    int i;
    long long ti_end;

    if (ThisTask == 0)
    {
        printf("\nAll.TimeMax has been changed in the parameterfile\n");
        printf("Need to adjust integer timeline\n\n\n");
    }

    if (TimeMax_new < TimeMax_old)
    {
        if (ThisTask == 0)
            printf("\nIt is not allowed to reduce All.TimeMax\n\n");
        endrun(556);
    }
    ti_end = (TimeMax_new - All.TimeBegin) / All.Timebase_interval;
    while (ti_end > TIMEBASE)
    {
        All.Timebase_interval *= 2.0;

        ti_end /= 2;
        All.Ti_Current /= 2;
        for (i = 0; i < NumPart; i++)
        {
            P[i].Ti_begstep /= 2;
            P[i].Ti_endstep /= 2;
        }
    }

    All.TimeMax = TimeMax_new;
}
