#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

/*! \file timestep.c
 *  \brief routines for 'kicking' particles in momentum space and assigning new timesteps
 */

static double dt_displacement = 0;


/*! This function advances the system in momentum space, i.e. it does apply
 *  the 'kick' operation after the forces have been computed. Additionally, it
 *  assigns new timesteps to particles. At start-up, a half-timestep is
 *  carried out, as well as at the end of the simulation. In between, the
 *  half-step kick that ends the previous timestep and the half-step kick for
 *  the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
    int i, j, no, ti_step, ti_min, tend, tstart;
    double dt_entr2, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2;
    double aphys;
#ifndef SAVETEMP_LTR
    double minentropy, dt_entr;
#endif
    FLOAT dv[3];
    dt_displacement = All.MaxSizeTimestep;
    for (i = 0; i < NumPart; i++)
    {
        if (P[i].Ti_endstep == All.Ti_Current)
        {
            ti_step = get_timestep(i, &aphys, 0);

            /* make it a power 2 subdivision */
            ti_min = TIMEBASE;
            while (ti_min > ti_step)
                ti_min >>= 1;
            ti_step = ti_min;
            if (ti_step > (P[i].Ti_endstep - P[i].Ti_begstep)) /* timestep wants to increase */
            {
                if (((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
                    ti_step = P[i].Ti_endstep - P[i].Ti_begstep; /* leave at old step */
            }
            if (All.Ti_Current == TIMEBASE) /* we here finish the last timestep. */
                ti_step = 0;

            if ((TIMEBASE - All.Ti_Current) < ti_step) /* check that we don't run beyond the end */
                ti_step = TIMEBASE - All.Ti_Current;

            tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2; /* midpoint of old step */
            tend = P[i].Ti_endstep + ti_step / 2; /* midpoint of new step */
            dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
#ifndef SAVETEMP_LTR
            dt_entr = (tend - tstart) * All.Timebase_interval;
#endif
            dt_gravkick2 = dt_hydrokick2 = dt_entr2
                = (tend - P[i].Ti_endstep) * All.Timebase_interval;
            P[i].Ti_begstep = P[i].Ti_endstep;
            P[i].Ti_endstep = P[i].Ti_begstep + ti_step;
            /* do the kick */

            for (j = 0; j < 3; j++)
            {
                dv[j] = P[i].GravAccel[j] * dt_gravkick;
                P[i].Vel[j] += dv[j];
            }

            if (P[i].Type < 2) /* SPH stuff */
            {
                for (j = 0; j < 3; j++)
                {
                    dv[j] = SphP[i].DragAccel[j] * dt_hydrokick;
                    P[i].Vel[j] += dv[j];
                    SphP[i].VelPred[j] = P[i].Vel[j] - dt_gravkick2 * P[i].GravAccel[j]
                        - SphP[i].DragAccel[j] * dt_hydrokick2;
                }
                if (P[i].Type == 0) /* SPH stuff */
                {
                    for (j = 0; j < 3; j++)
                    {
                        dv[j] = SphP[i].HydroAccel[j] * dt_hydrokick;
                        P[i].Vel[j] += dv[j];
                        SphP[i].VelPred[j] = P[i].Vel[j] - dt_gravkick2 * P[i].GravAccel[j]
                            - dt_hydrokick2 * SphP[i].HydroAccel[j]
                            - SphP[i].DragAccel[j] * dt_hydrokick2;
                    }
/* In case of cooling, we prevent that the entropy (and
   hence temperature decreases by more than a factor 0.5 */
#ifndef SAVETEMP_LTR
                    if (SphP[i].DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
                        SphP[i].Entropy += SphP[i].DtEntropy * dt_entr;
                    else
                        SphP[i].Entropy *= 0.5;

                    if (All.MinEgySpec)
                    {
                        minentropy
                            = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].Density, GAMMA_MINUS1);
                        if (SphP[i].Entropy < minentropy)
                        {
                            SphP[i].Entropy = minentropy;
                            SphP[i].DtEntropy = 0;
                        }
                    }

                    /* In case the timestep increases in the new step, we
                       make sure that we do not 'overcool' when deriving
                       predicted temperatures. The maximum timespan over
                       which prediction can occur is ti_step/2, i.e. from
                       the middle to the end of the current step */

                    dt_entr = ti_step / 2 * All.Timebase_interval;
                    if (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
                        SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
#else
                    double meanweight, R, Tmid, Tz;
                    meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);
                    R = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1]);
                    Tmid = Tcoef * All.CentralTemperature
                        * pow(All.CentralRadius * SOLAR_RADIUS * AU / All.UnitLength_in_cm / R,
                               Pindex);
                    Tz = Tmid;

                    if (Aindex > 0)
                    {
                        double hz, Tatm;
                        hz = heightScale * sqrt(BOLTZMANN * GAMMA * Tmid * R * R * R / meanweight
                                               / PROTONMASS / All.G / All.CentralMass)
                            / All.UnitVelocity_in_cm_per_s;
                        Tatm = All.CentralTemperature * pow(All.CentralRadius * SOLAR_RADIUS * AU
                                                                / All.UnitLength_in_cm / (2. * R),
                                                            Aindex);
                        Tz += (Tatm - Tmid) * pow(sin(M_PI * P[i].Pos[2] / 3. / hz), 4);
                        if(fabs(P[i].Pos[2])>(3.*hz))
                           Tz=Tatm;
                    }                    
                    SphP[i].Entropy =  BOLTZMANN * Tz / meanweight / PROTONMASS
                        * All.UnitMass_in_g / All.UnitEnergy_in_cgs;
                    SphP[i].DtEntropy = 0;
#endif
                }
            }

            /* if tree is not going to be reconstructed, kick parent nodes dynamically.
             */
            if (All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
            {
                no = Father[i];
                while (no >= 0)
                {
                    for (j = 0; j < 3; j++)
                        Extnodes[no].vs[j] += dv[j] * P[i].Mass / Nodes[no].u.d.mass;

                    no = Nodes[no].u.d.father;
                }
            }
        }
    }
}

/*! This function normally (for flag==0) returns the maximum allowed timestep
 *  of a particle, expressed in terms of the integer mapping that is used to
 *  represent the total simulated timespan. The physical acceleration is
 *  returned in `aphys'.
 */
int get_timestep(int p, /*!< particle index */
    double* aphys, /*!< acceleration (physical units) */
    int flag /*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */)
{
    double ax, ay, az, ac;
    double dt = 0, dt_courant = 0, dt_accel;
    int ti_step;
    ax = P[p].GravAccel[0];
    ay = P[p].GravAccel[1];
    az = P[p].GravAccel[2];
    if (P[p].Type < 2)
    {
        ax += SphP[p].DragAccel[0];
        ay += SphP[p].DragAccel[1];
        az += SphP[p].DragAccel[2];
        if (P[p].Type == 0)
        {
            ax += SphP[p].HydroAccel[0];
            ay += SphP[p].HydroAccel[1];
            az += SphP[p].HydroAccel[2];
        }
    }
    ac = sqrt(ax * ax + ay * ay + az * az); /* this is now the physical acceleration */
    *aphys = ac;
    if (ac == 0)
        ac = 1.0e-30;
    dt = dt_accel = All.MinSizeTimestep;
    if (P[p].Type < 2)
        dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * SphP[p].Hsml / ac);
    if (P[p].Type == 0)
    {
        dt_courant = 2 * All.CourantFac * SphP[p].Hsml / SphP[p].MaxSignalVel;
        if (dt_courant < dt)
            dt = dt_courant;
    }
    if (dt >= All.MaxSizeTimestep)
        dt = All.MaxSizeTimestep;

    if (dt >= dt_displacement)
        dt = dt_displacement;

    if (dt < All.MinSizeTimestep)
        dt = All.MinSizeTimestep;
    ti_step = dt / All.Timebase_interval;
    if (!(ti_step > 0 && ti_step < TIMEBASE))
    {
        printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
               "We better stop.\n"
               "Task=%d Part-ID=%d type=%d dt=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) "
               "tree=(%g|%g|%g)\n\n",
            ThisTask, (int)P[p].ID, (int)P[p].Type, dt, All.Timebase_interval, ti_step, ac,
            P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1],
            P[p].GravAccel[2]);
        if (P[p].Type < 2)
        {
            printf("drag-frc=(%g|%g|%g)\n", SphP[p].DragAccel[0], SphP[p].DragAccel[1],
                SphP[p].DragAccel[2]);
            if (P[p].Type == 0)
                printf("hydro-frc=(%g|%g|%g) entropy=%g\n", SphP[p].HydroAccel[0], SphP[p].HydroAccel[1],
                    SphP[p].HydroAccel[2],SphP[p].Entropy);
        }
        fflush(stdout);
        endrun(818);
    }

    return ti_step;
}
