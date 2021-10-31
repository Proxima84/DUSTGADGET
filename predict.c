#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"


/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly.
 */
void move_particles(int time0, int time1)
{
    int i, j;
    double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;
    dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    for (i = 0; i < NumPart; i++)
    {
        for (j = 0; j < 3; j++)
            P[i].Pos[j] += P[i].Vel[j] * dt_drift;

        if (P[i].Type == 0)
        {
            for (j = 0; j < 3; j++)
                SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick
                    + SphP[i].HydroAccel[j] * dt_hydrokick + SphP[i].DragAccel[j] * dt_hydrokick;
#ifndef HSMLCONSTANT
            if (SphP[i].Hsml < All.MinGasHsml)
                SphP[i].Hsml = All.MinGasHsml;
#endif
            dt_entr = (time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

            SphP[i].Pressure
                = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
        }
        else if (P[i].Type == 1)
        {
            for (j = 0; j < 3; j++)
                SphP[i].VelPred[j]
                    += P[i].GravAccel[j] * dt_gravkick + SphP[i].DragAccel[j] * dt_hydrokick;
        }
    }
#ifdef SHEREBORDER
 double Rgas,Vgas,Rdust,Vdust;
 Vgas=0;
 Rgas=0;
 Rdust=0;
 Vdust=0;
 for (i = 0; i < NumPart; i++)
 {
  if(P[i].Type==0)
  {
   if(P[i].ID!=0)
   {
    Vgas+=P[i].Pos[0]*SphP[i].VelPred[0]+P[i].Pos[1]*SphP[i].VelPred[1]+P[i].Pos[2]*SphP[i].VelPred[2];
    Rgas+=P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2];
   }
  }
  else if(P[i].Type==1)
  {
   if(P[i].ID!=0)
   {
    Vdust+=P[i].Pos[0]*SphP[i].VelPred[0]+P[i].Pos[1]*SphP[i].VelPred[1]+P[i].Pos[2]*SphP[i].VelPred[2];
    Rdust+=P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2];
   }
  }
 }
 for (i = 0; i < NumPart; i++)
 {
  if(P[i].Type==0)
  { 
   if(P[i].ID==0)
     for (j = 0; j < 3; j++)     
       SphP[i].VelPred[j] = Vgas/Rgas*P[i].Pos[j];                       
  }
  else if(P[i].Type==1)
  {
    if(P[i].ID==0)
     for (j = 0; j < 3; j++)     
       SphP[i].VelPred[j] = Vdust/Rdust*P[i].Pos[j];       
  }
 }
#endif

    /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
    if (All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
        for (i = 0; i < Numnodestree; i++)
            for (j = 0; j < 3; j++)
                Nodes[All.MaxPart + i].u.d.s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

        force_update_len();

        force_update_pseudoparticles();
    }
}

