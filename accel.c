#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file accel.c
 *  \brief driver routine to carry out force computation
 */

/*! This routine computes the accelerations for all active particles.
 *  First, the gravitational tree forces
 *  are computed. This constructs the tree, if needed.
 *
 *  If gas particles are present, the density-loop for active SPH particles
 *  is carried out. This includes an iteration on the correct number of
 *  neighbours.  Finally, the hydrodynamical forces are added. If dust
 *  particles are present the drag forces are calculated
 */
void compute_accelerations(int mode)
{
    gravity_tree(); /* computes gravity accel. */
    if (All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
        gravity_tree(); /* For the first timestep, we redo it
                         * to allow usage of relative opening
                         * criterion for consistent accuracy.
                         */
    density(); /* computes density, and pressure */
    force_update_hmax(); /* tell the tree nodes the new SPH smoothing length such that they are
                            guaranteed to hold the correct max(Hsml) */
    hydro_force(); /* adds hydrodynamical accelerations and computes viscous entropy injection  */
    if (All.Time > 0)
        interaction(); /* adds drag accelerations of gas and dust particles  */
}

