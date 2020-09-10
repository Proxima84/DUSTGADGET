#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file main.c
 *  \brief start of the program
 */

/*!
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0.  Then begrun() is called, which sets up
 *  the simulation either from IC's or from restart files.  Finally,
 *  run() is started, the main simulation loop, which iterates over
 *  the timesteps.
 */
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    for (PTask = 0; NTask > (1 << PTask); PTask++)
        ;

    if (argc < 2)
    {
        if (ThisTask == 0)
        {
            printf("Parameters are missing.\n");
            printf("Call with <ParameterFile> [<RestartFlag>]\n");
        }
        endrun(0);
    }

    strcpy(ParameterFile, argv[1]);

    if (argc >= 3)
        RestartFlag = atoi(argv[2]);
    else
        RestartFlag = 0;
    begrun(); /* set-up run  */
    run(); /* main simulation loop */

    MPI_Finalize(); /* clean up & finalize MPI */

    return 0;
}

/* ----------------------------------------------------------------------
   The rest of this file contains documentation for compiling and
   running the code, in a format appropriate for 'doxygen'.
   ----------------------------------------------------------------------
 */

/*! \mainpage Reference documentation for GADGET-2

\author Volker Springel \n
        Max-Planck-Institute for Astrophysics \n
        Karl-Schwarzschild-Str. 1 \n
        85740 Garching \n 
        Germany \n
        volker@mpa-garching.mpg.de \n

\n

\section prelim Getting started

GADGET-2 is a massively parallel code for hydrodynamical cosmological
simulations. It is a flexible code that can be applied to a variety of
different types of simulations, offering a number of sophisticated
simulation algorithms.

A full account of the numerical algorithms employed by the code is given
in the accompanying code paper, and detailed instructions for usage of
the code are given in the included code documentation.

This html-document serves as a cross-referenced documentation of the
source code itself - in fact, using the doxygen tool, the html-pages
have been produced from comments inlined in the source code. Apart from
the source-code documentation, a brief guide to code compilation is
given below, and under <b>Related Pages (see link on top)</b> you can find an
explanation of GADGET's parameterfile and a short guide to compile-time
options of the code.


\section install Compilation 

GADGET-2 needs the following non-standard libraries for compilation:

- \b MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g.  MPICH
  (http://www-unix.mcs.anl.gov/mpi/mpich/) or LAM
  (http://www.lam-mpi.org/).

Note that if any of the above libraries is not installed in standard
locations on your system, the \ref Gadget-Makefile "Makefile" provided with
the code may need slight adjustments. Similarly, compiler options,
particularly with respect to optimisations, may need adjustment to the
C-compiler that is used. Finally, the \ref Gadget-Makefile "Makefile"
contains a number of compile-time options that need to be set appropriately
for the type of simulation that is simulated.

The provided makefile is compatible with GNU-make, i.e. typing \b make or
\b gmake should then build the executable <b>Gadget2</b>.  If your site
does not have GNU-make, get it, or write your own makefile.

\section howtorun Running the code

In order to start the simulation code, a \ref parameterfile "parameterfile"
needs to be specified. An additional optional numerical parameter can be
used to signal whether a continuation from a set of restart files, or from
a snapshot file, is desired. A typical command to start the code looks like
the following: \n \n

 <b> mpirun -np 8 ./Gadget2 <parameterfile> [restartflag]</b> \n \n 

This would start the code using 8 processors, assuming that the parallel
environment uses the <em>mpirun</em> command to start MPI
applications. Depending on the operating system, other commands may be
required for this task, e.g. <em>poe</em> on IBM/AIX machines. Note that
the code can in principle be started using an arbitrary number of
processors, but the communication algorithms will be most efficient for
powers of 2. It is also possible to use a single processor only, in
which case the code behaves like a serial code, except that GADGET-2
will still go through some of the overhead induced by the
parallelization algorithms, so the code will not quite reach the same
performance as an optimum serial solution in this case.


The optional <em>restartflag</em> can have the values 0, 1, or 2, only. "1"
signals a continuation from restart files, while "2" can be used to restart
from a snapshot file produced by the code. If omitted (equivalent to the
default of "0"), the code starts from initial conditions.

*/

/*! \page parameterfile  Parameterfile of GADGET-2

The parameterfile for GADGET-2 is a simple text file, consisting of pairs of
tags and values. For each parameter, a separate line needs to be specified,
first listing the name (tag) of the parameter, and then the assigned value,
separated by whitespace. It is allowed to add further text behind the
assigned parameter value. The order of the parameters is arbitrary, but
each one needs to occur exactly one time, otherwise an error message will
be produced. Empty lines, or lines beginning with a \%-sign, are ignored and
treated as comments.

 
- \b InitCondFile \n The filename of the initial conditions file. If a
             restart from a snapshot with the "2" option is desired, one
             needs to specify the snapshot file here.

- \b OutputDir \n    Pathname of the output directory of the code.

- \b RestartFile \n Basename of restart-files produced by the code.

- \b SnapshotFileBase \n Basename of snapshot files produced by the code.

- \b TimeBegin \n This sets the starting time of a simulation when the code
             is started from initial conditions. For cosmological
             integrations, the value specified here is taken as the initial
             scale factor.

- \b TimeMax \n This sets the final time for the simulation. The code
             normally tries to run until this time is reached. For
             cosmological integrations, the value given here is the final
             scale factor.

- \b TimeBetSnapshot \n The time interval between two subsequent snapshot
             files in case a file with output times is not specified.  For
             cosmological simulations, this is a multiplicative factor
             applied to the time of the last snapshot, such that the
             snapshots will have a constant logarithmic spacing in the
             scale factor. Otherwise, the parameter is an additive constant
             that gives the linear spacing between snapshot times.

- \b TimeOfFirstSnapshot \n The time of the first desired snapshot file in
             case a file with output times is not specified. For
             cosmological simulations, the value given here is the scale factor
             of the first desired output.

- \b ErrTolIntAccuracy \n This dimensionless parameter controls the
             accuracy of the timestep criterion selected by
             <em>TypeOfTimestepCriterion</em>.

- \b CourantFac \n This sets the value of the Courant parameter used in the
             determination of the hydrodynamical timestep of SPH particles.

- \b MaxSizeTimestep \n This gives the maximum timestep a particle may
             take. This should be set to a sensible value in order to
             protect against too large timesteps for particles with very
             small acceleration. For cosmological simulations, the
             parameter given here is the maximum allowed step in the
             logarithm of the expansion factor. Note that the definition
             of MaxSizeTimestep has <em>changed</em> compared to Gadget-1.1 for cosmological simulations.


- \b MinSizeTimestep \n If a particle requests a timestep smaller than the
             value specified here, the code will normally terminate with a
             warning message. If compiled with the
             <em>NOSTOP_WHEN_BELOW_MINTIMESTEP</em> option, the code will
             instead force the timesteps to be at least as large as
             <em>MinSizeTimestep</em>.

- \b TypeOfOpeningCriterion \n This selects the type of cell-opening
             criterion used in the tree walks. A value of `0' results in
             standard Barnes & Hut, while `1' selects the relative opening
             criterion of GADGET-2.

- \b ErrTolTheta \n This gives the maximum opening angle if the BH
             criterion is used for the tree walk.  If the relative opening
             criterion is used instead, a first force estimate is computed
             using the BH algorithm, which is then recomputed with the
             relative opening criterion.

- \b ErrTolForceAcc \n The accuracy parameter for the relative opening
             criterion for the tree walk.

- \b TreeDomainUpdateFrequency \n The domain decomposition and tree
             construction need not necessarily be done every single
             timestep. Instead, tree nodes can be dynamically updated,
             which is faster. However, the tree walk will become more
             expensive since the tree nodes have to "grow" to keep
             accomodating all particles they enclose. The parameter
             <em>TreeDomainUpdateFrequency</em> controls how often the
             domain decomposition is carried out and the tree is
             reconstructed from scratch. For example, a value of 0.1 means
             that the domain decomposition and the tree are reconstructed
             whenever there have been more than 0.1*N force computations
             since the last reconstruction, where N is the total particle
             number. A value of 0 will reconstruct the tree every timestep.

- \b DesNumNgb \n This sets the desired number of SPH smoothing neighbours. 
                  It is ignored if HSMLCONSTANT preprocessor directive used.

- \b MaxNumNgbDeviation \n This sets the allowed variation of the number of
              neighbours around the target value <em>DesNumNgb</em>. It is 
              ignored if HSMLCONSTANT preprocessor directive used.

- \b ArtBulkViscConst \n This sets the value of the artificial viscosity
             parameter used by GADGET-2.


- \b InitGasTemp \n This sets the initial gas temperature (assuming either
             a mean molecular weight corresponding to full ionization or
             full neutrality, depending on whether the temperature is above
             or below 10^4 K) in Kelvin when initial conditions are
             read. However, the gas temperature is only set to a certain
             temperature if <em>InitGasTemp</em>>0, and if the temperature
             of the gas particles in the initial conditions file is zero,
             otherwise the initial gas temperature is left at the value
             stored in the IC file.

- \b MinGasTemp \n A minimum temperature floor imposed by the code. This
             may be set to zero.

- \b PartAllocFactor \n Each processor allocates space for
             <em>PartAllocFactor</em> times the average number of particles
             per processor. This number needs to be larger than 1 to allow
             the simulation to achieve a good work-load balancing, which
             requires to trade particle-load balance for work-load
             balance. It is good to make <em>PartAllocFactor</em> quite a
             bit larger than 1, but values in excess of 3 will typically
             not improve performance any more. For a value that is too
             small, the code may not be able to succeed in the domain
             decomposition and terminate.

- \b TreeAllocFactor \n To construct the BH-tree for N particles, somewhat
             less than N internal tree-nodes are necessary for `normal'
             particle distributions.  <em>TreeAllocFactor</em> sets the
             number of internal tree-nodes allocated in units of the
             particle number.  By experience, space for 0.65 N internal
             nodes is usually fully sufficient, so a value of 0.7 should
             put you on the safe side.

- \b BufferSize \n This specifies the size (in MByte per processor) of a
             communication buffer used by the code.

- \b UnitLength_in_cm \n This sets the internal length unit in cm/h, where
             H_0 = 100 h km/sec/Mpc. For example, a choice of 3.085678e21
             sets the length unit to 1.0 kpc/h.

- \b UnitMass_in_g \n This sets the internal mass unit in g/h, where H_0 =
             100 h km/sec/Mpc.  For example, a choice of 1.989e43 sets the
             mass unit to 10^10 M_sun/h.

- \b UnitVelocity_in_cm_per_s \n This sets the internal velocity unit in
             cm/sec. For example, a choice of 1e5 sets the velocity unit to
             km/sec.  Note that the specification of
             <em>UnitLength_in_cm</em>, <em>UnitMass_in_g</em>, and
             <em>UnitVelocity_in_cm_per_s</em> also determines the internal
             unit of time.

- \b GravityConstantInternal \n The numerical value of the gravitational
             constant G in internal units depends on the system of units
             you choose. For example, for the choices above, G=43007.1 in
             internal units.  For <em>GravityConstantInternal</em>=0, the
             code calculates the value corresponding to the physical value
             of G automatically.  However, you might want to set G
             yourself.  For example, by specifying
             <em>GravityConstantInternal</em>=1,
             <em>UnitLength_in_cm</em>=1, <em>UnitMass_in_g</em>=1, and
             <em>UnitVelocity_in_cm_per_s</em>=1, one obtains a `natural'
             system of units. Note that the code will nevertheless try to
             use the `correct' value of the Hubble constant in this case,
             so you should not set <em>GravityConstantInternal</em> in
             cosmological integrations.

- \b MinGasHsmlFractional \n This parameter sets the minimum allowed SPH
             smoothing length in units of the gravitational softening
             length of the gas particles. The smoothing length will be
             prevented from falling below this value. When this bound is
             actually reached, the number of smoothing neighbors will
             instead be increased above <em>DesNumNgb</em>. It is 
              ignored if HSMLCONSTANT preprocessor directive used.

- \b HsmlConstant \n The constant smoothing length for the gas and dust particles. It is 
              ignored if HSMLCONSTANT preprocessor directive not used.

- \b SofteningStars \n The Plummer equivalent gravitational softening
             length for particle type 2

//TODO
DustGasMechStep    0.005
DustSize	   0.0001 %cm
DustRho			1 %g/cm^3

AccretionRadius  1.0
BarrierDistance  15.0
StarMassIndicator 0.1
CentralTemperature 6000 % in Kelvin
CentralRadius 1.0 %in solar radius


*/

//TODO

/*! \page Gadget-Makefile  Makefile of GADGET-2

A number of features of GADGET-2 are controlled with compile-time options
in the makefile rather than by the parameterfile. This has been done in
order to allow the generation of highly optimised binaries by the compiler,
even when the underlying source code allows for many different ways to run the
code.

The makefile contains a dummy list of all available compile-time options,
with most of them commented out by default. To activate a certain feature,
the corresponding parameter should be commented in, and given the desired
value, where appropriate. Below, a brief guide to these options is
included.

<b>Important Note:</b> Whenever one of the compile-time options
described below is modified, a full recompilation of the code may be
necessary. To guarantee that this is done when a simple <b>make</b> is
specified, all source files have been specified in the Makefile as being
dependent on the Makefile itself. Alternatively, one can also issue the
command <b>make clean</b>, which will erase all object files, followed
by <b>make</b>.

Note that the above technique has the disadvantage that different
simulations may require different binaries of GADGET-2. If several
simulations are run concurrently, there is hence the danger that a
simulation is started/resumed with the `wrong' binary. Note that while
GADGET-2 checks the plausibility of some of the most important code
options, this is not done for all of them. To minimise the risk of using
the wrong executable for a simulation, it is recommended to produce a
separate executable for each simulation that is run. For example, a good
strategy is to make a copy of the whole code together with its makefile in
the output directory of each simulation run, and then to use this copy to
compile the code and to run the simulation.


\n
\section secmake1 Basic operation mode of code
- \b PERIODIC \n Set this if you want to have periodic boundary conditions.

- \b UNEQUALSOFTENINGS \n Set this if you use particles with different
     gravitational softening lengths.

\n
\section secmake2 Things that are always recommended
- \b PEANOHILBERT \n This is a tuning option. When set, the code will bring
     the particles into Peano-Hilbert order after each domain
     decomposition. This improves cache utilisation and performance.

- \b WALLCLOCK \n If set, a wallclock timer is used by the code to measure
     internal time consumption (see cpu-log file).  Otherwise, a timer that
     measures consumed processor ticks is used.

\n
\section secmake3 TreePM options
- \b PMGRID=128 \n This enables the TreePM method, i.e. the long-range
     force is computed with a PM-algorithm, and the short range force with
     the tree. The parameter has to be set to the size of the mesh that
     should be used, e.g.~64, 96, 128, etc. The mesh dimensions need not
     necessarily be a power of two, but the FFT is fastest for such a
     choice.  Note: If the simulation is not in a periodic box, then a FFT
     method for vacuum boundaries is employed, using a mesh with dimension
     twice that specified by <b>PMGRID</b>.

- \b PLACEHIGHRESREGION=1+8 \n If this option is set (will only work
     together with \b PMGRID), then the long range force is computed in two
     stages: One Fourier-grid is used to cover the whole simulation volume,
     allowing the computation of the large-scale force.  A second Fourier
     mesh is placed on the region occupied by `high-resolution' particles,
     allowing the computation of an intermediate-scale force. Finally, the
     force on very small scales is computed by the tree. This procedure can
     be useful for `zoom-simulations', where the majority of particles (the
     high-res particles) are occupying only a small fraction of the
     volume. To activate this option, the parameter needs to be set to an
     integer that encodes the particle types that make up the high-res
     particles in the form of a bit mask. For example, if types 0, 1, and 4
     are the high-res particles, then the parameter should be set to
     <b>PLACEHIGHRESREGION=1+2+16</b>, i.e. to the sum
     \f$2^0+2^1+2^4\f$. The spatial region covered by the high-res grid is
     determined automatically from the initial conditions. Note: If a
     periodic box is used, the high-res zone is not allowed to intersect the box
     boundaries.

- <b> ENLARGEREGION=1.1</b> \n The spatial region covered by the high-res zone
     normally has a fixed size during the simulation, which initially is
     set to the smallest region that encompasses all high-res
     particles. Normally, the simulation will be interrupted if high-res
     particles leave this region in the course of the run. However, by
     setting this parameter to a value larger than one, the high-res region
     can be expanded on the fly.  For example, setting it to 1.4 will enlarge its
     side-length by 40% in such an event (it remains centred on the high-res
     particles). Hence, with such a setting, the high-res region may expand
     or move by a limited amount. If in addition \b SYNCHRONIZATION is
     activated, then the code will be able to continue even if high-res
     particles leave the initial high-res grid. In this case, the code will
     update the size and position of the grid that is placed onto the
     high-resolution region automatically. To prevent that this potentially
     happens every single PM step, one should nevertheless assign a value
     slightly larger than 1 to \b ENLARGEREGION.

- <b> ASMTH=1.25</b> \n This can be used to override the value assumed for the
     scale that defines the long-range/short-range force-split in the
     TreePM algorithm. The default value is 1.25, in mesh-cells.

- <b> RCUT=4.5</b> \n This can be used to override the maximum radius in which
     the short-range tree-force is evaluated (in case the TreePM algorithm
     is used). The default value is 4.5, given in mesh-cells.

\n 
\section secmake4 Single or double precision
- \b DOUBLEPRECISION \n This makes the code store and compute internal
     particle data in double precision. Note that output files are
     nevertheless written by converting the values that are saved to single
     precision.

- \b DOUBLEPRECISION_FFTW \n If this is set, the code will use the
     double-precision version of FTTW, provided the latter has been
     explicitly installed with a "d" prefix, and NOTYPEPREFIX_FFTW is not
     set. Otherwise the single precision version ("s" prefix) is used.


\n
\section secmake5 Time integration options
- \b SYNCHRONIZATION \n When this is set, particles may only increase their
     timestep if the new timestep will put them into synchronisation with
     the higher time level. This typically means that only on half of the
     timesteps of a particle an increase of its step may occur. Especially
     for TreePM runs, it is usually advisable to set this option.

- \b FLEXSTEPS \n This is an alternative to SYNCHRONIZATION. Particle
     timesteps are here allowed to be integer multiples of the minimum
     timestep that occurs among the particles, which in turn is rounded
     down to the nearest power-of-two devision of the total simulated
     timespan. This option distributes particles more evenly over
     individual system timesteps, particularly once a simulation has run
     for a while, and may then result in a reduction of work-load imbalance
     losses.

- \b PSEUDOSYMMETRIC \n When this option is set, the code will try to
     `anticipate' timestep changes by extrapolating the change of the
     acceleration into the future. This in general improves the long-term
     integration behaviour of periodic orbits, because then the adaptive
     integration becomes more akin to a strictly time reversible
     integrator. Note: This option has no effect if FLEXSTEPS is set.

- \b NOSTOP_WHEN_BELOW_MINTIMESTEP \n If this is activated, the code will
     not terminate when the timestep falls below the value of \b
     MinSizeTimestep specified in the parameterfile. This is useful for
     runs where one wants to enforce a constant timestep for all
     particles. This can be done by activating this option, and by setting
     \b MinSizeTimestep and \b MaxSizeTimestep to an equal value.

- \b NOPMSTEPADJUSTMENT \n When this is set, the long-range timestep for
     the PM force computation is always determined by \b MaxSizeTimeStep.
     Otherwise, it is set to the minimum of \b MaxSizeTimeStep and the
     timestep obtained for the maximum long-range force with an effective
     softening scale equal to the PM smoothing-scale.

\n
\section secmake6  Output options
- \b HAVE_HDF5 \n If this is set, the code will be compiled with support
     for input and output in the HDF5 format. You need to have the HDF5
     libraries and headers installed on your computer for this option to
     work. The HDF5 format can then be selected as format "3" in Gadget's
     parameterfile.

- \b OUTPUTPOTENTIAL \n This will force the code to compute gravitational
     potentials for all particles each time a snapshot file is
     generated. These values are then included in the snapshot files. Note
     that the computation of the values of the potential costs additional
     time.

- \b OUTPUTACCELERATION \n This will include the physical acceleration of
     each particle in snapshot files.

- \b OUTPUTCHANGEOFENTROPY \n This will include the rate of change of
     entropy of gas particles in snapshot files.

- \b OUTPUTTIMESTEP \n This will include the timesteps actually taken by
     each particle in the snapshot files.

\n
\section secmake7 Things for special behaviour
- \b NOGRAVITY \n This switches off gravity. Makes only sense for pure SPH
     simulations in non-expanding space.

- \b NOTREERND \n If this is not set, the tree construction will succeed
     even when there are a few particles at identical locations. This is
     done by `rerouting' particles once the node-size has fallen below
     \f$10^{-3}\f$ of the softening length. When this option is activated,
     this will be suppressed and the tree construction will always fail if
     there are particles at extremely close or identical coordinates.

- \b NOTYPEPREFIX_FFTW \n If this is set, the fftw-header/libraries are
     accessed without type prefix (adopting whatever was chosen as default
     at compile-time of fftw). Otherwise, the type prefix 'd' for
     double-precision is used.
 
- \b LONG_X/Y/Z \n These options can be used together with PERIODIC and
     NOGRAVITY only.  When set, the options define numerical factors that
     can be used to distort the periodic simulation cube into a
     parallelepiped of arbitrary aspect ratio. This can be useful for
     idealized SPH tests.

- \b TWODIMS \n This effectively switches of one dimension in SPH,
     i.e. the code follows only 2d hydrodynamics in the xy-, yz-, or
     xz-plane. This only works with NOGRAVITY, and if all coordinates of
     the third axis are exactly equal. Can be useful for idealized SPH
     tests.

- \b SPH_BND_PARTICLES \n If this is set, particles with a particle-ID
     equal to zero do not receive any SPH acceleration. This can be useful
     for idealized SPH tests, where these particles represent fixed
     "walls".

- \b NOVISCOSITYLIMITER \n If this is set, there is no explicit upper
     limit on the viscosity.  In the default version, this limiter will
     try to protect against possible particle `reflections', which could
     in principle occur if very poor timestepping is used in the
     presence of strong shocks.

- \b COMPUTE_POTENTIAL_ENERGY \n When this option is set, the code will
     compute the gravitational potential energy each time a global
     statistics is computed. This can be useful for testing global energy
     conservation.

- \b ISOTHERM_EQS \n This special option makes the gas behave like an
     isothermal gas with equation of state \f$ P = c_s^2 \rho \f$. The
     sound-speed \f$ c_s \f$ is set by the thermal energy per unit mass in the
     intial conditions, i.e. \f$ c_s^2=u \f$. If the value for \f$ u \f$ is
     zero, then the initial gas temperature in the parameter file is used to
     define the sound speed according to \f$ c_s^2= k\,T/m_p \f$ , where \f$
     m_p \f$ is the proton mass.

- \b ADAPTIVE_GRAVSOFT_FORGAS \n When this option is set, the gravitational
     softening lengths used for gas particles is tied to their SPH smoothing
     length. This can be useful for dissipative collapse simulations. The
     option requires the setting of UNEQUALSOFTENINGS.

- \b SELECTIVE_NO_GRAVITY \n This can be used for special computations where
     one wants to exclude certain particle types from receiving gravitational
     forces. The particle types that are excluded in this fashion are specified
     by a bit mask, in the same as for the PLACEHIGHRESREGION option.

- \b LONGIDS \n If this is set, the code assumes that particle-IDs are
     stored as 64-bit long integers. This is only really needed if you want
     to go beyond ~2 billion particles.
     
\n
\section secmake8 Testing and Debugging options
- \b FORCETEST=0.01 \n This can be set to check the force accuracy of the
     code, and is only included as a debugging option. The option needs to
     be set to a number between 0 and 1 (e.g. 0.01), which specifies the
     fraction of randomly chosen particles for which at each timestep
     forces by direct summation are computed. The normal tree-forces and
     the `correct' direct summation forces are then collected in a file \b
     forcetest.txt for later inspection. Note that the simulation itself is
     unaffected by this option, but it will of course run much(!)  slower,
     particularly if <b> FORCETEST*NumPart*NumPart>>NumPart</b>
     Note: Particle IDs must be set to numbers >=1 for this
     option to work.

\n
\section secmake9 Glass making
- \b MAKEGLASS=262144 \n This option can be used to generate a glass-like
     particle configuration. The value assigned gives the particle load,
     which is initially generated as a Poisson sample and then evolved
     towards a glass with the sign of gravity reversed

*/
