/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's, typedef's, and enum's
 *     - add #include "allvars.h", delete the #ifndef ALLVARS_H conditional
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */


#define ALLVARS_H

#include <stdio.h>
#include "tags.h"

#define  GADGETVERSION   "2.0"   /*!< code version string */

#define  TIMEBASE        (1<<28) /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                  *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds to 2^29
                                  */

#define  MAXTOPNODES     200000   /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

typedef  long long  peanokey;    /*!< defines the variable type used for Peano-Hilbert keys */

#define  BITS_PER_DIMENSION 18	 /*!< Bits per dimension available for Peano-Hilbert order. 
				      Note: If peanokey is defined as type int, the allowed maximum is 10.
				      If 64-bit integers are used, the maximum is 21 */

#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))  /*!< The number of different Peano-Hilbert cells */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  MAXLEN_FILENAME  100    /*!< Maximum number of characters for filenames (including the full path) */

#ifndef   ISOTHERM_EQS 
 #ifdef   SAVETEMP_LTR
 #define  GAMMA       1.0 /*!< index for LTR gas */
 #define  Pindex      0.5 /*rate of decreasing of the midplane temperature with the distance from the center mass T=Tcoef*Tstar*(Rstar/r)^Pindex*/
 #define  Aindex      0  /*rate of decreasing of the atmosphere temperature with the distance from the disk plane T=Tstar*(Rstar/2/r)^Aindex*/
 #define  Tcoef       0.33437 /* for optical depth = 1 in optical spectra from Chiang & Goldreich (1997)*/
 #define  heightScale 3.0 /*the scale of half-thickness of the disk*/
 #else
 #define  GAMMA     1.4//(5.0/3)   /*!< adiabatic index of simulated gas */
 #endif
#else
 #define  GAMMA         1.0     /*!< index for isothermal gas */
#endif

#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.23  /*!< mass fraction of hydrogen, relevant only for radiative cooling */

#ifdef  MKSKHEMA
 #ifdef DIM1
 #define SIGMA 1
 #else
 #define SIGMA 0.333333333
 #endif
#endif

/* Some physical constants in cgs units */

#define  GRAVITY           6.672e-8   /*!< Gravitational constant (in cgs units) */
#define  SOLAR_MASS        1.989e33
#define  SOLAR_LUM         3.826e33
#define  SOLAR_RADIUS      0.004652 //(АU)
#define  AU                1.496e13
#define  RAD_CONST         7.565e-15
#define  AVOGADRO          6.0222e23
#define  BOLTZMANN         1.3806e-16
#define  GAS_CONST         8.31425e7
#define  C                 2.9979e10
#define  PLANCK            6.6262e-27
#define  CM_PER_MPC        3.085678e24
#define  PROTONMASS        1.6726e-24
#define  ELECTRONMASS      9.10953e-28
#define  THOMPSON          6.65245e-25
#define  ELECTRONCHARGE    4.8032e-10
#define  HUBBLE            3.2407789e-18	/* in h/sec */

/* Some conversion factors */
#define  SEC_PER_YEAR      3.155e7

#ifndef ASMTH
#define ASMTH 1.25  /*!< ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#endif

#ifndef RCUT
#define RCUT  4.5   /*!< RCUT gives the maximum distance (in units of the scale used for the force split) out to 
                         which short-range forces are evaluated in the short-range tree walk. */
#endif

#define MAX_NGB             20000  /*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST   500	   /*!< maxmimum number of entries in list of snapshot output times */

#define DRIFT_TABLE_LENGTH  1000   /*!< length of the lookup table used to hold the drift and kick factors */ 

#define MAXITER             150    /*!< maxmimum number of steps for SPH neighbour iteration */

            /*!< If defined, the variable type FLOAT is set to "double", otherwise to FLOAT */
#define FLOAT double

#ifdef DIM1
#define  NUMDIMS 1                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  1.333333333                 /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  8.0
#define  KERNEL_COEFF_3  24.0
#define  KERNEL_COEFF_4  16.0
#define  KERNEL_COEFF_5  2.666666667
#define  KERNEL_COEFF_6  (-8) 
#define  NORM_COEFF      1                /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */ 
#else
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470                 /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)  
#define  NORM_COEFF      4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */ 
#endif

extern int ThisTask;		/*!< the rank of the local processor */
extern int NTask;               /*!< number of processors */
extern int PTask;	        /*!< smallest integer such that NTask <= 2^PTask */

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern long long Ntype[3];      /*!< total number of particles of each type */
extern int NtypeLocal[3];       /*!< local number of particles of each type */

extern int NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
extern int NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

extern int RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                     initial conditions, 1 is resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

extern char *Exportflag;        /*!< Buffer used for flagging whether a particle needs to be exported to another process */

extern int  *Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */

extern int TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */
extern int DomainDecompositionFlag;

extern int Flag_FullStep;       /*!< This flag signals that the current step involves all particles */
extern double DomainCorner[3];    /*!< gives the lower left corner of simulation volume */
extern double DomainCenter[3];    /*!< gives the center of simulation volume */
extern double DomainLen;          /*!< gives the (maximum) side-length of simulation volume */
extern double DomainFac;          /*!< factor used for converting particle coordinates to a Peano-Hilbert mesh covering the simulation volume */
extern int    DomainMyStart;      /*!< first domain mesh cell that resides on the local processor */
extern int    DomainMyLast;       /*!< last domain mesh cell that resides on the local processor */
extern int    *DomainStartList;   /*!< a table that lists the first domain mesh cell for all processors */
extern int    *DomainEndList;     /*!< a table that lists the last domain mesh cell for all processors */
extern double *DomainWork;        /*!< a table that gives the total "work" due to the particles stored by each processor */
extern int    *DomainCount;       /*!< a table that gives the total number of particles held by each processor */
extern int    *DomainCountSph;    /*!< a table that gives the total number of SPH particles held by each processor */

extern int    *DomainTask;        /*!< this table gives for each leaf of the top-level tree the processor it was assigned to */
extern int    *DomainNodeIndex;   /*!< this table gives for each leaf of the top-level tree the corresponding node of the gravitational tree */
extern FLOAT  *DomainTreeNodeLen; /*!< this table gives for each leaf of the top-level tree the side-length of the corresponding node of the gravitational tree */
extern FLOAT  *DomainHmax;        /*!< this table gives for each leaf of the top-level tree the maximum SPH smoothing length among the particles of the corresponding node of the gravitational tree */

extern struct DomainNODE
{
  FLOAT s[3];                     /*!< center-of-mass coordinates */
  FLOAT vs[3];                    /*!< center-of-mass velocities */
  FLOAT mass;                     /*!< mass of node */
  FLOAT maxsoft;                  /*!< hold the maximum gravitational softening of particles in the 
                                       node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
}
 *DomainMoment;                   /*!< this table stores for each node of the top-level tree corresponding node data from the gravitational tree */

extern peanokey *DomainKeyBuf;    /*!< this points to a buffer used during the exchange of particle data */

extern peanokey *Key;             /*!< a table used for storing Peano-Hilbert keys for particles */
extern peanokey *KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for all particles, used to construct top-level tree */


extern int NTopnodes;             /*!< total number of nodes in top-level tree */
extern int NTopleaves;            /*!< number of leaves in top-level tree. Each leaf can be assigned to a different processor */

extern struct topnode_data
{
  int Daughter;                   /*!< index of first daughter cell (out of 8) of top-level node */
  int Pstart;                     /*!< for the present top-level node, this gives the index of the first node in the concatenated list of topnodes collected from all processors */
  int Blocks;                     /*!< for the present top-level node, this gives the number of corresponding nodes in the concatenated list of topnodes collected from all processors */
  int Leaf;                       /*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  peanokey Size;                  /*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey;              /*!< first Peano-Hilbert key in top-level node */
  long long Count;                /*!< counts the number of particles in this top-level node */
}
 *TopNodes;                       /*!< points to the root node of the top-level tree */

extern double TimeOfLastTreeConstruction; /*!< holds what it says, only used in connection with FORCETEST */

/* variables for input/output, usually only used on process 0 */

extern char ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for starting the simulation */
extern void *CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                                parallel tree-force computation, the SPH routines, etc. */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!< total particle numbers (global value) */
  long long TotN_gas;		/*!< total gas particle number (global value) */

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */
  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density computation */
  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH hydrodynamical force computation */
  int BunchSizeDrag;           /*!< number of particles fitting into the communication buffer in the SPH drag force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain decomposition */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				     NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				     the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				     the maximum(!) number of particles.  Note: A typical local tree for N
				     particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  double DesNumNgb;             /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */


  /* some force counters  */
  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */


  /* system of units  */

  double G;                        /*!< Gravity-constant in internal units */
  double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h */
  double UnitMass_in_g;            /*!< factor to convert internal mass unit to grams/h */
  double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to cm/sec */
  double UnitLength_in_cm;         /*!< factor to convert internal length unit to cm/h */
  double UnitPressure_in_cgs;      /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
  double UnitDensity_in_cgs;       /*!< factor to convert internal length unit to g/cm^3*h^2 */
  double UnitCoolingRate_in_cgs;   /*!< factor to convert internal cooling rate to cgs units */
  double UnitEnergy_in_cgs;        /*!< factor to convert internal energy to cgs units */
    /*!< factor to convert internal time to megayears/h */
  double GravityConstantInternal;  /*!< If set to zero in the parameterfile, the internal value of the
                                        gravitational constant is set to the Newtonian value based on the system of
                                        units specified. Otherwise the value provided is taken as internal gravity constant G. */

  /* Code options */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */

  /* Parameters determining output frequency */

  int SnapshotFileCount;        /*!< number of snapshot that is written next */
  double TimeBetSnapshot;       /*!< simulation time interval between snapshot files */
  double TimeOfFirstSnapshot;   /*!< simulation time of first snapshot files */ 
  int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this point */


  /* Current time of the simulation, global step, and end of simulation */

  double Time;                  /*!< current time of the simulation */
  double TimeBegin;             /*!< time of initial conditions of the simulation */
  double TimeStep;              /*!< difference between current times of previous and current timestep */
  double TimeMax;	        /*!< marks the point of time until the simulation is to be evolved */


  /* variables for organizing discrete timeline */

  double Timebase_interval;     /*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;               /*!< current time on integer timeline */ 
  int Ti_nextoutput;            /*!< next output time on integer timeline */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                     timestep is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                     timestep determined by the timestep criteria falls below this limit. */ 
  double MaxSizeTimestep;       /*!< maximum allowed timestep */

  double CourantFac;		/*!< SPH-Courant factor */
  /* frequency of tree reconstruction/domain decomposition */

  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */

  /* Gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length).
   * Five groups of particles are supported 0="gas", 1="dust", 2="stars"
   */ 
  double MinGasHsml;            /*!< minimum allowed SPH smoothing length */
  FLOAT HsmlConstant;          /*!< constant smoothing length if it is determined*/

  double DustGasMechStep;
  double DustSize;
  double DustRho;
  double ForceSoftening[3];     /*!< current (comoving) gravitational softening lengths for each particle type */
  double MassTable[3];          /*!< Table with particle masses for particle types with equal mass.
                                     If particle masses are all equal for one type, the corresponding entry in MassTable 
                                     is set to this value, allowing the size of the snapshot files to be reduced. */  
  double BarrierDistance;          /*!< radius of the calculated area */  
  double AccretionRadius;          /*!< radius of the accretion in stellar radius */  
  double StarMassIndicator;        /*!< minimum mass of the star to exclude disk selfgravity*/
  double CentralTemperature;      /*!< temperature in the centre of mass*/
  double CentralRadius;          /*!< radius for the disk temperature*/
#ifdef   SAVETEMP_LTR
  double CentralMass;          /*!< mass of the center of the mass*/
#endif
  /* some filenames */

  char InitCondFile[MAXLEN_FILENAME];          /*!< filename of initial conditions */
  char OutputDir[MAXLEN_FILENAME];             /*!< output directory of the code */
  char SnapshotFileBase[MAXLEN_FILENAME];      /*!< basename to construct the names of snapshotf files */
  char RestartFile[MAXLEN_FILENAME];           /*!< basename of restart-files */
}
 All;                                          /*!< a container variable for global variables that are equal on all processors */



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  FLOAT Pos[3];			/*!< particle position at its current time */
  FLOAT Mass;			/*!< particle mass */
  FLOAT Vel[3];			/*!< particle velocity at its current time */
  FLOAT GravAccel[3];		/*!< particle acceleration due to gravity */
  FLOAT OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
  int ID;		/*!< particle identifier */
  float GravCost;		/*!< weight factor used for balancing the work-load */
  int Type;		        /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  int Ti_endstep;               /*!< marks start of current timestep of particle on integer timeline */ 
  int Ti_begstep;               /*!< marks end of current timestep of particle on integer timeline */
  int Accretion;
  FLOAT Radius;
#ifdef NOSTARSELFGRAV
  int ItsStar;	
#endif
}
 *P,              /*!< holds particle data on local processor */
 *DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  FLOAT Entropy;                /*!< current value of entropy (actually entropic function) of particle */
  FLOAT Density;		/*!< current baryonic mass density of particle */
  FLOAT Hsml;			/*!< current smoothing length */
  FLOAT Left;                   /*!< lower bound in iterative smoothing length search */  
  FLOAT Right;                  /*!< upper bound in iterative smoothing length search */ 
  FLOAT NumNgb;                 /*!< weighted number of neighbours found */
  FLOAT Pressure;		/*!< current pressure */
  FLOAT DtEntropy;              /*!< rate of change of entropy */
  FLOAT DtDragEntropy;          /*!< rate of change of drag entropy */
  FLOAT HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
  FLOAT DragAccel[3];             /*!< acceleration due to drag force */
  FLOAT VelPred[3];		/*!< predicted SPH particle velocity at the current time */
  FLOAT DivVel;			/*!< local velocity divergence */
  FLOAT CurlVel;		/*!< local velocity curl */
  FLOAT Rot[3];		        /*!< local velocity curl */
#ifndef HSMLCONSTANT
  FLOAT DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative entropy formulation of SPH */
#endif
  FLOAT MaxSignalVel;           /*!< maximum "signal velocity" occuring for this particle */

  FLOAT GasVelMid[3];	        /*Middle gas velocity in cell for calculation of drag force*/
  FLOAT DustVelMid[3];          /*Middle dust velocity in cell for calculation of drag force*/
  FLOAT HydroAccelMid[3];       /*Middle hydro acceleration in cell for calculation of drag force*/
  FLOAT GasDensity;             /*Middle gas density in cell for calculation of drag force*/
  FLOAT DustDensity;            /*Middle dust density in cell for calculation of drag force*/
  FLOAT SoundSpeedMid;          /*Middle sound speed in cell for calculation of drag force*/
  int DustNumber;                /*Middle dust number in cell for calculation of drag force*/
  int GasNumber;                  /*Middle gas number in cell for calculation of drag force*/
#ifdef HSMLCONSTANT
  FLOAT   DustGasNgb;
#endif
}
 *SphP,                        	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;                 /*!< buffer for SPH particle data in domain decomposition */

/*  Variables for Tree
 */

extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */

extern struct NODE
{
  FLOAT len;			/*!< sidelength of treenode */
  FLOAT center[3];		/*!< geometrical center of node */
  FLOAT maxsoft;                /*!< hold the maximum gravitational softening of particles in the 
                                     node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      FLOAT s[3];               /*!< center of mass of node */
      FLOAT mass;               /*!< mass of node */
      int bitflags;             /*!< a bit-field with various information on the node */
      int sibling;              /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;             /*!< this gives the next node in case the current node needs to be opened */
      int father;               /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
 *Nodes_base,                   /*!< points to the actual memory allocted for the nodes */
 *Nodes;                        /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
 				     gives the first allocated node */
extern int *Nextnode;	        /*!< gives next node in tree walk */
extern int *Father;	        /*!< gives parent node in tree    */


extern struct extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
{
  FLOAT hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  FLOAT vs[3];			/*!< center-of-mass velocity */
}
 *Extnodes_base,                /*!< points to the actual memory allocted for the extended node information */
 *Extnodes;                     /*!< provides shifted access to extended node information, parallel to Nodes/Nodes_base */

/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[3];                        /*!< number of particles of each type in this file */
  double mass[3];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
//  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[3];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  //unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[184];	               /*!< fills to 256 Bytes */
}
 header;                               /*!< holds header for snapshot files */


#define IO_NBLOCKS 11   /*!< total number of defined information blocks for snapshot files.
                             Must be equal to the number of entries in "enum iofields" */

enum iofields           /*!< this enumeration lists the defined output blocks in snapshot files. Not all of them need to be present. */
{ 
  IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_RHO,
  IO_HSML,
  IO_POT,
  IO_ACCEL,
  IO_DTENTR,
  IO_TSTP,
};

/* Various structures for communication
 */
extern struct gravdata_in
{
  union
  {
    FLOAT Pos[3];
    FLOAT Acc[3];
    FLOAT Potential;
  }
  u;
  int Type;
  int id;
  FLOAT Mass;
  FLOAT Soft;
  union
  {
    FLOAT OldAcc;
    int Ninteractions;
  }
  w;
}
 *GravDataIn,                   /*!< holds particle data to be exported to other processors */
 *GravDataGet,                  /*!< holds particle data imported from other processors */
 *GravDataResult,               /*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;                  /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

extern struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 *GravDataIndexTable;           /*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */

extern struct densdata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  int Index;
  int Task;
  int Type;
  //  int id;
}
 *DensDataIn,                   /*!< holds particle data for SPH density computation to be exported to other processors */
 *DensDataGet;                  /*!< holds imported particle data for SPH density computation */

extern struct densdata_out
{
  FLOAT Rho;
  FLOAT Div, Rot[3];
#ifndef HSMLCONSTANT
  FLOAT DhsmlDensity;
#endif
  FLOAT Ngb;
  
}
 *DensDataResult,               /*!< stores the locally computed SPH density results for imported particles */
 *DensDataPartialResult;        /*!< imported partial SPH density results from other processors */


extern struct hydrodata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT Density;
  FLOAT Pressure;
  FLOAT F1;
#ifndef HSMLCONSTANT
  FLOAT DhsmlDensityFactor;
#endif
  int Timestep;
  int   Task;
  int   Index;
 // int id;

}
 *HydroDataIn,                  /*!< holds particle data for SPH hydro-force computation to be exported to other processors */
 *HydroDataGet;                 /*!< holds imported particle data for SPH hydro-force computation */

extern struct hydrodata_out
{
  FLOAT Acc[3];
  FLOAT DtEntropy;
  FLOAT MaxSignalVel;
  FLOAT AccDrag[3];  

}
 *HydroDataResult,              /*!< stores the locally computed SPH hydro results for imported particles */
 *HydroDataPartialResult;       /*!< imported partial SPH hydro-force results from other processors */

extern struct dragdata_in
{
  FLOAT Pos[3];
#ifdef HSMLCONSTANT
 int   Type;
#endif
  int   Timestep;
  int   Task;
  int   Index; 
}
 *DragDataIn,                  /*!< holds particle data for SPH drag-force computation to be exported to other processors */
 *DragDataGet;                 /*!< holds imported particle data for SPH drag-force computation */

extern struct dragdata_out
{
  FLOAT GasVelMid[3];
  FLOAT DustVelMid[3];
  FLOAT HydroAccelMid[3];
  FLOAT GasDensity;
  FLOAT DustDensity;
  FLOAT SoundSpeedMid;
  int DustNumber;
  int GasNumber;
#ifdef HSMLCONSTANT
  FLOAT   DustGasNgb;
#endif
 
}
 *DragDataResult,              /*!< stores the locally computed SPH drag results for imported particles */
 *DragDataPartialResult;       /*!< imported partial SPH drag-force results from other processors */


