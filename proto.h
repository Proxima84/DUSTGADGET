/*! \file proto.h
 *  \brief this file contains all function prototypes of the code
 */

#ifndef ALLVARS_H
#include "allvars.h"
#endif

void   accretion(void);
void   accretion_evaluate(int target, int mode);
void   advance_and_find_timesteps(void);
void   allocate_commbuffers(void);
void   allocate_memory(void);
void   begrun(void);
int    blockpresent(enum iofields blocknr);
void   catch_abort(int sig);
void   catch_fatal(int sig);
int    compare_key(const void *a, const void *b);
void   compute_accelerations(int mode);

void   compute_potential(void);
int    dens_compare_key(const void *a, const void *b);
void   density(void);
void   density_decouple(void);
void   density_evaluate(int i, int mode);

void   distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);
double dmax(double, double);
double dmin(double, double);

void   domain_Decomposition(void); 
int    domain_compare_key(const void *a, const void *b);
int    domain_compare_key(const void *a, const void *b);
int    domain_compare_toplist(const void *a, const void *b);
void   domain_countToGo(void);
void   domain_decompose(void);
void   domain_determineTopTree(void);
void   domain_exchangeParticles(int partner, int sphflag, int send_count, int recv_count);
void   domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void   domain_findExtent(void);
int    domain_findSplit(int cpustart, int ncpu, int first, int last);
void   domain_shiftSplit(void);
void   domain_sumCost(void);
void   domain_topsplit(int node, peanokey startkey);
void   domain_topsplit_local(int node, peanokey startkey);

double drift_integ(double a, void *param);
void   dump_particles(void);
void dustGas(void);
void dustGas_evaluate(int target, int mode);

void   empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);
void   endrun(int);

void   fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);

int    find_next_outputtime(int time);
void   find_next_sync_point_and_drift(void);

void   force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount, int *nextfree);
void   force_exchange_pseudodata(void);
void   force_flag_localnodes(void);
void   force_insert_pseudo_particles(void);
void   force_setupnonrecursive(int no);
void   force_treeallocate(int maxnodes, int maxpart); 
int    force_treebuild(int npart);
int    force_treebuild_single(int npart);
int    force_treeevaluate(int target, int mode);
void   force_treefree(void);
void   force_treeupdate_pseudos(void);
void   force_update_hmax(void);
void   force_update_len(void);
void   force_update_node(int no, int flag);
void   force_update_node_hmax_local(void);
void   force_update_node_hmax_toptree(void);
void   force_update_node_len_local(void);
void   force_update_node_len_toptree(void);
void   force_update_node_recursive(int no, int sib, int father);
void   force_update_pseudoparticles(void);
void   force_update_size_of_parent_node(int no);

void   free_memory(void);

int    get_bytes_per_blockelement(enum iofields blocknr);

int    get_datatype_in_block(enum iofields blocknr);
double get_drift_factor(int time0, int time1);
double get_gravkick_factor(int time0, int time1);
double get_hydrokick_factor(int time0, int time1);
int    get_particles_in_block(enum iofields blocknr, int *typelist);
int    get_timestep(int p, double *a, int flag);
int    get_values_per_blockelement(enum iofields blocknr);

int    grav_tree_compare_key(const void *a, const void *b);

void   gravity_tree(void);
double gravkick_integ(double a, void *param);

int    hydro_compare_key(const void *a, const void *b);
void   hydro_evaluate(int target, int mode);
void   hydro_force(void);
double hydrokick_integ(double a, void *param);

int    imax(int, int);
int    imin(int, int);

void   init(void);
void   init_peano_map(void);

void   interaction();
void   interactionGas();
void   interactGas_evaluate(int target, int mode);
void   interactionDust();
void   interactDust_evaluate(int target, int mode);

void   long_range_force(void);
void   long_range_init(void);
void   long_range_init_regionsize(void);
void   move_particles(int time0, int time1);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);

int    ngb_clear_buf(FLOAT searchcenter[3], FLOAT hguess, int numngb);
void   ngb_treeallocate(int npart);
void   ngb_treebuild(void);
int    ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode, int ptype);
int    ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hguess, int *startnode, int ptype);
int    ngb_treefind_cell(FLOAT searchcenter[3], FLOAT hcell, int *startnode, int ptype);
void   ngb_treefree(void);
void   ngb_treesearch(int);
void   ngb_treesearch_pairs(int);
void   ngb_update_nodes(void);
int    ngb_treefind_rad(FLOAT searchcenter[3],int *startnode, int ptype);

void   open_outputfiles(void);

peanokey peano_hilbert_key(int x, int y, int z, int bits);
void   peano_hilbert_order(void);
double pow(double, double);  /* on some old DEC Alphas, the correct prototype for pow() is missing, even when math.h is included */

void   read_file(char *fname, int readTask, int lastTask);
void   read_ic(char *fname);

void   read_parameter_file(char *fname);
void   readjust_timebase(double TimeMax_old, double TimeMax_new);

void   remove_particles(int curtime);
void   reorder_gas(void);
void   reorder_particles(void);
void   restart(int mod);
void   run(void);
void   star_gravity(void);
void   saveaccretion(void);
void   savepositions(int num);
void   set_softenings(void);
void   set_units(void);

void   setup_smoothinglengths(void);
void   statistics(void);
void   terminate_processes(void);

void   write_file(char *fname, int readTask, int lastTask);
void   write_pid_file(void);
