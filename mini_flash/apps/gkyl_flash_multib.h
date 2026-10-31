#pragma once

#include <gkyl_flash_block_geom.h>
#include <gkyl_flash.h>

typedef struct gkyl_flash_multib_app gkyl_flash_multib_app;

// Simulation statistics for the multi-block app, reduced over all ranks. 
// The gkyl_flash_stat members store the maximum and average over all ranks.
// Note `init_tm` includes communication time because of an internal barrier.
struct gkyl_flash_multib_stat {
  int nup; // Number of update calls (identical on all ranks).
  int nfeuler; // Number of forward Euler steps taken (identical on all ranks).
  int nfail; // Number of failed SSP-RK stages due to a smaller dt (identical on all ranks).

  struct gkyl_flash_stat max; // Maximum over ranks of each statistic.
  struct gkyl_flash_stat avg; // Average over ranks of each statistic.
};

// Domain decomposition scheme. Either:
// 1. Round-robin assignment of ranks to blocks
// 2. Row-major ordered distribution
enum gkyl_flash_multib_decomp_type {
  GKYL_FLASH_MULTIB_DECOMP_RR = 0,
  GKYL_FLASH_MULTIB_DECOMP_LINEAR
};

struct gkyl_flash_multib_species_inp_per_block {
  int block_id; // *GLOBAL* block ID.

  void *ctx; // Context struct for initial conditions.
  // Distribution function initialization function.
  void (*init_dist)(double t, const double *xn, double *fout, void *ctx);
  // Fluid quantities initialization function.
  void (*init_fluid)(double t, const double *xn, double *fout, void *ctx);

  bool init_dist_from_arrays; // Flag to initialize the distribution function as a Maxwellian from moment
                              // arrays [n, T/m] instead of the init_dist function. The moments are either
                              // given in init_dist_moms, or set after app creation with the interface
                              // functions (gkyl_flash_multib_app_init_maxwellian_from_moms).
  struct gkyl_array *init_dist_moms; // Optional DG moments [n, T/m] on the extended kinetic configuration-space
                                     // range of this block (2*num_conf_basis components). Copied at app creation and
                                     // projected onto a Maxwellian when initial conditions are applied.

  bool skip_init_dist; // Flag to skip initializing the distribution function.
};

struct gkyl_flash_multib_field_inp_per_block {
  int block_id; // *GLOBAL* block ID.

  void *ctx; // Context structure for initial condition init function.
  void (*init)(double t, const double *xn, double *fout, void *ctx); // Initialization function pointer

  void *ext_em_ctx; // Context structure for external field initial condition init function.
  void (*ext_em_init)(double t, const double *xn, double *fout, void *ctx); // Initialization function pointer for external field

  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2]; // Boundary conditions.
};

// BC for blocks
struct gkyl_flash_block_physical_bcs {
  int block_id; // *GLOBAL* block ID.
  int dir;  // Direction in which BC is specified.
  enum gkyl_edge_loc edge; // Edge for this BC.
  enum gkyl_species_bc_type bc_type; // BC type flag.
};

struct gkyl_flash_multib_species {
  char name_kin[128], name_flu[128]; // Species name.

  double mass; // Species mass.
  int cells[3]; // Number of velocity space cells for Vlasov species.
  double lower[3], upper[3]; // Velocity space bounds.

  double gas_gamma; // Adiabatic index of the Euler fluid (defaults to 5/3).
  enum gkyl_wave_limiter limiter; // Fluid slope limiter (defaults to monotonized centered).

  int num_physical_bcs; // Number of domain boundary conditions.
  struct gkyl_flash_block_physical_bcs *bcs; // Physical boundary conditions.

  bool write_pkpm_diags; // Optional flag to write PKPM diagnostic moments:
                         // (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp).

  bool duplicate_ics; // Flag for duplicating initial conditions across all blocks.
  struct gkyl_flash_multib_species_inp_per_block *blocks; // Per-block initial conditions.

  struct gkyl_flash_collisions collisions; // Collisions input.
};

struct gkyl_flash_multib_field {
  bool is_static;
  double epsilon0, mu0;

  double limiter_fac; // Limiter factor for adjusting diffusion in slope limiter.

  bool duplicate_ics; // Flag for duplicating initial conditions across all blocks.
  struct gkyl_flash_multib_field_inp_per_block *blocks; // Per-block initial conditions.
};


struct gkyl_flash_multib_inp {
  char name[128]; // App name.
  int cdim, vdim; // Configuration and velocity space dimensionality.
  int num_blocks_global; // Total number of blocks.

  int cells_flu[GKYL_MAX_CDIM]; // Cells per block.
  double lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM]; // Global configuration space bounds.
  int nghost_kin[GKYL_MAX_CDIM], nghost_flu[GKYL_MAX_CDIM]; // Number of ghost cells for kinetic/fluid vars (defaults to 1/2).
  struct gkyl_flash_block_geom *block_geom; // Block geometry object.

  double cfl_frac; // CFL scaling factor (default 1.0).

  enum gkyl_flash_multib_decomp_type decomp_type; // Flag for method of decomposing the global grid into blocks and 
                                                  // assigning ranks per block.

  // Periodic BCs
  int num_periodic_dir;
  int periodic_dirs[GKYL_MAX_CDIM];

  int num_species;
  struct gkyl_flash_multib_species species[GKYL_MAX_SPECIES]; // Species input structures.

  struct gkyl_flash_multib_field field; // Field input struct.

  bool skip_fluid_species; // Flag to not initialize fluid species.
  
  struct gkyl_comm *comm; // Communicator to use.
};

/**
 * Initialize a new multi-block app on the current rank.
 * 
 * @param mbinp 
 * @return gkyl_flash_multib_app* 
 */
gkyl_flash_multib_app*
gkyl_flash_multib_app_new(struct gkyl_flash_multib_inp *mbinp);

/**
 * Apply initial conditions for each block.
 * 
 * @param app 
 * @param t0 
 */
void
gkyl_flash_multib_app_apply_ic(gkyl_flash_multib_app *app, double t0);

/**
 * Apply initial conditions to species sidx on all local blocks.
 * Must be called on all ranks.
 *
 * @param app Multiblock app object.
 * @param sidx Species index.
 * @param t0 Time for initial conditions.
 */
void
gkyl_flash_multib_app_apply_ic_species(gkyl_flash_multib_app *app, int sidx, double t0);

/**
 * Apply field initial conditions on all local blocks.
 * Must be called on all ranks.
 *
 * @param app Multiblock app object.
 * @param t0 Time for initial conditions.
 */
void
gkyl_flash_multib_app_apply_ic_field(gkyl_flash_multib_app *app, double t0);

/**
 * Initialize the distribution functions on all local blocks as Maxwellians
 * matching DG moments [n, T/m] on each block's extended kinetic
 * configuration space range, then apply physical and inter-block boundary
 * conditions. Input init_dist_moms[b] holds the moments for local block b.
 * Only needed if the moments were not given in the species input at app creation.
 * Must be called once on all ranks.
 *
 * @param app Multiblock app object.
 * @param init_dist_moms DG moment arrays, one per local block (copied).
 */
void
gkyl_flash_multib_app_init_maxwellian_from_moms(gkyl_flash_multib_app *app,
  struct gkyl_array *init_dist_moms[]);

/**
 * Apply boundary conditions to distribution function and fluid vars for each block.
 * 
 * @param app 
 * @param tcurr 
 * @param distf Pointer to array of distribution functions for each species.
 * @param fluid 
 */
void
gkyl_flash_multib_app_apply_bc(gkyl_flash_multib_app *app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *fluid[]);

/**
 * Compute the maximum stable timestep over all blocks on all ranks.
 *
 * @param app
 * @return double
 */
double
gkyl_flash_multib_app_compute_dt_max(gkyl_flash_multib_app *app);

/**
 * Advance the hybrid system on every block by a single SSP-RK3 update.
 * The actual dt may be smaller than the suggested timestep dt0 but never larger.
 * The returned status struct contains the actual and suggested timesteps, which
 * are the minimum over all blocks on all ranks.
 *
 * @param app
 * @param dt0 Suggested timestep
 * @return struct gkyl_update_status
 */
struct gkyl_update_status
gkyl_flash_multib_app_update(gkyl_flash_multib_app *app, double dt0);

/**
 * Take a single forward Euler step of ONLY the kinetic equation on every block,
 * using the global stable timestep (never larger than dt0). 
 * The returned status struct contains the actual and suggested timesteps, which
 * are the minimum over all blocks on all ranks.
 *
 * @param app
 * @param dt0 Suggested timestep
 * @return struct gkyl_update_status
 */
struct gkyl_update_status
gkyl_flash_multib_app_kinetic_update(gkyl_flash_multib_app *app, double dt0);

/**
 * Write out diagnostic files for each block of the multib app.
 * 
 * @param app 
 * @param tm 
 * @param frame 
 */
void
gkyl_flash_multib_app_write(struct gkyl_flash_multib_app *app, double tm, int frame);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by only writing from rank 0 of the global communicator.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void
gkyl_flash_multib_app_cout(const gkyl_flash_multib_app *app, FILE *fp, const char *fmt, ...);

/**
 * Number of blocks handled by this rank.
 *
 * @param app Multiblock app object.
 * @return Number of local blocks.
 */
int
gkyl_flash_multib_app_num_local_blocks(const gkyl_flash_multib_app *app);

/**
 * Global block ID of local block b, for 0 <= b < number of local blocks.
 *
 * @param app Multiblock app object.
 * @param b Local block index.
 * @return Global block ID.
 */
int
gkyl_flash_multib_app_local_block_id(const gkyl_flash_multib_app *app, int b);

/**
 * Return the extended configuration space range over the fluid grid of
 * local block blockID.
 *
 * @param app Multiblock app object.
 * @param blockID Local block index.
 * @return Extended local fluid range.
 */
struct gkyl_range
gkyl_flash_multib_app_get_fluid_conf_ext_range(gkyl_flash_multib_app *app, int blockID);

/**
 * Return the configuration space range over the fluid grid of local block
 * blockID.
 *
 * @param app Multiblock app object.
 * @param blockID Local block index.
 * @return Local fluid range.
 */
struct gkyl_range
gkyl_flash_multib_app_get_fluid_conf_range(gkyl_flash_multib_app *app, int blockID);

/**
 * Compute simulation statistics reduced over all ranks: the maximum and
 * average of each timing, where a rank's timings are summed over its local
 * blocks. Must be called on all ranks.
 *
 * @param app Multiblock app object.
 * @return Reduced statistics.
 */
struct gkyl_flash_multib_stat
gkyl_flash_multib_app_stat(gkyl_flash_multib_app *app);

/**
 * Print timing statistics for multiblock app. Includes same statistics as
 * single-block apps, both the maxima and averages over all ranks.
 * 
 * @param app 
 * @param fp 
 */
void
gkyl_flash_multib_app_print_timings(gkyl_flash_multib_app *app, FILE *fp);

/**
 * Release memory for multiblock app object.
 * 
 * @param mbapp App to free
 */
void
gkyl_flash_multib_app_release(gkyl_flash_multib_app *mbapp);
