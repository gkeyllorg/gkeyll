#pragma once

#include <gkyl_flash_block_geom.h>
#include <gkyl_flash.h>

typedef struct gkyl_flash_multib_app gkyl_flash_multib_app;

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

  bool init_dist_from_arrays; // Flag to initialize distribution function from input fluid moment arrays.
  struct gkyl_array *init_dist_moms; // Input fluid moment arrays for initializing the distibution function with
                                     // the `init_dist_from_arrays` flag. Moments should be ordered as:
                                     // [rho, ux, uy, uz, temp].

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
 * Compute the maximum stable timestep over all blocks on all ranks.
 *
 * @param app
 * @return double
 */
double
gkyl_flash_multib_app_compute_dt_max(gkyl_flash_multib_app *app);

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
gkyl_flash_multib_kinetic_update(gkyl_flash_multib_app *app, double dt0);

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
gkyl_flash_multib_update(gkyl_flash_multib_app *app, double dt0);

/**
 * Release memory for multiblock app object.
 * 
 * @param mbapp App to free
 */
void
gkyl_flash_multib_app_release(gkyl_flash_multib_app *mbapp);