#pragma once

#include <gkyl_app.h>
#include <gkyl_wave_prop.h>

// Object type.
typedef struct gkyl_flash_app gkyl_flash_app;

// Simulation statistics for mini-FLASH hybrid simulation.
struct gkyl_flash_stat {
  int nup; // Number of update calls.
  int nfeuler; // Number of forward Euler steps taken.
  int nfail; // Number of failed SSP-RK stages due to a smaller dt.

  double total_tm; // Time spent in update calls (excludes initial conditions and I/O).
  double fwd_euler_tm; // Time spent in forward Euler steps.

  double init_tm; // Time to apply initial and boundary conditions to all species and the field.
  double init_species_tm; // Time to apply initial conditions to kinetic species.
  double init_fluid_species_tm; // Time to apply initial conditions to fluid species.
  double init_field_tm; // Time to apply initial conditions to the field.

  double field_vars_tm; // Time to compute b_hat, div(b) and max(b).
  double species_mom_tm; // Time to compute PKPM moments and variables for the kinetic update.
  double species_coll_mom_tm; // Time to compute collision operator moments.
  double coupling_tm; // Time to interpolate between kinetic and fluid grids and compute coupling fluxes.

  double species_rhs_tm; // Time to compute collisionless kinetic RHS.
  double species_coll_tm; // Time to compute collision operator RHS.
  double fluid_rhs_tm; // Time to compute fluid update (wave propagation and kinetic closure terms).

  double species_bc_tm; // Time to apply kinetic species boundary conditions.
  double fluid_species_bc_tm; // Time to apply fluid species boundary conditions.

  double io_tm; // Time to compute and write output frames (includes diagnostic moments).
};

// Input object specifying details for collisions in mini-FLASH updater.
// User specifies colliding species and collision frequency profiles
// for each colliding pair.
struct gkyl_flash_collisions {
  void *ctx;
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);

  int num_cross_collisions;
  char collide_with[GKYL_MAX_SPECIES][128];
};

// Input parameters for species in mini-FLASH hybrid updater.
// User specifies velocity space grid information and initial conditions 
// for the PKPM distribution functions (F0, G) and 
// Euler fluid variables (rho, rho u_x, rho u_y, rho u_z, energy_dens).
struct gkyl_flash_species {
  char name_kinetic[128], name_fluid[128]; // Species name.
  double mass; // Species mass.

  int cells[3]; // Number of velocity space cells for Vlasov species.
  double lower[3], upper[3]; // Velocity space bounds.

  void *ctx; // Context struct for initial conditions.
  // Distribution function initialization function.
  void (*init_dist)(double t, const double *xn, double *fout, void *ctx);
  // Fluid quantities initialization function.
  void (*init_fluid)(double t, const double *xn, double *fout, void *ctx);

  double gas_gamma;

  enum gkyl_wave_limiter limiter;
  double limiter_fac; // Optional input parameter for adjusting diffusion in slope limiter.
  bool limit_fluid; // Optional input parameter for applying limiters to fluid variables.

  struct gkyl_flash_collisions collisions;

  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];

  bool write_pkpm_diags; // Optional flag to write PKPM diagnostic moments:
                         // (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp).

  bool init_dist_from_arrays; // Flag to initialize distribution function from input fluid moment arrays.
  struct gkyl_array *init_dist_moms; // Input fluid moment arrays for initializing the distibution function with
                                     // the `init_dist_from_arrays` flag. Moments should be ordered as:
                                     // [rho, ux, uy, uz, temp].

  bool skip_init_dist; // Flag to skip initializing the distribution function.
};

// Input parameters for electromagnetic fields in mini-FLASH hybrid updater.
// User specifies field initial conditions and any static external fields.
struct gkyl_flash_field {
  bool is_static;
  double epsilon0, mu0;

  double limiter_fac; // Limiter factor for adjusting diffusion in slope limiter.

  void *ctx; // Context structure for initial condition init function.
  void (*init)(double t, const double *xn, double *fout, void *ctx); // Initialization function pointer

  void *ext_em_ctx; // Context structure for external field initial condition init function.
  void (*ext_em_init)(double t, const double *xn, double *fout, void *ctx); // Initialization function pointer for external field

  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2]; // Boundary conditions.
};

// Top-level input parameters for `gkyl_flash_app`.
// User specifies configuration space grid for finite volume fluid quantities.
// The kinetic grid is constructed from the fluid grid with half the resolution 
// in configuration space.
struct gkyl_flash_inp {
  char name[128];
  int cdim, vdim; // Configuration and velocity space dimensionality.
  double lower_conf[GKYL_MAX_CDIM], upper_conf[GKYL_MAX_CDIM]; // Configuration space bounds.
  int cells_fluid[GKYL_MAX_CDIM]; // Number of configuration space cells for finite volume fluid species.
  int nghost_kinetic[GKYL_MAX_CDIM], nghost_fluid[GKYL_MAX_CDIM]; // Number of ghost cells for kinetic/fluid vars (defaults to 1/2).
  double cfl_frac;

  int num_species;
  struct gkyl_flash_species species[GKYL_MAX_SPECIES]; // Species input object.

  struct gkyl_flash_field field; // Field input object.

  // Axes with periodic boundary conditions.
  int num_periodic_dir;
  int periodic_dirs[3];

  bool skip_fluid_species; // Flag to not initialize fluid species.

  struct gkyl_app_parallelism_inp parallelism; // Parallelism-related inputs.
};

/**
 * Initialize a new updater for coupled Euler-flash system.
 * 
 * @param inp 
 * @return gkyl_flash_app* 
 */
gkyl_flash_app*
gkyl_flash_app_new(struct gkyl_flash_inp *inp);

/**
 * Compute the maximum timestep for the kinetic update without actually 
 * computing the RHS of the kinetic equations.
 * 
 * @param app 
 * @return double 
 */
double
gkyl_flash_app_compute_dt_max(gkyl_flash_app *app);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void 
gkyl_flash_app_apply_ic(gkyl_flash_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_vm object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void 
gkyl_flash_app_apply_ic_species(gkyl_flash_app* app, int sidx, double t0);

/**
 * Initialize field by projecting initial conditions on basis functions.
 * 
 * @param app 
 * @param t0 
*/
void 
gkyl_flash_app_apply_ic_field(gkyl_flash_app* app, double t0);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_flash_app_write(gkyl_flash_app* app, double tm, int frame);

/**
 * Write species data to file: distribution function and fluid moments.
 * 
 * @param app 
 * @param sidx
 * @param tm 
 * @param frame 
 */
void gkyl_flash_app_write_species(gkyl_flash_app* app, int sidx, double tm, int frame);

void
gkyl_flash_app_write_fluid_species(gkyl_flash_app *app, int sidx, double tm, int frame);

/**
 * Compute and write PKPM diagnostic moments to file:
 * (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
 * 
 * @param app 
 * @param sidx
 * @param tm 
 * @param frame 
 */
void 
gkyl_flash_app_write_pkpm_moms(gkyl_flash_app* app, int sidx, double tm, int frame);

/**
 * Write electromagnetic field data to file.
 * 
 * @param app 
 * @param tm 
 * @param frame 
 */
void
gkyl_flash_app_write_field(gkyl_flash_app *app, double tm, int frame);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by only writing from rank 0 of the kinetic communicator.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void
gkyl_flash_app_cout(const gkyl_flash_app *app, FILE *fp, const char *fmt, ...);

/**
 * Advance the mini-FLASH hybrid system by a single timestep with suggested
 * interval dt. The actual dt used in the update may be smaller than the 
 * suggested one but never larger, and is returned in the `gkyl_update_status`
 * struct.  
 * 
 * @param app App object.
 * @param dt Suggested timestep.
 * @return struct gkyl_update_status 
 */
struct gkyl_update_status 
gkyl_flash_update(gkyl_flash_app* app, double dt);


/**
 * Get struct containing simulations statistics.
 * 
 * @param app 
 * @return struct gkyl_flash_stat 
 */
struct gkyl_flash_stat
gkyl_flash_app_stat(gkyl_flash_app *app);

/**
 * Print timing statistics to console (from rank 0 only).
 *
 * @param app App object.
 * @param fp File pointer for output.
 */
void
gkyl_flash_app_print_timings(gkyl_flash_app *app, FILE *fp);

/**
 * Release memory for app and associated pointer members.
 * 
 * @param app App object.
 */
void
gkyl_flash_app_release(gkyl_flash_app *app);