#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_vlasov_comms.h>
#include <gkyl_mom_type.h>

#include <stdbool.h>

// Parameters for projection
struct gkyl_vlasov_projection {
  enum gkyl_projection_id proj_id; // type of projection (see gkyl_eqn_type.h)
  enum gkyl_quad_type quad_type; // quadrature scheme to use: defaults to Gaussian
  
  union {
    struct {
      // pointer and context to initialization function 
      void *ctx_func; 
      void (*func)(double t, const double *xn, double *fout, void *ctx); 
    };
    struct {
      // pointers and contexts to initialization functions for LTE distribution projection
      // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
      void *ctx_density;
      void (*density)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_V_drift;
      void (*V_drift)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_temp;
      void (*temp)(double t, const double *xn, double *fout, void *ctx);

      // boolean if we are correcting all the moments or only density
      bool correct_all_moms;  
      double iter_eps; // error tolerance for moment fixes (density is always exact)
      int max_iter; // maximum number of iteration
      bool use_last_converged; // use last iteration value regardless of convergence?
    };
  };
};

// Parameters for species collisions
struct gkyl_vlasov_collisions {
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)
  bool write_diagnostics; // Whether to write diagnostics out.

  void *ctx; // context for collision function
  // function for computing self-collision frequency
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);

  // inputs for Spitzer collisionality
  bool normNu; // Set to true if you want to rescale collision frequency
  double nuFrac; // Parameter for rescaling collision frequency from SI values
  double hbar; // Planck's constant/2 pi 

  // BGK collisions specific inputs
  bool correct_all_moms; // boolean if we are correcting all the moments or only density
  double iter_eps; // error tolerance for moment fixes (density is always exact)
  int max_iter; // maximum number of iteration
  bool fixed_temp_relax; // Are BGK collisions relaxing to a fixed input temperature?
  bool use_last_converged; // use last iteration value regardless of convergence?

  // Boolean for using implicit BGK collisions (replaces rk3)   
  bool has_implicit_coll_scheme; 

  int num_cross_collisions; // number of species to cross-collide with
  char collide_with[GKYL_MAX_SPECIES][128]; // names of species to cross collide with

  char collide_with_fluid[128]; // name of fluid species to cross collide with
};

// Parameters for species radiation
struct gkyl_vlasov_radiation {
  enum gkyl_radiation_id radiation_id; // type of radiation (see gkyl_eqn_type.h)

  void *ctx_nu; // context for collision frequency
  // function for computing collision frequency
  void (*nu)(double t, const double *xn, double *fout, void *ctx);

  void *ctx_nu_rad_drag; // context for collision frequency * drag
  // function for computing collision frequency * drag
  void (*nu_rad_drag)(double t, const double *xn, double *fout, void *ctx);  
};

// Parameters for species source
struct gkyl_vlasov_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // optional parameter to write out source
  int num_sources;

  double source_length; // required for boundary flux source
  char source_species[128];
  
  // sources using projection routine
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];
};

// Parameters for boundary conditions
struct gkyl_vlasov_bc {
  enum gkyl_species_bc_type type;
  void *aux_ctx;
  void (*aux_profile)(double t, const double *xn, double *fout, void *ctx);  
  double aux_parameter;
};

struct gkyl_vlasov_bcs {
  struct gkyl_vlasov_bc lower, upper;
};

// Parameters for fluid species source
struct gkyl_vlasov_fluid_source {
  enum gkyl_source_id source_id; // type of source

  void *ctx; // context for source function
  // function for computing source profile
  void (*profile)(double t, const double *xn, double *aout, void *ctx);
};

// Parameters for fluid species advection
struct gkyl_vlasov_fluid_advection {
  void *velocity_ctx; // context for applied advection function
  // pointer to applied advection velocity function
  void (*velocity)(double t, const double *xn, double *aout, void *ctx);
};

// Parameters for fluid species diffusion
struct gkyl_vlasov_fluid_diffusion {
  double D; // constant diffusion coefficient
  int order; // integer for order of the diffusion (4 for grad^4, 6 for grad^6, default is grad^2)
  void* Dij_ctx; // context for applied diffusion function if using general diffusion tensor
  // pointer to applied diffusion function is using general diffusion tensor 
  void (*Dij)(double t, const double* xn, double* Dout, void* ctx);
};

// Parameters for Vlasov species
struct gkyl_vlasov_species {
  char name[128]; // species name

  enum gkyl_model_id model_id; // type of model 
                               // (e.g., SR, general geometry, see gkyl_eqn_type.h)

  double charge, mass; // charge and mass
  double lower[3], upper[3]; // lower, upper bounds of velocity-space
  int cells[3]; // velocity-space cells

  // initial conditions using projection routine
  int num_init; // number of initial condition functions
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];

  int num_diag_moments; // number of diagnostic moments
  enum gkyl_distribution_moments diag_moments[16]; // list of diagnostic moments

  // collisions to include
  struct gkyl_vlasov_collisions collisions;

  // radiation to include
  struct gkyl_vlasov_radiation radiation;

  // source to include
  struct gkyl_vlasov_source source;

  void *app_accel_ctx; // context for applied acceleration function
  // pointer to applied acceleration function
  void (*app_accel)(double t, const double *xn, double *aout, void *ctx);
  bool app_accel_evolve; // set to true if applied acceleration function is time dependent

  void *hamil_ctx; // context for hamiltonian function
  // pointer to hamilonian function
  void (*hamil)(double t, const double *xn, double *aout, void *ctx);

  void *h_ij_ctx; // context for spatial metric function (covariant)
  // pointer to metric (covariant components) function
  void (*h_ij)(double t, const double *xn, double *aout, void *ctx);

  void *h_ij_inv_ctx; // context for spatial metric function (contravariant)
  // pointer to metric (contravaraint components) function
  void (*h_ij_inv)(double t, const double *xn, double *aout, void *ctx);

  void *det_h_ctx; // context for determinant of the spatial metric
  // pointer to the determinant of the spatial metric
  void (*det_h)(double t, const double *xn, double *aout, void *ctx);

  bool output_f_lte; // Boolean for writing out f_lte (used for calculating transport coeff.)
  double iter_eps; // error tolerance for moment fixes of f_lte (density is always exact)
  int max_iter; // maximum number of iterations for correction output f_lte
  bool use_last_converged; // use last iteration value regardless of convergence for f_lte?

  // boundary conditions
  struct gkyl_vlasov_bcs bcx, bcy, bcz;
};

// Parameter for EM field
struct gkyl_vlasov_field {
  enum gkyl_field_id field_id; // type of field 
                               // (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
  bool is_static; // set to true if field does not change in time

  double epsilon0, mu0;
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  void *ext_em_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
  void (*ext_em)(double t, const double *xn, double *ext_em_out, void *ctx);
  bool ext_em_evolve; // set to true if external electromagnetic field function is time dependent

  void *app_current_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
  void (*app_current)(double t, const double *xn, double *app_current_out, void *ctx);
  bool app_current_evolve; // set to true if applied current function is time dependent
  
  double limiter_fac; // Optional input parameter for adjusting diffusion in slope limiter
  bool limit_em; // Optional input parameter for applying limiters to EM fields

  bool use_ghost_current; // Are we using ghost currents to correct dE/dt = -J in 1x
  
  // boundary conditions
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];

  // Options for Vlasov-Poisson.
  struct gkyl_poisson_bc poisson_bcs; // Boundary conditions for Poisson eqn.

  void *external_potentials_ctx; // Context for external (phi,A) potentials.
  // Pointer to function defining external potentials (phi,A).
  void (*external_potentials)(double t, const double *xn, double *ext_pot, void *ctx);
  bool external_potentials_evolve; // True if external potentials are time dependent.
};

// Parameter for Vlasov fluid species
struct gkyl_vlasov_fluid_species {
  char name[128]; // species name

  double charge, mass; // charge and mass
  
  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // Hyper-diffusion coefficient
  double nuHyp;

  struct gkyl_wv_eqn *equation; // equation object
  double limiter_fac; // Optional input parameter for adjusting diffusion in slope limiter

  // advection coupling (if scalar advection)
  // gkyl_eqn_type eqn_id = GKYL_EQN_ADVECTION
  struct gkyl_vlasov_fluid_advection advection;

  // source term
  struct gkyl_vlasov_fluid_source source;
  
  // diffusion coupling to include
  struct gkyl_vlasov_fluid_diffusion diffusion;

  void *can_pb_n0_ctx; // context for background density function in canonical PB fluid systems
  // pointer to background density function in canonical PB fluid systems
  void (*can_pb_n0)(double t, const double *xn, double *aout, void *ctx);

  void *app_accel_ctx; // context for applied acceleration function
  // pointer to applied acceleration function
  void (*app_accel)(double t, const double *xn, double *aout, void *ctx);
  bool app_accel_evolve; // set to true if applied acceleration function is time dependent
  
  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// Top-level app parameters
struct gkyl_vm {
  char name[128]; // name of app: used as output prefix

  int cdim, vdim; // conf, velocity space dimensions
  double lower[3], upper[3]; // lower, upper bounds of config-space
  int cells[3]; // config-space cells
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use

  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  double cfl_frac; // CFL fraction to use (default 1.0)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_vlasov_species species[GKYL_MAX_SPECIES]; // species objects

  int num_fluid_species; // number of fluid species
  struct gkyl_vlasov_fluid_species fluid_species[GKYL_MAX_SPECIES]; // fluid species objects
  
  bool skip_field; // Skip field update or no field specified
  struct gkyl_vlasov_field field; // field object
  bool is_electrostatic; // Indicate whether to use Vlasov-Poisson.

  struct gkyl_app_parallelism_inp parallelism; // Parallelism-related inputs.
};

// Simulation statistics
struct gkyl_vlasov_stat {
  bool use_gpu; // did this sim use GPU?
  
  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method
    
  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure
    
  double total_tm; // time for simulation (not including ICs)
  double rk3_tm; // time for SSP RK3 step
  double fl_em_tm; // time for implicit fluid-EM coupling step
  double init_species_tm; // time to initialize all species
  double init_fluid_species_tm; // time to initialize all fluid species
  double init_field_tm; // time to initialize fields

  double species_rhs_tm; // time to compute species collisionless RHS
  double fluid_species_rhs_tm; // time to compute fluid species RHS
  double fluid_species_vars_tm; // time to compute fluid variables (flow velocity and pressure)
    
  double species_coll_mom_tm; // time needed to compute various moments needed in LBO
  double species_lbo_coll_drag_tm[GKYL_MAX_SPECIES]; // time to compute LBO drag terms
  double species_lbo_coll_diff_tm[GKYL_MAX_SPECIES]; // time to compute LBO diffusion terms
  double species_coll_tm; // total time for collision updater (excluded moments)

  double species_rad_tm; // total time for radiation updater 

  double species_lte_tm; // time needed to compute the lte equilibrium

  long niter_self_bgk_corr[GKYL_MAX_SPECIES]; // number of iterations used to correct self collisions in BGK

  double species_bc_tm; // time to compute species BCs
  double fluid_species_bc_tm; // time to compute fluid species BCs
  double field_bc_tm; // time to compute field
  
  double field_rhs_tm; // time to compute field RHS
  double current_tm; // time to compute currents and accumulation

  long n_species_omega_cfl; // number of times CFL-omega all-reduce is called
  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl

  long n_field_omega_cfl; // number of times CFL-omega for field all-reduce is called
  double field_omega_cfl_tm; // time spent in all-reduce for omega-cfl for field

  long n_mom; // calls to moment calculation
  double mom_tm; // time to compute moments

  long n_diag; // calls to diagnostics
  double diag_tm; // time to compute diagnostics

  long n_io; // number of calls to IO
  double io_tm; // time to perform IO
};

// Object representing Vlasov app
typedef struct gkyl_vlasov_app gkyl_vlasov_app;

/**
 * Construct a new Vlasov app.
 *
 * @param vm App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New vlasov app object.
 */
gkyl_vlasov_app* gkyl_vlasov_app_new(struct gkyl_vm *vm);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0);

/**
 * Initialize field by projecting initial conditions on basis
 * functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_vm object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0);

/**
 * Initialize fluid species by projecting initial conditions on basis
 * functions. Fluid species index (sidx) is the same index used to specify
 * the species in the gkyl_vm object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of fluid species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_fluid_species(gkyl_vlasov_app* app, int sidx, double t0);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_app_from_file_field(gkyl_vlasov_app *app, const char *fname);

/**
 * Initialize Vlasov species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_vlasov_app_from_file_species(gkyl_vlasov_app *app, int sidx,
  const char *fname);

/**
 * Initialize fluid species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_vlasov_app_from_file_fluid_species(gkyl_vlasov_app *app, int sidx,
  const char *fname);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_field(gkyl_vlasov_app *app, int frame);

/**
 * Initialize Vlasov species from frame
 *
 * @param app App object
 * @param sidx gk species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_species(gkyl_vlasov_app *app, int sidx, int frame);

/**
 * Initialize fluid species from frame
 *
 * @param app App object
 * @param sidx gk species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_fluid_species(gkyl_vlasov_app *app, int sidx, int frame);

/**
 * Initialize the Vlasov app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_app_read_from_frame(gkyl_vlasov_app *app, int frame);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_calc_mom(gkyl_vlasov_app *app);

/**
 * Calculate integrated diagnostic moments.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app* app, double tm);

/**
 * Calculate integrated L2 norm of the distribution function, f^2.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_integrated_L2_f(gkyl_vlasov_app* app, double tm);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app* app, double tm);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame);

/**
 * Write species data to file - for the local equilbrium.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species_lte(gkyl_vlasov_app* app, int sidx, double tm, int frame);

/**
 * Write species data to file for H, G, drag coeff, diff coeff.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species_fpo(gkyl_vlasov_app* app, int sidx, double tm, int frame);

/**
 * Write fluid species data to file. 
 * 
 * @param app App object.
 * @param sidx Index of fluid species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_fluid_species(gkyl_vlasov_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_mom(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write integrated diagnostic moments for species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app);

/**
 * Write integrated diagnostic quantities for fluid species to file. Integrated
 * quantities are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_app_write_fluid_integrated_mom(gkyl_vlasov_app *app);

/**
 * Write integrated L2 norm of the species distribution function to file. Integrated
 * L2 norm is appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_app_write_integrated_L2_f(gkyl_vlasov_app *app);

/**
 * Write integrated correct lte status of the species distribution function to file. Correct
 * lte status is appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_app_write_lte_corr_status(gkyl_vlasov_app *app);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app* app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_stat_write(gkyl_vlasov_app* app);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_vlasov_app_cout(const gkyl_vlasov_app* app, FILE *fp, const char *fmt, ...);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status gkyl_vlasov_update(gkyl_vlasov_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_vlasov_stat gkyl_vlasov_app_stat(gkyl_vlasov_app* app);

/**
 * Run the RHS for the species update. This is used to compute kernel
 * timers and is not otherwise a useful function for a full
 * simulation.
 *
 * @param app App object.
 * @param update_vol_term Set to 1 to update vol term also, 0 otherwise
 */
void gkyl_vlasov_app_species_ktm_rhs(gkyl_vlasov_app* app, int update_vol_term);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_app_release(gkyl_vlasov_app* app);
