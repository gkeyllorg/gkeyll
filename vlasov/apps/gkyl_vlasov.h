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
#include <gkyl_vlasov_triad_geom.h>

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
  double nu_frac; // Rescales collision frequencies (default = 1).
  bool write_coll_diagnostics; // Whether to output diagnostics.

  // Function for computing self-collision frequency.
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);
  void *self_nu_ctx; // Context for self_nu.

  int num_cross_collisions; // Number of species to collide with.
  char collide_with[GKYL_MAX_SPECIES][128]; // Names of species to collide with.
  // Functions for computing cross-collision frequencies (one for each num_cross_collisions).
  // Either all partners specify a function or none do (Spitzer cross frequencies).
  evalf_t cross_nu[GKYL_MAX_SPECIES];
  void *cross_nu_ctx[GKYL_MAX_SPECIES]; // Context for cross_nu.

  // Parameters used to compute the Coulomb Logarithm.
  double den_ref; // Reference density.
  double temp_ref; // Reference temperature.
  double hbar, eps0, eV; // Planck's constant/2 pi, vacuum permittivity, elementary charge.

  // BGK collisions specific inputs
  bool fixed_temp_relax; // Are BGK collisions relaxing to a fixed input temperature?
  bool is_implicit; // Boolean for using implicit BGK collisions (replaces rk3).
};

// Parameters for species radiation
struct gkyl_vlasov_radiation {
  enum gkyl_vlasov_radiation_id radiation_id; // type of radiation (see gkyl_eqn_type.h)
  double t_cool; // Cooling time in radiation operator rad_force ~ -1/t_cool*drag
  double p0; // (four-) velocity to cool to.
};

// Parameters for species source
struct gkyl_vlasov_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // Are we writing out the source?
  bool evolve_source; // Are our sources time-dependent?
  int num_sources;

  double source_length; // required for boundary flux source
  char source_species[128];

  // Sources using projection routine.
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];

  // Adaptive source input parameters.
  int num_cross_source; // Number of species that we are sourcing with.
  char source_with[GKYL_MAX_SPECIES][128]; // Names of species that we are using for cross sources.
  double source_with_v_thresh
    [GKYL_MAX_SPECIES]; // Threshold velocity if re-scaling density based on partial moments.
  double source_with_f_thresh
    [GKYL_MAX_SPECIES]; // Threshold f for accumulating partial moments for re-scaling density.
  bool source_with_upper_half
    [GKYL_MAX_SPECIES]; // Are you using the upper-half or lower-half plane for partial moments?
  int source_with_proj
    [GKYL_MAX_SPECIES]; // Which projection function is being used with this adaptive source?
  bool filter; // Are we filtering the rescaled density?
  int num_filters; // Are we filtering repeatedly?
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
  void *Dij_ctx; // context for applied diffusion function if using general diffusion tensor
  // pointer to applied diffusion function is using general diffusion tensor
  void (*Dij)(double t, const double *xn, double *Dout, void *ctx);
};

struct gkyl_vlasov_correct_inp {
  bool correct_all_moms; // Boolean if we are correcting all the moments or only density.
  double iter_eps; // Error tolerance for moment fixes (density is always exact).
  int max_iter; // Maximum number of iteration
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
    // *even if* the scheme fails to converge.
  bool output_f_lte; // Boolean for writing out f_lte (used for calculating transport coeff.).
};

struct vlasov_mapc2p_vel {
  void *mapc2p_vel_ctx; // context for mapc2p function for velocity space
  // pointer to mapc2p function for velocity space:
  // xc are the computational space coordinates and on output
  // xp are the corresponding physical space coordinates.
  void (*mapc2p_vel_func)(double t, const double *xc, double *xp, void *ctx);
};

struct vlasov_mapc2p_pos {
  void *mapc2p_pos_ctx; // context for mapc2p function for configuration space
  // pointer to (per-direction) mapc2p function for configuration space:
  // xc is the computational space coordinate and on output xp is the
  // corresponding physical space coordinate. NULL => identity map.
  void (*mapc2p_pos_func)(double t, const double *xc, double *xp, void *ctx);
};

// Parameters for Vlasov geometry.
struct gkyl_vlasov_geom {
  double spin_bh, mass_bh; // Black hole spin and mass.
  bool use_preset_geom; // bool to determine if we are using triad input geom
  enum gkyl_triad_preset_geom_type
    triad_preset_geom_type; // geom type for preset geometries for triads
};

// Which evolved quantities a species owns: a kinetic distribution function, a
// fluid moment vector, or both.
enum gkyl_species_type {
  GKYL_SPECIES_VLASOV, // kinetic distribution only
  GKYL_SPECIES_FLUID, // fluid moments only
  GKYL_SPECIES_PKPM // both (reserved; not supported by the Vlasov app)
};

// Parameters for the kinetic block of a Vlasov species: the velocity-space
// grid, distribution-function initial conditions, collisions, sources, and
// kinetic BCs. Declared as the 'kinetic' block of struct gkyl_vlasov_species.
struct gkyl_vlasov_kinetic_species {
  enum gkyl_model_id model_id; // Type of model
    // (e.g., SR, general geometry, see gkyl_eqn_type.h).

  double lower[3], upper[3]; // Lower, upper bounds of velocity-space.
  int cells[3]; // Velocity-space cells.

  // Velocity-space mapping in each velocity-space dimension
  struct vlasov_mapc2p_vel mapc2p_vel[GKYL_MAX_CDIM];

  // Initial conditions using projection routine.
  int num_init; // Number of initial condition functions.
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];

  bool is_static; // Set to true if species does not change in time.
  bool no_collisionless_terms; // Set to true to turn off collisionles terms.
  bool write_omega_cfl; // Whether to output dt diagnostic for the CFL constraint.
  bool write_cell_avg; // Boolean for only writing cell average of f.
  bool use_lo; // bool to determine if using low-order kernels for non-canonical Hamiltonian models.

  bool
    use_vierbein; // bool to determine if we are using vierbein input or by default using tangent vectors/ triads
  bool
    use_extended_hamil_def; // bool to determine if we are using the extended hamil definitions which includes potentials

  // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped.
  double skip_cell_thresh;

  int num_diag_moments; // number of diagnostic moments
  enum gkyl_distribution_moments diag_moments[16]; // list of diagnostic moments

  void *app_accel_ctx; // Context for applied acceleration function.
  // Pointer to applied acceleration function.
  void (*app_accel)(double t, const double *xn, double *aout, void *ctx);
  bool app_accel_evolve; // Set to true if applied acceleration function is time dependent.

  void *cov_tangent_basis_ctx; // Context for covariant tangent basis function.
  // Pointer to covariant tangent basis function.
  void (*cov_tangent_basis)(double t, const double *xn, double *aout, void *ctx);

  void *triad_basis_ctx; // Context for triad basis function.
  // Pointer to triad basis function.
  void (*triad_basis)(double t, const double *xn, double *aout, void *ctx);

  void *triad_basis_gradient_ctx; // Context for triad basis gradient function.
  // Pointer to triad basis gradient function.
  void (*triad_basis_gradient)(double t, const double *xn, double *aout, void *ctx);

  void *vierbein_ctx; // Context for vierbein function.
  // Pointer to vierbein function.
  void (*vierbein)(double t, const double *xn, double *aout, void *ctx);

  void *vierbein_gradient_ctx; // Context for vierbein gradient function.
  // Pointer to vierbein gradient function.
  void (*vierbein_gradient)(double t, const double *xn, double *aout, void *ctx);

  void *hamil_ctx; // Context for general (phase space) Hamiltonian function.
  // Pointer to general (phase space) Hamilonian function.
  void (*hamil)(double t, const double *xn, double *aout, void *ctx);

  void *h_ij_ctx; // Context for spatial metric function (covariant).
  // Pointer to metric (covariant components) function.
  void (*h_ij)(double t, const double *xn, double *aout, void *ctx);

  void *h_ij_inv_ctx; // Context for spatial metric function (contravariant).
  // Pointer to metric (contravariant components) function.
  void (*h_ij_inv)(double t, const double *xn, double *aout, void *ctx);

  void *det_h_ctx; // Context for determinant of the spatial metric.
  // Pointer to the determinant of the spatial metric.
  void (*det_h)(double t, const double *xn, double *aout, void *ctx);

  void *background_flows_ctx; // Context for background flows in can-pb
  // Pointer to the background flows.
  void (*background_flows)(double t, const double *xn, double *aout, void *ctx);

  void *effective_potential_ctx; // Context for effective potential in can-pb
  // Pointer to the effective potential
  void (*effective_potential)(double t, const double *xn, double *aout, void *ctx);

  // Input quantities used by LTE (local thermodynamic equilibrium) projection
  // This projection operator is used by BGK collisions.
  struct gkyl_vlasov_correct_inp correct;

  // Collisions to include.
  struct gkyl_vlasov_collisions collisions;

  // Radiation to include.
  struct gkyl_vlasov_radiation radiation;

  // Source to include.
  struct gkyl_vlasov_source source;

  // Boundary conditions.
  struct gkyl_vlasov_bcs bcx, bcy, bcz;
};

// Parameters for the fluid block of a Vlasov species: the equation object,
// fluid initial condition, advection/diffusion, and fluid BCs. Declared as the
// 'fluid' block of struct gkyl_vlasov_species.
struct gkyl_vlasov_fluid_species {
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

// Parameters for one Vlasov species. The type selects which block is read
// (GKYL_SPECIES_VLASOV: 'kinetic', GKYL_SPECIES_FLUID: 'fluid'); the identity
// (name, charge, mass) is declared here and not in the blocks. Names must be
// unique. For example, a kinetic species and a fluid species are declared as
//   struct gkyl_vlasov_kinetic_species elc = { .lower = ..., .cells = ..., ... };
//   struct gkyl_vlasov_fluid_species ion = { .equation = eqn, .init = ..., ... };
//   .num_species = 2,
//   .species = {
//     { .name = "elc", .charge = -1.0, .mass = 1.0, .type = GKYL_SPECIES_VLASOV, .kinetic = elc },
//     { .name = "ion", .charge = 1.0, .mass = 1836.0, .type = GKYL_SPECIES_FLUID, .fluid = ion },
//   },
struct gkyl_vlasov_species {
  char name[128]; // Species name.
  double charge, mass; // Charge and mass.

  enum gkyl_species_type type; // Which aspect(s) this species owns.

  struct gkyl_vlasov_kinetic_species
    kinetic; // Kinetic block (velocity grid, projections, collisions, ...).
  struct gkyl_vlasov_fluid_species fluid; // Fluid block (equation, init, diffusion, ...).
};

// Parameter for electromagnetic fields.
struct gkyl_vlasov_field {
  enum gkyl_field_id field_id; // Type of field
    // (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h).
  bool is_static; // Set to true if field does not change in time.

  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  // Correction speeds as a fraction of the speed of light for div(E)/div(B) errors.
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;
  double K_phi, K_psi; // Geometric source coefficients for electric/magnetic cleaning fields.

  void *ctx; // Context for initial condition init function for Vlasov-Maxwell.
  // Pointer to initialization function for Vlasov-Maxwell fields.
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  double limiter_fac; // Optional input parameter for adjusting diffusion in slope limiter.
  bool limit_em; // Optional input parameter for applying limiters to EM fields.

  bool use_ghost_current; // Are we using ghost currents to correct dE/dt = -J in 1x?
  bool use_geom_sources; // Are we using geometric sources to correct dE/dt = -J in 1x?

  // Vlasov-Maxwell boundary conditions.
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];

  // Vlasov-Poisson boundary conditions.
  struct gkyl_poisson_bc poisson_bcs;

  void *sigma_ctx; // context for resistive layer to damp EM fields
  // pointer to resistive layer function to damp EM fields
  void (*sigma)(double t, const double *xn, double *fout, void *ctx);

  void *ext_em_ctx; // Context for external electromagnetic fields function (E,B).
  // Pointer to external electromagnetic fields function.
  void (*ext_em)(double t, const double *xn, double *ext_em_out, void *ctx);
  bool ext_em_evolve; // Set to true if external electromagnetic fields are time dependent.

  void *app_current_ctx; // Context for applied currents function.
  // Pointer to applied currents function.
  void (*app_current)(double t, const double *xn, double *app_current_out, void *ctx);
  bool app_current_evolve; // Set to true if applied currents are time dependent.

  void *external_potentials_ctx; // Context for external (phi,A) potentials.
  // Pointer to function defining external potentials (phi,A).
  void (*external_potentials)(double t, const double *xn, double *ext_pot, void *ctx);
  bool external_potentials_evolve; // Set to true if external potentials are time dependent.
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

  // Per-direction non-uniform configuration-space mapping (C^0 piecewise
  // linear). Each direction with a NULL func is the identity map. This is a
  // diagonal coordinate stretch and is independent of the general curvilinear
  // mapc2p above (which feeds the wave_geom object).
  struct vlasov_mapc2p_pos mapc2p_pos[GKYL_MAX_CDIM];

  double cfl_frac; // CFL fraction to use (default 1.0)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  struct gkyl_vlasov_geom geom; // geom object

  // Species list: kinetic and fluid species in any order, at most
  // GKYL_MAX_SPECIES in total. Declaration order sets the index used by the
  // per-species app functions (files are named by the species name). Kinetic-
  // only diagnostics (distribution moments, L2 norm) skip fluid species.
  int num_species; // total number of species (kinetic + fluid)
  struct gkyl_vlasov_species species[GKYL_MAX_SPECIES]; // species objects

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

  long
    n_iter_corr[GKYL_MAX_SPECIES]; // total number of iterations used to correct species LTE projection
  long num_corr
    [GKYL_MAX_SPECIES]; // total number of times correction updater for species LTE projection is called
  double species_coll_mom_tm; // time needed to compute various moments needed in LBO
  double species_lbo_coll_drag_tm[GKYL_MAX_SPECIES]; // time to compute LBO drag terms
  double species_lbo_coll_diff_tm[GKYL_MAX_SPECIES]; // time to compute LBO diffusion terms
  double species_coll_tm; // total time for collision updater (excluded moments)

  double species_rad_tm; // total time for radiation updater

  double species_lte_tm; // time needed to compute the lte equilibrium

  double species_bc_tm; // time to compute species BCs
  double fluid_species_bc_tm; // time to compute fluid species BCs
  double field_bc_tm; // time to compute field

  double field_rhs_tm; // time to compute field RHS
  double current_tm; // time to compute currents and accumulation

  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl for species
  double species_io_tm; // Time to write the species distribution.
  double species_diag_calc_tm; // Time to compute species diagnostics.
  double species_diag_io_tm; // Time to write species diagnostics.

  double field_omega_cfl_tm; // time spent in all-reduce for omega-cfl for fields
  double field_io_tm; // Time to write the fields.
  double field_diag_calc_tm; // Time to compute field diagnostics.
  double field_diag_io_tm; // Time to write field diagnostics.

  long n_species_omega_cfl; // number of times CFL-omega all-reduce is called for species
  long n_field_omega_cfl; // number of times CFL-omega all-reduce is called for fields
  long n_mom; // total number of calls to moment updater routines
  long n_diag; // total number of calls to diagnostics
  long n_io; // number of calls to IO
  long n_diag_io; // number of calls to IO for diagnostics
  long n_field_diag; // total number of calls to diagnostics for field
  long n_field_io; // number of calls to IO for field
  long n_field_diag_io; // number of calls to IO for field diagnostics
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
gkyl_vlasov_app *gkyl_vlasov_app_new(struct gkyl_vm *vm);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_vlasov_app_apply_ic(gkyl_vlasov_app *app, double t0);

/**
 * Initialize field by projecting initial conditions on basis
 * functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app *app, double t0);

/**
 * Initialize a species by projecting initial conditions on basis functions.
 * Species are indexed in the order they appear in gkyl_vm's species[] list
 * (kinetic and fluid alike); see gkyl_vlasov_app_find_species to resolve a
 * species by name.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app *app, int sidx, double t0);

/**
 * Find a species by name.
 *
 * @param app App object.
 * @param nm Species name.
 * @return Index of the species (as used by the per-species app functions), -1 if not found.
 */
int gkyl_vlasov_app_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_file_field(
  gkyl_vlasov_app *app, const char *fname
);

/**
 * Initialize Vlasov species from file
 *
 * @param app App object
 * @param sidx Index of species (see gkyl_vlasov_app_apply_ic_species).
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_file_species(
  gkyl_vlasov_app *app, int sidx, const char *fname
);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_frame_field(gkyl_vlasov_app *app, int frame);

/**
 * Initialize Vlasov species from frame
 *
 * @param app App object
 * @param sidx Index of species (see gkyl_vlasov_app_apply_ic_species).
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_frame_species(
  gkyl_vlasov_app *app, int sidx, int frame
);

/**
 * Initialize the Vlasov app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_read_from_frame(gkyl_vlasov_app *app, int frame);

/**
 * Calculate integrated diagnostic moments.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app *app, double tm);

/**
 * Calculate integrated L2 norm of the distribution function, f^2.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_integrated_L2_f(gkyl_vlasov_app *app, double tm);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app *app, double tm);

/**
 * Write field and species data to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write field data to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_field(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write species data to file.
 *
 * @param app App object.
 * @param sidx Index of species to write (see gkyl_vlasov_app_apply_ic_species).
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species(gkyl_vlasov_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_mom(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write integrated diagnostic moments for all species (kinetic and fluid) to file.
 * Integrated moments are appended to the same file.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app);

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
void gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app *app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_stat_write(gkyl_vlasov_app *app);

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
void gkyl_vlasov_app_cout(const gkyl_vlasov_app *app, FILE *fp, const char *fmt, ...);

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
struct gkyl_update_status gkyl_vlasov_update(gkyl_vlasov_app *app, double dt);

/**
 * Return simulation statistics.
 *
 * @return Return statistics object.
 */
struct gkyl_vlasov_stat gkyl_vlasov_app_stat(gkyl_vlasov_app *app);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_app_release(gkyl_vlasov_app *app);
