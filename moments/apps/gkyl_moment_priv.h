// Private header for use in moment app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fv_proj.h>
#include <gkyl_kep_scheme.h>
#include <gkyl_mhd_src.h>
#include <gkyl_moment.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_mp_scheme.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_ten_moment_nn_closure.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_gr_spacetime.h>

// number of components that various applied functions should return
enum {
  GKYL_MOM_APP_NUM_APPLIED_CURRENT = 3,
  GKYL_MOM_APP_NUM_EXT_EM = 6,
  GKYL_MOM_APP_NUM_APPLIED_ACCELERATION = 3,
  GKYL_MOM_APP_NUM_NT_SOURCE = 2
};

// Species data
struct moment_species {
  int ndim;
  char name[128]; // species name
  double charge, mass;

  bool is_static; // is the fluid static?

  // Does this fluid require source terms?
  // Set to true if any term which requires source terms is detected 
  // including app->has_field (because the fluids couple to EM fields via sources)
  // and various moment_species inputs such as has_friction, has_volume_sources, 
  // has_reactivity, has_app_accel, etc. 
  bool update_sources;
  
  double k0; // Closure parameter (default is 0.0, used by 10 moment).
  bool has_grad_closure; // Has gradient-based closure (only for 10 moment).
  bool has_nn_closure; // Has neural network-based closure (only for 10 moment).
  struct gkyl_kann_net* ann; // Neural network architecture.
  int poly_order; // Polynomial order of learned DG coefficients.
  enum gkyl_braginskii_type type_brag; // which Braginskii equations

  bool has_friction; // Run with frictional sources.
  bool use_explicit_friction; // Use an explicit (SSP-RK3) solver for integrating frictional sources.
  double friction_Z; // Ionization number for frictional sources.
  double friction_T_elc; // Electron temperature for frictional sources.
  double friction_Lambda_ee; // Electron-electron collisional term for frictional sources.

  bool has_volume_sources; // Run with volume-based geometrical sources.
  double volume_gas_gamma; // Adiabatic index for volume-based geometrical sources.
  double volume_U0; // Initial comoving plasma velocity for volume-based geometrical sources.
  double volume_R0; // Initial radial distance from expansion/contraction center for volume-based geometrical sources.

  bool has_reactivity; // Run with reactive sources.
  double reactivity_gas_gamma; // Adiabatic index for reactive sources.
  double reactivity_specific_heat_capacity; // Specific heat capacity for reactive sources.
  double reactivity_energy_of_formation; // Energy of formation for reactive sources.
  double reactivity_ignition_temperature; // Ignition temperature for reactive sources.
  double reactivity_reaction_rate; // Reaction rate for reactive sources.

  bool has_einstein_medium; // Run with coupled fluid-Einstein sources in plane-symmetric spacetimes.
  double medium_gas_gamma; // Adiabatic index for coupled fluid-Einstein sources in plane-symmetric spacetimes.
  double medium_kappa; // Stress-energy prefactor for coupled fluid-Einstein sources in plane-symmetric spacetimes.

  bool has_gr_ultra_rel; // Run with general relativistic source terms (Euler equations, ultra-relativistic equation of state).
  double gr_ultra_rel_gas_gamma; // Adiabatic index for general relativistic Euler equations (ultra-relativistic equation of state).

  bool has_gr_euler; // Run with general relativistic source terms (Euler equations, ideal gas equation of state).
  double gr_euler_gas_gamma; // Adiabatic index for general relativistic Euler equations (ideal gas equation of state).

  bool has_gr_twofluid; // Run with general relativistic two-fluid source terms.
  double gr_twofluid_mass_elc; // Electron mass for general relativistic two-fluid equations.
  double gr_twofluid_mass_ion; // Ion mass for general relativistic two-fluid equations.
  double gr_twofluid_charge_elc; // Electron charge for general relativistic two-fluid equations.
  double gr_twofluid_charge_ion; // Ion charge for general relativistic two-fluid equations.
  double gr_twofluid_gas_gamma_elc; // Adiabatic index for electrons in general relativistic two-fluid equations.
  double gr_twofluid_gas_gamma_ion; // Adiabatic index for ions in general relativistic two-fluid equations.
  double gr_twofluid_e_fact; // Electric field divergence error propagation speed for general relativistic two-fluid equations.

  bool has_vacuum_einstein; // Run with vacuum Einstein sources using the Bona-Masso formalism.
  double vacuum_einstein_excision_threshold; // Excision threshold (lapse) for vacuum Einstein equations using the Bona-Masso formalism.
  enum gkyl_spacetime_slicing vacuum_einstein_spacetime_slicing; // Spacetime slicing condition for vacuum Einstein equations using the Bona-Masso formalism.
  enum gkyl_spacetime_evolution vacuum_einstein_spacetime_evolution; // Spacetime evolution system for vacuum Einstein equations using the Bona-Masso formalism.

  bool has_vacuum_einstein_conformal; // Run with vacuum Einstein sources using the conformal Bona-Masso formalism.
  double vacuum_einstein_conformal_excision_threshold; // Excision threshold (lapse) for vacuum Einstein equations using the conformal Bona-Masso formalism.
  enum gkyl_spacetime_slicing vacuum_einstein_conformal_spacetime_slicing; // Spacetime slicing condition for vacuum Einstein equations using the conformal Bona-Masso formalism.
  enum gkyl_spacetime_evolution vacuum_einstein_conformal_spacetime_evolution; // Spacetime evolution system for vacuum Einstein equations using the conformal Bona-Masso formalism.
  
  bool has_gr_mhd; // Run with general relativistic source terms (general relativistic magnetohydrodynamics equations).
  double gr_mhd_gas_gamma; // Adiabatic index for general relativistic magnetohydrodynamics equations.

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  bool has_app_accel; // flag to indicate there is applied acceleration
  bool app_accel_evolve; // flag to indicate applied acceleration is time-dependent
  struct gkyl_array *app_accel; // applied acceleration
  gkyl_fv_proj *app_accel_proj; // projector for acceleration

  struct gkyl_array *nT_source; // array for num density and temperature sources
  // projection func for num density and temperature sources
  gkyl_fv_proj *proj_nT_source;
  bool nT_source_set_only_once; // set by user
  bool nT_source_is_set; // to be set at run time

  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  struct gkyl_array *embed_mask;

  enum gkyl_eqn_type eqn_type;  // type ID of equation
  int num_equations;            // number of equations in species
  struct gkyl_wv_eqn *equation; // equation object

  enum gkyl_moment_scheme scheme_type; // scheme to update equations

  // 
  // solvers and data to update fluid equations
  union {
    struct {
      gkyl_wave_prop *slvr[3];        // wave-prop solver in each direction
      struct gkyl_array *fdup, *f[4]; // arrays for updates
    };
    struct {
      union {
        gkyl_mp_scheme *mp_slvr;   // monotonicity-preserving scheme
        gkyl_kep_scheme *kep_slvr; // KEP scheme
      };
      struct gkyl_array *f0, *f1, *fnew; // arrays for updates
      struct gkyl_array *cflrate;        // CFL rate in each cell
      struct gkyl_array *alpha;          // for shock detector
    };
  };
  struct gkyl_array *fcurr; // points to current solution (depends on scheme)

  // boundary condition type
  enum gkyl_species_bc_type lower_bct[3], upper_bct[3];
  // boundary condition solvers on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_q; // integrated conserved quantities
  bool is_first_q_write_call; // flag for dynvec written first time
};

// Field data
struct moment_field {
  int ndim;
  double epsilon0, mu0;
  bool is_static; // is the field static?

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  struct gkyl_wv_eqn *maxwell; // pointer to Maxwell eqn obj

  bool has_ext_em; // flag to indicate there is external electromagnetic field
  bool ext_em_evolve; // flag to indicate external electromagnetic field is time dependent
  struct gkyl_array *ext_em; // external electromagnetic field
  gkyl_fv_proj *ext_em_proj; // projector for external electromagnetic field 
  double t_ramp_E; // linear ramp for turning on external E field

  bool has_app_current; // flag to indicate there is an applied current 
  bool app_current_evolve; // flag to indicate applied current is time dependent  
  struct gkyl_array *app_current; // applied current
  gkyl_fv_proj *app_current_proj;  // projector for applied current 
  double t_ramp_curr; // linear ramp for turning on applied currents

  bool use_explicit_em_coupling; // flag to indicate if em coupling should be explicit, defaults implicit
  struct gkyl_array *app_current1; // additional array for applied currents (for use_explicit_em_coupling stages)
  struct gkyl_array *app_current2; // additional array for applied currents (for use_explicit_em_coupling stages)

  bool has_volume_sources; // Run with volume-based geometrical sources.
  double volume_gas_gamma; // Adiabatic index for volume-based geometrical sources.
  double volume_U0; // Initial comoving plasma velocity for volume-based geometrical sources.
  double volume_R0; // Initial radial distance from expansion/contraction center for volume-based geometrical sources.

  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  struct gkyl_array *embed_mask;

 // scheme to update equations solvers and data to update fluid
 // equations
  enum gkyl_moment_scheme scheme_type;
  union {
    struct {
      gkyl_wave_prop *slvr[3]; // wave-prop solver in each direction
      struct gkyl_array *fdup, *f[4]; // arrays for updates
    };
    struct {
      gkyl_mp_scheme *mp_slvr; // monotonicity-preserving scheme
      struct gkyl_array *f0, *f1, *fnew; // arrays for updates
      struct gkyl_array *cflrate; // CFL rate in each cell
    };
  };
  struct gkyl_array *fcurr; // points to current solution (depends on scheme)

  // boundary condition type
  enum gkyl_field_bc_type lower_bct[3], upper_bct[3];
  // boundary conditions on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for dynvec written first time
};

// Source data
struct moment_coupling {
  // grid for braginskii variables (braginskii variables located at cell nodes)  
  struct gkyl_rect_grid non_ideal_grid;
  // local, local-ext ranges for braginskii variables (loop over nodes)  
  struct gkyl_range non_ideal_local, non_ideal_local_ext;

  // Gradient-based closure solver (if present).
  struct gkyl_ten_moment_grad_closure *grad_closure_slvr[GKYL_MAX_SPECIES];
  // Neural network-based closure solver (if present).
  struct gkyl_ten_moment_nn_closure *nn_closure_slvr[GKYL_MAX_SPECIES];
  // Braginskii solver (if present).
  struct gkyl_moment_braginskii *brag_slvr; 

  // array for stable time-step from non-ideal terms  
  struct gkyl_array *non_ideal_cflrate[GKYL_MAX_SPECIES];
  // array for non-ideal variables 
  // Braginskii variables, viscous stress tensor and heat-flux vector for Euler/Isothermal Euler
  // heat-flux tensor for ten-moment  
  struct gkyl_array *non_ideal_vars[GKYL_MAX_SPECIES];
  // array for storing RHS of each species from non-ideal term updates 
  // Braginskii tranport for Euler/Isothermal Euler
  // Gradient-based closure for ten-moment
  struct gkyl_array *pr_rhs[GKYL_MAX_SPECIES];
  // array for storing RHS of number density and temperature source terms
  struct gkyl_array *nT_sources[GKYL_MAX_SPECIES];

  gkyl_moment_em_coupling *slvr; // source solver function
};

struct mhd_src {
  gkyl_mhd_src *slvr; // source solver function
};

// Moment app object: used as opaque pointer in user code
struct gkyl_moment_app {
  char name[128]; // name of app
  int ndim; // space dimensions
  double tcurr; // current time
  double cfl; // CFL number

  enum gkyl_moment_scheme scheme_type; // scheme to use
  enum gkyl_wave_split_type split_type; // edge splitting to use  
  // 
  enum gkyl_mp_recon mp_recon; // reconstruction scheme to use
 // should shock-hybrid scheme be used when using KEP?  
  bool use_hybrid_flux_kep;

  bool has_braginskii; // has Braginskii transport
  double coll_fac; // multiplicative collisionality factor for Braginskii  

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
  int nghost[3]; // number of ghost-cells in each direction

  int is_dir_skipped[3]; // flags to tell if update in direction are skipped

  struct gkyl_rect_grid grid; // grid
  struct gkyl_range local, local_ext; // local, local-ext ranges
  struct gkyl_range global, global_ext; // global, global-ext ranges

  struct gkyl_rect_decomp *decomp; // decomposition object
  struct gkyl_comm *comm;   // communicator object

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx;   // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom *geom; // geometry needed for species and field solvers

  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  int has_field; // flag to indicate if we have a field
  struct moment_field field; // field data

  // species data
  int num_species;
  struct moment_species *species; // species data

  // work arrays for use in the KEP and MP scheme: these are stored
  // here so they can be reused
  struct {
    struct gkyl_array *ql, *qr;     // expansions on left/right edge of cell
    struct gkyl_array *amdq, *apdq; // minus/plus fluctuations
  };
  
  int update_sources; // flag to indicate if sources are to be updated
  struct moment_coupling sources; // sources

  int update_mhd_source;
  struct mhd_src mhd_source;

  struct gkyl_moment_stat stat; // statistics

  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_moment_app *app, double dt0);

  bool has_collision; // has collisions
  // scaling factors for collision frequencies so that nu_sr=nu_base_sr/rho_s
  // nu_rs=nu_base_rs/rho_r, and nu_base_sr=nu_base_rs
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
};

// Meta-data for IO
struct moment_output_meta {
  int frame; // frame number
  double stime; // output time
};

/** Some common functions to species and fields */

// functions for use in integrated quantities calculation
static inline void
integ_unit(int nc, const double *qin, double *integ_out)
{
  for (int i = 0; i < nc; ++i)
    integ_out[i] = qin[i];
}
static inline void
integ_sq(int nc, const double *qin, double *integ_out)
{
  for (int i = 0; i < nc; ++i)
    integ_out[i] = qin[i] * qin[i];
}

// function for copy BC
static inline void
bc_copy(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin,
  double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int c = 0; c < nc; ++c)
    ghost[c] = skin[c];
}

// function for skip BCs
static inline void
bc_skip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin,
  double *GKYL_RESTRICT ghost, void *ctx)
{
}

// Compute integrated quantities specified by i_func
void calc_integ_quant(const struct gkyl_wv_eqn *eqn, double vol,
  const struct gkyl_array *q,
  const struct gkyl_wave_geom *geom,
  struct gkyl_range update_rng, double *integ_q);

// Check array "q" for nans
bool check_for_nans(const struct gkyl_array *q, struct gkyl_range update_rng);

// Apply periodic BCs to corner cells of "f" (ONLY WORKS IN 2D)
void moment_apply_periodic_corner_sync_2d(const gkyl_moment_app *app,
  struct gkyl_array *f);

// Apply wedge-periodic BCs to array "f"
void moment_apply_wedge_bc(const gkyl_moment_app *app, double tcurr,
  const struct gkyl_range *update_rng,
  struct gkyl_array *bc_buffer, int dir,
  const struct gkyl_wv_apply_bc *lo,
  const struct gkyl_wv_apply_bc *up,
  struct gkyl_array *f);

/**
 * Return ghost cell layout for grid.
 *
 * @param app App object.
 * @param nghost On output, ghost-cells used for grid.
 *
 */
void gkyl_moment_app_nghost(gkyl_moment_app *app, int nghost[3]);

/** moment_species API */

// Initialize the moment species object
void moment_species_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app,
  struct moment_species *sp);

// Apply BCs to species data "f"
void moment_species_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_species *sp,
  struct gkyl_array *f);

// Maximum stable time-step from species
double moment_species_max_dt(const gkyl_moment_app *app,
  const struct moment_species *sp);

// Advance solution of species by time-step dt to tcurr+dt
struct gkyl_update_status moment_species_update(gkyl_moment_app *app,
  struct moment_species *sp,
  double tcurr, double dt);

// Compute RHS of moment equations
double moment_species_rhs(gkyl_moment_app *app, struct moment_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

// Free memory allocated by species
void moment_species_release(const struct moment_species *sp);

/** moment_field API */

// Initialize EM field
void moment_field_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_field *mom_fld,
  struct gkyl_moment_app *app, struct moment_field *fld);

// Apply BCs to EM field
void moment_field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field,
  struct gkyl_array *f);

// Maximum stable time-step due to EM fields
double moment_field_max_dt(const gkyl_moment_app *app,
  const struct moment_field *fld);

// Update EM field from tcurr to tcurr+dt
struct gkyl_update_status moment_field_update(gkyl_moment_app *app,
  const struct moment_field *fld,
  double tcurr, double dt);

// Compute RHS of EM equations
double moment_field_rhs(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

// Release the EM field object
void moment_field_release(const struct moment_field *fld);

/** moment_coupling API */

// initialize source solver: this should be called after all species
// and fields are initialized
void moment_coupling_init(const struct gkyl_moment_app *app,
                          struct moment_coupling *src);

/** mhd_src functions */

void mhd_src_init(const struct gkyl_moment_app *app,
  const struct gkyl_moment_species *sp, struct mhd_src *src);

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
void mhd_src_update(gkyl_moment_app *app, struct mhd_src *src, int nstrang,
  double tcurr, double dt);

void mhd_src_release(const struct mhd_src *src);

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
struct gkyl_update_status moment_coupling_update(gkyl_moment_app *app,
  struct moment_coupling *src, int nstrang, double tcurr, double dt);

// Release coupling sources
void moment_coupling_release(const struct gkyl_moment_app *app,
  const struct moment_coupling *src);

/** Top-level app API */

// Take a single time-step using a single-step time-stepper
struct gkyl_update_status moment_update_one_step(gkyl_moment_app *app,
  double dt0);

// Take a single time-step using a SSP-RK3 stepper
struct gkyl_update_status moment_update_ssp_rk3(gkyl_moment_app *app,
  double dt0);

/**
 * Create new array meta header from input struct. Free returned
 * object with moment_array_meta_release.
 *
 * @param meta Meta-data for output.
 * @return New meta object to pass to write method.
 */
struct gkyl_msgpack_data* moment_array_meta_new(struct moment_output_meta meta);

/**
 * Release meta struct
 *
 * @param mt Meta object to free
 */
void moment_array_meta_release(struct gkyl_msgpack_data *mt);

/**
 * Read meta-data from mpack formated binary input
 *
 * @param mt Mpack encoded meta-data
 * @return Meta-data for simulation
 */
struct moment_output_meta moment_meta_from_mpack(struct gkyl_msgpack_data *mt);
