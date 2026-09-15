// Private header for use in PKPM app: do not include in user-facing
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
#include <gkyl_bc_basic.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_pkpm_em_coupling.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_dist_vars.h>
#include <gkyl_dg_euler_pkpm.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_updater_diffusion_fluid.h>
#include <gkyl_dg_updater_lbo_pkpm.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_dg_updater_pkpm.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_dynvec.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_pkpm.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_pkpm.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_pkpm.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_pkpm.h>

// Definitions of private structs and APIs attached to these objects
// for use in PKPM app.

// Meta-data for IO
struct pkpm_output_meta {
  int frame; // frame number
  double stime; // output time
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
  char basis_type_nm[64]; // used during read
};

// data for moments
struct pkpm_species_moment {
  struct gkyl_dg_updater_moment *mcalc; // moment update

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

// forward declare species struct
struct pkpm_species;

struct pkpm_lbo_collisions {  
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *prim_moms, *nu_prim_moms; // LBO primitive moments
  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value
  struct gkyl_array *maxwellian_moms; // Maxwellian moments n, t, v_t^2.

  double betaGreenep1; // value of Greene's factor beta + 1
  double other_m[GKYL_MAX_SPECIES]; // masses of species being collided with
  struct gkyl_array *other_prim_moms[GKYL_MAX_SPECIES]; // self-primitive moments of species being collided with
  struct gkyl_array *cross_prim_moms[GKYL_MAX_SPECIES]; // LBO cross-primitive moments
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // LBO cross-species collision frequencies
  struct gkyl_array *other_nu[GKYL_MAX_SPECIES];
  struct gkyl_array *cross_nu_prim_moms; // weak multiplication of collision frequency and primitive moments
  
  struct gkyl_array *self_nu, *self_nu_prim_moms; // LBO self-primitive moments

  struct pkpm_species_moment moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)
  struct gkyl_array *m0;
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  gkyl_dg_bin_op_mem *greene_factor_mem; // memory needed in computing Greene factor
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct pkpm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

// species data
struct pkpm_species {
  struct gkyl_pkpm_species info; // data for species
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1, *fnew; // arrays for distribution function updates
  struct gkyl_array *fluid, *fluid1, *fluidnew; // arrays for momentum updates

  struct gkyl_array *cflrate_f; // CFL rate in each cell for distribution function update
  struct gkyl_array *cflrate_fluid; // CFL rate in each cell for momentum update

  struct gkyl_array *bc_buffer_dist; // buffer for BCs for distribution functions (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed_dist, *bc_buffer_up_fixed_dist; // fixed buffers for time independent BCs for distribution functions
  struct gkyl_array *bc_buffer_fluid; // buffer for BCs for momentum (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed_fluid, *bc_buffer_up_fixed_fluid; // fixed buffers for time independent BCs for momentum

  struct gkyl_array *f_host; // host copy of distribution function for use IO and initialization
  struct gkyl_array *fluid_host; // host copy of momentum for use IO and initialization

  // Duplicate copy of fluid data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  struct gkyl_array *fluid_dup;  

  struct gkyl_wv_eqn *equation; // For storing 10 moment equation object for upwinding fluid equations with Roe solve

  struct gkyl_array *qmem; // array for q/m*(E,B) for use in *explicit* update
                           // Note: this array is *only* used if the PKPM self-consistent EM fields are static
                           // If PKPM self-consistent EM fields are dynamics we utilize an implicit source update
                           // for the momentum equations and Ampere's Law

  struct pkpm_species_moment pkpm_moms; // for computing pkpm moments needed in update
  struct pkpm_species_moment pkpm_moms_diag; // for computing pkpm moments diagnostics

  // PKPM distribution function variables
  struct gkyl_array *g_dist_source; // g_dist_source = [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
                                    //                 (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]
  struct gkyl_array *F_k_p_1; // k+1 distribution function (first NP components are F_2) 
  struct gkyl_array *F_k_m_1; // k-1 distribution function (first NP components are F_1)

  // PKPM variables
  struct gkyl_array *pkpm_div_ppar; // div(p_parallel b_hat) used for computing self-consistent total pressure force 
  struct gkyl_array *pkpm_prim; // [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
  struct gkyl_array *pkpm_prim_surf; // Surface primitive variables. Ordered as:
                                     // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 3.0*Txx_xl/m, 3.0*Txx_xr/m, 
                                     //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 
                                     //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
  struct gkyl_array *pkpm_u; // [ux, uy, uz]
  struct gkyl_array *pkpm_p_ij; // (p_par - p_perp) b_i b_j + p_perp g_ij
  struct gkyl_array *cell_avg_prim; // Integer array for whether rho, p_par, or p_perp < 0.0 at control points
                                    // *only* currently used for diagnostic purposes
  struct gkyl_array *pkpm_accel; // Acceleration variables for PKPM, pkpm_accel:
                                 // 0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
                                 // 1: bb_grad_u (bb : grad(u))
                                 // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
                                 // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2*nu)
  struct gkyl_array *integ_pkpm_mom; // integrated PKPM variables [rho, rho ux, rho uy, rho uz, rho ux^2, rho uy^2, rho uz^2, p_par, p_perp]

  // Note: Each cell owns the *lower* surface expansions for pkpm_lax and pkpm_penalization 
  struct gkyl_array *pkpm_lax; // Surface expansion of Lax penalization lambda_i = |u_i| + sqrt(3.0*T_ii/m)
  struct gkyl_array *pkpm_penalization; // Surface expansion of total momentum penalization

  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars; // Updater to compute PKPM variables (primitive and acceleration variables)
  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars_ext; // Updater to compute PKPM variables (primitive and acceleration variables)
                                                     // over extended range (used when BCs are not absorbing to minimize apply BCs calls)
  struct gkyl_dg_calc_pkpm_dist_vars *calc_pkpm_dist_vars; // Updater to compute PKPM distribution function variables 
                                                           // div(p_parallel b_hat) and distribution function sources

  bool limit_fluid; // boolean for whether or not we are limiting fluid variables

  // Pointers for io for PKPM fluid variables, handled by kinetic species because of fluid-kinetic coupling.
  // For PKPM we construct the 10 moment conserved variables for ease of analysis 
  // along with an array of the various update variables, primitive and acceleration
  struct gkyl_array *fluid_io;
  struct gkyl_array *fluid_io_host;
  struct gkyl_array *pkpm_vars_io;
  struct gkyl_array *pkpm_vars_io_host;

  struct gkyl_array *L2_f; // L2 norm f^2
  double *red_L2_f; // for reduction of integrated L^2 norm on GPU
  double *red_integ_diag; // for reduction of integrated moments on GPU
  gkyl_dynvec integ_L2_f; // integrated L^2 norm reduced across grid
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_L2_write_call; // flag for integrated L^2 norm dynvec written first time
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  gkyl_dg_updater_pkpm *slvr; // PKPM solver for both momentum and kinetic equation

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC on distribution functions and momentum.
  struct gkyl_bc_basic *bc_lo_dist[3];
  struct gkyl_bc_basic *bc_up_dist[3];
  struct gkyl_bc_basic *bc_lo_fluid[3];
  struct gkyl_bc_basic *bc_up_fluid[3];
  // To simplify BC application, store local skin and ghost ranges for distribution functions
  struct gkyl_range lower_skin_dist[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost_dist[GKYL_MAX_DIM];
  struct gkyl_range upper_skin_dist[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost_dist[GKYL_MAX_DIM];
  bool bc_is_absorb; // boolean for absorbing BCs since 1/rho is undefined in absorbing BCs
                     // If BCs are *not* absorbing, primitive variables can be calculated on *extended* range 

  bool has_app_accel; // flag to indicate there is applied acceleration
  bool app_accel_evolve; // flag to indicate applied acceleration is time-dependent
  struct gkyl_array *app_accel; // applied acceleration
  struct gkyl_array *app_accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_accel_proj; // projector for acceleration

  enum gkyl_collision_id collision_id; // type of collisions
  struct pkpm_lbo_collisions lbo; // collisions object

  // fluid diffusion
  bool has_diffusion; // flag to indicate there is applied diffusion
  struct gkyl_array *diffD; // array for diffusion tensor
  struct gkyl_dg_updater_diffusion_fluid *diff_slvr; // Fluid diffusion equation solver

  double *omegaCfl_ptr_dist;
  double *omegaCfl_ptr_fluid;
};

// field data
struct pkpm_field {
  struct gkyl_pkpm_field info; // data for field

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *em, *em1, *emnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *em_host;  // host copy for use IO and initialization

  // Duplicate copy of EM data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  struct gkyl_array *em_dup;

  // Identity position map (PKPM runs on uniform grids); the Maxwell curl
  // kernels index its configuration-space Jacobian per conf cell.
  struct gkyl_vlasov_position_map *pos_map;

  bool has_ext_em; // flag to indicate there is external electromagnetic field
  bool ext_em_evolve; // flag to indicate external electromagnetic field is time dependent
  struct gkyl_array *ext_em; // external electromagnetic field
  struct gkyl_array *ext_em_host; // host copy for use in IO and projecting
  struct gkyl_array *tot_em; // total electromagnetic field
  gkyl_proj_on_basis *ext_em_proj; // projector for external electromagnetic field 

  bool has_app_current; // flag to indicate there is an applied current 
  bool app_current_evolve; // flag to indicate applied current is time dependent
  struct gkyl_array *app_current; // applied current
  struct gkyl_array *app_current_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_current_proj; // projector for applied current 

  struct gkyl_array *cell_avg_magB2; // Integer array for whether |B|^2 *only* uses cell averages for weak division
                                     // Determined when constructing the matrix if |B|^2 < 0.0 at control points
  struct gkyl_array *bvar; // magnetic field unit vector and tensor (diagnostic and for use in pkpm model)
  struct gkyl_array *bvar_surf; // Surface expansion magnetic field unit vector and tensor (for use in pkpm model)
  struct gkyl_array *div_b; // Volume expansion of div(b) (for use in pkpm model)
  struct gkyl_array *max_b; // max(|b_i|) penalization (for use in pkpm model)
  struct gkyl_dg_calc_em_vars *calc_bvar; // Updater to compute magnetic field unit vector and tensor

  bool limit_em; // boolean for whether or not we are limiting EM fields
  struct gkyl_dg_calc_em_vars *calc_em_vars; // Updater to limit EM fields 

  gkyl_hyper_dg *slvr; // Maxwell solver

  struct gkyl_array *em_energy; // EM energy components in each cell
  double *em_energy_red; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components

  bool is_first_energy_write_call; // flag for energy dynvec written first time

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_field_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];

  double* omegaCfl_ptr;
};

// fluid-EM coupling data in PKPM system
struct pkpm_fluid_em_coupling {
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each species
  struct gkyl_dg_calc_pkpm_em_coupling* slvr; // fluid-EM coupling solver
};

// PKPM object: used as opaque pointer in user code
struct gkyl_pkpm_app {
  char name[128]; // name of app
  struct gkyl_job_pool *job_pool; // Job pool
  
  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  double tcurr; // current time
  double cfl; // CFL number

  bool use_gpu; // should we use GPU (if present)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  struct gkyl_rect_grid grid; // config-space grid
  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges  
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  struct gkyl_basis basis, confBasis, velBasis; // phase-space, conf-space basis, vel-space basis

  struct gkyl_rect_decomp *decomp; // decomposition object
  struct gkyl_comm *comm;   // communicator object for conf-space arrays

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx;   // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom *geom; // geometry needed for species and field solvers (*only* p=1 right now JJ: 05/03/24)
  
  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *confBasis;
  } basis_on_dev;


  struct pkpm_field *field; // pointer to field object
  // species data
  int num_species;
  struct pkpm_species *species; // data for each species

  bool use_explicit_source; // Boolean to turn on explicit fluid-EM coupling
  struct pkpm_fluid_em_coupling *pkpm_em; // fluid-EM coupling data
  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_pkpm_app *app, double dt0);
  
  struct gkyl_pkpm_stat stat; // statistics
};

// Take a single forward Euler step of the PKPM system with the suggested time-step dt. 
void pkpm_forward_euler(gkyl_pkpm_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, 
  struct gkyl_update_status *st);

// Take a single time-step using a first-order operator split implicit fluid-EM coupling + SSP RK3
struct gkyl_update_status pkpm_update_op_split(gkyl_pkpm_app *app,
  double dt0);

// Take a fully explicit single time-step using a SSP RK3 (including explicit fluid-EM coupling)
struct gkyl_update_status pkpm_update_explicit_ssp_rk3(gkyl_pkpm_app *app,
  double dt0);

/** gkyl_pkpm_app private API */

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct pkpm_species* pkpm_find_species(const gkyl_pkpm_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of species, -1 if not found
 */
int pkpm_find_species_idx(const gkyl_pkpm_app *app, const char *nm);

/** pkpm_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app PKPM app object
 * @param s Species object 
 * @param sm Species moment object
 * @param nm Name string indicating moment type
 */
void pkpm_species_moment_init(struct gkyl_pkpm_app *app, struct pkpm_species *s,
  struct pkpm_species_moment *sm, bool is_diag);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void pkpm_species_moment_calc(const struct pkpm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release species moment object.
 *
 * @param app PKPM app object
 * @param sm Species moment object to release
 */
void pkpm_species_moment_release(const struct gkyl_pkpm_app *app,
  const struct pkpm_species_moment *sm);

/** pkpm_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app PKPM app object
 * @param s Species object 
 * @param lbo Species LBO object
 * @param collides_with_fluid Boolean for if kinetic species collides with a fluid species
 */
void pkpm_species_lbo_init(struct gkyl_pkpm_app *app, struct pkpm_species *s,
  struct pkpm_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app PKPM app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void pkpm_species_lbo_cross_init(struct gkyl_pkpm_app *app, struct pkpm_species *s,
  struct pkpm_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void pkpm_species_lbo_moms(gkyl_pkpm_app *app,
  const struct pkpm_species *species,
  struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param collides_with_fluid Boolean for if kinetic species collides with a fluid species
 * @param fluidin Input fluid array (size: num_fluid_species)
 */
void pkpm_species_lbo_cross_moms(gkyl_pkpm_app *app,
  const struct pkpm_species *species,
  struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute RHS from LBO collisions
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 * @return Maximum stable time-step
 */
void pkpm_species_lbo_rhs(gkyl_pkpm_app *app,
  const struct pkpm_species *species,
  struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species LBO object.
 *
 * @param app PKPM app object
 * @param sm Species LBO object to release
 */
void pkpm_species_lbo_release(const struct gkyl_pkpm_app *app, const struct pkpm_lbo_collisions *lbo);

/** pkpm_species API */

/**
 * Initialize species.
 *
 * @param pkpm Input PKPM data
 * @param app PKPM app object
 * @param s On output, initialized species object
 */
void pkpm_species_init(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app, struct pkpm_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app PKPM app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void pkpm_species_apply_ic(gkyl_pkpm_app *app, struct pkpm_species *species, double t0);

/**
 * Compute species applied acceleration term
 *
 * @param app PKPM app object
 * @param species Species object
 * @param tm Time for use in acceleration
 */
void pkpm_species_calc_app_accel(gkyl_pkpm_app *app, struct pkpm_species *species, double tm);

/**
 * Compute parallel-kinetic-perpendicular-moment (pkpm) model variables
 * These are the coupling moments [rho, p_par, p_perp], the self-consistent
 * pressure force (div(p_par b_hat)), and the primitive variables
 *
 * @param app PKPM app object
 * @param species Species object
 * @param fin Input distribution function
 * @param fluidin Input fluid species array (size: num_fluid_species)
 */
void pkpm_species_calc_pkpm_vars(gkyl_pkpm_app *app, struct pkpm_species *species, 
  const struct gkyl_array *fin, const struct gkyl_array *fluidin);

/**
 * Compute parallel-kinetic-perpendicular-moment (pkpm) model update variables
 * These are the acceleration variables in the kinetic equation and 
 * the source distribution functions for Laguerre couplings.
 *
 * @param app PKPM app object
 * @param species Species object
 * @param fin Input distribution function
 */
void pkpm_species_calc_pkpm_update_vars(gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin);

/**
 * Limit slopes of solution of fluid variables
 *
 * @param app PKPM app object
 * @param species Species object
 * @param fin Input distribution function 
 * @param fluid Input (and Output after limiting) array fluid species
 */
void pkpm_fluid_species_limiter(gkyl_pkpm_app *app, struct pkpm_species *species,
  struct gkyl_array *fin, struct gkyl_array *fluid);

/**
 * Compute RHS from species distribution function
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param em EM field
 * @param rhs On output, the RHS from the species object
 * @param fluidin Input fluid array for potential fluid force (size: num_fluid_species)
 * @return Maximum stable time-step
 */
double pkpm_species_rhs(gkyl_pkpm_app *app, struct pkpm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *em, 
  struct gkyl_array *rhs_f, struct gkyl_array *rhs_fluid);

/**
 * Apply BCs to species distribution functions 
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param f Distribution function to apply BCs to
 */
void pkpm_species_apply_bc(gkyl_pkpm_app *app, const struct pkpm_species *species, 
  struct gkyl_array *f);

/**
 * Apply BCs to species momentum
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param fluid momentum to apply BCs to
 */
void pkpm_fluid_species_apply_bc(gkyl_pkpm_app *app, const struct pkpm_species *species, 
  struct gkyl_array *fluid);

/**
 * Compute L2 norm (f^2) of the distribution function diagnostic
 *
 * @param app PKPM app object
 * @param tm Time at which diagnostic is computed
 * @param species Pointer to species
 */
void pkpm_species_calc_L2(gkyl_pkpm_app *app, double tm, const struct pkpm_species *species);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers
 */
void pkpm_species_coll_tm(gkyl_pkpm_app *app);

/**
 * Fill stat object in app with collisionless timers.
 *
 * @param app App object to update stat timers
 */
void pkpm_species_tm(gkyl_pkpm_app *app);

/**
 * Delete resources used in species.
 *
 * @param app PKPM app object
 * @param species Species object to delete
 */
void pkpm_species_release(const gkyl_pkpm_app* app, const struct pkpm_species *s);

/** pkpm_field API */

/**
 * Create new field object
 *
 * @param pkpm Input pkpm data
 * @param app PKPM app object
 * @return Newly created field
 */
struct pkpm_field* pkpm_field_new(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app);

/**
 * Compute field initial conditions.
 *
 * @param app PKPM app object
 * @param field Field object
 * @param t0 Time for use in ICs
 */
void pkpm_field_apply_ic(gkyl_pkpm_app *app, struct pkpm_field *field, double t0);

/**
 * Compute external electromagnetic fields
 *
 * @param app PKPM app object
 * @param field Field object
 * @param tm Time for use in external electromagnetic fields computation
 */
void pkpm_field_calc_ext_em(gkyl_pkpm_app *app, struct pkpm_field *field, double tm);

/**
 * Compute applied currents
 *
 * @param app PKPM app object
 * @param field Field object
 * @param tm Time for use in applied current computation
 */
void pkpm_field_calc_app_current(gkyl_pkpm_app *app, struct pkpm_field *field, double tm);

/**
 * Compute magnetic field unit vector and unit tensor
 *
 * @param app PKPM app object
 * @param field Field object (output bvar is stored in field object)
 * @param em Input electromagnetic fields
 */
void pkpm_field_calc_bvar(gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *em);

/**
 * Compute E x B velocity
 *
 * @param app PKPM app object
 * @param field Field object (output ExB is stored in field object)
 * @param em Input electromagnetic fields
 */
void pkpm_field_calc_ExB(gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *em);

/**
 * Accumulate current density onto RHS from field equations
 *
 * @param app PKPM app object
 * @param fluidin[] Input fluid array (num_species size)
 * @param emout On output, the RHS from the field solver *with* accumulated current density
 */
void pkpm_field_accumulate_current(gkyl_pkpm_app *app, 
  const struct gkyl_array *fluidin[], struct gkyl_array *emout);

/**
 * Limit slopes of solution of EM variables
 *
 * @param app PKPM app object
 * @param field Pointer to field 
 * @param em Input (and Output after limiting) EM fields
 */
void pkpm_field_limiter(gkyl_pkpm_app *app, struct pkpm_field *field, struct gkyl_array *em);

/**
 * Accumulate current density onto RHS from field equations in fully explicit update
 *
 * @param app PKPM app object
 * @param field Pointer to field 
 * @param fluidin[] Input fluid array (num_species size)
 * @param emout On output, the RHS from the field solver *with* accumulated current density
 */
void pkpm_field_explicit_accumulate_current(gkyl_pkpm_app *app, struct pkpm_field *field, 
  const struct gkyl_array *fluidin[], struct gkyl_array *emout);

/**
 * Compute RHS from field equations
 *
 * @param app PKPM app object
 * @param field Pointer to field
 * @param em Input field
 * @param rhs On output, the RHS from the field solver
 * @return Maximum stable time-step
 */
double pkpm_field_rhs(gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *em, struct gkyl_array *rhs);

/**
 * Apply BCs to field
 *
 * @param app PKPM app object
 * @param field Pointer to field
 * @param f Field to apply BCs
 */
void pkpm_field_apply_bc(gkyl_pkpm_app *app, const struct pkpm_field *field,
  struct gkyl_array *f);

/**
 * Compute field energy diagnostic
 *
 * @param app PKPM app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 */
void pkpm_field_calc_energy(gkyl_pkpm_app *app, double tm, const struct pkpm_field *field);

/**
 * Release resources allocated by field
 *
 * @param app PKPM app object
 * @param f Field object to release
 */
void pkpm_field_release(const gkyl_pkpm_app* app, struct pkpm_field *f);

/**
 * Create new fluid-EM coupling updater for the PKPM system
 *
 * @param app PKPM app object
 * @return Newly created fluid-EM coupling updater
 */
struct pkpm_fluid_em_coupling* pkpm_fluid_em_coupling_init(struct gkyl_pkpm_app *app);

/**
 * Compute implicit update of fluid-EM coupling for the PKPM system
 *
 * @param app PKPM app object
 * @param pkpm_em fluid-EM coupling updater
 * @param tcurr Current time
 * @param dt Time step size
 */
void pkpm_fluid_em_coupling_update(struct gkyl_pkpm_app *app, 
  struct pkpm_fluid_em_coupling *pkpm_em, double tcurr, double dt);

/**
 * Release resources allocated by fluid-EM coupling object for the PKPM system
 *
 * @param app PKPM app object
 * @param pkpm_em fluid-EM coupling updater to release
 */
void pkpm_fluid_em_coupling_release(struct gkyl_pkpm_app *app, 
  struct pkpm_fluid_em_coupling *pkpm_em);
