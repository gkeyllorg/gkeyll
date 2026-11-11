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
#include <gkyl_bc_pkpm_species_reflect.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_pkpm_em_coupling.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_flux_surf.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_pkpm_mom_diffusion.h>
#include <gkyl_dg_updater_lbo_pkpm.h>
#include <gkyl_dg_updater_moment_pkpm.h>
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
#include <gkyl_pkpm_bgk_collisions.h>
#include <gkyl_pkpm_coll_nu.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_dg_updater_lbo_pkpm_nodal.h>
#include <gkyl_pkpm_lbo_moms.h>
#include <gkyl_pkpm_lte.h>
#include <gkyl_pkpm_temp_source.h>
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
  int num_cdim, num_vdim; // phase-space split of a hybrid basis (written when num_vdim > 0)
  // Staggered operator split (0: none): the step of the last update and the step it suggested, so a
  // restart continues as the uninterrupted run would.
  double dt_prev, dt_suggested;
};

// data for moments
struct pkpm_species_moment {
  struct gkyl_dg_updater_moment *mcalc; // moment update

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

// forward declare species struct
struct pkpm_species;

// Context of the host-side projection of the source distributions: the p=1 configuration
// fields are evaluated at the computational coordinates of the phase-space quadrature points.
struct pkpm_source_eval_ctx {
  int cdim;
  double mass;
  const struct gkyl_rect_grid *conf_grid;
  const struct gkyl_basis *conf_basis;
  const struct gkyl_range *local_ext;
  const struct gkyl_pkpm_source *info;
  const struct gkyl_array *prof; // [n_s, T_s, u_s] of the source (physical)
  const struct gkyl_array *rho, *u, *xi, *bvar; // nodal state of the stage (host)
};

struct pkpm_source {
  enum gkyl_source_id source_id;
  bool calc_bflux; // rate from the wall losses of source_species
  double source_length;
  struct pkpm_species *source_species;
  double scale; // rate factor of the stage (1 for a projected source)
  struct gkyl_array *prof_host; // [n_s, T_s, u_s] (host)
  struct gkyl_array *fluid_src,
    *fluid_src_host; // J [rho, M, E] of the profile; physical copy for IO
  struct gkyl_array *f_src, *f_src_host; // J [S_F, S_G] of the stage
  struct gkyl_array *n_src; // m n_s (physical)
  struct gkyl_array *n_ratio; // m n_s / rho
  struct gkyl_array *rho_host, *u_host, *xi_host, *bvar_host; // host mirrors (GPU only)
  gkyl_dg_bin_op_mem *div_mem;
  gkyl_proj_on_basis *proj_f;
  struct pkpm_source_eval_ctx eval_ctx;
};

struct pkpm_bgk_collisions;

// Accumulate the kinetic BGK term of a species onto rhs.
typedef void (*pkpm_bgk_kinetic_rhs_t)(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

struct pkpm_bgk_collisions {
  // The operator at the configuration-space nodes.
  bool use_nodal; // false for GPU runs
  bool is_implicit; // backward Euler after each SSP-RK3 step (nodal only)
  struct gkyl_array *rate; // relaxation rate at the nodes
  struct gkyl_array *shape; // implicit shape factor S at the nodes
  struct gkyl_array *f_lte; // corrected unit Maxwellian in w, 3 modes per velocity cell
  gkyl_pkpm_bgk_collisions *up_bgk; // updater
  pkpm_bgk_kinetic_rhs_t kinetic_rhs; // the nodal or modal kinetic term

  // The modal form of the operator (GPU runs).
  struct gkyl_array *nu; // collision frequency
  struct gkyl_array *f_target; // J [F_t, G_t]
};

struct pkpm_source_bgk;

// Accumulate the temperature source of a species onto the kinetic and fluid right-hand sides.
typedef void (*pkpm_source_bgk_rhs_t)(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
);

// BGK-form temperature source at the configuration-space nodes.
struct pkpm_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // type of source
  bool is_implicit; // backward Euler after each SSP-RK3 step
  struct gkyl_array *rate; // source rate nu_Q at the nodes
  struct gkyl_array *vtsq; // source temperature over the mass, T_Q/m, at the nodes
  struct gkyl_array *shape; // rate (explicit) or shape factor (implicit) of the distributions
  struct gkyl_array *f_lte; // corrected unit Maxwellian in w, 3 modes per velocity cell
  gkyl_pkpm_temp_source *up_temp; // thermal-energy relaxation of the fluid state
  gkyl_pkpm_bgk_collisions *up_bgk; // relaxation of the distributions
  pkpm_source_bgk_rhs_t rhs; // the explicit terms, or nothing
};

struct pkpm_lbo_collisions;

// Compute the primitive moments and the collision frequency of the self LBO of a species.
typedef void (*pkpm_lbo_moms_t)(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

// Accumulate the self LBO of a species onto rhs.
typedef void (*pkpm_lbo_kinetic_rhs_t)(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

struct pkpm_lbo_collisions {
  // The operator at the configuration-space nodes.
  bool use_nodal; // false for GPU runs
  struct gkyl_array *nu_ref; // reference collision frequency at the nodes
  struct gkyl_array *nu; // collision frequency at the nodes
  struct gkyl_array *n_nodal; // number density at the nodes (Spitzer)
  struct gkyl_array *xi_nodal; // [xi, xi^2, p] at the nodes (Spitzer)
  struct gkyl_array *vtsq_nodal; // physical thermal speed squared at the nodes (Spitzer)
  gkyl_pkpm_coll_nu *up_nu; // collision frequency at the nodes
  gkyl_pkpm_lbo_moms *up_moms; // primitive moments at the nodes
  gkyl_dg_updater_lbo_pkpm_nodal *up_coll; // the operator at the nodes
  pkpm_lbo_moms_t calc_moms; // the nodal (constant or Spitzer frequency) or modal moments
  pkpm_lbo_kinetic_rhs_t kinetic_rhs; // the nodal or modal operator

  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *nu_prim_moms; // expansions of nu and nu [u, vtsq] (GPU runs)
  struct gkyl_array *prim_moms; // [u, vtsq]: node values, or expansions on GPU runs
  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array
    *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq
    *spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value
  struct gkyl_array *maxwellian_moms; // Maxwellian moments n, t, v_t^2.

  double betaGreenep1; // value of Greene's factor beta + 1
  double other_m[GKYL_MAX_SPECIES]; // masses of species being collided with
  struct gkyl_array
    *other_prim_moms[GKYL_MAX_SPECIES]; // self-primitive moments of species being collided with
  struct gkyl_array *cross_prim_moms[GKYL_MAX_SPECIES]; // LBO cross-primitive moments
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // LBO cross-species collision frequencies
  struct gkyl_array *other_nu[GKYL_MAX_SPECIES];
  struct gkyl_array
    *cross_nu_prim_moms; // weak multiplication of collision frequency and primitive moments

  struct gkyl_array *self_nu, *self_nu_prim_moms; // LBO self-primitive moments

  struct pkpm_species_moment
    moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)
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

// Blocks of the evolved normalized distributions, one phase-space expansion each:
// [F, G, M0_perp (x, y, z), N_perp (x, y, z)]. M0_perp and N_perp are the first gyroangle
// harmonic as Cartesian vectors perpendicular to b (normalized_first_harmonic.tex).
enum pkpm_dist_block {
  PKPM_DIST_F = 0,
  PKPM_DIST_G = 1,
  PKPM_DIST_M0 = 2, // first of the three Cartesian components of M0_perp
  PKPM_DIST_N = 5, // first of the three Cartesian components of N_perp
  PKPM_DIST_NUM_BLOCKS = 8
};

// Records of the fluid positivity limiter (gkyl_dg_calc_pkpm_vars_positivity): the per-cell
// count [cells scaled for the density, for the internal energy, for the thermal speed,
// inadmissible averages], and the statistics per integrated-moment interval, which append the
// number of limiter applications.
enum pkpm_positivity_record { PKPM_POS_NUM_COUNT = 4, PKPM_POS_NUM_STAT = 5 };

// species data
struct pkpm_species {
  struct gkyl_pkpm_species info; // data for species

  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local,
    local_ext; // local, local-ext phase-space ranges (local_ext has one velocity ghost cell above w_max)
  struct gkyl_range
    local_ext_conf; // local range extended in configuration space only (moments, distribution BCs, sync)
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges

  struct gkyl_comm *comm; // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1,
    *fnew; // normalized distributions [F, G, M0_perp, N_perp] on the w = v_par/xi grid
  struct gkyl_array *blk_in,
    *blk_out; // one-block (F or G) scratch for configuration x phase products
  struct gkyl_array
    *f_no_J; // physical distributions (position-map Jacobian divided out; output and L2)
  struct gkyl_array *f_no_J_host; // host copy of f_no_J for I/O
  struct gkyl_array *fluid, *fluid1, *fluidnew; // fluid state [rho, rho ux, rho uy, rho uz, E]

  struct gkyl_array *cflrate_f; // CFL rate in each cell for distribution function update
  // Operators that act in the velocity alone at each configuration-space node (collisions, the
  // temperature source) share one conversion of the distributions per stage and one conversion
  // of their summed increments back.
  bool has_nodal_ops; // true if the species has such operators (CPU runs)
  struct gkyl_array *f_nodal; // distributions in the configuration-nodal, velocity-modal form
  struct gkyl_array *rhs_nodal; // increments of the nodal operators, in the same form
  struct gkyl_array *cflrate_fluid; // CFL rate in each cell for momentum update

  struct gkyl_array *bc_buffer_dist; // buffer for BCs for distribution functions (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed_dist,
    *bc_buffer_up_fixed_dist; // fixed buffers for time independent BCs for distribution functions
  struct gkyl_array *bc_buffer_fluid; // buffer for BCs for momentum (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed_fluid,
    *bc_buffer_up_fixed_fluid; // fixed buffers for time independent BCs for momentum

  struct gkyl_array *f_host; // host copy of distribution function for use IO and initialization
  struct gkyl_array *fluid_host; // host copy of momentum for use IO and initialization

  // Duplicate copy of fluid data in case time step fails.
  // Needed because of implicit source split which modifies solution and
  // is always successful, so if a time step fails due to the SSP RK3
  // we must restore the old solution before restarting the time step
  struct gkyl_array *fluid_dup;

  struct gkyl_wv_eqn
    *equation; // For storing 10 moment equation object for upwinding fluid equations with Roe solve

  struct gkyl_array *qmem; // array for q/m*(E,B) for use in *explicit* update
    // Note: this array is *only* used if the PKPM self-consistent EM fields are static
    // If PKPM self-consistent EM fields are dynamics we utilize an implicit source update
    // for the momentum equations and Ampere's Law

  struct pkpm_species_moment pkpm_moms_diag; // for computing pkpm moments diagnostics
  struct pkpm_species_moment
    pkpm_nmoms; // raw moments [M0, M2, M0_G, M1, M3, M1_G] of the normalized [F, G] (mass 1)
  struct gkyl_array *
    pkpm_moms; // kinetic coupling moments [p_par, p_perp, q_par, q_perp, s (3), r (3)] = [p M2, p M0_G, p xi M3, p xi M1_G, p m1, p xi (m2/2 + n0)]
  struct gkyl_array *
    pkpm_cmoms; // the constraints [M0(F), M1(F), (M2(F) + 2 M0(G))/3, int M0_perp dw (x, y, z)], [1, 0, 1, 0, 0, 0] when the closures hold
  struct gkyl_array *pkpm_cmoms_host; // host copy for output

  // Refactored (p=1) PKPM update variables.
  struct gkyl_array
    *pkpm_u; // flow velocity [ux, uy, uz], nodal quotient M/rho (the weak division at p=1)
  struct gkyl_array *
    pkpm_grad_u; // gradient of u for the sources: in-cell plus the blended face lift, cdim x 3 expansions
  struct gkyl_array *pkpm_theta; // T_perp/T = M0(G)
  struct gkyl_array
    *pkpm_xi; // thermal speed and scalar pressure [xi, xi^2, p] of the fluid state (nodal)
  struct gkyl_array
    *pkpm_xi3; // copy of pkpm_xi for the weak-gradient kernels (only the first block is used)
  struct gkyl_array
    *pkpm_grad_xi; // gradient of [xi, ...] (pkpm_grad_u layout; blocks (d, 0) hold d_d xi)
  struct gkyl_array *pkpm_scratch; // one-expansion scratch for weak products
  struct gkyl_array *
    fluid_no_J; // PHYSICAL fluid state [rho, M, E] of the stage (position-map Jacobian divided out)
  struct gkyl_array *pkpm_moms_no_J; // physical coupling moments (diagnostics)
  struct gkyl_array *pkpm_rho; // physical mass density (leading block of fluid_no_J)
  struct gkyl_array *pkpm_rho3; // [rho, 0, 0] for the weak-gradient kernels
  struct gkyl_array
    *pkpm_mom3; // momentum block of the fluid array (scratch for the 3-component consumers)
  struct gkyl_array *pkpm_grad_rho; // gradient of [rho, 0, 0] (scratch)
  struct gkyl_array *pkpm_grad_ln_n; // grad(ln n), cdim expansions
  struct gkyl_array *pkpm_c_0; // normalization closure of the factored update
  struct gkyl_array *pkpm_alpha_0; // pressure force a_P/xi closing the first w moment
  struct gkyl_array *pkpm_c_h; // heating rate closing the scalar second moment M2(F) + 2 M0(G)
  struct gkyl_array *pkpm_a_perp; // perpendicular force closing int M0_perp dw [x, y, z]
  struct gkyl_array
    *pkpm_rhs_moms; // moments of the configuration-pass right-hand side (separate scratch)
  struct gkyl_array *pkpm_p_ij; // p=1 projection of the pressure tensor (limiter scratch only)
  struct gkyl_array
    *conf_flux_surf; // stored modal configuration-space kinetic fluxes, [dir][side][block][mode]
  struct gkyl_array *flux_blend; // weight of the label-paired flux of each lower face, [dir]
  struct gkyl_array *vel_flux_surf; // stored modal v_par kinetic fluxes
  struct gkyl_dg_pkpm_harm_vol
    *harm_vol; // modal volume integrals of the first-harmonic terms (stage A)
  struct gkyl_array *harm_coef; // stored operator of the stage (gkyl_pkpm_harm_ops.h)
  struct gkyl_array *harm_wnode; // raw w-flux values of the stage at the w-face nodes
  struct gkyl_dg_pkpm_harm_exp
    *harm_exp; // gyrophase-mixing integrating factor of the first harmonic (plan 17.44)
  struct gkyl_array *harm_exp_coef; // its field (q/m) B frozen at the start of the step
  struct gkyl_array *harm_exp_nu; // gyrophase-mixing rate at the phase nodes, frozen per step
  struct gkyl_array *harm_exp_act; // configuration cells with a damped node (momentum restoration)
  struct gkyl_array *harm_exp_col; // w integrals of F and M0_perp of the restoration
  struct gkyl_array *fluid_flux_surf; // stored modal fluid fluxes [F_rho, F_x, F_y, F_z, F_E]
  struct gkyl_array *
    integ_pkpm_mom; // integrated PKPM variables [rho, rho ux, rho uy, rho uz, rho ux^2, rho uy^2, rho uz^2, p_par, p_perp]

  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars; // Updater for the configuration-space variables
  struct gkyl_dg_pkpm_remap_flux_surf
    *calc_face_flux; // configuration-space face flux updater (kinetic and fluid)
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vel_flux; // v_par flux updater
  struct gkyl_dg_eqn *eqn_kinetic_conf; // configuration pass (stage A) equation object
  struct gkyl_dg_eqn *eqn_kinetic_vel; // velocity pass (stage B) equation object
  struct gkyl_dg_eqn *eqn_fluid; // momentum equation object
  gkyl_hyper_dg *slvr_kinetic_conf; // hyper_dg over the configuration directions
  gkyl_hyper_dg *slvr_kinetic_vel; // hyper_dg over v_par
  gkyl_hyper_dg *slvr_fluid; // hyper_dg for the momentum equation
  double kinetic_tm; // time spent in the kinetic right-hand side
  double fluid_tm; // time spent in the momentum right-hand side
  long num_null_rho; // accumulated count of nodes with an inadmissible density in the closure
  long
    num_null_closure; // accumulated count of closure nodes nulled (M0 <= 0 or a non-positive determinant)
  double closure_cond_min; // smallest det/det0 of the closure solves over the run
  bool closure_null_warned; // the velocity-extent warning was printed
  struct gkyl_array
    *pkpm_pos_count; // per-cell record of the positivity limiter [rho, e, xi, avg] (cumulative)
  struct gkyl_array *pkpm_pos_io,
    *pkpm_pos_io_host; // the record as cell averages of a p=1 field for IO
  long pos_count_since
    [PKPM_POS_NUM_STAT]; // limiter counts [rho, e, xi, avg, applications] since the last integrated-moment call
  gkyl_dynvec pos_stat; // the counts per integrated-moment interval, summed over ranks
  bool is_first_pos_stat_write_call; // flag for the positivity statistics dynvec written first time

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

  // boundary conditions on lower/upper edges in each direction
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC on distribution functions and momentum.
  struct gkyl_bc_basic *bc_lo_dist[3];
  struct gkyl_bc_basic *bc_up_dist[3];
  // reflecting walls of the distributions (0 at other walls): the mirror image in the ghost
  // cell's field-aligned coordinates, from b in the skin and ghost cells
  struct gkyl_bc_pkpm_species_reflect *bc_lo_reflect[3];
  struct gkyl_bc_pkpm_species_reflect *bc_up_reflect[3];
  struct gkyl_bc_basic *bc_lo_fluid[3];
  struct gkyl_bc_basic *bc_up_fluid[3];
  // absorbing walls (0 at other walls): theta mirrored (even image) into the ghost cells, whose
  // empty G has no T_perp/T of its own, so that its face averages are the interior wall traces
  struct gkyl_bc_basic *bc_lo_theta[3], *bc_up_theta[3];
  // To simplify BC application, store local skin and ghost ranges for distribution functions
  struct gkyl_range lower_skin_dist[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost_dist[GKYL_MAX_DIM];
  struct gkyl_range upper_skin_dist[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost_dist[GKYL_MAX_DIM];
  bool bc_is_absorb; // whether any boundary is absorbing (wall fluxes from the kinetic upwind flux)
    // If BCs are *not* absorbing, primitive variables can be calculated on *extended* range

  bool has_app_accel; // flag to indicate there is applied acceleration
  bool app_accel_evolve; // flag to indicate applied acceleration is time-dependent
  struct gkyl_array *app_accel; // applied acceleration
  struct gkyl_array *app_accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_accel_proj; // projector for acceleration

  enum gkyl_collision_id collision_id; // type of collisions
  struct pkpm_lbo_collisions lbo; // collisions object
  struct pkpm_bgk_collisions bgk; // BGK collisions object

  struct pkpm_source src; // source object
  struct pkpm_source_bgk src_bgk; // temperature source object
  double bflux
    [2 *
     GKYL_MAX_CDIM]; // particle loss rate through the absorbing walls [lower, upper] per direction (stage)
  gkyl_dynvec bflux_diag; // the loss rates on the cadence of the integrated moments
  bool is_first_bflux_write_call;
  struct gkyl_array *fluid_flux_surf_host; // host mirror of the stored fluid fluxes (GPU only)

  // fluid diffusion
  bool has_diffusion; // flag to indicate there is applied momentum diffusion
  struct gkyl_dg_pkpm_mom_diffusion *mom_diff; // momentum diffusion on the position map

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

  struct gkyl_array *em_host; // host copy for use IO and initialization

  // Duplicate copy of EM data in case time step fails.
  // Needed because of implicit source split which modifies solution and
  // is always successful, so if a time step fails due to the SSP RK3
  // we must restore the old solution before restarting the time step
  struct gkyl_array *em_dup;

  // Identity position map (PKPM runs on uniform grids); the Maxwell curl

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
  struct gkyl_array
    *app_current_no_J; // physical applied current (the implicit coupling); app_current carries J
  struct gkyl_array
    *em_no_J; // physical E, B of the stage (the evolved em carries the position-map Jacobian)
  struct gkyl_array *em_no_J_host; // host copy of em_no_J for I/O
  struct gkyl_vlasov_position_map *pos_map; // acquired reference to the app's position map
  gkyl_proj_on_basis *app_current_proj; // projector for applied current

  bool use_ghost_current; // 1D periodic: subtract the domain-mean x current from Ampere's law
  struct gkyl_array *ghost_current; // one component: the mean current (coefficient 0), constant
  double *red_ghost_current; // device reduction scratch
  struct gkyl_array *app_current_tot; // applied current plus the ghost current (implicit coupling)

  // Operator split (plan 17.29): Stormer-Verlet in the curl terms. E is advanced only by the
  // implicit source half steps, with curl B (at the half step's start, known) as an effective
  // current; B drifts over the transport step with the E of the first half kick, so the field
  // right-hand side of the transport step is that of RK stage 1 with its E rows zeroed.
  struct gkyl_array *rhs_split; // field right-hand side of the transport step (stage 1, E rows zero)
  double omega_cfl_split; // its time-step rate
  struct gkyl_array *curl_rhs,
    *curl_rhs_no_J; // field right-hand side at a source half step (E rows: c^2 curl B)
  struct gkyl_array
    *curl_current; // applied current - eps0 c^2 curl B (physical), for the implicit solve

  struct gkyl_array *bvar; // magnetic field unit vector b = B/|B| [bx, by, bz] (p=1)
  struct gkyl_array *grad_b; // weak (surface-averaged) gradient of b, cdim x 3 expansions
  // physical time derivative [E, B, phi, psi] of the stage fields, the Maxwell right-hand side
  // divided by the position-map Jacobian (computed first in each stage; zero for a static field;
  // its E part lacks the explicit current, which is added after the species update)
  struct gkyl_array *dt_em_no_J;
  // d_t B/|B| of the stage at the Gauss-Legendre nodes of b, interpolated to p=1: the force
  // kernels of the first harmonic form d_t b = (I - b b) . (d_t B/|B|) from it
  struct gkyl_array *dtB_over_magB;
  struct gkyl_dg_calc_pkpm_bvar *calc_bvar; // Updater computing b by nodal normalization

  bool limit_em; // boolean for whether or not we are limiting EM fields
  struct gkyl_array *em_lim; // physical E, B being limited (scratch; only with limit_em)
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

  double *omegaCfl_ptr;
};

// fluid-EM coupling data in PKPM system
struct pkpm_fluid_em_coupling {
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each species
  struct gkyl_dg_calc_pkpm_em_coupling *slvr; // fluid-EM coupling solver
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

  double op_split_dt_prev; // step taken by the previous operator-split update (0: none yet)
  double op_split_dt_suggested; // step suggested by the previous operator-split update (0: none yet)
  // Staggered split: between updates the source-updated state (fluid, E) is half a source step
  // behind tcurr. Output and diagnostics advance a copy of it by that half step (synchronized), and
  // the next update restores the staggered state.
  bool op_split_synced; // the fluid and field hold the synchronized state
  struct gkyl_array *op_split_fluid[GKYL_MAX_SPECIES]; // staggered fluid states while synchronized
  struct gkyl_array *op_split_em; // staggered field while synchronized
  // per-cell limiter record held across the synchronizing source step
  struct gkyl_array *op_split_pos_count[GKYL_MAX_SPECIES];
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
  struct gkyl_comm *comm; // communicator object for conf-space arrays

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom
    *geom; // geometry needed for species and field solvers (*only* p=1 right now JJ: 05/03/24)
  struct gkyl_vlasov_position_map *
    pos_map; // configuration-space position map (always created; identity by default); shared by the field and species
  struct gkyl_dg_conf_nodal
    *conf_nodal; // conversions to and from the configuration-nodal, velocity-modal form (CPU runs)

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
  double fluid_positivity_eps; // Relative floor of the vertex positivity limiter (always applied)
  int num_face_nodes; // Transverse Gauss-Legendre nodes of a configuration face, per direction
  struct pkpm_fluid_em_coupling *pkpm_em; // fluid-EM coupling data
  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_pkpm_app *app, double dt0);

  struct gkyl_pkpm_stat stat; // statistics
};

// Take a single forward Euler step of the PKPM system with the suggested time-step dt. The
// first-harmonic blocks of the stage state fin are projected onto the plane perpendicular to the
// stage's b first (in place). stage (1, 2 or 3) is the SSP-RK3 stage: where the gyrophase mixing
// damps, the first-harmonic blocks of fout are those of its integrating-factor stage in the Shu-Osher
// form of the driver.
void pkpm_forward_euler(
  gkyl_pkpm_app *app, int stage, double tcurr, double dt, struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *fout[],
  struct gkyl_array *fluidout[], struct gkyl_array *emout, struct gkyl_update_status *st
);

// Take a single time-step of the staggered operator split: one implicit fluid-EM source step over the
// interval between the midpoints of consecutive transport steps, then the SSP-RK3 transport step.
struct gkyl_update_status pkpm_update_op_split(gkyl_pkpm_app *app, double dt0);

// Staggered split: advance the source-updated state (fluid, E) by the half source step it lags tcurr,
// keeping the staggered state for the next update (no-op if already synchronized or not staggered).
void pkpm_op_split_sync(gkyl_pkpm_app *app);

// Staggered split: restore the staggered state saved by pkpm_op_split_sync (no-op if not synchronized).
void pkpm_op_split_unsync(gkyl_pkpm_app *app);

// Take a fully explicit single time-step using a SSP RK3 (including explicit fluid-EM coupling)
struct gkyl_update_status pkpm_update_explicit_ssp_rk3(gkyl_pkpm_app *app, double dt0);

// The SSP-RK3 step followed by one backward-Euler step of the implicit collisions over the step
// RK3 took (the transport step T of both update paths).
struct gkyl_update_status pkpm_update_ssp_rk3_implicit_coll(gkyl_pkpm_app *app, double dt0);

/** gkyl_pkpm_app private API */

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct pkpm_species *pkpm_find_species(const gkyl_pkpm_app *app, const char *nm);

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
void pkpm_species_moment_diag_calc(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin, const struct gkyl_array *fluid
);
void pkpm_species_moment_init_mass(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_species_moment *sm, bool is_diag,
  double mass
);
void pkpm_species_moment_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_species_moment *sm, bool is_diag
);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void pkpm_species_moment_calc(
  const struct pkpm_species_moment *sm, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin
);

/**
 * Calculate moments of an arbitrary phase-space array into a separate output array.
 *
 * @param sm Species moment object
 * @param phase_rng Phase-space range
 * @param conf_rng Configuration-space range
 * @param fin Input phase-space array
 * @param mout Output moment array
 */
void pkpm_species_moment_calc_into(
  const struct pkpm_species_moment *sm, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin, struct gkyl_array *mout
);

/**
 * Release species moment object.
 *
 * @param app PKPM app object
 * @param sm Species moment object to release
 */
void pkpm_species_moment_release(
  const struct gkyl_pkpm_app *app, const struct pkpm_species_moment *sm
);

/** pkpm_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app PKPM app object
 * @param s Species object 
 * @param lbo Species LBO object
 * @param collides_with_fluid Boolean for if kinetic species collides with a fluid species
 */
void pkpm_species_lbo_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app PKPM app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void pkpm_species_lbo_cross_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void pkpm_species_lbo_moms(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

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
void pkpm_species_lbo_cross_moms(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

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
void pkpm_species_lbo_rhs(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

/**
 * Release species LBO object.
 *
 * @param app PKPM app object
 * @param sm Species LBO object to release
 */
void pkpm_species_lbo_release(
  const struct gkyl_pkpm_app *app, const struct pkpm_lbo_collisions *lbo
);

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
/**
 * Convert the distributions of the stage to the configuration-nodal, velocity-modal form and
 * clear the increments of the nodal operators (species with such operators only).
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param fin Input distribution functions
 */
void pkpm_species_nodal_dist(
  gkyl_pkpm_app *app, const struct pkpm_species *species, const struct gkyl_array *fin
);

/**
 * Accumulate the increments of the nodal operators onto modal distributions (species with
 * such operators only).
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param out Modal distributions to increment
 */
void pkpm_species_nodal_accumulate(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *out
);

/**
 * Computational-to-physical coordinate transform of the position map, for projecting
 * initial conditions and sources on a non-uniform configuration grid.
 *
 * @param xcomp Computational coordinates (the leading cdim entries are mapped in place)
 * @param xphys Physical coordinates on output
 * @param ctx Position map (const struct gkyl_vlasov_position_map *)
 */
void pkpm_pos_c2p(const double *xcomp, double *xphys, void *ctx);

void pkpm_species_calc_pkpm_vars(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Configuration-space face fluxes of the stage, the kinetic and the fluid flux and both rates in
 * one pass (the species' rates must be cleared by the caller), then the frame gradients of the
 * volume terms, which need the blend of the faces.
 */
void pkpm_species_face_flux(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin
);

/** Particle loss rates through the absorbing walls from the stored fluxes (into species->bflux). */
void pkpm_species_calc_bflux(gkyl_pkpm_app *app, struct pkpm_species *s);

// Source object: see gkyl_pkpm_source in gkyl_pkpm.h.
void pkpm_species_source_init(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source *src);
void pkpm_species_source_calc_rate(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source *src
);
void pkpm_species_source_kinetic_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source *src, const struct gkyl_array *fin,
  struct gkyl_array *rhs
);
void pkpm_species_source_fluid_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source *src, struct gkyl_array *rhs
);
void pkpm_species_source_write(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source *src, double tm, int frame
);
void pkpm_species_source_release(const gkyl_pkpm_app *app, const struct pkpm_source *src);

// Project a configuration-space profile at the physical coordinates of the position map and
// scale it. Returns a new host array on the extended range: the values at the 2^cdim interior
// Gauss-Legendre nodes if is_nodal is set, the p=1 expansion otherwise.
struct gkyl_array *pkpm_species_coll_profile_new(
  gkyl_pkpm_app *app, void (*eval)(double t, const double *xn, double *fout, void *ctx), void *ctx,
  double scale, bool is_nodal
);
// The target of the BGK-form operators: a new host array with the unit Maxwellian in w on the
// velocity grid of the species (3 modes per cell), its discrete moments corrected to [1, 0, 1].
struct gkyl_array *pkpm_species_bgk_target_new(gkyl_pkpm_app *app, struct pkpm_species *s);

// BGK collisions: nu (f_target - f) with the normalized unit Maxwellian target. They conserve
// the density, momentum and energy, so they never touch the fluid.
void pkpm_species_bgk_init(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk
);

// Temperature source (pkpm_species_source_bgk.c). The explicit form enters the right-hand
// sides before the closures; the implicit form is one backward-Euler step over dt on the
// distributions and the fluid state of the species (interior only).
void pkpm_species_source_bgk_init(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src
);
void pkpm_species_source_bgk_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
);
void pkpm_species_source_bgk_implicit_advance(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src, double dt
);
void pkpm_species_source_bgk_release(const gkyl_pkpm_app *app, const struct pkpm_source_bgk *src);
// Project M0_perp and N_perp of the distributions f onto the plane perpendicular to the current b
// at the Gauss-Legendre nodes of b, in place, over the extended ranges (plan 17.11).
void pkpm_species_project_transverse(
  gkyl_pkpm_app *app, struct pkpm_species *species, struct gkyl_array *f
);
// Gyrophase-mixing integrating factor of the first-harmonic blocks (plan 17.44, 17.24), after the
// forward Euler accumulate of stage s (fout = fin + h rhs; the gyration is explicit in rhs). At
// s = 1 it freezes the mixing of the step (the field and the rate at the phase nodes). At damped
// nodes it writes the first-harmonic blocks of fout as the value whose SSP-RK3 combination with the
// state at the start of the step is the integrating-factor stage, and restores int M0_perp dw = 0
// where the mixing acted.
void pkpm_species_harm_exp_advance(
  gkyl_pkpm_app *app, struct pkpm_species *species, int stage, double h, struct gkyl_array *fout
);

// out += a conf * phase for every block of a phase array (the bin-op kernels act on one
// block of num_basis components); phase may hold fewer blocks than out (e.g. [F, G] only).
void pkpm_species_accumulate_conf_phase(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct gkyl_array *out, double a,
  const struct gkyl_array *conf, const struct gkyl_array *phase
);

void pkpm_species_bgk_kinetic_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);
// One backward-Euler step of the BGK operator over dt on the species distribution (interior only).
void pkpm_species_bgk_implicit_advance(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double dt
);
void pkpm_species_bgk_release(const gkyl_pkpm_app *app, const struct pkpm_bgk_collisions *bgk);

// Output metadata helpers (pkpm.c).
struct gkyl_msgpack_data *pkpm_array_meta_new(struct pkpm_output_meta meta);
void pkpm_array_meta_release(struct gkyl_msgpack_data *mt);

/**
 * Apply the vertex positivity limiter to a fluid state over the local range (no-op when
 * disabled), recording the cells it acted on and the run totals.
 *
 * @param app PKPM app object
 * @param species Species object
 * @param fluid Fluid state [rho, M, E] to limit in place
 */
void pkpm_species_apply_positivity(
  gkyl_pkpm_app *app, struct pkpm_species *species, struct gkyl_array *fluid
);

/**
 * Compute the diagnostic variables that need the configuration pass (the pressure
 * force a_P), for output. Uses the species scratch arrays.
 *
 * @param app PKPM app object
 * @param species Species object
 * @param fin Input distribution function
 * @param fluidin Input fluid species array
 */
void pkpm_species_calc_diag_vars(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

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
double pkpm_species_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, const struct gkyl_array *em, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
);

/**
 * Apply BCs to species distribution functions 
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param f Distribution function to apply BCs to
 */
void pkpm_species_apply_bc(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *f
);

/**
 * Apply BCs to species momentum
 *
 * @param app PKPM app object
 * @param species Pointer to species
 * @param fluid momentum to apply BCs to
 */
void pkpm_fluid_species_apply_bc(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *fluid
);

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
void pkpm_species_release(const gkyl_pkpm_app *app, const struct pkpm_species *s);

/** pkpm_field API */

/**
 * Create new field object
 *
 * @param pkpm Input pkpm data
 * @param app PKPM app object
 * @return Newly created field
 */
struct pkpm_field *pkpm_field_new(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app);

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

// Domain mean of the species x current (coefficient 0 of sum_s q_s M_x/m_s), stored as a
// constant (physical) in field->ghost_current. fluidin are the J-weighted (evolved) fluid states:
// their cell mean is the physical domain mean on a position map that fixes the end points.
void pkpm_field_calc_ghost_current(
  gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *fluidin[]
);

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
void pkpm_field_accumulate_current(
  gkyl_pkpm_app *app, const struct gkyl_array *fluidin[], struct gkyl_array *emout
);

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
void pkpm_field_explicit_accumulate_current(
  gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *fluidin[],
  struct gkyl_array *emout
);

/**
 * Compute RHS from field equations
 *
 * @param app PKPM app object
 * @param field Pointer to field
 * @param stage SSP-RK3 stage (1, 2 or 3). With the operator split the right-hand side of stage 1
 *   (E rows zeroed: B drifts with the E of the source half kick) is kept for stages 2 and 3.
 * @param em Input field
 * @param rhs On output, the RHS from the field solver
 * @return Maximum stable time-step
 */
double pkpm_field_rhs(
  gkyl_pkpm_app *app, struct pkpm_field *field, int stage, const struct gkyl_array *em,
  struct gkyl_array *rhs
);

/**
 * Apply BCs to field
 *
 * @param app PKPM app object
 * @param field Pointer to field
 * @param f Field to apply BCs
 */
void pkpm_field_apply_bc(gkyl_pkpm_app *app, const struct pkpm_field *field, struct gkyl_array *f);

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
void pkpm_field_release(const gkyl_pkpm_app *app, struct pkpm_field *f);

/**
 * Create new fluid-EM coupling updater for the PKPM system
 *
 * @param app PKPM app object
 * @return Newly created fluid-EM coupling updater
 */
struct pkpm_fluid_em_coupling *pkpm_fluid_em_coupling_init(struct gkyl_pkpm_app *app);

/**
 * Compute implicit update of fluid-EM coupling for the PKPM system
 *
 * @param app PKPM app object
 * @param pkpm_em fluid-EM coupling updater
 * @param tcurr Current time
 * @param dt Time step size
 */
void pkpm_fluid_em_coupling_update(
  struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em, double tcurr, double dt
);

/**
 * Release resources allocated by fluid-EM coupling object for the PKPM system
 *
 * @param app PKPM app object
 * @param pkpm_em fluid-EM coupling updater to release
 */
void pkpm_fluid_em_coupling_release(
  struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em
);
