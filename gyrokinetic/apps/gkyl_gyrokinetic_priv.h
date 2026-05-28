// Private header for use in Gyrokinetic app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_dg_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_basic.h>
#include <gkyl_bc_basic_gyrokinetic.h>
#include <gkyl_bc_emission.h>
#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_boundary_flux.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_gk_neut_hamil.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_gk_collisionless_flux.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_recomb.h>
#include <gkyl_dg_updater_gk_anomalous_diffusion.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_dg_updater_lbo_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_deflated_fem_poisson.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_gk_maxwellian_correct.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_gk_maxwellian_moments.h>
#include <gkyl_gk_neut_fluid_prim_vars.h>
#include <gkyl_tok_geo.h>
#include <gkyl_velocity_map.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_null_pool.h>
#include <gkyl_position_map.h>
#include <gkyl_prim_cross_m0deltas.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_powsqrt_on_basis.h>
#include <gkyl_proj_exp_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_radiation_read.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_tok_geo.h>
#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_positivity_shift_vlasov.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>

// Definitions of private structs and APIs attached to these objects
// for use in Gyrokinetic app.

// struct for holding moment correction inputs
struct correct_all_moms_inp {
  bool correct_all_moms; // boolean if we are correcting all the moments or only density
  double iter_eps; // error tolerance for moment fixes (density is always exact)
  int max_iter; // maximum number of iterations
  bool use_last_converged; // use last iteration value regardless of convergence?
};

// data for gyrokinetic moments
struct gk_species_moment {
  struct gk_geometry *gk_geom; // geometry struct for dividing moments by Jacobian
  struct gkyl_dg_bin_op_mem *mem_geo; // memory needed in dividing moments by Jacobian
  bool is_integrated; // boolean for if computing integrated moments 
                      // integrated moments do not need to divide by Jacobian since
                      // the inverse Jacobian is already included in the computation
  int num_mom; // number of moments 

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)

  bool is_maxwellian_moms; // =True for Maxwellian moments (n, u, T/m).
  bool is_bimaxwellian_moms; // =True for BiMaxwellian moments (n, u, Tpar/m, Tperp/m).

  // MF 2025/07/18: These two objects shouldn't be here. They are needed by
  // fluid neutral moments, and we should get them directly from app->basis or
  // species->info, but these are not passes to some of the methods. Should
  // change that.
  int num_basis_conf; // Number of conf-space basis functions.
  double mass; // Species mass.
  // This method should probably be replaced by an updater.
  void (*fluid_calc_M)(const struct gk_species_moment *sm, const struct gkyl_range phase_rng,
    const struct gkyl_range conf_rng, const struct gkyl_array *fin);

  union {
    // Kinetic species .............................................. //
    struct {
      // Options for moment calculation: 
      // 1. Compute the moment directly with dg_updater_moment_gyrokinetic
      // 2. Compute the moments of the equivalent Maxwellian (n, u_par, T/m)
      // 3. Compute the moments of the equivalent Bi-Maxwellian (n, u_par, T_par/m, T_perp/m)
      //    Latter two options use specialized gkyl_gyrokinetic_maxwellian_moments updater
      union {
        struct {
          struct gkyl_gk_maxwellian_moments *gyrokinetic_maxwellian_moms; 
        };
        struct {
          struct gkyl_vlasov_lte_moments *vlasov_lte_moms; // Updater for computing LTE moments
        };
        struct {
          struct gkyl_dg_updater_moment *mcalc; 
        };
      };
    };

    // Fluid species .............................................. //
    struct {
      struct gkyl_gk_neut_fluid_prim_vars *nf_prim_vars; 
    };
  };

  // Methods chosen at runtime.
  void (*calc_func)(const struct gk_species_moment *sm, const struct gkyl_range phase_rng,
    const struct gkyl_range conf_rng, const struct gkyl_array *fin);
  void (*release_func)(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm);
  void (*diag_jacobgeo_div_func)(const struct gkyl_gyrokinetic_app *app, struct gk_species_moment *sm,
    struct gkyl_array *Jmom_in, struct gkyl_array *mom_out);
};

// Forward declare species struct.
struct gk_species;
struct gk_neut_species;

struct gk_lte {  
  struct gkyl_array *f_lte; // Maxwellian distribution.

  struct gk_species_moment moms; // Computes Maxwellian moments.

  // LTE distribution function projection object
  // also corrects the density of projected distribution function
  union {
    struct {
      struct gkyl_gk_maxwellian_proj_on_basis *proj_max; 
      struct gkyl_gk_maxwellian_correct *corr_max; 
    };
    struct {
      struct gkyl_vlasov_lte_proj_on_basis *proj_lte; 
      struct gkyl_vlasov_lte_correct *corr_lte; 
    };
  };

  long n_iter; // Total number of iterations from correcting moments.
  long num_corr; // Total number of times the correction updater is called.
  bool correct_all_moms; // Boolean if we are correcting all the moments.
  gkyl_dynvec corr_stat;
  bool is_first_corr_status_write_call;

  // Methods chosen at runtime.
  void (*from_moms_func)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
    struct gk_lte *lte, const struct gkyl_array *moms_lte);
  void (*from_f_func)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
    struct gk_lte *lte, const struct gkyl_array *fin);
  void (*write_max_corr_status_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gk_ns);
  void (*release_func)(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte);
};

struct gk_rad_drag {  
  enum gkyl_radiation_id radiation_id; // Type of radiation.
  bool write_diagnostics; // Whether to write diagnostics out.
  int num_cross_collisions; // Number of species we cross-collide with.
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // Cross-species we collide with.
  struct gk_neut_species *collide_with_neut[GKYL_MAX_SPECIES]; // Neutral cross-species we collide with.
  int collide_with_idx[2*GKYL_MAX_SPECIES]; // Index of species we collide with.
  bool is_neut_species[2*GKYL_MAX_SPECIES]; // Flag of whether neutral or Gk species.
  
  // Drag coefficients in vparallel and mu for each species being collided with.
  struct gkyl_gk_rad_drag *vnu_surf;
  struct gkyl_gk_rad_drag *vnu;
  struct gkyl_gk_rad_drag *vsqnu_surf;
  struct gkyl_gk_rad_drag *vsqnu;
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars;
  struct gkyl_array *rad_fit_ne[2*GKYL_MAX_SPECIES];

  struct gk_species_moment moms[2*GKYL_MAX_SPECIES]; // Moments needed in radiation update (need number density).

  struct gk_species_moment m2; // M2 of radiation update (needed for emissivity).
  struct gkyl_array *emissivity[2*GKYL_MAX_SPECIES];
  struct gkyl_array *emissivity_host[2*GKYL_MAX_SPECIES];
  struct gkyl_array *emissivity_rhs;
  struct gkyl_array *emissivity_denominator;

  struct gkyl_array *nvnu_surf; // Total vparallel radiation drag surface expansion including density scaling.
  struct gkyl_array *nvnu; // Total vparallel radiation drag volume expansion including density scaling.
  struct gkyl_array *nvsqnu_surf; // Total mu radiation drag surface expansion including density scaling.
  struct gkyl_array *nvsqnu; // Total mu radiation drag volume expansion including density scaling.

  struct gkyl_array **vtsq_min_per_species;  // Smallest vtsq that radiation is calculated (one for each
                                             // fit), divided by configuration space normalization
  struct gk_species_moment prim_moms; // Primitive moments.
  struct gkyl_array *vtsq; // Thermal speed squared.
  struct gkyl_array *m0; // M0 moment.

  // Host-side copies for I/O.
  struct gkyl_array *nvnu_surf_host; 
  struct gkyl_array *nvnu_host; 
  struct gkyl_array *nvsqnu_surf_host; 
  struct gkyl_array *nvsqnu_host; 

  gkyl_dg_updater_collisions *drag_slvr; // Radiation solver.

  struct gk_species_moment integ_moms; // Integrated moments.
  double *red_integ_diag, *red_integ_diag_global; // Reduced integrated moments.
  gkyl_dynvec integ_diag; // Integrated moments in time.
  bool is_first_integ_write_call; // Whether dynvec is being written for the 1st time.
  
  struct gkyl_msgpack_map_elem* io_meta_surfmu; // Metadata for I/O of mu surf array.
  int io_meta_surfmu_len; // Number of elements in io_meta_surfmu.
  struct gkyl_msgpack_map_elem* io_meta_surfvpar; // Metadata for I/O of vpar surf array.
  int io_meta_surfvpar_len; // Number of elements in io_meta_surfvpar.

  // Methods chosen at runtime:
  void (*moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);
  void (*rhs_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gk_rad_drag *rad, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_drag_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_emissivity_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);
  void (*write_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks);
};

struct gk_collisionless {
  enum gkyl_gk_collisionless_type collisionless_id; // Type of collisionless terms.
  bool write_diagnostics; // Whether to write diagnostics out.
  double scale_fac; // Factor multiplying collisionless terms.

  // Organization of the different equation objects and the required data and solvers
  union {
    // Charged (gyrokinetic) species ............................................ //
    struct {
      struct gkyl_array *flux_surf; // Array for surface phase space flux
      struct gkyl_array *flux_surf_ho; // Host array for surface phase space flux
      struct gkyl_array *apar; // A_parallel.
      struct gkyl_array *apardot; // d/dt A_parallel.

      struct gkyl_gk_collisionless_flux *surf_flux_op; // Collisionless fluxes.
      gkyl_dg_updater_gyrokinetic *slvr; // Collisionless solver.
 
      // Methods chosen at runtime.
      void (*flux_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
        struct gk_collisionless *gkcls, const struct gkyl_array *fin);
    };
    // Neutral (Vlasov) species ............................................ //
    struct {
      struct gkyl_array *alpha_surf; // array for surface phase space flux (v^i = v . e^i)
      struct gkyl_array *sgn_alpha_surf; // array for the sign of the surface phase space flux at quadrature points
                                         // utilized for numerical flux function
                                         // F = alpha_surf/2 ( (f^+ + f^-) - sign_alpha_surf*(f^+ - f^-) )
      struct gkyl_array *const_sgn_alpha; // boolean array for if the surface phase space flux is single signed
                                          // if true, numerical flux function inside kernels simplifies to
                                          // F = alpha_surf*f^- (if sign_alpha_surf = 1), 
                                          // F = alpha_surf*f^+ (if sign_alpha_surf = -1)

      gkyl_dg_updater_vlasov *vlasov_slvr; // Vlasov solver.

    };
  };
  void (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*fdot_scaling)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_collisionless *gkcls, struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_range *rng);
  void (*write_diags_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_collisionless *gkcls, double tm, int frame);
  // Methods for neutral species.
  void (*rhs_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
    struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_diags_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
    struct gk_collisionless *gkcls, double tm, int frame);
};

struct gk_lbo_collisions {  
  enum gkyl_collision_id collision_id; // type of collisions
  bool write_diagnostics; // Whether to write diagnostics out.
  bool not_in_dfdt; // Whether to not add collision contribution to df/dt 

  struct gkyl_array *self_nu; // Self-collision frequency.
  struct gkyl_array *boundary_corrections; // LBO boundary corrections.
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator.
  struct gkyl_array *nu_sum; // Sum of collision frequencies.
  struct gkyl_array *prim_moms, *nu_prim_moms; // Primitive moments.
  struct gkyl_array *nu_sum_host, *nu_prim_moms_host; // Host arrays for I/O.
  bool norm_nu_self; // Whether to compute self-species collision frequency in space and time.
  double norm_nu_fac_self; // Self collision frequency without factor of n_s/(2*v_ts^2)^(3/2).
  double vtsq_min; // Minimum vtsq.
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value.
  struct gkyl_array *nu_boundary_corrections; // Boundary corrections multiplied by nu.
  struct gkyl_array *nu_moms; // Moments multiplied by nu.
  struct gkyl_array *m2self; // Second velocity moment of self species.
  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator

  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  bool norm_nu_cross; // Whether to compute cross-species collision frequency in space and time.
  double norm_nu_fac_cross[GKYL_MAX_SPECIES]; // Cross collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2).
  double alpha_E_fac[GKYL_MAX_SPECIES]; // Time-independent factor in alpha_E.
  double betaGreenep1; // Galue of Greene's factor beta + 1.
  double delta_sr; // Free parameter in relationship between alpha_E and nu_sr.
  double other_m[GKYL_MAX_SPECIES]; // Masses of species colliding with.
  struct gkyl_array *other_prim_moms[GKYL_MAX_SPECIES]; // Self-primitive moments of species colliding with.
  struct gkyl_array *cross_prim_moms[GKYL_MAX_SPECIES]; // Cross-primitive moments.
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // Cross-species collision frequencies.
  struct gkyl_array *cross_nu_prim_moms; // Weak multiplication of collision frequency and primitive moments.
  struct gkyl_array *alpha_E; // Morse's alpha_E factor.
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  
  struct gk_species_moment moms; // Moments needed in LBO (M0, M1, M2).

  gkyl_dg_updater_collisions *coll_slvr; // collision solver

  // Pointers to methods chosen at runtime.
  void (*moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gk_lbo_collisions *lbo, const struct gkyl_array *fin);
  void (*self_nu_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gk_lbo_collisions *lbo, const struct gkyl_array *fin);
  void (*cross_nu_func)(gkyl_gyrokinetic_app *app, const struct gk_species *s,
    struct gk_lbo_collisions *lbo, int coll_idx);
  void (*alpha_E_func)(gkyl_gyrokinetic_app *app, const struct gk_species *s,
    struct gk_lbo_collisions *lbo, int coll_idx);
  void (*cross_moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gk_lbo_collisions *lbo);
  void (*rhs_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
};

struct gk_bgk_collisions {  
  enum gkyl_collision_id collision_id; // Type of collisions.
  bool write_diagnostics; // Whether to write diagnostics out.

  struct gkyl_array *self_nu; // Self-collision frequency.
  struct gkyl_array *nu_sum; // Sum of collision frequencies.
  struct gkyl_array *nu_sum_host; // Host arrays for I/O.
  bool norm_nu_self; // Whether to compute self-species collision frequency in space and time.
  double norm_nu_fac_self; // Self collision frequency without factor of n_s/(2*v_ts^2)^(3/2).
  double vtsq_min; // Minimum vtsq.
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  bool norm_nu_cross; // Whether to compute cross-species collision frequency in space and time.
  double norm_nu_fac_cross[GKYL_MAX_SPECIES]; // Cross collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2).
  double alpha_E_fac[GKYL_MAX_SPECIES]; // Time-independent factor in alpha_E.
  double betaGreenep1; // Galue of Greene's factor beta + 1.
  double delta_sr; // Free parameter in relationship between alpha_E and nu_sr.
  double other_m[GKYL_MAX_SPECIES]; // Masses of species colliding with.
  struct gkyl_array *other_prim_moms[GKYL_MAX_SPECIES]; // Self-primitive moments of species colliding with.
  struct gkyl_array *cross_prim_moms; // Cross-primitive moments.
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // Cross-species collision frequencies.
  struct gkyl_array *alpha_E; // Morse's alpha_E factor.
  struct gkyl_gyrokinetic_cross_prim_moms_bgk *cross_calc; // Cross-species moment computation.
  
  struct gkyl_array *nu_fmax; // Collision frequency times Maxwellian.
  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep).

  bool implicit_step; // Whether or not to take an implcit BGK step.
  double dt_implicit; // Timestep used by the implicit collisions.

  // Pointers to methods chosen at runtime.
  void (*moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*moms_func_implicit)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*self_nu_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_bgk_collisions *bgk);
  void (*cross_nu_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_bgk_collisions *bgk, int coll_idx);
  void (*alpha_E_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_bgk_collisions *bgk, int coll_idx);
  void (*cross_moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_bgk_collisions *bgk, int coll_idx);
  void (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);
  void (*rhs_func_implicit)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  // For neutral species (will hopefully remove when we unify species types).
  void (*moms_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*moms_func_implicit_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*self_nu_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk);
  void (*cross_nu_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, int coll_idx);
  void (*alpha_E_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, int coll_idx);
  void (*cross_moms_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, int coll_idx);
  void (*rhs_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);
  void (*rhs_func_implicit_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
    struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);
  void (*write_mom_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);
};

struct gk_anomalous_diff {
  enum gkyl_gk_anomalous_diff_id anom_diff_id; // Type of diffusion.
  bool write_diagnostics; // Whether to write diagnostics out.
  struct gkyl_array *diffD; // Diffusivity.
  struct gkyl_dg_updater_gk_anomalous_diffusion *slvr; // Anomalous diffusion equation solver.
  // Methods chosen at runtime.
  void (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_diags_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_anomalous_diff *gkad, double tm, int frame);
};

enum gkyl_species_bflux_type {
  GK_SPECIES_BFLUX_NONE = 0, // Do nothing.
  GK_SPECIES_BFLUX_CALC_FLUX, // Put boundary fluxes in ghost of rhs.
  GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS, // Also step the moments of bflux.
  GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS, // Also compute/write diagnostics.
};

#define BFLUX_MAX_MOM_NAMES 12

struct gk_boundary_fluxes {
  bool allocated_solver, allocated_moms, allocated_diags;
  // Objects used to compute boundary fluxes.
  int num_boundaries; // Number of boundaries to compute bfluxes at.
  int boundaries_dir[2*GKYL_MAX_CDIM]; // Direction of bflux boundaries.
  enum gkyl_edge_loc boundaries_edge[2*GKYL_MAX_CDIM]; // Edge of bflux boundaries.
  struct gkyl_range *boundaries_conf_skin[2*GKYL_MAX_CDIM]; // Conf-space ghost range of boundaries.
  struct gkyl_range *boundaries_conf_ghost[2*GKYL_MAX_CDIM]; // Conf-space ghost range of boundaries.
  struct gkyl_range *boundaries_phase_skin[2*GKYL_MAX_CDIM]; // Phase-space skin range of boundaries.
  struct gkyl_range *boundaries_phase_ghost[2*GKYL_MAX_CDIM]; // Phase-space ghost range of boundaries.
  struct gkyl_range *boundaries_phase_ghost_nosub; // Not a sub range of local_ext.
  struct gkyl_range *boundaries_conf_skin_fullx[2*GKYL_MAX_CDIM]; // Whole x range (for IWL).
  int num_eqns; // Number of equations with boundary fluxes.
  const struct gkyl_dg_eqn **eqns; // Array of equation objects with boundary fluxes.
  gkyl_boundary_flux *flux_slvr[2*GKYL_MAX_CDIM]; // boundary flux solver
  struct gkyl_array **flux; // Array storing boundary fluxes.
  // Objects used for calculating moments.
  int num_calc_moms; // Number of moments of boundary fluxes to compute.
  enum gkyl_distribution_moments calc_mom_names[BFLUX_MAX_MOM_NAMES]; // Names of moments calculated.
  bool *is_hamiltonian_mom; // True if need Hamiltonian moments.
  bool a_hamiltonian_mom; // There is one Hamiltonian moment.
  struct gkyl_bc_basic_gyrokinetic *gfss_bc_op[2*GKYL_MAX_CDIM]; // Applies BCs to bmag and phi.
  struct gkyl_array *bc_buffer; // Buffer used by gfss_bc_op;
  struct gkyl_array **f, **f1, **fnew; // Boundary flux through each boundary (one for each RK stage).
  struct gk_species_moment *moms_op; // Moments calculator.
  // Objects used for calculating diagnostics.
  int *diag_mom_idx; // Index of each diag mom in the array of calc moms.
  int *diag_int_mom_idx; // Index of each diag int mom in the array of calc moms.
  int *add_mom_idx; // Index of each additional mom in the array of calc moms.
  struct gkyl_translate_dim **transdim; // Updater that projects to lower dim.
  struct gkyl_rect_decomp *decomp_surf[GKYL_MAX_CDIM]; // Surface decomposition object.
  struct gkyl_comm *comm_surf[GKYL_MAX_CDIM]; // Surface communicator object.
  struct gkyl_range surf_local[GKYL_MAX_CDIM], surf_local_ext[GKYL_MAX_CDIM]; // local, local-ext surface ranges.
  struct gkyl_rect_grid grid_surf[GKYL_MAX_CDIM]; // Surface grid.
  struct gkyl_array **mom_surf, **mom_surf_ho; // Moments on the surface.
  struct gkyl_array_integrate **integ_op; // Operator that integrates over volume.
  double *int_moms_local, *int_moms_global; // Integrated moments in this time step.
  gkyl_dynvec *intmom; // Integrated moments of the boundary fluxes.
  double **intmom_cumm_buff; // Cummulative (in time) integrated moments of the boundary fluxes.
  bool is_first_intmom_write_call[2*GKYL_MAX_CDIM]; // Flag 1st writing of blux_intmom.
  // Function pointers to various methods.
  void (*bflux_rhs_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
    const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*bflux_calc_moms_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
    const struct gkyl_array *rhs, struct gkyl_array **bflux_moms);
  void (*bflux_get_flux_func)(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
    struct gkyl_array *out, const struct gkyl_range *out_rng);
  void (*bflux_get_flux_mom_func)(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
    enum gkyl_distribution_moments mom_type, struct gkyl_array **bflux_moms, struct gkyl_array *out, const struct gkyl_range *out_rng);
  void (*bflux_clear_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fin, double val);
  void (*bflux_scale_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fin, double val);
  void (*bflux_step_f_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fout,
    double dt, struct gkyl_array **fin);
  void (*bflux_set_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fout,
    double fac, struct gkyl_array **fin);
  void (*bflux_accumulate_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fout,
    double fac, struct gkyl_array **fin);
  void (*bflux_combine_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
    struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2);
  void (*bflux_copy_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
    struct gkyl_array **fout, struct gkyl_array **fin);
  void (*bflux_calc_voltime_int_mom_func)(gkyl_gyrokinetic_app* app, void *species,
    struct gk_boundary_fluxes *bflux, double tm);
  void (*bflux_calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, void *species,
    struct gk_boundary_fluxes *bflux, double tm);
  void (*bflux_write_integrated_mom_func)(gkyl_gyrokinetic_app *app, void *species,
    struct gk_boundary_fluxes *bflux);
  void (*bflux_write_mom_func)(gkyl_gyrokinetic_app* app, void *species,
    struct gk_boundary_fluxes *bflux, double tm, int frame);
};

struct gk_recycle_wall {
  // Recycling wall boundary conditions.
  // Input struct must include number of impact species,
  // name of the species, and recycling fraction.
  struct gkyl_gyrokinetic_emission_inp *params;
  int num_species;

  int dir; // Direction orthogonal to the recycling boundary..
  enum gkyl_edge_loc edge; // Edge of the recycling boundary.
  gkyl_boundary_flux *f0_flux_slvr; // Operator that computes the neutral flux.
  struct gkyl_dg_bin_op_mem *mem_geo; // Memory for calculating flux ratio.

  struct gkyl_rect_grid emit_grid; // Ghost grid of the emitting boundary.
  struct gkyl_range *emit_ghost_r; // Ghost range of the emitting boundary (sub of local_ext).
  struct gkyl_range *emit_skin_r; // Skin range of the emitting boundary (sub of local_ext).
  struct gkyl_range emit_buff_r; // Ghost range of the emitting boundary.
  struct gkyl_range emit_cbuff_r; // Conf-space ghost range of the emitting boundary.
  struct gkyl_range emit_normal_r;

  struct gk_species *impact_species[GKYL_MAX_SPECIES]; // Pointers to impacting species

  struct gkyl_rect_grid impact_grid[GKYL_MAX_SPECIES]; // Ghost grid of the impacting boundary.
  struct gkyl_range *impact_ghost_r[GKYL_MAX_SPECIES];  // Ghost range of the impacting boundary (sub of local_ext).
  struct gkyl_range impact_buff_r[GKYL_MAX_SPECIES];  // Ghost range of the impacting boundary.
  struct gkyl_range impact_cbuff_r[GKYL_MAX_SPECIES]; // Conf-space ghost range of the impacting boundary.

  struct gkyl_array *bc_buffer; // Fixed buffers for recycle BCs.
  struct gkyl_array *phase_flux_gk[GKYL_MAX_SPECIES]; // Array to put phase-flux into.
  struct gkyl_array *m0_flux_gk[GKYL_MAX_SPECIES]; // M0 moment of ion flux.
  struct gkyl_dg_updater_moment *m0op_gk[GKYL_MAX_SPECIES]; // M0 moment solver for ion flux.

  struct gkyl_array *unit_phase_flux_neut; // Unit-density neutral flux.  
  struct gkyl_array *unit_m0_flux_neut; // MO moment of unit-density flux.
  struct gkyl_dg_updater_moment *m0op_neut; // M0 moment solver for unit-density flux.

  struct gkyl_array *spectrum[GKYL_MAX_SPECIES]; // Unit-density Maxwellian is copied and scaled here.
  struct gkyl_range impact_normal_r[GKYL_MAX_SPECIES]; // Phase-space range w/ only velocities towards the boundary.
  struct gkyl_array *f_emit; // Array to fill neutral ghost cell with scaled Maxwellian distf.

  // Objects used for diagnostics.
  bool write_diagnostics;
  struct gkyl_array *f_diag;
  struct gkyl_array *emit_flux;
  struct gkyl_array *diag_out; 
  struct gkyl_array *diag_out_ho;
  // Function pointers to diagnostic functions.
  void (*write_flux_func)(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
    struct gk_recycle_wall *recyc, double tm, int frame);
};

struct gk_react {
  int num_react; // Number of reactions.
  bool all_gk; // Boolean for if reactions are only between gyrokinetic species.
  bool write_diagnostics; // Whether to write diagnostics out.
  struct gkyl_gyrokinetic_react_type react_type[GKYL_MAX_REACT]; // Input struct for type of reactions.

  struct gkyl_array *f_react; // Distribution function array which holds update for each reaction
                              // form depends on type_self, e.g., for ionization and type_self == GKYL_SELF_ELC
                              // f_react = n_donor*(fmax1(n_elc, upar_elc, vtiz1^2) + fmax2(n_elc, upar_donor, vtiz2^2) - f_elc)
                              // RHS update is then obtained by incrementing rhs += coeff_react*f_react

  enum gkyl_react_id react_id[GKYL_MAX_REACT]; // What type of reaction (ionization, charge exchange, recombination).
  enum gkyl_react_self_type type_self[GKYL_MAX_REACT]; // What is the role of species in this reaction
  struct gk_species *species_elc[GKYL_MAX_REACT]; // Electron species being reacted with.
  struct gk_species *species_ion[GKYL_MAX_REACT]; // Ion species being reacted with.
  int elc_idx[GKYL_MAX_REACT]; // Integer index of electron species being reacted with.
  int ion_idx[GKYL_MAX_REACT]; // Integer index of ion species being reacted with.
  int donor_idx[GKYL_MAX_REACT]; // Integer index of donor species being reacted with.
  int partner_idx[GKYL_MAX_REACT]; // Integer index of neut species in cx reaction.

  struct gkyl_array *coeff_react[GKYL_MAX_REACT]; // Reaction rate.
  struct gkyl_array *coeff_react_host[GKYL_MAX_REACT]; // Reaction rate.
  struct gkyl_array *vt_sq_iz1[GKYL_MAX_REACT]; // Ionization temperature.
  struct gkyl_array *vt_sq_iz2[GKYL_MAX_REACT]; // Ionization temperature.
  struct gkyl_array *Jm0_elc[GKYL_MAX_REACT]; // J*electron density.
  struct gkyl_array *Jm0_ion[GKYL_MAX_REACT]; // J*ion density.
  struct gkyl_array *Jm0_donor[GKYL_MAX_REACT]; // J*donor density.
  struct gkyl_array *Jm0_partner[GKYL_MAX_REACT]; // J*partner density.

  struct gkyl_array *react_lte_moms[GKYL_MAX_REACT]; // LTE/Maxwellian moments used for Maxwellian projection.
  struct gkyl_array *u_i[GKYL_MAX_REACT]; // Vector flow velocity for neutrals (ux, uy, uz).
  struct gkyl_array *u_i_dot_b_i[GKYL_MAX_REACT]; // u_i . b_i (Cartesian compoments of magnetic field unit vector).
  struct gkyl_array *vt_sq_donor[GKYL_MAX_REACT]; // Donor thermal velocity.
  struct gkyl_array *upar_ion[GKYL_MAX_REACT]; // Ion vector parallel flow velocity upar b_i.
  struct gkyl_array *vt_sq_ion[GKYL_MAX_REACT]; // Ion thermal velocity.
  struct gkyl_array *vt_sq_partner[GKYL_MAX_REACT]; // Neutral (partner) thermal velocity.
  
  union {
    // Ionization.
    struct {
      struct gkyl_dg_iz *iz[GKYL_MAX_REACT];
    };
    // Recombination.
    struct {
      struct gkyl_dg_recomb *recomb[GKYL_MAX_REACT];
    };
    // Charge exchange.
    struct {
      struct gkyl_dg_cx *cx[GKYL_MAX_REACT];
    };
  };

  // Methods chosen at runtime:
  void (*cross_moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);
  void (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *s,
    struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, struct gk_react *gkr,
    int ridx, double tm, int frame);
  // Neutral methods (MF 2025/10/28: to get rid of by unifying species types).
  void (*cross_moms_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
    struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);
  void (*rhs_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
    struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
    int ridx, double tm, int frame);
};

struct gk_scaling {
  bool write_diagnostics; // Whether to write diagnostics.
  enum gkyl_gk_species_scaling_type type; // Type of scaling operation.

  union {
    struct {
      // Objects for GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE.
      int num_boundaries; // Number of boundaries.
      int boundaries_dir[GKYL_MAX_CDIM*2]; // Direction of boundaries.
      enum gkyl_edge_loc boundaries_edge[GKYL_MAX_CDIM*2]; // Edge of boundaries.
      int elc_idx; // Index of electron species.
      int ion_idx; // Index of ion species.
      double recycling_coeff; // Recycling coefficient.
     
      struct gkyl_array *Jm0_init; // Initial number density times jacobgeo.
      struct gkyl_array *reactivity; // Reactivity.
      struct gkyl_array *dndt_react; // Reaction contribution.
      struct gkyl_dg_iz *iz_react_calc; // Operator to compute ionization rate.
      struct gkyl_array_integrate *integrate_op; // Operator that integrates an array.
      double react_vol_integ, *react_vol_integ_local, *react_vol_integ_global; // Volume integral of reaction contribution.
    
      struct gkyl_range boundaries_conf_ghost[2*GKYL_MAX_CDIM]; // Ghost cells to where integrated M0 is computed.
     
      double *bflux_m0_vol_integ_local; // Local volume integral of boundary flux M0.
    };
    struct {
      // Info for GKYL_GK_SPECIES_SCALING_FIXED_FRACTION.
      int ref_species_idx; // Index of reference species species.
      double fixed_fraction; // Fraction of reference species density.
    };
    struct {
      // Objects for GKYL_GK_SPECIES_SCALING_BOLTZMANN.
      struct gkyl_proj_exp_on_basis *proj_exp; // Operator to project exp(A).
      struct gkyl_array *sheath_val; // Value at the sheath entrance, copied to interior cells.
      struct gkyl_array *buffer_conf; // Conf-space buffer array.
    };
  };

  // Methods chosen at runtime.
  void (*cross_moms_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);
  void (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*apply_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[]);
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_scaling *sca, int ridx, double tm, int frame);
  // Neutral species function pointers (MF 2025/11/03: to be removed when we unify species types).
  void (*cross_moms_func_neut)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *ns,
    struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);
  void (*rhs_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
    struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*apply_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
    struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[]);
  void (*write_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *ns,
    struct gk_scaling *sca, int ridx, double tm, int frame);
};

// Context for c2p function passed to proj_on_basis.
struct gk_proj_on_basis_c2p_func_ctx {
  int cdim, vdim;
  struct gkyl_position_map *pos_map;
  struct gkyl_velocity_map *vel_map;
};

struct gk_proj {
  enum gkyl_projection_id proj_id; // Type of projection.
  struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
  union {
    // Function projection.
    struct {
      struct gkyl_proj_on_basis *proj_func; // Projection operator for specified function.
      struct gkyl_array *proj_host; // Array for projection on host-side if running on GPUs.
    };
    // Maxwellian and Bi-Maxwellian projection from primitive moments.
    struct {
      struct gkyl_array *dens; // Density (on device for GK Maxwellian Gaussian).

      struct gkyl_array *prim_moms_host; // Host-side prim_moms for initialization with proj_on_basis
      struct gkyl_array *prim_moms; // Prim_moms we pass to Maxwellian projection object (potentially on device)
      bool maxwellian_moms_from_file;
      bool bimaxwellian_moms_from_file;
      bool dens_from_file;
      bool upar_from_file;
      bool temp_from_file;
      bool temppar_from_file;
      bool tempperp_from_file;

      bool correct_all_moms; // Boolean if we are correcting all the moments.

      struct gkyl_proj_on_basis *proj_dens; // Projection operator for density
      struct gkyl_array *vtsq; // vt^2 = T/m (temperature/mass) (on device for GK Maxwellian Gaussian).
      struct gkyl_proj_on_basis *proj_temp; // Projection operator for temperature.
      struct gkyl_array *vtsqpar; // Host-side vt_par^2 = Tpar/m (parallel temperature/mass).
      struct gkyl_array *vtsqperp; // Host-side vt_perp^2 = Tperp/m (perpendicular temperature/mass).
      struct gkyl_proj_on_basis *proj_temppar; // Projection operator for parallel temperature.
      struct gkyl_proj_on_basis *proj_tempperp; // Projection operator for perpendicular temperature.

      struct gkyl_array *gaussian_profile; // Shape of the source in configuration space (on device).
      struct gkyl_proj_on_basis *proj_gaussian; // Projection operator for the shape of the source in config space.

      union {
        struct { 
          struct gkyl_array *upar; // upar for GK Maxwellian/Bi-Maxwellian projection (on device for GK Maxwellian Gaussian)
          struct gkyl_proj_on_basis *proj_upar; // Projection operator for upar for GK Maxwellian/Bi-Maxwellian projection.
          struct gkyl_gk_maxwellian_proj_on_basis *proj_max; // Maxwellian projection object for GK.
          struct gkyl_gk_maxwellian_correct *corr_max; // Maxwellian correction object for GK.
        };
        struct { 
          struct gkyl_array *udrift; // Host-side udrift for Vlasov neutrals LTE projection.
          struct gkyl_proj_on_basis *proj_udrift; // Projection operator for udrift for Vlasov neutrals LTE projection.
          struct gkyl_vlasov_lte_proj_on_basis *proj_lte; // Maxwellian projection object for Vlasov neutrals.
          struct gkyl_vlasov_lte_correct *corr_lte; // Maxwellian correction object for Vlasov neutrals.
        };     
      };  
    };
  };

  // Functions chosen at runtime.
  void (*projection_calc)(gkyl_gyrokinetic_app *app, struct gk_species *s, 
    struct gk_proj *proj, struct gkyl_array *f, double tm);
  void (*neut_calc_func)(gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
    struct gk_proj *proj, struct gkyl_array *f, double tm);
  void (*moms_correct)(gkyl_gyrokinetic_app *app, struct gk_species *s, 
    struct gk_proj *proj, struct gkyl_array *f, double tm);
  void (*release_func)(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj);
};

struct gk_adapt_source {
  bool adapt_particle, adapt_energy; // Adaptation flags.
  struct gk_species *adapt_species; // Non-owning pointer to the species whose particle loss is being adapted; lifetime managed elsewhere.
  int adapt_species_idx; // Index of the species whose particle loss is being adapted.

  int num_boundaries; // Number of boundaries to adapt to.
  int dir[2*GKYL_MAX_CDIM]; // Direction to adapt.
  enum gkyl_edge_loc edge[2*GKYL_MAX_CDIM]; // Edge to adapt.
  struct gkyl_range boundaries_phase_ghost[2*GKYL_MAX_CDIM]; // Range of computation of the bflux and moment (ALL phase space ghost).
  struct gkyl_range boundaries_conf_ghost[2*GKYL_MAX_CDIM]; // Range of integration in each boundary (SOL config space ghost).

  struct gk_species_moment integ_threemoms; // Integrated moment updater.
  double *red_integ_mom, *red_integ_mom_global; // For reduction.

  struct gkyl_array_integrate *integrate_op; // Operator that integrates an array.
  struct gkyl_array *bflux_m0, *bflux_m2; // Boundary flux moments.
  double *integ_m0, *integ_m2; // Integrated boundary flux moments for one boundary and one mpi process. (on device)

  double particle_src_curr, energy_src_curr; // current injection rates of the source.
  double particle_rate_loss, energy_rate_loss; // Loss rates we adapt to.
  double temperature_curr; // Current density and temperature.
};

struct gk_source {
  enum gkyl_source_id source_id; // Type of source.
  bool evolve; // Whether the source is time dependent.
  struct gkyl_array *source; // Applied source.
  struct gkyl_array *source_host; // Host copy for use in IO and projecting.
  struct gk_proj proj_source[GKYL_MAX_SOURCES]; // Projector for source.
  int num_sources; // Number of sources.

  int num_diag_mom; // Number of diagnostics moments.
  struct gk_species_moment *moms; // Diagnostic moments.
  int num_diag_int_mom; // Number of volume-integrated diagnostic moments.
  struct gk_species_moment integ_moms; // Integrated moments.
  double *red_integ_diag, *red_integ_diag_global; // For reduction of integrated moments.
  gkyl_dynvec integ_diag; // Integrated moments reduced across grid.
  gkyl_dynvec temp_diag, part_diag; // Src temperature and particle count diags.
  bool is_first_integ_write_call; // Flag for integrated moments dynvec written first time.
  struct gk_adapt_source adapt[GKYL_MAX_SOURCES]; // Adaptation source.
  int num_adapt_sources; // Number of adaptive sources.
  // Functions chosen at runtime.
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);
  void (*write_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks);
  void (*adapt_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, struct gk_source *src, 
    struct gkyl_array *f_buffer, struct gkyl_array **bflux_moms[], double tm);
};

struct gk_damping {
  enum gkyl_gyrokinetic_damping_type type; // Type of damping term.
  bool evolve; // Whether the source is time dependent.
  bool write_rate; // Whether to write damping-rate diagnostics.
  bool write_fbar; // Whether to write low-pass-filter fbar diagnostics.
  bool cellwise_const; // Whether low-pass-filter fbar is stored as one value per cell.
  bool do_not_reset_fbar; // Whether runtime reset should preserve existing fbar.
  struct gkyl_array *rate; // Damping rate.
  struct gkyl_array *rate_host; // Host copy for use in IO and projecting.
  // Low-pass filter variables
  bool fbar_initialized; // Whether fbar is allocated/owned and should be released.
  struct gkyl_array *fbar; // Filtered/averaged distribution function.
  struct gkyl_array *fbar1, *fbarnew; // SSPRK3 stage arrays for filtered distribution.
  struct gkyl_array *fbar_host; // Host copy of fbar for use in IO.
  // Functions chosen at runtime.
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_rate_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_fbar_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*advance_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_damping *damp, 
    const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
    struct gkyl_array *f_buffer,
    struct gkyl_array *rhs, struct gkyl_array *cflrate);
  void (*set_fbar_to_f_func)(const struct gk_species *gks, struct gk_damping *damp,
    const struct gkyl_array *f);
  void (*calc_fbar_rhs_func)(const struct gk_damping *damp,
    const struct gkyl_array *fin, const struct gkyl_array *fbar_in, struct gkyl_array *rhs_fbar);
  void (*forward_euler_func)(struct gk_species *gks, const struct gkyl_array *fin,
    const struct gkyl_array *fbar_in, struct gkyl_array *fbar_out, double dt);
  void (*combine_func)(struct gk_species *gks, struct gkyl_array *fout, double c1,
    const struct gkyl_array *f1, double c2, const struct gkyl_array *f2,
    const struct gkyl_range *rng);
  void (*copy_func)(struct gk_species *gks, struct gkyl_array *fout,
    const struct gkyl_array *fin, const struct gkyl_range *range);
};

struct gk_fdot_multiplier {
  enum gkyl_gyrokinetic_fdot_multiplier_type type; // Type of multiplicative function term.
  bool write_diagnostics; // Whether to write diagnostics out.
  bool evolve; // Whether the multiplicative function is time dependent.
  struct gkyl_array *multiplier; // Damping rate.
  struct gkyl_array *multiplier_host; // Host copy for use in IO and projecting.
  struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
  struct gkyl_loss_cone_mask_gyrokinetic *lcm_proj_op; // Operator that projects the loss cone mask.
  double *bmag_max; // Maximum magnetic field amplitude.
  double *bmag_max_coord; // Location of bmag_max.
  double *phi_m, *phi_m_global; // Electrostatic potential at bmag_max.
  // Functions chosen at runtime.
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*advance_times_rate_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out);
  void (*advance_times_cfl_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out);
};

struct gk_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // Type of source_bgk source.
  bool write_diagnostics; // Whether to write diagnostics out.
  double norm_power; // Normalized power, 2*P/(vdim_phys*m).
  struct gkyl_array *rate; // Sourcing rate.
  struct gkyl_array *Jrate; // Sourcing rate times the conf-space Jacobian.
  struct gkyl_array *vtsq_shape; // Spatial profile of the Maxwellian's v_t^2.
  struct gkyl_array *Jrate_vtsq_shape; // Jrate times vtsq_shape.
  struct gkyl_array *Jrate_fmax; // Jrate times the Maxwellian.
  struct gkyl_array *Jrate_mom; // Jrate times a velocity moment.
  struct gkyl_array *Jrate_cap; // Max value for Jrate_mom.
  struct gkyl_array_integrate *vol_integ_op; // Volume integrator.
  double *volint_local, *volint_global; // Local and global volume integrals.
  struct gkyl_array *M0dot, *M1dot, *M2dot; // Source rates for external sourcing
  struct gkyl_array *M0dot_host, *M1dot_host, *M2dot_host; // Host source rates for external sourcing
  double injection_time; // Injection time for external heating model
  double damping_factor; // Damping factor used to prevent negative densities
  double vtsq_amplitude; // Amplitude of squared thermal speed.
  struct gkyl_bgk_collisions *bgk_op; // BGK operator.
  bool implicit_step; // Whether or not to take an implcit BGK step.
  double dt_implicit; // Timestep used by the implicit collisions.
  struct gk_species_moment integ_mom_op; // Integrated moments.
  struct gk_species_moment correct_mom_op; // Correction moments.
  gkyl_dynvec vtsq_amp_diag; // Stores vtsq_amplitude for diagnostics.
  double *red_integ_diag, *red_integ_diag_global; // Reduced integrated moments.
  gkyl_dynvec integ_diag; // Integrated moments of the source.
  bool is_first_diag_dynvec_write_call; // Whether dynvec is being written for the first time.
  // Methods chosen at runtime.
  void (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_diags_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_source_bgk *src, double tm, int frame);
  void (*calc_integrated_diags_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_source_bgk *src, double tm);
  void (*write_integrated_diags_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_source_bgk *src);
};

struct gk_positivity {
  // Updater that enforces positivity by shifting f.
  enum gkyl_gyrokinetic_positivity_type type; // Type of positivity enforcement algorithm.
  bool quasineut_rescale; // Whether to rescale this species to enforce quasineutrality in the simulation.
  bool write_diagnostics; // Whether to output diagnostics.

  struct gkyl_array *fbuffer_ptr; // Pointer to an array were we store delta f.
  struct gkyl_array *delta_m0; // Number density of the positivity shift.
  union {
    struct {
      struct gkyl_positivity_shift_gyrokinetic *shift_op_gk;
      struct gkyl_array *delta_m0s_tot; // Density of total positivity shift (like-species).
      struct gkyl_array *delta_m0r_tot; // Density of total positivity shift (other species).
    };
    struct {
      struct gkyl_positivity_shift_vlasov *shift_op_vlasov;
    };
  };
  
  struct gk_species_moment moms; // Positivity shift diagnostic moments.
  struct gk_species_moment integ_moms; // Integrated moments.
  double *red_integ_diag, *red_integ_diag_global; // Reduced integrated moments.
  gkyl_dynvec integ_diag; // Integrated moments of the positivity shift.
  bool is_first_integ_write_call; // Flag first time writing integ_diag.
  bool integ_diag_file_exists; // Whether the integrated diagnostics file exists.
 
  // Methods chosen at runtime.
  void (*apply_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout);
  void (*deltaf_moms_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_positivity *pos);
  void (*deltaf_integ_moms_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_positivity *pos);
  void (*write_diags_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_positivity *pos, double tm, int frame);
  void (*calc_integrated_diags_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks,
    struct gk_positivity *pos, double tm);
  void (*write_integrated_diags_func)(gkyl_gyrokinetic_app *app, struct gk_species *gks,
    struct gk_positivity *pos);
  // Neutral species methods (MF 2025/10/29: to get rid of when we unify species types).
  void (*apply_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
    struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout);
  void (*deltaf_moms_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
    struct gk_positivity *pos);
  void (*deltaf_integ_moms_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
    struct gk_positivity *pos);
  void (*write_diags_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
    struct gk_positivity *pos, double tm, int frame);
  void (*calc_integrated_diags_func_neut)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
    struct gk_positivity *pos, double tm);
  void (*write_integrated_diags_func_neut)(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
    struct gk_positivity *pos);
};

// Species data.
struct gk_species {
  struct gkyl_gyrokinetic_species info; // Input data.

  struct gkyl_basis basis; // Phase-space basis.

  // Basis on device (points to host basis if running w/o GPU).
  struct gkyl_basis *basis_on_dev; 
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  struct gkyl_range global, global_ext; // Global, global-ext conf-space ranges.

  struct gkyl_comm *comm;   // Communicator object for phase-space arrays.
  int nghost[GKYL_MAX_DIM]; // Number of ghost-cells in each direction.

  struct gkyl_rect_grid grid_vel; // Velocity space grid.
  struct gkyl_range local_vel, local_ext_vel; // Local, local-ext velocity-space ranges.

  struct gkyl_velocity_map *vel_map; // Velocity mapping objects.

  struct gkyl_msgpack_map_elem* io_meta; // Metadata for I/O.
  int io_meta_len; // Number of elements in io_meta.

  struct gkyl_array *f, *f1, *fnew; // Arrays for updates.
  struct gkyl_array *cflrate; // CFL rate in each cell.
  struct gkyl_array *cflrate_ho; // CFL rate in each cell on host-side.

  struct gkyl_array *bc_buffer; // Buffer for BCs (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // Buffers for time independent BCs.

  struct gkyl_array *f_host; // Host copy for IO and initialization.

  struct gkyl_array *gyro_phi; // Gyroaveraged electrostatic potential.
  
  struct gk_species_moment m0; // Computes charge density.
  struct gk_species_moment integ_moms; // Integrated moments.
  struct gk_species_moment *moms; // Diagnostic moments
  double *red_integ_diag, *red_integ_diag_global; // Reduced integrated moments.
  gkyl_dynvec integ_diag; // Integrated moments reduced across grid
  bool is_first_integ_write_call;  // Whether dynvec is being written for the first time.

  struct gkyl_array *fdot_mom_old, *fdot_mom_new; // Moments of f_old and f_new.
  gkyl_dynvec fdot_integ_diag; // Integrated moments of Delta f=f_new - f_old..
  bool is_first_fdot_integ_write_call; // Whether dynvec is being written for the first time.

  struct gkyl_array_integrate* integ_wfsq_op; // Operator to integrate w*f^2.
  double *L2norm_local, *L2norm_global; // L2norm in local MPI process and across the communicator.
  gkyl_dynvec L2norm; // L2 norm.
  bool is_first_L2norm_write_call; // Whether dynvec is being written for the first time.

  int num_periodic_dir; // Number of periodic directions.
  int periodic_dirs[GKYL_MAX_CDIM]; // List of periodic directions.
  bool bc_is_np[GKYL_MAX_CDIM]; // Whether BC is nonperiodic.

  // Boundary conditions on lower/upper edges in each direction.
  struct gkyl_gyrokinetic_bc lower_bc[GKYL_MAX_CDIM], upper_bc[GKYL_MAX_CDIM];
  // gyrokinetic sheath boundary conditions
  struct gkyl_bc_sheath_gyrokinetic *bc_sheath_lo;
  struct gkyl_bc_sheath_gyrokinetic *bc_sheath_up;
  // Pointers to updaters that apply (non-sheath) BC.
  struct gkyl_bc_basic_gyrokinetic *bc_lo[GKYL_MAX_CDIM];
  struct gkyl_bc_basic_gyrokinetic *bc_up[GKYL_MAX_CDIM];
  // Local skin/ghost ranges.
  struct gkyl_range local_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range local_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range local_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range local_upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];
  // Core and SOL ranges for IWL sims.
  struct gkyl_range global_core, global_ext_core, global_sol, global_ext_sol;
  struct gkyl_range local_core, local_ext_core, local_sol, local_ext_sol;
  struct gkyl_range local_lower_skin_par_core, local_lower_ghost_par_core;
  struct gkyl_range local_upper_skin_par_core, local_upper_ghost_par_core;
  struct gkyl_range local_lower_skin_par_sol , local_lower_ghost_par_sol;
  struct gkyl_range local_upper_skin_par_sol , local_upper_ghost_par_sol;
  // GK IWL sims need a core range extended in z, and a TS BC updater.
  struct gkyl_range local_par_ext_core; // Core range extended in parallel direction.
  struct gkyl_bc_twistshift *bc_ts_lo, *bc_ts_up;

  struct gk_proj proj_init; // Projector for initial conditions.

  struct gk_collisionless collisionless; // Collisionless terms.

  struct gk_source src; // Plasma source.

  struct gk_damping damping; // Damping term: -nu(z)*f.

  struct gk_fdot_multiplier fdot_mult; // Function multiplying df/dt.

  struct gk_anomalous_diff anom_diff; // Anomalous diffusion.
  
  struct gk_source_bgk bgk_src; // BGK source.

  // Boundary fluxes used for other solvers and diagnostics.
  struct gk_boundary_fluxes bflux;

  struct gk_lte lte; // Object constructing LTE distributions.

  struct gk_scaling sca; // Recycle react scale.

  struct gk_lbo_collisions lbo; // LBO collisions object.
  struct gk_bgk_collisions bgk; // BGK collisions object.

  struct gk_react react; // Object for reactions with charged species.
  struct gk_react react_neut; // Object for reactions with neutral species.

  struct gk_rad_drag rad; // Radiation object.

  struct gk_positivity positivity; // Positivity enforcing operator.

  // Pointer to various functions selected at runtime.
  double (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
    struct gkyl_array *rhs, struct gkyl_array **bflux_moms);
  double (*rhs_implicit_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt);
  void (*bc_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gkyl_array *f);
  void (*release_func)(const gkyl_gyrokinetic_app* app, const struct gk_species *s);
  void (*step_f_func)(struct gkyl_array* out, double dt, const struct gkyl_array* inp); 
  void (*combine_func)(struct gkyl_array *out, double c1,
    const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
    const struct gkyl_range *rng);
  void (*copy_func)(struct gkyl_array *out, const struct gkyl_array *inp,
    const struct gkyl_range *range);
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_cfl_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);
  void (*write_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks);
  void (*calc_L2norm_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);
  void (*write_L2norm_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks);
  void (*calc_int_mom_dt_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom);

  // Quantities used for FLR model:
  struct gkyl_array *m0_gyroavg; // Gyroaveraged particle density.
  struct gkyl_array *flr_rhoSqD2; // Laplacian weight in FLR operator.
  struct gkyl_array *flr_kSq; // Field multiplying phi in FLR operator.
  struct gkyl_deflated_fem_poisson *flr_op; // Helmholtz solver to invert FLR operator.
  // Pointer to function that performs the gyroaverage.
  void (*gyroaverage)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    struct gkyl_array *field_in, struct gkyl_array *field_gyroavg);

  double *omega_cfl; // Maximum Omega_CFL in this MPI process.
  double *m0_max; // Maximum number density in this MPI process.
};

// Neutral species data.
struct gk_neut_species {
  struct gkyl_gyrokinetic_neut_species info; // Data for neutral species.

  bool is_fluid; // Whether this is a fluid species.
  enum gkyl_field_id field_id; // Type of Field equation (always GKYL_FIELD_NULL).
  enum gkyl_model_id model_id; // Type of Vlasov equation (always GKYL_MODEL_CANONICAL_PB)

  struct gkyl_basis basis; // Species basis.
  struct gkyl_basis *basis_on_dev; // Species basis on device (or host if use_gpu=false).

  struct gkyl_job_pool *job_pool; // Job pool.
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // Local, local-ext species ranges.
  struct gkyl_range global, global_ext; // Global, global-ext species ranges.

  struct gkyl_comm *comm;   // Communicator object for this species.
  int nghost[GKYL_MAX_DIM]; // Number of ghost-cells in each direction

  struct gkyl_msgpack_map_elem* io_meta; // Metadata for I/O.
  int io_meta_len; // Number of elements in io_meta.

  struct gkyl_array *f, *f1, *fnew; // Arrays for updates.
  struct gkyl_array *f_host; // Host array for initialization and I/O.

  struct gkyl_array *cflrate; // CFL rate in each cell.

  struct gk_species_moment *moms; // Compute diagnostic moments.
  struct gk_species_moment integ_moms; // Computes integrated moments.

  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  // Local skin/ghost ranges.
  struct gkyl_range local_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range local_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range local_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range local_upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];

  struct gk_proj proj_init; // Projector for initial conditions.

  struct gk_source src; // External source.

  struct gk_lte lte; // Object needed for the lte equilibrium.
 
  struct gk_scaling sca; // Recycle react scale.

  int num_periodic_dir; // Number of periodic directions.
  int periodic_dirs[GKYL_MAX_CDIM]; // List of periodic directions.
  bool bc_is_np[GKYL_MAX_CDIM]; // Whether BC is nonperiodic.
    
  bool enforce_positivity; // Enforces positivity of f or mass/energy density.

  double *omega_cfl; // CFL frequency.

  struct gk_bgk_collisions bgk; // BGK collisions object.

  struct gk_positivity positivity; // Positivity enforcing operator.

  struct gk_react react_neut; // Reaction object.

  // Boundary fluxes used for other solvers and diagnostics.
  struct gk_boundary_fluxes bflux;
        
  union {
    // Kinetic neutrals ............................................ //
    struct {
      struct gkyl_rect_grid grid_vel; // Velocity space grid.
      struct gkyl_range local_vel, local_ext_vel; // Local, local-ext velocity-space ranges.

      struct gkyl_velocity_map *vel_map; // Velocity mapping objects.

      struct gkyl_array *g_ij, *gij; // Metric tensor and its conjugate.
      struct gkyl_array *hamil; // Specified hamiltonian function for canonical poisson bracket
      struct gkyl_array *hamil_host; // Host side hamiltonian array for intial projection
      struct gk_species_moment m0; // Computes density.

      struct gkyl_array *bc_buffer; // Buffer for BCs (used by bc_basic)
      struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // Buffers for ghost-buffer BCs.

      // Recycling wall boundaries.
      struct gk_recycle_wall bc_recycle_lo, bc_recycle_up;
      bool recyc_lo, recyc_up;
      
      // Boundary conditions on lower/upper edges in each direction.
      struct gkyl_gyrokinetic_bc lower_bc[GKYL_MAX_CDIM], upper_bc[GKYL_MAX_CDIM];
      // Pointers to updaters that apply BC.
      struct gkyl_bc_basic *bc_lo[GKYL_MAX_CDIM];
      struct gkyl_bc_basic *bc_up[GKYL_MAX_CDIM];

      struct gk_collisionless collisionless; // Collisionless terms.

      // Updater that enforces positivity by shifting f.
      struct gkyl_positivity_shift_vlasov *pos_shift_op;
      struct gkyl_array *ps_delta_m0; // Number density of the positivity shift.
      struct gk_species_moment ps_moms; // Positivity shift diagnostic moments.
      gkyl_dynvec ps_integ_diag; // Integrated moments of the positivity shift.
      bool is_first_ps_integ_write_call; // Flag first time writing ps_integ_diag.
    };

    // Fluid neutrals .............................................. //
    struct {
      int num_moments; // Number of fluid moments.
      struct gkyl_array *prim_var; // Primitive variables (udrift, pressure).
      struct gkyl_array *prim_var_host; // prim_var on the host for I/O.
    };
  };

  // Pointer to various functions selected at runtime.
  double (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms);
  double (*rhs_implicit_func)(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt);
  void (*bc_func)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
    struct gkyl_array *f);
  void (*apply_pos_shift_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);
  void (*step_f_func)(struct gkyl_array* out, double dt, const struct gkyl_array* inp); 
  void (*combine_func)(struct gkyl_array *out, double c1,
    const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
    const struct gkyl_range *rng);
  void (*copy_func)(struct gkyl_array *out, const struct gkyl_array *inp,
    const struct gkyl_range *range);
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);
  void (*write_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);
  void (*report_n_iter_corr_func)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx);
  void (*release_is_static_func)(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s);
  void (*release_func)(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s);
};

// Field data.
struct gk_field {
  struct gkyl_gyrokinetic_field info; // Data for field.

  enum gkyl_gkfield_id gkfield_id;

  bool update_field; // Are we updating the field?.
  bool calc_init_field; // Whether to compute the t=0 field.

  // Arrays for local charge density, global charge density, and global smoothed (in z) charge density.
  struct gkyl_array *rho_c;
  struct gkyl_array *rho_c_global_dg;
  struct gkyl_array *rho_c_global_smooth; 
  struct gkyl_array *phi_fem, *phi_smooth; // Arrays for updates.

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  bool init_phi_pol; // Whether to use the initial user polarization phi.
  struct gkyl_array *phi_pol; // Initial polarization density potential.

  struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
                                      // for solving subset of Poisson solves with parallelization in z

  // organization of the different equation objects and the required data and solvers
  union {
    struct {
      struct gkyl_ambi_bolt_potential *ambi_pot;
      struct gkyl_array *sheath_vals[2*GKYL_MAX_CDIM];
    };
    // EM GK model
    struct {
      struct gkyl_array *apar_fem; // array for A_parallel
      struct gkyl_array *apardot_fem; // array for d/dt A_parallel
    };
  };

  double es_energy_fac_1d; 
  struct gkyl_array *es_energy_fac; 
  bool is_dirichletvar; // Whether user provided spatially varying phi BCs.
  struct gkyl_array *phi_bc; // Spatially varying BC.
  struct gkyl_array *epsilon; // Polarization weight including geometric factors.
  struct gkyl_array *kSq; // Weight for Poisson/Helmholtz solver. 

  struct gkyl_poisson_bias_line_list fem_parproj_bias_line_list; // Biased lines constraining the solution.

  struct gkyl_fem_parproj *fem_parproj; // Projects DG function onto continuous  FEM basis:  weight*phi_{fem} = phi_{dg} 
  struct gkyl_fem_parproj *fem_parproj_rho_sol , *fem_parproj_phi_sol ; // FEM projection in the SOL.
  struct gkyl_fem_parproj *fem_parproj_rho_core, *fem_parproj_phi_core; // FEM projection in the core.

  struct gkyl_fem_poisson_perp *fem_poisson_perp; // Solves - nabla . (epsilon * nabla phi) - kSq * phi = rho.

  // Objects needed for FLR effects.
  bool use_flr; // Whether to apply FLR effects.
  void (*invert_flr)(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi); // Function inverting  FLR  operator.
  struct gkyl_array *flr_rhoSq_sum; // Laplacian weight in FLR operator.
  struct gkyl_array *flr_kSq; // Field multiplying phi in FLR operator.
  struct gkyl_deflated_fem_poisson *flr_op; // Helmholtz solver to invert FLR operator.

  struct gkyl_array_integrate *calc_em_energy; // Operator computing EM energy.
  double *em_energy_red, *em_energy_red_global; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for energy dynvec written first time

  double *em_energy_red_old, *em_energy_red_new; // memory for use in GPU reduction of old EM energy.
  gkyl_dynvec integ_energy_dot; // d/dt of integrated energy components.
  bool is_first_energy_dot_write_call; // flag for d(energy)/dt dynvec written first time

  bool has_phi_wall_lo; // flag to indicate there is biased wall potential on lower wall
  bool phi_wall_lo_evolve; // flag to indicate biased wall potential on lower wall is time dependent
  struct gkyl_array *phi_wall_lo; // biased wall potential on lower wall
  struct gkyl_array *phi_wall_lo_host; // host copy for use in IO and projecting
  gkyl_eval_on_nodes *phi_wall_lo_proj; // projector for biased wall potential on lower wall 

  bool has_phi_wall_up; // flag to indicate there is biased wall potential on upper wall
  bool phi_wall_up_evolve; // flag to indicate biased wall potential on upper wall is time dependent
  struct gkyl_array *phi_wall_up; // biased wall potential on upper wall
  struct gkyl_array *phi_wall_up_host; // host copy for use in IO and projecting
  gkyl_eval_on_nodes *phi_wall_up_proj; // projector for biased wall potential on upper wall 

  // Pointer to function that computes the time rate of change of the energy.
  void (*calc_energy_dt_func)(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced);

  // Objects used in IWL simulations and TS BCs.
  struct gkyl_bc_twistshift *bc_ts_lo; // Fills lower core z-ghost with TS BC.
  struct gkyl_bc_basic_gyrokinetic *gfss_bc_op_core_up; // Fills upper core  z-ghost with skin  boundary value.
  struct gkyl_array *bc_buffer; // Buffer for bc_basic.
  
  // Pointer to functions that make phi continuous along z.
  void (*fem_projection_par_rho_func)(gkyl_gyrokinetic_app *app, struct gk_field *field,
    struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);
  void (*fem_projection_par_phi_func)(gkyl_gyrokinetic_app *app, struct gk_field *field,
    struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);
  // Pointer to function to calculate the potential.
  void (*rhs_phi_func)(struct gkyl_gyrokinetic_app *app, struct gk_field *field);
  // Pointer to function that calculates the field energy.
  void (*calc_energy_func)(gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm);
  // Pointer to function that accumulates the charge density. 
  void (*accumulate_rhoc_func)(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);
  // Pointer to function that frees memory.
  void (*release_func)(const struct gkyl_gyrokinetic_app* app, struct gk_field *field);
};

// Eirene data
struct gk_eirene {
  struct gkyl_gyrokinetic_eirene info; // data for coupling
  struct gk_species *coupling_species[GKYL_MAX_SPECIES]; // pointers to species to couple
  struct gk_source_bgk bgk_src[GKYL_MAX_SPECIES];
};

// Gyrokinetic object: used as opaque pointer in user code.
struct gkyl_gyrokinetic_app {
  char name[128]; // Name of app.
  struct gkyl_job_pool *job_pool; // Job pool.
  
  int cdim; // Configuration space dimensions.
  int poly_order; // Polynomial order.
  double tcurr; // Current time.
  double cfl; // CFL number.
  double cfl_omegaH; // CFL number used for omega_H.
  double bmag_ref; // Reference magnetic field.

  bool use_gpu; // Should we use GPU (if present).

  int num_periodic_dir; // Number of periodic directions.
  int periodic_dirs[3]; // List of periodic directions.
    
  struct gkyl_rect_grid grid; // Config-space grid.
  struct gkyl_range local, local_ext; // Local, local-ext conf-space ranges.
  struct gkyl_range global, global_ext; // Global, global-ext conf-space ranges.
  // To simplify BC application, store local skin and ghost ranges.
  struct gkyl_range local_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range local_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range local_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range local_upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges, valid (i.e. volume>0) in ranks abutting boundaries.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];

  // Core and SOL ranges for sims with a LCFS.
  struct gkyl_range global_core, global_ext_core, global_sol, global_ext_sol;
  struct gkyl_range global_par_ext_core; // Core range extended in parallel direction.
  struct gkyl_range global_lower_skin_par_core, global_lower_ghost_par_core;
  struct gkyl_range global_upper_skin_par_core, global_upper_ghost_par_core;
  struct gkyl_range global_lower_skin_par_sol , global_lower_ghost_par_sol;
  struct gkyl_range global_upper_skin_par_sol , global_upper_ghost_par_sol;
  struct gkyl_range local_core, local_ext_core, local_sol, local_ext_sol;
  struct gkyl_range local_par_ext_core; // Core range extended in parallel direction.
  struct gkyl_range local_lower_skin_par_core, local_lower_ghost_par_core;
  struct gkyl_range local_upper_skin_par_core, local_upper_ghost_par_core;
  struct gkyl_range local_lower_skin_par_sol , local_lower_ghost_par_sol;
  struct gkyl_range local_upper_skin_par_sol , local_upper_ghost_par_sol;

  struct gkyl_basis basis; // conf-space basis
  
  struct gkyl_rect_decomp *decomp; // Volume decomposition object.
  struct gkyl_comm *comm; // Volume communicator object for conf-space arrays.

  struct gkyl_rect_decomp *decomp_plane[GKYL_MAX_CDIM]; // Plane/slab decomposition object.
  struct gkyl_comm *comm_plane[GKYL_MAX_CDIM]; // Plane/slab communicator object for conf-space arrays.

  // pointer to basis on device
  // (points to host structs if not on GPU)
  struct gkyl_basis *basis_on_dev; 

  struct gk_geometry *gk_geom;
  struct gkyl_dg_geom *dg_geom;
  struct gkyl_gk_dg_geom *gk_dg_geom;
  struct gkyl_array *jacobtot_inv_weak; // 1/(J.B) computed via weak mul and div.
  double omegaH_gf; // Geometry and field model dependent part of omega_H.
  
  struct gkyl_position_map *position_map; // Position mapping object.

  struct gk_field *field; // pointer to field object
  // Pointer to function that computes the fields.
  void (*calc_field_func)(gkyl_gyrokinetic_app* app, double tcurr,
    const struct gkyl_array *fin[], struct gkyl_array **bflux[]);

  struct gk_eirene *eirene; // EIRENE coupling object.

  int num_species; // Number of charged species.
  struct gk_species *species; // Data for each charged species.

  int num_neut_species; // Number of neutral species.
  struct gk_neut_species *neut_species; // Data for each neutral species.

  // Enforce quasineutrality of sum_s delta f_s after enforcing positivity of charged species.
  bool post_positivity_quasineut;
  struct gkyl_array *ps_delta_m0_ions; // Number density of the total ion positivity shift.
  struct gkyl_array *ps_delta_m0_elcs; // Number density of the total elc positivity shift.
  void (*post_pos_quasineut_func)(gkyl_gyrokinetic_app *app, struct gkyl_array *fout[]);
  
  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_gyrokinetic_app *app, double dt0);

  struct gkyl_gyrokinetic_stat stat; // Statistics.

  struct gkyl_msgpack_map_elem* io_meta_basic; // Basic metadata for I/O.
  int io_meta_basic_len; // Number of elements in io_meta_basic.
  struct gkyl_msgpack_map_elem* io_meta; // Metadata for I/O.
  int io_meta_len; // Number of elements in io_meta.

  gkyl_dynvec dts; // Record time step over time.
  bool is_first_dt_write_call; // Flag for integrated moments dynvec written first time.
};

/** gkyl_gyrokinetic_app private API */

static inline struct gkyl_gyrokinetic_bc *
gk_fetch_bc_with_dir_edge(struct gkyl_gyrokinetic_bc *bc_list, int num_bcs,
  int dir, enum gkyl_edge_loc edge)
{
  // Given a list of gyrokinetic BCs `bc_list` of length `num_bcs`, return the
  // entry corresponding to the BC in the `dir` direction and the `edge` edge.
  // It omits GKYL_BC_GK_SKIP BCs.
  struct gkyl_gyrokinetic_bc *out = 0;
  for (int i=0; i<num_bcs; ++i) {
    struct gkyl_gyrokinetic_bc *bc_curr = &bc_list[i];
    if (bc_curr->dir == dir && bc_curr->edge == edge && bc_curr->type != GKYL_BC_GK_SKIP) {
      out = bc_curr;
      break;
    }
  }
  return out;
}

/**
 * Allocate a new gyrokinetic app and initialize its conf-space grid and
 * geometry. This method needs to be complemented by
 * gkyl_gyrokinetic_app_new_solver below.
 *
 * @param gk Gyrokinetic input struct.
 * @return A gyrokinetic app object.
 */
gkyl_gyrokinetic_app*
gkyl_gyrokinetic_app_new_geom(struct gkyl_gk *gk);

/**
 * Initialize the rest of the gyrokinetic app solver, after having called
 * the gkyl_gyrokinetic_app_new_geom method.
 *
 * @param gk Gyrokinetic input struct.
 * @param app Gyrokinetic app.
 */
void
gkyl_gyrokinetic_app_new_solver(struct gkyl_gk *gk, gkyl_gyrokinetic_app *app);

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct gk_species* gk_find_species(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into.
 * @param nm Name of species.
 * @return Index of species, -1 if not found.
 */
int gk_find_species_idx(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Find neutral species with given name.
 *
 * @param app Top-level app to look into.
 * @param nm Name of neutral species.
 * @return Pointer to neutral species with given name. NULL if not found.
 */
struct gk_neut_species* gk_find_neut_species(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Return index of neutral species in the order it appears in the input.
 *
 * @param app Top-level app to look into.
 * @param nm Name of neutral species.
 * @return Index of neutral species, -1 if not found.
 */
int gk_find_neut_species_idx(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_app_apply_ic_species(gkyl_gyrokinetic_app* app, int sidx, double t0);

/**
 * Initialize neutral species by projecting initial conditions on basis
 * functions. Neutral species index (sidx) is the same index used to specify
 * the neutral species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_app_apply_ic_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0);

/**
 * Perform part of initialization that depends on the other species being
 * (partially) initialized.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_app_apply_ic_cross_species(gkyl_gyrokinetic_app* app, int sidx, double t0);

/**
 * Perform part of neutral species initialization that depends on the other species being
 * (partially) initialized.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_app_apply_ic_cross_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0);

/** gk_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param sm Species moment object.
 * @param nm Name string indicating moment type
 * @param is_integrated Whether to compute the volume integrated moment.
 */
void gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated);

/**
 * Calculate moment, given distribution function @a fin.
 * 
 * @param sm Species moment object.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param fin Input distribution function array
 */
void gk_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Divide the moment by the configuration-space Jacobian for diagnostic
 * purposes. Not all components of the diagnostic moment are divided, depending
 * on the type.
 * 
 * @param app Gyrokinetic app object.
 * @param sm Species moment object.
 * @param Jmom_in Moment(s) to be divided by J.
 * @param mom_out Array in which to place the output.
 */
void gk_species_moment_diag_jacobgeo_div(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out);

/**
 * Release species moment object.
 *
 * @param app gyrokinetic app object
 * @param sm Species moment object to release.
 */
void gk_species_moment_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species_moment *sm);

/** gk_positivity API */

/**
 * Initialize species positivity operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param pos Positivity object.
 */
void gk_species_positivity_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_positivity *pos);

void gk_neut_species_positivity_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_positivity *pos);

/**
 * Determine how many species participate in the quasineutrality enforcement
 * after applying the positivity operator.
 *
 * @param app Gyrokinetic app object.
 */
int gk_species_positivity_num_species_in_quasineut(gkyl_gyrokinetic_app* app);

/**
 * Apply the positivity operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param pos Positivity object.
 * @param fbuffer Buffer array where delta f will be held.
 * @param fout State array to apply the positivity operation to.
 */
void gk_species_positivity_apply(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout);

void gk_neut_species_positivity_apply(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout);

/**
 * Write conf-space diagnostics of the positivity operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param pos Positivity object.
 * @param tm Current simulation time.
 * @param frame Current I/O frame.
 */
void gk_species_positivity_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_positivity *pos, double tm, int frame);

void gk_neut_species_positivity_write_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, double tm, int frame);

/**
 * Calculate integrated diagnostics of the positivity operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param pos Positivity object.
 * @param tm Current simulation time.
 */
void gk_species_positivity_calc_integrated_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_positivity *pos, double tm);

void gk_neut_species_positivity_calc_integrated_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, double tm);

/**
 * Write integrated diagnostics of the positivity operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param pos Positivity object.
 */
void gk_species_positivity_write_integrated_diags(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_positivity *pos);

void gk_neut_species_positivity_write_integrated_diags(gkyl_gyrokinetic_app *app,
  struct gk_neut_species *gkns, struct gk_positivity *pos);

/**
 * Release species positivity operator.
 *
 * @param app Gyrokinetic app object.
 * @param pos Positivity object.
 */
void gk_species_positivity_release(const struct gkyl_gyrokinetic_app *app, const struct gk_positivity *pos);

void gk_neut_species_positivity_release(const struct gkyl_gyrokinetic_app *app, const struct gk_positivity *pos);

/**
 * Reset the positivity operator.
 *
 * @param app gyrokinetic app object.
 * @param tm Time-stamp.
 * @param gks Species object to delete.
 * @param pos Positivity object.
 * @param pos_inp Positivity input parameters.
 */
void gk_species_positivity_reset(gkyl_gyrokinetic_app* app, double tm,
  struct gk_species *gks, struct gk_positivity *pos, struct gkyl_gyrokinetic_positivity pos_inp);

void gk_neut_species_positivity_reset(gkyl_gyrokinetic_app* app, double tm,
  struct gk_neut_species *gkns, struct gk_positivity *pos, struct gkyl_gyrokinetic_positivity pos_inp);

/** gk_species_lte API */

/**
 * Initialize species lte object.
 *
 * @param app Gyrokinetic app object
 * @param s Species object 
 * @param lte Species lte object
 * @param corr_inp Input struct with moment correction inputs
 */
void gk_species_lte_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_lte *lte, struct correct_all_moms_inp corr_inp);

/**
 * Compute LTE distribution from input moments
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to species
 * @param lte Pointer to lte object
 * @param moms_lte Input LTE moments
 */
void gk_species_lte_from_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lte *lte, const struct gkyl_array *moms_lte);

/**
 * Compute equivalent LTE distribution from input distribution function. 
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to species
 * @param lte Pointer to lte
 * @param fin Input distribution function
 */
void gk_species_lte(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lte *lte, const struct gkyl_array *fin);

/**
 * Write the LTE correction status. 
 *
 * @param app Gyrokinetic app object
 * @param gks Pointer to species
 */
void gk_species_lte_write_max_corr_status(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Release species lte object.
 *
 * @param app gyrokinetic app object
 * @param lte Species lte object to release
 */
void gk_species_lte_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte);

/** gk_species_collisionless API */

/**
 * Initialize species collisionless object.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param gkcls Species collisionless object.
 */
void gk_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls);

/**
 * Compute the collisionless phase-space flux.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param gkcls Species collisionless object.
 * @param fin Input distribution function.
 */
void gk_species_collisionless_flux(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin);

/**
 * Compute RHS contribution from collisionless terms.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param gkcls Species collisionless object.
 * @param fin Input distribution function.
 * @param rhs collisionless contribution to df/dt.
 */
void gk_species_collisionless_rhs(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write out diagnostics from the collisionless terms.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param gkcls Species collisionless object.
 * @param tm Current simulation time.
 * @param frame Current I/O frame.
 */
void gk_species_collisionless_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame);

/**
 * Release species collisionless object.
 *
 * @param app Gyrokinetic app object.
 * @param gkcls Species collisionless object.
 */
void gk_species_collisionless_release(const struct gkyl_gyrokinetic_app *app, const struct gk_collisionless *gkcls);

/**
 * Reset the collisionless solver.
 *
 * @param app gyrokinetic app object.
 * @param tm Time-stamp.
 * @param gks Species object to delete.
 * @param gkcls Input parameters for collisionless solver.
 */
void gk_species_collisionless_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks, 
  struct gk_collisionless *gkcls, struct gkyl_gyrokinetic_collisionless gkcls_inp);

/** gk_neut_species_collisionless API */

/**
 * Initialize species collisionless object.
 *
 * @param app Gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param gkcls Species collisionless object.
 */
void gk_neut_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls);

/**
 * Compute the collisionless phase-space flux.
 *
 * @param app Gyrokinetic app object.
 * @param gkns Pointer to neutral species.
 * @param gkcls Species collisionless object.
 * @param fin Input distribution function.
 */
void gk_neut_species_collisionless_flux(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin);

/**
 * Compute RHS contribution from collisionless terms.
 *
 * @param app Gyrokinetic app object.
 * @param gkns Pointer to neutral species.
 * @param gkcls Species collisionless object.
 * @param fin Input distribution function.
 * @param rhs collisionless contribution to df/dt.
 */
void gk_neut_species_collisionless_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write out diagnostics from the collisionless terms.
 *
 * @param app Gyrokinetic app object.
 * @param gkns Pointer to neutral species.
 * @param gkcls Species collisionless object.
 * @param tm Current simulation time.
 * @param frame Current I/O frame.
 */
void gk_neut_species_collisionless_write_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame);

/**
 * Release species collisionless object.
 *
 * @param app Gyrokinetic app object.
 * @param gkcls Species collisionless object.
 */
void gk_neut_species_collisionless_release(const struct gkyl_gyrokinetic_app *app, const struct gk_collisionless *gkcls);

/** gk_species_radiation API */

/**
 * Initialize species radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param rad Species radiation drag object
 */
void gk_species_radiation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_rad_drag *rad);

/**
 * Compute necessary moments for radiation drag object
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution functions (size num_species)
 * @param fin_neut Input neutral distribution functions (size num_species)
 */
void gk_species_radiation_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute RHS from radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void gk_species_radiation_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write species radiation drag.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_radiation_write_drag(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Write species radiation emissivity.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_radiation_write_emissivity(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Calculate species radiation integrated moments.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 */
void gk_species_radiation_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Write species radiation integrated moments.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 */
void gk_species_radiation_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Release species radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param rad Species radiation drag object to release
 */
void gk_species_radiation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_rad_drag *rad);

/** gk_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param lbo Species LBO object.
 */
void gk_species_lbo_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param lbo Species LBO object.
 */
void gk_species_lbo_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param lbo Pointer to LBO.
 * @param fin Input distribution function.
 */
void gk_species_lbo_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param lbo Pointer to LBO.
 * @param fin Input distribution function.
 */
void gk_species_lbo_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin);

/**
 * Compute RHS from LBO collisions.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param lbo Pointer to LBO.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from LBO.
 */
void gk_species_lbo_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write moments from LBO object.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void gk_species_lbo_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Release species LBO object.
 *
 * @param app gyrokinetic app object.
 * @param lbo Species LBO object to release.
 */
void gk_species_lbo_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lbo_collisions *lbo);

/** gk_species_bgk API */

/**
 * Initialize species BGK collisions object.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param bgk Species BGK object.
 */
void gk_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk);

/**
 * Initialize species BGK cross-collisions object.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param bgk Species BGK object.
 */
void gk_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk);

/**
 * Compute necessary moments and boundary
 * corrections for BGK collisions.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 */
void gk_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

void gk_species_bgk_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species BGK collisions.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 */
void gk_species_bgk_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

void gk_species_bgk_cross_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

/**
 * Compute RHS from BGK collisions (explicit integrator).
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from BGK.
 */
void gk_species_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Compute RHS from BGK collisions (implicit integrator).
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 * @param dt Time step.
 * @param rhs On output, the RHS from BGK.
 */
void gk_species_bgk_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);

/**
 * Write moments from BGK object.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void gk_species_bgk_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Release species BGK object.
 *
 * @param app gyrokinetic app object.
 * @param bgk Species BGK object to release.
 */
void gk_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk);

/** gk_species_react API */

/**
 * Initialize species reactions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param inp Input reaction struct for determining types of reactions
 * @param react Species reaction object
 * @param all_gk Boolean for if the reactions are between only GK species
 */
void gk_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gkyl_gyrokinetic_react inp, struct gk_react *react, bool all_gk);

/**
 * Initialize species reactions "cross-collisions" object
 * for who is reacting with whom
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param react Species react object
 */
void gk_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_react *react);

/**
 * Compute necessary rates and moments for reactions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param react Pointer to react
 * @param f_self Input self distribution function
 * @param fin Input distribution functions (size: num_species)
 * @param fin_neut Input neutral distribution functions (size: num_neut_species)
 */
void gk_species_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute RHS from reactions 
 * (e.g., ionization, charge exchange, recombination, or radiation)
 *
 * @param app gyrokinetic app object
 * @param s Pointer to species
 * @param react Pointer to react
 * @param fin Input distribution function
 * @param rhs On output, the RHS from react (df/dt)
 */
void gk_species_react_rhs(gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write reaction rate.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param gkr Pointer to react object
 * @param ridx Index for reaction species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_react_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, struct gk_react *gkr,
  int ridx, double tm, int frame);

/**
 * Release species react object.
 *
 * @param app gyrokinetic app object
 * @param react Species react object to release
 */
void gk_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react);

/** gk_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object. 
 * @param bflux Species boundary flux object.
 * @param bflux_type Indicate whether to compute fluxes, flux moments or diagnostics.
 * @param add_moms_inp Additional moments to step in time if bflux_type=GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS.
 */
void gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, void *species,
  struct gk_boundary_fluxes *bflux, enum gkyl_species_bflux_type bflux_type,
  struct gkyl_phase_diagnostics_inp add_moms_inp);

/**
 * In case time-integrated volume-integrated boundary flux diagnostics are
 * requested, read the final volume-time integrated diagnostics from the previous simulation.
 *
 * @param app Gyrokinetic app object.
 * @param spec_in Species object. 
 * @param bflux Species boundary flux object.
 */
void
gk_species_bflux_read_voltime_integrated_mom(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux);

/**
 * Compute boundary flux, either for another solver or for diagnostics.
 * Note: stores the boundary flux in the ghost cells of rhs.
 * The ghost cells are overwritten by apply_bc.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param fin Input distribution function.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 */
void gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Like gk_species_bflux_rhs but actually calculates boundary fluxes, unlike
 * gk_species_bflux_rhs which may be empty.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param fin Input distribution function.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 */
void gk_species_bflux_rhs_calc(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Copy the boundary fluxes into a given range of a given phase-space array.
 *
 * @param bflux Species boundary flux object.
 * @param dir Direction of the boundary.
 * @param edge Edge of the boundary.
 * @param out Array to copy the boundary flux into.
 * @param out_rng Range to copy the boundary flux into.
 */
void
gk_species_bflux_get_flux(struct gk_boundary_fluxes *bflux, int dir,
  enum gkyl_edge_loc edge, struct gkyl_array *out, const struct gkyl_range *out_rng);

/**
 * Copy a boundary fluxe moment into a given range of a given conf-space array.
 *
 * @param bflux Species boundary flux object.
 * @param dir Direction of the boundary.
 * @param edge Edge of the boundary.
 * @param mom_type Name of the moment desired.
 * @param bflux_moms Array of moments of boundary fluxes through every boundary.
 * @param out Array to copy the boundary flux moment into.
 * @param out_rng Range to copy the boundary flux moment into.
 */
void
gk_species_bflux_get_flux_mom(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  enum gkyl_distribution_moments mom_type, struct gkyl_array **bflux_moms, struct gkyl_array *out,
  const struct gkyl_range *out_rng);

/**
 * Compute moments of the boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 * @param bflux_moms Array of moments of boundary fluxes through every boundary.
 */
void gk_species_bflux_calc_moms(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *rhs, struct gkyl_array **bflux_moms);

/**
 * Clear the boundary fluxes at each boundary.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_in Array of boundary fluxes to clear.
 * @param val Value to set array to.
 */
void
gk_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val);

/**
 * Scale the boundary fluxes at each boundary.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_in Array of boundary fluxes to clear.
 * @param val Value to scale the fluxes by.
 */
void
gk_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val);

/**
 * Step the diagnotic boundary fluxes forward once.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param dt Time step.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_species_bflux_step_f(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double dt, struct gkyl_array **bflux_in);

/**
 * Set a stage of the boundary fluxes equal to another.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac Multiplicative factor.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_species_bflux_set(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double fac, struct gkyl_array **bflux_in);

/**
 * Accumulate a stage of the boundary fluxes onto another.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac Multiplicative factor.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_species_bflux_accumulate(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double fac, struct gkyl_array **bflux_in);

/**
 * Combine the diagnotic boundary fluxes for multi-stage RK stepper.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac1 Factor to multiply bflux_in1 by.
 * @param bflux_in1 Array of input boundary fluxes.
 * @param fac2 Factor to multiply bflux_in2 by.
 * @param bflux_in2 Array of input boundary fluxes.
 */
void
gk_species_bflux_combine(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2);

/**
 * Copy diagnotic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_species_bflux_copy(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_fout, struct gkyl_array **bflux_in);

/**
 * Calculate the integrated moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object
 * @param species Species object.
 * @param tm Current simulation time.
 */
void
gk_species_bflux_calc_integrated_mom(gkyl_gyrokinetic_app* app, void *species,
  struct gk_boundary_fluxes *bflux, double tm);

/**
 * Calculate the time integrated, integrated moments
 * of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object
 * @param species Species object.
 * @param tm Current simulation time.
 */
void
gk_species_bflux_calc_voltime_integrated_mom(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm);

/**
 * Write the integrated moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 */
void
gk_species_bflux_write_integrated_mom(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux);

/**
 * Write the moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 * @param tm Current simulation time.
 * @param frame Current frame.
 */
void
gk_species_bflux_write_mom(gkyl_gyrokinetic_app *app, void *species,
  struct gk_boundary_fluxes *bflux, double tm, int frame);

/**
 * Release species boundary flux object.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 * @param bflux Species boundary flux object to release.
 */
void gk_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const void *species, 
  const struct gk_boundary_fluxes *bflux);

/** gk_neut_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object. 
 * @param bflux Species boundary flux object.
 * @param bflux_type Indicate whether to compute fluxes, flux moments or diagnostics.
 * @param add_moms_inp Additional moments to step in time if bflux_type=gk_neut_species_BFLUX_CALC_FLUX_STEP_MOMS.
 */
void gk_neut_species_bflux_init(struct gkyl_gyrokinetic_app *app, void *species,
  struct gk_boundary_fluxes *bflux, enum gkyl_species_bflux_type bflux_type,
  struct gkyl_phase_diagnostics_inp add_moms_inp);

/**
 * In case time-integrated volume-integrated boundary flux diagnostics are
 * requested, read the final volume-time integrated diagnostics from the previous simulation.
 *
 * @param app Gyrokinetic app object.
 * @param spec_in Species object. 
 * @param bflux Species boundary flux object.
 */
void
gk_neut_species_bflux_read_voltime_integrated_mom(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux);

/**
 * Compute boundary flux, either for another solver or for diagnostics.
 * Note: stores the boundary flux in the ghost cells of rhs.
 * The ghost cells are overwritten by apply_bc.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param fin Input distribution function.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 */
void gk_neut_species_bflux_rhs(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Like gk_neut_species_bflux_rhs but actually calculates boundary fluxes, unlike
 * gk_neut_species_bflux_rhs which may be empty.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param fin Input distribution function.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 */
void gk_neut_species_bflux_rhs_calc(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Copy the boundary fluxes into a given range of a given phase-space array.
 *
 * @param bflux Species boundary flux object.
 * @param dir Direction of the boundary.
 * @param edge Edge of the boundary.
 * @param out Array to copy the boundary flux into.
 * @param out_rng Range to copy the boundary flux into.
 */
void
gk_neut_species_bflux_get_flux(struct gk_boundary_fluxes *bflux, int dir,
  enum gkyl_edge_loc edge, struct gkyl_array *out, const struct gkyl_range *out_rng);

/**
 * Copy a boundary fluxe moment into a given range of a given conf-space array.
 *
 * @param bflux Species boundary flux object.
 * @param dir Direction of the boundary.
 * @param edge Edge of the boundary.
 * @param mom_name Name of the moment desired.
 * @param bflux_moms Array of moments of boundary fluxes through every boundary.
 * @param out Array to copy the boundary flux moment into.
 * @param out_rng Range to copy the boundary flux moment into.
 */
void
gk_neut_species_bflux_get_flux_mom(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  enum gkyl_distribution_moments mom_type, struct gkyl_array **bflux_moms, struct gkyl_array *out,
  const struct gkyl_range *out_rng);

/**
 * Compute moments of the boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 * @param bflux_out Array of moments of boundary fluxes through every boundary.
 */
void gk_neut_species_bflux_calc_moms(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *rhs, struct gkyl_array **bflux_moms);

/**
 * Clear the boundary fluxes at each boundary.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_in Array of boundary fluxes to clear.
 * @param val Value to set array to.
 */
void
gk_neut_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val);

/**
 * Scale the boundary fluxes at each boundary.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_in Array of boundary fluxes to clear.
 * @param val Value to scale the fluxes by.
 */
void
gk_neut_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val);

/**
 * Step the diagnotic boundary fluxes forward once.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param dt Time step.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_neut_species_bflux_step_f(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double dt, struct gkyl_array **bflux_in);

/**
 * Set a stage of the boundary fluxes equal to another.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Neutral species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac Multiplicative factor.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_neut_species_bflux_set(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double fac, struct gkyl_array **bflux_in);

/**
 * Accumulate a stage of the boundary fluxes onto another.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Neutral species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac Multiplicative factor.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_neut_species_bflux_accumulate(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double fac, struct gkyl_array **bflux_in);

/**
 * Combine the diagnotic boundary fluxes for multi-stage RK stepper.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac1 Factor to multiply bflux_in1 by.
 * @param bflux_in1 Array of input boundary fluxes.
 * @param fac2 Factor to multiply bflux_in2 by.
 * @param bflux_in2 Array of input boundary fluxes.
 */
void
gk_neut_species_bflux_combine(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2);

/**
 * Copy diagnotic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_neut_species_bflux_copy(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_fout, struct gkyl_array **bflux_in);

/**
 * Calculate the integrated moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object
 * @param species Species object.
 * @param tm Current simulation time.
 */
void
gk_neut_species_bflux_calc_integrated_mom(gkyl_gyrokinetic_app* app, void *species,
  struct gk_boundary_fluxes *bflux, double tm);

/**
 * Calculate the time integrated, integrated moments
 * of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object
 * @param species Species object.
 * @param tm Current simulation time.
 */
void
gk_neut_species_bflux_calc_voltime_integrated_mom(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm);

/**
 * Write the integrated moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 */
void
gk_neut_species_bflux_write_integrated_mom(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux);

/**
 * Write the moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 * @param tm Current simulation time.
 * @param frame Current frame.
 */
void
gk_neut_species_bflux_write_mom(gkyl_gyrokinetic_app *app, void *species,
  struct gk_boundary_fluxes *bflux, double tm, int frame);

/**
 * Release species boundary flux object.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 * @param bflux Species boundary flux object to release.
 */
void gk_neut_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const void *species, 
  const struct gk_boundary_fluxes *bflux);

/** gk_species_projection API */

/**
 * Initialize species projection object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param inp Input struct for projection (contains functions pointers for type of projection).
 * @param proj Species projection object.
 */
void gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj);

/**
 * Compute species projection.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param proj Species projection object.
 * @param f Output distribution function from projection.
 * @param tm Time for use in projection.
 */
void gk_species_projection_calc(gkyl_gyrokinetic_app *app, struct gk_species *species, 
  struct gk_proj *proj, struct gkyl_array *f, double tm);

/**
 * Release species projection object.
 *
 * @param app gyrokinetic app object.
 * @param proj Species projection object to release.
 */
void gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj);

/** gk_species_source API */

/**
 * Initialize species source object.
 *
 * @param app Gyrokinetic app object.
 * @param s Species object.
 * @param src Species source object.
 */
void gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_source *src);

/**
 * Compute species applied source term.
 *
 * @param app Gyrokinetic app object.
 * @param species Species object.
 * @param src Species source object.
 * @param f_buffer Phase-space buffer used to project the source.
 * @param tm Time for use in source.
 */
void gk_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_species *species, 
  struct gk_source *src, struct gkyl_array *f_buffer, double tm);

/**
 * Adapt source to maintain input power and/or particle content constant.
 * 
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param src Species source object.
 * @param f_buffer Phase-space buffer used to accumulate the source.
 * @param bflux_moms Moments of boundary fluxes passed in through the gyrokinetic rhs routine.
 * @param tm Time for use in source.
 */
void
gk_species_source_adapt(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, 
  struct gkyl_array **bflux_moms[], double tm);

/**
 * Compute RHS contribution from source.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param src Pointer to source.
 * @param fin Input distribution function.
 * @param rhs Source contribution to df/dt.
 */
void gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write source diagnostics.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void gk_species_source_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Write source moment diagnostics.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void gk_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Calc source integrated moment diagnostics.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for source diagnostic.
 */
void gk_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Write source integrated moment diagnostics.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 */
void gk_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Release species source object.
 *
 * @param app Gyrokinetic app object.
 * @param src Species source object to release.
 */
void gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_species_damping API */

/**
 * Initialize species damping object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param damp Species damping object.
 */
void gk_species_damping_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_damping *damp);

/**
 * Initialize filtered distribution function from initial distribution.
 * Should be called after initial conditions are applied to f.
 *
 * @param gks Pointer to species.
 * @param damp Species damping object.
 * @param f Initial distribution function.
 */
void gk_species_damping_set_fbar_to_f(const struct gk_species *gks, struct gk_damping *damp, 
  const struct gkyl_array *f);


/**
 * Compute species applied source term.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param damp Species damping object.
 * @param phi Current electrostatic potential.
 * @param fin Current distribution function.
 * @param fbar_in Filtered distribution corresponding to fin RK stage.
 * @param f_buffer Phase-space buffer.
 * @param rhs df/dt damping term gets added to.
 * @param cflrate CFL frequency in phase space.
 */
void gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate);

/**
 * Write damping diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for damping diagnostic.
 * @param frame Output frame.
 */
void gk_species_damping_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Compute filtered distribution RHS for low-pass filter damping:
 *   d(fbar)/dt = rate * (f - fbar).
 *
 * @param damp Species damping object.
 * @param fin Current distribution function f.
 * @param fbar_in Current filtered distribution function fbar.
 * @param rhs_fbar RHS buffer for fbar.
 */
void gk_species_damping_calc_fbar_rhs(const struct gk_damping *damp,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in, struct gkyl_array *rhs_fbar);

/**
 * Take a forward-Euler substep for the filtered distribution.
 *
 * This is a runtime-dispatched no-op unless the damping type enables
 * filtered-state evolution (currently low-pass filter).
 *
 * @param gks Species object.
 * @param fin Current distribution function f.
 * @param fbar_in Current filtered distribution state.
 * @param fbar_out Output filtered distribution state after FE substep.
 * @param dt Time-step for substep.
 */
void gk_species_damping_forward_euler(struct gk_species *gks,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *fbar_out, double dt);

/**
 * Combine filtered distributions like gk_species_combine.
 *
 * This is runtime-dispatched and a no-op unless the damping type enables
 * filtered-state evolution (currently low-pass filter).
 */
void gk_species_damping_combine(struct gk_species *gks, struct gkyl_array *fout, double c1,
  const struct gkyl_array *f1, double c2, const struct gkyl_array *f2,
  const struct gkyl_range *rng);

/**
 * Copy filtered distribution ranges like gk_species_copy_range.
 *
 * This is runtime-dispatched and a no-op unless the damping type enables
 * filtered-state evolution (currently low-pass filter).
 */
void gk_species_damping_copy_range(struct gk_species *gks, struct gkyl_array *fout,
  const struct gkyl_array *fin, const struct gkyl_range *range);

/**
 * Release species damping object.
 *
 * @param app gyrokinetic app object.
 * @param damp Species damping object to release.
 */
void gk_species_damping_release(const struct gkyl_gyrokinetic_app *app, const struct gk_damping *damp);

/**
 * Reset species damping object.
 *
 * @param app gyrokinetic app object.
 * @param tm Time.
 * @param gks Species object.
 * @param damp Species damping object.
 * @param damp_inp New damping input.
 */
void gk_species_damping_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks,
  struct gk_damping *damp, struct gkyl_gyrokinetic_damping damp_inp);

/** gk_species_fdot_multiplier API */

/**
 * Initialize species df/dt multiplier object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param fdmul Species df/dt multiplier object.
 */
void gk_species_fdot_multiplier_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_fdot_multiplier *fdmul);

/**
 * Multiply the CFL rate.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param fdmul Species df/dt multiplier object.
 * @param phi Current electrostatic potential.
 * @param fin Current distribution function.
 * @param f_buffer Phase-space buffer.
 * @param out CFL rate to multiply.
 */
void gk_species_fdot_multiplier_advance_times_cfl(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out);

/**
 * Multiply df/dt.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param fdmul Species df/dt multiplier object.
 * @param phi Current electrostatic potential.
 * @param fin Current distribution function.
 * @param f_buffer Phase-space buffer.
 * @param out df/dt to multiply.
 */
void gk_species_fdot_multiplier_advance_times_rate(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out);

/**
 * Write damping diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for damping diagnostic.
 * @param frame Output frame.
 */
void gk_species_fdot_multiplier_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Release species damping object.
 *
 * @param app gyrokinetic app object.
 * @param fdmul Species df/dt multiplier object.
 */
void gk_species_fdot_multiplier_release(const struct gkyl_gyrokinetic_app *app, const struct gk_fdot_multiplier *fdmul);

/**
 * Reset the df/dt multiplier operator.
 *
 * @param app gyrokinetic app object.
 * @param tm Current simulation time.
 * @param gks Pointer to species.
 * @param fdmul Species df/dt multiplier object.
 * @param fdot_mult_inp New input struct for the fdot_multiplier.
 */
void gk_species_fdot_multiplier_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp);

/** gk_anomalous_diff API */

/**
 * Initialize species anomalous diffusion object.
 *
 * @param app Gyrokinetic app object.
 * @param s Species object.
 * @param gkad Species anomalous diffusion object.
 */
void gk_species_anomalous_diff_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_anomalous_diff *gkad);

/**
 * Compute RHS contribution from anomalous diffusion.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param gkad Species anomalous diffusion object.
 * @param fin Input distribution function.
 * @param rhs Anomalous diffusion contribution to df/dt.
 */
void gk_species_anomalous_diff_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write out diagnostics from the anomalous diffusion term.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param gkad Species anomalous diffusion object.
 * @param tm Current simulation time.
 * @param frame Current I/O frame.
 */
void gk_species_anomalous_diff_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame);

/**
 * Release species anomalous diffusion object.
 *
 * @param app Gyrokinetic app object.
 * @param gkad Species anomalous diffusion object.
 */
void gk_species_anomalous_diff_release(const struct gkyl_gyrokinetic_app *app, const struct gk_anomalous_diff *gkad);

/** gk_source_bgk API */

/**
 * Initialize species source_bgk object.
 *
 * @param app Gyrokinetic app object.
 * @param s Species object.
 * @param src Species source_bgk object.
 */
void gk_species_source_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source_bgk *src);

/**
 * Compute RHS contribution from source_bgk source.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param src Pointer to source.
 * @param fin Input distribution function.
 * @param rhs source_bgk contribution to df/dt.
 */
void gk_species_source_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write out diagnostics from the source_bgk source.
 *
 * @param app Gyrokinetic app object.
 * @param gks Pointer to species.
 * @param src Pointer to source.
 * @param tm Current simulation time.
 * @param frame Current I/O frame.
 */
void gk_species_source_bgk_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm, int frame);

/**
 * Calculate integrated diagnostics of the BGK source operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param src BGK source object.
 * @param tm Current simulation time.
 */
void gk_species_source_bgk_calc_integrated_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm);

/**
 * Write integrated diagnostics of the BGK source operator.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 * @param src BGK source object.
 */
void gk_species_source_bgk_write_integrated_diags(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_source_bgk *src);

/**
 * Release species source_bgk object.
 *
 * @param app Gyrokinetic app object.
 * @param src Species source_bgk object to release.
 */
void gk_species_source_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source_bgk *src);

/** gk_species API */

/**
 * Initialize species.
 *
 * @param gk Input gk data.
 * @param app gyrokinetic app object.
 * @param s On output, initialized species object.
 */
void gk_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param t0 Time for use in ICs.
 */
void gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute the part of the species initial conditions that depends on other
 * species.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param t0 Time for use in ICs.
 */
void gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute RHS from species distribution function
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param fin Input distribution function.
 * @param fbar_in Input filtered distribution for the current RK stage.
 * @param rhs On output, the RHS from the species object.
 * @param bflux_moms Output boundary flux moments.
 * @return Maximum stable time-step.
 */
double gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms);

/**
 * Compute the *implicit* RHS from species distribution function
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from the species object.
 * @param bflux_moms Output boundary flux moments.
 * @param dt timestep size (used in the implcit coef.).
 * @return Maximum stable time-step.
 */
double gk_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt);

/**
 * Scale and accumulate for forward euler method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param dt Timestep.
 * @param inp Input array.
 */
void gk_species_step_f(struct gk_species *species, struct gkyl_array* out, double dt,
  const struct gkyl_array* inp);

/**
 * Combine for rk3 method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param c1 Scaling factor.
 * @param arr1 Input array.
 * @param c2 Scaling factor.
 * @param arr2 Input array.
 * @param rng Range.
 */
void gk_species_combine(struct gk_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng);

/**
 * Copy for rk3 method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param inp Input array.
 * @param range Range.
 */
void gk_species_copy_range(struct gk_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Apply the positivity shift (to enforce f>=0) to a charged species.
 *
 * @param app Gyrokinetic app object.
 * @param gks Gyrokinetic species object.
 */
void gk_species_apply_pos_shift(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Apply BCs to dynamic species distribution function.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param f Field to apply BCs.
 */
void gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers.
 */
void gk_species_coll_tm(gkyl_gyrokinetic_app *app);

/**
 * Fill stat object in app with total number of iterations
 * used to correct moments in LTE projection object.
 * Also fills stat object with number of times correction object called. 
 *
 * @param app App object to update stat timers.
 * @param s Species object.
 * @param sidx Index of current species.
 */
void gk_species_n_iter_corr(gkyl_gyrokinetic_app *app, const struct gk_species *gks, int sidx);

/**
 * Species write function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_species_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Species moment write function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Species calc integrated moment function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 */
void gk_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Species write integrated moment function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 */
void gk_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Species calc L2norm function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 */
void gk_species_calc_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Species write L2norm function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 */
void gk_species_write_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Calculate the integrated moments divided by dt, for particle/energy balance
 * diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param dt Time step.
 * @param fdot_int_mom Integrated moment divided by dt (not yet reduced over comm).
 */
void
gk_species_calc_int_mom_dt(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom);

/**
 * Delete resources used in species.
 *
 * @param app gyrokinetic app object.
 * @param species Species object to delete.
 */
void gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *s);

/** gk_neut_species_moment API */

/**
 * Initialize neutral species moment object.
 *
 * @param app gyrokinetic app object.
 * @param s Neutral species object.
 * @param sm Neutral species moment object.
 * @param nm Name string indicating moment type.
 * @param is_integrated Whether to compute the volume integrated moment.
 */
void gk_neut_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated);

/**
 * Calculate neutral species moment, given input neutral distribution function fin.
 * 
 * @param sm Neutral species moment object.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param fin Input neutral distribution function array.
 */
void gk_neut_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Divide the moment by the configuration-space Jacobian for diagnostic
 * purposes. Not all components of the diagnostic moment are divided, depending
 * on the type.
 * 
 * @param app Gyrokinetic app object.
 * @param sm Species moment object.
 * @param Jmom_in Moment(s) to be divided by J.
 * @param mom_out Array in which to place the output.
 */
void gk_neut_species_moment_diag_jacobgeo_div(const struct gkyl_gyrokinetic_app *app,
  struct gk_species_moment *sm, struct gkyl_array *Jmom_in, struct gkyl_array *mom_out);

/**
 * Release neutral species moment object.
 *
 * @param app gyrokinetic app object.
 * @param sm Neutral species moment object to release.
 */
void gk_neut_species_moment_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species_moment *sm);

/** gk_neut_species_lte API */

/**
 * Initialize species lte object.
 *
 * @param app Gyrokinetic app object.
 * @param s Neutral species object.
 * @param lte Neutral species lte object.
 * @param corr_inp Input struct with moment correction inputs.
 */
void gk_neut_species_lte_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_lte *lte, struct correct_all_moms_inp corr_inp);

/**
 * Compute LTE distribution from input moments.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to neutral species.
 * @param lte Pointer to lte object.
 * @param moms_lte Input LTE moments.
 */
void gk_neut_species_lte_from_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_lte *lte, const struct gkyl_array *moms_lte);

/**
 * Compute equivalent LTE distribution from input distribution function. 
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to neutral species.
 * @param lte Pointer to lte.
 * @param fin Input distribution function.
 */
void gk_neut_species_lte(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_lte *lte, const struct gkyl_array *fin);

/**
 * Write the LTE correction status for the neutral species. 
 *
 * @param app Gyrokinetic app object.
 * @param gk_ns Pointer to neutral species.
 */
void gk_neut_species_lte_write_max_corr_status(gkyl_gyrokinetic_app* app, struct gk_neut_species *gk_ns);

/**
 * Release species lte object.
 *
 * @param app gyrokinetic app object.
 * @param lte Neutral species lte object to release.
 */
void gk_neut_species_lte_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte);

/** gk_neut_species_bgk API */

/**
 * Initialize neut species BGK collisions object.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param bgk Species BGK object.
 */
void gk_neut_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk);

/**
 * Initialize species BGK cross-collisions object.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param bgk Species BGK object.
 */
void gk_neut_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk);

/**
 * Compute necessary moments and boundary
 * corrections for BGK collisions.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 */
void gk_neut_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

void gk_neut_species_bgk_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species BGK collisions.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 */
void gk_neut_species_bgk_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

void gk_neut_species_bgk_cross_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin);

/**
 * Compute RHS from BGK collisions (explicit integrator).
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from BGK.
 */
void gk_neut_species_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Compute RHS from BGK collisions (implicit integrator).
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 * @param dt Time step.
 * @param rhs On output, the RHS from BGK.
 */
void gk_neut_species_bgk_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);

/**
 * Write moments from BGK object.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void gk_neut_species_bgk_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Release species BGK object.
 *
 * @param app gyrokinetic app object.
 * @param bgk Species BGK object to release.
 */
void gk_neut_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk);

/** gk_neut_species_recycle API **/

/**
 * Initialize recycling object.
 *
 * @param app Gyrokinetic app object.
 * @param recyc Recycling bc object.
 * @param dir Direction for BC (x, y, or z).
 * @param edge Edge for BC (lower/upper).
 * @param params Input params for recycling BCs.
 * @param s Gk_neut_species to apply BCs for.
 * @param use_gpu Boolean for using GPUs.
 */
void gk_neut_species_recycle_init(struct gkyl_gyrokinetic_app *app, struct gk_recycle_wall *recyc,
  int dir, enum gkyl_edge_loc edge, struct gkyl_gyrokinetic_emission_inp *params,
  struct gk_neut_species *s, bool use_gpu);

/**
 * Initialize recycling cross moments.
 *
 * @param app Gyrokinetic app object.
 * @param s Gk_neut_species to apply BCs for.
 * @param recyc Recycling bc object.
 */
void gk_neut_species_recycle_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc);

/**
 * Apply recycling BCs.
 *
 * @param app Gyrokinetic app object.
 * @param recyc Recycling bc object.
 * @param s Gk_neut_species to apply BCs for.
 * @param fout Gk_neut_species distf.
 */
void gk_neut_species_recycle_apply_bc(struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc,
  const struct gk_neut_species *s, struct gkyl_array *fout);

/**
 * Write recycle flux diagnostics.
 *
 * @param app Gyrokinetic app object.
 * @param s Gk_neut_species to apply BCs for.
 * @param recyc Recycling bc object.
 * @param tm Simulation time.
 * @param frame Simulation frame.
 */
void
gk_neut_species_recycle_write_flux(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc, double tm, int frame);

/**
 * Release recycle BC object.
 *
 * @param app Gyrokinetic app object
 * @param recyc Recycling bc object
 */
void gk_neut_species_recycle_release(const struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc);

/** gk_neut_species_react API */

/**
 * Initialize neutral species reactions object.
 *
 * @param app gyrokinetic app object.
 * @param s Neutral species object.
 * @param inp Input reaction struct for determining types of reactions.
 * @param react Neutral species reaction object.
 */
void gk_neut_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gkyl_gyrokinetic_react inp, struct gk_react *react);

/**
 * Initialize neutral species reactions "cross-collisions" object
 * for who is reacting with whom.
 *
 * @param app gyrokinetic app object.
 * @param s Neutral species object.
 * @param react Neutral species react object.
 */
void gk_neut_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_react *react);

/**
 * Compute necessary rates and moments for reactions.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to neutral species.
 * @param react Pointer to react.
 * @param fin Input distribution functions (size: num_species).
 * @param fin_neut Input neutral distribution functions (size: num_neut_species).
 */
void gk_neut_species_react_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species,
  struct gk_react *react,
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute RHS from reactions for neutrals
 * (e.g., ionization, charge exchange, recombination, or radiation).
 *
 * @param app gyrokinetic app object.
 * @param s Pointer to neutral species.
 * @param react Pointer to react.
 * @param fin Input neutral distribution function.
 * @param rhs On output, the neutral RHS from react (df/dt).
 */
void gk_neut_species_react_rhs(gkyl_gyrokinetic_app *app,
  struct gk_neut_species *s, struct gk_react *react,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write neutral reaction rate.
 *
 * @param app gyrokinetic app object.
 * @param gkns Pointer to neutral species.
 * @param gkr Pointer to react object.
 * @param ridx Index for reaction species.
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void gk_neut_species_react_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
  int ridx, double tm, int frame);

/**
 * Release neutral species react object.
 *
 * @param app gyrokinetic app object.
 * @param react Neutral species react object to release.
 */
void gk_neut_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react);

/** gk_species_scaling API */

/**
 * Initialize operator that scales the species according to some criteria.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 */
void gk_species_scaling_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_scaling *sca);

/**
 * Initialize the part of scaling that depends on other species.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 */
void  gk_species_scaling_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca);

/**
 * Store initial condition of the species.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 */
void gk_species_scaling_apply_ic_cross(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca);

/**
 * Compute the cross-species moments needed.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 * @param fin Input distribution for charged species.
 * @param fin_neut Input distribution/moments for neutral species.
 */
void gk_species_scaling_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute factor to scale species by.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 * @param fin Input distribution/moments for species.
 * @param rhs df/dt to add to.
 */
void gk_species_scaling_rhs(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Add up the integrated boundary fluxes and scale the species.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 * @param fin Distribution/moments to scale.
 * @param bflux Boundary fluxes of charged species.
 */
void gk_species_scaling_apply(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[]);

/**
 * Write scaling diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param sca Recycle react scale object.
 * @param fin Input distribution/moments for species.
 * @param rhs df/dt to add to.
 */
void gk_species_scaling_write(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame);

/**
 * Free memory associated with the scaling operator.
 *
 * @param app gyrokinetic app object.
 * @param sca Recycle react scale object.
 */
void gk_species_scaling_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_scaling *sca);

/** gk_neut_species_scaling API */

/**
 * Initialize operator that scales the neutral species according to some criteria.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 */
void gk_neut_species_scaling_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns, 
  struct gk_scaling *sca);

/**
 * Initialize the part of scaling that depends on other species.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 */
void  gk_neut_species_scaling_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca);

/**
 * Store initial condition of the neutral species.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 */
void gk_neut_species_scaling_apply_ic_cross(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca);

/**
 * Compute the cross-species moments needed.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 * @param fin Input distribution for charged species.
 * @param fin_neut Input distribution/moments for neutral species.
 */
void gk_neut_species_scaling_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *ns,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute factor to scale neutrals by.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 * @param fin Input distribution/moments for neutral species.
 * @param rhs df/dt to add to.
 */
void gk_neut_species_scaling_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Add up the integrated boundary fluxes and scale the neutral species.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 * @param fin Distribution/moments to scale.
 * @param bflux Boundary fluxes of charged species.
 */
void gk_neut_species_scaling_apply(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[]);

/**
 * Write scaling diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param ns Neutral species object.
 * @param sca Recycle react scale object.
 * @param fin Input distribution/moments for neutral species.
 * @param rhs df/dt to add to.
 */
void gk_neut_species_scaling_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *ns,
  struct gk_scaling *sca, int ridx, double tm, int frame);

/**
 * Free memory associated with the scaling operator.
 *
 * @param app gyrokinetic app object.
 * @param sca Recycle react scale object.
 */
void gk_neut_species_scaling_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_scaling *sca);

/** gk_neut_species_projection API */

/**
 * Initialize neutral species projection object.
 *
 * @param app gyrokinetic app object.
 * @param s Neutral species object.
 * @param inp Input struct for projection (contains functions pointers for type of projection).
 * @param proj Neutral species projection object.
 */
void gk_neut_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj);

/**
 * Compute neutral species projection
 *
 * @param app gyrokinetic app object.
 * @param species Neutral species object.
 * @param proj Neutral species projection object.
 * @param f Output Neutral distribution function from projection.
 * @param tm Time for use in projection.
 */
void gk_neut_species_projection_calc(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, 
  struct gk_proj *proj, struct gkyl_array *f, double tm);

/**
 * Release neutral species projection object.
 *
 * @param app gyrokinetic app object.
 * @param proj Neutral species projection object to release.
 */
void gk_neut_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj);

/** gk_neut_species_source API */

/**
 * Initialize neutral species source object.
 *
 * @param app gyrokinetic app object.
 * @param s Neutral species object.
 * @param src Neutral species source object.
 */
void gk_neut_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_source *src);

/**
 * Compute Neutral species applied source term
 *
 * @param app gyrokinetic app object.
 * @param species Neutral species object.
 * @param src Neutral species source object.
 * @param f_buffer Phase-space buffer used to project the source.
 * @param tm Time for use in source.
 */
void gk_neut_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, 
  struct gk_source *src, struct gkyl_array *f_buffer, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to Neutral species.
 * @param src Pointer to source.
 * @param fin Input neutral distribution function.
 * @param rhs On output, the incremented rhs (df/dt).
 */
void gk_neut_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write neutral source diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gkns Pointer to species.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void gk_neut_species_source_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Write neutral source moment diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gkns Pointer to species.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void gk_neut_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Calc neutral source integrated moment diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gkns Pointer to species.
 * @param tm Time for source diagnostic.
 */
void gk_neut_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);

/**
 * Write neutral source integrated moment diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gkns Pointer to species.
 * @param tm Time for source diagnostic.
 */
void gk_neut_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

/**
 * Release Neutral species source object.
 *
 * @param app gyrokinetic app object.
 * @param src Neutral species source object to release.
 */
void gk_neut_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_neut_species API */

/**
 * Initialize neutral species.
 *
 * @param gk Input gk data.
 * @param app gyrokinetic app object.
 * @param s On output, initialized neutral species object.
 */
void gk_neut_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s);

/**
 * Compute neutral species initial conditions.
 *
 * @param app gyrokinetic app object.
 * @param species Neutral species object.
 * @param t0 Time for use in ICs.
 */
void gk_neut_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0);

/**
 * Compute the part of the neutral species initial conditions
 * that depends on other species.
 *
 * @param app gyrokinetic app object.
 * @param species Neutral species object.
 * @param t0 Time for use in ICs.
 */
void gk_neut_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0);

/**
 * Compute RHS from neutral species distribution function.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to neutral species.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from the neutral species object (df/dt).
 * @param bflux_moms Output boundary flux moments.
 * @return Maximum stable time-step.
 */
double gk_neut_species_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms);

/**
 * Compute the *implicit* RHS from neutral species distribution function.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to neutral species.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from the species object.
 * @param bflux_moms Output boundary flux moments.
 * @param dt timestep size (used in the implcit coef).
 * @return Maximum stable time-step.
 */
double gk_neut_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt);

/**
 * Apply BCs to neutral species distribution function.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to neutral species.
 * @param f Field to apply BCs.
 */
void gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with total number of iterations
 * used to correct moments in LTE projection object for neutral species.
 * Also fills stat object with number of times correction object called. 
 *
 * @param app App object to update stat timers.
 * @param s Species object.
 * @param sidx Index of current species.
 */
void gk_neut_species_n_iter_corr(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx);

/**
 * Scale and accumulate for forward euler method.
 *
 * @param species Pointer to neutral species.
 * @param out Output array.
 * @param dt Timestep.
 * @param inp Input array.
 */
void gk_neut_species_step_f(struct gk_neut_species *species, struct gkyl_array* out, double dt,
  const struct gkyl_array* inp);

/**
 * Combine for rk3 method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param c1 Scaling factor.
 * @param arr1 Input array.
 * @param c2 Scaling factor.
 * @param arr2 Input array.
 * @param rng Range.
 */
void gk_neut_species_combine(struct gk_neut_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng);

/**
 * Copy for rk3 method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param inp Input array.
 * @param range Range.
 */
void gk_neut_species_copy_range(struct gk_neut_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Apply the positivity shift (to enforce f>=0) to a neutral species.
 *
 * @param app Gyrokinetic app object.
 * @param gkns Neutral species object.
 */
void gk_neut_species_apply_pos_shift(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

/**
 * Species write function.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_neut_species_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Species moment write function.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_neut_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Species calc integrated moment function.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 * @param tm simulation time.
 */
void gk_neut_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);

/**
 * Species write integrated moment function.
 *
 * @param app gyrokinetic app object.
 * @param gkns Neutral species object.
 */
void gk_neut_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

/**
 * Delete resources used in neutral species.
 *
 * @param app gyrokinetic app object.
 * @param species Neutral species object to delete.
 */
void gk_neut_species_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s);

/** gk_field API */

/**
 * Create new field object.
 *
 * @param gk Input gk data.
 * @param app gyrokinetic app object.
 * @return Newly created field.
 */
struct gk_field* gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app);

/**
 * Compute biased wall potentials.
 *
 * @param app gyrokinetic app object.
 * @param field Pointer to field.
 * @param tm Time to compute biased wall potentials at.
 */
void gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm);

/**
 * Accumulate charge density for Poisson solve.
 *
 * @param app gyrokinetic app object.
 * @param field Pointer to field.
 * @param fin[] Input distribution function (num_species size).
 * @param bflux Boundary fluxes of charged species.
 */
void gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[], struct gkyl_array **bflux[]);

/**
 * Compute EM field.
 *
 * @param app gyrokinetic app object.
 * @param field Pointer to field.
 * @param em Output field.
 */
void gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Read the field from a file.
 *
 * @param app Gyrokinetic app.
 * @param inp Input struct with importing parameters (and file name).
 */
void gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp);

/**
 * Project the initial field using a user provided function.
 *
 * @param app Gyrokinetic app.
 */
void gk_field_project_init(struct gkyl_gyrokinetic_app *app);

/**
 * Compute field energy diagnostic.
 *
 * @param app gyrokinetic app object.
 * @param tm Time at which diagnostic is computed.
 * @param field Pointer to field.
 */
void gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field);

/**
 * Compute field energy divided by dt for energy balance diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param field Pointer to field.
 * @param dt Time step.
 * @param energy_reduced Integrated field energy (single element double array).
 */
void gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced);

/**
 * Release resources allocated by field.
 *
 * @param app gyrokinetic app object.
 * @param f Field object to release.
 */
void gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f);

/** gk_eirene API */

/**
 * Initialize Eirene coupling object.
 *
 * @param app Gyrokinetic app object.
 * @param gk Input gk data.
 * @returns Eirene object
 */
struct gk_eirene* gk_eirene_init(gkyl_gyrokinetic_app *app, struct gkyl_gk *gk);

/**
 * Compute RHS for sources from eirene
 *
 * @param app Gyrokinetic app object.
 * @param fin Array of distribution functions (one for each species) .
 * @param fout Output array of distribution functions (one for each species).
 */
void gk_eirene_rhs(gkyl_gyrokinetic_app *app, const struct gkyl_array *fin[], struct gkyl_array *rhs[]);

/**
 * Write moments of eirene sources
 *
 * @param app Gyrokinetic app object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_eirene_write(gkyl_gyrokinetic_app *app, double tm, int frame);

/**
* Compute integrated moments of eirene sources
 *
 * @param app Gyrokinetic app object.
 * @param tm simulation time.
 */
void gk_eirene_calc_integrated_diagnostics(struct gkyl_gyrokinetic_app *app, double tm);

/**
* Write integrated moments of eirene sources
 *
 * @param app Gyrokinetic app object.
 */
void gk_eirene_write_integrated_diagnostics(struct gkyl_gyrokinetic_app *app);

/**
* Release gk_eirene object
 *
 * @param app Gyrokinetic app object.
 * @param eirene Eirene object.
 */
void gk_eirene_release(gkyl_gyrokinetic_app *app, struct gk_eirene *eirene);

/** Time stepping API */

/**
 * Compute the gyrokinetic fields.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param fin Array of distribution functions (one for each species) .
 * @param bflux Boundary fluxes of charged species.
 */
void gyrokinetic_calc_field(gkyl_gyrokinetic_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[]);

/**
 * Compute the gyrokinetic fields and apply boundary conditions.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param distf Array of distribution functions (for each charged species).
 * @param bflux Boundary fluxes of charged species.
 * @param distf_neut Array of distribution functions (for each neutral species).
 */
void gyrokinetic_calc_field_and_apply_bc(gkyl_gyrokinetic_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array **bflux[], struct gkyl_array *distf_neut[]);

/**
 * Compute the RHS of the gyrokinetic equation (df/dt) and the minimum time
 * step it requires for stability based on the CFL constraint.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param dt Suggested time step.
 * @param fin Input array of charged-species distribution functions.
 * @param fbar_in Input array of charged-species filtered distributions.
 * @param fout Output array of charged-species distribution functions.
 * @param bflux_out Output array of charged-species boundary fluxes.
 * @param fin_neut Input array of neutral-species distribution functions.
 * @param fout_neut Output array of neutral-species distribution functions.
 * @param bflux_out_neut Output array of neutral-species boundary fluxes.
 * @param st Time stepping status object.
 */
void gyrokinetic_rhs(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fbar_in[],
  struct gkyl_array *fout[], struct gkyl_array **bflux_out[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], struct gkyl_array **bflux_out_neut[],
  struct gkyl_update_status *st); 

/**
 * Compute the RHS of the gyrokinetic equation (df/dt) due to implicit terms.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param dt Suggested time step.
 * @param fin Input array of charged-species distribution functions.
 * @param fout Output array of charged-species distribution functions.
 * @param bflux_out Output array of charged-species boundary fluxes.
 * @param fin_neut Input array of neutral-species distribution functions.
 * @param fout_neut Output array of neutral-species distribution functions.
 * @param bflux_out_neut Output array of neutral-species boundary fluxes.
 * @param st Time stepping status object.
 */
void gyrokinetic_rhs_implicit(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_array *fin[], struct gkyl_array *fout[], struct gkyl_array **bflux_out[], 
  struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], struct gkyl_array **bflux_out_neut[],
  struct gkyl_update_status *st); 

/**
 * Take time-step using the RK3 method. Also sets the status object
 * which has the actual and suggested dts used. These can be different
 * from the actual time-step.
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
struct gkyl_update_status gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app* app, double dt0);

/**
 * Take time-step of the (BGK) collision operator using a first order implicit method. 
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
void gyrokinetic_update_implicit_coll(gkyl_gyrokinetic_app *app,  double dt0);

/**
 * Take time-step using a first order operator split combining 
 * the RK3 method for the collisionless advection with a first order implicit
 * method for BGK collisions. The first order implicit step occurs *after* the 
 * RK3 step and utilizes the stable RK3 time step. Also sets the status object
 * which has the actual and suggested dts used. These can be different
 * from the actual time-step.
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
struct gkyl_update_status gyrokinetic_update_op_split(gkyl_gyrokinetic_app *app,  double dt0);

/**
 * Initialize enforcement of quasineutrality of the guiding centers after
 * applying the positivity operator to enforce f>=0 of each charged species.
 *
 * @param app Gyrokinetic app object.
 */
void gyrokinetic_post_positivity_quasineut_init(gkyl_gyrokinetic_app* app);

/**
 * Enforce quasineutrality of the guiding centers after applying the positivity
 * operator to enforce f>=0 of each charged species.
 *
 * @param app Gyrokinetic app object.
 * @param fout State field of each species.
 */
void gyrokinetic_post_positivity_quasineut(gkyl_gyrokinetic_app *app, struct gkyl_array *fout[]);

/**
 * Release resources used to enforce quasineutrality of the guiding centers after
 * applying the positivity operator to enforce f>=0 of each charged species.
 *
 * @param app Gyrokinetic app object.
 */
void gyrokinetic_post_positivity_quasineut_release(gkyl_gyrokinetic_app* app);
