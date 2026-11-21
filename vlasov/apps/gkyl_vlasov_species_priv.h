#pragma once

// Private header for the Vlasov species objects: the species container (struct
// vlasov_species), its kinetic (struct vm_species) and fluid (struct
// vm_fluid_species) aspects, their sub-structs, and the species private API.
// struct vm_geom is defined in gkyl_vlasov_priv.h.

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_basic.h>
#include <gkyl_bc_emission.h>
#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_prim_vars.h>
#include <gkyl_dg_calc_fluid_vars.h>
#include <gkyl_dg_calc_fluid_em_coupling.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_canonical_pb_fluid.h>
#include <gkyl_dg_gr_maxwell_conf_flux_surf.h>
#include <gkyl_dg_gr_maxwell_current_deposition.h>
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_lorentz_conf.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_gaussian_filter.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_updater_fluid.h>
#include <gkyl_dg_updater_diffusion_fluid.h>
#include <gkyl_dg_updater_diffusion_gen.h>
#include <gkyl_dg_updater_lbo_vlasov.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov_calc_radiation.h>
#include <gkyl_dg_vlasov_conf_flux_surf.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_util.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_cross_prim_moms_bgk.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>
#include <gkyl_bc_pkpm_species_reflect.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_pkpm_mom_diffusion.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_mom_pkpm.h>
#include <gkyl_pkpm_bgk_collisions.h>
#include <gkyl_pkpm_coll_nu.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_dg_updater_lbo_pkpm_nodal.h>
#include <gkyl_pkpm_coll_cross.h>
#include <gkyl_pkpm_stage_weights.h>
#include <gkyl_pkpm_lbo_moms.h>
#include <gkyl_pkpm_lte.h>
#include <gkyl_pkpm_temp_source.h>
#include <gkyl_wv_canonical_pb_fluid.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_maxwell.h>

struct vm_geom; // geometry data, defined in gkyl_vlasov_priv.h, owned by the app.

// ---- species/fluid building-block sub-structs ----
// struct for holding moment correction inputs
struct correct_all_moms_inp {
  bool correct_all_moms; // boolean if we are correcting all the moments or only density
  double iter_eps; // error tolerance for moment fixes (density is always exact)
  int max_iter; // maximum number of iterations
  bool use_last_converged; // use last iteration value regardless of convergence?
};

// data for moments
struct vm_species_moment {
  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
  // Options for moment calculation:
  // 1. Compute the moment directly with gkyl_mom_calc
  // 2. Compute the moments of the equivalent LTE (local thermodynamic equilibrium)
  //    distribution (n, V_drift, T/m) with specialized updater
  //    Note: in relativity V_drift is the bulk four-velocity (GammaV, GammaV*V_drift)
  union {
    struct {
      struct gkyl_vlasov_lte_moments *vlasov_lte_moms; // Updater for computing LTE moments.
    };
    struct {
      struct gkyl_mom_type *mom_type; // Moment type.
      struct gkyl_mom_calc *mom_calc; // Moment calculator.
    };
  };

  // Moment specific pointers for the moment hamiltonian and it's range
  struct gkyl_array *mom_hamil;
  struct gkyl_range *mom_hamil_range;

  bool is_vlasov_lte_moms;
  bool is_integrated; // =True means volume integrated moment.
};

// forward declare species struct
struct vm_species;

// Context for the computational->physical coordinate transform (c2p) used when
// projecting on non-uniform meshes: configuration coords are mapped by the
// position map, velocity coords (for phase-space projections) by the velocity map.
struct vm_proj_c2p_ctx {
  int cdim; // number of configuration-space dimensions
  const struct gkyl_vlasov_position_map *pos_map; // configuration-space map
  const struct gkyl_vlasov_velocity_map *vel_map; // velocity-space map
};

// Phase-space computational-to-physical map (position map for the conf coordinates, velocity
// map for the velocity coordinates), shared by the projections and the emission BCs.
void vm_proj_c2p_phase(const double *xcomp, double *xphys, void *ctx);

struct vm_proj {
  enum gkyl_projection_id proj_id; // type of projection
  enum gkyl_model_id model_id;
  struct vm_proj_c2p_ctx
    c2p_ctx; // coordinate-map context for projections (captured by proj_on_basis)
  // organization of the different projection objects and the required data and solvers
  union {
    // function projection
    struct {
      struct gkyl_proj_on_basis *proj_func; // projection operator for specified function
      struct gkyl_array *proj_host; // array for projection on host-side if running on GPUs
    };
    // LTE (Local thermodynamic equilibrium) distribution function project with moment correction
    // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
    struct {
      struct gkyl_array *dens; // host-side density
      struct gkyl_array *V_drift; // host-side V_drift
      struct gkyl_array *T_over_m; // host-side T/m (temperature/mass)

      struct gkyl_array *vlasov_lte_moms_host; // host-side LTE moms (n, V_drift, T/m)
      struct gkyl_array *vlasov_lte_moms; // LTE moms (n, V_drift, T/m) for passing to updaters

      struct gkyl_proj_on_basis *proj_dens; // projection operator for density
      struct gkyl_proj_on_basis *proj_V_drift; // projection operator for V_drift
      struct gkyl_proj_on_basis *proj_temp; // projection operator for temperature

      // LTE distribution function projection object
      // also corrects the density of projected distribution function
      struct gkyl_vlasov_lte_proj_on_basis *proj_lte;

      // Correction updater for ensuring LTE distribution has desired LTE (n, V_drift, T/m) moments
      bool correct_all_moms; // boolean if we are correcting all the moments
      struct gkyl_vlasov_lte_correct *corr_lte;
    };
  };
};

struct vm_lte {
  struct gkyl_array *f_lte;

  // LTE moments (n, V_drift, T/m) of the stored J_x J_v f, with the configuration-space
  // Jacobian divided out of the density (see vm_species_lte_moms), so they are physical.
  struct vm_species_moment moms;

  // LTE distribution function projection object
  // also corrects the density of projected distribution function
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte;
  // LTE correction object.
  struct gkyl_vlasov_lte_correct *corr_lte;

  long n_iter; // total number of iterations from correcting moments
  long num_corr; // total number of times the correction updater is called
  bool correct_all_moms; // boolean if we are correcting all the moments
  gkyl_dynvec corr_stat;
  bool is_first_corr_status_write_call;
};

struct vm_collisionless {
  bool has_app_accel; // flag to indicate there is applied acceleration
  bool app_accel_evolve; // flag to indicate applied acceleration is time dependent
  struct gkyl_array *app_accel; // applied acceleration
  struct gkyl_array *app_accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_accel_proj; // projector for acceleration

  bool use_lo; // bool to determine if using low-order kernels for non-canonical Hamiltonian models.
  bool use_preset_geom; // bool to determine if we are using triad input geom
  bool use_vierbein; // bool to determine if using vierbein inputs for triads.
  bool
    use_extended_hamil_def; // bool to determine if we are using the extended hamil definitions which includes potentials

  double qbym; // Charge (q) divided by mass (m).
  struct gkyl_array *qmem; // array for q/m*(E,B)
  struct gkyl_array *
    em_no_J; // physical E,B (stored field J*E,J*B divided by the conf Jacobian) for the Lorentz force on a mapped grid
  struct gkyl_array *pot_tot; // array for total potentials (q/m*(phi + phi_ext) + phi_g, q/m*A_ext)
  bool has_E; // Do we have electric fields?
  bool has_phi; // Do we have scalar potentials (electrostatic/gravitational)?
  bool has_B; // Do we have magnetic fields?
  bool has_gr_em_triad_coupling; // Do we need GR-Maxwell Lorentz-force fields for triad species?
  int num_surf_conf_basis; // Surface modal basis size per configuration direction.
  int num_surf_vel_basis; // Surface modal basis size per velocity direction.

  struct gkyl_array *conf_flux_surf; // Modal expansion of surface fluxes at conf-space surfaces.
  struct gkyl_array *vel_flux_surf; // Modal expansion of surface fluxes at velocity-space surfaces.
  struct gkyl_dg_vlasov_conf_flux_surf
    *calc_conf_flux; // Updater for computing modal expansion of surface fluxes (conf).
  struct gkyl_dg_vlasov_vel_flux_surf
    *calc_vel_flux; // Updater for computing modal expansion of surface fluxes (vel).
  struct gkyl_dg_gr_maxwell_lorentz_conf *calc_lorentz; // Updater for local GR Lorentz-force fields.
  struct gkyl_dg_gr_maxwell_current_deposition
    *calc_current_dep; // Updater for GR current deposition.

  struct gkyl_dg_eqn *eqn; // Vlasov equation object.
  struct gkyl_hyper_dg *slvr; // Vlasov solver.

  void (*rhs_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, struct vm_collisionless *cls,
    const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs
  );
};

struct vm_lbo_collisions {
  enum gkyl_collision_id collision_id; // type of collisions
  bool write_coll_diagnostics; // Whether to write diagnostics out.

  struct gkyl_array *self_nu; // Self-collision frequency.
  struct gkyl_array *boundary_corrections; // LBO boundary corrections.
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator.
  struct gkyl_array *nu_sum; // Sum of collision frequencies.
  struct gkyl_array *prim_moms, *nu_prim_moms; // Primitive moments.
  struct gkyl_array *nu_sum_host, *nu_prim_moms_host; // Host arrays for I/O.
  bool norm_nu_self; // Whether to compute self-species collision frequency in space and time.
  double norm_nu_fac_self; // Self collision frequency without factor of n_s/(2*v_ts^2)^(3/2).
  double vtsq_min; // Minimum vtsq.
  struct gkyl_spitzer_coll_freq
    *spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value.
  struct gkyl_array *nu_boundary_corrections; // Boundary corrections multiplied by nu.
  struct gkyl_array *nu_moms; // Moments multiplied by nu.
  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator

  int num_cross_collisions; // number of species we cross-collide with
  struct vm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  bool norm_nu_cross; // Whether to compute cross-species collision frequency in space and time.
  bool needs_lte_moms; // Whether the per-stage LTE moments (n, V_drift, T/m) must be computed:
    // set if this species' self or cross nu is computed (Spitzer), or if a
    // partner's computed cross nu reads this species' LTE moments.
  double norm_nu_fac_cross
    [GKYL_MAX_SPECIES]; // Cross collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2).
  double alpha_E_fac[GKYL_MAX_SPECIES]; // Time-independent factor in alpha_E.
  double betaGreenep1; // Value of Greene's factor beta + 1.
  double delta_sr; // Free parameter in relationship between alpha_E and nu_sr.
  double other_m[GKYL_MAX_SPECIES]; // Masses of species colliding with.
  struct gkyl_array
    *other_prim_moms[GKYL_MAX_SPECIES]; // Self-primitive moments of species colliding with.
  struct gkyl_array *cross_prim_moms[GKYL_MAX_SPECIES]; // Cross-primitive moments.
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // Cross-species collision frequencies.
  struct gkyl_array
    *cross_nu_prim_moms; // Weak multiplication of collision frequency and primitive moments.
  struct gkyl_array *alpha_E; // Morse's alpha_E factor.
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  struct gkyl_array *other_moms[GKYL_MAX_SPECIES]; // Moments (M0, M1, M2) of species colliding with.

  // PKPM partners (three velocity dimensions, CPU). A kinetic partner gives the cross operator
  // its own moments and primitive moments; for a PKPM partner this species forms both from the
  // node values of the partner's fluid state. The PKPM species has no collision operator that
  // moves momentum or energy, so its fluid state takes the opposite of what the cross operator
  // gives this species: exch, read by the right-hand side of the partner's fluid state.
  int num_pkpm_partners; // number of PKPM species among the partners
  const struct vlasov_species *pkpm_partner[GKYL_MAX_SPECIES]; // a PKPM partner (NULL: kinetic)
  struct gkyl_array *exch[GKYL_MAX_SPECIES]; // [d(m M1)/dt (3), d(m M2)/dt] with a PKPM partner
  struct gkyl_array *exch_term; // one term of the above
  struct gkyl_array *prim_nodal; // node values [n, u (3), vtsq] of this species, for the partners
  gkyl_pkpm_coll_cross *pkpm_conv; // conversions between node values and moments
  // Moments of partner coll_idx for the cross operator (nothing to do for a kinetic partner).
  void (*other_moms_func[GKYL_MAX_SPECIES])(
    gkyl_vlasov_app *app, struct vm_lbo_collisions *lbo, int coll_idx
  );
  // Momentum and energy the cross operator with partner coll_idx gives this species, from the
  // cross primitive moments times the frequency (nothing to do for a kinetic partner).
  void (*exch_func[GKYL_MAX_SPECIES])(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo, int coll_idx,
    const struct gkyl_array *nu_cross_prim_moms
  );
  // Node values of this species for its PKPM partners (nothing to do without one).
  void (*prim_nodal_func)(gkyl_vlasov_app *app, struct vm_lbo_collisions *lbo);

  struct vm_species_moment moms; // Moments needed in LBO (M0, M1, M2).

  gkyl_dg_updater_collisions *coll_slvr; // collision solver

  // Pointers to methods chosen at runtime.
  void (*moms_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo,
    const struct gkyl_array *fin
  );
  void (*self_nu_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo,
    const struct gkyl_array *fin
  );
  void (*cross_nu_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo, int coll_idx
  );
  void (*cross_moms_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo
  );
  void (*rhs_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo,
    const struct gkyl_array *fin, struct gkyl_array *rhs
  );
  void (*write_mom_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
};

struct vm_bgk_collisions {
  enum gkyl_collision_id collision_id; // Type of collisions.
  bool write_coll_diagnostics; // Whether to write diagnostics out.

  struct gkyl_array *self_nu; // Self-collision frequency.
  struct gkyl_array
    *ref_self_nu; // Reference self-collision frequency needed for robustness checks n, T < 0.
  struct gkyl_array *nu_sum; // Sum of collision frequencies.
  struct gkyl_array *nu_sum_host; // Host arrays for I/O.
  bool norm_nu_self; // Whether to compute self-species collision frequency in space and time.
  double norm_nu_fac_self; // Self collision frequency without factor of n_s/(2*v_ts^2)^(3/2).
  double vtsq_min; // Minimum vtsq.
  struct gkyl_spitzer_coll_freq
    *spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  int num_cross_collisions; // number of species we cross-collide with
  struct vm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  bool norm_nu_cross; // Whether to compute cross-species collision frequency in space and time.
  double norm_nu_fac_cross
    [GKYL_MAX_SPECIES]; // Cross collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2).
  double alpha_E_fac[GKYL_MAX_SPECIES]; // Time-independent factor in alpha_E.
  double betaGreenep1; // Value of Greene's factor beta + 1.
  double delta_sr; // Free parameter in relationship between alpha_E and nu_sr.
  double other_m[GKYL_MAX_SPECIES]; // Masses of species colliding with.
  struct gkyl_array
    *other_prim_moms[GKYL_MAX_SPECIES]; // Self-primitive moments of species colliding with.
  struct gkyl_array *cross_prim_moms; // Cross-primitive moments.
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // Cross-species collision frequencies.
  struct gkyl_array *ref_cross_nu
    [GKYL_MAX_SPECIES]; // Reference cross-collision frequencies needed for robustness checks n, T < 0.
  struct gkyl_array *alpha_E; // Morse's alpha_E factor.
  struct gkyl_vlasov_cross_prim_moms_bgk *cross_calc; // Cross-species moment computation.

  struct gkyl_array *nu_f_lte; // Collision frequency times Maxwellian.
  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep).

  bool fixed_temp_relax; // Boolean for whether the temperature being relaxed to is fixed in time.
  struct gkyl_array *fixed_temp; // Array of fixed temperature BGK collisions are relaxing to.

  bool implicit_step; // Whether or not to take an implicit BGK step.

  // Pointers to methods chosen at runtime.
  void (*moms_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk,
    const struct gkyl_array *fin
  );
  void (*moms_func_implicit)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk,
    const struct gkyl_array *fin
  );
  void (*self_nu_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk
  );
  void (*cross_nu_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk, int coll_idx
  );
  void (*alpha_E_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk, int coll_idx
  );
  void (*cross_moms_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk, int coll_idx
  );
  void (*rhs_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, struct vm_bgk_collisions *bgk,
    const struct gkyl_array *fin, double dt, struct gkyl_array *rhs
  );
  void (*rhs_func_implicit)(
    gkyl_vlasov_app *app, struct vm_species *vms, struct vm_bgk_collisions *bgk,
    const struct gkyl_array *fin, double dt, struct gkyl_array *rhs
  );
  void (*target_func_implicit)(
    gkyl_vlasov_app *app, struct vm_species *vms, struct vm_bgk_collisions *bgk,
    const struct gkyl_array *fin
  );
  void (*fixed_temp_calc_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk,
    const struct gkyl_array *fin
  );
  void (*write_mom_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
};

// BGK-form temperature source nu_Q (f_M(n, u, T_Q) - f) (see struct gkyl_vlasov_source_bgk)
struct vm_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // type of temperature source
  bool is_implicit; // backward Euler in each RK stage instead of the explicit update

  struct gkyl_array *nu; // source rate nu_Q, zeroed where the cell-average density or T_Q <= 0
  struct gkyl_array *ref_nu; // source rate nu_Q as projected
  struct gkyl_array *vtsq; // source temperature over the mass, T_Q/m
  struct gkyl_array *moms; // LTE moments (n, V_drift, T_Q/m) of the target
  struct gkyl_array *nu_f_lte; // source rate times the target Maxwellian
  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep)

  // Implicit source only.
  bool has_implicit_coll; // the species also has implicit BGK collisions (taken in the same step)
  struct gkyl_array *nu_tot; // source rate plus the collision frequencies (with collisions only)
  struct gkyl_array *vtsq_new; // T'/m, the temperature the source gives the species over a step
  struct gkyl_array *div_num, *div_den; // numerator and denominator of T'/m
  gkyl_dg_bin_op_mem *div_mem; // memory of the weak division

  // Methods chosen at runtime. Explicit update: moms_func, rhs_func. Implicit update:
  // moms_implicit_func, implicit_func. The methods of the other update, and all of them without
  // a source or for a static species, are empty.
  void (*moms_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
    const struct gkyl_array *fin
  );
  void (*rhs_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
    const struct gkyl_array *fin, struct gkyl_array *rhs
  );
  void (*moms_implicit_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
    const struct gkyl_array *fin, double dt
  );
  void (*implicit_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
    const struct gkyl_array *fin, double dt, struct gkyl_array *out
  );
};

struct vm_boundary_fluxes {
  struct gkyl_rect_grid boundary_grid[2 * GKYL_MAX_CDIM];
  struct gkyl_array *flux_arr[2 * GKYL_MAX_CDIM];
  struct gkyl_array *mom_arr[2 * GKYL_MAX_CDIM];
  struct gkyl_range flux_r[2 * GKYL_MAX_CDIM];
  struct gkyl_range conf_r[2 * GKYL_MAX_CDIM];
  struct gkyl_mom_type *mom_type;
  struct gkyl_mom_calc *integ_moms[2 * GKYL_MAX_CDIM];
  gkyl_ghost_surf_calc *flux_slvr; // boundary flux solver
  double *red_mom; // reduction of the integrated moments of a boundary (on the device in GPU runs)
};

struct vm_emitting_wall {
  // emitting wall sheath boundary conditions
  int num_species;
  int dir;
  enum gkyl_edge_loc edge;
  double *scale_ptr;
  double t_bound;
  bool elastic;
  bool write;

  struct gkyl_spectrum_model *spectrum_model[GKYL_MAX_SPECIES];
  struct gkyl_yield_model *yield_model[GKYL_MAX_SPECIES];
  struct gkyl_elastic_model *elastic_model;
  struct gkyl_bc_emission_ctx *params;

  struct gkyl_bc_emission_spectrum *update[GKYL_MAX_SPECIES];
  struct gkyl_bc_emission_elastic *elastic_update;
  struct gkyl_array *f_emit;
  struct gkyl_array *f_emit_host;
  struct gkyl_array *buffer;
  struct gkyl_array *elastic_yield;
  struct gkyl_array *yield[GKYL_MAX_SPECIES]; // projected secondary electron yield
  struct gkyl_array *spectrum[GKYL_MAX_SPECIES]; // projected secondary electron spectrum
  struct gkyl_array *weight[GKYL_MAX_SPECIES];
  struct gkyl_array *flux[GKYL_MAX_SPECIES];
  struct gkyl_array *bflux_arr[GKYL_MAX_SPECIES];
  struct gkyl_array *k[GKYL_MAX_SPECIES];
  // Yield diagnostic: per impact species, the flux-weighted effective yield and the incoming
  // flux the emitted spectrum is normalized with, appended at every frame.
  struct gkyl_array *weight_host[GKYL_MAX_SPECIES];
  struct gkyl_array *flux_host[GKYL_MAX_SPECIES];
  gkyl_dynvec yield_diag;
  bool is_first_yield_write_call;
  struct vm_species *impact_species[GKYL_MAX_SPECIES]; // pointers to impacting species
  // Coordinate maps: the yield is evaluated on the impacting species' phase space, the spectrum
  // and elastic yield on the emitting species' phase space.
  struct vm_proj_c2p_ctx c2p_emit;
  struct vm_proj_c2p_ctx c2p_impact[GKYL_MAX_SPECIES];
  struct gkyl_range impact_normal_r[GKYL_MAX_SPECIES];
  // Moment type of each impact species (its basis, velocity map and Hamiltonian),
  // for the boundary-flux moments.
  struct gkyl_mom_type *mom_type[GKYL_MAX_SPECIES];
  struct gkyl_mom_calc *flux_slvr[GKYL_MAX_SPECIES]; // integrated moments

  struct gkyl_rect_grid *impact_grid[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_ghost_r[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_skin_r[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_buff_r[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_cbuff_r[GKYL_MAX_SPECIES];

  struct gkyl_rect_grid *emit_grid;
  struct gkyl_range *emit_buff_r;
  struct gkyl_range *emit_ghost_r;
  struct gkyl_range *emit_skin_r;
};

struct vm_source {
  bool write_source; // Are we writing out the source?
  bool evolve_source; // Are our sources time-dependent?

  bool calc_bflux; // boolean for if we are using boundary fluxes to rescale sources
  double scale_factor; // factor to scale source function
  double source_length; // length used to scale the source function
  struct vlasov_species *source_species; // species whose wall losses set the rate (of any type)

  bool rescale_m0; // boolean for if we are rescaling M0
  int num_cross_source; // how many other species are we obtaining sources from?
  struct vm_species *adapt_source_species[GKYL_MAX_SPECIES]; // list of species to use for the source
  int adapt_source_slot[GKYL_MAX_SPECIES]; // this species' slot in each partner's source_with list
    // (indexes the partner's scale_m0 in the cross rescale)
  struct gkyl_array
    *scale_m0[GKYL_MAX_SPECIES]; // Time-dependent re-scaling of the density of the source.
  struct gkyl_array *scale_m0_host[GKYL_MAX_SPECIES]; // host copy for use in IO
  struct gkyl_mom_calc *m0_reduced[GKYL_MAX_SPECIES]; // Reduced density update for rescaling source.
  struct gkyl_array *adapt_source[GKYL_MAX_SPECIES]; // adaptive source array
  int adapt_proj_source[GKYL_MAX_SPECIES]; // Index of projection function to use for adaptive source.

  bool filter; // boolean for if we are filtering rescaled M0
  int num_filters; // number of times to apply filter
  gkyl_dg_gaussian_filter *gauss_filter; // updater for filtering rescaled M0
  // Ghost cells of the rescaled M0 are filled before every filter pass: periodic directions are
  // synchronized and non-periodic directions copy the skin value into the ghosts.
  struct gkyl_bc_basic *filter_bc_lo[GKYL_MAX_CDIM]; // copy BC updaters, lower edges
  struct gkyl_bc_basic *filter_bc_up[GKYL_MAX_CDIM]; // copy BC updaters, upper edges
  struct gkyl_array *filter_bc_buffer; // buffer for the copy BC updaters

  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO
  struct gkyl_array *source_tmp; // temporary array for sources for accumulation if num_sources>1
  struct vm_proj proj_source[GKYL_MAX_PROJ]; // projector for source
  int num_sources; // Number of sources.

  int num_diag_moments; // number of diagnostics moments
  struct vm_species_moment *moms; // diagnostic moments
  struct vm_species_moment integ_moms; // integrated moments
  double *red_integ_diag; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time
};

// ---- PKPM species: sub-objects and the PKPM-only state of the two aspects ----
struct vlasov_species;
struct vm_species;
struct vm_fluid_species;

// data for moments
struct vm_pkpm_species_moment {
  struct gkyl_dg_updater_moment *mcalc; // moment update

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

// Context of the host-side projection of the source distributions: the p=1 configuration
// fields are evaluated at the computational coordinates of the phase-space quadrature points.
struct vm_pkpm_source_eval_ctx {
  int cdim;
  double mass;
  const struct gkyl_rect_grid *conf_grid;
  const struct gkyl_basis *conf_basis;
  const struct gkyl_range *local_ext;
  const struct gkyl_vlasov_projection *info; // LTE profile functions of the source
  const struct gkyl_array *prof; // [n_s, T_s, u_s] of the source (physical)
  const struct gkyl_array *rho, *u, *xi, *bvar; // nodal state of the stage (host)
};

struct vm_pkpm_source {
  enum gkyl_source_id source_id;
  bool calc_bflux; // rate from the wall losses of source_species
  double source_length;
  struct vlasov_species *source_species;
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
  struct vm_pkpm_source_eval_ctx eval_ctx;
};

struct vm_pkpm_bgk_collisions;

// Compute the node values [n, u (3), vtsq] of a species from a fluid state (with the Jacobian).
typedef void (*vm_coll_cross_prim_t)(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
);

// Compute the rate of the kinetic operator of a species with its partners at the nodes.
typedef void (*vm_coll_cross_moms_t)(gkyl_vlasov_app *app, struct vlasov_species *s);

// Collision frequencies and cross-species collisions of a species with a fluid aspect at the
// configuration-space nodes (CPU runs): the node values of its fluid state, and for each partner
// the frequencies and the coefficient of the momentum and energy exchange (struct
// vm_coll_exchange). A PKPM species (LBO or BGK) also has its self-collision frequency when it is
// computed from the local moments and the rate each partner adds to its kinetic operator; an
// Euler fluid has no operator.
struct vm_coll_cross {
  bool is_active; // the species computes its self frequency or has partners
  int num_cross_collisions; // number of species we cross-collide with
  struct vlasov_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  int other_idx[GKYL_MAX_SPECIES]; // position of this species in each partner's collide_with
  struct gkyl_array *prim_fluid; // [n, u (3), vtsq] of the species at the nodes
  struct gkyl_array *self_nu_ref; // constant of the computed self frequency
  struct gkyl_array *nu_ref[GKYL_MAX_SPECIES]; // reference cross frequencies at the nodes
  struct gkyl_array *nu_self; // frequency of this species on the current partner
  struct gkyl_array *nu_other; // frequency of the current partner on this species
  struct gkyl_array *coef[GKYL_MAX_SPECIES]; // exchange coefficient with each partner at the nodes
  struct gkyl_array *temp_ratio; // temperature before the exchange over that after it
  // Kinetic partners of a PKPM species (three velocity dimensions). Their cross operator carries
  // the momentum and energy of the pair, so they are not in the exchange: the species adds the
  // rate they give its kinetic operator, and its fluid state takes the opposite of what their
  // operator receives (kin_exch of the fluid aspect).
  int num_kin_partners; // number of kinetic species we cross-collide with
  struct vlasov_species *kin_partner[GKYL_MAX_SPECIES]; // those species
  struct gkyl_array *kin_nu_self[GKYL_MAX_SPECIES]; // frequency of this species on each, at the nodes
  struct gkyl_array
    *kin_nu_other[GKYL_MAX_SPECIES]; // frequency of each on this species, at the nodes
  const struct gkyl_array *kin_prim[GKYL_MAX_SPECIES]; // node values [n, u (3), vtsq] of each
  struct gkyl_array *nu_kin_self; // rate of the kinetic operator from self collisions
  struct gkyl_array *nu_kin; // rate of the kinetic operator (self plus partners)
  gkyl_pkpm_coll_nu *up_self_nu; // computed self frequency
  gkyl_pkpm_coll_nu *up_cross_nu; // cross frequencies: given, or computed
  gkyl_pkpm_coll_cross *up; // kinetic rate and exchange coefficient with a partner
  vm_coll_cross_prim_t calc_prim; // node values and self frequency of the stage
  vm_coll_cross_prim_t calc_step_prim; // node values and self frequency of the implicit stage
  vm_coll_cross_moms_t calc_moms; // rate of the kinetic operator of the stage
};

// Momentum and energy exchange of the species with a fluid aspect that have partners (PKPM
// species and Euler fluids): one backward-Euler step over all of them in each SSP-RK3 stage.
struct vm_coll_exchange {
  int num_active; // number of species with collision data (partners or a computed self frequency)
  struct vlasov_species *active[GKYL_MAX_SPECIES]; // those species
  int num_species; // number of species with partners
  struct vlasov_species *species[GKYL_MAX_SPECIES]; // those species
  struct gkyl_pkpm_coll_cross_species inp[GKYL_MAX_SPECIES]; // their data in the exchange
  gkyl_pkpm_coll_cross *up; // the backward-Euler step
  // the step on the given fluid states, or nothing without such species
  void (*exchange)(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[]);
};

// Accumulate the kinetic BGK term of a species onto rhs.
typedef void (*vm_pkpm_bgk_kinetic_rhs_t)(
  gkyl_vlasov_app *app, struct vm_species *kin, struct vm_pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

// Compute the shape factor of the implicit BGK step of a species.
typedef void (*vm_pkpm_bgk_shape_t)(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_bgk_collisions *bgk, double dt
);

struct vm_pkpm_bgk_collisions {
  // The operator at the configuration-space nodes.
  bool use_nodal; // false for GPU runs
  bool is_implicit; // backward Euler in each SSP-RK3 stage (nodal only)
  struct gkyl_array *rate_self; // self-collision frequency at the nodes (given or computed)
  struct gkyl_array *rate; // relaxation rate at the nodes (rate_self itself without partners)
  struct gkyl_array *shape; // implicit shape factor S at the nodes
  vm_pkpm_bgk_shape_t calc_shape; // shape factor of the implicit step (self, or with a partner)
  struct gkyl_array *f_lte; // corrected unit Maxwellian in w, 3 modes per velocity cell
  gkyl_pkpm_bgk_collisions *up_bgk; // updater
  vm_pkpm_bgk_kinetic_rhs_t kinetic_rhs; // the nodal or modal kinetic term

  // The modal form of the operator (GPU runs).
  struct gkyl_array *nu; // collision frequency
  struct gkyl_array *f_target; // J [F_t, G_t]
};

struct vm_pkpm_source_bgk;

// Compute the energy relaxation of the temperature source of a species for its fluid aspect, and
// the rate of the distributions.
typedef void (*vm_pkpm_source_bgk_calc_fluid_t)(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source_bgk *src,
  const struct gkyl_array *fluidin
);

// Accumulate the temperature source of a species onto the kinetic right-hand side.
typedef void (*vm_pkpm_source_bgk_kinetic_rhs_t)(
  gkyl_vlasov_app *app, struct vm_species *kin, struct vm_pkpm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs_f
);

// BGK-form temperature source, with its rate and energy relaxation at the configuration-space
// nodes.
struct vm_pkpm_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // type of source
  bool is_implicit; // backward Euler in each SSP-RK3 stage
  bool has_implicit_coll; // the implicit step includes the implicit BGK collisions of the species
  struct gkyl_array *rate; // source rate nu_Q at the nodes
  struct gkyl_array *vtsq; // source temperature over the mass, T_Q/m, at the nodes
  struct gkyl_array *shape; // rate (explicit) or shape factor (implicit) of the distributions
  struct gkyl_array *temp_ratio; // temperature before the implicit step over that after it
  struct gkyl_array *f_lte; // corrected unit Maxwellian in w, 3 modes per velocity cell
  gkyl_pkpm_temp_source *up_temp; // thermal-energy relaxation of the fluid state
  gkyl_pkpm_bgk_collisions *up_bgk; // relaxation of the distributions
  vm_pkpm_source_bgk_calc_fluid_t calc_fluid; // the explicit fluid term of the stage, or nothing
  vm_pkpm_source_bgk_kinetic_rhs_t kinetic_rhs; // the explicit kinetic term, or nothing
};

struct vm_pkpm_lbo_collisions;

// Compute the primitive moments of the self LBO of a species.
typedef void (*vm_pkpm_lbo_moms_t)(
  gkyl_vlasov_app *app, const struct vlasov_species *species, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

// Accumulate the self LBO of a species onto rhs.
typedef void (*vm_pkpm_lbo_kinetic_rhs_t)(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

struct vm_pkpm_lbo_collisions {
  // The operator at the configuration-space nodes.
  struct gkyl_array *nu; // collision frequency of the operator at the nodes
  gkyl_pkpm_lbo_moms *up_moms; // primitive moments at the nodes
  gkyl_dg_updater_lbo_pkpm_nodal *up_coll; // the operator at the nodes
  vm_pkpm_lbo_moms_t calc_moms; // the moments (constant or Spitzer frequency)
  vm_pkpm_lbo_kinetic_rhs_t kinetic_rhs; // the operator

  struct gkyl_array *prim_moms; // [u, vtsq] at the nodes
  struct gkyl_array *nu_self; // self-collision frequency at the nodes (given or computed)
};

// Records of the fluid positivity limiter (gkyl_dg_calc_pkpm_vars_positivity): the per-cell
// count [cells scaled for the density, for the internal energy, inadmissible averages], and
// the statistics per integrated-moment interval, which append the number of limiter
// applications.
enum vm_pkpm_positivity_record { PKPM_POS_NUM_COUNT = 3, PKPM_POS_NUM_STAT = 4 };

// PKPM-only state of the kinetic aspect of a PKPM species (struct vm_species): the distributions [F, G, M0_perp, N_perp] on the grid of
// the normalized velocity and everything the kinetic right-hand side reads
struct vm_pkpm_kinetic {
  struct gkyl_range
    local_ext_conf; // local range extended in configuration space only (moments, distribution BCs, sync)
  struct gkyl_array *blk_in,
    *blk_out; // one-block (F or G) scratch for configuration x phase products
  struct gkyl_array *f_no_J_host; // host copy of f_no_J for I/O
  // Operators that act in the velocity alone at each configuration-space node (collisions, the
  // temperature source) share one conversion of the distributions per stage and one conversion
  // of their summed increments back.
  bool has_nodal_ops; // true if the species has an operator at the nodes: the LBO (CPU runs)
  struct gkyl_array *f_nodal; // distributions in the configuration-nodal, velocity-modal form
  struct gkyl_array *rhs_nodal; // increment of the nodal operator, in the same form

  struct vm_pkpm_species_moment moms_diag; // for computing pkpm moments diagnostics
  // raw moments [M0, M2, M0_G, M1, M3, M1_G, mu (3), m1 (3), m2 (3), n0 (3)] of the normalized
  // distributions (mass 1): mu, m1, m2 = int (1, w, w^2) M0_perp dw and n0 = int N_perp dw
  struct vm_pkpm_species_moment nmoms;
  struct gkyl_array *
    cmoms; // the constraints [M0(F), M1(F), (M2(F) + 2 M0(G))/3, int M0_perp dw (x, y, z)], [1, 0, 1, 0, 0, 0] when the closures hold
  struct gkyl_array *cmoms_host; // host copy for output
  struct gkyl_array *
    grad_u; // gradient of u for the sources: in-cell plus the blended face lift, cdim x 3 expansions
  struct gkyl_array *theta; // T_perp/T = M0(G)
  // Copies of the flow velocity and of [xi, xi^2, p] of the fluid aspect, set with them
  // (vm_pkpm_species_calc_pkpm_vars).
  struct gkyl_array *u, *xi;
  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars; // Updater for the closures and the constraints
  struct gkyl_array *grad_xi; // gradient of [xi, ...] (grad_u layout; blocks (d, 0) hold d_d xi)
  struct gkyl_array *grad_ln_n; // grad(ln n), cdim expansions
  struct gkyl_array *c_0; // normalization closure (zeroth w moment of F)
  struct gkyl_array *alpha_0; // pressure force a_P/xi closing the first w moment
  struct gkyl_array *c_h; // heating rate closing the scalar second moment M2(F) + 2 M0(G)
  struct gkyl_array *a_perp; // perpendicular force closing int M0_perp dw [x, y, z]
  // [alpha_0, c] of the last stage of the last completed step: the closures of the first stage
  // of a step take the speeds of the velocity-domain edges from them. Written with the species.
  struct gkyl_array *closure_step;
  bool closure_from_step; // the stage is the first of a step (set with the moments of the stage)
  struct gkyl_array *closure_step_host; // host copy for output
  struct gkyl_array *rhs_moms; // moments of the configuration-pass right-hand side (separate scratch)
  struct gkyl_array
    *conf_flux_surf; // stored modal configuration-space kinetic fluxes, [dir][side][block][mode]
  struct gkyl_array *flux_blend; // weight of the label-paired flux of each lower face, [dir]
  struct gkyl_array *vel_flux_surf; // stored modal v_par kinetic fluxes
  struct gkyl_dg_pkpm_harm_vol
    *harm_vol; // modal volume integrals of the first-harmonic terms (configuration pass)
  struct gkyl_array *harm_coef; // stored operator of the stage (gkyl_pkpm_harm_ops.h)
  struct gkyl_array *harm_wnode; // raw w-flux values of the stage at the w-face nodes
  struct gkyl_dg_pkpm_harm_exp *harm_exp; // gyrophase-mixing integrating factor of the first harmonic
  struct gkyl_array *harm_exp_coef; // its field (q/m) B frozen at the start of the step
  struct gkyl_array *harm_exp_nu; // gyrophase-mixing rate at the phase nodes, frozen per step
  struct gkyl_array *harm_exp_act; // configuration cells with a damped node (momentum restoration)
  struct gkyl_array *harm_exp_col; // w integrals of F and M0_perp of the restoration
  struct gkyl_dg_pkpm_remap_flux_surf
    *calc_face_flux; // configuration-space face flux updater (kinetic and fluid)
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vel_flux; // v_par flux updater
  struct gkyl_dg_eqn *eqn_kinetic_conf; // equation object of the configuration pass
  struct gkyl_dg_eqn *eqn_kinetic_vel; // equation object of the velocity pass
  gkyl_hyper_dg *slvr_kinetic_conf; // hyper_dg over the configuration directions
  gkyl_hyper_dg *slvr_kinetic_vel; // hyper_dg over v_par
  long
    num_null_closure; // accumulated count of closure nodes nulled (M0 <= 0 or a non-positive determinant)
  double closure_cond_min; // smallest det/det0 of the closure solves over the run
  bool closure_null_warned; // the velocity-extent warning was printed
  // Diagnostic output: the variables
  // [u (3), T_perp/m, a_P, div b, bb:grad u, T_perp/m div b, xi, c]
  struct gkyl_array *vars_io;
  struct gkyl_array *vars_io_host;
  int num_face_nodes; // Transverse Gauss-Legendre nodes of a configuration face, per direction
  // reflecting walls of the distributions (0 at other walls): the mirror image in the ghost
  // cell's field-aligned coordinates, from b in the skin and ghost cells
  struct gkyl_bc_pkpm_species_reflect *bc_lo_reflect[3];
  struct gkyl_bc_pkpm_species_reflect *bc_up_reflect[3];
  // absorbing walls (0 at other walls): theta mirrored (even image) into the ghost cells, whose
  // empty G has no T_perp/T of its own, so that its face averages are the interior wall traces
  struct gkyl_bc_basic *bc_lo_theta[3], *bc_up_theta[3];
  bool bc_is_absorb; // whether any boundary is absorbing (wall fluxes from the kinetic upwind flux)

  enum gkyl_collision_id collision_id; // type of collisions
  struct vm_pkpm_lbo_collisions lbo; // collisions object
  struct vm_pkpm_bgk_collisions bgk; // BGK collisions object
  struct vm_coll_cross cross; // cross-species collisions (LBO or BGK)
  // Stage combination of a species whose implicit kinetic relaxation carries a temperature change
  // (implicit BGK with partners, implicit temperature source): the stages combine the
  // distributions with weights that are linear in the physical second moments.
  bool has_stage_weights;
  struct gkyl_array *stage_lam; // temperature factor of the implicit steps of this step, at the nodes
  gkyl_pkpm_stage_weights *up_stage_weights; // the weighted combination

  struct vm_pkpm_source src; // source object
  struct vm_pkpm_source_bgk src_bgk; // temperature source object
  double bflux
    [2 *
     GKYL_MAX_CDIM]; // particle loss rate through the absorbing walls [lower, upper] per direction (stage)
  gkyl_dynvec bflux_diag; // the loss rates on the cadence of the integrated moments
  bool is_first_bflux_write_call;
};

// fluid aspect of a PKPM species: the fluid state [rho, rho ux, rho uy, rho uz, E] and
// everything the fluid right-hand side reads
struct vm_pkpm_fluid {
  struct gkyl_array *
    moms; // kinetic coupling moments [p_par, p_perp, q_par, q_perp, s (3), r (3)] = [p M2, p M0_G, p xi M3, p xi M1_G, p m1, p xi (m2/2 + n0)]

  // Configuration-space variables of the update (p=1).
  struct gkyl_array *u; // flow velocity [ux, uy, uz], nodal quotient M/rho (the weak division at p=1)
  struct gkyl_array *xi; // thermal speed and scalar pressure [xi, xi^2, p] of the fluid state (nodal)
  struct gkyl_array *xi3; // copy of xi for the weak-gradient kernels (only the first block is used)
  struct gkyl_array *scratch; // one-expansion scratch for weak products
  struct gkyl_array *
    fluid_no_J; // PHYSICAL fluid state [rho, M, E] of the stage (position-map Jacobian divided out)
  struct gkyl_array *fluid_no_J_host; // host copy of fluid_no_J for I/O
  struct gkyl_array *moms_no_J; // physical coupling moments (diagnostics)
  struct gkyl_array *rho; // physical mass density (leading block of fluid_no_J)
  struct gkyl_array *rho3; // [rho, 0, 0] for the weak-gradient kernels
  struct gkyl_array
    *mom3; // momentum block of the fluid array (scratch for the 3-component consumers)
  struct gkyl_array *grad_rho; // gradient of [rho, 0, 0] (scratch)
  struct gkyl_array *fluid_flux_surf; // stored modal fluid fluxes [F_rho, F_x, F_y, F_z, F_E]

  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars; // Updater for the configuration-space variables
  struct gkyl_dg_eqn *eqn_fluid; // fluid equation object
  gkyl_hyper_dg *slvr_fluid; // hyper_dg for the fluid equations
  struct gkyl_array
    *pos_count; // per-cell record of the positivity limiter [rho, e, avg] (cumulative)
  struct gkyl_array *pos_io,
    *pos_io_host; // the record as cell averages of a p=1 field for IO
  long pos_count_since
    [PKPM_POS_NUM_STAT]; // limiter counts [rho, e, avg, applications] since the last integrated-moment call
  // both records as saved by vm_pkpm_species_positivity_record_save
  struct gkyl_array *pos_count_saved;
  long pos_count_since_saved[PKPM_POS_NUM_STAT];
  gkyl_dynvec pos_stat; // the counts per integrated-moment interval, summed over ranks
  bool is_first_pos_stat_write_call; // flag for the positivity statistics dynvec written first time

  // Diagnostic output: the pressure tensor and heat flux [P_ij (6), q_i (3)] of the fluid equations
  struct gkyl_array *prim_io;
  struct gkyl_array *prim_io_host;
  struct gkyl_array *fluid_flux_surf_host; // host mirror of the stored fluid fluxes (GPU only)
  struct gkyl_dg_pkpm_mom_diffusion *mom_diff; // momentum diffusion on the position map

  // What the fluid right-hand side needs from the sources of the kinetic aspect, set by the
  // source objects and the staging phases: the fluid image J [rho, M, E] of the source profile
  // with the rate of the stage, and the energy relaxation of the explicit temperature source.
  const struct gkyl_array *src_fluid; // fluid image of the source profile (0 without a source)
  double src_scale; // rate factor of the stage (0 without a source)
  bool has_src_bgk_rhs; // flag to indicate there is an explicit temperature source
  struct gkyl_array *src_bgk_rhs; // its energy relaxation of the stage
  // Kinetic collision partners: the momentum and energy their cross operators receive from this
  // species per unit time, [d(m M1)/dt (3), d(m M2)/dt]; the right-hand side takes the opposite.
  int num_kin_exch;
  const struct gkyl_array *kin_exch[GKYL_MAX_SPECIES];
};

// ---- kinetic and fluid aspect objects ----
// species data
struct vm_species {
  struct gkyl_vlasov_kinetic_species info; // data for species

  char name[128]; // Species name.
  double charge, mass; // Charge and mass.

  struct gkyl_basis basis; // Phase-space basis.
  struct gkyl_basis basis_vel; // Velocity-space basis.
  struct gkyl_basis basis_surf; // Surface basis at velocity-space surfaces for velocity-space fluxes.
  struct gkyl_basis *basis_on_dev; // Pointer to phase-space basis on device.

  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext phase-space ranges

  struct gkyl_comm *comm; // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)
  // buffers of the fixed-function BCs per direction (allocated for the edges that have one)
  struct gkyl_array *bc_buffer_lo_fixed[3], *bc_buffer_up_fixed[3];

  struct gkyl_array *f_host; // Host-side distribution function for I/O on GPUs.
  struct gkyl_array *cflrate_host; // Host-side cflrate for I/O on GPUs.
  bool write_cell_avg; // Boolean for only writing cell average of f.

  enum gkyl_field_id field_id; // Type of field equation.
  enum gkyl_model_id model_id; // type of Vlasov equation (e.g., non-relativistic vs. relativistic).
  enum gkyl_hamil_id hamil_id; // Representation of the equation-object Hamiltonian (hamil):
    // sparse/dense velocity-space expansion or phase-space expansion.
  enum gkyl_hamil_id mom_hamil_id; // Representation of the moment Hamiltonian (mom_hamil).
    // Differs from hamil_id for triad-GR: the equation Hamiltonian
    // is a phase-space expansion while moments use the dense
    // velocity-space (SR) Hamiltonian.
  enum gkyl_triad_preset_geom_type
    triad_preset_geom_type; // geom type for preset geometries for triads

  struct gkyl_vlasov_velocity_map
    *vel_map; // Velocity-space mapping object (owns all velocity map arrays).
  struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space mapping object (acquired reference to the app's pos_map).

  struct gkyl_array *f_no_J; // Distribution function without velocity-space Jacobian.
    // When using uniform velocity-space mesh, just stores the distribution function at that RK stage.

  // Organization of the different equation objects and the required data.
  bool has_rad; // Do we have a radiation drag force?
  struct gkyl_array *rad; // array for radiation drag force.

  // Geometry
  struct vm_geom *geom; // Geometry structure for vm

  // Organization of the different equation objects and the required data.
  struct gkyl_range
    mom_hamil_range; // Range Hamiltonian (for moments) is defined over (only velocity-space or all phase-space).
  struct gkyl_array
    *mom_hamil; // Specified Hamiltonian (for moments) function for canonical poisson bracket.
  struct gkyl_range
    hamil_range; // Range Hamiltonian is defined over (only velocity-space or all phase-space).
  struct gkyl_array *hamil; // Specified Hamiltonian function for canonical poisson bracket.
  struct gkyl_array *hamil_host; // Host-side Hamiltonian array for initial projection.
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation
  // Poisson tensor read by the configuration-space surface flux: the face-node storage
  // of the nodal triad kernels, otherwise the same array as conf_poisson_tensor.
  struct gkyl_array *conf_poisson_tensor_surf;
  struct gkyl_array
    *conf_poisson_tensor_host; // Host-side configuration space Poisson tensor representation
  union {
    // Special relativistic Vlasov-Maxwell model.
    struct {
      struct gkyl_array *gamma_inv; // array for 1/gamma = 1.0/sqrt(1 + p^2)
      struct gkyl_dg_calc_sr_vars *sr_vars; // updater for computing SR variables
    };
    // Canonical Poisson Bracket using specified Hamiltonian in phase space.
    struct {
      struct gkyl_array *h_ij; // Specified metric (covariant) for canonical poisson bracket
      struct gkyl_array *h_ij_host; // Host side metric array for initial projection
      struct gkyl_array *h_ij_inv; // Specified metric inverse for canonical poisson bracket
      struct gkyl_array *h_ij_inv_host; // Host side metric inverse array for initial projection
      struct gkyl_array *det_h; // Specified metric determinant
      struct gkyl_array *det_h_host; // Host side metric determinant
      struct gkyl_array *background_flows; // Specified background flows
      struct gkyl_array *background_flows_host; // Host side background flows
      struct gkyl_array *effective_potential; // Specified effective potential
      struct gkyl_array *effective_potential_host; // Host side effective potential
    };
  };

  struct vm_species_moment m1i; // for computing currents
  struct vm_species_moment m0; // for computing charge density
  struct vm_species_moment integ_moms; // integrated moments
  struct vm_species_moment *moms; // diagnostic moments
  struct gkyl_array *L2_f; // L2 energy f^2
  double *red_L2_f; // for reduction of integrated L^2 energy on GPU
  double *red_integ_diag; // for reduction of integrated moments on GPU
  gkyl_dynvec integ_L2_f; // integrated L^2 energy reduced across grid
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_L2_write_call; // flag for integrated L^2 energy dynvec written first time
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  // boundary conditions on lower/upper edges in each direction
  struct gkyl_vlasov_bc lower_bc[3], upper_bc[3];
  // emitting wall sheath boundary conditions
  struct vm_emitting_wall bc_emission_lo;
  struct vm_emitting_wall bc_emission_up;
  bool emit_lo; // flag to indicate if there emission BCs
  bool emit_up;
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges, valid (i.e. volume>0) in ranks abutting boundaries.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];

  int num_init; // Number of initial condition functions.
  struct vm_proj proj_init[GKYL_MAX_PROJ]; // projectors for initial conditions

  struct vm_lte lte; // object needed for the lte equilibrium

  struct vm_collisionless collisionless; // Collisionless terms.

  struct vm_lbo_collisions lbo; // LBO collisions object
  struct vm_bgk_collisions bgk; // BGK collisions object
  struct vm_source_bgk src_bgk; // temperature source object

  bool calc_bflux; // are we computing boundary fluxes?
  struct vm_boundary_fluxes bflux; // boundary flux object

  enum gkyl_source_id source_id; // type of source
  struct vm_source src; // applied source

  struct vm_pkpm_kinetic pkpm; // PKPM-only state (PKPM species)

  double *omega_cfl;

  // Pointer to various functions selected at runtime.
  void (*apply_ic_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm);
  double (*rhs_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, const struct gkyl_array *fin,
    const struct gkyl_array *em, struct gkyl_array *rhs
  );
  double (*rhs_implicit_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, const struct gkyl_array *fin,
    struct gkyl_array *rhs, double dt
  );
  // Moments of the implicit update over dt (computed for every species before
  // any rhs_implicit_func, which reads the moments of the collision partners).
  void (*calc_implicit_moms_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, const struct gkyl_array *fin, double dt
  );
  // Source phases of a forward-Euler stage, run after every species' RHS:
  // moments that rescale adaptive sources, then the source onto the RHS.
  void (*calc_source_moms_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, const struct gkyl_array *fin
  );
  void (*rhs_source_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, double tcurr, struct gkyl_array *rhs
  );
  void (*bc_func)(
    gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f, double tcurr
  );
  // Read the distribution from a file and rebuild what the file does not carry.
  struct gkyl_app_restart_status (*from_file_func)(
    gkyl_vlasov_app *app, struct vm_species *vms, const char *fname
  );
  void (*release_func)(const gkyl_vlasov_app *app, const struct vm_species *vms);
  void (*step_f_func)(struct gkyl_array *out, double dt, const struct gkyl_array *inp);
  void (*combine_func)(
    struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
    const struct gkyl_array *arr2, const struct gkyl_range *rng
  );
  void (*copy_func)(
    struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
  );
  void (*write_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
  void (*write_lte_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
  void (*write_cfl_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
  void (*write_cell_avg_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
  void (*write_mom_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm);
  void (*write_integrated_mom_func)(gkyl_vlasov_app *app, struct vm_species *vms);
  void (*calc_L2_func)(gkyl_vlasov_app *app, struct vm_species *vms, double tm);
  void (*write_L2_func)(gkyl_vlasov_app *app, struct vm_species *vms);
};

struct vm_fluid_source {
  struct vm_species_moment moms; // source moments

  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *source_proj; // projector for source
};

// fluid species data
struct vm_fluid_species {
  struct gkyl_vlasov_fluid_species info; // data for fluid

  char name[128]; // Species name.
  double charge, mass; // Charge and mass.

  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_array *fluid, *fluid1, *fluidnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)

  struct gkyl_array *fluid_host; // host copy for use IO and initialization

  enum gkyl_eqn_type eqn_type; // type ID of equation
  int num_equations; // number of equations in species
  struct gkyl_wv_eqn *equation; // equation object
  bool has_poisson;
  // organization of the different equation objects and the required data and solvers
  union {
    // Applied advection
    struct {
      struct gkyl_array *app_advect; // applied advection
      struct gkyl_array *app_advect_host; // host copy for use in IO and projecting
    };
    // Euler/Isothermal Euler
    struct {
      // u: flow velocity (ux, uy, uz). p: pressure, vth^2*rho for isothermal Euler and
      // (gamma - 1)*(E - 1/2 rho u^2) for Euler.
      // Also a prim_vars and prim_vars_host array for I/O of (u,p)
      struct gkyl_array *u;
      struct gkyl_array *p;
      struct gkyl_array *prim_vars;
      struct gkyl_array *prim_vars_host;
      struct gkyl_array *
        cell_avg_prim; // Integer array for whether e.g., rho *only* uses cell averages for weak division
      // Determined when constructing the matrix if rho < 0.0 at control points

      struct gkyl_array *u_surf;
      struct gkyl_array *p_surf;
      struct gkyl_dg_calc_fluid_vars
        *calc_fluid_vars; // Updater to compute fluid variables (flow velocity and pressure)
      struct gkyl_dg_calc_fluid_vars
        *calc_fluid_vars_ext; // Updater to compute fluid variables (flow velocity and pressure)
      // over extended range (used when BCs are not absorbing to minimize apply BCs calls)
    };
    // Canonical PB Fluid such as incompressible Euler or Hasegawa-Wakatani
    struct {
      struct gkyl_array *
        phi; // potential determined by canonical PB Poisson equation on local range used by updater
      struct gkyl_array *
        phi_global; // potential determined by canonical PB Poisson equation on global range given by Poisson solver
      struct gkyl_array *
        poisson_rhs_global; // global RHS of Poisson equation, simply an all-gather of, e.g., the vorticity
      struct gkyl_array *phi_host; // host copy for use IO
      struct gkyl_array
        *can_pb_n0; // background density gradient for driving turbulence in some fluid systems.
      // Rate of the adiabatic (alpha) coupling for the slowest box mode, alpha (1 + k_min^2)/k_min^2,
      // seeding the CFL rate so the explicit RK stepping of this stiff mode stays stable.
      double can_pb_alpha_cfl_rate;
      struct gkyl_array *
        epsilon; // Permittivity in Poisson equation, set to -1.0 for canonical PB Poisson equations.
      struct gkyl_array *
        kSq; // k^2 factor in Helmholtz equation needed for Hasegawa-Mima where we solve (grad^2 - 1) phi = RHS

      struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
        // for solving Poisson equation on each MPI process in parallel

      struct gkyl_fem_poisson
        *fem_poisson; // Poisson solver for - nabla . (epsilon * nabla phi) - kSq * phi = rho.

      struct gkyl_array
        *alpha_surf; // Surface configuration space velocity (derivatives of potential, phi)
      struct gkyl_array *sgn_alpha_surf; // sign(alpha_surf) at quadrature points
      struct gkyl_array
        *const_sgn_alpha; // boolean for if sign(alpha_surf) is a constant, either +1 or -1
      struct gkyl_dg_calc_canonical_pb_fluid_vars
        *calc_can_pb_fluid_vars; // Updater for computing surface alpha and sources.
      struct gkyl_array
        *can_pb_energy_fac; // Factor in calculation of canonical PB energy diagnostic.
      struct gkyl_array_integrate *calc_can_pb_energy;
      double *red_can_pb_energy,
        *red_can_pb_energy_global; // Memory for use in GPU reduction of canonical PB energy.
    };
  };

  struct gkyl_dg_updater_fluid *advect_slvr; // Fluid equation solver

  // fluid diffusion
  bool has_diffusion; // flag to indicate there is applied diffusion
  struct gkyl_array *diffD; // array for diffusion tensor
  struct gkyl_dg_updater_diffusion_fluid *diff_slvr; // Fluid equation solver
  struct gkyl_dg_updater_diffusion_gen *diff_slvr_gen;

  // boundary conditions on lower/upper edges in each direction
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // buffers of the fixed-function BCs per direction (allocated for the edges that have one)
  struct gkyl_array *bc_buffer_lo_fixed[3], *bc_buffer_up_fixed[3];
  bool bc_is_absorb; // boolean for absorbing BCs since 1/rho is undefined in absorbing BCs
    // If BCs are *not* absorbing, primitive variables can be calculated on *extended* range

  struct gkyl_array *integ_mom; // Integrated moments
  double *red_integ_diag; // for reduction on GPU
  gkyl_dynvec integ_diag; // Integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for int-moments dynvec written first time

  bool has_app_accel; // flag to indicate there is applied acceleration
  bool app_accel_evolve; // flag to indicate applied acceleration is time dependent
  struct gkyl_array *app_accel; // applied acceleration
  struct gkyl_array *app_accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_accel_proj; // projector for acceleration

  struct vm_coll_cross coll_cross; // collisions of a fluid species with other species (exchange)

  // Vertex positivity limiter of the fluid state (Euler and isothermal Euler at polynomial order
  // 1; a no-op otherwise).
  double positivity_eps; // relative floor of the vertex density (and internal energy for Euler)
  void (*positivity_func)(
    gkyl_vlasov_app *app, struct vm_fluid_species *f, struct gkyl_array *fluid
  );

  // Forces in the fluid's own RHS: the applied acceleration without a field, and with it the
  // Lorentz force of a Maxwell field coupled explicitly (use_explicit_sources). None with the
  // implicit fluid-EM coupling, which applies both.
  bool has_explicit_force; // the RHS of the fluid applies its forces
  struct gkyl_array *qmem; // [a, q/m B]: acceleration (applied + q/m E) and q/m times B (NULL
    // unless the field is coupled explicitly)
  struct gkyl_array *mom; // momentum density, for the current of the explicit coupling (same)
  void (*rhs_force_func)(
    gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid,
    const struct gkyl_array *em, struct gkyl_array *rhs
  );

  // fluid source
  enum gkyl_source_id source_id; // type of source
  struct vm_fluid_source src; // applied source

  struct vm_pkpm_fluid pkpm; // PKPM-only state (PKPM species)

  double *omegaCfl_ptr;

  // Pointers to functions selected at runtime.
  void (*apply_ic_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, double t0);
  double (*rhs_func)(
    gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid,
    const struct gkyl_array *em, struct gkyl_array *rhs
  );
  // Source onto the RHS, run after every species' RHS.
  void (*rhs_source_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, struct gkyl_array *rhs);
  void (*step_f_func)(struct gkyl_array *out, double dt, const struct gkyl_array *inp);
  void (*combine_func)(
    struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
    const struct gkyl_array *arr2, const struct gkyl_range *rng
  );
  void (*copy_func)(
    struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
  );
  void (*bc_func)(gkyl_vlasov_app *app, const struct vm_fluid_species *f, struct gkyl_array *fluid);
  // Slope and positivity limiters of a stepped state, then the BCs (set by the equation type).
  void (*limiter_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, struct gkyl_array *fluid);
  // Read the fluid state from a file and rebuild what the file does not carry.
  struct gkyl_app_restart_status (*from_file_func)(
    gkyl_vlasov_app *app, struct vm_fluid_species *f, const char *fname
  );
  // Function pointers for computing primitive/auxiliary variables,
  // and also write method, release method, and method for calculating integrated quantities.
  void (*prim_vars_func)(
    gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid
  );
  void (*calc_integrated_mom_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm);
  void (*write_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm, int frame);
  void (*release_func)(const gkyl_vlasov_app *app, struct vm_fluid_species *f);
};

// RK buffers of a species: the solution (f, fluid), the first-stage buffer
// (f1, fluid1) and the new-stage buffer (fnew, fluidnew).
enum vm_rk_buf { VM_RK_F, VM_RK_F1, VM_RK_FNEW };

// ---- species container ----
// A species owns a kinetic aspect (kinetic), a fluid aspect (fluid) or, a PKPM
// species, both; an absent one is NULL. Per-aspect operations apply to the
// aspects present; the staging phases, the explicit field coupling, the
// diagnostics, output and restart of the species as a whole and the release are
// dispatched through the function pointers below, set by the constructors.
struct vlasov_species {
  enum gkyl_species_type type; // Which aspects this species owns.

  // Identity, copied into each aspect at construction and never changed after.
  char name[128]; // Species name.
  double charge, mass; // Charge and mass.

  struct vm_species *kinetic; // kinetic aspect (NULL if absent)
  struct vm_fluid_species *fluid; // fluid aspect (NULL if absent)

  // Staging phases that fill the pre-RHS auxiliary arrays. calc_self_moms:
  // moments of the species' own stage state (kinetic: LBO moments and boundary
  // corrections, explicit BGK moments, target moments of the explicit temperature
  // source; fluid: nothing). The stage distribution is passed non-const; no
  // species modifies it. calc_coupled_vars: quantities that read other species or
  // the fluid state, i.e. LBO cross-collision moments (kinetic) and primitive
  // variables (fluid); it runs after every species' self moments. BGK cross
  // moments are formed inside the BGK RHS.
  void (*calc_self_moms_func)(
    gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
    const struct gkyl_array *fluidin
  );
  void (*calc_coupled_vars_func)(
    gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
    const struct gkyl_array *fluidin
  );
  // Update of the stepped state at the end of a forward-Euler stage, after
  // step_f and before the boundary conditions (empty for kinetic and fluid
  // species).
  void (*finish_stage_func)(
    gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, double dt, struct gkyl_array *fout,
    struct gkyl_array *fluidout
  );
  // Accumulate the species' coupling onto the field's target array. Maxwell: the
  // current onto the EM RHS (kinetic species; Euler fluids with
  // use_explicit_sources). Poisson: the weighted number density onto the source
  // of each potential (target unused). Empty for the null field and all other
  // fluid species.
  void (*accumulate_field_coupling_func)(
    gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
    const struct gkyl_array *fluidin, struct gkyl_array *target
  );
  // Completion of a time step, once its last stage has been accepted (empty for
  // kinetic and fluid species).
  void (*complete_step_func)(gkyl_vlasov_app *app, struct vlasov_species *sp);
  // Diagnostics, output and restart of the species as a whole (those of the aspect
  // of a kinetic or of a fluid species; a PKPM species forms its diagnostics from
  // both aspects and names the files of both).
  void (*calc_integrated_mom_func)(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm);
  void (*write_func)(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame);
  void (*write_mom_func)(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame);
  void (*write_integrated_mom_func)(gkyl_vlasov_app *app, struct vlasov_species *sp);
  struct gkyl_app_restart_status (*read_from_frame_func)(
    gkyl_vlasov_app *app, struct vlasov_species *sp, int frame
  );
  // Collision data of a species with a fluid aspect (that of the kinetic aspect of a PKPM species,
  // of the fluid aspect of a fluid species; NULL for a kinetic species), for its partners.
  struct vm_coll_cross *coll_cross;
  // Number of particles per unit time the species loses through its walls, from the boundary
  // fluxes of the stage (kinetic species: known once its right-hand side is; PKPM species: once
  // its self moments are; zero for a fluid species).
  double (*wall_loss_rate_func)(gkyl_vlasov_app *app, const struct vlasov_species *sp);
  // Release the aspects of the species and the container's allocations.
  void (*release_func)(const gkyl_vlasov_app *app, struct vlasov_species *sp);
};

/** vlasov_species API: type-agnostic operations on the species container (vlasov_species.c). */

/**
 * Construct a species from its input: validate the declared type against the
 * input blocks (a kinetic species needs a velocity grid and no fluid equation;
 * a fluid species needs an equation and no velocity grid; a PKPM species needs
 * a velocity grid and no fluid equation) and dispatch to the typed constructor.
 * The field must already be constructed.
 *
 * @param app Vlasov app object
 * @param inp Species input
 * @param sp Species container to construct
 */
void vlasov_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
);

/**
 * Construct a kinetic species: allocate the kinetic aspect, store its input
 * block and identity, and set the container's methods. The aspect is then
 * initialized by vm_species_init.
 *
 * @param app Vlasov app object
 * @param inp Species input
 * @param sp Species container to construct
 */
void vlasov_kinetic_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
);

/**
 * Construct a fluid species: allocate the fluid aspect, store its input block
 * and identity, and set the container's methods. The aspect is then initialized
 * by vm_fluid_species_init.
 *
 * @param app Vlasov app object
 * @param inp Species input
 * @param sp Species container to construct
 */
void vlasov_fluid_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
);

/**
 * Construct a PKPM species: allocate both aspects, store their input blocks and
 * the identity, and set the container's methods. The aspects are then
 * initialized together by vm_pkpm_species_init (vlasov_species_init_pkpm).
 *
 * @param app Vlasov app object
 * @param inp Species input
 * @param sp Species container to construct
 */
void vlasov_pkpm_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
);

/**
 * Compute the moments of the species' own stage state needed by its RHS: LBO
 * moments and boundary corrections, explicit BGK moments, target moments of the
 * explicit temperature source (no-op for fluid species; staging phase 1).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param stage SSP-RK3 stage (1, 2 or 3)
 * @param fin Stage distribution function (NULL if no kinetic aspect)
 * @param fluidin Stage fluid moments (NULL if no fluid aspect)
 */
void vlasov_species_calc_self_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Compute the coupled variables (staging phase 2): LBO cross-collision moments
 * (kinetic) or primitive variables (fluid). Must run after phase 1 for all
 * species.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fin Input distribution function (NULL if no kinetic aspect)
 * @param fluidin Input fluid moments (NULL if no fluid aspect)
 */
void vlasov_species_calc_coupled_vars(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Accumulate the species' coupling onto the field's target array. Maxwell: the
 * current onto emout (kinetic species; Euler fluids with use_explicit_sources).
 * Poisson: the weighted number density onto the source of each potential the
 * field solves for. Empty for the null field and all other fluid species.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fin Input distribution function (NULL if no kinetic aspect)
 * @param fluidin Input fluid moments (NULL if no fluid aspect)
 * @param target Field array to accumulate onto (unused for Poisson)
 */
void vlasov_species_accumulate_field_coupling(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
);

/**
 * Compute the time-dependent applied acceleration of each present aspect at a
 * stage time (no-op for an aspect whose acceleration does not evolve, and for
 * a fluid aspect whose acceleration the implicit fluid-EM coupling applies).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param tcurr Current time
 */
void vlasov_species_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr);

/**
 * Compute the RHS of each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fin Input distribution function (NULL if no kinetic aspect)
 * @param fluidin Input fluid moments (NULL if no fluid aspect)
 * @param emin Input EM field
 * @param fout On output, the distribution RHS
 * @param fluidout On output, the fluid RHS
 * @return Maximum stable time-step across aspects
 */
double vlasov_species_rhs(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, const struct gkyl_array *emin, struct gkyl_array *fout,
  struct gkyl_array *fluidout
);

/**
 * Compute the moments needed to rescale adaptive sources (kinetic aspect).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fin Input distribution function (NULL if no kinetic aspect)
 */
void vlasov_species_calc_source_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin
);

/**
 * Accumulate the source onto the RHS of each present aspect. Run after the RHS
 * of every species: a boundary-flux source reads the boundary fluxes its
 * source species computed in its RHS.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param tcurr Current time
 * @param fout Distribution RHS (NULL if no kinetic aspect)
 * @param fluidout Fluid RHS (NULL if no fluid aspect)
 */
void vlasov_species_source_rhs(
  gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr, struct gkyl_array *fout,
  struct gkyl_array *fluidout
);

/**
 * Forward-Euler accumulate, out = dt*out + in, for each present aspect.
 *
 * @param sp Species container
 * @param dt Time-step
 * @param fin Input distribution function (NULL if no kinetic aspect)
 * @param fluidin Input fluid moments (NULL if no fluid aspect)
 * @param fout Distribution to update
 * @param fluidout Fluid moments to update
 */
void vlasov_species_step_f(
  struct vlasov_species *sp, double dt, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *fout, struct gkyl_array *fluidout
);

/**
 * Finish a forward-Euler stage of a species: update of its stepped state after
 * step_f and before the boundary conditions (empty for kinetic and fluid
 * species).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param stage SSP-RK3 stage (1, 2 or 3)
 * @param dt Time-step taken
 * @param fout Stepped distribution function (NULL if no kinetic aspect)
 * @param fluidout Stepped fluid moments (NULL if no fluid aspect)
 */
void vlasov_species_finish_stage(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, double dt, struct gkyl_array *fout,
  struct gkyl_array *fluidout
);

/**
 * Complete a time step of a species, once its last stage has been accepted and
 * before the stage state is copied into the solution (empty for kinetic and
 * fluid species).
 *
 * @param app Vlasov app object
 * @param sp Species container
 */
void vlasov_species_complete_step(gkyl_vlasov_app *app, struct vlasov_species *sp);

/**
 * Combine RK buffers, out = c1*b1 + c2*b2, for each present aspect. The
 * stepper selects the buffers.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param out Output buffer
 * @param c1 Coefficient of b1
 * @param b1 First input buffer
 * @param c2 Coefficient of b2
 * @param b2 Second input buffer
 */
void vlasov_species_combine(
  gkyl_vlasov_app *app, struct vlasov_species *sp, enum vm_rk_buf out, double c1, enum vm_rk_buf b1,
  double c2, enum vm_rk_buf b2
);

/**
 * Copy one RK buffer into another, out = inp, for each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param out Output buffer
 * @param inp Input buffer
 */
void vlasov_species_copy_range(
  gkyl_vlasov_app *app, struct vlasov_species *sp, enum vm_rk_buf out, enum vm_rk_buf inp
);

/**
 * Gather the RK-stage input/output arrays of every species, indexed by species
 * (NULL where a species lacks that aspect).
 *
 * @param app Vlasov app object
 * @param in Buffer this stage reads
 * @param out Buffer this stage writes
 * @param fin On output, the input distributions
 * @param fout On output, the output distributions
 * @param fluidin On output, the input fluid moments
 * @param fluidout On output, the output fluid moments
 */
void vlasov_species_gather_rk_state(
  gkyl_vlasov_app *app, enum vm_rk_buf in, enum vm_rk_buf out, struct gkyl_array *fin[],
  struct gkyl_array *fout[], const struct gkyl_array *fluidin[], struct gkyl_array *fluidout[]
);

/**
 * Apply boundary conditions to each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param f Distribution function (NULL if no kinetic aspect)
 * @param fluid Fluid moments (NULL if no fluid aspect)
 * @param tcurr Current time
 */
void vlasov_species_apply_bc(
  gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *f, struct gkyl_array *fluid,
  double tcurr
);

/**
 * Apply the limiters to the stepped fluid state (fluid aspect).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fluid Fluid moments to limit (NULL if no fluid aspect)
 */
void vlasov_species_limiter(
  gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fluid
);

/**
 * Implicit update of a stage (BGK collisions and the temperature source), phase 1:
 * moments of the input distribution.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fin Input distribution function (NULL if no kinetic aspect)
 * @param dt Time-step of the implicit update
 */
void vlasov_species_calc_implicit_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin, double dt
);

/**
 * Implicit update, phase 2: the backward-Euler step, fout = fin + dt*(implicit
 * RHS). Must run after phase 1 for all species (cross collisions read the
 * partners' moments).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fin Input distribution function (NULL if no kinetic aspect)
 * @param fout On output, the updated distribution
 * @param dt Time-step
 */
void vlasov_species_rhs_implicit(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  struct gkyl_array *fout, double dt
);

/**
 * Implicit update, phase 3: apply BCs and copy the update into the stage
 * state.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fout Updated distribution (NULL if no kinetic aspect)
 * @param fstage Stage distribution the update replaces (NULL if no kinetic aspect)
 * @param tcurr Current time
 */
void vlasov_species_finish_implicit_update(
  gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fout,
  struct gkyl_array *fstage, double tcurr
);

/**
 * Project initial conditions for each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param t0 Time for initial conditions
 */
void vlasov_species_apply_ic(gkyl_vlasov_app *app, struct vlasov_species *sp, double t0);

/**
 * Compute integrated diagnostic moments for each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param tm Time at which the diagnostic is computed
 */
void vlasov_species_calc_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm);

/**
 * Compute the integrated L2 norm of the distribution (kinetic aspect).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param tm Time at which the diagnostic is computed
 */
void vlasov_species_calc_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm);

/**
 * Write the evolved state and attendant diagnostics (sources, emission spectra,
 * fluid primitive variables) of each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param tm Time-stamp
 * @param frame Frame number
 */
void vlasov_species_write(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame);

/**
 * Write diagnostic moments (kinetic aspect).
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param tm Time-stamp
 * @param frame Frame number
 */
void vlasov_species_write_mom(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame);

/**
 * Append integrated diagnostic moments for each present aspect.
 *
 * @param app Vlasov app object
 * @param sp Species container
 */
void vlasov_species_write_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *sp);

/**
 * Append the integrated L2 norm of the distribution (kinetic aspect).
 *
 * @param app Vlasov app object
 * @param sp Species container
 */
void vlasov_species_write_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp);

/**
 * Append the LTE-correction iteration status (kinetic aspect).
 *
 * @param app Vlasov app object
 * @param sp Species container
 */
void vlasov_species_write_lte_corr_status(gkyl_vlasov_app *app, struct vlasov_species *sp);

/**
 * Read each present aspect's evolved state from the named file and rebuild what
 * a restart does not carry (Jacobian rescales, boundary fluxes, BCs, sources,
 * applied accelerations at the time of the file), through the aspect's own method.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param fname File to read
 * @return Restart status
 */
struct gkyl_app_restart_status vlasov_species_from_file(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const char *fname
);

/**
 * Restart a species from a frame: project the initial conditions of each
 * present aspect (fills the fixed-function BC buffers), read each aspect's
 * frame file, mark its diagnostics to append, and set the app clock to the
 * frame time.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @param frame Frame number
 * @return Restart status
 */
struct gkyl_app_restart_status vlasov_species_read_from_frame(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int frame
);

/**
 * Initialize the kinetic aspect of every kinetic species. Species are
 * looked up by name during initialization, so every container must already
 * be constructed.
 *
 * @param vm App inputs
 * @param app Vlasov app object
 */
void vlasov_species_init_kinetic_aspects(struct gkyl_vm *vm, gkyl_vlasov_app *app);

/**
 * Link the kinetic aspects to each other: emission walls to their impact
 * species, cross-collision partners, and sources (which may create boundary-
 * flux updaters on their source species). Runs after all kinetic aspects are
 * initialized.
 *
 * @param app Vlasov app object
 */
void vlasov_species_link_kinetic_aspects(gkyl_vlasov_app *app);

/**
 * Initialize the fluid aspect of every fluid species, then its source.
 *
 * @param vm App inputs
 * @param app Vlasov app object
 */
void vlasov_species_init_fluid_aspects(struct gkyl_vm *vm, gkyl_vlasov_app *app);

/**
 * Write the record of the closure solves of a PKPM species to the statistics
 * (the nodes where the closures were nulled and the smallest conditioning proxy,
 * reduced over ranks).
 *
 * @param app Vlasov app object
 * @param sp Species container (a PKPM species)
 * @param fp File of the statistics (NULL on ranks other than 0)
 */
void vlasov_pkpm_species_stat_write(gkyl_vlasov_app *app, struct vlasov_species *sp, FILE *fp);

/**
 * Initialize both aspects of every PKPM species, then what refers to other
 * species: the sources (the boundary-flux rate) and the cross-species collision
 * terms.
 *
 * @param vm App inputs
 * @param app Vlasov app object
 */
void vlasov_species_init_pkpm(struct gkyl_vm *vm, gkyl_vlasov_app *app);

/**
 * Number of particles per unit time the species loses through its walls, from the
 * boundary fluxes of the stage. For the boundary-flux source of a species of any type.
 *
 * @param app Vlasov app object
 * @param sp Species container
 * @return Loss rate summed over the walls and the ranks
 */
double vlasov_species_wall_loss_rate(gkyl_vlasov_app *app, const struct vlasov_species *sp);

/**
 * Release the aspects of the species and the container's allocations.
 *
 * @param app Vlasov app object
 * @param sp Species container
 */
void vlasov_species_release(const gkyl_vlasov_app *app, struct vlasov_species *sp);

/**
 * Gather the current evolved state of every species into arrays indexed by
 * species (NULL where a species lacks that aspect).
 *
 * @param app Vlasov app object
 * @param distf On output, distf[i] is the distribution of species i (may be NULL to skip)
 * @param fluid On output, fluid[i] is the fluid moments of species i (may be NULL to skip)
 */
void vlasov_species_gather_state(
  gkyl_vlasov_app *app, struct gkyl_array *distf[], struct gkyl_array *fluid[]
);

/**
 * Gather the distribution function of every species as read-only inputs for the
 * field initial conditions, energy diagnostic and output.
 *
 * @param app Vlasov app object
 * @param fin On output, fin[i] is the distribution of species i (NULL if it has no kinetic aspect)
 */
void vlasov_species_gather_dist(gkyl_vlasov_app *app, const struct gkyl_array *fin[]);

/** vm_species / vm_fluid_species private API (vm_species*.c, vm_fluid_species*.c). */

/** vm_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param sm Species moment object
 * @param mom_type Type of moment to compute.
 * @param is_integrated Whether to compute volume-integrated moment.
 */
void vm_species_moment_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_species_moment *sm,
  enum gkyl_distribution_moments mom_type, bool is_integrated
);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param sm vm species moment
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void vm_species_moment_calc(
  const struct vm_species_moment *sm, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin
);

/**
 * Divide the configuration-space Jacobian J_x of the position map out of a moment of the stored
 * J_x*J_v*f, giving the physical moment (no-op on an identity map). Such moments carry J_x in every
 * component, except the LTE moments (n, V_drift, T/m), where only the density does.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param sm Species moment object @a mom holds a moment of
 * @param mom Moment array (sm->marr or sm->marr_host), divided in place over the local range
 */
void vm_species_moment_divide_Jx(
  const struct gkyl_vlasov_app *app, const struct vm_species *s, const struct vm_species_moment *sm,
  struct gkyl_array *mom
);

/**
 * Release species moment object.
 *
 * @param app Vlasov app object
 * @param sm Species moment object to release
 */
void vm_species_moment_release(
  const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm
);

/** vm_species_projection API */

/**
 * Project a configuration-space function (e.g. a collision frequency) onto the
 * configuration-space basis, sampling it at physical coordinates so that it is
 * consistent with initial conditions and sources when a position map is used.
 *
 * @param app vlasov app object
 * @param vms Species object (provides the position map)
 * @param func Function to project
 * @param ctx Context for func
 * @param fout_ho Output (host) array
 */
void vm_species_proj_conf_func(
  struct gkyl_vlasov_app *app, struct vm_species *vms, evalf_t func, void *ctx,
  struct gkyl_array *fout_ho
);

// Reference parameters for Spitzer collision frequencies and Morse's alpha_E.
struct vm_coll_ref_params {
  double den, temp; // Reference density and temperature.
  double hbar, eps0; // Planck's constant/2 pi, vacuum permittivity.
};

/**
 * Reference parameters of a collisions input block, with unset (zero) entries
 * defaulting to 1 (normalized units), matching the Lua interface.
 *
 * @param coll Collisions input block
 * @return Reference parameters
 */
struct vm_coll_ref_params vm_species_coll_ref_params(const struct gkyl_vlasov_collisions *coll);

/**
 * Initialize species projection object.
 *
 * @param app vlasov app object
 * @param s Species object
 * @param inp Input struct for projection (contains functions pointers for type of projection)
 * @param proj Species projection object
 */
void vm_species_projection_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct gkyl_vlasov_projection inp,
  struct vm_proj *proj
);

/**
 * Compute species projection
 *
 * @param app vlasov app object
 * @param species Species object
 * @param proj Species projection object
 * @param f Output distribution function from projection
 * @param tm Time for use in projection
 */
void vm_species_projection_calc(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_proj *proj,
  struct gkyl_array *f, double tm
);

/**
 * Release species projection object.
 *
 * @param app vlasov app object
 * @param proj Species projection object to release
 */
void vm_species_projection_release(const struct gkyl_vlasov_app *app, const struct vm_proj *proj);

/** vm_species_lte API */

/**
 * Initialize species lte object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param lte Species lte object
 * @param corr_inp Input struct with moment correction inputs
 */
void vm_species_lte_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lte *lte,
  struct correct_all_moms_inp corr_inp
);

/**
 * Compute the LTE moments (n, V_drift, T/m) of fin into lte->moms, dividing the
 * configuration-space Jacobian out of the density in place so the moments are
 * physical (the stored fin is J_x J_v f).
 *
 * @param app Vlasov app object
 * @param vms Pointer to species
 * @param lte Pointer to lte object
 * @param fin Input distribution function
 */
void vm_species_lte_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, const struct vm_lte *lte,
  const struct gkyl_array *fin
);

/**
 * Compute LTE distribution from input moments
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lte Pointer to lte object
 * @param moms_lte Input LTE moments
 */
void vm_species_lte_from_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lte *lte,
  const struct gkyl_array *moms_lte
);

/**
 * Compute equivalent LTE distribution from input distribution function.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lte Pointer to lte
 * @param fin Input distribution function
 */
void vm_species_lte(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lte *lte,
  const struct gkyl_array *fin
);

/**
 * Write the LTE correction status.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 */
void vm_species_lte_write_max_corr_status(gkyl_vlasov_app *app, struct vm_species *vms);

/**
 * Release species lte object.
 *
 * @param app Vlasov app object
 * @param lte Species lte object to release
 */
void vm_species_lte_release(const struct gkyl_vlasov_app *app, const struct vm_lte *lte);

/** vm_species_collisionless API */

/**
 * Whether a kinetic species couples to a GR-Maxwell field through the triad
 * Lorentz force and current deposition: requires a GR-Maxwell field, a triad
 * model, and a geometry object whose GR coupling is enabled.
 *
 * @param app Vlasov app object
 * @param field_id Field type
 * @param model_id Species model type
 * @return True if the GR triad-Maxwell coupling is active for this species
 */
bool vm_species_has_gr_em_triad_coupling(
  const struct gkyl_vlasov_app *app, enum gkyl_field_id field_id, enum gkyl_model_id model_id
);

/**
 * Initialize species collisionless object.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param cls Species collisionless object.
 */
void vm_species_collisionless_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_collisionless *cls
);

/**
 * Compute species applied acceleration term
 *
 * @param app Vlasov app object
 * @param cls Species collisionless object
 * @param tm Time for use in acceleration
 */
void vm_species_collisionless_app_accel(
  gkyl_vlasov_app *app, struct vm_collisionless *cls, double tm
);

/**
 * Compute RHS contribution from collisionless terms.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param cls Species collisionless object.
 * @param fin Input distribution function.
 * @param em Input electromagnetic fields.
 * @param rhs Collisionless contribution to df/dt.
 */
void vm_species_collisionless_rhs(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_collisionless *cls,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs
);

/**
 * Release species collisionless object.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param cls Species collisionless object.
 */
void vm_species_collisionless_release(
  const struct gkyl_vlasov_app *app, const struct vm_species *vms,
  const struct vm_collisionless *cls
);

/** vm_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param lbo Species LBO object
 */
void vm_species_lbo_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo
);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param lbo Species LBO object
 */
void vm_species_lbo_cross_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo
);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vm_species_lbo_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vm_species_lbo_cross_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

/**
 * Compute RHS from LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void vm_species_lbo_rhs(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

/**
 * Write moments from LBO object.
 *
 * @param app Vlasov app object.
 * @param species Pointer to species.
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void vm_species_lbo_write_mom(
  gkyl_vlasov_app *app, struct vm_species *species, double tm, int frame
);

/**
 * Release species LBO object.
 *
 * @param app Vlasov app object
 * @param lbo Species LBO object to release
 */
void vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo);

/** vm_species_bgk API */

/**
 * Initialize species BGK collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param bgk Species BGK object
 */
void vm_species_bgk_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_bgk_collisions *bgk
);

/**
 * Initialize species BGK cross-collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param bgk Species BGK object
 */
void vm_species_bgk_cross_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_bgk_collisions *bgk
);

/**
 * Compute necessary moments for BGK collisions.
 *
 * @param app Vlasov app object
 * @param vms Pointer to species.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 */
void vm_species_bgk_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin
);

/**
 * Compute the moments needed by implicit BGK collisions (no-op without them).
 *
 * @param app Vlasov app object
 * @param vms Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void vm_species_bgk_moms_implicit(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin
);

/**
 * Compute and store a fixed temperature for BGK collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void vm_species_bgk_moms_fixed_temp(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin
);

/**
 * Compute RHS from BGK collisions.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 * @param rhs On output, the RHS from bgk
 */
void vm_species_bgk_rhs(
  gkyl_vlasov_app *app, struct vm_species *species, struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

/**
 * Build the target of implicit BGK collisions without applying the operator (no-op without
 * them): the sum over self and cross collisions of the collision frequency times its
 * Maxwellian in bgk->nu_f_lte and the sum of the frequencies in bgk->nu_sum. For an implicit
 * update that takes the collisions together with another operator of the same form.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void vm_species_bgk_target_implicit(
  gkyl_vlasov_app *app, struct vm_species *species, struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin
);

/**
 * Compute RHS from BGK collisions (implicit integrator).
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 * @param dt Time step.
 * @param rhs On output, the implicit RHS from bgk
 */
void vm_species_bgk_rhs_implicit(
  gkyl_vlasov_app *app, struct vm_species *species, struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin, double dt, struct gkyl_array *rhs
);

/**
 * Write moments from BGK object.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void vm_species_bgk_write_mom(
  gkyl_vlasov_app *app, struct vm_species *species, double tm, int frame
);

/**
 * Release species BGK object.
 *
 * @param app Vlasov app object
 * @param bgk Species BGK object to release
 */
void vm_species_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_bgk_collisions *bgk);

/** vm_species_source_bgk API */

/**
 * Initialize the temperature source nu_Q (f_M(n, u, T_Q) - f) of a species: project the
 * rate and the temperature, allocate the target and choose the methods of the explicit or
 * the implicit update (all empty without a source and for a static species). The collision
 * objects of the species must be initialized first.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param src Temperature source object
 */
void vm_species_source_bgk_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source_bgk *src
);

/**
 * Compute the moments of the target of the explicit temperature source: the LTE moments of
 * the input distribution with the temperature replaced by T_Q.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to temperature source
 * @param fin Input distribution function
 */
void vm_species_source_bgk_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_source_bgk *src,
  const struct gkyl_array *fin
);

/**
 * Accumulate the explicit temperature source onto the RHS of the species and its rate onto
 * the CFL frequency. Run after vm_species_source_bgk_moms.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to temperature source
 * @param fin Input distribution function
 * @param rhs RHS of the species to accumulate onto
 */
void vm_species_source_bgk_rhs(
  gkyl_vlasov_app *app, struct vm_species *species, struct vm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

/**
 * Compute the moments of the implicit temperature source over dt (no-op without one):
 * those of its target, and the temperature it gives the species over the step,
 * T' = (T + dt nu_Q T_Q)/(1 + dt nu_Q), written into the LTE moments of the species, where
 * implicit BGK collisions (of the species and of its partners) read their targets. Run after
 * vm_species_bgk_moms_implicit.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to temperature source
 * @param fin Input distribution function
 * @param dt Time-step of the backward-Euler step
 */
void vm_species_source_bgk_moms_implicit(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt
);

/**
 * Accumulate the implicit RHS of the temperature source, together with the implicit BGK
 * collisions of the species if it has them, onto out (no-op without an implicit source):
 * (nu_c f_M(n, u, T') + nu_Q f_M(n, u, T_Q) - (nu_c + nu_Q) fin)/(1 + dt (nu_c + nu_Q)),
 * so that fin + dt out is the backward-Euler step of both. Run after the moments of every
 * species.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to temperature source
 * @param fin Input distribution function
 * @param dt Time-step of the backward-Euler step
 * @param out Array to accumulate the implicit RHS onto
 */
void vm_species_source_bgk_implicit_advance(
  gkyl_vlasov_app *app, struct vm_species *species, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt, struct gkyl_array *out
);

/**
 * Release the temperature source object.
 *
 * @param app Vlasov app object
 * @param src Temperature source object to release
 */
void vm_species_source_bgk_release(
  const struct gkyl_vlasov_app *app, const struct vm_source_bgk *src
);

/** vm_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param bflux Species boundary flux object
 */
void vm_species_bflux_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_boundary_fluxes *bflux
);

/**
 * Compute the boundary fluxes of the input distribution and their integrated
 * moments (used by boundary-flux sources).
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bflux Species boundary flux object
 * @param fin Input distribution function
 * @param rhs Work array whose ghost cells receive the boundary fluxes (overwritten by the BCs)
 */
void vm_species_bflux_rhs(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

/**
 * Number of particles per unit time through the boundaries, from the integrated
 * moments of the boundary fluxes vm_species_bflux_rhs last computed.
 *
 * @param app Vlasov app object
 * @param bflux Species boundary flux object
 * @return Loss rate summed over the boundaries and the ranks
 */
double vm_species_bflux_loss_rate(gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux);

/**
 * Release species boundary flux object.
 *
 * @param app Vlasov app object
 * @param bflux Species boundary flux object to release
 */
void vm_species_bflux_release(
  const struct gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux
);

/** vm_species_emission API */

/**
 * Initialize species emission object
 *
 * @param app Vlasov app object
 * @param emit Species emission object
 * @param dir Direction of BC
 * @param edge Edge of configuration space
 * @param ctx Emission context
 */
void vm_species_emission_init(
  struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit, int dir, enum gkyl_edge_loc edge,
  void *ctx
);

/**
 * Initialize emission BC cross-species object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param emit Species emission object
 */
void vm_species_emission_cross_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_emitting_wall *emit
);

/**
 * Apply emission BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param emit Pointer to emission object
 * @param fout Field to apply BCs
 * @param tcurr Current time
 */
void vm_species_emission_apply_bc(
  struct gkyl_vlasov_app *app, const struct vm_species *s, const struct vm_emitting_wall *emit,
  struct gkyl_array *fout, double tcurr
);

/**
 * Write emission spectrum distribution function
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param emit Pointer to emission object
 * @param mt Write meta
 * @param frame Current frame
 */
void vm_species_emission_write(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_emitting_wall *emit,
  struct gkyl_msgpack_data *mt, int frame
);

/**
 * Release species emission object.
 *
 * @param emit Species emission object to release
 */
void vm_species_emission_release(const struct vm_emitting_wall *emit);

/** vm_species_source API */

/**
 * Initialize species source object.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param src Species source object
 */
void vm_species_source_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src
);

/**
 * Compute species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 * @param tm Time for use in source
 */
void vm_species_source_calc(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_source *src, double tm
);

/**
 * Compute the moments needed to re-scale the density of adaptive sources.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 * @param fin Input distribution function
 */
void vm_species_source_adapt_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_source *src,
  const struct gkyl_array *fin
);

/**
 * Re-scale the adaptive source from the moments computed by
 * vm_species_source_adapt_moms.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 */
void vm_species_source_adapt(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_source *src
);

/**
 * Accumulate the source onto the RHS of the species (a boundary-flux source is
 * scaled by the boundary fluxes of its source species first).
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param rhs RHS of the species to accumulate onto
 */
void vm_species_source_rhs(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_source *src,
  struct gkyl_array *rhs
);

/**
 * Write source diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void vm_species_source_write(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source *src, double tm, int frame
);

/**
 * Write source moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void vm_species_source_write_mom(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source *src, double tm, int frame
);

/**
 * Calc source integrated moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 * @param tm Time for source diagnostic.
 */
void vm_species_source_calc_integrated_mom(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source *src, double tm
);

/**
 * Write source integrated moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 */
void vm_species_source_write_integrated_mom(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source *src
);

/**
 * Release species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src);

/** vm_species API */

/**
 * Initialize species.
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param s On output, initialized species object
 */
void vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0);

/**
 * Compute RHS from species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param em EM field
 * @param rhs On output, the RHS from the species object
 * @return Maximum stable time-step
 */
double vm_species_rhs(
  gkyl_vlasov_app *app, struct vm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *em, struct gkyl_array *rhs
);

/**
 * Take the implicit step of the species (implicit BGK collisions and/or
 * temperature source): one backward-Euler step over dt.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param rhs On output, the updated distribution fin + dt*(implicit RHS)
 * @param dt Time-step of the backward-Euler step
 * @return DBL_MAX (the implicit step does not limit the time-step)
 */
double vm_species_rhs_implicit(
  gkyl_vlasov_app *app, struct vm_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, double dt
);

/**
 * Compute the moments of the implicit update of the species (BGK collisions and
 * the temperature source). Run for every species before any
 * vm_species_rhs_implicit, which reads the partners' moments.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param dt Time-step of the implicit update
 */
void vm_species_calc_implicit_moms(
  gkyl_vlasov_app *app, struct vm_species *species, const struct gkyl_array *fin, double dt
);

/**
 * Compute the time-dependent applied acceleration (no-op when it does not
 * evolve in time).
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param tcurr Current time
 */
void vm_species_calc_app_accel(gkyl_vlasov_app *app, struct vm_species *species, double tcurr);

/**
 * Compute the moments that rescale adaptive sources (no-op for other sources).
 * Run for every species before any source is accumulated: an adaptive source
 * reads its partners' moments.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 */
void vm_species_calc_source_moms(
  gkyl_vlasov_app *app, struct vm_species *species, const struct gkyl_array *fin
);

/**
 * Accumulate the source of the species onto its RHS (no-op without a source):
 * project a time-dependent source, rescale an adaptive one, and accumulate (a
 * static species only projects a time-dependent source). Run after
 * vm_species_calc_source_moms for every species.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param tcurr Current time
 * @param rhs RHS of the species to accumulate onto
 */
void vm_species_rhs_source(
  gkyl_vlasov_app *app, struct vm_species *species, double tcurr, struct gkyl_array *rhs
);

/**
 * Scale and accumulate for the forward Euler method: out = a*out + inp.
 *
 * @param vms Pointer to species.
 * @param out Output array.
 * @param a Scale factor (the time-step).
 * @param inp Input array.
 */
void vm_species_step_f(
  struct vm_species *vms, struct gkyl_array *out, double a, const struct gkyl_array *inp
);

/**
 * Combine for rk3 method.
 *
 * @param vms Pointer to species.
 * @param out Output array.
 * @param c1 Scaling factor.
 * @param arr1 Input array.
 * @param c2 Scaling factor.
 * @param arr2 Input array.
 * @param rng Range.
 */
void vm_species_combine(
  struct vm_species *vms, struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range *rng
);

/**
 * Copy for rk3 method.
 *
 * @param vms Pointer to species.
 * @param out Output array.
 * @param inp Input array.
 * @param range Range.
 */
void vm_species_copy_range(
  struct vm_species *vms, struct gkyl_array *out, const struct gkyl_array *inp,
  const struct gkyl_range *range
);

/**
 * Apply BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param f Field to apply BCs
 * @param tcurr Current time
 */
void vm_species_apply_bc(
  gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f, double tcurr
);

/**
 * Species write function: the distribution function and, when requested, the
 * source and the emission spectra of the frame.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void vm_species_write(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);

/**
 * Species moment write function: the diagnostic moments of the distribution
 * and, when requested, of the source.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void vm_species_write_mom(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame);

/**
 * Species calc integrated moment function (of the distribution and, when
 * requested, of the source).
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm Time at which diagnostic is computed.
 */
void vm_species_calc_integrated_mom(gkyl_vlasov_app *app, struct vm_species *vms, double tm);

/**
 * Species write integrated moment function (of the distribution and, when
 * requested, of the source).
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 */
void vm_species_write_integrated_mom(gkyl_vlasov_app *app, struct vm_species *vms);

/**
 * Read the distribution function from a file and rebuild what the file does
 * not carry: the velocity- and configuration-space Jacobian rescale, the
 * boundary fluxes, the BCs, and the source and the applied acceleration at the
 * time of the file.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param fname File to read.
 * @return Restart status.
 */
struct gkyl_app_restart_status vm_species_from_file(
  gkyl_vlasov_app *app, struct vm_species *vms, const char *fname
);

/**
 * Restart the species from its file of a frame (see vm_species_from_file) and
 * mark its diagnostics to append to their existing files.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param frame Frame number.
 * @return Restart status.
 */
struct gkyl_app_restart_status vm_species_read_from_frame(
  gkyl_vlasov_app *app, struct vm_species *vms, int frame
);

/**
 * Compute L2 energy (f^2) of the distribution function diagnostic
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm Time at which diagnostic is computed.
 */
void vm_species_calc_L2(gkyl_vlasov_app *app, struct vm_species *vms, double tm);

/**
 * Species write integrated L^2 energy (f^2) function.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 */
void vm_species_write_L2(gkyl_vlasov_app *app, struct vm_species *vms);

/**
 * Fill stat object in app with total number of iterations
 * used to correct moments in LTE projection object.
 * Also fills stat object with number of times correction object called.
 *
 * @param app App object to update stat timers.
 */
void vm_species_n_iter_corr(gkyl_vlasov_app *app);

/**
 * Delete resources used in species.
 *
 * @param app Vlasov app object
 * @param s Species object to delete
 */
void vm_species_release(const gkyl_vlasov_app *app, const struct vm_species *s);

/** The vm_field/vp_field API is declared in gkyl_vlasov_field_priv.h. */

/** vm_fluid_species_source API */

/**
 * Initialize fluid species source object.
 *
 * @param app Vlasov app object
 * @param fluid_species Species object
 * @param src Species source object
 */
void vm_fluid_species_source_init(
  struct gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct vm_fluid_source *src
);

/**
 * Compute fluid species applied source term
 *
 * @param app Vlasov app object
 * @param fluid_species Species object
 * @param tm Time for use in source
 */
void vm_fluid_species_source_calc(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm
);

/**
 * Accumulate the fluid source onto the fluid RHS.
 *
 * @param app Vlasov app object
 * @param species Pointer to fluid species
 * @param src Pointer to source
 * @param rhs Fluid RHS of the species to accumulate onto
 */
void vm_fluid_species_source_rhs(
  gkyl_vlasov_app *app, const struct vm_fluid_species *species, struct vm_fluid_source *src,
  struct gkyl_array *rhs
);

/**
 * Release fluid species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vm_fluid_species_source_release(
  const struct gkyl_vlasov_app *app, const struct vm_fluid_source *src
);

/** vm_fluid_species API */

/**
 * Create new fluid species object
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param f On output, initialized fluid species object
 */
void vm_fluid_species_init(
  struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f
);

/**
 * Compute fluid species initial conditions.
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param t0 Time for use in ICs
 */
void vm_fluid_species_apply_ic(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0
);

/**
 * Compute fluid species applied acceleration term
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param tm Time for use in acceleration
 */
void vm_fluid_species_calc_app_accel(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm
);

/**
 * Recompute a time-dependent applied acceleration at a stage time, when the
 * RHS of the fluid species applies it (no field, or a field coupled
 * explicitly). Otherwise the implicit fluid-EM coupling applies the
 * acceleration and evaluates it at its own time.
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param tm Stage time
 */
void vm_fluid_species_calc_stage_app_accel(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm
);

/**
 * Apply the vertex positivity limiter to a fluid state over the local range
 * (Euler and isothermal Euler at polynomial order 1; a no-op otherwise). The
 * boundary conditions are not reapplied.
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param fluid In/out fluid state
 */
void vm_fluid_species_apply_positivity(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct gkyl_array *fluid
);

/**
 * Accumulate the current of a fluid species coupled explicitly to the Maxwell
 * field onto the RHS of the electric field: -q/(m eps0) times its momentum.
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param fluid Input array fluid species
 * @param emout On output, the RHS of the field with the current accumulated
 */
void vm_fluid_species_accumulate_current(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, const struct gkyl_array *fluid,
  struct gkyl_array *emout
);

/**
 * Compute primitive variables (bulk velocity, u, and pressure, p, if pressure present)
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object (where primitive variables are stored)
 * @param fluid Input array fluid species
 */
void vm_fluid_species_prim_vars(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, const struct gkyl_array *fluid
);

/**
 * Limit the slopes of the fluid state, apply the vertex positivity limiter and
 * reapply the BCs (Euler and isothermal Euler; no-op for other equations).
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param fluid Fluid state, limited in place
 */
void vm_fluid_species_limiter(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct gkyl_array *fluid
);

/**
 * Compute RHS from fluid species equations
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param fluid Input fluid species
 * @param em EM field
 * @param rhs On output, the RHS from the fluid species solver
 * @return Maximum stable time-step
 */
double vm_fluid_species_rhs(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, const struct gkyl_array *fluid,
  const struct gkyl_array *em, struct gkyl_array *rhs
);

/**
 * Accumulate the source of the fluid species onto its RHS (no-op without a
 * source). Run after the RHS of every species.
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param rhs Fluid RHS of the species to accumulate onto
 */
void vm_fluid_species_rhs_source(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct gkyl_array *rhs
);

/**
 * Forward-Euler accumulate for the fluid state: out = dt*out + inp.
 */
void vm_fluid_species_step_f(
  struct vm_fluid_species *fluid_species, struct gkyl_array *out, double dt,
  const struct gkyl_array *inp
);

/**
 * Combine fluid RK stages: out = c1*arr1 + c2*arr2 over rng.
 */
void vm_fluid_species_combine(
  struct vm_fluid_species *fluid_species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng
);

/**
 * Copy the fluid state: out = inp over range.
 */
void vm_fluid_species_copy_range(
  struct vm_fluid_species *fluid_species, struct gkyl_array *out, const struct gkyl_array *inp,
  const struct gkyl_range *range
);

/**
 * Apply BCs to fluid species
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param f Fluid Species to apply BCs
 */
void vm_fluid_species_apply_bc(
  gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f
);

/**
 * Computed the integrated quantities for the fluid system.
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param tm Time integrated quantities are being computed at.
 */
void vm_fluid_species_calc_integrated_mom(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm
);

/**
 * Write out the evolved fluid species and other potential primitive/auxiliary variables.
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param tm Time fluid quantities are being written at.
 * @param frame Frame number for I/O.
 */
void vm_fluid_species_write(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm, int frame
);

/**
 * Write fluid species integrated moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param fluid_species Pointer to fluid species.
 */
void vm_fluid_species_write_integrated_mom(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species
);

/**
 * Read the fluid state from a file and rebuild what the file does not carry:
 * the BCs, and the source and the applied acceleration at the time of the file.
 *
 * @param app Vlasov app object.
 * @param fluid_species Pointer to fluid species.
 * @param fname File to read.
 * @return Restart status.
 */
struct gkyl_app_restart_status vm_fluid_species_from_file(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, const char *fname
);

/**
 * Restart the fluid species from its file of a frame (see
 * vm_fluid_species_from_file) and mark its diagnostics to append to their
 * existing files.
 *
 * @param app Vlasov app object.
 * @param fluid_species Pointer to fluid species.
 * @param frame Frame number.
 * @return Restart status.
 */
struct gkyl_app_restart_status vm_fluid_species_read_from_frame(
  gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, int frame
);

/**
 * Release resources allocated by fluid species
 *
 * @param app Vlasov app object
 * @param f Fluid_Species object to release
 */
void vm_fluid_species_release(const gkyl_vlasov_app *app, struct vm_fluid_species *f);

// ---- PKPM species route (vm_pkpm_species*.c) ----

/** vm_pkpm_species_moment API */

/**
 * Compute the diagnostic moments [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
 * into s->kinetic->pkpm.moms_diag: the moments of the normalized distributions times n xi^k,
 * k = [0, 1, 2, 2, 3, 3, 4, 4]. Uses s->fluid->pkpm.prim_io as scratch; s->fluid->pkpm.xi must be
 * current.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param phase_rng Phase-space range
 * @param conf_rng Configuration-space range
 * @param fin Input distribution functions
 * @param fluid Input fluid state (the density gives n)
 */
void vm_pkpm_species_moment_diag_calc(
  struct gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin, const struct gkyl_array *fluid
);

/**
 * Initialize species moment object with a given mass factor.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param sm Species moment object
 * @param is_diag True for the diagnostic moment set, false for the moments of the update
 * @param mass Mass factor of the moments (1.0 for the moments of the normalized distributions)
 */
void vm_pkpm_species_moment_init_mass(
  struct gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_species_moment *sm,
  bool is_diag, double mass
);

/**
 * Initialize species moment object with the species mass.
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param sm Species moment object
 * @param is_diag True for the diagnostic moment set, false for the moments of the update
 */
void vm_pkpm_species_moment_init(
  struct gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_species_moment *sm,
  bool is_diag
);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void vm_pkpm_species_moment_calc(
  const struct vm_pkpm_species_moment *sm, const struct gkyl_range phase_rng,
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
void vm_pkpm_species_moment_calc_into(
  const struct vm_pkpm_species_moment *sm, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin, struct gkyl_array *mout
);

/**
 * Release species moment object.
 *
 * @param app Vlasov app object
 * @param sm Species moment object to release
 */
void vm_pkpm_species_moment_release(
  const struct gkyl_vlasov_app *app, const struct vm_pkpm_species_moment *sm
);

/** vm_pkpm_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void vm_pkpm_species_lbo_init(
  struct gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_lbo_collisions *lbo
);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vm_pkpm_species_lbo_moms(
  gkyl_vlasov_app *app, const struct vlasov_species *species, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
);

/**
 * Compute RHS from LBO collisions
 *
 * @param app Vlasov app object
 * @param kin Kinetic aspect of the species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void vm_pkpm_species_lbo_rhs(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);

/**
 * Release species LBO object.
 *
 * @param app Vlasov app object
 * @param sm Species LBO object to release
 */
void vm_pkpm_species_lbo_release(
  const struct gkyl_vlasov_app *app, const struct vm_pkpm_lbo_collisions *lbo
);

/** vm_pkpm_species API */

/**
 * Initialize species.
 *
 * @param pkpm Input PKPM data
 * @param app Vlasov app object
 * @param s On output, initialized species object
 */
void vm_pkpm_species_init(
  struct gkyl_vm *pkpm, struct gkyl_vlasov_app *app, struct vlasov_species *s
);

/**
 * Convert the distributions of the stage to the configuration-nodal, velocity-modal form and
 * clear the increments of the nodal operators (species with such operators only).
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution functions
 */
void vm_pkpm_species_nodal_dist(
  gkyl_vlasov_app *app, const struct vlasov_species *species, const struct gkyl_array *fin
);

/**
 * Accumulate the increments of the nodal operators onto modal distributions (species with
 * such operators only).
 *
 * @param app Vlasov app object
 * @param kin Kinetic aspect of the species
 * @param out Modal distributions to increment
 */
void vm_pkpm_species_nodal_accumulate(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct gkyl_array *out
);

/**
 * Compute the configuration-space variables of a stage over the extended range: the raw
 * moments of the normalized distributions, the flow u and the thermal speed xi (nodal
 * quotients of the fluid state), the kinetic coupling moments
 * [p_par, p_perp, q_par, q_perp, s (3), r (3)], the constraints and theta = T_perp/T.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param fin Input distribution functions
 * @param fluidin Input fluid state [rho, M (3), E]
 */
void vm_pkpm_species_calc_pkpm_vars(
  gkyl_vlasov_app *app, struct vlasov_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Configuration-space face fluxes of the stage, the kinetic and the fluid flux and both rates in
 * one pass (the species' rates must be cleared by the caller), then the frame gradients of the
 * volume terms, which need the blend of the faces.
 */
void vm_pkpm_species_face_flux(
  gkyl_vlasov_app *app, struct vlasov_species *species, const struct gkyl_array *fin
);

/** Particle loss rates through the absorbing walls from the stored fluxes (into the bflux of the kinetic aspect). */
void vm_pkpm_species_calc_bflux(gkyl_vlasov_app *app, struct vlasov_species *s);

// Source object: see the source member of gkyl_vlasov_kinetic_species in gkyl_vlasov.h.
void vm_pkpm_species_source_init(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src
);
void vm_pkpm_species_source_calc_rate(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src
);
void vm_pkpm_species_source_calc_dist(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src
);
void vm_pkpm_species_source_kinetic_rhs(
  gkyl_vlasov_app *app, struct vm_species *kin, struct vm_pkpm_source *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);
void vm_pkpm_species_source_write(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src, double tm, int frame
);
void vm_pkpm_species_source_release(const gkyl_vlasov_app *app, const struct vm_pkpm_source *src);

// Project a configuration-space profile at the physical coordinates of the position map and
// scale it. Returns a new host array on the extended range: the values at the 2^cdim interior
// Gauss-Legendre nodes if is_nodal is set, the p=1 expansion otherwise.
struct gkyl_array *vm_pkpm_species_coll_profile_new(
  gkyl_vlasov_app *app, void (*eval)(double t, const double *xn, double *fout, void *ctx),
  void *ctx, double scale, bool is_nodal
);
// The self-collision frequency of the input at the nodes, scaled by nu_frac: a new host array on
// the extended range. Without an input function it is zero and the species' collision-frequency
// object (vm_pkpm_species_coll_cross.c) fills it every stage.
struct gkyl_array *vm_pkpm_species_coll_self_nu_new(
  gkyl_vlasov_app *app, const struct vlasov_species *s, double nu_frac
);
// The target of the BGK-form operators: a new host array with the unit Maxwellian in w on the
// velocity grid of the species (3 modes per cell), its discrete moments corrected to [1, 0, 1].
struct gkyl_array *vm_pkpm_species_bgk_target_new(gkyl_vlasov_app *app, struct vlasov_species *s);

// BGK collisions: nu (f_target - f) with the normalized unit Maxwellian target. Self collisions
// conserve the density, momentum and energy, so they never touch the fluid.
void vm_pkpm_species_bgk_init(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_bgk_collisions *bgk
);

// Cross-species collisions of the species with a fluid aspect (vm_coll_exchange.c).
// Operations of a species that needs none.
void vm_coll_cross_prim_none(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
);
void vm_coll_cross_moms_none(gkyl_vlasov_app *app, struct vlasov_species *s);
// Node values [n, u (3), vtsq] of a fluid state of the species.
void vm_coll_cross_prim_nodal(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
);
// The two frequencies of the pair (species, partner i) at the nodes, into nu_self and nu_other of
// the species. The partner forms the same two from the same node values.
void vm_coll_cross_pair_nu(gkyl_vlasov_app *app, struct vlasov_species *s, int i);
// Reference density, temperature, Planck constant and permittivity of the Coulomb logarithm of
// a computed frequency (1 where the input does not set them).
void vm_coll_cross_ref_params(const struct gkyl_vlasov_collisions *coll, double ref[4]);
// The partners of a species among the first num_partners of its input: those with a fluid
// aspect, with their reference frequencies and exchange coefficients, and the kinetic ones (of
// a PKPM species); and the node values of its fluid state. Call once the containers of all
// species exist and the kinetic species are linked.
void vm_coll_cross_partners_init(gkyl_vlasov_app *app, struct vlasov_species *s, int num_partners);
void vm_coll_cross_partners_release(const struct vm_coll_cross *cross);
// Collision data of a fluid species; call for every fluid species. One without partners gets
// operations that do nothing.
void vm_fluid_species_coll_cross_init(gkyl_vlasov_app *app, struct vlasov_species *s);
void vm_fluid_species_coll_cross_release(
  const gkyl_vlasov_app *app, const struct vm_coll_cross *cross
);
// Initialize the exchange of the species that have partners; call after the collision data of
// every species (vm_pkpm_species_coll_cross_init, vm_fluid_species_coll_cross_init).
void vm_coll_exchange_init(gkyl_vlasov_app *app);
// The momentum and energy exchange between the species with partners over dt, on the given
// fluid states (indexed by the species of the app), which it updates in place. It first forms
// the node values of those states (and the self frequencies of implicit operators).
void vm_coll_exchange_advance(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[]);
void vm_coll_exchange_release(const gkyl_vlasov_app *app);

// Computed self-collision frequency and the kinetic part of the cross-species collisions of a
// PKPM species (vm_pkpm_species_coll_cross.c), for LBO and BGK alike.
// Initialize them for a species; call for every PKPM species, after the collision objects of all
// species exist. A species with a given self frequency and no partners gets operations that do
// nothing.
void vm_pkpm_species_coll_cross_init(gkyl_vlasov_app *app, struct vlasov_species *s);
// Node values [n, u (3), vtsq] of the species from the fluid state of the stage, and its
// self-collision frequency when it is computed.
void vm_pkpm_species_coll_cross_prim(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluidin
);
// The rate of the kinetic operator of the stage with every partner at the nodes. Call after the
// node values (and the LBO moments) of every species.
void vm_pkpm_species_coll_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *s);
// Shape factor of the implicit BGK step of a species with partners, after vm_coll_exchange_advance.
void vm_pkpm_species_coll_cross_shape(
  gkyl_vlasov_app *app, struct vlasov_species *s, double dt, struct gkyl_array *shape
);
void vm_pkpm_species_coll_cross_release(
  const gkyl_vlasov_app *app, const struct vm_coll_cross *cross
);

// Temperature source (vm_pkpm_species_source_bgk.c). The explicit form enters the right-hand
// sides before the closures: its energy relaxation is computed with the species' staging and
// added by the fluid right-hand side, the relaxation of the distributions by the kinetic one. The
// implicit form is one backward-Euler step over dt on a distribution and a fluid state of the
// species (interior only), together with its implicit BGK collisions if it has them: _moms steps
// the energy of the fluid and forms the shape factor (after vm_coll_exchange_advance), _advance
// steps the distribution.
void vm_pkpm_species_source_bgk_init(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source_bgk *src
);
void vm_pkpm_species_source_bgk_calc_fluid(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source_bgk *src,
  const struct gkyl_array *fluidin
);
void vm_pkpm_species_source_bgk_kinetic_rhs(
  gkyl_vlasov_app *app, struct vm_species *kin, struct vm_pkpm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs_f
);
void vm_pkpm_species_source_bgk_implicit_moms(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source_bgk *src, double dt,
  struct gkyl_array *fluid
);
void vm_pkpm_species_source_bgk_implicit_advance(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source_bgk *src,
  struct gkyl_array *f
);
void vm_pkpm_species_source_bgk_release(
  const gkyl_vlasov_app *app, const struct vm_pkpm_source_bgk *src
);
// Project M0_perp and N_perp of the distributions f onto the plane perpendicular to the current b
// at the Gauss-Legendre nodes of b, in place, over the extended ranges.
void vm_pkpm_species_project_transverse(
  gkyl_vlasov_app *app, struct vlasov_species *species, struct gkyl_array *f
);
// Gyrophase-mixing integrating factor of the first-harmonic blocks, after the
// forward Euler accumulate of stage s (fout = fin + h rhs; the gyration is explicit in rhs). At
// s = 1 it freezes the mixing of the step (the field and the rate at the phase nodes). At damped
// nodes it writes the first-harmonic blocks of fout as the value whose SSP-RK3 combination with the
// state at the start of the step is the integrating-factor stage, and restores int M0_perp dw = 0
// where the mixing acted.
void vm_pkpm_species_harm_exp_advance(
  gkyl_vlasov_app *app, struct vlasov_species *species, int stage, double h, struct gkyl_array *fout
);

// out += a conf * phase for every block of a phase array (the bin-op kernels act on one
// block of num_basis components); phase may hold fewer blocks than out (e.g. [F, G] only).
void vm_pkpm_species_accumulate_conf_phase(
  gkyl_vlasov_app *app, struct vm_species *kin, struct gkyl_array *out, double a,
  const struct gkyl_array *conf, const struct gkyl_array *phase
);

void vm_pkpm_species_bgk_kinetic_rhs(
  gkyl_vlasov_app *app, struct vm_species *kin, struct vm_pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
);
// One backward-Euler step of the BGK operator over dt on a distribution of the species (interior
// only): _moms forms the shape factor of the step (after vm_coll_exchange_advance), _advance
// steps the distribution.
void vm_pkpm_species_bgk_implicit_moms(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_bgk_collisions *bgk, double dt
);
void vm_pkpm_species_bgk_implicit_advance(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_bgk_collisions *bgk,
  struct gkyl_array *f
);
void vm_pkpm_species_bgk_release(
  const gkyl_vlasov_app *app, const struct vm_pkpm_bgk_collisions *bgk
);

// Save and restore the records of the positivity limiter of a fluid aspect (the per-cell counts
// and the counts since the last integrated-moment call), around an application of the limiter
// that is not part of the evolution (the synchronization of the staggered split for output).
void vm_pkpm_species_positivity_record_save(struct vm_fluid_species *fl);
void vm_pkpm_species_positivity_record_restore(struct vm_fluid_species *fl);

/**
 * Apply the vertex positivity limiter to a fluid state over the local range and record the
 * cells it acted on and the counts of the current integrated-moment interval.
 *
 * @param app Vlasov app object
 * @param fl Fluid aspect of the species
 * @param fluid Fluid state [rho, M, E] to limit in place
 */
void vm_pkpm_species_apply_positivity(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, struct gkyl_array *fluid
);

/**
 * Compute the diagnostic variables that need the configuration pass (the pressure
 * force a_P), for output. Uses the species scratch arrays.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param fin Input distribution function
 * @param fluidin Input fluid species array
 */
void vm_pkpm_species_calc_diag_vars(
  gkyl_vlasov_app *app, struct vlasov_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Compute everything of the stage that reads this species and the field alone (staging phase
 * 1): the transverse projection of the first harmonic, the moments, the primitive variables,
 * the coupling moments and T_perp/T, the distributions at the nodes, the LBO moments, the node
 * values of the fluid state, the face fluxes with their gradients and wall losses, and the
 * energy relaxation of the explicit temperature source.
 *
 * @param app Vlasov app object
 * @param sp Species object
 * @param stage SSP-RK3 stage (1, 2 or 3): the closures take the speeds of the velocity-domain
 *              edges from the closures of the stage before, those of stage 1 from the last stage
 *              of the last completed step
 * @param fin Distribution functions of the stage (the first harmonic is projected in place)
 * @param fluidin Fluid state of the stage
 */
void vm_pkpm_species_calc_self_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Compute what reads other species (staging phase 2): the rate of the kinetic operator with
 * the collision partners, the rate of the boundary-flux source and the source distributions.
 * Must run after phase 1 for all species.
 *
 * @param app Vlasov app object
 * @param sp Species object
 * @param fin Distribution functions of the stage
 * @param fluidin Fluid state of the stage
 */
void vm_pkpm_species_calc_coupled_vars(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
);

/**
 * Compute the right-hand side of the kinetic equations of a species: the configuration pass,
 * the first-harmonic terms, sources, collisions, the closures and the velocity pass.
 *
 * @param app Vlasov app object
 * @param kin Kinetic aspect of the species
 * @param fin Input distribution functions
 * @param em EM field (unused: the species reads the physical fields of the stage on the field)
 * @param rhs On output, the right-hand side of the kinetic equations
 * @return Maximum stable time-step
 */
double vm_pkpm_kinetic_species_rhs(
  gkyl_vlasov_app *app, struct vm_species *kin, const struct gkyl_array *fin,
  const struct gkyl_array *em, struct gkyl_array *rhs
);

/**
 * Compute the right-hand side of the fluid equations of a species: the fluid solve on the
 * stored face fluxes, the source and the energy relaxation of the temperature source, momentum
 * diffusion and, with explicit sources, the electromagnetic force.
 *
 * @param app Vlasov app object
 * @param fl Fluid aspect of the species
 * @param fluidin Input fluid state
 * @param em EM field
 * @param rhs On output, the right-hand side of the fluid equations
 * @return Maximum stable time-step
 */
double vm_pkpm_fluid_species_rhs(
  gkyl_vlasov_app *app, struct vm_fluid_species *fl, const struct gkyl_array *fluidin,
  const struct gkyl_array *em, struct gkyl_array *rhs
);

/**
 * Complete a time step of a species: keep the closures of its last stage for the first stage of
 * the next step.
 *
 * @param app Vlasov app object
 * @param sp Species object
 */
void vm_pkpm_species_complete_step(gkyl_vlasov_app *app, struct vlasov_species *sp);

/**
 * Complete the stage of a species after the forward-Euler accumulate: the gyrophase-mixing
 * integrating factor of the first harmonic and the positivity limiter of the fluid state.
 *
 * @param app Vlasov app object
 * @param sp Species object
 * @param stage Stage of the SSP-RK3 step (1, 2 or 3)
 * @param dt Time-step
 * @param fout Distribution functions of the completed stage
 * @param fluidout Fluid state of the completed stage
 */
void vm_pkpm_species_finish_stage(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, double dt, struct gkyl_array *fout,
  struct gkyl_array *fluidout
);

/**
 * Apply BCs to species distribution functions 
 *
 * @param app Vlasov app object
 * @param kin Kinetic aspect of the species
 * @param f Distribution function to apply BCs to
 */
void vm_pkpm_species_apply_bc(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct gkyl_array *f
);

/**
 * Apply BCs to the fluid state of a species
 *
 * @param app Vlasov app object
 * @param fl Fluid aspect of the species
 * @param fluid Fluid state [rho, M, E] to apply BCs to
 */
void vm_pkpm_fluid_species_apply_bc(
  gkyl_vlasov_app *app, const struct vm_fluid_species *fl, struct gkyl_array *fluid
);

/**
 * Delete resources used in species.
 *
 * @param app Vlasov app object
 * @param species Species object to delete
 */
void vm_pkpm_species_release(const gkyl_vlasov_app *app, struct vlasov_species *s);
