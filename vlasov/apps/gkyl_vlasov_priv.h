// Private header for use in Vlasov app: do not include in user-facing
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
#include <gkyl_dg_calc_canonical_pb_fluid_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
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
#include <gkyl_dg_gr_maxwell_geom_source.h>
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
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_maxwell.h>

// The field object (struct vm_field) and its API. Included here so the app
// struct below can hold a struct vm_field *; struct vm_geom (used only as a
// pointer in struct vm_field) is forward-declared there and defined below.
#include <gkyl_vlasov_field_priv.h>

// Definitions of private structs and APIs attached to these objects
// for use in Vlasov app.

// Meta-data for IO
struct vlasov_output_meta {
  int frame; // frame number
  double stime; // output time
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
  char basis_type_nm[64]; // used during read
};

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
      struct gkyl_vlasov_lte_moments *vlasov_lte_moms; // Epdater for computing LTE moments.
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

struct vm_proj {
  enum gkyl_projection_id proj_id; // type of projection
  enum gkyl_model_id model_id;
  struct vm_proj_c2p_ctx c2p_ctx; // coordinate-map context for projections (captured by proj_on_basis)
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

      // Correction updater for insuring LTE distribution has desired LTE (n, V_drift, T/m) moments
      bool correct_all_moms; // boolean if we are correcting all the moments
      struct gkyl_vlasov_lte_correct *corr_lte;    
    };
  };
};

struct vm_lte {  
  struct gkyl_array *f_lte;

  struct vm_species_moment moms; // moments needed in the equilibrium

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
  bool use_extended_hamil_def; // bool to determine if we are using the extended hamil defintions which includes potentials

  double qbym; // Charge (q) divided by mass (m).
  struct gkyl_array *qmem; // array for q/m*(E,B)
  struct gkyl_array *em_no_J; // physical E,B (stored field J*E,J*B divided by the conf Jacobian) for the Lorentz force on a mapped grid
  struct gkyl_array *pot_tot; // array for total potentials (q/m*phi + m*phi_g, q/m*A)
  bool has_E; // Do we have electric fields? 
  bool has_phi; // Do we have scalar potentials (electrostatic/gravitational)?
  bool has_B; // Do we have magnetic fields? 
  bool has_gr_em_triad_coupling; // Do we need GR-Maxwell Lorentz-force fields for triad species?
  int num_surf_conf_nodes; // number of surface nodes at configuration-space surfaces
  int num_surf_vel_nodes; // number of surface nodes at velocity-space surfaces

  struct gkyl_array *conf_flux_surf; // Modal expansion of surface fluxes at conf-space surfaces. 
  struct gkyl_array *vel_flux_surf; // Modal expansion of surface fluxes at velocity-space surfaces. 
  struct gkyl_dg_vlasov_conf_flux_surf *calc_conf_flux; // Updater for computing modal expansion of surface fluxes (conf). 
  struct gkyl_dg_vlasov_vel_flux_surf *calc_vel_flux; // Updater for computing modal expansion of surface fluxes (vel).   
  struct gkyl_dg_gr_maxwell_lorentz_conf *calc_lorentz; // Updater for local GR Lorentz-force fields.
  struct gkyl_dg_gr_maxwell_current_deposition *calc_current_dep; // Updater for GR current deposition.

  struct gkyl_dg_eqn *eqn; // Vlasov equation object.
  struct gkyl_hyper_dg *slvr; // Vlasov solver.  
  
  void (*rhs_func)(gkyl_vlasov_app *app, struct vm_species *vms,
    struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
    struct gkyl_array *rhs);  
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
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value.
  struct gkyl_array *nu_boundary_corrections; // Boundary corrections multiplied by nu.
  struct gkyl_array *nu_moms; // Moments multiplied by nu.
  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator

  int num_cross_collisions; // number of species we cross-collide with
  struct vm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
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
  
  struct vm_species_moment moms; // Moments needed in LBO (M0, M1, M2).

  gkyl_dg_updater_collisions *coll_slvr; // collision solver

  // Pointers to methods chosen at runtime.
  void (*moms_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_lbo_collisions *lbo, const struct gkyl_array *fin);
  void (*self_nu_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_lbo_collisions *lbo, const struct gkyl_array *fin);
  void (*cross_nu_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_lbo_collisions *lbo, int coll_idx);
  void (*alpha_E_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_lbo_collisions *lbo, int coll_idx);
  void (*cross_moms_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_lbo_collisions *lbo);
  void (*rhs_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs);
  void (*write_mom_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
};

struct vm_bgk_collisions {  
  enum gkyl_collision_id collision_id; // Type of collisions.
  bool write_coll_diagnostics; // Whether to write diagnostics out.

  struct gkyl_array *self_nu; // Self-collision frequency.
  struct gkyl_array *ref_self_nu; // Reference self-collision frequency needed for robustness checks n, T < 0.
  struct gkyl_array *nu_sum; // Sum of collision frequencies.
  struct gkyl_array *nu_sum_host; // Host arrays for I/O.
  bool norm_nu_self; // Whether to compute self-species collision frequency in space and time.
  double norm_nu_fac_self; // Self collision frequency without factor of n_s/(2*v_ts^2)^(3/2).
  double vtsq_min; // Minimum vtsq.
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  int num_cross_collisions; // number of species we cross-collide with
  struct vm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  bool norm_nu_cross; // Whether to compute cross-species collision frequency in space and time.
  double norm_nu_fac_cross[GKYL_MAX_SPECIES]; // Cross collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2).
  double alpha_E_fac[GKYL_MAX_SPECIES]; // Time-independent factor in alpha_E.
  double betaGreenep1; // Galue of Greene's factor beta + 1.
  double delta_sr; // Free parameter in relationship between alpha_E and nu_sr.
  double other_m[GKYL_MAX_SPECIES]; // Masses of species colliding with.
  struct gkyl_array *other_prim_moms[GKYL_MAX_SPECIES]; // Self-primitive moments of species colliding with.
  struct gkyl_array *cross_prim_moms; // Cross-primitive moments.
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // Cross-species collision frequencies.
  struct gkyl_array *ref_cross_nu[GKYL_MAX_SPECIES]; // Reference cross-collision frequencies needed for robustness checks n, T < 0.
  struct gkyl_array *alpha_E; // Morse's alpha_E factor.
  struct gkyl_vlasov_cross_prim_moms_bgk *cross_calc; // Cross-species moment computation.
  
  struct gkyl_array *nu_f_lte; // Collision frequency times Maxwellian.
  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep).

  bool fixed_temp_relax; // Boolean for whether the temperature being relaxed to is fixed in time.
  struct gkyl_array *fixed_temp; // Array of fixed temperature BGK collisions are relaxing to.

  bool implicit_step; // Whether or not to take an implcit BGK step.

  // Pointers to methods chosen at runtime.
  void (*moms_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*moms_func_implicit)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*self_nu_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk);
  void (*cross_nu_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk, int coll_idx);
  void (*alpha_E_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk, int coll_idx);
  void (*cross_moms_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk, int coll_idx);
  void (*rhs_func)(gkyl_vlasov_app *app, struct vm_species *vms,
    struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);
  void (*rhs_func_implicit)(gkyl_vlasov_app *app, struct vm_species *vms,
    struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);
  void (*fixed_temp_calc_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct vm_bgk_collisions *bgk, const struct gkyl_array *fin);
  void (*write_mom_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
};

struct vm_boundary_fluxes {
  struct gkyl_rect_grid boundary_grid[2*GKYL_MAX_CDIM];
  struct gkyl_array *flux_arr[2*GKYL_MAX_CDIM];
  struct gkyl_array *mom_arr[2*GKYL_MAX_CDIM];
  struct gkyl_range flux_r[2*GKYL_MAX_CDIM];
  struct gkyl_range conf_r[2*GKYL_MAX_CDIM];
  struct gkyl_mom_type *mom_type;
  struct gkyl_mom_calc *integ_moms[2*GKYL_MAX_CDIM];
  gkyl_ghost_surf_calc *flux_slvr; // boundary flux solver
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
  struct vm_species *impact_species[GKYL_MAX_SPECIES]; // pointers to impacting species
  struct gkyl_range impact_normal_r[GKYL_MAX_SPECIES];
  struct gkyl_mom_type *mom_type;
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
  double *scale_ptr;
  struct vm_species *source_species; // species to use for the source
  int source_species_idx; // index of source species

  bool rescale_m0; // boolean for if we are rescaling M0 
  int num_cross_source; // how many other species are we obtaining sources from?
  struct vm_species *adapt_source_species[GKYL_MAX_SPECIES]; // list of species to use for the source
  int adapt_source_species_idx[GKYL_MAX_SPECIES]; // list of indices of source species
  struct gkyl_array *scale_m0[GKYL_MAX_SPECIES]; // Time-dependent re-scaling of the density of the source. 
  struct gkyl_array *scale_m0_host[GKYL_MAX_SPECIES]; // host copy for use in IO
  struct gkyl_mom_calc *m0_reduced[GKYL_MAX_SPECIES]; // Reduced density update for rescaling source. 
  struct gkyl_array *adapt_source[GKYL_MAX_SPECIES]; // adaptive source array
  int adapt_proj_source[GKYL_MAX_SPECIES]; // Index of projection function to use for adaptive source. 

  bool filter; // boolean for if we are filtering recaled M0
  int num_filters; // number of times to apply filter
  gkyl_dg_gaussian_filter *gauss_filter; // updater for filtering rescaled M0
  
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

// geometry data
struct vm_geom {
  struct gkyl_vlasov_geom info; // data for vlasov geometry
  double spin_bh, mass_bh; // Charge and mass.
  bool use_preset_geom; // bool to determine if we are using triad input geom
  enum gkyl_triad_preset_geom_type triad_preset_geom_type; // geom type for preset geometries for triads
  int theta_pole_lo[GKYL_MAX_CDIM]; // (lower bound) Determines if the theta pole BC is being used
  int theta_pole_up[GKYL_MAX_CDIM]; // (lower bound) Determines if the theta pole BC is being used

  // Geometry needed for GR-DG-Maxwells
  bool has_gr_fields; // Boolean for determining if we have fields for GR-DG-Maxwells
  bool has_gr_em_triad_coupling; // Boolean for GR-DG-Maxwell coupled to triad species
  struct gkyl_surf_and_vol_node_arrays *lapse; // lapse scalar (ADM \alpha)
  struct gkyl_surf_and_vol_node_arrays *shift; // shift vector - contravaraint radial component (ADM \beta^r)
  struct gkyl_surf_and_vol_node_arrays *geom_factor_con; // contravariant geometric source factors
  struct gkyl_surf_and_vol_node_arrays *h_ij; // Spatial metric, covaraint components, h_ij
  struct gkyl_surf_and_vol_node_arrays *h_ij_inv; // Spatial metric, contravariant components, h^ij
  struct gkyl_surf_and_vol_node_arrays *det_h; // Squareroot of the spatial determinant from Jc = sqrt(det(h_ij))
  struct gkyl_surf_and_vol_node_arrays *vierb_cov; // Covariant vierbein components
  struct gkyl_surf_and_vol_node_arrays *vierb_con; // Contravariant vierbein components

  // Geometry copy for initalization (for GPU only)
  struct gkyl_surf_and_vol_node_arrays *lapse_init; // lapse scalar (ADM \alpha)
  struct gkyl_surf_and_vol_node_arrays *shift_init; // shift vector - contravaraint radial component (ADM \beta^r)
  struct gkyl_surf_and_vol_node_arrays *geom_factor_con_init; // contravariant geometric source factors
  struct gkyl_surf_and_vol_node_arrays *h_ij_init; // Spatial metric, covaraint components, h_ij
  struct gkyl_surf_and_vol_node_arrays *h_ij_inv_init; // Spatial metric, contravariant components, h^ij
  struct gkyl_surf_and_vol_node_arrays *det_h_init; // Squareroot of the spatial determinant from Jc = sqrt(det(h_ij))
  struct gkyl_surf_and_vol_node_arrays *vierb_cov_init; // Covariant vierbein components
  struct gkyl_surf_and_vol_node_arrays *vierb_con_init; // Contravariant vierbein components

};

// species data
struct vm_species {
  struct gkyl_vlasov_species info; // data for species

  struct gkyl_basis basis; // Phase-space basis. 
  struct gkyl_basis basis_vel; // Velocity-space basis. 
  struct gkyl_basis basis_surf; // Surface basis at velocity-space surfaces for velocity-space fluxes. 
  struct gkyl_basis *basis_on_dev; // Pointer to phase-space basis on device.   

  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // fixed buffers for time independent BCs 

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
  enum gkyl_triad_preset_geom_type triad_preset_geom_type; // geom type for preset geometries for triads

  struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object (owns all velocity map arrays).
  struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (acquired reference to the app's pos_map).

  struct gkyl_array *f_no_J; // Distribution function without velocity-space Jacobian.
                             // When using uniform velocity-space mesh, just stores the distribution function at that RK stage. 

  // Organization of the different equation objects and the required data.
  bool has_rad; // Do we have a radiation drag force?
  struct gkyl_array *rad; // array for radiation drag force. 

  // Geometry
  struct vm_geom *geom; // Geometry structure for vm

  // Organization of the different equation objects and the required data.
  struct gkyl_range mom_hamil_range; // Range Hamiltonian (for moments) is defined over (only velocity-space or all phase-space).
  struct gkyl_array *mom_hamil; // Specified Hamiltonian (for moments) function for canonical poisson bracket.
  struct gkyl_range hamil_range; // Range Hamiltonian is defined over (only velocity-space or all phase-space).
  struct gkyl_array *hamil; // Specified Hamiltonian function for canonical poisson bracket.
  struct gkyl_array *hamil_host; // Host-side Hamiltonian array for initial projection.
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation
  struct gkyl_array *conf_poisson_tensor_host; // Host-side configuration space Poisson tensor representation  
  union {
    // Special relativistic Vlasov-Maxwell model.
    struct {
      struct gkyl_array *gamma_inv; // array for 1/gamma = 1.0/sqrt(1 + p^2) 
      struct gkyl_dg_calc_sr_vars *sr_vars; // updater for computing SR variables
    };
    // Canonical Poisson Bracket using specified Hamiltonian in phase space. 
    struct {
      struct gkyl_array *h_ij; // Specified metric inverse for canonical poisson bracket
      struct gkyl_array *h_ij_host; // Host side metric inverse array for intial projection
      struct gkyl_array *h_ij_inv; // Specified metric inverse for canonical poisson bracket
      struct gkyl_array *h_ij_inv_host; // Host side metric inverse array for intial projection
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

  bool calc_bflux; // are we computing boundary fluxes?
  struct vm_boundary_fluxes bflux; // boundary flux object

  enum gkyl_source_id source_id; // type of source
  struct vm_source src; // applied source

  double *omega_cfl;

  // Pointer to various functions selected at runtime.
  void (*apply_ic_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm);
  double (*rhs_func)(gkyl_vlasov_app *app, struct vm_species *vms,
    const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs);
  double (*rhs_implicit_func)(gkyl_vlasov_app *app, struct vm_species *vms,
    const struct gkyl_array *fin, struct gkyl_array *rhs, double dt);
  void (*bc_func)(gkyl_vlasov_app *app, const struct vm_species *vms,
    struct gkyl_array *f, double tcurr);
  void (*release_func)(const gkyl_vlasov_app* app, const struct vm_species *vms);
  void (*step_f_func)(struct gkyl_array* out, double dt, const struct gkyl_array* inp); 
  void (*combine_func)(struct gkyl_array *out, double c1,
    const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
    const struct gkyl_range *rng);
  void (*copy_func)(struct gkyl_array *out, const struct gkyl_array *inp,
    const struct gkyl_range *range);
  void (*write_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
  void (*write_lte_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
  void (*write_cfl_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
  void (*write_cell_avg_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
  void (*write_mom_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm);
  void (*write_integrated_mom_func)(gkyl_vlasov_app* app, struct vm_species *vms);
  void (*calc_L2_func)(gkyl_vlasov_app* app, struct vm_species *vms, double tm);
  void (*write_L2_func)(gkyl_vlasov_app* app, struct vm_species *vms);  
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

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *fluid, *fluid1, *fluidnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)

  struct gkyl_array *fluid_host;  // host copy for use IO and initialization

  // Duplicate copy of fluid data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  struct gkyl_array *fluid_dup;  

  enum gkyl_eqn_type eqn_type;  // type ID of equation
  int num_equations;            // number of equations in species
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
      // For isothermal Euler, u : (ux, uy, uz), p : (vth*rho)
      // For Euler, u : (ux, uy, uz, T/m), p : (gamma - 1)*(E - 1/2 rho u^2)
      // Also a prim_vars and prim_vars_host array for I/O of (u,p)
      struct gkyl_array *u; 
      struct gkyl_array *p; 
      struct gkyl_array *prim_vars; 
      struct gkyl_array *prim_vars_host; 
      struct gkyl_array *cell_avg_prim; // Integer array for whether e.g., rho *only* uses cell averages for weak division
                                        // Determined when constructing the matrix if rho < 0.0 at control points

      // Arrays for kinetic energy at old and new time steps.
      // These are used because implicit source solve updates momentum but does not affect 
      // the pressure, so we can construct the updated energy from the updated momentum.
      struct gkyl_array *ke_old; 
      struct gkyl_array *ke_new; 

      struct gkyl_array *u_surf; 
      struct gkyl_array *p_surf;
      struct gkyl_dg_calc_fluid_vars *calc_fluid_vars; // Updater to compute fluid variables (flow velocity and pressure)
      struct gkyl_dg_calc_fluid_vars *calc_fluid_vars_ext; // Updater to compute fluid variables (flow velocity and pressure)
                                                           // over extended range (used when BCs are not absorbing to minimize apply BCs calls) 
    };
    // Canonical PB Fluid such as incompressible Euler or Hasegawa-Wakatani
    struct {
      struct gkyl_array *phi; // potential determined by canonical PB Poisson equation on local range used by updater
      struct gkyl_array *phi_global; // potential determined by canonical PB Poisson equation on global range given by Poisson solver
      struct gkyl_array *poisson_rhs_global; // global RHS of Poisson equation, simply an all-gather of, e.g., the vorticity
      struct gkyl_array *phi_host; // host copy for use IO
      struct gkyl_array *can_pb_n0; // background density gradient for driving turbulence in some fluid systems. 
      struct gkyl_array *epsilon; // Permittivity in Poisson equation, set to -1.0 for canonical PB Poisson equations. 
      struct gkyl_array *kSq; // k^2 factor in Helmholtz equation needed for Hasegawa-Mima where we solve (grad^2 - 1) phi = RHS

      struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
                                          // for solving Poisson equation on each MPI process in parallel
    
      struct gkyl_fem_poisson *fem_poisson; // Poisson solver for - nabla . (epsilon * nabla phi) - kSq * phi = rho.

      struct gkyl_array *alpha_surf; // Surface configuration space velocity (derivatives of potential, phi)
      struct gkyl_array *sgn_alpha_surf; // sign(alpha_surf) at quadrature points
      struct gkyl_array *const_sgn_alpha; // boolean for if sign(alpha_surf) is a constant, either +1 or -1
      struct gkyl_dg_calc_canonical_pb_fluid_vars *calc_can_pb_fluid_vars; // Updater for computing surface alpha and sources. 
      struct gkyl_array *can_pb_energy_fac; // Factor in calculation of canonical PB energy diagnostic.
      struct gkyl_array_integrate *calc_can_pb_energy;
      double *red_can_pb_energy, *red_can_pb_energy_global; // Memory for use in GPU reduction of canonical PB energy.
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

  // fluid source
  enum gkyl_source_id source_id; // type of source
  struct vm_fluid_source src; // applied source

  double* omegaCfl_ptr;

  // Function pointers for computing primitive/auxiliary variables, 
  // and also write method, release method, and method for calculating integrated quantities. 
  void (*prim_vars_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid); 
  void (*calc_integrated_mom_func)(gkyl_vlasov_app* app, struct vm_fluid_species *f, double tm); 
  void (*write_func)(gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm, int frame); 
  void (*release_func)(const gkyl_vlasov_app *app, struct vm_fluid_species *f);   
};

// fluid-EM coupling data
struct vm_fluid_em_coupling {
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each species
  struct gkyl_dg_calc_fluid_em_coupling* slvr; // fluid-EM coupling solver
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_app {
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

  struct gkyl_basis basis; // Configuration-space basis. 
  struct gkyl_basis *basis_on_dev; // Pointer to configuration-space basis on device. 

  struct gkyl_rect_decomp *decomp; // decomposition object
  struct gkyl_comm *comm;   // communicator object for conf-space arrays

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx;   // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom *geom; // geometry needed for species and field solvers (*only* p=1 right now JJ: 11/24/23)

  struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (owns all position map arrays); shared by all species.

  bool has_field; // has field
  struct vm_field *field; // pointer to field object (its dispatch methods are
                          // set by vlasov_field_new; callers use vlasov_field_*).

  // geometry data
  struct vm_geom *vm_geom;

  // species data
  int num_species;
  struct vm_species *species; // data for each species
  
  // fluid data
  int num_fluid_species;
  struct vm_fluid_species *fluid_species; // data for each fluid species

  bool has_fluid_em_coupling; // Boolean for if there is implicit fluid-EM coupling
  struct vm_fluid_em_coupling *fl_em; // fluid-EM coupling data

  bool has_implicit_coll_scheme; // Boolean for using implicit bgk scheme (over explicit rk3)

  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_vlasov_app *app, double dt0);

  struct gkyl_vlasov_stat stat; // statistics
};

// Take a single forward Euler step of the Vlasov-Maxwell system 
// with the suggested time-step dt. Also supports just Maxwell's equations
// and fluid equations (Euler's) with potential Vlasov-fluid coupling. 
void vlasov_forward_euler(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, 
  struct gkyl_update_status *st);

// Calls the vlasov implicit contribution for all vm species
void vlasov_update_implicit_coll(gkyl_vlasov_app *app,  double dt0);

// Take a single time-step using a first-order operator split 
// implicit fluid-EM coupling and/or implicit BGK collisions + SSP RK3
struct gkyl_update_status vlasov_update_op_split(gkyl_vlasov_app *app,  double dt0);

// Take a single time-step using a SSP-RK3 stepper
struct gkyl_update_status vlasov_update_ssp_rk3(gkyl_vlasov_app *app,
  double dt0);

/** gkyl_vlasov_app private API */

/**
 * Create a new array metadata object. It must be freed using
 * vlasov_array_meta_release.
 *
 * @param meta Vlasov metadata object.
 * @return Array metadata object.
 */
struct gkyl_msgpack_data*
vlasov_array_meta_new(struct vlasov_output_meta meta);

/**
 * Free memory for array metadata object.
 *
 * @param mt Array metadata object.
 */
void
vlasov_array_meta_release(struct gkyl_msgpack_data *mt);

/**
 * Return the metadata for outputing vlasov data.
 *
 * @param mt Array metadata object.
 * @return A vlasov metadata object.
 */
struct vlasov_output_meta
vlasov_meta_from_mpack(struct gkyl_msgpack_data *mt);

/**
 * Apply BCs to kinetic species, fluid species and EM fields.
 *
 * @param app Top-level Vlasov app.
 * @param tcurr Current simulation time.
 * @param distf Array of distribution functions (for each species).
 * @param fluid Array of moments (for each species).
 * @param emfield Electromagnetic fields.
 */
void vm_apply_bc(gkyl_vlasov_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *fluid[], struct gkyl_array *emfield);

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct vm_species* vm_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of species, -1 if not found
 */
int vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm);

/**
 * Find fluid species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of fluid species
 * @return Pointer to fluid species with given name. NULL if not found.o
 */
struct vm_fluid_species *vm_find_fluid_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Return index fluid species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of fluid species
 * @return Index of species, -1 if not found
 */
int vm_find_fluid_species_idx(const gkyl_vlasov_app *app, const char *nm);


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
void vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param sm vm species moment
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release species moment object.
 *
 * @param app Vlasov app object
 * @param sm Species moment object to release
 */
void vm_species_moment_release(const struct gkyl_vlasov_app *app,
  const struct vm_species_moment *sm);

/** vm_species_projection API */

/**
 * Initialize species projection object.
 *
 * @param app vlasov app object
 * @param s Species object 
 * @param inp Input struct for projection (contains functions pointers for type of projection)
 * @param proj Species projection object
 */
void vm_species_projection_init(struct gkyl_vlasov_app *app, struct vm_species *s, 
  struct gkyl_vlasov_projection inp, struct vm_proj *proj);

/**
 * Compute species projection
 *
 * @param app vlasov app object
 * @param species Species object
 * @param proj Species projection object
 * @param f Output distribution function from projection
 * @param tm Time for use in projection
 */
void vm_species_projection_calc(gkyl_vlasov_app *app, const struct vm_species *species, 
  struct vm_proj *proj, struct gkyl_array *f, double tm);

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
void vm_species_lte_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_lte *lte, struct correct_all_moms_inp corr_inp);

/**
 * Compute LTE distribution from input moments
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lte Pointer to lte object
 * @param moms_lte Input LTE moments
 */
void vm_species_lte_from_moms(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lte *lte,
  const struct gkyl_array *moms_lte);

/**
 * Compute equivalent LTE distribution from input distribution function. 
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lte Pointer to lte
 * @param fin Input distribution function
 */
void vm_species_lte(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lte *lte,
  const struct gkyl_array *fin);

/**
 * Write the LTE correction status. 
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 */
void vm_species_lte_write_max_corr_status(gkyl_vlasov_app* app, struct vm_species *vms);

/**
 * Release species lte object.
 *
 * @param app Vlasov app object
 * @param lte Species lte object to release
 */
void vm_species_lte_release(const struct gkyl_vlasov_app *app, const struct vm_lte *lte);

/** vm_species_collisionless API */

/**
 * Initialize species collisionless object.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param cls Species collisionless object.
 */
void vm_species_collisionless_init(struct gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls);

/**
 * Compute species applied acceleration term
 *
 * @param app Vlasov app object
 * @param cls Species collisionless object
 * @param tm Time for use in acceleration
 */
void vm_species_collisionless_app_accel(gkyl_vlasov_app *app, struct vm_collisionless *cls, double tm);

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
void vm_species_collisionless_rhs(gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs);

/**
 * Release species collisionless object.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param cls Species collisionless object.
 */  
void vm_species_collisionless_release(const struct gkyl_vlasov_app *app, 
  const struct vm_species *vms, const struct vm_collisionless *cls);

/** vm_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void vm_species_lbo_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vm_species_lbo_moms(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vm_species_lbo_cross_moms(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute RHS from LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 * @return Maximum stable time-step
 */
void vm_species_lbo_rhs(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write moments from LBO object.
 *
 * @param app Vlasov app object.
 * @param species Pointer to species.
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void vm_species_lbo_write_mom(gkyl_vlasov_app* app, struct vm_species *species, double tm, int frame);

/**
 * Release species LBO object.
 *
 * @param app Vlasov app object
 * @param sm Species LBO object to release
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
void vm_species_bgk_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_bgk_collisions *bgk);

/**
 * Initialize species BGK cross-collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param bgk Species BGK object
 */
void vm_species_bgk_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_bgk_collisions *bgk);

/**
 * Compute necessary moments for BGK collisions.
 *
 * @param app Vlasov app object
 * @param vms Pointer to species.
 * @param bgk Pointer to BGK.
 * @param fin Input distribution function.
 */
void vm_species_bgk_moms(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin);

void vm_species_bgk_moms_implicit(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin);

/**
 * Compute and store a fixed temperature for BGK collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void vm_species_bgk_moms_fixed_temp(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin);

/**
 * Compute RHS from BGK collisions.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 * @param rhs On output, the RHS from bgk
 */
void vm_species_bgk_rhs(gkyl_vlasov_app *app,
  struct vm_species *species,
  struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

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
void vm_species_bgk_rhs_implicit(gkyl_vlasov_app *app,
  struct vm_species *species,
  struct vm_bgk_collisions *bgk,
  const struct gkyl_array *fin, double dt, struct gkyl_array *rhs);

/**
 * Write moments from BGK object.
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param tm Simulation time.
 * @param frame Simulation output frame.
 */
void vm_species_bgk_write_mom(gkyl_vlasov_app* app, struct vm_species *species, double tm, int frame);

/**
 * Release species BGK object.
 *
 * @param app Vlasov app object
 * @param bgk Species BGK object to release
 */
void vm_species_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_bgk_collisions *bgk);

/** vm_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param bflux Species boundary flux object
 */
void vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_boundary_fluxes *bflux);

/**
 * Compute boundary flux from rhs
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bflux Species boundary flux object
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void vm_species_bflux_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species boundary flux object.
 *
 * @param app Vlasov app object
 * @param bflux Species boundary flux object to release
 */
void vm_species_bflux_release(const struct gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux);

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
void vm_species_emission_init(struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit,
  int dir, enum gkyl_edge_loc edge, void *ctx);

/**
 * Initialize emission BC cross-species object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param emit Species emission object
 */
void vm_species_emission_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_emitting_wall *emit);

/**
 * Apply emission BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param emit Pointer to emission object
 * @param fout Field to apply BCs
 * @param tcurr Current time
 */
void vm_species_emission_apply_bc(struct gkyl_vlasov_app *app, const struct vm_species *s, 
  const struct vm_emitting_wall *emit, struct gkyl_array *fout, double tcurr);

/**
 * Write emission spectrum distribution function
 *
 * @param app Vlasov app object
 * @param s Species object
 * @param emit Pointer to emission object
 * @param mt Write meta
 * @param frame Current frame
 */
void vm_species_emission_write(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_emitting_wall *emit, struct gkyl_msgpack_data *mt, int frame);

/**
 * Release species emission object.
 *
 * @param app Vlasov app object
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
void vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src);

/**
 * Compute species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 * @param tm Time for use in source
 */
void vm_species_source_calc(gkyl_vlasov_app *app, const struct vm_species *species, 
  struct vm_source *src, double tm);

/**
 * Compute density re-scaling for adaptive sourcing. 
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param tm Time for use in source
 */
void vm_species_source_adapt_moms(gkyl_vlasov_app *app, const struct vm_species *species, 
  struct vm_source *src, const struct gkyl_array *fin);

/**
 * Adapt source based on density re-scaling. 
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param tm Time for use in source
 */
void vm_species_source_adapt(gkyl_vlasov_app *app, const struct vm_species *species, 
  struct vm_source *src);

/**
 * Compute RHS contribution from source
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function
 */
void vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[]);

/**
 * Write source diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void vm_species_source_write(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src, double tm, int frame);

/**
 * Write source moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void vm_species_source_write_mom(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src, double tm, int frame);

/**
 * Calc source integrated moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 * @param tm Time for source diagnostic.
 */
void vm_species_source_calc_integrated_mom(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src, double tm);

/**
 * Write source integrated moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param vms Pointer to species.
 * @param src Pointer to source.
 */
void vm_species_source_write_integrated_mom(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src);

/**
 * Release species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src);

/** vm_species API */

/**
 * Initialize geom.
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param s On output, initialized geom object
 */
void vm_geom_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_geom *s);


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
double vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs);

/**
 * Compute the *implicit* RHS from species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param rhs On output, the RHS from the species object
 * @param dt timestep size (used in the implcit coef.)
 * @return Maximum stable time-step
 */
double vm_species_rhs_implicit(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt);

/**
 * Scale and accumulate for forward euler method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param dt Timestep.
 * @param inp Input array.
 */
void vm_species_step_f(struct vm_species *vms, struct gkyl_array* out, double a,
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
void vm_species_combine(struct vm_species *vms, struct gkyl_array *out, double c1,
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
void vm_species_copy_range(struct vm_species *vms, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Apply BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param f Field to apply BCs
 * @param tcurr Current time
 */
void vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f, double tcurr);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers
 */
void vm_species_coll_tm(gkyl_vlasov_app *app);

/**
 * Species write function.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void vm_species_write(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);

/**
 * Species moment write function.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void vm_species_write_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame);

/**
 * Species calc integrated moment function.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 * @param tm Time at which diagnostic is computed.
 */
void vm_species_calc_integrated_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm);

/**
 * Species write integrated moment function.
 *
 * @param app Vlasov app object.
 * @param vms Species object.
 */
void vm_species_write_integrated_mom(gkyl_vlasov_app* app, struct vm_species *vms);

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
void vm_species_write_L2(gkyl_vlasov_app* app, struct vm_species *vms);

/**
 * Fill stat object in app with total number of iterations
 * used to correct moments in LTE projection object.
 * Also fills stat object with number of times correction object called. 
 *
 * @param app App object to update stat timers.
 */
void vm_species_n_iter_corr(gkyl_vlasov_app *app);

/**
 * Delete resources used in geom.
 *
 * @param app Vlasov app object
 * @param species Geom object to delete
 */
void vm_geom_release(const gkyl_vlasov_app* app, const struct vm_geom *s);

/**
 * Delete resources used in species.
 *
 * @param app Vlasov app object
 * @param species Species object to delete
 */
void vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s);


/**
 * Initialize fluid species source object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param src Species source object
 */
void vm_fluid_species_source_init(struct gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct vm_fluid_source *src);

/**
 * Compute fluid species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param tm Time for use in source
 */
void vm_fluid_species_source_calc(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function RHS
 */
void vm_fluid_species_source_rhs(gkyl_vlasov_app *app, const struct vm_fluid_species *species,
  struct vm_fluid_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[]);

/**
 * Release fluid species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vm_fluid_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_fluid_source *src);

/** vm_fluid_species API */

/**
 * Create new fluid species object
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param f On output, initialized fluid species object
 */
void vm_fluid_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f);

/**
 * Compute fluid species initial conditions.
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param t0 Time for use in ICs
 */
void vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0);

/**
 * Compute fluid species applied acceleration term
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param tm Time for use in acceleration
 */
void vm_fluid_species_calc_app_accel(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm);

/**
 * Compute primitive variables (bulk velocity, u, and pressure, p, if pressure present)
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object (where primitive variables are stored)
 * @param fluid Input array fluid species
 */
void vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid);

/**
 * Limit slopes of solution of fluid variables
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species (where primitive variables are stored)
 * @param fluid Input (and Output after limiting) array fluid species
 */
void vm_fluid_species_limiter(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  struct gkyl_array *fluid);

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
double vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, 
  const struct gkyl_array *fluid, const struct gkyl_array *em, 
  struct gkyl_array *rhs);

/**
 * Apply BCs to fluid species
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param f Fluid Species to apply BCs
 */
void vm_fluid_species_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f);

/**
 * Computed the integrated quantities for the fluid system.
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param tm Time integrated quantities are being computed at. 
 */
void vm_fluid_species_calc_integrated_mom(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm);

/**
 * Write out the evolved fluid species and other potential primitive/auxiliary variables.
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param tm Time fluid quantities are being written at.
 * @param frame Frame number for I/O.  
 */
void vm_fluid_species_write(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm, int frame);

/**
 * Write fluid species integrated moment diagnostics.
 *
 * @param app Vlasov app object.
 * @param fluid_species Pointer to fluid species.
 */
void vm_fluid_species_write_integrated_mom(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species);

/**
 * Release resources allocated by fluid species
 *
 * @param app Vlasov app object
 * @param f Fluid_Species object to release
 */
void vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f);

/** vm_fluid_em_coupling API */

/**
 * Create new fluid-EM coupling updater
 *
 * @param app Vlasov app object
 * @return Newly created fluid-EM coupling updater
 */
struct vm_fluid_em_coupling* vm_fluid_em_coupling_init(struct gkyl_vlasov_app *app);

/**
 * Compute implicit update of fluid-EM coupling 
 *
 * @param app Vlasov app object
 * @param fl_em fluid-EM coupling updater
 * @param tcurr Current time
 * @param dt Time step size
 */
void vm_fluid_em_coupling_update(struct gkyl_vlasov_app *app, 
  struct vm_fluid_em_coupling *fl_em, double tcurr, double dt);

/**
 * Release resources allocated by fluid-EM coupling object
 *
 * @param app Vlasov app object
 * @param fl_em fluid-EM coupling updater to release
 */
void vm_fluid_em_coupling_release(struct gkyl_vlasov_app *app, 
  struct vm_fluid_em_coupling *fl_em);
