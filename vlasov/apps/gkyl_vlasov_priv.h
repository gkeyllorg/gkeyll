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

// The species objects (struct vm_species, struct vm_fluid_species) and the
// unified container (struct vlasov_species), with their building-block structs.
// Included here so the app struct below can hold a struct vlasov_species *;
// struct vm_geom (held only as a pointer) is forward-declared there and
// defined below.
#include <gkyl_vlasov_species_priv.h>

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

// Implicit fluid-EM coupling data: an app-level assembly object. Participation
// is a property of the species (it has a fluid aspect, under a dynamic Maxwell
// field), so construction gathers the fluid-bearing species from app->species[]
// by scanning for that property -- in declaration order, giving the grid-local
// kernel the dense, stably-indexed lists it needs -- rather than assuming any
// layout of the species array.
struct vm_fluid_em_coupling {
  int num_fluid; // number of gathered fluid-bearing species
  struct vlasov_species *species[GKYL_MAX_SPECIES]; // gathered species, in declaration order
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each gathered species
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

  // species data. One backing array of unified containers: kinetic species
  // occupy indices [0, num_species), fluid species [num_species, num_species +
  // num_fluid_species). 'species' points at the head, 'fluid_species' is a view
  // into the fluid-typed tail (so the fluid view stays dense/contiguous, which
  // restart IO relies on; the fluid-EM coupling gathers its own species list by
  // property scan). 'species' owns the allocation.
  int num_species;
  struct vlasov_species *species; // unified container array (owns the backing storage)

  // fluid data
  int num_fluid_species;
  struct vlasov_species *fluid_species; // view into the fluid tail of 'species'

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

// The implicit half of the op-split step: per-species implicit collisions
// (BGK) followed by the holistic implicit fluid-EM coupling, both taken with
// the explicit step's actual dt. Checks its own participation flags; a no-op
// when neither implicit scheme is active.
void vlasov_update_implicit(gkyl_vlasov_app *app, double dt0);

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
 * Read the grid/array header of a restart file.
 *
 * @param app App object.
 * @param fname File to read the header from.
 * @return Restart status: IO status plus the frame number and simulation time
 *   from the file's embedded metadata.
 */
struct gkyl_app_restart_status vlasov_header_from_file(gkyl_vlasov_app *app,
  const char *fname);

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

/** The species/fluid private API (vm_species_* / vm_fluid_species_*) and the
 *  unified vlasov_species_* dispatch API now live in gkyl_vlasov_species_priv.h. */

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
