// Private header for use in Vlasov app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

// App-level dependencies; the field and species private headers below carry
// the dependencies of their own layers.
#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_fluid_em_coupling.h>
#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>

// The field and species objects and their private APIs.
#include <gkyl_vlasov_field_priv.h>
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
  double spin_bh, mass_bh; // Black hole spin and mass.
  bool use_preset_geom; // bool to determine if we are using triad input geom
  enum gkyl_triad_preset_geom_type triad_preset_geom_type; // geom type for preset geometries for triads
  int theta_pole_lo[GKYL_MAX_CDIM]; // (lower bound) Determines if the theta pole BC is being used
  int theta_pole_up[GKYL_MAX_CDIM]; // (upper bound) Determines if the theta pole BC is being used

  // Geometry needed for GR-DG-Maxwells
  bool has_gr_fields; // Boolean for determining if we have fields for GR-DG-Maxwells
  bool has_gr_em_triad_coupling; // Boolean for GR-DG-Maxwell coupled to triad species
  struct gkyl_surf_and_vol_node_arrays *lapse; // lapse scalar (ADM \alpha)
  struct gkyl_surf_and_vol_node_arrays *shift; // shift vector - contravariant radial component (ADM \beta^r)
  struct gkyl_surf_and_vol_node_arrays *geom_factor_con; // contravariant geometric source factors
  struct gkyl_surf_and_vol_node_arrays *h_ij; // Spatial metric, covariant components, h_ij
  struct gkyl_surf_and_vol_node_arrays *h_ij_inv; // Spatial metric, contravariant components, h^ij
  struct gkyl_surf_and_vol_node_arrays *det_h; // Square root of the spatial determinant from Jc = sqrt(det(h_ij))
  struct gkyl_surf_and_vol_node_arrays *vierb_cov; // Covariant vierbein components
  struct gkyl_surf_and_vol_node_arrays *vierb_con; // Contravariant vierbein components

  // Geometry copy for initialization (for GPU only)
  struct gkyl_surf_and_vol_node_arrays *lapse_init; // lapse scalar (ADM \alpha)
  struct gkyl_surf_and_vol_node_arrays *shift_init; // shift vector - contravariant radial component (ADM \beta^r)
  struct gkyl_surf_and_vol_node_arrays *geom_factor_con_init; // contravariant geometric source factors
  struct gkyl_surf_and_vol_node_arrays *h_ij_init; // Spatial metric, covariant components, h_ij
  struct gkyl_surf_and_vol_node_arrays *h_ij_inv_init; // Spatial metric, contravariant components, h^ij
  struct gkyl_surf_and_vol_node_arrays *det_h_init; // Square root of the spatial determinant from Jc = sqrt(det(h_ij))
  struct gkyl_surf_and_vol_node_arrays *vierb_cov_init; // Covariant vierbein components
  struct gkyl_surf_and_vol_node_arrays *vierb_con_init; // Contravariant vierbein components

};

// Implicit fluid-EM coupling data: the fluid-bearing species, gathered in
// declaration order at construction.
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

  // Species data: one array of species containers in input order.
  int num_species; // Total number of species (kinetic + fluid).
  struct vlasov_species *species; // Species containers.

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

// The implicit half of the operator-split step: implicit BGK collisions for each
// species, then the implicit fluid-EM coupling, both with the explicit step's
// actual dt. A no-op when neither implicit scheme is active.
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
 * Find the species with the given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of the species in app->species, -1 if not found
 */
int vlasov_find_species_idx(const gkyl_vlasov_app *app, const char *nm);

/**
 * Find the kinetic species with the given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to the kinetic aspect of the species. NULL if not found.
 */
struct vm_species* vm_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Find the kinetic species with the given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of the species in app->species, -1 if not found
 */
int vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm);

/** vm_fluid_em_coupling API */

/**
 * Create new fluid-EM coupling updater
 *
 * @param app Vlasov app object
 * @return Newly created fluid-EM coupling updater
 */
struct vm_fluid_em_coupling* vm_fluid_em_coupling_new(struct gkyl_vlasov_app *app);

/**
 * Whether a fluid species takes part in the implicit fluid-EM coupling: only
 * momentum-carrying fluids (Euler) do.
 *
 * @param f Fluid species
 * @return True if the fluid species couples to the EM field
 */
bool vm_fluid_em_coupling_supported(const struct vm_fluid_species *f);

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

/** vm_geom API */

/**
 * Initialize geom.
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param s On output, initialized geom object
 */
void vm_geom_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_geom *s);



/**
 * Delete resources used in geom.
 *
 * @param app Vlasov app object
 * @param s Geom object to delete
 */
void vm_geom_release(const gkyl_vlasov_app* app, const struct vm_geom *s);
