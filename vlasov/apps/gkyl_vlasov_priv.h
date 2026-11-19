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
  // Staggered operator split (0: none): the step of the last update and the step it suggested, so a
  // restart continues as the uninterrupted run. Written when nonzero (field frames only).
  double dt_prev, dt_suggested;
};

// geometry data
struct vm_geom {
  struct gkyl_vlasov_geom info; // data for vlasov geometry
  double spin_bh, mass_bh; // Black hole spin and mass.
  bool use_preset_geom; // bool to determine if we are using triad input geom
  enum gkyl_triad_preset_geom_type
    triad_preset_geom_type; // geom type for preset geometries for triads
  int theta_pole_lo[GKYL_MAX_CDIM]; // (lower bound) Determines if the theta pole BC is being used
  int theta_pole_up[GKYL_MAX_CDIM]; // (upper bound) Determines if the theta pole BC is being used

  // Geometry needed for GR-DG-Maxwells
  bool has_gr_fields; // Boolean for determining if we have fields for GR-DG-Maxwells
  bool has_gr_em_triad_coupling; // Boolean for GR-DG-Maxwell coupled to triad species
  struct gkyl_surf_and_vol_node_arrays *lapse; // lapse scalar (ADM \alpha)
  struct gkyl_surf_and_vol_node_arrays
    *shift; // shift vector - contravariant radial component (ADM \beta^r)
  struct gkyl_surf_and_vol_node_arrays *geom_factor_con; // contravariant geometric source factors
  struct gkyl_surf_and_vol_node_arrays *h_ij; // Spatial metric, covariant components, h_ij
  struct gkyl_surf_and_vol_node_arrays *h_ij_inv; // Spatial metric, contravariant components, h^ij
  struct gkyl_surf_and_vol_node_arrays
    *det_h; // Square root of the spatial determinant from Jc = sqrt(det(h_ij))
  struct gkyl_surf_and_vol_node_arrays *vierb_cov; // Covariant vierbein components
  struct gkyl_surf_and_vol_node_arrays *vierb_con; // Contravariant vierbein components

  // Geometry copy for initialization (for GPU only)
  struct gkyl_surf_and_vol_node_arrays *lapse_init; // lapse scalar (ADM \alpha)
  struct gkyl_surf_and_vol_node_arrays
    *shift_init; // shift vector - contravariant radial component (ADM \beta^r)
  struct gkyl_surf_and_vol_node_arrays
    *geom_factor_con_init; // contravariant geometric source factors
  struct gkyl_surf_and_vol_node_arrays *h_ij_init; // Spatial metric, covariant components, h_ij
  struct gkyl_surf_and_vol_node_arrays
    *h_ij_inv_init; // Spatial metric, contravariant components, h^ij
  struct gkyl_surf_and_vol_node_arrays
    *det_h_init; // Square root of the spatial determinant from Jc = sqrt(det(h_ij))
  struct gkyl_surf_and_vol_node_arrays *vierb_cov_init; // Covariant vierbein components
  struct gkyl_surf_and_vol_node_arrays *vierb_con_init; // Contravariant vierbein components
};

// Implicit fluid-EM coupling data: the fluid-bearing species, gathered in
// declaration order at construction.
struct vm_fluid_em_coupling {
  int num_fluid; // number of gathered fluid-bearing species
  struct vlasov_species *species[GKYL_MAX_SPECIES]; // gathered species, in declaration order
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each gathered species
  struct gkyl_dg_calc_fluid_em_coupling *slvr; // fluid-EM coupling solver

  // Currents of the solve besides the species': the applied current and the effective current
  // -eps0 c^2 curl B of Ampere's law, whose E advances only in the source step.
  struct gkyl_array *curl_rhs; // Maxwell right-hand side of the field at the start of the source step
  struct gkyl_array *curl_current; // applied current - eps0 * (E rows of curl_rhs)

  // Staggered states (fluids, field): those at the start of an update while it runs (for a
  // transport step that takes less than planned), and those kept while the app holds the
  // synchronized ones.
  struct gkyl_array *fluid_stag[GKYL_MAX_SPECIES];
  struct gkyl_array *em_stag;

  // Ghost current (1x periodic; allocated with it): the fluids and the field at the start of a
  // source step, and the solve of the step (with or without the ghost current of the fluids).
  struct gkyl_array *fluid_ghost[GKYL_MAX_SPECIES];
  struct gkyl_array *em_ghost;
  void (*advance_func)(
    struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double dt,
    const struct gkyl_array *app_accels[GKYL_MAX_SPECIES],
    struct gkyl_array *fluids[GKYL_MAX_SPECIES]
  );

  // Kinetic species of a hybrid run: their current at the start of a source step is a known
  // current of the step (a leapfrog with the transport step, which moves them with the E of the
  // source step), so the step has to resolve their plasma frequency.
  int num_kinetic; // number of kinetic species
  struct vlasov_species *kinetic[GKYL_MAX_SPECIES]; // kinetic species, in declaration order
  struct gkyl_array *omega_p_sq; // per cell: sum over the kinetic species of q^2 n / (eps0 m)
  struct gkyl_array *m0_avg; // per cell: average density of a kinetic species (scratch)
  double *omega_p_sq_max; // maximum of omega_p_sq over the local range
  // Largest step of the leapfrog of the kinetic currents (DBL_MAX without kinetic species).
  double (*src_dt_func)(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em);
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
  struct gkyl_comm *comm; // communicator object for conf-space arrays

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom
    *geom; // geometry needed for species and field solvers (*only* p=1 right now JJ: 11/24/23)

  struct gkyl_vlasov_position_map *
    pos_map; // Configuration-space mapping object (owns all position map arrays); shared by all species.

  bool has_field; // has field
  struct vm_field *field; // pointer to field object (its dispatch methods are
    // set by vlasov_field_new; callers use vlasov_field_*).

  // geometry data
  struct vm_geom *vm_geom;

  // Species data: one array of species containers in input order.
  int num_species; // Total number of species (kinetic + fluid).
  struct vlasov_species *species; // Species containers.

  // Momentum-carrying fluids with the Maxwell field: fully explicit coupling in the SSP-RK3
  // stages instead of the implicit fluid-EM coupling of the staggered operator split.
  bool use_explicit_sources;

  bool has_fluid_em_coupling; // Boolean for if there is implicit fluid-EM coupling
  struct vm_fluid_em_coupling *fl_em; // fluid-EM coupling data

  // Staggered operator split: between updates the source-updated state (fluids, E) lags tcurr by
  // half the last step. Output and diagnostics synchronize it (vlasov_op_split_sync) and the next
  // update restores the staggered state.
  double op_split_dt_prev; // step taken by the previous operator-split update (0: none yet)
  double op_split_dt_suggested; // step suggested by the previous operator-split update (0: none yet)
  bool op_split_synced; // the fluids and the field hold the synchronized state

  // Boolean for implicit operators of the species in the RK3 stages (BGK collisions,
  // temperature sources)
  bool has_implicit_species_ops;

  // Explicit fluid-EM coupling (use_explicit_sources with momentum-carrying fluids and the
  // Maxwell field): the plasma frequency squared of the coupled fluids in each cell, for the rate
  // of the source in the time-step of each fluid; computed every stage by src_omega_p_func (a
  // no-op otherwise, and the arrays are NULL).
  void (*src_omega_p_func)(gkyl_vlasov_app *app, const struct gkyl_array *fluidin[]);
  struct gkyl_array *src_omega_p_sq; // per cell: sum over the fluids of (q/m)^2 rho / eps0
  struct gkyl_array *src_rho_avg; // per cell: average density of a fluid (scratch)

  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_vlasov_app *app, double dt0);

  struct gkyl_vlasov_stat stat; // statistics
};

// Take a single forward Euler step of the Vlasov-Maxwell system
// with the suggested time-step dt. Also supports just Maxwell's equations
// and fluid equations (Euler's) with potential Vlasov-fluid coupling.
// 'stage' is the SSP-RK3 stage (1, 2 or 3) the step belongs to. The stage
// distributions fin[] are not const (the signature of calc_self_moms_func); no
// species modifies them.
void vlasov_forward_euler(
  gkyl_vlasov_app *app, int stage, double tcurr, double dt, struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *fout[],
  struct gkyl_array *fluidout[], struct gkyl_array *emout, struct gkyl_update_status *st
);

// Explicit fluid-EM coupling: nothing to prepare (assigned to the app's src_omega_p_func).
void vlasov_src_omega_p_none(gkyl_vlasov_app *app, const struct gkyl_array *fluidin[]);

// Explicit fluid-EM coupling: the plasma frequency squared of the coupled fluids in each cell,
// the sum of (q/m)^2 rho / eps0 with the cell averages of the stage densities, into
// app->src_omega_p_sq. With the cyclotron frequency of a fluid it bounds the frequencies of its
// source (see gkyl_dg_calc_fluid_vars_em_source), which enters the time-step of the fluid. Zero
// for a static field, which does not respond to the currents. Kinetic species are not included.
void vlasov_src_omega_p_explicit(gkyl_vlasov_app *app, const struct gkyl_array *fluidin[]);

// The implicit operators of SSP-RK3 stage 'stage' (1, 2 or 3), applied to the
// stage state after the explicit part of the stage, its BCs and the stage
// combination: implicit BGK collisions and implicit temperature sources, one
// backward-Euler step of both over the share of dt the stage adds to its state
// (dt, dt/4, 2 dt/3). A no-op without implicit operators.
void vlasov_stage_implicit(gkyl_vlasov_app *app, int stage, double dt);

// Take a single time-step of the staggered operator split: one implicit fluid-EM source step over
// the interval between the midpoints of consecutive transport steps, and the SSP-RK3 transport step.
struct gkyl_update_status vlasov_update_op_split(gkyl_vlasov_app *app, double dt0);

// Staggered split: advance the source-updated state (fluids, E) by the half source step it lags
// tcurr, keeping the staggered state for the next update. A no-op when there is no implicit
// fluid-EM coupling, when already synchronized, and before the first update.
void vlasov_op_split_sync(gkyl_vlasov_app *app);

// Staggered split: restore the staggered state saved by vlasov_op_split_sync (no-op if not
// synchronized).
void vlasov_op_split_unsync(gkyl_vlasov_app *app);

// Staggered split: write the staggered field of a frame (<app>-field_restart_<frame>.gkyl), as the
// stepper holds it, for a restart. A no-op when there is no implicit fluid-EM coupling.
void vlasov_op_split_write_field(gkyl_vlasov_app *app, double tm, int frame);

// Staggered split: write the staggered fluid state of a species of a frame
// (<app>-<species>_restart_<frame>.gkyl). A no-op for a species the source step does not advance.
void vlasov_op_split_write_species(
  gkyl_vlasov_app *app, const struct vlasov_species *sp, double tm, int frame
);

// Staggered split: after the frames of a restart are read, replace the synchronized fluids and
// field by the staggered ones of the frame and restore the steps of the split.
struct gkyl_app_restart_status vlasov_op_split_restart(gkyl_vlasov_app *app, int frame);

// Take a single time-step using a SSP-RK3 stepper
struct gkyl_update_status vlasov_update_ssp_rk3(gkyl_vlasov_app *app, double dt0);

// Attempt a single SSP-RK3 step of size dt0: as vlasov_update_ssp_rk3, but when a stage allows
// only a smaller step it returns with the solution untouched and dt_actual < dt0 in the status,
// for the caller to prepare the smaller step (the staggered operator split repeats its source
// step).
struct gkyl_update_status vlasov_update_ssp_rk3_attempt(gkyl_vlasov_app *app, double dt0);

/** gkyl_vlasov_app private API */

/**
 * Create a new array metadata object. It must be freed using
 * vlasov_array_meta_release.
 *
 * @param meta Vlasov metadata object.
 * @return Array metadata object.
 */
struct gkyl_msgpack_data *vlasov_array_meta_new(struct vlasov_output_meta meta);

/**
 * Read the grid/array header of a restart file.
 *
 * @param app App object.
 * @param fname File to read the header from.
 * @return Restart status: IO status plus the frame number and simulation time
 *   from the file's embedded metadata.
 */
struct gkyl_app_restart_status vlasov_header_from_file(gkyl_vlasov_app *app, const char *fname);

/**
 * Read the metadata of an output file.
 *
 * @param fname File to read the metadata from.
 * @return Metadata of the file (zeros if the file cannot be read or has none).
 */
struct vlasov_output_meta vlasov_meta_from_file(const char *fname);

/**
 * Free memory for array metadata object.
 *
 * @param mt Array metadata object.
 */
void vlasov_array_meta_release(struct gkyl_msgpack_data *mt);

/**
 * Return the metadata for outputing vlasov data.
 *
 * @param mt Array metadata object.
 * @return A vlasov metadata object.
 */
struct vlasov_output_meta vlasov_meta_from_mpack(struct gkyl_msgpack_data *mt);

/**
 * Apply BCs to kinetic species, fluid species and EM fields.
 *
 * @param app Top-level Vlasov app.
 * @param tcurr Current simulation time.
 * @param distf Array of distribution functions (for each species).
 * @param fluid Array of moments (for each species).
 * @param emfield Electromagnetic fields.
 */
void vm_apply_bc(
  gkyl_vlasov_app *app, double tcurr, struct gkyl_array *distf[], struct gkyl_array *fluid[],
  struct gkyl_array *emfield
);

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
struct vm_species *vm_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Find the kinetic species with the given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of the species in app->species, -1 if not found
 */
int vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm);

/**
 * Configuration-space computational-to-physical map for projections
 * (proj_on_basis_c2p_t / eval_on_nodes_c2p_t): maps the computational
 * coordinates of the projection's nodes to physical ones with the position map,
 * so user functions of physical coordinates mean the same thing with and without
 * a map (the identity map returns xcomp unchanged). Phase-space projections use
 * vm_proj_c2p_phase.
 *
 * @param xcomp Computational coordinates
 * @param xphys On output, physical coordinates
 * @param ctx Position map (const struct gkyl_vlasov_position_map *)
 */
void vm_pos_c2p(const double *xcomp, double *xphys, void *ctx);

/** vm_fluid_em_coupling API */

/**
 * Create new fluid-EM coupling updater
 *
 * @param app Vlasov app object
 * @return Newly created fluid-EM coupling updater
 */
struct vm_fluid_em_coupling *vm_fluid_em_coupling_new(struct gkyl_vlasov_app *app);

/**
 * Whether a fluid species takes part in the fluid-EM coupling: only
 * momentum-carrying fluids (Euler and isothermal Euler) do.
 *
 * @param f Fluid species
 * @return True if the fluid species couples to the EM field
 */
bool vm_fluid_em_coupling_supported(const struct vm_fluid_species *f);

/**
 * Implicit source step of the fluid-EM coupling (implicit midpoint at the nodes): advances the
 * fluid momenta and energies and the electric field in place. E advances only here, with
 * c^2 curl B of the field and the current of the kinetic species at the start of the step as
 * known currents (Stormer-Verlet in the curl terms and a leapfrog in the kinetic currents; the
 * transport step drifts B and moves the kinetic species with this E).
 *
 * @param app Vlasov app object
 * @param fl_em fluid-EM coupling updater
 * @param tcurr Time at which time-dependent applied fields, currents and accelerations are evaluated
 * @param dt Time step size
 */
void vm_fluid_em_coupling_update(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double tcurr, double dt
);

/**
 * Largest step the leapfrog of the kinetic currents allows: cfl_frac / omega_p, with omega_p^2
 * the largest sum over the kinetic species of the plasma frequency squared of the cell-average
 * density (the leapfrog is stable for omega_p dt < 2 at a fixed step). DBL_MAX when the run has
 * no kinetic species.
 *
 * @param app Vlasov app object
 * @param fl_em fluid-EM coupling updater
 * @return Largest step of the source update
 */
double vm_fluid_em_coupling_dt(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em);

/**
 * Limiters after a source step of an update: the vertex positivity limiter of the fluids and
 * the slope limiter of the field. Not part of the source step itself: the synchronization for
 * output applies the source step alone.
 *
 * @param app Vlasov app object
 * @param fl_em fluid-EM coupling updater
 */
void vm_fluid_em_coupling_limit(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em);

/**
 * Release resources allocated by fluid-EM coupling object
 *
 * @param app Vlasov app object
 * @param fl_em fluid-EM coupling updater to release
 */
void vm_fluid_em_coupling_release(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em);

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
void vm_geom_release(const gkyl_vlasov_app *app, const struct vm_geom *s);
