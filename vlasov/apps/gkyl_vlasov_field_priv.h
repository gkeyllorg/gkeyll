#pragma once

// Private header for the Vlasov field object: the field struct shared by
// Vlasov-Maxwell and Vlasov-Poisson (with its dispatch function pointers), the
// type-agnostic vlasov_field_* API, and the vm_field_* / vp_field_* APIs.

#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_gr_maxwell_conf_flux_surf.h>
#include <gkyl_dg_gr_maxwell_geom_source.h>
#include <gkyl_dynvec.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_job_pool.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_vlasov.h>

struct vm_geom; // geometry data, defined in gkyl_vlasov_priv.h, owned by the app.
struct gkyl_vlasov_position_map; // configuration-space position map, defined in zero/, owned by the app.

// Coordinate-map context for projecting external fields/potentials on a
// non-uniform configuration mesh: the projection's computational quadrature
// coordinates are mapped to physical ones via the position map before the
// user-supplied function is evaluated (the identity map is transparent).
struct vm_field_proj_c2p_ctx {
  const struct gkyl_vlasov_position_map *pos_map; // configuration-space position map.
};

// field data
struct vm_field {
  struct gkyl_vlasov_field info; // data for field
  enum gkyl_field_id field_id; // Type of field.

  // Pointers to various functions selected at runtime (by field type).
  double (*update_func)(gkyl_vlasov_app *app, double tcurr,
    const struct gkyl_array *fin[], const struct gkyl_array *emin, struct gkyl_array *emout);
  void (*combine_func)(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out,
    double c1, const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2);
  void (*copy_range_func)(gkyl_vlasov_app *app, struct vm_field *field,
    struct gkyl_array *out, const struct gkyl_array *inp);
  void (*apply_ic_func)(gkyl_vlasov_app *app, struct vm_field *field,
    const struct gkyl_array *fin[], double t0);
  void (*apply_bc_func)(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *em);
  void (*limiter_func)(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *em);
  void (*complete_update_func)(gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
    const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout);
  void (*calc_ext_em_func)(gkyl_vlasov_app *app, struct vm_field *field, double tm);
  void (*calc_app_current_func)(gkyl_vlasov_app *app, struct vm_field *field, double tm);
  void (*calc_ext_pot_func)(gkyl_vlasov_app *app, struct vm_field *field, double tm);
  void (*calc_energy_func)(gkyl_vlasov_app *app, double tm, struct vm_field *field,
    const struct gkyl_array *fin[]);
  void (*write_func)(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[]);
  void (*write_energy_func)(gkyl_vlasov_app *app);
  // Restart the field from the named file. Vlasov-Maxwell reads the EM field;
  // Vlasov-Poisson evaluates its static external potentials/fields (the
  // potential is solved from the distribution whenever it is needed); the null
  // field is a no-op.
  struct gkyl_app_restart_status (*from_file_func)(gkyl_vlasov_app *app, struct vm_field *field,
    const char *fname);
  void (*release_func)(const gkyl_vlasov_app *app, struct vm_field *field);

  union {
    // Vlasov-Maxwell.
    struct {
      struct gkyl_job_pool *job_pool; // Job pool
      struct gkyl_array *em, *em1, *emnew; // arrays for updates
      struct gkyl_array *cflrate; // CFL rate in each cell
      struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)
      struct gkyl_array *bc_buffer_lo_fixed[3], *bc_buffer_up_fixed[3]; // fixed buffers for time independent BCs
      struct gkyl_array *bc_buffer_lo_fixed_no_J[3], *bc_buffer_up_fixed_no_J[3]; // fixed buffers for GR fields without Jc

      struct gkyl_array *em_host;  // host copy for use IO and initialization

      // Duplicate copy of EM data in case time step fails.
      // Needed because of implicit source split which modifies solution and
      // is always successful, so if a time step fails due to the SSP RK3
      // we must restore the old solution before restarting the time step
      struct gkyl_array *em_dup;

      bool has_sigma; // flag to indicate there is a resistive layer
      struct gkyl_array *sigma; // resistive layer for damping EM fields
      struct gkyl_array *sigmaEM; // resistivity*E/B for incrementing onto RHS

      gkyl_hyper_dg *slvr; // Maxwell solver

      bool limit_em; // boolean for whether or not we are limiting EM fields
      struct gkyl_dg_calc_em_vars *calc_em_vars; // Updater to limit EM fields

      struct gkyl_array *em_energy; // EM energy components in each cell
      double *em_energy_red; // memory for use in GPU reduction of EM energy

      bool use_ghost_current; // Are we using ghost currents to correct dE/dt = -J in 1x
      struct gkyl_array *ghost_current; // Array for storying global average of current density
      double *red_ghost_current; // memory for use in GPU reduction of average of current density

      bool use_geom_sources; // Are we using geometric sources to correct dE/dt = -J in 1x?
      struct gkyl_array *geom_source; // Geometric source contribution to field RHS.
      struct gkyl_dg_gr_maxwell_geom_source *calc_geom_source; // Updater for geometric source contribution.

      // boundary conditions on lower/upper edges in each direction
      enum gkyl_field_bc_type lower_bc[3], upper_bc[3];
      // Pointers to updaters that apply BC.
      struct gkyl_bc_basic *bc_lo[3];
      struct gkyl_bc_basic *bc_up[3];

      double* omega_cfl;
    };

    // Vlasov-Poisson.
    struct {
      struct gkyl_array *epsilon;  // Permittivity in Poisson equation.

      struct gkyl_array *rho_c, *rho_c_global; // Local and global charge density.
      struct gkyl_array *phi, *phi_global; // Local and global potential.

      struct gkyl_array *phi_host;  // host copy for use IO and initialization

      struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
                                          // for solving subset of Poisson solves with parallelization in z

      struct gkyl_fem_poisson *fem_poisson; // Poisson solver for - nabla . (epsilon * nabla phi) - kSq * phi = rho.

      struct gkyl_array *es_energy_fac; // Factor in calculation of ES energy diagnostic.
      struct gkyl_array_integrate *calc_es_energy;
      double *es_energy_red, *es_energy_red_global; // Memory for use in GPU reduction of ES energy.
    };
  };

  struct vm_geom *geom; // Geometry data for GR-DG-Maxwell (owned by app as app->vm_geom)
  bool weight_by_pos_jacob; // True for the standard E_B Maxwell field on a non-identity position map:
                            // em stores J*E, J*B; em_no_J holds the physical E, B for force/I/O.
  struct gkyl_array *em_no_J; // arrays for storing em field without Jc
  struct gkyl_array *em_no_J_host; // host copy of primitive GR fields for I/O
  int num_surf_conf_nodes; // number of surface nodes at configuration-space surfaces
  struct gkyl_array *conf_flux_surf; // Modal expansion of surface fluxes at conf-space surfaces.
  struct gkyl_dg_gr_maxwell_conf_flux_surf *calc_conf_flux; // Updater for computing modal expansion of surface fluxes (conf).

  bool has_ext_em; // flag to indicate there are external electromagnetic fields (E, B)
  bool ext_em_evolve; // flag to indicate external electromagnetic fields are time dependent
  struct gkyl_array *ext_em; // external electromagnetic fields
  struct gkyl_array *ext_em_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *ext_em_proj; // projector for external electromagnetic field

  bool has_app_current; // flag to indicate there is an applied current
  bool app_current_evolve; // flag to indicate applied current is time dependent
  struct gkyl_array *app_current; // applied current
  struct gkyl_array *app_current_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_current_proj; // projector for applied current

  bool has_ext_pot; // flag to indicate there are external potentials (phi, A)
  bool ext_pot_evolve; // flag to indicate external external potentials are time dependent
  struct gkyl_array *ext_pot; // external potentials
  struct gkyl_array *ext_pot_host; // host copy for use in IO and projecting
  gkyl_eval_on_nodes *ext_pot_proj; // projector for external potentials
  struct vm_field_proj_c2p_ctx ext_c2p_ctx; // comp->phys map for external-field projection on mapped grids

  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for energy dynvec written first time
};

/** vlasov_field API: type-agnostic operations on the field object (vlasov_field.c). */

/**
 * Create the field object, dispatching on the field type (Maxwell, Poisson, or
 * the null field when no field is present).
 *
 * @param vm App inputs
 * @param app Vlasov app object
 * @return New field object
 */
struct vm_field* vlasov_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app);

/**
 * Update the field at the current time: the RHS of Maxwell's equations, or the
 * Poisson solve for the potential from the charge density.
 *
 * @param app Vlasov app object
 * @param tcurr Current time
 * @param fin[] Input distribution functions (overall species count)
 * @param emin Input EM field
 * @param emout On output, the field RHS
 * @return Maximum stable time-step (DBL_MAX for the Poisson solve)
 */
double vlasov_field_update(gkyl_vlasov_app *app, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *emin, struct gkyl_array *emout);

/**
 * Combine RK stages of the field state, out = c1*arr1 + c2*arr2 (no-op for
 * Vlasov-Poisson).
 *
 * @param app Vlasov app object
 * @param out Output array
 * @param c1 Coefficient of arr1
 * @param arr1 First input array
 * @param c2 Coefficient of arr2
 * @param arr2 Second input array
 */
void vlasov_field_combine(gkyl_vlasov_app *app, struct gkyl_array *out,
  double c1, const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2);

/**
 * Copy the field state, out = inp (no-op for Vlasov-Poisson).
 *
 * @param app Vlasov app object
 * @param out Output array
 * @param inp Input array
 */
void vlasov_field_copy_range(gkyl_vlasov_app *app, struct gkyl_array *out, const struct gkyl_array *inp);

/**
 * Compute field initial conditions.
 *
 * @param app Vlasov app object
 * @param fin[] Input distribution functions (overall species count)
 * @param t0 Time for initial conditions
 */
void vlasov_field_apply_ic(gkyl_vlasov_app *app, const struct gkyl_array *fin[], double t0);

/**
 * Apply boundary conditions to the field (no-op for Vlasov-Poisson).
 *
 * @param app Vlasov app object
 * @param em Field to apply BCs to
 */
void vlasov_field_apply_bc(gkyl_vlasov_app *app, struct gkyl_array *em);

/**
 * Limit slopes of the field (no-op for Vlasov-Poisson).
 *
 * @param app Vlasov app object
 * @param em Field to limit
 */
void vlasov_field_limiter(gkyl_vlasov_app *app, struct gkyl_array *em);

/**
 * Complete the field update: for Vlasov-Maxwell, accumulate the species current
 * onto the RHS and finalize emout = emin + dt*RHS (no-op for Vlasov-Poisson).
 *
 * @param app Vlasov app object
 * @param dt Time-step
 * @param fin[] Input distribution functions (overall species count)
 * @param fluidin[] Input fluid moments (overall species count)
 * @param emin Input EM field
 * @param emout On input the field RHS, on output the updated field
 */
void vlasov_field_complete_update(gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout);

/**
 * Compute the external electromagnetic field.
 *
 * @param app Vlasov app object
 * @param tm Time at which to evaluate the field
 */
void vlasov_field_calc_ext_em(gkyl_vlasov_app *app, double tm);

/**
 * Compute the applied current.
 *
 * @param app Vlasov app object
 * @param tm Time at which to evaluate the current
 */
void vlasov_field_calc_app_current(gkyl_vlasov_app *app, double tm);

/**
 * Compute the external potentials.
 *
 * @param app Vlasov app object
 * @param tm Time at which to evaluate the potentials
 */
void vlasov_field_calc_ext_pot(gkyl_vlasov_app *app, double tm);

/**
 * Compute the field energy diagnostic. Vlasov-Poisson solves for the potential
 * from fin[] so the diagnostic is at time tm.
 *
 * @param app Vlasov app object
 * @param tm Time at which the diagnostic is computed
 * @param fin[] Input distribution functions (overall species count)
 */
void vlasov_field_calc_energy(gkyl_vlasov_app *app, double tm, const struct gkyl_array *fin[]);

/**
 * Write the field. Vlasov-Poisson solves for the potential from fin[] so the
 * written potential is at time tm.
 *
 * @param app Vlasov app object
 * @param tm Time-stamp
 * @param frame Frame number
 * @param fin[] Input distribution functions (overall species count)
 */
void vlasov_field_write(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[]);

/**
 * Write the field energy diagnostic.
 *
 * @param app Vlasov app object
 */
void vlasov_field_write_energy(gkyl_vlasov_app *app);

/**
 * Restart the field from the named file (dispatches by field type).
 *
 * @param app Vlasov app object
 * @param fname File to read
 * @return Restart status
 */
struct gkyl_app_restart_status vlasov_field_from_file(gkyl_vlasov_app *app, const char *fname);

/**
 * Restart the field from a frame: name the frame file, restart from it, and
 * mark the energy diagnostic to append.
 *
 * @param app Vlasov app object
 * @param frame Frame number
 * @return Restart status
 */
struct gkyl_app_restart_status vlasov_field_read_from_frame(gkyl_vlasov_app *app, int frame);

/**
 * Release the field object.
 *
 * @param app Vlasov app object
 */
void vlasov_field_release(gkyl_vlasov_app *app);

/** vm_field API (Vlasov-Maxwell): implemented in vm_field.c. */

/**
 * Create new field object
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @return Newly created field
 */
struct vm_field* vm_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app);

/**
 * Compute field initial conditions.
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param fin[] Input distribution functions; unused for Vlasov-Maxwell
 * @param t0 Time for use in ICs
 */
void vm_field_apply_ic(gkyl_vlasov_app *app, struct vm_field *field,
  const struct gkyl_array *fin[], double t0);

/**
 * Compute external electromagnetic fields
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param tm Time for use in external electromagnetic fields computation
 */
void vm_field_calc_ext_em(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Compute applied currents
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param tm Time for use in applied current computation
 */
void vm_field_calc_app_current(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Compute external potentials
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param tm Time for use in external potentials computation
 */
void vm_field_calc_ext_pot(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Accumulate current density onto RHS from field equations. Loops over the
 * species and dispatches each species' explicit field coupling
 * (kinetic species accumulate -q/eps0*m1i; implicitly-coupled fluid species
 * are a no-op), then adds the applied current.
 *
 * @param app Vlasov app object
 * @param fin[] Input distribution functions, indexed over the overall species count (NULL for fluid-only species)
 * @param fluidin[] Input fluid arrays, indexed over the overall species count (NULL for kinetic-only species)
 * @param emout On output, the RHS from the field solver *with* accumulated current density
 */
void vm_field_accumulate_current(gkyl_vlasov_app *app,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], struct gkyl_array *emout);

/**
 * Accumulate geometric source terms onto RHS from field equations.
 *
 * @param app Vlasov app object
 * @param emin Input field at the start of the step
 * @param vm_geom Geometry data
 * @param emout On output, the RHS from the field solver *with* geometric sources
 */
void vm_field_accumulate_geom_sources(gkyl_vlasov_app *app,
  const struct gkyl_array *emin, const struct vm_geom *vm_geom, struct gkyl_array *emout);

/**
 * Limit slopes of solution of EM variables
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param em Input (and Output after limiting) EM fields
 */
void vm_field_limiter(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *em);

/**
 * Seed fixed-function field BC buffers from the stored field state.
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 */
void vm_field_buffer_fixed_func_bc(gkyl_vlasov_app *app, struct vm_field *field);

/**
 * Compute RHS from field equations
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param em Input field
 * @param rhs On output, the RHS from the field solver
 * @return Maximum stable time-step
 */
double vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em, struct gkyl_array *rhs);

/**
 * Update the field at the current time: compute the RHS of Maxwell's equations.
 * Assigned to the field's update_func; returns the max stable time-step.
 *
 * @param app Vlasov app object
 * @param tcurr Current time
 * @param fin[] Input distribution function (num_species size); unused for Maxwell
 * @param emin Input field
 * @param emout On output, the RHS from the field solver
 * @return Maximum stable time-step
 */
double vm_field_update(gkyl_vlasov_app *app, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *emin, struct gkyl_array *emout);

/**
 * Complete the field update: accumulate the species current onto the RHS
 * (unless the field is static) and finalize emout = emin + dt*RHS.
 *
 * @param app Vlasov app object
 * @param dt Time-step taken
 * @param fin[] Input distribution function (num_species size)
 * @param fluidin[] Input fluid array (num_fluid_species size)
 * @param emin Input field at the start of the step
 * @param emout On output, the completed field emin + dt*RHS
 */
void vm_field_complete_update(gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout);

/**
 * Combine RK stages of the field state (out = c1*arr1 + c2*arr2).
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param out On output, c1*arr1 + c2*arr2
 * @param c1 Coefficient of arr1
 * @param arr1 First field array
 * @param c2 Coefficient of arr2
 * @param arr2 Second field array
 */
void vm_field_combine(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out,
  double c1, const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2);

/**
 * Copy the field state (out = inp).
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param out On output, a copy of inp
 * @param inp Field array to copy
 */
void vm_field_copy_range(gkyl_vlasov_app *app, struct vm_field *field,
  struct gkyl_array *out, const struct gkyl_array *inp);

/**
 * Apply BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param f Field to apply BCs
 */
void vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  struct gkyl_array *f);

/**
 * Write out electromagnetic fields.
 *
 * @param app Vlasov app object
 * @param tm Time-stamp
 * @param frame Frame number
 * @param fin[] Input distribution functions; unused for Vlasov-Maxwell
 */
void vm_field_write(gkyl_vlasov_app* app, double tm, int frame, const struct gkyl_array *fin[]);

/**
 * Compute field energy diagnostic.
 *
 * @param app Vlasov app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 * @param fin[] Input distribution functions; unused for Vlasov-Maxwell
 */
void vm_field_calc_energy(gkyl_vlasov_app *app, double tm, struct vm_field *field,
  const struct gkyl_array *fin[]);

/**
 * Write out electromagnetic field energy.
 *
 * @param app Vlasov app object
 */
void vm_field_write_energy(gkyl_vlasov_app* app);

/**
 * Restart the Vlasov-Maxwell EM field from the named file: seed the fixed-
 * function BC buffers from the initial conditions, read the field, rescale it
 * to the evolved (Jacobian-weighted) representation, apply BCs, and recompute
 * the external EM field and applied current.
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param fname File to read
 * @return Restart status
 */
struct gkyl_app_restart_status vm_field_from_file(gkyl_vlasov_app *app,
  struct vm_field *field, const char *fname);

/**
 * Release resources allocated by field
 *
 * @param app Vlasov app object
 * @param f Field object to release
 */
void vm_field_release(const gkyl_vlasov_app* app, struct vm_field *f);

/** vp_field API (Vlasov-Poisson): implemented in vp_field.c. */

/**
 * Create new field object.
 *
 * @param vm Input Vlasov data.
 * @param app Vlasov app object.
 * @return Newly created field.
 */
struct vm_field* vp_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app);

/**
 * Compute external electromagnetic fields.
 *
 * @param app Vlasov app object.
 * @param field Field object.
 * @param tm Time for use in external fields computation.
 */
void vp_field_calc_ext_em(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Compute external applied currents.
 *
 * @param app Vlasov app object.
 * @param field Field object.
 * @param tm Time for use in external applied currents computation.
 */
void vp_field_calc_app_current(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Compute external potentials.
 *
 * @param app Vlasov app object.
 * @param field Field object.
 * @param tm Time for use in external potentials computation.
 */
void vp_field_calc_ext_pot(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Vlasov-Poisson initial conditions: evaluate the static external potentials
 * and fields. The potential itself is not solved here; it is solved from the
 * distribution wherever it is needed (each forward Euler step and the field
 * diagnostics), so initial conditions and restart share this code path.
 *
 * @param app Vlasov app object.
 * @param field Field object.
 * @param fin[] Input distribution functions; unused.
 * @param t0 Time at which to evaluate the externals.
 */
void vp_field_apply_ic(gkyl_vlasov_app *app, struct vm_field *field,
  const struct gkyl_array *fin[], double t0);

/**
 * Restart for Vlasov-Poisson: no field file is read; evaluates the static
 * external potentials and fields exactly as the initial conditions do.
 *
 * @param app Vlasov app object.
 * @param field Field object.
 * @param fname Unused.
 * @return Restart status (always success; frame/time come from the species files).
 */
struct gkyl_app_restart_status vp_field_from_file(gkyl_vlasov_app *app,
  struct vm_field *field, const char *fname);

/**
 * Solve for the electrostatic potential from the distribution functions:
 * accumulate the charge density over all species, then solve the Poisson
 * equation. The single definition of the potential, used by the field update
 * (each forward Euler stage) and by the field diagnostics.
 *
 * @param app Vlasov app object.
 * @param field Pointer to field.
 * @param fin[] Input distribution functions, indexed over the overall species count.
 */
void vp_field_solve(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[]);

/**
 * Update the field at the current time: solve for the potential from the charge
 * density. Assigned to the field's update_func; the elliptic solve imposes no
 * CFL constraint, so it returns DBL_MAX.
 *
 * @param app Vlasov app object.
 * @param tcurr Current time.
 * @param fin[] Input distribution function (num_species size).
 * @param emin Unused (no EM RK state for Vlasov-Poisson).
 * @param emout Unused (no EM RK state for Vlasov-Poisson).
 * @return DBL_MAX (no stability constraint from the elliptic solve).
 */
double vp_field_update(gkyl_vlasov_app *app, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *emin, struct gkyl_array *emout);

/**
 * Complete the field update. No-op for Vlasov-Poisson (the potential is solved
 * from the distribution wherever it is needed, not carried in the RK state).
 */
void vp_field_complete_update(gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout);

/**
 * Combine RK stages of the field state. No-op for Vlasov-Poisson.
 */
void vp_field_combine(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out,
  double c1, const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2);

/**
 * Copy the field state. No-op for Vlasov-Poisson.
 */
void vp_field_copy_range(gkyl_vlasov_app *app, struct vm_field *field,
  struct gkyl_array *out, const struct gkyl_array *inp);

/**
 * Apply BCs to field. No-op for Vlasov-Poisson (the potential has no EM BCs).
 */
void vp_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *em);

/**
 * Limit slopes of EM variables. No-op for Vlasov-Poisson.
 */
void vp_field_limiter(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *em);

/**
 * Write out the potential (solved here from fin[] so it is at time tm) and the
 * external fields/potentials.
 *
 * @param app Vlasov app object
 * @param tm Time-stamp
 * @param frame Frame number
 * @param fin[] Input distribution functions, indexed over the overall species count
 */
void vp_field_write(gkyl_vlasov_app* app, double tm, int frame, const struct gkyl_array *fin[]);

/**
 * Compute the electrostatic field energy diagnostic (the potential is solved
 * here from fin[] so it is at time tm).
 *
 * @param app Vlasov app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 * @param fin[] Input distribution functions, indexed over the overall species count
 */
void vp_field_calc_energy(gkyl_vlasov_app *app, double tm, struct vm_field *field,
  const struct gkyl_array *fin[]);

/**
 * Write out potential field energy.
 *
 * @param app Vlasov app object
 */
void vp_field_write_energy(gkyl_vlasov_app* app);

/**
 * Release resources allocated by field.
 *
 * @param app Vlasov app object.
 * @param f Field object to release.
 */
void vp_field_release(const gkyl_vlasov_app* app, struct vm_field *f);
