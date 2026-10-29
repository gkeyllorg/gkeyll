#pragma once

#include <assert.h>
#include <float.h>
#include <string.h>

#include <gkyl_bc_basic.h>
#include <gkyl_comm.h>
#include <gkyl_dg_calc_flash_coupling.h>
#include <gkyl_dg_calc_pkpm_dist_vars.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_updater_flash.h>
#include <gkyl_dg_updater_lbo_pkpm.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_dynvec.h>
#include <gkyl_flash.h>
#include <gkyl_flash_interp_vars.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>

// Metadata struct for I/O
struct flash_io_meta {
  int frame; // Frame number.
  double time; // Output time of frame.
  const char *basis_type; // Name of basis set, or NULL for finite-volume data.
  int poly_order; // Polynomial order (ignored if basis_type is NULL).
  int num_cdim, num_vdim; // Phase-space split, written if num_vdim > 0.
};

// LBO collisions object.
struct flash_lbo_collisions {
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *prim_moms, *nu_prim_moms; // LBO primitive moments

  struct gkyl_array *self_nu; // LBO self collision frequency
  struct gkyl_array *m0;

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

// Fluid species object.
struct flash_fluid_species {
  struct gkyl_flash_species info;

  struct gkyl_rect_grid grid; // Configuration space grid.
  struct gkyl_range local, local_ext; // Local, local-ext configuration space ranges.
  struct gkyl_range global, global_ext; // Global, global-ext configuration space ranges.

  struct gkyl_array *fluid, *fluid1, *fluidnew; // Fluid variable arrays.
  struct gkyl_array *fluid_prim; // (rho, ux, uy, uz, p) for diagnositcs.
  struct gkyl_array *fluid_phi;

  struct gkyl_array *flash_coupling_fluxes; // Flux contribution terms, [div(P)-grad(p), (P - pI) : grad(u)].

  double gas_gamma;
  struct gkyl_wave_geom *wv_geom; // Geometry object for Euler equation.
  struct gkyl_wv_eqn *wv_eqn; // Euler wave equation object.
  struct gkyl_wave_prop *up_wave_prop; // Euler equation wave propagation updater.

  // Boundary conditions on lower/upper edges in each direction.
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  struct gkyl_wv_apply_bc *bc_lo[3];
  struct gkyl_wv_apply_bc *bc_up[3];
  bool bc_is_absorb;

  // Skin and ghost ranges for BC application.
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Species object.
struct flash_species {
  struct gkyl_flash_species info; // Storage for input parameters struct.

  double mass;

  struct gkyl_rect_grid grid; // Phase space grid.
  struct gkyl_rect_grid grid_vel; // Velocity space grid.

  struct gkyl_range local, local_ext; // Local, local-ext phase space ranges.
  struct gkyl_range global, global_ext; // Global, global-ext phase space ranges.
  struct gkyl_range local_vel, local_ext_vel; // Local, local-ext velocity space ranges.

  struct gkyl_rect_decomp *decomp; // Decomposition object.
  struct gkyl_comm *comm; // Communicator object for phase space.

  struct gkyl_array *f, *f1, *fnew; // Distribution function arrays.
  struct gkyl_array *L2_f; // Array to store L2 norm of the distribution function.
  gkyl_dynvec integ_L2_f; // Storage for integrated L^2 norm reduced across grid.
  struct gkyl_array *pkpm_div_ppar; // dif(p_par b_hat) for selv-consistent total pressure force.
  struct gkyl_array *pkpm_prim, *pkpm_p_ij; // Primitive moments.
  struct gkyl_array *cell_avg_prim; // Integer array for whether rho, p_par, or p_perp < 0.0 at control points
  struct gkyl_array *pkpm_prim_surf; // Surface expansions of primitive variables.
  struct gkyl_array *pkpm_lax, *pkpm_penalization; // Expansions of Lax penalization and momentum penalization.
  struct gkyl_array *pkpm_accel; // PKPM acceleration variables.
  struct gkyl_array *pkpm_moms, *pkpm_moms_diag; // PKPM moments (rho, p_par, p_perp, M1).

  struct gkyl_array *flash_coupling_fluxes; // Flux contribution terms, [div(P)-grad(p), (P - pI) : grad(u)].
  struct gkyl_array *fluid_rhou; // Momentum density vector (rho u_x, rho u_y, rho u_z).

  struct gkyl_array *g_dist_source, *F_k_p_1, *F_k_m_1; // Distribution function arrays for coupling different Laguerre moments.

  struct gkyl_array *cflrate_f; // CFL rate per cell.

  struct flash_lbo_collisions lbo;

  struct gkyl_wave_geom *ten_mom_geom; // Geometry object for Euler equation.
  struct gkyl_wv_eqn *eqn_ten_mom; // 10-moment equation object for PKPM moment calc.

  struct gkyl_dg_calc_pkpm_dist_vars *calc_pkpm_dist_vars; // Updater to compute div(p_par b_hat) and source terms

  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars, *calc_pkpm_vars_ext; // Updater to compute accel vars, pressure,
                                                                      // primitive variables, source terms, and integrated variables

  struct gkyl_dg_calc_flash_coupling *calc_flash_coupling; // Updater to compute flux contributions to Euler equations from
                                                           // kinetic equations: [div(P)-grad(p), (P - pI) : grad(u)]

  struct gkyl_dg_updater_moment *pkpm_mom_calc, *pkpm_mom_calc_diag; // Updaters for PKPM moment calculation.
  struct gkyl_dg_updater_flash *updater; // Updater for flow-frame kinetic solve.

  double *omegaCFL_ptr; // Pointer for `gkyl_array_reduce` operation when computing CFL frequency.

  // Boundary conditions on lower/upper edges in each direction.
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC on distribution functions and momentum.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  struct gkyl_bc_basic *bc_conf_lo[3];
  struct gkyl_bc_basic *bc_conf_up[3];

  // Buffer arrays used for applying BCs
  struct gkyl_array *bc_buffer, *bc_buffer_lo_fixed, *bc_buffer_up_fixed;
  bool bc_is_absorb;

  // To simplify BC application, store local skin and ghost ranges for distribution functions.
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  struct gkyl_range lower_skin_conf[GKYL_MAX_DIM];
  struct gkyl_range upper_skin_conf[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost_conf[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost_conf[GKYL_MAX_DIM];

  bool init_dist_from_arrays; // Flag to initialize distribution function from input fluid moment arrays.
  struct gkyl_array *init_dist_moms; // Input fluid moment arrays for initializing the distibution function with
                                     // the `init_dist_from_arrays` flag. Moments should be ordered as:
                                     // [rho, ux, uy, uz, temp].

  bool skip_init_dist; // Flag to skip initializing the distribution function.
};

// Electromagnetic field object.
struct flash_field {
  struct gkyl_flash_field info;
  double mu0, epsilon0;

  struct gkyl_array *em; // Electromagnetic field variables (Ex, Ey, Ez, Bx, By, Bz, phi, psi).
  struct gkyl_array *tot_em; // Total electromagnetic field.
  struct gkyl_array *bvar; // Magnetic field unit vector and tensor.
  struct gkyl_array *bvar_surf; // Surface expansion magnetic field unit vector and tensor (for use in pkpm model).
  struct gkyl_array *div_b; // Volume expansion of div(b) (for use in pkpm model).
  struct gkyl_array *max_b; // max(|b_i|) penalization (for use in pkpm model).
  struct gkyl_array *cell_avg_magB2; // Integer array for whether |B|^2 *only* uses cell averages for weak division.
                                     // Determined when constructing the matrix if |B|^2 < 0.0 at control points.

  bool has_ext_em; // Flag for external electromagnetic fields.
  struct gkyl_array *ext_em; // External electromagnetic field variables.
  gkyl_proj_on_basis *proj_ext_em; // Projection updater for external EM fields.

  struct gkyl_wave_geom *geom; // Geometry object for fields.
  struct gkyl_dg_calc_em_vars *calc_bvar; // Updater to compute magnetic field unit vector and tensor.

  // Boundary conditions on lower/upper edges in each direction.
  enum gkyl_field_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
};

// Top-level app object, containing species and fields.
struct gkyl_flash_app {
  char name[128]; // Name used for output file prefix.
  bool use_gpu;

  double cfl;

  struct gkyl_flash_stat stat;

  double tcurr; // Current time.

  int poly_order; // Polynomial order of kinetic species basis.
  int cdim, vdim; // Configuration and velocity space dimensionality.
  struct gkyl_rect_grid grid_kinetic, grid_fluid; // Config. space grids for Vlasov and fluid species.
  int nghost_fluid[GKYL_MAX_CDIM], nghost_kinetic[GKYL_MAX_CDIM]; // Number of config. space ghost cells.
  int num_periodic_dir; // Number of directions with periodic BCs.
  int periodic_dirs[GKYL_MAX_CDIM]; // Array of directions with periodic BCs.
  
  struct gkyl_range local_fluid, local_fluid_ext; // local, local-ext config. space ranges.
  struct gkyl_range global_fluid, global_fluid_ext; // global, global-ext config. space ranges.
  struct gkyl_range local_kinetic, local_kinetic_ext; // local, local-ext config. space ranges.
  struct gkyl_range global_kinetic, global_kinetic_ext; // global, global-ext config. space ranges.
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin_fluid[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost_fluid[GKYL_MAX_DIM];
  struct gkyl_range upper_skin_fluid[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost_fluid[GKYL_MAX_DIM];

  struct gkyl_rect_decomp *decomp_fluid, *decomp_kinetic; // Grid decomposition objects.
  struct gkyl_comm *comm_fluid, *comm_kinetic; // Communicator objects for configuration space.

  struct gkyl_basis basis, confBasis; // Config. and phase space basis sets.

  struct flash_field *field; // Pointer to field object.

  gkyl_flash_interp_vars *interp_vars; // Updater to interpolate between modal DG and finite volume representations.

  bool skip_fluid_species; // Flag to not initialize fluid species.
  
  int num_species;
  struct flash_species *species; // Pointer to species objects.
  struct flash_fluid_species *fluid_species; // Pointer to fluid species objects.
};

static inline void
bc_copy(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin,
  double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int c = 0; c < nc; ++c)
    ghost[c] = skin[c];
}

/**
 * Initialize kinetic species object. 
 * 
 * @param inp 
 * @param app 
 * @param spec 
 */
void
flash_species_init(struct gkyl_flash_inp *inp, 
  struct gkyl_flash_app *app, struct flash_species *spec);

/**
 * Interpolate finite volume momentum density rho * u onto kinetic grid.
 * 
 * @param app 
 * @param species 
 * @param fv_fluidin 
 */
void
flash_species_interp_rhou_to_dg(struct gkyl_flash_app *app,
  struct flash_species *species, const struct gkyl_array *fv_fluidin);

// Compute pkpm moments, [rho, p_par, p_perp].
void
flash_species_calc_pkpm_moms(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin);

/**
 * Compute PKPM model coupling moments (rho, p_par, p_perp), self-consistent
 * pressure force (div(p_par b_hat)), and primitive variables 
 * (ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m).
 * 
 * @param app 
 * @param species 
 * @param fin 
*/
void
flash_species_calc_pkpm_vars(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin);

/**
  * Compute PKPM update variables: acceleration variables in kinetic equation
  * and source distribution functions for Laguerre couplings.
  * 
  * @param app 
  * @param species 
  * @param fin 
*/
void
flash_species_calc_pkpm_update_vars(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin);

/**
 * Compute PKPM diagnostic moments, (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp).
 * 
 * @param app 
 * @param species 
 * @param fin 
 */
void
flash_species_calc_pkpm_diags(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin);

/**
 * Compute RHS of hybrid PKPM update, returning maximum stable timestep.
 * 
 * @param app 
 * @param species 
 * @param fin 
 * @param fluidin 
 * @param rhs_f 
 * @param rhs_fluid 
 * @return double 
 */
double flash_species_rhs(gkyl_flash_app *app, 
  struct flash_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs_f);

/**
 * Compute maximum stable timestep for kinetic equation without computing the RHS.
 * 
 * @param app 
 * @param s 
 * @return double 
 */
double
flash_species_compute_dt(gkyl_flash_app *app, struct flash_species *s);

/**
 * Compute L2 norm of the distribution function, storing result in species->integ_L2_f dynvec.
 * 
 * @param app 
 * @param tm 
 * @param species 
 */
void
flash_species_calc_L2(gkyl_flash_app *app, double tm, const struct flash_species *species);

/**
 * Apply initial conditions to PKPM distribution functions.
 * 
 * @param app 
 * @param t0 
 */
void
flash_species_apply_ic(gkyl_flash_app *app, struct flash_species *s, double t0);

/**
 * Free memory allocated for species object and associated members.
 * 
 * @param app
 * @param species 
 */
void
flash_species_release(gkyl_flash_app *app, struct flash_species *species);

/**
 * Apply boundary conditions to distribution function and fluid variables.
 * 
 * @param app 
 * @param species 
 * @param f 
 * @param fv_fluid 
 */
void
flash_species_apply_bc(gkyl_flash_app *app, const struct flash_species *species, 
  struct gkyl_array *f);

/**
 * Evaluate PKPM pressure at cell boundaries and center for use in 
 * computing contributions to Euler fluxes.
 * 
 * @param app  
 * @param species 
 * @param fin 
 * @param pkpm_flux 
 */
void
flash_species_calc_coupling_fluxes(struct gkyl_flash_app *app,
  struct flash_species *species, const struct gkyl_array *fin);


/**
 * Initialize a Maxwellian distribution using the given input moment array [n, ux, uy, uz, T/m]
 * 
 * @param app 
 * @param s 
 * @param init_dist_moms
 */
void
flash_species_calc_distf_from_moms(gkyl_flash_app *app, struct flash_species *s, 
  struct gkyl_array *init_dist_moms);

/**
 * Adjust the velocity space extents of the input species' internal phase space grid.
 * 
 * @param s 
 * @param lower 
 * @param upper 
 */
void
flash_species_set_vel_space_extents(gkyl_flash_app *app, struct flash_species *s, 
  double lower_vel[GKYL_MAX_VDIM], double upper_vel[GKYL_MAX_VDIM]);

/**
 * Initialize LBO collisions object.
 * 
 * @param app 
 * @param s 
 * @param lbo 
 */
void
flash_species_lbo_init(gkyl_flash_app *app, struct flash_species *s, struct flash_lbo_collisions *lbo);

/**
 * Compute required moments, primitive moments, and boundary corrections for LBO updater.
 * 
 * @param app 
 * @param species 
 * @param lbo 
 * @param fin 
 */
void
flash_species_lbo_moms(gkyl_flash_app *app, const struct flash_species *species,
  struct flash_lbo_collisions *lbo, const struct gkyl_array *fin);

/**
 * Compute RHS of PKPM LBO collision operator.
 * 
 * @param app 
 * @param species 
 * @param lbo 
 * @param fin 
 * @param rhs 
 */
void
flash_species_lbo_rhs(gkyl_flash_app *app, const struct flash_species *species, 
  struct flash_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Compute CFL frequency for LBO update without computing the RHS, storing it in species->cflrate_f.
 * 
 * @param app 
 * @param species 
 * @param lbo 
 */
void
flash_species_lbo_compute_cflrate(gkyl_flash_app *app, const struct flash_species *species,
  struct flash_lbo_collisions *lbo);

/**
 * Release memory for LBO updater.
 * 
 * @param app 
 * @param lbo 
 */
void 
flash_species_lbo_release(const struct flash_lbo_collisions *lbo);

/**
 * Initialize finite volume fluid species object. 
 * 
 * @param app 
 * @param spec 
 */
void
flash_fluid_species_init(struct gkyl_flash_app *app, struct flash_fluid_species *spec);

/**
 * Apply initial conditions to finite volume fluid quantities.
 * 
 * @param app 
 * @param t0 
 */
void
flash_fluid_species_apply_ic(gkyl_flash_app *app, struct flash_fluid_species *fs, double t0);

/**
 * Apply boundary conditions to distribution function and fluid variables.
 * 
 * @param app 
 * @param species 
 * @param f 
 * @param fv_fluid 
 */
void
flash_fluid_species_apply_bc(gkyl_flash_app *app, const struct flash_fluid_species *species, 
  struct gkyl_array *fluid);
  
/**
 * Advance fluid solution vector to next time step using wave propagation algorithm, 
 * and accumulate fluxes from kinetic closure.
 * 
 * @param app 
 * @param tcurr 
 * @param dt0 
 * @param s 
 * @param fs 
 * @param fluidin 
 * @param fluidout 
 * @return double 
 */
double
flash_fluid_species_update(gkyl_flash_app *app, double tcurr, double dt0,
  struct flash_species *s, struct flash_fluid_species *fs, const struct gkyl_array *fluidin,
  struct gkyl_array *fluidout);

/**
 * Compute primitive variables (rho, ux, uy, uz, p) from conservative fluid quantities.
 * 
 * @param app 
 * @param s 
 * @param fluid 
 * @param fluid_prim 
*/
void
flash_fluid_species_prim_vars(gkyl_flash_app *app, 
  struct flash_fluid_species *s, const struct gkyl_array *fluid, struct gkyl_array *fluid_prim);

/**
 * Free memory allocated for species object and associated members.
 * 
 * @param species 
 */
void
flash_fluid_species_release(struct flash_fluid_species *species);

/**
 * Create a new field object for mini-FLASH app. 
 * Used to compute magnetic field information for PKPM updates.
 * 
 * @param inp 
 * @param app 
 * @return struct flash_field* 
 */
struct flash_field*
flash_field_new(struct gkyl_flash_inp *inp, struct gkyl_flash_app *app);

/**
 * Initialize field object by using user-specified initial conditions
 * to compute electromagnetic fields.
 * 
 * @param app 
 * @param field 
 * @param t0 
 */
void
flash_field_apply_ic(struct gkyl_flash_app *app, struct flash_field *field, double t0);

/**
 * Compute magnetic field unit vector and unit tensor and div(b).
 * 
 * @param app 
 * @param field 
 * @param em 
 */
void
flash_field_calc_bvar(gkyl_flash_app *app, struct flash_field *field, struct gkyl_array *em);

/**
 * Release memory for field object and associated arrays and updaters.
 * 
 * @param f 
 */
void
flash_field_release(struct flash_field *f);

/**
 * Take a single forward Euler timestep for the Euler-PKPM hybrid system.
 * Timestep is never larger than the suggested timestep dt0, but may be smaller.
 * The true dt used in the update is returned in the `gkyl_update_status` struct.
 * 
 * @param app 
 * @param tcurr 
 * @param dt0 
 * @param fin 
 * @param fv_fluidin 
 * @param fout 
 * @param fv_fluidout 
 * @param st 
 */
void
flash_forward_euler(gkyl_flash_app *app, double tcurr, double dt0,
  const struct gkyl_array *fin[], const struct gkyl_array *fv_fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fv_fluidout[],
  struct gkyl_update_status *st);

/**
 * Take a single forward Euler timestep of ONLY the kinetic piece of the hybrid system.
 * Timestep is never larger than the suggested timestep dt0, but may be smaller.
 * The true dt used in the update is returned in the `gkyl_update_status` struct.
 * 
 * @param app 
 * @param dt0 
 * @param fin 
 * @param fout 
 * @param st 
*/
void
flash_vlasov_forward_euler(gkyl_flash_app *app, double dt0,
  const struct gkyl_array *fin[], struct gkyl_array *fout[],
  struct gkyl_update_status *st);

/**
 * Advance hybrid system by a single timestep using SSP-RK3 method.
 * Timestep is never larger than the suggested timestep dt, but may be smaller.
 * The true dt used in the update is returned in the `gkyl_update_status` struct.
 * 
 * @param app 
 * @param dt 
 * @return struct gkyl_update_status 
 */
struct gkyl_update_status
gkyl_flash_update_ssp_rk3(gkyl_flash_app *app, double dt);

/**
 * Advance hybrid system by a single timestep using two-stage RK method.
 * Timestep is never larger than the suggested timestep dt, but may be smaller.
 * The true dt used in the update is returned in the `gkyl_update_status` struct.
 * 
 * @param app 
 * @param dt 
 * @return struct gkyl_update_status 
 */
struct gkyl_update_status
gkyl_flash_update_ssp_rk2(gkyl_flash_app *app, double dt);