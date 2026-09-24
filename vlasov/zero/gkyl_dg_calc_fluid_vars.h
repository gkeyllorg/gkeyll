#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_dg_calc_fluid_vars gkyl_dg_calc_fluid_vars;

/**
 * Create new updater to compute fluid variables needed in 
 * update of DG fluid equations (e.g., isothermal Euler, Euler, 10 moment). Methods compute:
 * p : Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2)
 * p_surf : [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr]
 * u : [ux, uy, uz]
 * u_surf : [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
 *           ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
 *           ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
 * 
 * Updater also stores the kernels to compute fluid source terms and fluid integrated moments.
 * 
 * @param wv_eqn      Wave equation (stores function pointers for computing waves and limiting solution)
 * @param geom        Wave geometry object for computing waves in local coordinate system
 * @param cbasis      Configuration space basis functions
 * @param mem_range   Configuration space range that sets the size of the bin_op memory
 *                    for computing primitive moments. Note range is stored so 
 *                    updater loops over consistent range for primitive moments
 * @param limiter_fac Optional parameter for changing diffusion in sloper limiter 
 *                    by changing relationship between slopes and cell average differences.
 *                    By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
 *                    and a cell slope estimate from two adjacent cells is (for the x variation): 
 *                    integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
 *                    where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
 *                    This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
 * @param use_gpu     bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_fluid_vars* 
gkyl_dg_calc_fluid_vars_new(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  double limiter_fac, bool use_gpu);

/**
 * Create new updater to compute fluid variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_fluid_vars* 
gkyl_dg_calc_fluid_vars_cu_dev_new(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  double limiter_fac);

/**
 * Compute flow velocity from mass density and momentum density.
 *
 * @param up            Updater for computing fluid variables 
 * @param fluid         Input array of fluid variables [rho, rho ux, rho uy, rho uz, ...]
 * @param cell_avg_prim Array for storing boolean value of whether rho is negative at corners 
 *                      Note: Only used for diagnostic purposes (not for adjusting solution)
 * @param u             Output array of volume expansion of flow velocity [ux, uy, uz]
 * @param u_surf        Output array of surface expansion of flow velocity
 *                      [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
 *                       ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
 *                       ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
 */
void gkyl_dg_calc_fluid_vars_advance(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_array* fluid, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* u, struct gkyl_array* u_surf);

/**
 * Compute pressure from fluid variables in the volume and at needed surfaces
 *
 * @param up Updater for computing fluid variables 
 * @param conf_range Configuration space range
 * @param fluid      Input array of fluid variables [rho, rho ux, rho uy, rho uz, ...]
 * @param u          Input array of volume expansion of flow velocity [ux, uy, uz]
 * @param p          Output array of volume expansion of pressure 
 * @param p_surf     Output array of surface expansion of pressure [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 
 */
void gkyl_dg_calc_fluid_vars_pressure(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf);

/**
 * Compute kinetic energy from fluid variables in the volume 
 *
 * @param up Updater for computing fluid variables 
 * @param conf_range Configuration space range
 * @param fluid      Input array of fluid variables [rho, rho ux, rho uy, rho uz, ...]
 * @param u          Input array of volume expansion of flow velocity [ux, uy, uz]
 * @param ke         Output array of volume expansion of kinetic energy
 */
void gkyl_dg_calc_fluid_vars_ke(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* ke);

/**
 * Limit slopes for fluid variables
 *
 * @param up         Updater for computing fluid variables 
 * @param conf_range Configuration space range
 * @param fluid      Input (and Output after limiting) array of fluid variables [rho, rho ux, rho uy, rho uz, ...]
 */
void gkyl_dg_calc_fluid_vars_limiter(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, struct gkyl_array* fluid);

/**
 * Compute integrated fluid variables (rho, rhoux, rhouy, rhouz, rhou^2, ...).
 * For isothermal Euler, integrated internal energy is just weighted integrated density (rho*vth)
 * For Euler, computes integrated internal energy p (without the 1/(gas_gamma - 1) factor) 
 * For 10 moment, computes integrated internal energy from sum over directions (Pxx + Pyy + Pzz)
 * Note: post-processing requires the addition of the needed factors (such as 1/2 in rhou^2)
 *
 * @param up Updater for computing fluid variables 
 * @param conf_range Configuration space range
 * @param fluid Input array of fluid variables [rho, rhoux, rhouy, rhouz, ...]
 * @param u_i Input array of flow velocity [ux, uy, uz]
 * @param p_ij Input array of pressure 
 * @param int_fluid_vars Output array of integrated variables (6 components)
 */
void gkyl_dg_calc_fluid_integrated_vars(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* fluid, 
  const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* fluid_int_vars);

/**
 * Compute fluid model source terms.
 *
 * @param up         Updater for computing fluid variables 
 * @param conf_range Configuration space range
 * @param app_accel  Input array of applied acceleration (external forces)
 * @param fluid      Input array of fluid variables [rho, rhoux, rhouy, rhouz, ...]
 * @param rhs        Output increment to fluid variables
 */
void gkyl_dg_calc_fluid_vars_source(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* app_accel, const struct gkyl_array* fluid, 
  struct gkyl_array* rhs);

/**
 * Delete pointer to updater to compute fluid variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_fluid_vars_release(struct gkyl_dg_calc_fluid_vars *up);

/**
 * Host-side wrappers for fluid vars operations on device
 */

void gkyl_dg_calc_fluid_vars_advance_cu(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_array* fluid, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* u, struct gkyl_array* u_surf);

void gkyl_dg_calc_fluid_vars_pressure_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf);

void gkyl_dg_calc_fluid_vars_ke_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* ke);

void gkyl_dg_calc_fluid_vars_limiter_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, struct gkyl_array* fluid);

void gkyl_dg_calc_fluid_integrated_vars_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* fluid, 
  const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* fluid_int_vars);

void gkyl_dg_calc_fluid_vars_source_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* app_accel, const struct gkyl_array* fluid, 
  struct gkyl_array* rhs);
