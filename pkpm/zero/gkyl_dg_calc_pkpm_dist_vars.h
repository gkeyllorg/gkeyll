#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_dist_vars gkyl_dg_calc_pkpm_dist_vars;

/**
 * Create new updater to compute pkpm distribution function variables needed in 
 * updates and used for diagnostics. 
 * 
 * Updater stores the kernels to compute 
 * 1. pkpm distribution functions for the mirror force and source terms
 * 2. pkpm div(p_par b) from the distribution function for self-consistency 
 *    and to avoid spurious development of first moment of kinetic equation.
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param cbasis Configuration space basis functions
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_dist_vars* 
gkyl_dg_calc_pkpm_dist_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis* cbasis, bool use_gpu);

/**
 * Create new updater to compute pkpm variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_pkpm_dist_vars* 
gkyl_dg_calc_pkpm_dist_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis* cbasis);

/**
 * Compute parallel-kinetic-perpendicular-moment model distribution function
 * in the mirror force for the T_perp/m*G = T_perp/m*(F_0 - F_1) kinetic equation,
 * along with the vperp characteristics which are a pure source term in the first Laguerre moment update.
 * g_dist_source = [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
 *                  (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]
 * First output is mirror force source *distribution*, second output is *total* vperp characteristics source.
 * Note that T_perp/m*G is the evolved quantity for the first Laguerre moment. 
 * Also outputs F_1 from T_perp/m*G for the evolution of F_2 if F_2 is present. 
 * To simplify internal Gkeyll logic, assume F_2 is present and output F_1 even if F_2 = 0.0.
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Range to index configuration space fields
 * @param phase_range Range to calculate mirror force distribution function 
 * @param pkpm_prim Input array of primitive moments [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp]
 *        This updater needs T_perp/m and m/T_perp
 * @param nu_prim_moms_sum Input array of nu*primitive moments, summed over self and cross collisions (needed for nu*vth^2)
 * @param div_b Input array of div(b)
 * @param pkpm_accel_vars Input arrary of pkpm acceleration variables ordered as:
 *        0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
          This updater p_perp_source.
 * @param fIn Input array of pkpm distribution functions: [F_0, T_perp/m G = T_perp/m (F_0 - F_1)]
 * @param F_k_p_1 Input array of k+1 distribution function. F_2 expansion is the first NP coefficients.
 * @param g_dist_source Output array: 2.0*(T_perp/m*G) + T_perp/m*(F_2 - F_0)
 * @param F_k_m_1 Output array of k-1 distribution function. F_1 expansion is the first NP coefficients.
 */
void gkyl_dg_calc_pkpm_dist_vars_mirror_force(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* pkpm_prim, const struct gkyl_array* nu_prim_moms_sum, 
  const struct gkyl_array* div_b, const struct gkyl_array* pkpm_accel_vars, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1);

/**
 * Compute divergence of parallel pressure for use in pressure force in kinetic equation.
 * Needed for self-consistency and to avoid spurious development of first moment of kinetic equation.
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Range to index configuration space fields
 * @param phase_range Range to index phase space fields
 * @param bvar_surf Input array of surface expansion of magnetic field unit vector and unit tensor
 * @param bvar Input array of volume expansion of magnetic field unit vector and unit tensor
 * @param fIn Input array of pkpm distribution functions: [F_0, T_perp/m G = T_perp/m (F_0 - F_1)]
 * @param max_b Input array of surface expansion of max(|b_i|) penalization
 * @param pkpm_div_ppar Output array of divergence of p_parallel b_hat
 */
void gkyl_dg_calc_pkpm_dist_vars_div_ppar(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, const struct gkyl_array* fIn, 
  const struct gkyl_array* max_b, struct gkyl_array* pkpm_div_ppar);

/**
 * Delete pointer to updater to compute pkpm variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_dist_vars_release(struct gkyl_dg_calc_pkpm_dist_vars *up);

/**
 * Host-side wrappers for pkpm dist vars operations on device
 */

void gkyl_dg_calc_pkpm_dist_vars_mirror_force_cu(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* pkpm_prim, const struct gkyl_array* nu_prim_moms_sum, 
  const struct gkyl_array* div_b, const struct gkyl_array* pkpm_accel_vars, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1);

void gkyl_dg_calc_pkpm_dist_vars_div_ppar_cu(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, const struct gkyl_array* fIn, 
  const struct gkyl_array* max_b, struct gkyl_array* pkpm_div_ppar);

