#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Object type
typedef struct gkyl_dg_calc_sr_vars gkyl_dg_calc_sr_vars;

/**
 * Create new updater to compute relativistic variables needed in 
 * updates and used for diagnostics. Methods compute:
 * "p_vars" : the particle Lorentz boost factor gamma = sqrt(1 + p^2) and its inverse
 * n : the rest-frame plasma density n = GammaV_inv*M0 where GammaV_inv = sqrt(1 - |V_drift|^2)
 *     and V_drift is the bulk velocity computed via weak division V_drift = M1i/M0
 * u_i : the bulk four-velocity (GammaV, GammaV*V_drift) computed using weak division from M0 and M1i
 * p = n*T : the rest-frame pressure. The rest-frame pressure is computed as a velocity moment 
 *           with the weight = gamma*GammaV^2 - 2*GammaV*(v . p) + 1/gamma*((v . p)^2 - 1)
 *           where v is the spatial component of the bulk four-velocity: GammaV*V_drift, 
 *           GammaV is the bulk Lorentz boost factor: sqrt(1 + v^2), 
 *           p is the spatial component of the particle four-velocity, 
 *           and gamma = sqrt(1 + p^2) is the particle Lorentz boost factor.
 * 
 * @param phase_grid Phase-space grid (for getting cell spacing and cell center in pressure calculation) 
 * @param vel_grid   Momentum (four-velocity)-space grid 
 *                   (for getting cell spacing and cell center in gamma and 1/gamma calculation) 
 * @param conf_basis Configuration-space basis functions
 * @param vel_basis  Momentum (four-velocity)-space basis functions
 * @param mem_range  Configuration-space range that sets the size of the bin_op memory
 *                   for computing V_drift. Note range is stored so updater loops 
 *                   over consistent range solving linear systems since memory is pre-allocated.
 * @param vel_range  Momentum (four-velocity)-space
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_sr_vars* 
gkyl_dg_calc_sr_vars_new(const struct gkyl_rect_grid *phase_grid, const struct gkyl_rect_grid *vel_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *vel_basis, 
  const struct gkyl_range *mem_range, const struct gkyl_range *vel_range, bool use_gpu);

/**
 * Create new updater to compute relativistic variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_sr_vars* 
gkyl_dg_calc_sr_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, const struct gkyl_rect_grid *vel_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *vel_basis, 
  const struct gkyl_range *mem_range, const struct gkyl_range *vel_range);

/**
 * Compute the momentum grid variables for special relativistic simulations
 * Uses special kernels which convert between a Gauss-Lobatto nodal basis and
 * our modal basis to insure continuity of the momentum (four-velocity)-grid variables.
 *
 * @param up        Updater for computing sr variables 
 * @param gamma     Output array of particle Lorentz boost factor, gamma = sqrt(1 + p^2) 
 * @param gamma_inv Output array of inverse particle Lorentz boost factor, 1/gamma = 1/sqrt(1 + p^2) 
 */
void gkyl_calc_sr_vars_init_p_vars(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv);

/**
 * Compute the rest-frame density n = GammaV_inv*M0 where GammaV_inv = sqrt(1 - |V_drift|^2).
 * Updater first computes V_drift from M1i = GammaV*n*V_drift and M0 = GammaV*n using weak division. 
 * The updater then attempts to do the following operation
 * 1. Compute (GammaV_inv)^2 = 1 - |V_drift|^2 using basis_exp_sq 
 *    (see gkyl_basis_*_exp_sq.h in kernels/basis/)
 * 2. Project onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt (see gkyl_basis_*_sqrt.h in kernels/basis/)
 * If any quadrature points are assessed to be negative, the updater then switches to a more robust
 * evaluation of V_drift at p=1 Gauss-Lobatto quadrature points (the corners of the reference element)
 * and then evaluates 1 - (V_drift_lobatto)^2 (with V_drift_lobatto = 1.0e-16 if the result is still negative).
 * We can then guarantee the 1 - (V_drift_lobatto)^2 > 0.0 at quadrature points and GammaV_inv is well defined.
 * Finally, the rest-frame density is computed with weak multiplication with the final GammaV_inv basis expansion.
 *
 * @param up  Updater for computing sr variables 
 * @param M0  Input lab-frame density = GammaV*n
 * @param M1i Input lab-frame flux = GammaV*n*V_drift
 * @param n   Output rest-frame density.
 */
void gkyl_dg_calc_sr_vars_n(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_array* M0, const struct gkyl_array* M1i, struct gkyl_array* n);

/**
 * Compute derived quantities from spatial component of the bulk four-velocity u_i = GammaV*V_drift.
 * These include the components squared u_i^2, the Lorentz boost factor GammaV = sqrt(1 + |u_i|^2)
 * and the square of the Lorentz boost factor. 
 *
 * @param up         Updater for computing sr variables 
 * @param conf_range Configuration-space range
 * @param u_i        Input spatial component of the bulk four-velocity GammaV*V_drift
 * @param u_i_sq     Output squared spatial components of the bulk four-velocity 
 * @param GammaV     Output Lorentz boost factor for the bulk four-velocity sqrt(1 + |u_i|^2)
 * @param GammaV_sq  Output square of the Lorentz boost factor for the bulk four-velocity
 */
void gkyl_dg_calc_sr_vars_GammaV(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range,
  const struct gkyl_array* u_i, struct gkyl_array* u_i_sq, 
  struct gkyl_array* GammaV, struct gkyl_array* GammaV_sq);

/**
 * Compute the rest-frame pressure = n*T. The rest-frame pressure is computed as a velocity moment.
 * The weight = gamma*GammaV^2 - 2*GammaV*(v . p) + 1/gamma*((v . p)^2 - 1)
 * where v is the spatial component of the bulk four-velocity: GammaV*V_drift, 
 * GammaV is the bulk Lorentz boost factor: sqrt(1 + v^2), 
 * p is the spatial component of the particle four-velocity, 
 * and gamma = sqrt(1 + p^2) is the particle Lorentz boost factor.
 *
 * @param up          Updater for computing sr variables 
 * @param conf_range  Configuration-space range
 * @param vel_range   Momentum (four-velocity)-space range
 * @param phase_range Phase-space range
 * @param gamma       Input array of particle Lorentz boost factor, gamma = sqrt(1 + p^2) 
 * @param gamma_inv   Input array of inverse particle Lorentz boost factor, 1/gamma = 1/sqrt(1 + p^2) 
 * @param u_i         Input spatial component of the bulk four-velocity GammaV*V_drift
 * @param u_i_sq      Input squared spatial components of the bulk four-velocity 
 * @param GammaV      Input Lorentz boost factor for the bulk four-velocity sqrt(1 + |u_i|^2)
 * @param GammaV_sq   Input square of the Lorentz boost factor for the bulk four-velocity
 * @param f           Input distribution function
 * @param sr_pressure Output pressure
 */
void gkyl_dg_calc_sr_vars_pressure(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* gamma, const struct gkyl_array* gamma_inv, 
  const struct gkyl_array* u_i, const struct gkyl_array* u_i_sq, 
  const struct gkyl_array* GammaV, const struct gkyl_array* GammaV_sq, 
  const struct gkyl_array* f, struct gkyl_array* sr_pressure);

/**
 * Delete pointer to updater to compute sr variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_sr_vars_release(struct gkyl_dg_calc_sr_vars *up);

/**
 * Host-side wrappers for sr vars operations on device
 */

void gkyl_calc_sr_vars_init_p_vars_cu(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv);

void gkyl_dg_calc_sr_vars_n_cu(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_array* M0, const struct gkyl_array* M1i, struct gkyl_array* n);

void gkyl_dg_calc_sr_vars_GammaV_cu(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range,
  const struct gkyl_array* u_i, struct gkyl_array* u_i_sq, 
  struct gkyl_array* GammaV, struct gkyl_array* GammaV_sq);

void gkyl_dg_calc_sr_vars_pressure_cu(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* gamma, const struct gkyl_array* gamma_inv, 
  const struct gkyl_array* u_i, const struct gkyl_array* u_i_sq, 
  const struct gkyl_array* GammaV, const struct gkyl_array* GammaV_sq, 
  const struct gkyl_array* f, struct gkyl_array* sr_pressure);
