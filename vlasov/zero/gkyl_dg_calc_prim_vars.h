#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>

/**
 * Compute u from state vector which contains *both* density and momentum.
 * Data always organized such that density is 0th component of statevec
 * & momentum is 1st, 2nd, and 3rd component of statevec.
 * Examples include: 
 * Isothermal Euler state vector [rho, rhoux, rhouy, rhouz],
 * Euler state vector [rho, rhoux, rhouy, rhouz, E],
 * Nonrelativistic Vlasov Five Moments array [M0, M1i, M2] (u_i = M1i/M0),
 * Special Relativistic Vlasov N_i array [Gamma*n, Gamma*n*u] where u is the bulk flow
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param statevec Input state vector which contains *both* density and momentum
 * @param u_i Output array of bulk flow velocity
 */
void gkyl_calc_prim_vars_u_from_statevec(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* statevec, struct gkyl_array* u_i);

/**
 * Compute u from input density and momentum vectors.
 * Examples include: 
 * Nonrelativistic Vlasov M0 & M1i (u_i = M1i/M0),
 * Special Relativistic M0 (Gamma*n) and M1i (Gamma*n*u) where u is the bulk flow
 * Parallel-kinetic-perpendicular-moment model with vlasov_pkpm_moms and euler_pkpm
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param rho Input density 
 * @param rhou Input momentum
 * @param u_i Output array of bulk flow velocity
 */
void gkyl_calc_prim_vars_u_from_rhou(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* rho, const struct gkyl_array* rhou, struct gkyl_array* u_i);
