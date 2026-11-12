#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_temp_source gkyl_pkpm_temp_source;

/**
 * Create an updater for the fluid part of the PKPM temperature source
 * S = nu_Q (f_M(n, u, T_Q) - f), at the 2^cdim interior Gauss-Legendre nodes of configuration
 * space. The density and the momentum are unchanged. The thermal energy
 * eps = E - |M|^2/(2 rho) relaxes to eps_Q = (3/2) rho vtsq_Q at the rate nu_Q, and the
 * normalized distributions relax to the unit Maxwellian at the rate nu_Q eps_Q/eps, which the
 * updater returns for the BGK updater (gkyl_pkpm_bgk_collisions).
 *
 * The fluid state may carry the position-map Jacobian, which is constant in a cell: every
 * quantity formed here is homogeneous in it.
 *
 * At a node with a nonpositive density the energy is left unchanged, and at a node with a
 * nonpositive density or thermal energy the distributions relax at the rate nu_Q.
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param use_gpu Whether to run on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_pkpm_temp_source *gkyl_pkpm_temp_source_new(
  const struct gkyl_basis *conf_basis, bool use_gpu
);

/**
 * Explicit update: accumulate nu_Q (eps_Q - eps) at the nodes onto the energy block of rhs and
 * the largest node rate nu_Q of each cell onto cflrate, and compute the relaxation rate of the
 * normalized distributions, nu_Q eps_Q/eps, at the nodes.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param rate Source rate nu_Q at the nodes.
 * @param vtsq Source temperature over the mass, T_Q/m, at the nodes.
 * @param fluid Fluid state [rho, M (3), E] (p=1 expansions).
 * @param rhs Right-hand side of the fluid state to accumulate onto.
 * @param cflrate CFL rate of the fluid state to accumulate onto.
 * @param shape_rate Output relaxation rate of the normalized distributions at the nodes.
 */
void gkyl_pkpm_temp_source_advance(
  const struct gkyl_pkpm_temp_source *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rate, const struct gkyl_array *vtsq, const struct gkyl_array *fluid,
  struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_array *shape_rate
);

/**
 * One backward-Euler step of the thermal energy in place,
 * eps <- (eps + nu_Q dt eps_Q)/(1 + nu_Q dt) at each node, and the shape factor of the same step
 * for the normalized distributions, S = 1/(1 + nu_Q dt eps_Q/eps), with eps the thermal energy
 * before the step. Ghost cells are not updated.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param dt Time step.
 * @param rate Source rate nu_Q at the nodes.
 * @param vtsq Source temperature over the mass, T_Q/m, at the nodes.
 * @param fluid Fluid state [rho, M (3), E] (p=1 expansions), updated in place.
 * @param shape Output shape factor S at the nodes.
 * @param temp_ratio Output: thermal energy (temperature) before the step over that after it, at
 *   the nodes.
 */
void gkyl_pkpm_temp_source_implicit_advance(
  const struct gkyl_pkpm_temp_source *up, const struct gkyl_range *conf_range, double dt,
  const struct gkyl_array *rate, const struct gkyl_array *vtsq, struct gkyl_array *fluid,
  struct gkyl_array *shape, struct gkyl_array *temp_ratio
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_temp_source_release(struct gkyl_pkpm_temp_source *up);
