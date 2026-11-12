#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_bgk_collisions gkyl_pkpm_bgk_collisions;

/**
 * Create an updater for the PKPM BGK operator at the 2^cdim interior Gauss-Legendre nodes of
 * configuration space: all 8 blocks relax to [f_t, f_t, 0, 0]. The distributions and the
 * increments are in the configuration-nodal, velocity-modal form (gkyl_dg_conf_nodal).
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param phase_basis Phase-space basis (p=1 x p=2 hybrid).
 * @param use_gpu Whether to run on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_pkpm_bgk_collisions *gkyl_pkpm_bgk_collisions_new(
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu
);

/**
 * Accumulate the explicit BGK term onto the nodal right-hand side, rate_q (J f_t - f_q) at each
 * node, and the largest node rate of each cell onto cflrate.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param vel_range Velocity-space range indexing f_lte.
 * @param phase_range Phase-space range to update.
 * @param rate Relaxation rate at the nodes.
 * @param jacob Position-map Jacobian (first component, constant in a cell).
 * @param f_lte Target: 3 coefficients per velocity cell.
 * @param fnodal Input distribution in the nodal form (8 blocks).
 * @param rhs_nodal Right-hand side in the nodal form to accumulate onto (8 blocks).
 * @param cflrate CFL rate to accumulate onto.
 */
void gkyl_pkpm_bgk_collisions_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *rate, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  const struct gkyl_array *fnodal, struct gkyl_array *rhs_nodal, struct gkyl_array *cflrate
);

/**
 * Implicit relaxation of the distributions, in place: f <- S f + (1 - S) J f_target, with the
 * shape factor S given at the nodes (S = 1/(1 + nu dt) for a backward-Euler step), F and G
 * relaxing to the Maxwellian and the first-harmonic blocks to zero. It is the weak product of the
 * expansion of S with each block of f; at p=1 that equals the product at the nodes, so no nodal
 * copy of f is needed.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param vel_range Velocity-space range indexing f_lte.
 * @param phase_range Phase-space range to update.
 * @param shape Shape factor S at the nodes.
 * @param jacob Position-map Jacobian (first component, constant in a cell).
 * @param f_lte Target: 3 coefficients per velocity cell.
 * @param f Distributions (8 blocks of expansions), updated in place.
 */
void gkyl_pkpm_bgk_collisions_implicit_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *shape, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  struct gkyl_array *f
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_bgk_collisions_release(struct gkyl_pkpm_bgk_collisions *up);
