#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_bgk_collisions gkyl_pkpm_bgk_collisions;

/**
 * Create an updater for the PKPM BGK operator: all 8 blocks relax to [f_t, f_t, 0, 0], with a
 * rate (or an implicit shape factor) given at the 2^cdim interior Gauss-Legendre nodes of
 * configuration space. The target f_t is one velocity-space function, the unit Maxwellian in w,
 * the same in every cell, so the operator is a factor times (J f_t - f): it is applied as the
 * weak product of the expansion of the factor with each block of the distributions. At p=1 that
 * equals the product at the nodes, and no nodal copy of the distributions is needed.
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
 * Accumulate the explicit BGK term onto the right-hand side, rate (J f_t - f), and the largest
 * node rate of each cell onto cflrate.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param vel_range Velocity-space range indexing f_lte.
 * @param phase_range Phase-space range to update.
 * @param rate Relaxation rate at the nodes.
 * @param jacob Position-map Jacobian (first component, constant in a cell).
 * @param f_lte Target: 3 coefficients per velocity cell.
 * @param f Distributions (8 blocks of expansions).
 * @param rhs Right-hand side to accumulate onto (8 blocks of expansions).
 * @param cflrate CFL rate to accumulate onto.
 */
void gkyl_pkpm_bgk_collisions_advance(
  const struct gkyl_pkpm_bgk_collisions *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *rate, const struct gkyl_array *jacob, const struct gkyl_array *f_lte,
  const struct gkyl_array *f, struct gkyl_array *rhs, struct gkyl_array *cflrate
);

/**
 * Implicit relaxation of the distributions, in place: f <- S f + (1 - S) J f_target, with the
 * shape factor S given at the nodes (S = 1/(1 + nu dt) for a backward-Euler step), F and G
 * relaxing to the Maxwellian and the first-harmonic blocks to zero.
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
