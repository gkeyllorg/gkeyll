#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_stage_weights gkyl_pkpm_stage_weights;

/**
 * Create an updater for the stage combination of the normalized PKPM distributions of a species
 * whose implicit operators change its temperature inside the SSP-RK3 stages (implicit BGK with
 * cross-species collisions, implicit temperature source).
 *
 * The stages combine the state M at the start of the step with a forward-Euler output Y,
 * c_M M + c_Y Y. A physical second moment is the temperature times its normalized value, so if an
 * implicit step has multiplied the temperature of the state behind Y by lam (lam = 1 for M), the
 * combination that is linear in the physical second moments is
 *
 *   w_M M + w_Y Y,  w_Y = c_Y lam / (c_M + c_Y lam),  w_M = 1 - w_Y,
 *
 * and the temperature factor of the combined state is lam_c = c_M + c_Y lam. The weights are
 * non-negative and sum to one at every node of configuration space, so the combination is convex
 * and keeps the constraints; with lam = 1 it is the SSP-RK3 combination. lam is given at the
 * 2^cdim interior Gauss-Legendre nodes, and the weights multiply the distributions as weak
 * products (equal at p=1 to the products at the nodes).
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param phase_basis Phase-space basis (p=1 x p=2 hybrid).
 * @param use_gpu Whether to run on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_pkpm_stage_weights *gkyl_pkpm_stage_weights_new(
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu
);

/**
 * Turn the SSP-RK3 combination f = c_M f_old + c_Y f_fe into the weighted one, in place:
 * f += (w_Y - c_Y) (f_fe - f_old), for every block of the distributions, and replace lam by the
 * temperature factor of the combined state, c_M + c_Y lam.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range to update.
 * @param coef_fe SSP-RK3 weight c_Y of the forward-Euler output (c_M = 1 - c_Y).
 * @param f_old State at the start of the step.
 * @param f_fe Forward-Euler output of the stage.
 * @param lam Temperature factor of the state behind f_fe at the nodes; on output, that of the
 *   combined state.
 * @param f The SSP-RK3 combination on input, the weighted combination on output.
 */
void gkyl_pkpm_stage_weights_advance(
  const struct gkyl_pkpm_stage_weights *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, double coef_fe, const struct gkyl_array *f_old,
  const struct gkyl_array *f_fe, struct gkyl_array *lam, struct gkyl_array *f
);

/**
 * Account for an implicit step in the temperature factor: lam <- lam / temp_ratio at the nodes,
 * temp_ratio being the temperature before the step over that after it.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param temp_ratio Temperature before the implicit step over that after it, at the nodes.
 * @param lam Temperature factor at the nodes, updated in place.
 */
void gkyl_pkpm_stage_weights_scale(
  const struct gkyl_pkpm_stage_weights *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *temp_ratio, struct gkyl_array *lam
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_stage_weights_release(struct gkyl_pkpm_stage_weights *up);
