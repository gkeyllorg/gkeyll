#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_lbo_pkpm_nodal.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_lbo_pkpm_nodal gkyl_dg_updater_lbo_pkpm_nodal;

/**
 * Create a new updater for the Lenard-Bernstein (Dougherty) operator of the PKPM distributions
 * in the configuration-nodal, velocity-modal form: the equation object gkyl_dg_lbo_pkpm_nodal
 * with a hyper_dg solver over the velocity directions (zero flux at the ends of the velocity
 * grid).
 *
 * @param phase_grid Phase-space grid object.
 * @param conf_basis Configuration-space basis functions (p=1 serendipity).
 * @param phase_basis Phase-space basis functions (hybrid).
 * @param conf_range Configuration-space range.
 * @param auxin Auxiliary fields of the operator (see gkyl_dg_lbo_pkpm_nodal.h).
 * @param use_gpu Bool for whether the updater is on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_dg_updater_lbo_pkpm_nodal *gkyl_dg_updater_lbo_pkpm_nodal_new(
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  struct gkyl_dg_lbo_pkpm_nodal_auxfields *auxin, bool use_gpu
);

/**
 * Accumulate the operator onto rhs and its rate onto cflrate. The distributions and the
 * right-hand side are in the configuration-nodal, velocity-modal form.
 *
 * @param up Updater.
 * @param update_rng Phase-space range on which to compute the update.
 * @param fnodal Input distributions (8 blocks).
 * @param cflrate CFL rate in each cell, incremented.
 * @param rhs_nodal Right-hand side, incremented.
 */
void gkyl_dg_updater_lbo_pkpm_nodal_advance(
  struct gkyl_dg_updater_lbo_pkpm_nodal *up, const struct gkyl_range *update_rng,
  const struct gkyl_array *GKYL_RESTRICT fnodal, struct gkyl_array *GKYL_RESTRICT cflrate,
  struct gkyl_array *GKYL_RESTRICT rhs_nodal
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_updater_lbo_pkpm_nodal_release(struct gkyl_dg_updater_lbo_pkpm_nodal *up);
