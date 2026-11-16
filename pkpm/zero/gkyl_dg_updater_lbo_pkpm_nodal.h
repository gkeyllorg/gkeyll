#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_lbo_pkpm_nodal.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_lbo_pkpm_nodal gkyl_dg_updater_lbo_pkpm_nodal;

// Input packaged as a struct
struct gkyl_dg_updater_lbo_pkpm_nodal_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 x p=2 hybrid).
  const struct gkyl_range *conf_range; // Configuration-space range for indexing nu and prim_moms.
  const struct gkyl_array *nu; // Collision frequency at the configuration-space nodes.
  const struct gkyl_array *prim_moms; // Nodal drift and thermal speed squared [u, vtsq] in w units.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater for the Lenard-Bernstein (Dougherty) operator of the PKPM distributions
 * in the configuration-nodal, velocity-modal form: the equation object gkyl_dg_lbo_pkpm_nodal
 * with a hyper_dg solver over the velocity directions (zero flux at the ends of the velocity
 * grid).
 *
 * @param inp Input parameters defined in gkyl_dg_updater_lbo_pkpm_nodal_inp struct.
 * @return New updater pointer.
 */
struct gkyl_dg_updater_lbo_pkpm_nodal *gkyl_dg_updater_lbo_pkpm_nodal_inew(
  const struct gkyl_dg_updater_lbo_pkpm_nodal_inp *inp
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
