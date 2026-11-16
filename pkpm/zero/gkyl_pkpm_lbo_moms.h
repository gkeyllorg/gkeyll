#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_lbo_moms gkyl_pkpm_lbo_moms;

// Input packaged as a struct
struct gkyl_pkpm_lbo_moms_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 x p=2 hybrid).
  const double *vel_bounds; // Lower and upper bound of the velocity grid, [w_min, w_max].
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create an updater for the primitive moments [u, vtsq] of the PKPM Lenard-Bernstein operator,
 * the drift and the thermal speed squared of the normalized distributions in the w coordinate
 * (not the [n, u (x, y, z), vtsq] of gkyl_pkpm_layout.h), at the 2^cdim interior Gauss-Legendre
 * nodes of configuration space. At each node they solve
 *
 *   M0 u - cM vtsq = M1
 *   M1 u + (3 M0 - cE) vtsq = M2 + 2 M0_G
 *
 * with cM = F(w_max) - F(w_min) and cE = w_max F(w_max) - w_min F(w_min), the corrections for
 * the finite velocity extent, so that the operator has no [M0, M1, M2 + 2 M0_G] moments. At
 * p=1 the node values equal those of the Galerkin solve of the same system.
 *
 * The moments and the distributions may carry the position-map Jacobian, which is constant in a
 * cell: the primitive moments are homogeneous in it.
 *
 * The determinant, M0 (3 M0 - cE) + cM M1, is positive as long as F is small at the edges of the
 * velocity grid; the solve has no guard.
 *
 * @param inp Input parameters defined in gkyl_pkpm_lbo_moms_inp struct.
 * @return New updater pointer.
 */
struct gkyl_pkpm_lbo_moms *gkyl_pkpm_lbo_moms_inew(const struct gkyl_pkpm_lbo_moms_inp *inp);

/**
 * Compute the primitive moments at the nodes.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range (its velocity extent locates the edge cells).
 * @param moms Moments [M0, M2, M0_G, M1, ...] of the normalized distributions (p=1 expansions).
 * @param fnodal Distribution functions in the configuration-nodal, velocity-modal form
 *   (gkyl_dg_conf_nodal; the leading block is F).
 * @param prim_moms Output primitive moments [u, vtsq] in w at the nodes (one value per node each).
 */
void gkyl_pkpm_lbo_moms_advance(
  const struct gkyl_pkpm_lbo_moms *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *moms,
  const struct gkyl_array *fnodal, struct gkyl_array *prim_moms
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_lbo_moms_release(struct gkyl_pkpm_lbo_moms *up);
