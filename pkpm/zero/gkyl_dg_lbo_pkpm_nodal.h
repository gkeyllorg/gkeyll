#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_pkpm_nodal_auxfields {
  const struct gkyl_array *nu; // collision frequency at the configuration-space nodes
  const struct gkyl_array *prim_moms; // primitive moments [u, vtsq] at the nodes
};

/**
 * Create a new equation object for the Lenard-Bernstein (Dougherty) operator of the PKPM
 * distributions in the configuration-nodal, velocity-modal form (gkyl_dg_conf_nodal):
 *
 *   d_w[(nu w - nu u) h + nu vtsq d_w h]
 *
 * on each of the 8 blocks h = [F, G, M0_perp (3), N_perp (3)], with zero flux at the ends of the
 * velocity grid, plus the couplings 2 (nu vtsq F - nu G) on G, -nu M0_perp on M0_perp and
 * 4 nu vtsq M0_perp - 3 nu N_perp on N_perp. Drag and diffusion are one equation. Its kernels
 * act in the velocity alone and are applied at each node of a configuration-space cell with the
 * node's nu, u and vtsq, which at p=1 equals the Galerkin form of the operator. A node whose
 * vtsq is not between 0 and the square of the largest velocity of the grid is skipped.
 *
 * The CFL rate of a phase-space cell is the largest over its nodes of
 * 2.5 (2/dv) nu |u - w_c| + 3 nu + 9 (2/dv)^2 nu vtsq.
 *
 * @param cbasis Configuration-space basis functions (p=1 serendipity).
 * @param pbasis Phase-space basis functions (hybrid).
 * @param conf_range Configuration-space range for use in indexing the auxiliary fields.
 * @param pgrid Phase-space grid object.
 * @param use_gpu Bool to determine if equation object is on host or device (not available yet).
 * @return Pointer to the equation object.
 */
struct gkyl_dg_eqn *gkyl_dg_lbo_pkpm_nodal_new(
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_rect_grid *pgrid, bool use_gpu
);

/**
 * Set the auxiliary fields of the operator.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_lbo_pkpm_nodal_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_pkpm_nodal_auxfields auxin
);
