#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_fem_poisson_bctype.h>

// Object type
typedef struct gkyl_fem_parproj gkyl_fem_parproj;

// Boundary condition types.
enum gkyl_fem_parproj_bc_type {
  GKYL_FEM_PARPROJ_PERIODIC = 0,
  GKYL_FEM_PARPROJ_DIRICHLET_GHOST, // Solution = ghost evaluated at the boundary.
  GKYL_FEM_PARPROJ_DIRICHLET_SKIN, // Solution = skin evaluated at the boundary.
  GKYL_FEM_PARPROJ_NONE,      // does not enforce a BC.
};

/**
 * Create new updater to project a DG field onto the FEM (nodal) basis
 * in order to make the field continuous or, thanks to the option to pass
 * a multiplicative weight, solve 1D algebraic equations in which the output
 * field is continuous (but the input may not be). That is, we solve
 *    wgtL*phi_{fem} \doteq wgtR*rho_{dg}
 * where wgt is the weight field, phi_{fem} is the (continuous field)
 * we wish to compute, rho_{dg} is the (discontinuous) input source field,
 * and \doteq implies weak equality with respect to the FEM basis.
 * Free using gkyl_fem_parproj_release method.
 *
 * @param solve_range Range in which to perform the projection operation.
 * @param grid Grid object.
 * @param basis Basis functions of the DG field.
 * @param bctype Type of boundary condition (see gkyl_fem_parproj_bc_type).
 * @param bias_line_list List of points (in 2x solves) or lines (3x solves) to bias.
 * @param weight_left Weight on left-side of the operator (time-independent).
 * @param weight_right Weight on right-side of the operator (time-independent).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_parproj* gkyl_fem_parproj_new(const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype, struct gkyl_poisson_bias_line_list* bias_line_list,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right,
  bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs.
 */
void gkyl_fem_parproj_set_rhs(struct gkyl_fem_parproj* up,
  const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve(struct gkyl_fem_parproj* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_release(struct gkyl_fem_parproj *up);
