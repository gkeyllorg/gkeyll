#pragma once

#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_fem_poisson gkyl_fem_poisson;

/**
 * Create new updater to solve the Helmholtz problem
 *   - nabla . (epsilon * nabla phi) - kSq * phi = rho
 * using a FEM to ensure phi is continuous. This solver is also
 * used as a Poisson solver by passing a zero kSq. The input is the
 * DG field rho, which is translated to FEM. The output is the
 * DG field phi, after we've translated the FEM solution to DG.
 * Free using gkyl_fem_poisson_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param bcs Boundary conditions.
 * @param bias_plane_list List of points (1D), lines (2D), planes (3D) to bias.
 * @param epsilon_var Permittivity tensor. Defined over the extended range.
 * @param kSq Squared wave number (factor multiplying phi in Helmholtz eq).
 * @param is_epsilon_const =true if permittivity is constant in space.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_poisson* gkyl_fem_poisson_new(
  const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid, const struct gkyl_basis basis,
  struct gkyl_poisson_bc *bcs, struct gkyl_poisson_bias_plane_list* bias_plane_list, struct gkyl_array *epsilon_var,
  struct gkyl_array *kSq, bool is_epsilon_const, bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM poisson updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Spatially varying BC as a DG (volume) field, defined in the whole
                domain but really only applicable to and used in the skin cell.
 */
void gkyl_fem_poisson_set_rhs(gkyl_fem_poisson* up, struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_poisson_solve(gkyl_fem_poisson* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_release(gkyl_fem_poisson *up);
