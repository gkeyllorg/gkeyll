#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_proj_exp_on_basis gkyl_proj_exp_on_basis;

/**
 * Create new updater to project alpha*exp(beta*f) onto the basis via quadrature.
 *
 * @param basis Basis object (configuration space).
 * @param num_quad Number of quadrature nodes.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_proj_exp_on_basis *gkyl_proj_exp_on_basis_new(
  const struct gkyl_basis *basis, int num_quad, bool use_gpu);

/**
 * Compute alpha*exp(beta*fIn) via quadrature.
 *
 * @param up Projection updater.
 * @param range Config-space range
 * @param fIn Input scalar field.
 * @param fOut Ouput scalar field.
 */
void gkyl_proj_exp_on_basis_advance(const gkyl_proj_exp_on_basis *up,
  const struct gkyl_range *range, double alpha, double beta, const struct gkyl_array *fIn,
  struct gkyl_array *fOut);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_exp_on_basis_release(gkyl_proj_exp_on_basis *up);
