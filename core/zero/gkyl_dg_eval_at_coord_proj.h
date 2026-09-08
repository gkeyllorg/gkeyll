#pragma once

#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type.
typedef struct gkyl_dg_eval_at_coord_proj gkyl_dg_eval_at_coord_proj;

/**
 * Create an updater to down-project a higher-dimensional donor DG field
 * onto lower dimensions by evaluating the polynomial at specified
 * reference coordinates and projecting onto the target basis.
 *
 * Supported cases (serendipity basis):
 *   1x-3x : p1 and p2; any subset of directions
 *   4x    : p1 and p2; any subset of directions
 *   5x-6x : p1 only;   any subset of directions
 *
 * @param cdim_do       Number of configuration-space dimensions in the donor field.
 * @param basis_do      Basis of the donor field (ndim may exceed cdim_do for phase-space fields).
 * @param num_eval_dirs Number of directions to evaluate / eliminate.
 * @param eval_dirs     Sorted list of direction indices to evaluate (0-based).
 * @param use_gpu       Whether to run on GPU.
 * @return New updater object.
 */
struct gkyl_dg_eval_at_coord_proj *gkyl_dg_eval_at_coord_proj_new(int cdim_do,
  const struct gkyl_basis *basis_do, int num_eval_dirs, const int *eval_dirs, bool use_gpu);

/**
 * Apply the dg_eval_at_coord operator.
 *
 * @param up          Updater object.
 * @param eval_coords Computational coords for each evaluated direction.
 * @param grid        Computational grid object.
 * @param pick_lower  When coordinate is at cell boundary, this instructs to pick the lower cell.
 * @param known_index Cell index owning eval_coords, if known. Entry <0 means not known. May be NULL.
 * @param rng_do      Range of the donor field.
 * @param rng_tar     Range of the target field.
 * @param fdo         Donor DG field.
 * @param ncomp       Number of scalar DG fields packed in fdo/ftar.
 * @param ftar        Target DG field (written by this call).
 */
void gkyl_dg_eval_at_coord_proj_advance(struct gkyl_dg_eval_at_coord_proj *up,
  const double *eval_coords, const struct gkyl_rect_grid *grid, const bool *pick_lower,
  const int *known_index, const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Return the properties of the target basis.
 *
 * @param up         Updater object.
 * @param cdim       Conf-space dimensionality.
 * @param ndim       Dimensionality.
 * @param btype      Basis type.
 * @param poly_order Polynomial order.
 * @param num_basis  Number of basis functions.
 */
void gkyl_dg_eval_at_coord_proj_target_basis(struct gkyl_dg_eval_at_coord_proj *up, int *cdim,
  int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis);

/**
 * Release the memory associated with this updater.
 *
 * @param up Updater to release.
 */
void gkyl_dg_eval_at_coord_proj_release(struct gkyl_dg_eval_at_coord_proj *up);
