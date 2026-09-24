#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_array_integrate gkyl_array_integrate;

enum gkyl_array_integrate_op {
  GKYL_ARRAY_INTEGRATE_OP_NONE = 0,  // int dx f
  GKYL_ARRAY_INTEGRATE_OP_ABS,  // int dx |f|
  GKYL_ARRAY_INTEGRATE_OP_SQ,  // int dx f^2
  GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED,  // int dx w * f^2
  GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ,  // int dx |nabla f|^2
  GKYL_ARRAY_INTEGRATE_OP_GRADPERP_SQ,  // int dx |nabla_perp f|^2
  GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ,  // int dx epsilon*|nabla_perp f|^2
};

/**
 * Create a new updater that integrates a gkyl_array, with an option to perform
 * additional operations during the integration (e.g. square, abs, etc). *
 *
 * @param grid Grid array is defined on.
 * @param basis Basis array is defined on.
 * @param num_comp Number of (vector) components in the array.
 * @param op Additional operator to apply in very cell.
 * @param use_gpu Indicate whether to perform integral on the device.
 */
struct gkyl_array_integrate*
gkyl_array_integrate_new(const struct gkyl_rect_grid* grid, const struct gkyl_basis* basis,
  int num_comp, enum gkyl_array_integrate_op op, bool use_gpu);

/**
 * Compute the array integral.
 *
 * @param up array_integrate updater.
 * @param fin Input gkyl_array.
 * @param factor Factor to multiply by.
 * @param weight Weight field we multiply by inside the integral.
 * @param range Range we'll integrate over.
 * @param weight_range Range of the weight.
 * @return out Output integral result(s). On device memory if use_gpu=true.
 */
void gkyl_array_integrate_advance(gkyl_array_integrate *up, const struct gkyl_array *fin,
  double factor, const struct gkyl_array *weight, const struct gkyl_range *range,
  const struct gkyl_range *weight_range, double *out);

/**
 * Release memory associated with this updater.
 *
 * @param up array_integrate updater.
 */
void gkyl_array_integrate_release(gkyl_array_integrate *up);
