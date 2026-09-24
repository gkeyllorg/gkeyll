#pragma once

#include <float.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

/**
 * Perform a "reduce" operation of data in the array accounting for the DG
 * representation within the array. It evaluates the DG field at Gauss-Legendre
 * nodes in each cell, and reduces them before reducing over cells. It allows
 * specification of the vector component for arrays that contain multiple
 * scalar DG fields.
 *
 * @param out Reduced output value (size 1).
 * @param arr Array to perform reduction on.
 * @param comp Vector component in arr to reduce.
 * @param op Reduction operators.
 * @param basis Basis DG coefficients are expanded in (device pointer if use_gpu=true).
 */
void gkyl_array_dg_reducec(double *out, const struct gkyl_array *arr, int comp,
  enum gkyl_array_op op, const struct gkyl_basis *basis);

/**
 * Perform a "reduce" operation of data in the array accounting for the DG
 * representation within the array. It evaluates the DG field at Gauss-Legendre
 * nodes in each cell, and reduces them before reducing over cells. It allows
 * specification of the vector component for arrays that contain multiple
 * scalar DG fields.
 *
 * @param out Reduced output value (size 1).
 * @param arr Array to perform reduction on.
 * @param comp Vector component in arr to reduce.
 * @param op Reduction operators.
 * @param basis Basis DG coefficients are expanded in (device pointer if use_gpu=true).
 * @param range Range specifying region.
 */
void gkyl_array_dg_reducec_range(double *out, const struct gkyl_array *arr, int comp,
  enum gkyl_array_op op, const struct gkyl_basis *basis, const struct gkyl_range *range);


