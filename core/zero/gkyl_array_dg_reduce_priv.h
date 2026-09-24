#pragma once

#include <gkyl_array_dg_reduce.h>

#ifdef GKYL_HAVE_CUDA

/**
 * Max reduce a gkyl_array evaluating the DG field in each cell at Gauss-Legendre nodes
 * and reducing those values. Allows for reducing a single vector component for arrays
 * with multiple scalar fields within it.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_dg_reducec_max_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis);

/**
 * Min reduce a gkyl_array evaluating the DG field in each cell at Gauss-Legendre nodes
 * and reducing those values. Allows for reducing a single vector component for arrays
 * with multiple scalar fields within it.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_dg_reducec_min_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis);

/**
 * Sum reduce a gkyl_array evaluating the DG field in each cell at Gauss-Legendre nodes
 * and reducing those values. Allows for reducing a single vector component for arrays
 * with multiple scalar fields within it.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_dg_reducec_sum_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis);

/**
 * Max reduce a gkyl_array evaluating the DG field in each cell (within the input range)
 * at Gauss-Legendre nodes and reducing those values. Allows for reducing a single vector
 * component for arrays with multiple scalar fields within it.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_dg_reducec_range_max_cu(double *out_d, const struct gkyl_array* inp,
  int comp, const struct gkyl_basis *basis, const struct gkyl_range *range);

/**
 * Min reduce a gkyl_array evaluating the DG field in each cell (within the input range)
 * at Gauss-Legendre nodes and reducing those values. Allows for reducing a single vector
 * component for arrays with multiple scalar fields within it.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_dg_reducec_range_min_cu(double *out_d, const struct gkyl_array* inp,
  int comp, const struct gkyl_basis *basis, const struct gkyl_range *range);

/**
 * Sum reduce a gkyl_array evaluating the DG field in each cell (within the input range)
 * at Gauss-Legendre nodes and reducing those values. Allows for reducing a single vector
 * component for arrays with multiple scalar fields within it.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_dg_reducec_range_sum_cu(double *out_d, const struct gkyl_array* inp,
  int comp, const struct gkyl_basis *basis, const struct gkyl_range *range);

#endif
