#pragma once

#include <gkyl_array_reduce.h>

#ifdef GKYL_HAVE_CUDA

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_min_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_min_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_sum_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

#endif
