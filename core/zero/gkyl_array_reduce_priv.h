#pragma once

#include <gkyl_array_reduce.h>

#ifdef GKYL_HAVE_CUDA

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array *inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_min_cu(double *out_d, const struct gkyl_array *inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_min_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_sum_cu(double *out_d, const struct gkyl_array *inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 */
void gkyl_array_reduce_abs_max_cu(double *out_d, const struct gkyl_array *inp);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 */
void gkyl_array_reduce_sq_sum_cu(double *out_d, const struct gkyl_array *inp);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_range_abs_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_range_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 */
void gkyl_array_reduce_weighted_min_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 */
void gkyl_array_reduce_weighted_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 */
void gkyl_array_reduce_weighted_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 */
void gkyl_array_reduce_weighted_abs_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 */
void gkyl_array_reduce_weighted_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 */
void gkyl_array_reduce_weighted_cell_mean_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_weighted_range_min_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_weighted_range_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_weighted_range_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_weighted_range_abs_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_weighted_range_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
);

/**
 * Reduce device data componentwise, with CPU-wrapper validation already performed.
 *
 * @param out_d Device-accessible output of inp->ncomp doubles, not aliasing inputs.
 * @param inp Input device array.
 * @param wgt Weight device array with matching shape.
 * @param range Selection with the input arrays' parent indexing.
 */
void gkyl_array_reduce_weighted_range_cell_mean_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
);

#endif
