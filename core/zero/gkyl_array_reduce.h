#pragma once

#include <float.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

/**
 * Reduce a double array over all allocated cells, separately per component.
 * Supported operators are MIN, MAX, SUM, ABS_MAX, and SQ_SUM. The latter two
 * calculate max(abs(x)) and sum(x*x). No MPI reduction is performed.
 * Empty reductions return DBL_MAX for MIN, -DBL_MAX for MAX, zero otherwise.
 *
 * @param res Output ncomp doubles, not overlapping arr; host storage for CPU,
 * device-accessible storage (device or mapped pinned memory) for CUDA.
 * @param arr Input double array.
 * @param op Supported reduction operator, prefixed by GKYL_.
 */
void gkyl_array_reduce(double *res, const struct gkyl_array *arr, enum gkyl_array_op op);

/**
 * Reduce selected cells componentwise with gkyl_array_reduce's operators and
 * empty-range identities. This is a local reduction, excluding cells outside
 * the supplied range.
 *
 * @param res Output ncomp doubles, not overlapping arr and accessible to its backend.
 * @param arr Input double array.
 * @param op Supported reduction operator.
 * @param range Range with indexing inherited from arr's parent range.
 */
void gkyl_array_reduce_range(
  double *res, const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range
);

/**
 * Reduce a weighted double array over all allocated cells, separately per
 * component. arr and wgt must have matching sizes, positive component counts,
 * and backends. For component k, supported formulas over cells i are:
 * MIN/MAX/SUM: min/max/sum(w[i,k]*x[i,k]); ABS_MAX: max(abs(w[i,k]*x[i,k]));
 * SQ_SUM: sum(w[i,k]*x[i,k]^2); CELL_MEAN_SQ_SUM:
 * sum(w[i,k]*(sum_j x[i,j]^2 / ncomp)). The weight appears once in all formulas.
 * CELL_MEAN_SQ_SUM is supported only by the weighted API. This routine does not
 * take a final square root, normalize over cells, or perform an MPI reduction.
 * Empty reductions return DBL_MAX for MIN, -DBL_MAX for MAX, zero otherwise.
 *
 * @param res Output ncomp doubles, not overlapping either input; host storage for
 * CPU, device-accessible storage (device or mapped pinned memory) for CUDA.
 * @param arr Input double array x.
 * @param wgt Input double array w, with one weight per cell/component.
 * @param op Supported reduction operator, prefixed by GKYL_.
 */
void gkyl_array_reduce_weighted(
  double *res, const struct gkyl_array *arr, const struct gkyl_array *wgt, enum gkyl_array_op op
);

/**
 * Reduce selected cells using gkyl_array_reduce_weighted's formulas, matching
 * array requirements, and empty-range identities. No MPI reduction is performed.
 *
 * @param res Output ncomp doubles, not overlapping inputs and accessible to their backend.
 * @param arr Input double array.
 * @param wgt Weight array, with one weight per cell/component.
 * @param op Supported weighted reduction operator.
 * @param range Range with indexing inherited from the arrays' parent range.
 */
void gkyl_array_reduce_weighted_range(
  double *res, const struct gkyl_array *arr, const struct gkyl_array *wgt, enum gkyl_array_op op,
  const struct gkyl_range *range
);
