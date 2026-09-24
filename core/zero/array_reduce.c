#include <gkyl_array_reduce.h>
#include <gkyl_array_reduce_priv.h>
#include <assert.h>
#include <float.h>

void 
gkyl_array_reduce(double *out, const struct gkyl_array *arr, enum gkyl_array_op op)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_max_cu(out, arr);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_min_cu(out, arr);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_sum_cu(out, arr);
        break;
    }
    return;
  }
#endif

  long nc = arr->ncomp;
  double *arr_d = arr->data;

  switch (op) {
    case GKYL_MIN:
      for (long k=0; k<nc; ++k) out[k] = DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmin(out[k], d[k]);
      }
      break;

    case GKYL_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], d[k]);
      }
      break;

    case GKYL_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] += d[k];
      }
      break;
  }
}

void
gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_range_max_cu(res, arr, range);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_range_min_cu(res, arr, range);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_range_sum_cu(res, arr, range);
        break;
    }
    return;
  }
#endif

  long n = arr->ncomp;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  switch (op) {
    case GKYL_MIN:
      for (long i=0; i<n; ++i) res[i] = DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmin(res[i], d[i]);
      }
      break;
    case GKYL_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], d[i]);
      }
      break;
    case GKYL_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] += d[i];
      }
      break;
  }
}

