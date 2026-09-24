#include <gkyl_array_dg_reduce.h>
#include <gkyl_array_dg_reduce_priv.h>
#include <assert.h>
#include <float.h>

void 
gkyl_array_dg_reducec(double *out, const struct gkyl_array *arr, int comp,
  enum gkyl_array_op op, const struct gkyl_basis *basis)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_dg_reducec_max_cu(out, arr, comp, basis);
        break;                                            
      case GKYL_MIN:                                      
        gkyl_array_dg_reducec_min_cu(out, arr, comp, basis);
        break;                                            
      case GKYL_SUM:                                      
        gkyl_array_dg_reducec_sum_cu(out, arr, comp, basis);
        break;
    }
    return;
  }
#endif

  double *arr_d = arr->data;

  int num_nodes = basis->num_quad;

  switch (op) {
    case GKYL_MIN:
      out[0] = DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        double arr_nodal[num_nodes];
        for (int k=0; k<num_nodes; ++k) {
          basis->modal_to_quad_nodal(&d[comp*basis->num_basis], arr_nodal, k);
          out[0] = fmin(out[0], arr_nodal[k]);
        }
      }
      break;

    case GKYL_MAX:
      out[0] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        double arr_nodal[num_nodes];
        for (int k=0; k<num_nodes; ++k) {
          basis->modal_to_quad_nodal(&d[comp*basis->num_basis], arr_nodal, k);
          out[0] = fmax(out[0], arr_nodal[k]);
        }
      }
      break;

    case GKYL_SUM:
      out[0] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        double arr_nodal[num_nodes];
        for (int k=0; k<num_nodes; ++k) {
          basis->modal_to_quad_nodal(&d[comp*basis->num_basis], arr_nodal, k);
          out[0] += arr_nodal[k];
        }
      }
      break;
  }
}

void 
gkyl_array_dg_reducec_range(double *out, const struct gkyl_array *arr, int comp,
  enum gkyl_array_op op, const struct gkyl_basis *basis, const struct gkyl_range *range)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_dg_reducec_range_max_cu(out, arr, comp, basis, range);
        break;
      case GKYL_MIN:
        gkyl_array_dg_reducec_range_min_cu(out, arr, comp, basis, range);
        break;
      case GKYL_SUM:
        gkyl_array_dg_reducec_range_sum_cu(out, arr, comp, basis, range);
        break;
    }
    return;
  }
#endif

  double *arr_d = arr->data;

  int num_nodes = basis->num_quad;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  switch (op) {
    case GKYL_MIN:
      out[0] = DBL_MAX;
      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        double arr_nodal[num_nodes];
        for (int k=0; k<num_nodes; ++k) {
          basis->modal_to_quad_nodal(&d[comp*basis->num_basis], arr_nodal, k);
          out[0] = fmin(out[0], arr_nodal[k]);
        }
      }
      break;

    case GKYL_MAX:
      out[0] = -DBL_MAX;
      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        double arr_nodal[num_nodes];
        for (int k=0; k<num_nodes; ++k) {
          basis->modal_to_quad_nodal(&d[comp*basis->num_basis], arr_nodal, k);
          out[0] = fmax(out[0], arr_nodal[k]);
        }
      }
      break;

    case GKYL_SUM:
      out[0] = 0;
      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        double arr_nodal[num_nodes];
        for (int k=0; k<num_nodes; ++k) {
          basis->modal_to_quad_nodal(&d[comp*basis->num_basis], arr_nodal, k);
          out[0] += arr_nodal[k];
        }
      }
      break;
  }
}

