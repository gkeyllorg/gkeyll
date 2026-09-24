#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_gauss_quad_data.h>

// create range to loop over quadrature points.
static inline struct gkyl_range get_qrange(int dim, int num_quad) {
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<dim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

// Sets weights and basis functions at ords. Returns total
// number of quadrature nodes.
static int
init_quad_values(const struct gkyl_basis *basis, int num_quad,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  double ordinates1[num_quad], weights1[num_quad];

  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  }
  else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }

  struct gkyl_range qrange = get_qrange(ndim, num_quad);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  } else {
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);

    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];

    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu)
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  else
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));

  // copy host array to device array
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}
