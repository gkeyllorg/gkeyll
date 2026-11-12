#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_dg_gaussian_filter.h>
#include <gkyl_dg_gaussian_filter_priv.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <assert.h>

// Sets ordinates, weights and basis functions at ords.
// Returns the total number of quadrature nodes
static int
init_quad_values(
  int cdim, const struct gkyl_basis *basis, int num_quad, struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu
)
{
  double ordinates1[num_quad], weights1[num_quad];
  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }

  int qshape[GKYL_MAX_CDIM];
  for (int i = 0; i < cdim; ++i) {
    qshape[i] = num_quad;
  }
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, cdim, qshape);

  int tot_quad = qrange.volume;
  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, cdim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *ordinates = gkyl_array_cu_dev_new(GKYL_DOUBLE, cdim, tot_quad);
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  } else {
    *ordinates = gkyl_array_new(GKYL_DOUBLE, cdim, tot_quad);
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);

    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i = 0; i < cdim; ++i) {
      ord[i] = ordinates1[iter.idx[i] - qrange.lower[i]];
    }

    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i = 0; i < cdim; ++i) {
      wgt[0] *= weights1[iter.idx[i] - qrange.lower[i]];
    }
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu) {
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  } else {
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  }

  for (int n = 0; n < tot_quad; ++n) {
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));
  }

  // copy host array to device array
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

struct gkyl_dg_gaussian_filter *
gkyl_dg_gaussian_filter_inew(const struct gkyl_dg_gaussian_filter_inp *inp)
{
  gkyl_dg_gaussian_filter *up = gkyl_malloc(sizeof(*up));

  up->conf_grid = *inp->conf_grid;
  up->conf_basis = *inp->conf_basis;

  up->cdim = up->conf_basis.ndim;
  up->num_conf_basis = up->conf_basis.num_basis;
  up->use_gpu = inp->use_gpu;
  up->extend_filter = inp->extend_filter;

  // Initialize data needed for configuration-space quadrature on host.
  // We will use host-side initialization to construct the matrix for
  // the modal to nodal and nodal to modal calculations on host.
  int num_quad = up->conf_basis.poly_order + 1;
  up->tot_quad = init_quad_values(
    up->cdim, &up->conf_basis, num_quad, &up->ordinates, &up->weights, &up->basis_at_ords, false
  );

  // Allocate the memory for computing the specific modal to nodal
  // and nodal to modal calculations on host.
  struct gkyl_mat_mm_array_mem *n2m_mem_ho;
  n2m_mem_ho = gkyl_mat_mm_array_mem_new(
    up->num_conf_basis, up->tot_quad, 1.0, 0.0, GKYL_NO_TRANS, GKYL_NO_TRANS, false
  );
  struct gkyl_mat_mm_array_mem *m2n_mem_ho;
  m2n_mem_ho = gkyl_mat_mm_array_mem_new(
    up->tot_quad, up->num_conf_basis, 1.0, 0.0, GKYL_NO_TRANS, GKYL_NO_TRANS, false
  );

  // Compute the matrix A for the conf modal to nodal and nodal to modal memory on host.
  const double *conf_w = (const double *)up->weights->data;
  const double *confb_o = (const double *)up->basis_at_ords->data;
  for (int n = 0; n < up->tot_quad; ++n) {
    for (int k = 0; k < up->num_conf_basis; ++k) {
      gkyl_mat_set(n2m_mem_ho->A, k, n, conf_w[n] * confb_o[k + up->num_conf_basis * n]);
      gkyl_mat_set(m2n_mem_ho->A, n, k, confb_o[k + up->num_conf_basis * n]);
    }
  }
  // Allocate matrix multiplication memory on either host or device (depending on use_gpu)
  // for use in update. Copy results of host-side initializion to updater.
  up->n2m_mem = gkyl_mat_mm_array_mem_new(
    up->num_conf_basis, up->tot_quad, 1.0, 0.0, GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu
  );
  gkyl_mat_copy(up->n2m_mem->A, n2m_mem_ho->A);
  gkyl_mat_mm_array_mem_release(n2m_mem_ho);

  up->m2n_mem = gkyl_mat_mm_array_mem_new(
    up->num_conf_basis, up->tot_quad, 1.0, 0.0, GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu
  );
  gkyl_mat_copy(up->m2n_mem->A, m2n_mem_ho->A);
  gkyl_mat_mm_array_mem_release(m2n_mem_ho);

  // Configuration space field evaluation at quadrature points and filter array.
  if (up->use_gpu) {
    up->conf_arr_quad =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->conf_basis.num_basis, inp->conf_range_ext->volume);
    up->conf_arr_filter =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, up->conf_basis.num_basis, inp->conf_range_ext->volume);
    // Re-initialize data needed for conf-space quadrature on device.
    up->tot_quad = init_quad_values(
      up->cdim, &up->conf_basis, num_quad, &up->ordinates, &up->weights, &up->basis_at_ords,
      up->use_gpu
    );
  } else {
    up->conf_arr_quad =
      gkyl_array_new(GKYL_DOUBLE, up->conf_basis.num_basis, inp->conf_range_ext->volume);
    up->conf_arr_filter =
      gkyl_array_new(GKYL_DOUBLE, up->conf_basis.num_basis, inp->conf_range_ext->volume);
  }

  return up;
}

void
gkyl_dg_gaussian_filter_advance(
  gkyl_dg_gaussian_filter *up, const struct gkyl_range *conf_range, struct gkyl_array *conf_arr
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_dg_gaussian_filter_advance_cu(up, conf_range, conf_arr);
  }
#endif

  int cdim = up->cdim;
  int tot_quad = up->tot_quad;
  int num_conf_basis = up->num_conf_basis;
  gkyl_array_clear(up->conf_arr_quad, 0.0);
  gkyl_array_clear(up->conf_arr_filter, 0.0);

  // Convert from modal basis to nodal quadrature basis.
  gkyl_mat_mm_array(up->m2n_mem, conf_arr, up->conf_arr_quad);

  int cidx_l[GKYL_MAX_CDIM], cidx_c[GKYL_MAX_CDIM], cidx_r[GKYL_MAX_CDIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  double xmu[GKYL_MAX_DIM];
  double xmul[GKYL_MAX_DIM], xmuc[GKYL_MAX_DIM], xmur[GKYL_MAX_DIM];
  double dx[GKYL_MAX_CDIM];
  for (int d = 0; d < cdim; ++d) {
    dx[d] = up->conf_grid.dx[d];
  }
  // Loop over range and compute the Gaussian filter.
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    gkyl_copy_int_arr(cdim, conf_iter.idx, cidx_l);
    gkyl_copy_int_arr(cdim, conf_iter.idx, cidx_c);
    gkyl_copy_int_arr(cdim, conf_iter.idx, cidx_r);
    // Assumes filter is only applied in first dimension for now JJ 09/28/25
    cidx_l[0] = cidx_l[0] - 1;
    cidx_r[0] = cidx_r[0] + 1;

    long loc_l = gkyl_range_idx(conf_range, cidx_l);
    long loc_c = gkyl_range_idx(conf_range, cidx_c);
    long loc_r = gkyl_range_idx(conf_range, cidx_r);

    const double *conf_arr_quad_l = gkyl_array_cfetch(up->conf_arr_quad, loc_l);
    const double *conf_arr_quad_c = gkyl_array_cfetch(up->conf_arr_quad, loc_c);
    const double *conf_arr_quad_r = gkyl_array_cfetch(up->conf_arr_quad, loc_r);
    double *conf_arr_filter_d = gkyl_array_fetch(up->conf_arr_filter, loc_c);

    // Fetch the cell center in the lower, center, and upper cells.
    gkyl_rect_grid_cell_center(&up->conf_grid, cidx_l, xcl);
    gkyl_rect_grid_cell_center(&up->conf_grid, cidx_c, xcc);
    gkyl_rect_grid_cell_center(&up->conf_grid, cidx_r, xcr);

    for (int i = 0; i < tot_quad; ++i) {
      // Fetch the quadrature point location for the quadrature point we are updating.
      comp_to_phys(cdim, gkyl_array_cfetch(up->ordinates, i), dx, xcc, xmu);
      // First construct the filter normalization
      double denominator = filter_norm(
        cdim, tot_quad, up->ordinates, up->weights, dx, xmu, xcl, xmul, xcc, xmuc, xcr, xmur
      );
      for (int j = 0; j < tot_quad; ++j) {
        // Fetch the quadrature point locations for all other quadrature points.
        comp_to_phys(cdim, gkyl_array_cfetch(up->ordinates, j), dx, xcl, xmul);
        comp_to_phys(cdim, gkyl_array_cfetch(up->ordinates, j), dx, xcc, xmuc);
        comp_to_phys(cdim, gkyl_array_cfetch(up->ordinates, j), dx, xcr, xmur);
        const double *wq = gkyl_array_cfetch(up->weights, j);
        // Sum the contributions from each cell in the three-cell filter.
        double expamp_l = -pow((xmu[0] - xmul[0]) / (sqrt(2.0) * dx[0]), 2.0);
        double expamp_c = -pow((xmu[0] - xmuc[0]) / (sqrt(2.0) * dx[0]), 2.0);
        double expamp_r = -pow((xmu[0] - xmur[0]) / (sqrt(2.0) * dx[0]), 2.0);
        conf_arr_filter_d[i] +=
          wq[0] / denominator *
          (conf_arr_quad_l[j] * exp(expamp_l) + conf_arr_quad_c[j] * exp(expamp_c) +
           conf_arr_quad_r[j] * exp(expamp_r));
      }
    }
  }

  // Convert back to modal basis
  gkyl_mat_mm_array(up->n2m_mem, up->conf_arr_filter, conf_arr);
}

void
gkyl_dg_gaussian_filter_release(gkyl_dg_gaussian_filter *up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);

  gkyl_mat_mm_array_mem_release(up->n2m_mem);
  gkyl_mat_mm_array_mem_release(up->m2n_mem);

  gkyl_array_release(up->conf_arr_quad);
  gkyl_array_release(up->conf_arr_filter);

  gkyl_free(up);
}
