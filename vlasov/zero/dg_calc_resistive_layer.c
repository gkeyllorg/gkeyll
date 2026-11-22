#include <assert.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_calc_resistive_layer.h>
#include <gkyl_dg_calc_resistive_layer_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>

struct gkyl_dg_calc_resistive_layer *
gkyl_dg_calc_resistive_layer_inew(const struct gkyl_dg_calc_resistive_layer_inp *inp)
{
  struct gkyl_dg_calc_resistive_layer *up = gkyl_malloc(sizeof(*up));
  const struct gkyl_basis *basis = inp->basis;
  int cdim = basis->ndim, num_quad = basis->poly_order + 1;
  assert(num_quad <= gkyl_gauss_max);
  up->num_basis = basis->num_basis;
  up->use_gpu = inp->use_gpu;

  // Tensor Gauss-Legendre points, row-major (the last direction fastest).
  int qshape[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    qshape[d] = num_quad;
  }
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, cdim, qshape);
  up->tot_quad = qrange.volume;
  struct gkyl_array *ords_ho = gkyl_array_new(GKYL_DOUBLE, cdim, up->tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  struct gkyl_array *basis_ho = gkyl_array_new(GKYL_DOUBLE, up->num_basis, up->tot_quad);
  struct gkyl_range_iter qiter;
  gkyl_range_iter_init(&qiter, &qrange);
  while (gkyl_range_iter_next(&qiter)) {
    long q = gkyl_range_idx(&qrange, qiter.idx);
    double *ord = gkyl_array_fetch(ords_ho, q), *w = gkyl_array_fetch(weights_ho, q);
    w[0] = 1.0;
    for (int d = 0; d < cdim; ++d) {
      ord[d] = gkyl_gauss_ordinates[num_quad][qiter.idx[d]];
      w[0] *= gkyl_gauss_weights[num_quad][qiter.idx[d]];
    }
    basis->eval(ord, gkyl_array_fetch(basis_ho, q));
  }

  // sigma at the quadrature points of each cell, in physical coordinates.
  assert(inp->sigma->ncomp == up->tot_quad);
  struct gkyl_array *sigma_ho = gkyl_array_new(GKYL_DOUBLE, up->tot_quad, inp->sigma->size);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, inp->range);
  while (gkyl_range_iter_next(&iter)) {
    double xc[GKYL_MAX_DIM], xcomp[GKYL_MAX_DIM], xphys[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(inp->grid, iter.idx, xc);
    double *sig = gkyl_array_fetch(sigma_ho, gkyl_range_idx(inp->range, iter.idx));
    for (int q = 0; q < up->tot_quad; ++q) {
      const double *ord = gkyl_array_cfetch(ords_ho, q);
      for (int d = 0; d < cdim; ++d) {
        xcomp[d] = xc[d] + 0.5 * inp->grid->dx[d] * ord[d];
      }
      if (inp->c2p_func) {
        inp->c2p_func(xcomp, xphys, inp->c2p_func_ctx);
      } else {
        memcpy(xphys, xcomp, sizeof(double[cdim]));
      }
      inp->eval(0.0, xphys, &sig[q], inp->ctx);
    }
  }

  if (up->use_gpu) {
    up->weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, up->tot_quad);
    up->basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_basis, up->tot_quad);
  } else {
    up->weights = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
    up->basis_at_ords = gkyl_array_new(GKYL_DOUBLE, up->num_basis, up->tot_quad);
  }
  gkyl_array_copy(up->weights, weights_ho);
  gkyl_array_copy(up->basis_at_ords, basis_ho);
  gkyl_array_copy(inp->sigma, sigma_ho);
  up->sigma = gkyl_array_acquire(inp->sigma);

  gkyl_array_release(ords_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_ho);
  gkyl_array_release(sigma_ho);
  return up;
}

void
gkyl_dg_calc_resistive_layer_advance(
  const struct gkyl_dg_calc_resistive_layer *up, const struct gkyl_range *range,
  const struct gkyl_array *em, struct gkyl_array *rhs
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_dg_calc_resistive_layer_advance_cu(up, range, em, rhs);
  }
#endif
  const double *weights = up->weights->data, *basis_at_ords = up->basis_at_ords->data;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    dg_calc_resistive_layer_cell(
      up->num_basis, up->tot_quad, weights, basis_at_ords, gkyl_array_cfetch(up->sigma, loc),
      gkyl_array_cfetch(em, loc), gkyl_array_fetch(rhs, loc)
    );
  }
}

void
gkyl_dg_calc_resistive_layer_release(struct gkyl_dg_calc_resistive_layer *up)
{
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->sigma);
  gkyl_free(up);
}
