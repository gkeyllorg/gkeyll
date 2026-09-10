#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_on_basis_priv.h>
#include <gkyl_range.h>

// Identity comp to phys coord mapping, for when user doesn't provide a map.
static inline void
c2p_identity(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_rect_grid *grid = ctx;
  int ndim = grid->ndim;
  for (int d=0; d<ndim; d++) xphys[d] = xcomp[d];
}

struct gkyl_proj_on_basis*
gkyl_proj_on_basis_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  return gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = grid,
      .basis = basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = num_quad,
      .num_ret_vals = num_ret_vals,
      .eval = eval,
      .ctx = ctx,
      .c2p_func = 0,
      .c2p_func_ctx = NULL,
      .use_gpu = false,
    }
  );
}

struct gkyl_proj_on_basis*
gkyl_proj_on_basis_inew(const struct gkyl_proj_on_basis_inp *inp)
{
  struct gkyl_proj_on_basis *up = gkyl_malloc(sizeof(struct gkyl_proj_on_basis));

  up->grid = *inp->grid;
  int num_quad = up->num_quad = inp->num_quad == 0 ? inp->basis->poly_order+1 : inp->num_quad;
  int num_ret_vals = up->num_ret_vals = inp->num_ret_vals;
  up->eval = inp->eval;
  up->ctx = inp->ctx;
  up->num_basis = inp->basis->num_basis;
  up->use_gpu = inp->use_gpu;

  if (up->use_gpu) {
    // On the GPU eval and c2p_func must be device function pointers; a NULL
    // c2p_func selects an identity mapping assigned on the device.
    up->c2p = inp->c2p_func;
    up->c2p_ctx = inp->c2p_func_ctx;
  }
  else if (inp->c2p_func == 0) {
    up->c2p = c2p_identity;
    up->c2p_ctx = &up->grid; // Use grid as the context since all we need is ndim.
  }
  else {
    up->c2p = inp->c2p_func;
    up->c2p_ctx = inp->c2p_func_ctx;
  }

  double ordinates1[num_quad], weights1[num_quad];

  if (inp->qtype == GKYL_GAUSS_QUAD) {
    if (num_quad <= gkyl_gauss_max) {
      // use pre-computed values if possible (these are more accurate
      // than computing them on the fly)
      memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
      memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
    }
    else {
      gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
    }
  }
  else if (inp->qtype == GKYL_GAUSS_LOBATTO_QUAD) {
    assert( (num_quad > 1) && (num_quad <= gkyl_gauss_max) );

    // Gauss-Lobatto quadrature
    memcpy(ordinates1, gkyl_gauss_lobatto_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_lobatto_weights[num_quad], sizeof(double[num_quad]));
  }
  else {
    fprintf(stderr, "Quadrature type not available. Exiting... \n");
    assert(false);
  }

  // create range to loop over quadrature points
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<inp->grid->ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, inp->grid->ndim, qshape);

  int tot_quad = up->tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  up->ordinates = gkyl_array_new(GKYL_DOUBLE, inp->grid->ndim, tot_quad);
  up->weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    long node = gkyl_range_idx(&qrange, iter.idx);

    // set ordinates
    double *ord = gkyl_array_fetch(up->ordinates, node);
    for (int i=0; i<inp->grid->ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];

    // set weights
    double *wgt = gkyl_array_fetch(up->weights, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  up->basis_at_ords = gkyl_array_new(GKYL_DOUBLE, inp->basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    inp->basis->eval(gkyl_array_fetch(up->ordinates, n),
      gkyl_array_fetch(up->basis_at_ords, n));

  up->ordinates_cu = 0;
  up->weights_cu = 0;
  up->basis_at_ords_cu = 0;
  up->on_dev = up; // On the CPU the updater points to itself.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Mirror the quadrature data on the device and create the device clone.
    up->ordinates_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, inp->grid->ndim, tot_quad);
    up->weights_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
    up->basis_at_ords_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, inp->basis->num_basis, tot_quad);
    gkyl_array_copy(up->ordinates_cu, up->ordinates);
    gkyl_array_copy(up->weights_cu, up->weights);
    gkyl_array_copy(up->basis_at_ords_cu, up->basis_at_ords);
    up->on_dev = gkyl_proj_on_basis_cu_dev_new(up);
  }
#else
  assert(!up->use_gpu);
#endif

  return up;
}

int gkyl_proj_on_basis_get_tot_quad(const struct gkyl_proj_on_basis *up)
{
  return up->tot_quad;
}

double* gkyl_proj_on_basis_fetch_ordinate(const struct gkyl_proj_on_basis *up, long node)
{
  return gkyl_array_fetch(up->ordinates, node);
}

void
gkyl_proj_on_basis_quad(const struct gkyl_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis;
  int tot_quad = up->tot_quad;
  int num_ret_vals = up->num_ret_vals;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  // arrangement of f is as:
  // c0[0], c0[1], ... c1[0], c1[1], ....
  // where c0, c1, ... are components of f (num_ret_vals)
  int offset = 0;
  for (int n=0; n<num_ret_vals; ++n) {
    for (int k=0; k<num_basis; ++k) f[offset+k] = 0.0;

    for (int imu=0; imu<tot_quad; ++imu) {
      double tmp = weights[imu]*func_at_ords[n+num_ret_vals*imu];
      for (int k=0; k<num_basis; ++k)
        f[offset+k] += tmp*basis_at_ords[k+num_basis*imu];
    }
    offset += num_basis;
  }
}

void
gkyl_proj_on_basis_advance(const struct gkyl_proj_on_basis *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_proj_on_basis_advance_cu(up, tm, update_range, arr);
    return;
  }
#endif

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  int tot_quad = up->tot_quad;
  struct gkyl_array *fun_at_ords = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i=0; i<tot_quad; ++i) {
      proj_on_basis_log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->ordinates, i),
        up->grid.dx, xc, xmu);
      up->c2p(xmu, xmu, up->c2p_ctx);
      up->eval(tm, xmu, gkyl_array_fetch(fun_at_ords, i), up->ctx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    gkyl_proj_on_basis_quad(up, fun_at_ords, gkyl_array_fetch(arr, lidx));
  }

  gkyl_array_release(fun_at_ords);
}

void
gkyl_proj_on_basis_release(struct gkyl_proj_on_basis* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_release(up->ordinates_cu);
    gkyl_array_release(up->weights_cu);
    gkyl_array_release(up->basis_at_ords_cu);
    gkyl_cu_free(up->on_dev);
  }
#endif
  gkyl_free(up);
}
