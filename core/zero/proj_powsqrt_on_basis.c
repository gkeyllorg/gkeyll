#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_powsqrt_on_basis.h>
#include <gkyl_proj_powsqrt_on_basis_priv.h>
#include <gkyl_range.h>

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

gkyl_proj_powsqrt_on_basis*
gkyl_proj_powsqrt_on_basis_new(const struct gkyl_basis *basis, int num_quad, bool use_gpu)
{
  gkyl_proj_powsqrt_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_powsqrt_on_basis));

  up->ndim = basis->ndim;
  up->num_quad = num_quad;
  up->num_basis = basis->num_basis;
  up->use_gpu = use_gpu;

  // initialize data needed for quadrature
  up->tot_quad = init_quad_values(basis, num_quad, &up->weights,
    &up->basis_at_ords, use_gpu);

  if (up->use_gpu)
    up->fun_at_ords = NULL;
  else
    up->fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad); // Only used in CPU implementation.

  return up;
}

static void
proj_on_basis(const gkyl_proj_powsqrt_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;

  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

void
gkyl_proj_powsqrt_on_basis_advance(const gkyl_proj_powsqrt_on_basis *up,
  const struct gkyl_range *range, double expIn, const struct gkyl_array *fIn,
  struct gkyl_array *fOut)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_proj_powsqrt_on_basis_advance_cu(up, range, expIn, fIn, fOut);
#endif

  // Create range to loop over quadrature points.
  struct gkyl_range qrange = get_qrange(up->ndim, up->num_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(range, iter.idx);

    const double *fIn_d = gkyl_array_cfetch(fIn, linidx);

    struct gkyl_range_iter qiter;
    gkyl_range_iter_init(&qiter, &qrange);
    while (gkyl_range_iter_next(&qiter)) {
      
      int qidx = gkyl_range_idx(&qrange, qiter.idx);

      // Evaluate input function f at quad point.
      const double *b_ord = gkyl_array_cfetch(up->basis_at_ords, qidx);
      double fIn_q=0.;
      for (int k=0; k<up->num_basis; ++k)
        fIn_q += fIn_d[k]*b_ord[k];

      double *fOut_q = gkyl_array_fetch(up->fun_at_ords, qidx);

      // Evaluate pow(sqrt()) at quad point.
      fOut_q[0] = (fIn_q < 0.) ? 1.e-40 : pow(sqrt(fIn_q), expIn);
    }

    // Compute expansion coefficients.
    proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(fOut, linidx));
  }

}

void
gkyl_proj_powsqrt_on_basis_release(gkyl_proj_powsqrt_on_basis* up)
{
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  if (!up->use_gpu)
    gkyl_array_release(up->fun_at_ords);
  gkyl_free(up);
}
