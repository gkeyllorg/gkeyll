// Tests of the resistive layer updater (gkyl_dg_calc_resistive_layer):
//   * sigma is sampled at the (p+1)^cdim Gauss-Legendre points of each cell, row-major (the last
//     direction fastest), at the physical coordinates given by the position map;
//   * a uniform sigma gives -sigma E exactly (p = 1, 2, 3; 1x, 2x, 3x) and leaves B alone;
//   * at p = 1 the increment is the exact weak product of E and the p=1 interpolant of sigma;
//   * at p = 2 and 3 the operator E -> sigma E of each cell is the matrix
//     sum_q w_q sigma_q b_j(x_q) b_k(x_q), positive semidefinite for a sigma that steps inside a
//     cell, while the weak product with the L2 projection of the step has a negative eigenvalue
//     (it adds energy for some E).
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_resistive_layer.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

struct setup {
  int cdim;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
};

static void
setup_init(struct setup *s, int cdim, int poly_order, bool tensor)
{
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {4, 3, 2}, ghost[3] = {1, 1, 1};
  s->cdim = cdim;
  gkyl_rect_grid_init(&s->grid, cdim, lower, upper, cells);
  if (tensor) {
    gkyl_cart_modal_tensor(&s->basis, cdim, poly_order);
  } else {
    gkyl_cart_modal_serendip(&s->basis, cdim, poly_order);
  }
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
}

static int
num_quad_tot(const struct setup *s)
{
  int n = 1;
  for (int d = 0; d < s->cdim; ++d) {
    n *= s->basis.poly_order + 1;
  }
  return n;
}

static void
eval_sigma_smooth(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 2.0 + sin(3.0 * xn[0]) + 0.5 * cos(2.0 * y) * cos(z);
}

static void
eval_sigma_uniform(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = *(double *)ctx;
}

// Zero for x < -0.3 and 5 after it: the step is inside a cell of the 4-cell grid.
static void
eval_sigma_step(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = xn[0] < -0.3 ? 0.0 : 5.0;
}

// Smallest eigenvalue of a symmetric n x n matrix (power iteration on shift I - A).
static double
min_eig(int n, const double *A)
{
  double shift = 0.0, v[64], w[64];
  for (int i = 0; i < n; ++i) {
    double r = 0.0;
    for (int j = 0; j < n; ++j) {
      r += fabs(A[i * n + j]);
    }
    shift = fmax(shift, r);
    v[i] = 1.0 + 0.1 * i;
  }
  if (shift == 0.0) {
    return 0.0;
  }
  for (int it = 0; it < 20000; ++it) {
    double norm = 0.0;
    for (int i = 0; i < n; ++i) {
      w[i] = shift * v[i];
      for (int j = 0; j < n; ++j) {
        w[i] -= A[i * n + j] * v[j];
      }
      norm += w[i] * w[i];
    }
    for (int i = 0; i < n; ++i) {
      v[i] = w[i] / sqrt(norm);
    }
  }
  double ray = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ray += v[i] * A[i * n + j] * v[j];
    }
  }
  return ray;
}

static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = sin(2.0 * x) * cos(y) + 0.3 * z;
  fout[1] = cos(3.0 * x + y) - 0.2;
  fout[2] = 0.5 + x * y * z + exp(-x * x);
  fout[3] = 0.7;
  fout[4] = -0.4 * x;
  fout[5] = cos(y);
  fout[6] = 0.1;
  fout[7] = -0.2;
}

static void
c2p_stretch(const double *xcomp, double *xphys, void *ctx)
{
  int cdim = *(int *)ctx;
  for (int d = 0; d < cdim; ++d) {
    xphys[d] = xcomp[d] + 0.1 * xcomp[d] * xcomp[d];
  }
}

static struct gkyl_array *
project(const struct setup *s, int num_quad, int ncomp, evalf_t f, void *ctx)
{
  struct gkyl_array *arr =
    gkyl_array_new(GKYL_DOUBLE, ncomp * s->basis.num_basis, s->local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, num_quad, ncomp, f, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, arr);
  gkyl_proj_on_basis_release(proj);
  return arr;
}

static struct gkyl_dg_calc_resistive_layer *
make_layer(const struct setup *s, evalf_t f, void *ctx, void *c2p, struct gkyl_array *sigma)
{
  return gkyl_dg_calc_resistive_layer_inew(&(struct gkyl_dg_calc_resistive_layer_inp){
    .grid = &s->grid,
    .basis = &s->basis,
    .range = &s->local_ext,
    .eval = f,
    .ctx = ctx,
    .c2p_func = c2p ? c2p_stretch : 0,
    .c2p_func_ctx = c2p,
    .sigma = sigma,
    .use_gpu = false,
  });
}

// Sampling: values at the tensor Gauss-Legendre points, row-major, through the position map.
static void
test_sample(int cdim, int poly_order)
{
  struct setup s;
  setup_init(&s, cdim, poly_order, false);
  int nq = poly_order + 1, ntot = num_quad_tot(&s);
  struct gkyl_array *sigma = gkyl_array_new(GKYL_DOUBLE, ntot, s.local_ext.volume);
  struct gkyl_dg_calc_resistive_layer *up = make_layer(&s, eval_sigma_smooth, &cdim, &cdim, sigma);

  double max_err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    double xc[3];
    gkyl_rect_grid_cell_center(&s.grid, iter.idx, xc);
    const double *sg = gkyl_array_cfetch(sigma, gkyl_range_idx(&s.local_ext, iter.idx));
    for (int q = 0; q < ntot; ++q) {
      double xcomp[3], xphys[3], ref;
      for (int d = 0, rem = q; d < cdim; ++d) {
        int stride = (int)pow(nq, cdim - 1 - d), i = rem / stride;
        rem -= i * stride;
        xcomp[d] = xc[d] + 0.5 * s.grid.dx[d] * gkyl_gauss_ordinates[nq][i];
      }
      c2p_stretch(xcomp, xphys, &cdim);
      eval_sigma_smooth(0.0, xphys, &ref, &cdim);
      max_err = fmax(max_err, fabs(sg[q] - ref));
    }
  }
  TEST_CHECK(max_err < 1e-15);
  TEST_MSG("cdim=%d p=%d: sampling error %g", cdim, poly_order, max_err);
  gkyl_dg_calc_resistive_layer_release(up);
  gkyl_array_release(sigma);
}

// Uniform sigma: rhs gains -sigma E, and B and the other components are untouched.
static void
test_uniform(int cdim, int poly_order, bool tensor)
{
  struct setup s;
  setup_init(&s, cdim, poly_order, tensor);
  int nb = s.basis.num_basis;
  double sig0 = 2.7;
  struct gkyl_array *sigma = gkyl_array_new(GKYL_DOUBLE, num_quad_tot(&s), s.local_ext.volume);
  struct gkyl_dg_calc_resistive_layer *up = make_layer(&s, eval_sigma_uniform, &sig0, 0, sigma);
  struct gkyl_array *em = project(&s, poly_order + 1, 8, eval_em, &cdim);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, 8 * nb, s.local_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_calc_resistive_layer_advance(up, &s.local, em, rhs);

  double max_err = 0.0, max_other = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *e = gkyl_array_cfetch(em, loc), *r = gkyl_array_cfetch(rhs, loc);
    for (int k = 0; k < 3 * nb; ++k) {
      max_err = fmax(max_err, fabs(r[k] + sig0 * e[k]));
    }
    for (int k = 3 * nb; k < 8 * nb; ++k) {
      max_other = fmax(max_other, fabs(r[k]));
    }
  }
  TEST_CHECK(max_err < 1e-13);
  TEST_CHECK(max_other == 0.0);
  TEST_MSG("cdim=%d p=%d tensor=%d: error %g", cdim, poly_order, tensor, max_err);
  gkyl_dg_calc_resistive_layer_release(up);
  gkyl_array_release(sigma);
  gkyl_array_release(em);
  gkyl_array_release(rhs);
}

// p = 1: the increment is minus the exact weak product of E and the interpolant of sigma at the
// Gauss-Legendre points (the 2-point quadrature projection).
static void
test_p1_weak_product(int cdim, bool tensor)
{
  struct setup s;
  setup_init(&s, cdim, 1, tensor);
  int nb = s.basis.num_basis;
  struct gkyl_array *sigma = gkyl_array_new(GKYL_DOUBLE, num_quad_tot(&s), s.local_ext.volume);
  struct gkyl_dg_calc_resistive_layer *up = make_layer(&s, eval_sigma_smooth, &cdim, 0, sigma);
  struct gkyl_array *em = project(&s, 2, 8, eval_em, &cdim);
  struct gkyl_array *sig_modal = project(&s, 2, 1, eval_sigma_smooth, &cdim);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, 8 * nb, s.local_ext.volume);
  struct gkyl_array *prod = gkyl_array_new(GKYL_DOUBLE, 3 * nb, s.local_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_calc_resistive_layer_advance(up, &s.local, em, rhs);
  for (int c = 0; c < 3; ++c) {
    gkyl_dg_mul_op_range(&s.basis, c, prod, 0, sig_modal, c, em, &s.local);
  }

  double max_err = 0.0, max_amp = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *r = gkyl_array_cfetch(rhs, loc), *p = gkyl_array_cfetch(prod, loc);
    for (int k = 0; k < 3 * nb; ++k) {
      max_err = fmax(max_err, fabs(r[k] + p[k]));
      max_amp = fmax(max_amp, fabs(p[k]));
    }
  }
  TEST_CHECK(max_amp > 0.1);
  TEST_CHECK(max_err < 1e-13 * max_amp);
  TEST_MSG("cdim=%d tensor=%d: difference to the weak product %g", cdim, tensor, max_err);
  gkyl_dg_calc_resistive_layer_release(up);
  gkyl_array_release(sigma);
  gkyl_array_release(sig_modal);
  gkyl_array_release(em);
  gkyl_array_release(rhs);
  gkyl_array_release(prod);
}

// p >= 2, sigma stepping inside a cell: per cell, the matrix A_jk of E -> sigma E (column k from
// E = b_k) equals sum_q w_q sigma_q b_j(x_q) b_k(x_q) and is positive semidefinite. The weak product
// with the L2 projection of the step has a negative eigenvalue.
static void
test_dissipation(int poly_order, bool tensor)
{
  struct setup st;
  setup_init(&st, 2, poly_order, tensor);
  const struct setup *s = &st;
  evalf_t f = eval_sigma_step;
  void *ctx = 0;
  int cdim = s->cdim;
  int nb = s->basis.num_basis, nq = poly_order + 1, ntot = num_quad_tot(s);
  struct gkyl_array *sigma = gkyl_array_new(GKYL_DOUBLE, ntot, s->local_ext.volume);
  struct gkyl_dg_calc_resistive_layer *up = make_layer(s, f, ctx, 0, sigma);
  struct gkyl_array *sig_modal = project(s, 32, 1, f, ctx);
  struct gkyl_array *em = gkyl_array_new(GKYL_DOUBLE, 8 * nb, s->local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, 8 * nb, s->local_ext.volume);
  struct gkyl_array *prod = gkyl_array_new(GKYL_DOUBLE, nb, s->local_ext.volume);
  long ncell = s->local_ext.volume;
  double *A = gkyl_malloc(sizeof(double[ncell * nb * nb]));
  double *Aw = gkyl_malloc(sizeof(double[ncell * nb * nb]));
  for (int k = 0; k < nb; ++k) {
    gkyl_array_clear(em, 0.0);
    gkyl_array_shiftc(em, 1.0, k); // E_x = b_k in every cell
    gkyl_array_clear(rhs, 0.0);
    gkyl_dg_calc_resistive_layer_advance(up, &s->local, em, rhs);
    gkyl_dg_mul_op_range(&s->basis, 0, prod, 0, sig_modal, 0, em, &s->local);
    for (long c = 0; c < ncell; ++c) {
      const double *r = gkyl_array_cfetch(rhs, c), *pr = gkyl_array_cfetch(prod, c);
      for (int j = 0; j < nb; ++j) {
        A[(c * nb + j) * nb + k] = -r[j];
        Aw[(c * nb + j) * nb + k] = pr[j];
      }
    }
  }

  double ords[3], bq[64][64], wq[64];
  for (int q = 0; q < ntot; ++q) {
    wq[q] = 1.0;
    for (int d = 0, rem = q; d < cdim; ++d) {
      int stride = (int)pow(nq, cdim - 1 - d), i = rem / stride;
      rem -= i * stride;
      ords[d] = gkyl_gauss_ordinates[nq][i];
      wq[q] *= gkyl_gauss_weights[nq][i];
    }
    s->basis.eval(ords, bq[q]);
  }
  double min_lam = 1e300, min_lam_weak = 1e300, max_err = 0.0, max_sig = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->local);
  while (gkyl_range_iter_next(&iter)) {
    long c = gkyl_range_idx(&s->local_ext, iter.idx);
    const double *sg = gkyl_array_cfetch(sigma, c);
    for (int j = 0; j < nb; ++j) {
      for (int k = 0; k < nb; ++k) {
        double ref = 0.0;
        for (int q = 0; q < ntot; ++q) {
          ref += wq[q] * sg[q] * bq[q][j] * bq[q][k];
        }
        max_err = fmax(max_err, fabs(A[(c * nb + j) * nb + k] - ref));
      }
    }
    for (int q = 0; q < ntot; ++q) {
      max_sig = fmax(max_sig, sg[q]);
    }
    min_lam = fmin(min_lam, min_eig(nb, &A[c * nb * nb]));
    min_lam_weak = fmin(min_lam_weak, min_eig(nb, &Aw[c * nb * nb]));
  }
  TEST_CHECK(max_err < 1e-13);
  TEST_CHECK(min_lam > -1e-12);
  TEST_CHECK(min_lam_weak < -0.1);
  TEST_MSG(
    "p=%d %s: quadrature form error %g, smallest eigenvalue %g (largest sigma at a point %g); weak "
    "product smallest eigenvalue %g",
    poly_order, s->basis.id, max_err, min_lam, max_sig, min_lam_weak
  );
  gkyl_free(A);
  gkyl_free(Aw);
  gkyl_dg_calc_resistive_layer_release(up);
  gkyl_array_release(sigma);
  gkyl_array_release(sig_modal);
  gkyl_array_release(em);
  gkyl_array_release(rhs);
  gkyl_array_release(prod);
}

void
test_sample_all()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    for (int p = 1; p <= 3; ++p) {
      test_sample(cdim, p);
    }
  }
}
void
test_uniform_ser()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    for (int p = 1; p <= 3; ++p) {
      test_uniform(cdim, p, false);
    }
  }
}
void
test_uniform_tensor()
{
  for (int cdim = 1; cdim <= 2; ++cdim) {
    for (int p = 1; p <= 2; ++p) {
      test_uniform(cdim, p, true);
    }
  }
}
void
test_p1_weak_product_all()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    test_p1_weak_product(cdim, false);
    test_p1_weak_product(cdim, true);
  }
}
void
test_dissipation_p2()
{
  test_dissipation(2, false);
  test_dissipation(2, true);
}
void
test_dissipation_p3()
{
  test_dissipation(3, false);
}

TEST_LIST = {
  {"sample", test_sample_all},
  {"uniform_ser", test_uniform_ser},
  {"uniform_tensor", test_uniform_tensor},
  {"p1_weak_product", test_p1_weak_product_all},
  {"dissipation_p2", test_dissipation_p2},
  {"dissipation_p3", test_dissipation_p3},
  {NULL, NULL}
};
