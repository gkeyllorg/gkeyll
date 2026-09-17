/* Test the ingredients of the flux-surface average (FSA) used by the 2x/3x
 * adiabatic electron field solve: the (y,z) average of J*phi divided by that
 * of J (gkyl_array_average), and the extension of a 1D function of x to a
 * cdim field constant in the other directions (chained gkyl_translate_dim).
 * With J linear in z and a zonal function linear in x, all p=1 operations
 * are exact, so the checks are to round-off.
 */

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_average.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_translate_dim.h>
#include <gkyl_util.h>
#include <math.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

struct fsa_ctx {
  int cdim;
  double lower[3], upper[3];
  bool jac_linear; // J = 1+0.3z if true, J = 1 otherwise.
};

static void
eval_jac(double t, const double *xn, double *fout, void *ctx)
{
  struct fsa_ctx *c = ctx;
  fout[0] = c->jac_linear? 1.0 + 0.3*xn[c->cdim-1] : 1.0;
}

static double
zonal_func(double x)
{
  return 0.7 + 0.4*x;
}

static void
eval_zonal(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = zonal_func(xn[0]);
}

static void
eval_zonal_plus_nonzonal(double t, const double *xn, double *fout, void *ctx)
{
  // Add a function with zero average over a period in y (3x) or z (2x).
  struct fsa_ctx *c = ctx;
  int d = c->cdim == 3? 1 : 1; // Direction of the periodic variation.
  double L = c->upper[d]-c->lower[d];
  double x = xn[0], s = xn[d];
  fout[0] = zonal_func(x) + (1.0 + 0.5*x*x)*sin(2.0*M_PI*(s-c->lower[d])/L);
}

// Objects needed to compute the FSA and the 1D->cdim extension.
struct fsa_test {
  int cdim;
  bool use_gpu;
  struct gkyl_rect_grid grid, grid_x;
  struct gkyl_basis basis, basis_x, basis_xy;
  struct gkyl_range local, local_ext, local_x, local_x_ext, local_xy, local_xy_ext;
  struct gkyl_array_average *avg;
  struct gkyl_translate_dim *infl_lo, *infl_up;
  struct gkyl_array *jac, *jphi, *avg_jphi, *avg_jphi_ho, *avg_j_ho, *tmp_xy;
  gkyl_dg_bin_op_mem *div_mem;
};

static void
fsa_test_new(struct fsa_test *t, struct fsa_ctx *ctx, const int *cells, bool use_gpu)
{
  int cdim = ctx->cdim;
  int poly_order = 1;
  t->cdim = cdim;
  t->use_gpu = use_gpu;

  gkyl_rect_grid_init(&t->grid, cdim, ctx->lower, ctx->upper, cells);
  gkyl_rect_grid_init(&t->grid_x, 1, ctx->lower, ctx->upper, cells);
  gkyl_cart_modal_serendip(&t->basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&t->basis_x, 1, poly_order);
  gkyl_cart_modal_serendip(&t->basis_xy, 2, poly_order);

  int ghost[] = {1, 1, 1};
  gkyl_create_grid_ranges(&t->grid, ghost, &t->local_ext, &t->local);
  gkyl_range_init(&t->local_x_ext, 1, &t->local_ext.lower[0], &t->local_ext.upper[0]);
  gkyl_sub_range_init(&t->local_x, &t->local_x_ext, &t->local.lower[0], &t->local.upper[0]);
  gkyl_range_init(&t->local_xy_ext, 2, t->local_ext.lower, t->local_ext.upper);
  gkyl_sub_range_init(&t->local_xy, &t->local_xy_ext, t->local.lower, t->local.upper);

  int nb = t->basis.num_basis, nb_x = t->basis_x.num_basis;
  t->jac = mkarr(use_gpu, nb, t->local_ext.volume);
  t->jphi = mkarr(use_gpu, nb, t->local_ext.volume);
  t->avg_jphi = mkarr(use_gpu, nb_x, t->local_x_ext.volume);
  t->avg_jphi_ho = use_gpu? mkarr(false, nb_x, t->local_x_ext.volume) : gkyl_array_acquire(t->avg_jphi);
  t->avg_j_ho = mkarr(false, nb_x, t->local_x_ext.volume);
  t->div_mem = gkyl_dg_bin_op_mem_new(t->local_x.volume, nb_x);

  // Jacobian.
  struct gkyl_array *jac_ho = use_gpu? mkarr(false, nb, t->local_ext.volume) : gkyl_array_acquire(t->jac);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&t->grid, &t->basis, poly_order+1, 1, eval_jac, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &t->local_ext, jac_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(t->jac, jac_ho);
  gkyl_array_release(jac_ho);

  // Average over the directions other than x.
  int avg_dim[3] = {0};
  for (int d=1; d<cdim; d++) avg_dim[d] = 1;
  t->avg = gkyl_array_average_inew(&(struct gkyl_array_average_inp) {
    .grid = &t->grid,
    .basis = t->basis,
    .basis_avg = t->basis_x,
    .local = &t->local,
    .local_avg = &t->local_x,
    .local_avg_ext = &t->local_x_ext,
    .weight = NULL,
    .avg_dim = avg_dim,
    .use_gpu = use_gpu,
  });
  gkyl_array_average_advance(t->avg, t->jac, t->avg_jphi);
  gkyl_array_copy(t->avg_j_ho, t->avg_jphi);

  // 1D -> cdim extension.
  t->tmp_xy = 0;
  t->infl_up = 0;
  if (cdim == 2) {
    t->infl_lo = gkyl_translate_dim_new(1, t->basis_x, 2, t->basis, 0, GKYL_NO_EDGE, use_gpu);
  }
  else {
    t->tmp_xy = mkarr(use_gpu, t->basis_xy.num_basis, t->local_xy_ext.volume);
    t->infl_lo = gkyl_translate_dim_new(1, t->basis_x, 2, t->basis_xy, 0, GKYL_NO_EDGE, use_gpu);
    t->infl_up = gkyl_translate_dim_new(2, t->basis_xy, 3, t->basis, 0, GKYL_NO_EDGE, use_gpu);
  }
}

static void
fsa_test_release(struct fsa_test *t)
{
  gkyl_array_release(t->jac);
  gkyl_array_release(t->jphi);
  gkyl_array_release(t->avg_jphi);
  gkyl_array_release(t->avg_jphi_ho);
  gkyl_array_release(t->avg_j_ho);
  gkyl_dg_bin_op_mem_release(t->div_mem);
  gkyl_array_average_release(t->avg);
  gkyl_translate_dim_release(t->infl_lo);
  if (t->cdim == 3) {
    gkyl_translate_dim_release(t->infl_up);
    gkyl_array_release(t->tmp_xy);
  }
}

static void
fsa_test_fsa(struct fsa_test *t, const struct gkyl_array *phi, struct gkyl_array *psi_ho)
{
  // psi = int J phi dy dz / int J dy dz (host output).
  gkyl_dg_mul_op_range(&t->basis, 0, t->jphi, 0, t->jac, 0, phi, &t->local);
  gkyl_array_average_advance(t->avg, t->jphi, t->avg_jphi);
  gkyl_array_copy(t->avg_jphi_ho, t->avg_jphi);
  gkyl_dg_div_op_range(t->div_mem, &t->basis_x, 0, psi_ho, 0, t->avg_jphi_ho, 0, t->avg_j_ho, &t->local_x);
}

static void
fsa_test_inflate(struct fsa_test *t, const struct gkyl_array *psi, struct gkyl_array *out)
{
  if (t->cdim == 2) {
    gkyl_translate_dim_advance(t->infl_lo, &t->local_x, &t->local, psi, 1, out);
  }
  else {
    gkyl_translate_dim_advance(t->infl_lo, &t->local_x, &t->local_xy, psi, 1, t->tmp_xy);
    gkyl_translate_dim_advance(t->infl_up, &t->local_xy, &t->local, t->tmp_xy, 1, out);
  }
}

static double
max_abs_diff_1x(struct fsa_test *t, const struct gkyl_array *a, const struct gkyl_array *b)
{
  double err = 0.0;
  for (int i=t->local_x.lower[0]; i<=t->local_x.upper[0]; i++) {
    long lidx = gkyl_range_idx(&t->local_x, &i);
    const double *a_c = gkyl_array_cfetch(a, lidx);
    const double *b_c = gkyl_array_cfetch(b, lidx);
    for (int k=0; k<t->basis_x.num_basis; k++) err = GKYL_MAX2(err, fabs(a_c[k]-b_c[k]));
  }
  return err;
}

static void
test_fsa(int cdim, bool use_gpu)
{
  struct fsa_ctx ctx = {
    .cdim = cdim,
    .lower = {-1.0, 0.0, 0.0},
    .upper = { 1.0, 2.0, 1.0},
    .jac_linear = true,
  };
  int cells[] = {6, 8, 4};
  if (cdim == 2) {
    ctx.upper[1] = 1.0;
    cells[1] = 8;
  }
  int poly_order = 1;

  struct fsa_test t;
  fsa_test_new(&t, &ctx, cells, use_gpu);
  int nb = t.basis.num_basis, nb_x = t.basis_x.num_basis;

  struct gkyl_array *phi = mkarr(use_gpu, nb, t.local_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, nb, t.local_ext.volume) : gkyl_array_acquire(phi);
  struct gkyl_array *psi_ho = mkarr(false, nb_x, t.local_x_ext.volume);
  struct gkyl_array *psi_ref = mkarr(false, nb_x, t.local_x_ext.volume);

  // Reference: 1D projection of the zonal function.
  gkyl_proj_on_basis *proj_x = gkyl_proj_on_basis_new(&t.grid_x, &t.basis_x, poly_order+1, 1, eval_zonal, &ctx);
  gkyl_proj_on_basis_advance(proj_x, 0.0, &t.local_x, psi_ref);
  gkyl_proj_on_basis_release(proj_x);

  // 1) FSA of a function constant on flux surfaces returns it exactly.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&t.grid, &t.basis, poly_order+1, 1, eval_zonal, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &t.local, phi_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(phi, phi_ho);
  fsa_test_fsa(&t, phi, psi_ho);
  double err = max_abs_diff_1x(&t, psi_ho, psi_ref);
  TEST_CHECK( err < 1e-12 );
  TEST_MSG("FSA of zonal function: max error %g", err);

  // 2) The non-zonal part averages to zero. In 2x the periodic variation is
  // along the averaged direction, so use J=1 there for the integral to vanish.
  struct fsa_test t2 = t;
  bool use_t2 = cdim == 2;
  if (use_t2) {
    ctx.jac_linear = false;
    fsa_test_new(&t2, &ctx, cells, use_gpu);
  }
  proj = gkyl_proj_on_basis_new(&t.grid, &t.basis, poly_order+1, 1, eval_zonal_plus_nonzonal, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &t.local, phi_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(phi, phi_ho);
  fsa_test_fsa(use_t2? &t2 : &t, phi, psi_ho);
  err = max_abs_diff_1x(&t, psi_ho, psi_ref);
  TEST_CHECK( err < 1e-12 );
  TEST_MSG("FSA of zonal+nonzonal function: max error %g", err);
  if (use_t2)
    fsa_test_release(&t2);

  // 3) Extension of a 1D function, and FSA of the extension is the identity.
  struct gkyl_array *psi = mkarr(use_gpu, nb_x, t.local_x_ext.volume);
  struct gkyl_array *psi_in_ho = use_gpu? mkarr(false, nb_x, t.local_x_ext.volume) : gkyl_array_acquire(psi);
  gkyl_array_clear(psi_in_ho, 0.0);
  for (int i=t.local_x.lower[0]; i<=t.local_x.upper[0]; i++) {
    double *c = gkyl_array_fetch(psi_in_ho, gkyl_range_idx(&t.local_x, &i));
    c[0] = 1.0 + 0.1*i;
    c[1] = 0.3 - 0.05*i;
  }
  gkyl_array_copy(psi, psi_in_ho);
  fsa_test_inflate(&t, psi, phi);
  gkyl_array_copy(phi_ho, phi);

  // The constant extension scales the 1 and x coefficients by sqrt(2)^(cdim-1).
  double fac = pow(sqrt(2.0), cdim-1);
  err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &t.local);
  while (gkyl_range_iter_next(&iter)) {
    const double *p = gkyl_array_cfetch(phi_ho, gkyl_range_idx(&t.local, iter.idx));
    const double *c = gkyl_array_cfetch(psi_in_ho, gkyl_range_idx(&t.local_x, &iter.idx[0]));
    for (int k=0; k<nb; k++) {
      double ref = k<2? fac*c[k] : 0.0;
      err = GKYL_MAX2(err, fabs(p[k]-ref));
    }
  }
  TEST_CHECK( err < 1e-14 );
  TEST_MSG("Extension coefficients: max error %g", err);

  fsa_test_fsa(&t, phi, psi_ho);
  err = max_abs_diff_1x(&t, psi_ho, psi_in_ho);
  TEST_CHECK( err < 1e-12 );
  TEST_MSG("FSA of extension: max error %g", err);

  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(psi);
  gkyl_array_release(psi_in_ho);
  gkyl_array_release(psi_ho);
  gkyl_array_release(psi_ref);
  fsa_test_release(&t);
}

void test_2x_p1() { test_fsa(2, false); }
void test_3x_p1() { test_fsa(3, false); }

#ifdef GKYL_HAVE_CUDA
void gpu_test_2x_p1() { test_fsa(2, true); }
void gpu_test_3x_p1() { test_fsa(3, true); }
#endif

TEST_LIST = {
  { "test_2x_p1", test_2x_p1 },
  { "test_3x_p1", test_3x_p1 },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_2x_p1", gpu_test_2x_p1 },
  { "gpu_test_3x_p1", gpu_test_3x_p1 },
#endif
  { NULL, NULL },
};
