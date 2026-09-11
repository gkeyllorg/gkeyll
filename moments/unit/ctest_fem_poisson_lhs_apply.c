// Test gkyl_fem_poisson_lhs_apply, which applies the LHS operator
// M^{-1}*(M+K) = weak (1 - rho^2*Lap) using the assembled FEM matrices.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static struct gkyl_array*
mkarr_dev(bool on_gpu, long nc, long size)
{
#ifdef GKYL_HAVE_CUDA
  if (on_gpu) return gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
#endif
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static double
calc_l2(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_range range_ext,
  struct gkyl_basis basis, struct gkyl_array *f1, struct gkyl_array *f2)
{
  struct gkyl_array *diff = mkarr(basis.num_basis, range_ext.volume);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *f1p = gkyl_array_cfetch(f1, lidx);
    const double *f2p = f2? gkyl_array_cfetch(f2, lidx) : NULL;
    double *diffp = gkyl_array_fetch(diff, lidx);
    for (int i=0; i<basis.num_basis; i++) diffp[i] = f1p[i] - (f2p? f2p[i] : 0.0);
  }
  struct gkyl_array *l2 = mkarr(1, range_ext.volume);
  gkyl_dg_calc_l2_range(&basis, 0, l2, 0, diff, range);
  gkyl_array_scale_range(l2, grid.cellVolume, &range);
  double l2red[1];
  gkyl_array_reduce_range(l2red, l2, GKYL_SUM, &range);
  gkyl_array_release(diff);
  gkyl_array_release(l2);
  return sqrt(l2red[0]);
}

static void evalFunc_1x(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = cos(xn[0])*sin(xn[0]);
}
static void evalFunc_2x(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = cos(xn[0])*sin(xn[1])*cos(xn[1]);
}

static void
test_lhs_apply(int dim, int cells[], double rho, enum gkyl_poisson_bc_type bc_type, bool use_gpu)
{
  double lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
  for (int d=0; d<dim; d++) { lower[d] = -M_PI; upper[d] = M_PI; }
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int nghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_poisson_bc bcs = {0};
  for (int d=0; d<dim; d++) {
    bcs.lo_type[d] = bc_type;
    bcs.up_type[d] = bc_type;
  }

  // Permittivity epsilon = rho^2 (isotropic) and kSq = -1, projected source, all on host.
  double dg0norm = pow(sqrt(2.0), dim);
  int neps = 2*(dim/3)+1;
  struct gkyl_array *epsilon_ho = mkarr(neps*basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon_ho, rho*rho*dg0norm, 0*basis.num_basis);
  if (dim > 2) {
    gkyl_array_shiftc(epsilon_ho, rho*rho*dg0norm, 2*basis.num_basis);
  }
  struct gkyl_array *kSq_ho = mkarr(basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(kSq_ho, -dg0norm, 0);

  struct gkyl_array *src_ho = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1,
    dim==1? evalFunc_1x : evalFunc_2x, NULL);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, src_ho);
  gkyl_proj_on_basis_release(proj);

  // Working arrays (on device if use_gpu), populated from the host copies.
  struct gkyl_array *epsilon = mkarr_dev(use_gpu, neps*basis.num_basis, local_ext.volume);
  struct gkyl_array *kSq = mkarr_dev(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *src = mkarr_dev(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *w = mkarr_dev(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *g = mkarr_dev(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *w2 = mkarr_dev(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_array_copy(epsilon, epsilon_ho);
  gkyl_array_copy(kSq, kSq_ho);
  gkyl_array_copy(src, src_ho);

  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&local, &grid, basis, &bcs, NULL,
    epsilon, kSq, false, use_gpu);

  // Continuous field w = A^{-1}*src.
  gkyl_fem_poisson_set_rhs(poisson, src, NULL);
  gkyl_fem_poisson_solve(poisson, w);

  // Apply the LHS operator, g = (1 - rho^2*Lap) w.
  gkyl_fem_poisson_lhs_apply(poisson, w, g);

  // Exact identity: solving with the same operator inverts the apply.
  gkyl_fem_poisson_set_rhs(poisson, g, NULL);
  gkyl_fem_poisson_solve(poisson, w2);

  // Copy results to host for the L2 checks.
  struct gkyl_array *w_ho = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *g_ho = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *w2_ho = mkarr(basis.num_basis, local_ext.volume);
  gkyl_array_copy(w_ho, w);
  gkyl_array_copy(g_ho, g);
  gkyl_array_copy(w2_ho, w2);

  double err = calc_l2(grid, local, local_ext, basis, w2_ho, w_ho);
  double norm_w = calc_l2(grid, local, local_ext, basis, w_ho, NULL);
  TEST_CHECK(err < 1.0e-10*norm_w);
  TEST_MSG("round-trip L2 error = %g (|w| = %g)", err, norm_w);

  // Applying (1 - rho^2*Lap) amplifies the mode, so |g| > |w|.
  double norm_g = calc_l2(grid, local, local_ext, basis, g_ho, NULL);
  TEST_CHECK(norm_g > norm_w);
  TEST_MSG("|g| = %g, |w| = %g", norm_g, norm_w);

  gkyl_array_release(epsilon_ho);
  gkyl_array_release(kSq_ho);
  gkyl_array_release(src_ho);
  gkyl_array_release(epsilon);
  gkyl_array_release(kSq);
  gkyl_array_release(src);
  gkyl_array_release(w);
  gkyl_array_release(g);
  gkyl_array_release(w2);
  gkyl_array_release(w_ho);
  gkyl_array_release(g_ho);
  gkyl_array_release(w2_ho);
  gkyl_fem_poisson_release(poisson);
}

void test_1x_p1(void) { 
  int cells[] = {16};     
  test_lhs_apply(1, cells, 0.3, GKYL_POISSON_DIRICHLET, false);
  test_lhs_apply(1, cells, 0.3, GKYL_POISSON_PERIODIC, false);
}

void test_2x_p1(void) { int cells[] = {24, 16}; 
  test_lhs_apply(2, cells, 0.3, GKYL_POISSON_DIRICHLET, false);
  test_lhs_apply(2, cells, 0.3, GKYL_POISSON_PERIODIC, false);
 }

#ifdef GKYL_HAVE_CUDA

void gpu_test_1x_p1(void) { 
  int cells[] = {16};     
  test_lhs_apply(1, cells, 0.3, GKYL_POISSON_DIRICHLET, true);
  test_lhs_apply(1, cells, 0.3, GKYL_POISSON_PERIODIC, true);
}

void gpu_test_2x_p1(void) { 
  int cells[] = {24, 16};     
  test_lhs_apply(2, cells, 0.3, GKYL_POISSON_DIRICHLET, true);
  test_lhs_apply(2, cells, 0.3, GKYL_POISSON_PERIODIC, true);
}

#endif

TEST_LIST = {
  { "test_1x_p1", test_1x_p1 },
  { "test_2x_p1", test_2x_p1 },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_1x_p1", gpu_test_1x_p1 },
  { "gpu_test_2x_p1", gpu_test_2x_p1 },
#endif
  { NULL, NULL },
};
