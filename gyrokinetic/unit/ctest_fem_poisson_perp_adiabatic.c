/* Test the 2x/3x field solve with adiabatic electrons and a flux-surface
 * average (FSA),
 *   -nabla_perp . (eps nabla_perp phi) + K (phi - <phi>) = rho,
 * solved with the Woodbury identity used by the gyrokinetic app:
 *   phi = phi1 + H^{-1} K E psi,  H phi1 = rho,  (I - G) psi = R phi1,
 * with H = -nabla_perp.(eps nabla_perp) + K, R the FSA, E the extension of a
 * function of x to a field constant in the other directions and
 * G = R H^{-1} K E. Slab with J=1, eps=1, Dirichlet x, periodic y (3x) or
 * z (2x), and the manufactured solution
 *   phi = A(x) + B(x) cos(k s),  A = sin(pi x),  B = sin(2 pi x),
 * with s = y (3x) or s = z (2x), so that <phi> = A and
 *   rho = pi^2 A + (4 pi^2 + kperp^2 + K) B cos(k s),
 * where kperp^2 = k^2 in 3x and 0 in 2x. Without the FSA the zonal part is
 * damped to pi^2/(pi^2+K) A.
 */

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_average.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_mat.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_translate_dim.h>
#include <gkyl_util.h>
#include <math.h>
#include <stdlib.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static bool verbose(void) { return getenv("TEST_VERBOSE") != NULL; }

struct sol_ctx {
  int cdim;
  double K; // e^2 n0/Te.
  double k; // Wavenumber of the non-zonal part.
  bool damped; // Whether the zonal part is damped (no FSA).
};

static double zonal_func(double x) { return sin(M_PI*x); }
static double nonzonal_func(double x) { return sin(2.0*M_PI*x); }

static void
eval_phi(double t, const double *xn, double *fout, void *ctx)
{
  struct sol_ctx *c = ctx;
  double x = xn[0], s = xn[1];
  fout[0] = zonal_func(x) + nonzonal_func(x)*cos(c->k*s);
}

static void
eval_rho(double t, const double *xn, double *fout, void *ctx)
{
  struct sol_ctx *c = ctx;
  double x = xn[0], s = xn[1];
  double kperpSq = c->cdim == 3? c->k*c->k : 0.0;
  fout[0] = M_PI*M_PI*zonal_func(x) + (4.0*M_PI*M_PI + kperpSq + c->K)*nonzonal_func(x)*cos(c->k*s);
}

static void
eval_zonal(double t, const double *xn, double *fout, void *ctx)
{
  struct sol_ctx *c = ctx;
  double fac = c->damped? M_PI*M_PI/(M_PI*M_PI + c->K) : 1.0;
  fout[0] = fac*zonal_func(xn[0]);
}

static double
error_L2norm(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_basis basis,
  const struct gkyl_array *field1, const struct gkyl_array *field2)
{
  // L2 norm of the difference between 2 host fields (field2 may be NULL).
  struct gkyl_array *diff = gkyl_array_new(GKYL_DOUBLE, field1->ncomp, field1->size);
  gkyl_array_copy(diff, field1);
  if (field2)
    gkyl_array_accumulate(diff, -1.0, field2);

  struct gkyl_array *l2_cell = gkyl_array_new(GKYL_DOUBLE, 1, field1->size);
  gkyl_dg_calc_l2_range(&basis, 0, l2_cell, 0, diff, range);
  gkyl_array_scale_range(l2_cell, grid.cellVolume, &range);

  double l2[1];
  gkyl_array_reduce_range(l2, l2_cell, GKYL_SUM, &range);

  gkyl_array_release(diff);
  gkyl_array_release(l2_cell);
  return sqrt(l2[0]);
}

// Objects for the Woodbury solve.
struct woodbury {
  int cdim;
  bool use_gpu;
  double K;
  struct gkyl_rect_grid grid, grid_x;
  struct gkyl_basis basis, basis_x, basis_xy;
  struct gkyl_range local, local_ext, local_x, local_x_ext, local_xy, local_xy_ext;
  struct gkyl_fem_poisson_perp *solver;
  struct gkyl_array_average *avg;
  struct gkyl_translate_dim *infl_lo, *infl_up;
  struct gkyl_array *eps, *kSq, *tmp_xy, *avg_phi, *psi, *psi_ho, *phi2, *rhs2;
  int m;
  struct gkyl_mat *A, *A_lu, *rhs_m;
  gkyl_mem_buff ipiv;
};

static void
woodbury_fsa(struct woodbury *w, const struct gkyl_array *phi)
{
  // With J=1 the FSA is the plain average over the directions other than x.
  gkyl_array_average_advance(w->avg, phi, w->avg_phi);
  gkyl_array_copy(w->psi_ho, w->avg_phi);
}

static void
woodbury_inflate(struct woodbury *w, const struct gkyl_array *psi, struct gkyl_array *out)
{
  if (w->cdim == 2) {
    gkyl_translate_dim_advance(w->infl_lo, &w->local_x, &w->local, psi, 1, out);
  }
  else {
    gkyl_translate_dim_advance(w->infl_lo, &w->local_x, &w->local_xy, psi, 1, w->tmp_xy);
    gkyl_translate_dim_advance(w->infl_up, &w->local_xy, &w->local, w->tmp_xy, 1, out);
  }
}

static void
woodbury_response(struct woodbury *w, const struct gkyl_array *psi, struct gkyl_array *phi)
{
  // phi = H^{-1} K E psi.
  woodbury_inflate(w, psi, phi);
  gkyl_array_set(w->rhs2, w->K, phi);
  gkyl_fem_poisson_perp_set_rhs(w->solver, w->rhs2);
  gkyl_fem_poisson_perp_solve(w->solver, phi);
}

static void
woodbury_pack(struct woodbury *w, const struct gkyl_array *arr_ho, double *vec)
{
  int k = 0;
  for (int i=w->local_x.lower[0]; i<=w->local_x.upper[0]; i++) {
    const double *c = gkyl_array_cfetch(arr_ho, gkyl_range_idx(&w->local_x, &i));
    for (int b=0; b<w->basis_x.num_basis; b++) vec[k++] = c[b];
  }
}

static void
woodbury_unpack(struct woodbury *w, const double *vec, struct gkyl_array *arr_ho)
{
  int k = 0;
  for (int i=w->local_x.lower[0]; i<=w->local_x.upper[0]; i++) {
    double *c = gkyl_array_fetch(arr_ho, gkyl_range_idx(&w->local_x, &i));
    for (int b=0; b<w->basis_x.num_basis; b++) c[b] = vec[k++];
  }
}

static void
woodbury_new(struct woodbury *w, int cdim, const int *cells, double K, bool use_gpu)
{
  int poly_order = 1;
  w->cdim = cdim;
  w->use_gpu = use_gpu;
  w->K = K;

  double lower[] = {0.0, 0.0, -M_PI}, upper[] = {1.0, 1.0, M_PI};
  if (cdim == 2) {
    lower[1] = -M_PI;
    upper[1] = M_PI;
  }
  gkyl_rect_grid_init(&w->grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&w->grid_x, 1, lower, upper, cells);
  gkyl_cart_modal_serendip(&w->basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&w->basis_x, 1, poly_order);
  gkyl_cart_modal_serendip(&w->basis_xy, 2, poly_order);

  int ghost[] = {1, 1, 1};
  gkyl_create_grid_ranges(&w->grid, ghost, &w->local_ext, &w->local);
  gkyl_range_init(&w->local_x_ext, 1, &w->local_ext.lower[0], &w->local_ext.upper[0]);
  gkyl_sub_range_init(&w->local_x, &w->local_x_ext, &w->local.lower[0], &w->local.upper[0]);
  gkyl_range_init(&w->local_xy_ext, 2, w->local_ext.lower, w->local_ext.upper);
  gkyl_sub_range_init(&w->local_xy, &w->local_xy_ext, w->local.lower, w->local.upper);

  int nb = w->basis.num_basis, nb_x = w->basis_x.num_basis;
  double dg0norm = pow(sqrt(2.0), cdim);

  // Helmholtz solver: eps = identity, kSq = -K, Dirichlet x, periodic y.
  int epsnum = cdim == 3? 3 : 1;
  w->eps = mkarr(use_gpu, epsnum*nb, w->local_ext.volume);
  gkyl_array_shiftc(w->eps, dg0norm, 0*nb);
  if (cdim == 3)
    gkyl_array_shiftc(w->eps, dg0norm, 2*nb);
  w->kSq = mkarr(use_gpu, nb, w->local_ext.volume);
  gkyl_array_shiftc(w->kSq, -K*dg0norm, 0);

  struct gkyl_poisson_bc bcs = { };
  bcs.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bcs.up_type[0] = GKYL_POISSON_DIRICHLET;
  bcs.lo_value[0].v[0] = 0.0;
  bcs.up_value[0].v[0] = 0.0;
  if (cdim == 3) {
    bcs.lo_type[1] = GKYL_POISSON_PERIODIC;
    bcs.up_type[1] = GKYL_POISSON_PERIODIC;
  }
  w->solver = gkyl_fem_poisson_perp_new(&w->local, &w->grid, w->basis, &bcs, NULL, w->eps, w->kSq, use_gpu);

  // FSA.
  int avg_dim[3] = {0};
  for (int d=1; d<cdim; d++) avg_dim[d] = 1;
  w->avg = gkyl_array_average_inew(&(struct gkyl_array_average_inp) {
    .grid = &w->grid,
    .basis = w->basis,
    .basis_avg = w->basis_x,
    .local = &w->local,
    .local_avg = &w->local_x,
    .local_avg_ext = &w->local_x_ext,
    .weight = NULL,
    .avg_dim = avg_dim,
    .use_gpu = use_gpu,
  });
  w->avg_phi = mkarr(use_gpu, nb_x, w->local_x_ext.volume);
  w->psi = mkarr(use_gpu, nb_x, w->local_x_ext.volume);
  w->psi_ho = use_gpu? mkarr(false, nb_x, w->local_x_ext.volume) : gkyl_array_acquire(w->psi);

  // 1D -> cdim extension.
  w->tmp_xy = 0;
  w->infl_up = 0;
  if (cdim == 2) {
    w->infl_lo = gkyl_translate_dim_new(1, w->basis_x, 2, w->basis, 0, GKYL_NO_EDGE, use_gpu);
  }
  else {
    w->tmp_xy = mkarr(use_gpu, w->basis_xy.num_basis, w->local_xy_ext.volume);
    w->infl_lo = gkyl_translate_dim_new(1, w->basis_x, 2, w->basis_xy, 0, GKYL_NO_EDGE, use_gpu);
    w->infl_up = gkyl_translate_dim_new(2, w->basis_xy, 3, w->basis, 0, GKYL_NO_EDGE, use_gpu);
  }

  w->phi2 = mkarr(use_gpu, nb, w->local_ext.volume);
  w->rhs2 = mkarr(use_gpu, nb, w->local_ext.volume);

  // Zonal system A = I - G, one Helmholtz solve per column of G.
  int m = w->local_x.volume*nb_x;
  w->m = m;
  w->A = gkyl_mat_new(m, m, 0.0);
  w->A_lu = gkyl_mat_new(m, m, 0.0);
  w->rhs_m = gkyl_mat_new(m, 1, 0.0);
  w->ipiv = gkyl_mem_buff_new(sizeof(long[m]));
  double *ej = gkyl_malloc(sizeof(double[m]));
  for (int j=0; j<m; j++) {
    for (int i=0; i<m; i++) ej[i] = i==j? 1.0 : 0.0;
    woodbury_unpack(w, ej, w->psi_ho);
    gkyl_array_copy(w->psi, w->psi_ho);
    woodbury_response(w, w->psi, w->phi2);
    woodbury_fsa(w, w->phi2);
    woodbury_pack(w, w->psi_ho, gkyl_mat_get_col(w->A, j));
  }
  gkyl_free(ej);
  for (int i=0; i<m; i++)
    for (int j=0; j<m; j++)
      gkyl_mat_set(w->A, i, j, (i==j? 1.0 : 0.0) - gkyl_mat_get(w->A, i, j));
}

static void
woodbury_release(struct woodbury *w)
{
  gkyl_fem_poisson_perp_release(w->solver);
  gkyl_array_average_release(w->avg);
  gkyl_translate_dim_release(w->infl_lo);
  if (w->cdim == 3) {
    gkyl_translate_dim_release(w->infl_up);
    gkyl_array_release(w->tmp_xy);
  }
  gkyl_array_release(w->eps);
  gkyl_array_release(w->kSq);
  gkyl_array_release(w->avg_phi);
  gkyl_array_release(w->psi);
  gkyl_array_release(w->psi_ho);
  gkyl_array_release(w->phi2);
  gkyl_array_release(w->rhs2);
  gkyl_mat_release(w->A);
  gkyl_mat_release(w->A_lu);
  gkyl_mat_release(w->rhs_m);
  gkyl_mem_buff_release(w->ipiv);
}

static void
woodbury_solve(struct woodbury *w, struct gkyl_array *rho, struct gkyl_array *phi, bool use_fsa)
{
  // H phi1 = rho.
  gkyl_fem_poisson_perp_set_rhs(w->solver, rho);
  gkyl_fem_poisson_perp_solve(w->solver, phi);

  if (use_fsa) {
    // (I - G) psi = <phi1>, phi2 = H^{-1} K E psi.
    woodbury_fsa(w, phi);
    woodbury_pack(w, w->psi_ho, gkyl_mat_get_col(w->rhs_m, 0));
    gkyl_mat_copy(w->A_lu, w->A);
    bool status = gkyl_mat_linsolve_lu(w->A_lu, w->rhs_m, gkyl_mem_buff_data(w->ipiv));
    TEST_CHECK( status );
    woodbury_unpack(w, gkyl_mat_get_ccol(w->rhs_m, 0), w->psi_ho);
    gkyl_array_copy(w->psi, w->psi_ho);
    woodbury_response(w, w->psi, w->phi2);
    gkyl_array_accumulate_range(phi, 1.0, w->phi2, &w->local);
  }
}

static void
solve_case(int cdim, const int *cells, bool use_fsa, bool use_gpu, double *err_phi, double *err_zonal, double *norm_zonal)
{
  int poly_order = 1;
  double K = 20.0;
  struct woodbury w;
  woodbury_new(&w, cdim, cells, K, use_gpu);
  int nb = w.basis.num_basis, nb_x = w.basis_x.num_basis;

  struct sol_ctx ctx = {
    .cdim = cdim,
    .K = K,
    .k = cdim == 3? 2.0*M_PI : 1.0,
    .damped = !use_fsa,
  };

  struct gkyl_array *rho = mkarr(use_gpu, nb, w.local_ext.volume);
  struct gkyl_array *phi = mkarr(use_gpu, nb, w.local_ext.volume);
  struct gkyl_array *rho_ho = use_gpu? mkarr(false, nb, w.local_ext.volume) : gkyl_array_acquire(rho);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, nb, w.local_ext.volume) : gkyl_array_acquire(phi);
  struct gkyl_array *phi_sol = mkarr(false, nb, w.local_ext.volume);
  struct gkyl_array *psi_sol = mkarr(false, nb_x, w.local_x_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&w.grid, &w.basis, poly_order+1, 1, eval_rho, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &w.local, rho_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(rho, rho_ho);

  proj = gkyl_proj_on_basis_new(&w.grid, &w.basis, 2*(poly_order+1), 1, eval_phi, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &w.local, phi_sol);
  gkyl_proj_on_basis_release(proj);

  proj = gkyl_proj_on_basis_new(&w.grid_x, &w.basis_x, 2*(poly_order+1), 1, eval_zonal, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &w.local_x, psi_sol);
  gkyl_proj_on_basis_release(proj);

  woodbury_solve(&w, rho, phi, use_fsa);
  gkyl_array_copy(phi_ho, phi);

  // Error in phi (only meaningful with the FSA) and in its zonal part.
  *err_phi = error_L2norm(w.grid, w.local, w.basis, phi_ho, phi_sol);
  woodbury_fsa(&w, phi);
  *err_zonal = error_L2norm(w.grid_x, w.local_x, w.basis_x, w.psi_ho, psi_sol);
  *norm_zonal = error_L2norm(w.grid_x, w.local_x, w.basis_x, psi_sol, NULL);

  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phi_sol);
  gkyl_array_release(psi_sol);
  woodbury_release(&w);
}

static void
test_adiabatic(int cdim, bool use_gpu)
{
  int cells[3] = {16, 16, 2};
  int cells2[3] = {32, 32, 2};
  double err_phi[2], err_zonal[2], norm_zonal[2];

  // With the FSA: phi and its zonal part converge at 2nd order.
  solve_case(cdim, cells, true, use_gpu, &err_phi[0], &err_zonal[0], &norm_zonal[0]);
  solve_case(cdim, cells2, true, use_gpu, &err_phi[1], &err_zonal[1], &norm_zonal[1]);
  double ratio_phi = err_phi[0]/err_phi[1], ratio_zonal = err_zonal[0]/err_zonal[1];
  if (verbose())
    printf("\nFSA: phi err %.4e %.4e (ratio %.3f), zonal err %.4e %.4e (ratio %.3f), zonal norm %.4e\n",
      err_phi[0], err_phi[1], ratio_phi, err_zonal[0], err_zonal[1], ratio_zonal, norm_zonal[1]);
  TEST_CHECK( ratio_phi >= 3.5 );
  TEST_MSG("phi L2 error ratio (2x refinement) = %.4f", ratio_phi);
  TEST_CHECK( ratio_zonal >= 3.5 );
  TEST_MSG("zonal L2 error ratio (2x refinement) = %.4f", ratio_zonal);
  TEST_CHECK( err_zonal[1] < 1e-2*norm_zonal[1] );
  TEST_MSG("zonal relative L2 error = %.4e", err_zonal[1]/norm_zonal[1]);

  // Without the FSA the zonal part is damped by pi^2/(pi^2+K).
  solve_case(cdim, cells, false, use_gpu, &err_phi[0], &err_zonal[0], &norm_zonal[0]);
  solve_case(cdim, cells2, false, use_gpu, &err_phi[1], &err_zonal[1], &norm_zonal[1]);
  ratio_zonal = err_zonal[0]/err_zonal[1];
  if (verbose())
    printf("no FSA: zonal err %.4e %.4e (ratio %.3f), zonal norm %.4e\n",
      err_zonal[0], err_zonal[1], ratio_zonal, norm_zonal[1]);
  TEST_CHECK( ratio_zonal >= 3.5 );
  TEST_MSG("damped zonal L2 error ratio (2x refinement) = %.4f", ratio_zonal);
  TEST_CHECK( err_zonal[1] < 1e-2*norm_zonal[1] );
  TEST_MSG("damped zonal relative L2 error = %.4e", err_zonal[1]/norm_zonal[1]);
}

void test_2x_p1() { test_adiabatic(2, false); }
void test_3x_p1() { test_adiabatic(3, false); }

#ifdef GKYL_HAVE_CUDA
void gpu_test_2x_p1() { test_adiabatic(2, true); }
void gpu_test_3x_p1() { test_adiabatic(3, true); }
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
