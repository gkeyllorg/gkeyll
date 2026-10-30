/* Test the fem_poisson_perp_adiabatic updater, which solves the 2x/3x field equation with an
 * adiabatic response and a flux-surface average (FSA),
 *   -nabla_perp . (eps nabla_perp phi) + [K - alpha nabla_perp . (eps nabla_perp)] (phi - <phi>) = rho,
 * with the Woodbury identity (alpha = 0 without FLR effects). Slab with J=1, eps=1, Dirichlet
 * x, periodic y (3x) or z (2x), and the manufactured solution
 *   phi = A(x) + B(x) cos(k s),  A = sin(pi x) + phi_lo + (phi_up - phi_lo) x,  B = sin(2 pi x),
 * with s = y (3x) or s = z (2x), so that <phi> = A and
 *   rho = pi^2 sin(pi x) + ((1+alpha) (4 pi^2 + kperp^2) + K) B cos(k s),
 * where kperp^2 = k^2 in 3x and 0 in 2x. Nonzero Dirichlet values phi_lo, phi_up exercise the
 * removal of the rhs = 0 solution from the zonal response. For contrast, a Boltzmann response
 * without the FSA (a plain Helmholtz solve, alpha = 0) damps the zonal part to
 * pi^2/(pi^2+K) sin(pi x).
 */

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_fem_poisson_perp_adiabatic.h>
#include <gkyl_flux_surf_avg.h>
#include <gkyl_null_comm.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <math.h>
#include <stdlib.h>

static struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                   gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static bool
verbose(void)
{
  return getenv("TEST_VERBOSE") != NULL;
}

struct sol_ctx {
  int cdim;
  double K; // e^2 n0/Te.
  double lap_fac; // Factor alpha of the Laplacian term in the response.
  double k; // Wavenumber of the non-zonal part.
  double phi_lo, phi_up; // Dirichlet values in x.
  bool damped; // Whether the zonal part is damped (Boltzmann response, no FSA).
};

static double
zonal_func(double x, const struct sol_ctx *c)
{
  return sin(M_PI * x) + c->phi_lo + (c->phi_up - c->phi_lo) * x;
}

static double
nonzonal_func(double x)
{
  return sin(2.0 * M_PI * x);
}

static void
eval_phi(double t, const double *xn, double *fout, void *ctx)
{
  struct sol_ctx *c = ctx;
  double x = xn[0], s = xn[1];
  fout[0] = zonal_func(x, c) + nonzonal_func(x) * cos(c->k * s);
}

static void
eval_rho(double t, const double *xn, double *fout, void *ctx)
{
  struct sol_ctx *c = ctx;
  double x = xn[0], s = xn[1];
  double kperpSq = c->cdim == 3 ? c->k * c->k : 0.0;
  fout[0] =
    M_PI * M_PI * sin(M_PI * x) +
    ((1.0 + c->lap_fac) * (4.0 * M_PI * M_PI + kperpSq) + c->K) * nonzonal_func(x) * cos(c->k * s);
}

static void
eval_zonal(double t, const double *xn, double *fout, void *ctx)
{
  struct sol_ctx *c = ctx;
  fout[0] = c->damped ? M_PI * M_PI / (M_PI * M_PI + c->K) * sin(M_PI * xn[0]) :
                        zonal_func(xn[0], c);
}

static double
error_L2norm(
  struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_basis basis,
  const struct gkyl_array *field1, const struct gkyl_array *field2
)
{
  // L2 norm of the difference between 2 host fields (field2 may be NULL).
  struct gkyl_array *diff = gkyl_array_new(GKYL_DOUBLE, field1->ncomp, field1->size);
  gkyl_array_copy(diff, field1);
  if (field2) {
    gkyl_array_accumulate(diff, -1.0, field2);
  }

  struct gkyl_array *l2_cell = gkyl_array_new(GKYL_DOUBLE, 1, field1->size);
  gkyl_dg_calc_l2_range(&basis, 0, l2_cell, 0, diff, range);
  gkyl_array_scale_range(l2_cell, grid.cellVolume, &range);

  double l2[1];
  gkyl_array_reduce_range(l2, l2_cell, GKYL_SUM, &range);

  gkyl_array_release(diff);
  gkyl_array_release(l2_cell);
  return sqrt(l2[0]);
}

// Grids, ranges and coefficients of the slab problem, and the FSA updater.
struct slab {
  bool use_gpu;
  struct gkyl_rect_grid grid, grid_x;
  struct gkyl_basis basis, basis_x;
  struct gkyl_range local, local_ext, local_x, local_x_ext;
  struct gkyl_poisson_bc bcs;
  struct gkyl_array *jac, *eps, *kJ;
  struct gkyl_comm *comm;
  struct gkyl_flux_surf_avg *fsa;
};

static void
slab_new(struct slab *sl, const struct sol_ctx *ctx, const int *cells, bool use_gpu)
{
  int cdim = ctx->cdim, poly_order = 1;
  sl->use_gpu = use_gpu;

  double lower[] = {0.0, 0.0, -M_PI}, upper[] = {1.0, 1.0, M_PI};
  if (cdim == 2) {
    lower[1] = -M_PI;
    upper[1] = M_PI;
  }
  gkyl_rect_grid_init(&sl->grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&sl->grid_x, 1, lower, upper, cells);
  gkyl_cart_modal_serendip(&sl->basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&sl->basis_x, 1, poly_order);

  int ghost[] = {1, 1, 1};
  gkyl_create_grid_ranges(&sl->grid, ghost, &sl->local_ext, &sl->local);
  gkyl_range_init(&sl->local_x_ext, 1, &sl->local_ext.lower[0], &sl->local_ext.upper[0]);
  gkyl_sub_range_init(&sl->local_x, &sl->local_x_ext, &sl->local.lower[0], &sl->local.upper[0]);

  // J = 1, eps = identity and K constant.
  int nb = sl->basis.num_basis;
  double dg0norm = pow(sqrt(2.0), cdim);
  sl->jac = mkarr(use_gpu, nb, sl->local_ext.volume);
  gkyl_array_shiftc(sl->jac, dg0norm, 0);
  int epsnum = cdim == 3 ? 3 : 1;
  sl->eps = mkarr(use_gpu, epsnum * nb, sl->local_ext.volume);
  gkyl_array_shiftc(sl->eps, dg0norm, 0 * nb);
  if (cdim == 3) {
    gkyl_array_shiftc(sl->eps, dg0norm, 2 * nb);
  }
  sl->kJ = mkarr(use_gpu, nb, sl->local_ext.volume);
  gkyl_array_shiftc(sl->kJ, ctx->K * dg0norm, 0);

  // Dirichlet in x, periodic in y (3x).
  sl->bcs = (struct gkyl_poisson_bc){};
  sl->bcs.lo_type[0] = GKYL_POISSON_DIRICHLET;
  sl->bcs.up_type[0] = GKYL_POISSON_DIRICHLET;
  sl->bcs.lo_value[0].v[0] = ctx->phi_lo;
  sl->bcs.up_value[0].v[0] = ctx->phi_up;
  if (cdim == 3) {
    sl->bcs.lo_type[1] = GKYL_POISSON_PERIODIC;
    sl->bcs.up_type[1] = GKYL_POISSON_PERIODIC;
  }

  sl->comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = use_gpu});
  sl->fsa = gkyl_flux_surf_avg_inew(&(struct gkyl_flux_surf_avg_inp){
    .grid = &sl->grid,
    .basis = sl->basis,
    .basis_x = sl->basis_x,
    .local = &sl->local,
    .local_ext = &sl->local_ext,
    .local_x = &sl->local_x,
    .local_x_ext = &sl->local_x_ext,
    .jacobgeo = sl->jac,
    .is_deflated = true,
    .comm = sl->comm,
    .use_gpu = use_gpu,
  });
}

static void
slab_release(struct slab *sl)
{
  gkyl_flux_surf_avg_release(sl->fsa);
  gkyl_comm_release(sl->comm);
  gkyl_array_release(sl->jac);
  gkyl_array_release(sl->eps);
  gkyl_array_release(sl->kJ);
}

// Solve with the adiabatic response (use_fsa, Laplacian term lap_fac) or the Boltzmann response
// (lap_fac = 0), and compute the L2 errors of phi and of its zonal part, and the L2 norm of the
// zonal part.
static void
solve_case(
  int cdim, const int *cells, bool use_fsa, double lap_fac, double phi_lo, double phi_up,
  bool use_gpu, double *err_phi, double *err_zonal, double *norm_zonal
)
{
  int poly_order = 1;
  struct sol_ctx ctx = {
    .cdim = cdim,
    .K = 20.0,
    .lap_fac = lap_fac,
    .k = cdim == 3 ? 2.0 * M_PI : 1.0,
    .phi_lo = phi_lo,
    .phi_up = phi_up,
    .damped = !use_fsa,
  };
  struct slab sl;
  slab_new(&sl, &ctx, cells, use_gpu);
  int nb = sl.basis.num_basis, nb_x = sl.basis_x.num_basis;

  struct gkyl_array *rho = mkarr(use_gpu, nb, sl.local_ext.volume);
  struct gkyl_array *phi = mkarr(use_gpu, nb, sl.local_ext.volume);
  struct gkyl_array *psi = mkarr(use_gpu, nb_x, sl.local_x_ext.volume);
  struct gkyl_array *rho_ho = use_gpu ? mkarr(false, nb, sl.local_ext.volume) :
                                        gkyl_array_acquire(rho);
  struct gkyl_array *phi_ho = use_gpu ? mkarr(false, nb, sl.local_ext.volume) :
                                        gkyl_array_acquire(phi);
  struct gkyl_array *psi_ho = mkarr(false, nb_x, sl.local_x_ext.volume);
  struct gkyl_array *phi_sol = mkarr(false, nb, sl.local_ext.volume);
  struct gkyl_array *psi_sol = mkarr(false, nb_x, sl.local_x_ext.volume);

  gkyl_proj_on_basis *proj =
    gkyl_proj_on_basis_new(&sl.grid, &sl.basis, poly_order + 1, 1, eval_rho, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &sl.local, rho_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(rho, rho_ho);

  proj = gkyl_proj_on_basis_new(&sl.grid, &sl.basis, 2 * (poly_order + 1), 1, eval_phi, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &sl.local, phi_sol);
  gkyl_proj_on_basis_release(proj);

  proj = gkyl_proj_on_basis_new(&sl.grid_x, &sl.basis_x, 2 * (poly_order + 1), 1, eval_zonal, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &sl.local_x, psi_sol);
  gkyl_proj_on_basis_release(proj);

  if (use_fsa) {
    struct gkyl_fem_poisson_perp_adiabatic *solver =
      gkyl_fem_poisson_perp_adiabatic_inew(&(struct gkyl_fem_poisson_perp_adiabatic_inp){
        .grid = &sl.grid,
        .basis = sl.basis,
        .local = &sl.local,
        .local_ext = &sl.local_ext,
        .basis_x = sl.basis_x,
        .local_x = &sl.local_x,
        .local_x_ext = &sl.local_x_ext,
        .fsa = sl.fsa,
        .bcs = &sl.bcs,
        .bias_line_list = NULL,
        .epsilon = sl.eps,
        .kJ = sl.kJ,
        .lap_fac = lap_fac,
        .use_gpu = use_gpu,
      });
    double cond = gkyl_fem_poisson_perp_adiabatic_cond(solver);
    TEST_CHECK(isfinite(cond) && cond >= 1.0);
    TEST_MSG("cond_1(I-G) = %g", cond);
    gkyl_fem_poisson_perp_adiabatic_advance(solver, rho, phi);
    gkyl_fem_poisson_perp_adiabatic_release(solver);
  } else {
    // Boltzmann response: (L + K) phi = rho, with kSq = -K.
    struct gkyl_array *kSq = mkarr(use_gpu, nb, sl.local_ext.volume);
    gkyl_array_set(kSq, -1.0, sl.kJ);
    struct gkyl_fem_poisson_perp *helmholtz =
      gkyl_fem_poisson_perp_new(&sl.local, &sl.grid, sl.basis, &sl.bcs, NULL, sl.eps, kSq, use_gpu);
    gkyl_fem_poisson_perp_set_rhs(helmholtz, rho);
    gkyl_fem_poisson_perp_solve(helmholtz, phi);
    gkyl_fem_poisson_perp_release(helmholtz);
    gkyl_array_release(kSq);
  }
  gkyl_array_copy(phi_ho, phi);

  // Error in phi (only meaningful with the FSA) and in its zonal part.
  *err_phi = error_L2norm(sl.grid, sl.local, sl.basis, phi_ho, phi_sol);
  gkyl_flux_surf_avg_advance(sl.fsa, phi, psi);
  gkyl_array_copy(psi_ho, psi);
  *err_zonal = error_L2norm(sl.grid_x, sl.local_x, sl.basis_x, psi_ho, psi_sol);
  *norm_zonal = error_L2norm(sl.grid_x, sl.local_x, sl.basis_x, psi_sol, NULL);

  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(psi);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(psi_ho);
  gkyl_array_release(phi_sol);
  gkyl_array_release(psi_sol);
  slab_release(&sl);
}

// Check 2nd order convergence of phi and its zonal part with the FSA.
static void
check_adiabatic(int cdim, double lap_fac, double phi_lo, double phi_up, bool use_gpu)
{
  int cells[3] = {16, 16, 2};
  int cells2[3] = {32, 32, 2};
  double err_phi[2], err_zonal[2], norm_zonal[2];

  solve_case(
    cdim, cells, true, lap_fac, phi_lo, phi_up, use_gpu, &err_phi[0], &err_zonal[0], &norm_zonal[0]
  );
  solve_case(
    cdim, cells2, true, lap_fac, phi_lo, phi_up, use_gpu, &err_phi[1], &err_zonal[1], &norm_zonal[1]
  );
  double ratio_phi = err_phi[0] / err_phi[1], ratio_zonal = err_zonal[0] / err_zonal[1];
  if (verbose()) {
    printf(
      "\nFSA (lap_fac=%g, phi_lo=%g, phi_up=%g): phi err %.4e %.4e (ratio %.3f), zonal err %.4e "
      "%.4e (ratio %.3f), zonal norm %.4e\n",
      lap_fac, phi_lo, phi_up, err_phi[0], err_phi[1], ratio_phi, err_zonal[0], err_zonal[1],
      ratio_zonal, norm_zonal[1]
    );
  }
  TEST_CHECK(ratio_phi >= 3.5);
  TEST_MSG("phi L2 error ratio (2x refinement) = %.4f", ratio_phi);
  TEST_CHECK(ratio_zonal >= 3.5);
  TEST_MSG("zonal L2 error ratio (2x refinement) = %.4f", ratio_zonal);
  TEST_CHECK(err_zonal[1] < 1e-2 * norm_zonal[1]);
  TEST_MSG("zonal relative L2 error = %.4e", err_zonal[1] / norm_zonal[1]);
}

static void
test_posson_perp_adiabatic(int cdim, bool use_gpu)
{
  // Adiabatic response, with homogeneous and nonzero Dirichlet values, without and with the
  // Laplacian term of FLR effects.
  check_adiabatic(cdim, 0.0, 0.0, 0.0, use_gpu);
  check_adiabatic(cdim, 0.0, 0.5, -0.3, use_gpu);
  check_adiabatic(cdim, 0.6, 0.0, 0.0, use_gpu);
  check_adiabatic(cdim, 0.6, 0.5, -0.3, use_gpu);

  // Without the FSA the zonal part is damped by pi^2/(pi^2+K).
  int cells[3] = {16, 16, 2};
  int cells2[3] = {32, 32, 2};
  double err_phi[2], err_zonal[2], norm_zonal[2];
  solve_case(cdim, cells, false, 0.0, 0.0, 0.0, use_gpu, &err_phi[0], &err_zonal[0], &norm_zonal[0]);
  solve_case(
    cdim, cells2, false, 0.0, 0.0, 0.0, use_gpu, &err_phi[1], &err_zonal[1], &norm_zonal[1]
  );
  double ratio_zonal = err_zonal[0] / err_zonal[1];
  if (verbose()) {
    printf(
      "no FSA: zonal err %.4e %.4e (ratio %.3f), zonal norm %.4e\n", err_zonal[0], err_zonal[1],
      ratio_zonal, norm_zonal[1]
    );
  }
  TEST_CHECK(ratio_zonal >= 3.5);
  TEST_MSG("damped zonal L2 error ratio (2x refinement) = %.4f", ratio_zonal);
  TEST_CHECK(err_zonal[1] < 1e-2 * norm_zonal[1]);
  TEST_MSG("damped zonal relative L2 error = %.4e", err_zonal[1] / norm_zonal[1]);
}

void
test_posson_perp_adiabatic_2x_p1_ho()
{
  test_posson_perp_adiabatic(2, false);
}
void
test_posson_perp_adiabatic_3x_p1_ho()
{
  test_posson_perp_adiabatic(3, false);
}

#ifdef GKYL_HAVE_CUDA
void
test_posson_perp_adiabatic_2x_p1_dev()
{
  test_posson_perp_adiabatic(2, true);
}
void
test_posson_perp_adiabatic_3x_p1_dev()
{
  test_posson_perp_adiabatic(3, true);
}
#endif

TEST_LIST = {
  {"test_posson_perp_adiabatic_2x_p1_ho", test_posson_perp_adiabatic_2x_p1_ho},
  {"test_posson_perp_adiabatic_3x_p1_ho", test_posson_perp_adiabatic_3x_p1_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_posson_perp_adiabatic_2x_p1_dev", test_posson_perp_adiabatic_2x_p1_dev},
  {"test_posson_perp_adiabatic_3x_p1_dev", test_posson_perp_adiabatic_3x_p1_dev},
#endif
  {NULL, NULL}
};
