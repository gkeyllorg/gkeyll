// Test gkyl_gk_dens_mu_cells, which computes the density of a gyrokinetic
// species in each mu cell, against the analytic result for a distribution
// that is linear in each coordinate.
//
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_dens_mu_cells.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>

// Slope of the distribution function along mu.
static const double mu_slope = 0.4;

static double
conf_profile(int cdim, const double *xn)
{
  return 1.0 + 0.5 * xn[0] + (cdim > 1 ? 0.3 * xn[1] : 0.0);
}

static void
eval_distf(double t, const double *xn, double *fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double vpar = xn[cdim], mu = xn[cdim + 1];
  fout[0] = conf_profile(cdim, xn) * (1.0 + 0.2 * vpar) * (1.0 + mu_slope * mu);
}

static void
eval_conf_profile(double t, const double *xn, double *fout, void *ctx)
{
  int cdim = *(int *)ctx;
  fout[0] = conf_profile(cdim, xn);
}

static void
test_dens_mu_cells(int cdim, bool use_gpu)
{
  int pdim = cdim + 2;
  double mass = 1.7;
  double lower[] = {0.0, 0.0, 0.0, 0.0}, upper[] = {1.0, 2.0, 0.0, 0.0};
  int cells[] = {4, 3, 0, 0};
  // Velocity space: vpar symmetric about 0, and mu.
  lower[cdim] = -2.0;
  upper[cdim] = 2.0;
  cells[cdim] = 6;
  lower[cdim + 1] = 0.0;
  upper[cdim + 1] = 3.0;
  cells[cdim + 1] = 5;
  int num_mu_cells = cells[cdim + 1];

  struct gkyl_rect_grid conf_grid, phase_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_basis conf_basis, phase_basis;
  gkyl_cart_modal_serendip(&conf_basis, cdim, 1);
  gkyl_cart_modal_gkhybrid(&phase_basis, cdim, 2);
  int conf_ghost[] = {1, 1}, phase_ghost[] = {1, 1, 0, 0};
  phase_ghost[cdim] = 0;
  struct gkyl_range conf_local, conf_local_ext, phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  struct gkyl_array *distf =
    gkyl_array_new(GKYL_DOUBLE, phase_basis.num_basis, phase_local_ext.volume);
  struct gkyl_array *profile =
    gkyl_array_new(GKYL_DOUBLE, conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *dens =
    gkyl_array_new(GKYL_DOUBLE, num_mu_cells * conf_basis.num_basis, conf_local_ext.volume);

  gkyl_proj_on_basis *proj_distf =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, 3, 1, eval_distf, &cdim);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &phase_local, distf);
  gkyl_proj_on_basis_release(proj_distf);
  gkyl_eval_on_nodes *proj_profile =
    gkyl_eval_on_nodes_new(&conf_grid, &conf_basis, 1, eval_conf_profile, &cdim);
  gkyl_eval_on_nodes_advance(proj_profile, 0.0, &conf_local, profile);
  gkyl_eval_on_nodes_release(proj_profile);

  struct gkyl_gk_dens_mu_cells *dens_up =
    gkyl_gk_dens_mu_cells_new(&phase_grid, &conf_basis, &phase_basis, mass, use_gpu);
  if (use_gpu) {
    // Compute the density on the device and copy it to the host for the checks.
    struct gkyl_array *distf_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, distf->ncomp, distf->size);
    struct gkyl_array *dens_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, dens->ncomp, dens->size);
    gkyl_array_copy(distf_dev, distf);
    gkyl_gk_dens_mu_cells_advance(dens_up, &phase_local, &conf_local, distf_dev, dens_dev);
    gkyl_array_copy(dens, dens_dev);
    gkyl_array_release(distf_dev);
    gkyl_array_release(dens_dev);
  } else {
    gkyl_gk_dens_mu_cells_advance(dens_up, &phase_local, &conf_local, distf, dens);
  }

  // Integral of (1+mu_slope*mu) over each mu cell. The part linear in vpar
  // integrates to zero.
  double dmu = phase_grid.dx[cdim + 1];
  double vpar_length = upper[cdim] - lower[cdim];
  double cell_weight[16] = {0.0};
  for (int j = 0; j < num_mu_cells; j++) {
    double mu_lo = lower[cdim + 1] + j * dmu;
    cell_weight[j] = dmu + mu_slope * (mu_lo * dmu + 0.5 * dmu * dmu);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &conf_local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&conf_local, iter.idx);
    const double *profile_d = gkyl_array_cfetch(profile, linidx);
    const double *dens_d = gkyl_array_cfetch(dens, linidx);
    for (int n = 0; n < num_mu_cells; n++) {
      double fac = 2.0 * M_PI / mass * vpar_length * cell_weight[n];
      for (int k = 0; k < conf_basis.num_basis; k++) {
        double expected = fac * profile_d[k];
        TEST_CHECK(gkyl_compare_double(expected, dens_d[n * conf_basis.num_basis + k], 1e-12));
        TEST_MSG(
          "mu cell=%d cell=%ld coeff=%d: expected %.12e, got %.12e", n, linidx, k, expected,
          dens_d[n * conf_basis.num_basis + k]
        );
      }
    }
  }

  gkyl_gk_dens_mu_cells_release(dens_up);
  gkyl_array_release(distf);
  gkyl_array_release(profile);
  gkyl_array_release(dens);
}

void
test_1x2v(void)
{
  test_dens_mu_cells(1, false);
}

void
test_2x2v(void)
{
  test_dens_mu_cells(2, false);
}

#ifdef GKYL_HAVE_CUDA
void
test_1x2v_cu(void)
{
  test_dens_mu_cells(1, true);
}

void
test_2x2v_cu(void)
{
  test_dens_mu_cells(2, true);
}
#endif

TEST_LIST = {
  {"test_1x2v", test_1x2v},
  {"test_2x2v", test_2x2v},
#ifdef GKYL_HAVE_CUDA
  {"test_1x2v_cu", test_1x2v_cu},
  {"test_2x2v_cu", test_2x2v_cu},
#endif
  {NULL, NULL}
};
