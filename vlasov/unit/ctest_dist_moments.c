#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_dist_moments.h>
#include <gkyl_util.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Density
void 
eval_density(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;  
}

// Drift velocity
void
eval_drift_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 0.5;
  // double x = xn[0];
  // fout[0] = 3 * (1/cosh(x/0.2));
}

// T/m
void
eval_temperature_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

struct lte_test_ctx {double n; double drift; double T_over_m;};

void
eval_lte_1x1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct lte_test_ctx *lte = ctx;

  double v = xn[1];
  double dv = v-lte->drift;

  fout[0] = lte->n/sqrt(2.0*GKYL_PI*lte->T_over_m)*exp(-dv*dv/(2.0*lte->T_over_m));
}

void
test_1x1v_lte(int poly_order)
{
  double lower[] = {0.1, -8.0}, upper[] = {1.0, 8.0};
  int cells[] = {2, 32};
  int cdim = 1, vdim = 1;
  int pdim = cdim+vdim;

  double conf_lower[] = {lower[0]}, conf_upper[] = {upper[0]};
  double vel_lower[] = {lower[1]}, vel_upper[] = {upper[1]};
  int conf_cells[] = {cells[0]};
  int vel_cells[] = {cells[1]};

  // Grids
  struct gkyl_rect_grid phase_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, conf_lower, conf_upper, conf_cells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, vel_lower, vel_upper, vel_cells);

  // Basis
  struct gkyl_basis phase_basis, conf_basis, vel_basis;
  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&vel_basis, vdim, poly_order);

  // Ranges
  int conf_ghost[] = {1};
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[] = {0};
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int phase_ghost[] = {1, 0};
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  // Reference LTE distribution
  struct lte_test_ctx lte_ctx = {.n = 1.0, .drift = 0.5, .T_over_m = 1.0,};

  struct gkyl_array *distf = mkarr(phase_basis.num_basis, phase_local_ext.volume);

  gkyl_proj_on_basis *proj_dist = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    poly_order+1, 1, eval_lte_1x1v, &lte_ctx);

  gkyl_proj_on_basis_advance(proj_dist, 0.0, &phase_local, distf);

  // Generalized distribution moments updater
  struct gkyl_vlasov_dist_moments_inp inp = {
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &conf_basis,
    .vel_basis = &vel_basis,
    .phase_basis = &phase_basis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_local_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .dist_id = GKYL_DIST_TYPE_LTE,
    .model_id = GKYL_MODEL_DEFAULT,
    .mass = 1.0,
    .use_gpu = false,
  };

  gkyl_vlasov_dist_moments *dist_moms = gkyl_vlasov_dist_moments_inew(&inp);
  TEST_CHECK(dist_moms != NULL);

  // Output moments: (n, V_drift, T/m)
  struct gkyl_array *moms = mkarr((vdim+2)*conf_basis.num_basis,
    conf_local_ext.volume);

  gkyl_vlasov_dist_moments_advance(dist_moms, &phase_local, &conf_local,
    distf, moms);

  // Reference moments
  struct gkyl_array *density_ref = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *drift_ref = mkarr(vdim*conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *temperature_ref = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *moms_ref = mkarr((vdim+2)*conf_basis.num_basis,
    conf_local_ext.volume);

  gkyl_proj_on_basis *proj_density = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, 1, eval_density, NULL);
  gkyl_proj_on_basis *proj_drift = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, vdim, eval_drift_1v, NULL);
  gkyl_proj_on_basis *proj_temperature = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, 1, eval_temperature_1v, NULL);

  gkyl_proj_on_basis_advance(proj_density, 0.0, &conf_local, density_ref);
  gkyl_proj_on_basis_advance(proj_drift, 0.0, &conf_local, drift_ref);
  gkyl_proj_on_basis_advance(proj_temperature, 0.0, &conf_local, temperature_ref);

  gkyl_array_set_offset_range(moms_ref, 1.0, density_ref,
    0*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms_ref, 1.0, drift_ref,
    1*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms_ref, 1.0, temperature_ref,
    (vdim+1)*conf_basis.num_basis, &conf_local);

  // Compare calculated moments to reference moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &conf_local);

  while (gkyl_range_iter_next(&iter)) {

    long idx = gkyl_range_idx(&conf_local_ext, iter.idx);
    const double *m = gkyl_array_cfetch(moms, idx);
    const double *m_ref = gkyl_array_cfetch(moms_ref, idx);

    for (int k=0; k<(vdim+2)*conf_basis.num_basis; ++k) {
      TEST_CHECK(gkyl_compare_double(m[k], m_ref[k], 1e-10));
    }
  }

  // Test density only moment calculation
  struct gkyl_array *density = mkarr(conf_basis.num_basis, conf_local_ext.volume);

  gkyl_vlasov_dist_density_moment_advance(dist_moms, &phase_local, &conf_local,
    distf, density);

  gkyl_range_iter_init(&iter, &conf_local);

  while (gkyl_range_iter_next(&iter)) {

    long idx = gkyl_range_idx(&conf_local_ext, iter.idx);
    const double *n = gkyl_array_cfetch(density, idx);
    const double *n_ref = gkyl_array_cfetch(density_ref, idx);

    for (int k=0; k<conf_basis.num_basis; ++k) {
      TEST_CHECK(gkyl_compare_double(n[k], n_ref[k], 1e-10));
    }
  }

  char fname[1024];
  sprintf(fname, "ctest_vlasov_dist_moments_lte_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&conf_grid, &conf_local, 0, moms, fname);

  gkyl_array_release(distf);
  gkyl_array_release(moms);
  gkyl_array_release(moms_ref);
  gkyl_array_release(density);
  gkyl_array_release(density_ref);
  gkyl_array_release(drift_ref);
  gkyl_array_release(temperature_ref);
  gkyl_proj_on_basis_release(proj_dist);
  gkyl_proj_on_basis_release(proj_density);
  gkyl_proj_on_basis_release(proj_drift);
  gkyl_proj_on_basis_release(proj_temperature);
  gkyl_vlasov_dist_moments_release(dist_moms);
}

void
test_1x1v_lte_p2(void) {test_1x1v_lte(2);}

TEST_LIST = {
  {"test_1x1v_lte_p2", test_1x1v_lte_p2},
  {NULL, NULL},
};