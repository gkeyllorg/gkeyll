#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_dist_type_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_dist_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>
#include <stdio.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Distribution context and analytical callbacks.
struct dist_test_ctx {
  int vdim;
  int num_temp;
  double n;
  double drift[3];
  double T_over_m[3];
  double kappa;
};

static void
eval_density(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct dist_test_ctx *params = ctx;
  fout[0] = params->n;
}

static void
eval_drift(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct dist_test_ctx *params = ctx;
  for (int d = 0; d < params->vdim; ++d)
    fout[d] = params->drift[d];
}

static void
eval_temperature(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct dist_test_ctx *params = ctx;
  for (int d = 0; d < params->num_temp; ++d)
    fout[d] = params->T_over_m[d];
}

// LTE uses one temperature; multi-temperature Maxwellian uses one per direction.
static void
eval_maxwellian_dist(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct dist_test_ctx *params = ctx;
  double value = params->n;
  for (int d = 0; d < params->vdim; ++d) {
    double temperature = params->T_over_m[params->num_temp == 1 ? 0 : d];
    double relative_velocity = xn[d+1] - params->drift[d];
    value *= exp(-0.5*relative_velocity*relative_velocity/temperature)
      /sqrt(2.0*GKYL_PI*temperature);
  }
  fout[0] = value;
}

static void
eval_kappa_dist(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct dist_test_ctx *params = ctx;
  double alpha = (2.0*params->kappa-params->vdim)*params->T_over_m[0];
  double relative_speed_sq = 0.0;
  for (int d = 0; d < params->vdim; ++d) {
    double relative_velocity = xn[d+1] - params->drift[d];
    relative_speed_sq += relative_velocity*relative_velocity;
  }
  double norm = params->n*exp(
    lgamma(params->kappa+1.0) - lgamma(params->kappa+1.0-0.5*params->vdim)
    - 0.5*params->vdim*log(GKYL_PI*alpha));
  fout[0] = norm*pow(1.0+relative_speed_sq/alpha, -(params->kappa+1.0));
}

// Test configuration
typedef void (*test_evalf_t)(double t, const double *xn,
  double *restrict fout, void *ctx);

struct dist_test_case {
  enum gkyl_distribution_projection dist_id;

  int vdim;
  int num_mom;
  double lower[4];
  double upper[4];
  int cells[4];
  void *ctx;
  double *param;
  test_evalf_t eval_dist;
  test_evalf_t eval_density;
  test_evalf_t eval_drift;
  test_evalf_t eval_temperature;
  const char *name;
};

// Grid, basis, and range objects
struct dist_test_geometry {
  struct gkyl_rect_grid phase_grid;
  struct gkyl_rect_grid conf_grid;
  struct gkyl_rect_grid vel_grid;

  struct gkyl_basis phase_basis;
  struct gkyl_basis conf_basis;
  struct gkyl_basis vel_basis;

  struct gkyl_range phase_local;
  struct gkyl_range phase_local_ext;
  struct gkyl_range conf_local;
  struct gkyl_range conf_local_ext;
  struct gkyl_range vel_local;
  struct gkyl_range vel_local_ext;
};

static void
init_test_geometry(struct dist_test_geometry *g, const struct dist_test_case *tc, int poly_order)
{
  int cdim = 1;
  int vdim = tc->vdim;
  int pdim = cdim + vdim;

  double conf_lower[] = {tc->lower[0]};
  double conf_upper[] = {tc->upper[0]};
  int conf_cells[] = {tc->cells[0]};

  double vel_lower[3], vel_upper[3];
  int vel_cells[3];

  for (int d=0; d<vdim; ++d) {
    vel_lower[d] = tc->lower[d+1];
    vel_upper[d] = tc->upper[d+1];
    vel_cells[d] = tc->cells[d+1];
  }

  // Grids
  gkyl_rect_grid_init(&g->phase_grid, pdim, tc->lower, tc->upper, tc->cells);
  gkyl_rect_grid_init(&g->conf_grid, cdim, conf_lower, conf_upper, conf_cells);
  gkyl_rect_grid_init(&g->vel_grid, vdim, vel_lower, vel_upper, vel_cells);

  // Basis
  gkyl_cart_modal_serendip(&g->phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&g->conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&g->vel_basis, vdim, poly_order);

  // Ranges
  int conf_ghost[] = {1};
  int vel_ghost[] = {0, 0, 0};
  int phase_ghost[] = {1, 0, 0, 0};
  gkyl_create_grid_ranges(&g->conf_grid, conf_ghost, &g->conf_local_ext, &g->conf_local);
  gkyl_create_grid_ranges(&g->vel_grid, vel_ghost, &g->vel_local_ext, &g->vel_local);
  gkyl_create_grid_ranges(&g->phase_grid, phase_ghost, &g->phase_local_ext, &g->phase_local);
}

// Project analytical density, drift, and temperature into the target moment array.
static void
project_target_moments(const struct dist_test_geometry *g,
  const struct dist_test_case *tc, int poly_order,
  struct gkyl_array *moms_target)
{
  int vdim = tc->vdim;
  int num_temp = tc->num_mom - vdim - 1;
  int ncb = g->conf_basis.num_basis;

  struct gkyl_array *density = mkarr(ncb, g->conf_local_ext.volume);
  struct gkyl_array *drift = mkarr(vdim*ncb, g->conf_local_ext.volume);
  struct gkyl_array *temperature = mkarr(num_temp*ncb, g->conf_local_ext.volume);

  gkyl_proj_on_basis *proj_density = gkyl_proj_on_basis_new(&g->conf_grid, &g->conf_basis,
      poly_order+1, 1, tc->eval_density, tc->ctx);
  gkyl_proj_on_basis *proj_drift = gkyl_proj_on_basis_new(&g->conf_grid, &g->conf_basis,
      poly_order+1, vdim, tc->eval_drift, tc->ctx);
  gkyl_proj_on_basis *proj_temperature = gkyl_proj_on_basis_new(&g->conf_grid, &g->conf_basis,
      poly_order+1, num_temp, tc->eval_temperature, tc->ctx);

  gkyl_proj_on_basis_advance(proj_density, 0.0, &g->conf_local, density);
  gkyl_proj_on_basis_advance(proj_drift, 0.0, &g->conf_local, drift);
  gkyl_proj_on_basis_advance(proj_temperature, 0.0, &g->conf_local, temperature);

  gkyl_array_set_offset_range(moms_target, 1.0, density, 0, &g->conf_local);
  gkyl_array_set_offset_range(moms_target, 1.0, drift, ncb, &g->conf_local);
  gkyl_array_set_offset_range(moms_target, 1.0, temperature, (vdim+1)*ncb, &g->conf_local);

  gkyl_array_release(density);
  gkyl_array_release(drift);
  gkyl_array_release(temperature);
  gkyl_proj_on_basis_release(proj_density);
  gkyl_proj_on_basis_release(proj_drift);
  gkyl_proj_on_basis_release(proj_temperature);
}

// General distribution projection test.
static void
test_distribution_projection(int poly_order, struct dist_test_case *tc, double tol)
{
  struct dist_test_geometry g;
  init_test_geometry(&g, tc, poly_order);
  struct dist_test_ctx *ctx = tc->ctx;
  int vdim = tc->vdim;
  int ncb = g.conf_basis.num_basis;

  struct gkyl_array *moms = mkarr(tc->num_mom*ncb, g.conf_local_ext.volume);
  project_target_moments(&g, tc, poly_order, moms);

  struct gkyl_dist_proj_type *dist_proj = NULL;
  switch (tc->dist_id) {
    case GKYL_DIST_TYPE_LTE:
      dist_proj = gkyl_vlasov_lte_new(1, 1+vdim, poly_order,
        ncb, g.phase_basis.num_basis, GKYL_MODEL_DEFAULT);
      break;
    case GKYL_DIST_TYPE_MULTI_TEMP_MAX:
      dist_proj = gkyl_vlasov_multi_temp_max_new(1, 1+vdim, poly_order,
        ncb, g.phase_basis.num_basis, GKYL_MODEL_DEFAULT);
      break;
    case GKYL_DIST_TYPE_KAPPA:
      dist_proj = gkyl_vlasov_kappa_new(1, 1+vdim, poly_order,
        ncb, g.phase_basis.num_basis, GKYL_MODEL_DEFAULT);
      break;
    default:
      break;
  }

  TEST_CHECK(dist_proj != NULL);
  TEST_CHECK(dist_proj->dist_id == tc->dist_id);
  TEST_CHECK(gkyl_dist_proj_type_num_mom(dist_proj) == tc->num_mom);

  // Kernel tests at the drift velocity and at a thermal offset.
  double xc_test[4] = {0.5}, dx_test[4];
  int idx_test[4] = {1, 1, 1, 1};
  double moms_test[7] = {ctx->n};
  for (int d = 0; d < 1+vdim; ++d)
    dx_test[d] = g.phase_grid.dx[d];
  for (int d = 0; d < vdim; ++d) {
    xc_test[d+1] = ctx->drift[d];
    moms_test[d+1] = ctx->drift[d];
  }
  for (int d = 0; d < ctx->num_temp; ++d)
    moms_test[1+vdim+d] = ctx->T_over_m[d];

  double amplitude = ctx->n;
  double alpha = (2.0*ctx->kappa-vdim)*ctx->T_over_m[0];
  if (tc->dist_id == GKYL_DIST_TYPE_KAPPA) {
    amplitude *= exp(lgamma(ctx->kappa+1.0) - lgamma(ctx->kappa+1.0-0.5*vdim)
      - 0.5*vdim*log(GKYL_PI*alpha));
  }
  else {
    for (int d = 0; d < vdim; ++d)
      amplitude /= sqrt(2.0*GKYL_PI*ctx->T_over_m[ctx->num_temp == 1 ? 0 : d]);
  }

  double f_test[1];
  gkyl_dist_proj_type_calc(dist_proj, xc_test, dx_test, idx_test,
    moms_test, tc->param, f_test);
  TEST_CHECK(gkyl_compare_double(f_test[0], amplitude, 1e-14));

  double exponent = 0.0;
  double relative_speed_sq = 0.0;
  for (int d = 0; d < vdim; ++d) {
    double temperature = ctx->T_over_m[ctx->num_temp == 1 ? 0 : d];
    double offset_sq = (d+1)*temperature;
    xc_test[d+1] += sqrt(offset_sq);
    exponent += 0.5*(d+1);
    relative_speed_sq += offset_sq;
  }
  double f_expected = tc->dist_id == GKYL_DIST_TYPE_KAPPA
    ? amplitude*pow(1.0+relative_speed_sq/alpha, -(ctx->kappa+1.0))
    : amplitude*exp(-exponent);
  gkyl_dist_proj_type_calc(dist_proj, xc_test, dx_test, idx_test,
    moms_test, tc->param, f_test);
  TEST_CHECK(gkyl_compare_double(f_test[0], f_expected, 1e-14));

  struct gkyl_vlasov_dist_proj_on_basis_inp inp = {
    .dist_proj = dist_proj,
    .param = tc->param,
    .phase_grid = &g.phase_grid,
    .vel_grid = &g.vel_grid,
    .conf_basis = &g.conf_basis,
    .vel_basis = &g.vel_basis,
    .phase_basis = &g.phase_basis,
    .conf_range = &g.conf_local,
    .conf_range_ext = &g.conf_local_ext,
    .vel_range = &g.vel_local,
    .phase_range = &g.phase_local,
    .model_id = GKYL_MODEL_DEFAULT,
    .use_gpu = false,
    .quad_type = GKYL_GAUSS_QUAD,
  };
  gkyl_vlasov_dist_proj_on_basis *proj_dist = gkyl_vlasov_dist_proj_on_basis_inew(&inp);
  TEST_CHECK(proj_dist != NULL);

  struct gkyl_array *distf = mkarr(g.phase_basis.num_basis, g.phase_local_ext.volume);
  struct gkyl_array *distf_ref = mkarr(g.phase_basis.num_basis, g.phase_local_ext.volume);
  gkyl_vlasov_dist_proj_on_basis_moments_advance(proj_dist, &g.conf_local, moms);
  gkyl_vlasov_dist_proj_on_basis_advance(proj_dist, &g.phase_local, &g.conf_local, distf);

  gkyl_proj_on_basis *proj_ref = gkyl_proj_on_basis_new(&g.phase_grid, &g.phase_basis,
    poly_order+1, 1, tc->eval_dist, tc->ctx);
  gkyl_proj_on_basis_advance(proj_ref, 0.0, &g.phase_local, distf_ref);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &g.phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&g.phase_local_ext, iter.idx);
    const double *f = gkyl_array_cfetch(distf, idx);
    const double *f_ref = gkyl_array_cfetch(distf_ref, idx);
    for (int k = 0; k < g.phase_basis.num_basis; ++k)
      TEST_CHECK(gkyl_compare_double(f[k], f_ref[k], tol));
  }

  if (tc->dist_id == GKYL_DIST_TYPE_KAPPA) {
    // Test parameter updates using the same projection updater.
    ctx->kappa = 6.0;
    tc->param[0] = ctx->kappa;
    struct gkyl_array *distf_updated = mkarr(g.phase_basis.num_basis, g.phase_local_ext.volume);
    struct gkyl_array *distf_ref_updated = mkarr(g.phase_basis.num_basis, g.phase_local_ext.volume);
    gkyl_vlasov_dist_proj_on_basis_advance(proj_dist, &g.phase_local, &g.conf_local, distf_updated);
    gkyl_proj_on_basis_advance(proj_ref, 0.0, &g.phase_local, distf_ref_updated);

    bool changed = false;
    gkyl_range_iter_init(&iter, &g.phase_local);
    while (gkyl_range_iter_next(&iter)) {
      long idx = gkyl_range_idx(&g.phase_local_ext, iter.idx);
      const double *f_old = gkyl_array_cfetch(distf, idx);
      const double *f_new = gkyl_array_cfetch(distf_updated, idx);
      const double *f_ref = gkyl_array_cfetch(distf_ref_updated, idx);
      for (int k = 0; k < g.phase_basis.num_basis; ++k) {
        TEST_CHECK(gkyl_compare_double(f_new[k], f_ref[k], vdim == 3 ? tol : 1e-12));
        if (fabs(f_new[k]-f_old[k]) > 1e-10)
          changed = true;
      }
    }
    TEST_CHECK(changed);
    gkyl_array_release(distf_updated);
    gkyl_array_release(distf_ref_updated);
  }

  char fname[1024];
  sprintf(fname, "ctest_dist_proj_on_basis_%s_p%d.gkyl", tc->name, poly_order);
  gkyl_grid_sub_array_write(&g.phase_grid, &g.phase_local, 0, distf, fname);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ref);
  gkyl_proj_on_basis_release(proj_ref);
  gkyl_vlasov_dist_proj_on_basis_release(proj_dist);
  gkyl_dist_proj_type_release(dist_proj);
}

void
test_1x1v_lte(int poly_order)
{
  struct dist_test_ctx ctx = {
    .vdim = 1,
    .num_temp = 1,
    .n = 1.0,
    .drift = {0.5},
    .T_over_m = {1.0},
  };

  struct dist_test_case tc = {
    .dist_id = GKYL_DIST_TYPE_LTE,
    .vdim = 1,
    .num_mom = 3,
    .lower = {0.1, -8.0},
    .upper = {1.0, 8.0},
    .cells = {2, 32},
    .ctx = &ctx,
    .param = NULL,
    .eval_dist = eval_maxwellian_dist,
    .eval_density = eval_density,
    .eval_drift = eval_drift,
    .eval_temperature = eval_temperature,
    .name = "lte_1x1v",
  };

  test_distribution_projection(poly_order, &tc, 1e-12);
}

void
test_1x3v_lte(int poly_order)
{
  struct dist_test_ctx ctx = {
    .vdim = 3,
    .num_temp = 1,
    .n = 1.25,
    .drift = {0.5, -0.25, 0.75},
    .T_over_m = {1.0, 1.0, 1.0},
  };

  struct dist_test_case tc = {
    .dist_id = GKYL_DIST_TYPE_LTE,
    .vdim = 3,
    .num_mom = 5,
    .lower = {0.1, -8.0, -8.0, -8.0},
    .upper = {1.0, 8.0, 8.0, 8.0},
    .cells = {2, 16, 16, 16},
    .ctx = &ctx,
    .param = NULL,
    .eval_dist = eval_maxwellian_dist,
    .eval_density = eval_density,
    .eval_drift = eval_drift,
    .eval_temperature = eval_temperature,
    .name = "lte_1x3v",
  };

  test_distribution_projection(poly_order, &tc, 1e-10);
}

void
test_1x2v_multi_temp_max(int poly_order)
{
  struct dist_test_ctx ctx = {
    .vdim = 2,
    .num_temp = 2,
    .n = 1.25,
    .drift = {0.5, -0.25},
    .T_over_m = {1.0, 2.0},
  };

  struct dist_test_case tc = {
    .dist_id = GKYL_DIST_TYPE_MULTI_TEMP_MAX,
    .vdim = 2,
    .num_mom = 5,
    .lower = {0.1, -8.0, -8.0},
    .upper = {1.0, 8.0, 8.0},
    .cells = {2, 16, 16},
    .ctx = &ctx,
    .param = NULL,
    .eval_dist = eval_maxwellian_dist,
    .eval_density = eval_density,
    .eval_drift = eval_drift,
    .eval_temperature = eval_temperature,
    .name = "multi_temp_max_1x2v",
  };

  test_distribution_projection(poly_order, &tc, 1e-12);
}

void
test_1x3v_multi_temp_max(int poly_order)
{
  struct dist_test_ctx ctx = {
    .vdim = 3,
    .num_temp = 3,
    .n = 1.25,
    .drift = {0.5, -0.25, 0.75},
    .T_over_m = {1.0, 2.0, 0.75},
  };

  struct dist_test_case tc = {
    .dist_id = GKYL_DIST_TYPE_MULTI_TEMP_MAX,
    .vdim = 3,
    .num_mom = 7,
    .lower = {0.1, -8.0, -8.0, -8.0},
    .upper = {1.0, 8.0, 8.0, 8.0},
    .cells = {2, 16, 16, 16},
    .ctx = &ctx,
    .param = NULL,
    .eval_dist = eval_maxwellian_dist,
    .eval_density = eval_density,
    .eval_drift = eval_drift,
    .eval_temperature = eval_temperature,
    .name = "multi_temp_max_1x3v",
  };

  test_distribution_projection(poly_order, &tc, 1e-10);
}

void
test_1x1v_kappa(int poly_order)
{
  struct dist_test_ctx ctx = {
    .vdim = 1,
    .num_temp = 1,
    .n = 1.0,
    .drift = {0.5},
    .T_over_m = {1.0},
    .kappa = 7.0,
  };

  double param[] = {ctx.kappa};

  struct dist_test_case tc = {
    .dist_id = GKYL_DIST_TYPE_KAPPA,
    .vdim = 1,
    .num_mom = 3,
    .lower = {0.1, -8.0},
    .upper = {1.0, 8.0},
    .cells = {2, 32},
    .ctx = &ctx,
    .param = param,
    .eval_dist = eval_kappa_dist,
    .eval_density = eval_density,
    .eval_drift = eval_drift,
    .eval_temperature = eval_temperature,
    .name = "kappa_1x1v",
  };

  test_distribution_projection(poly_order, &tc, 1e-11);
}

void
test_1x3v_kappa(int poly_order)
{
  struct dist_test_ctx ctx = {
    .vdim = 3,
    .num_temp = 1,
    .n = 1.25,
    .drift = {0.5, -0.25, 0.75},
    .T_over_m = {1.0},
    .kappa = 10.0,
  };

  double param[] = {ctx.kappa};

  struct dist_test_case tc = {
    .dist_id = GKYL_DIST_TYPE_KAPPA,
    .vdim = 3,
    .num_mom = 5,
    .lower = {0.1, -8.0, -8.0, -8.0},
    .upper = {1.0, 8.0, 8.0, 8.0},
    .cells = {2, 16, 16, 16},
    .ctx = &ctx,
    .param = param,
    .eval_dist = eval_kappa_dist,
    .eval_density = eval_density,
    .eval_drift = eval_drift,
    .eval_temperature = eval_temperature,
    .name = "kappa_1x3v",
  };

  test_distribution_projection(poly_order, &tc, 1e-10);
}

void
test_1x1v_lte_p2(void){test_1x1v_lte(2);}

void
test_1x3v_lte_p2(void){test_1x3v_lte(2);}

void
test_1x2v_multi_temp_max_p2(void){test_1x2v_multi_temp_max(2);}

void
test_1x3v_multi_temp_max_p2(void){test_1x3v_multi_temp_max(2);}

void
test_1x1v_kappa_p2(void){test_1x1v_kappa(2);}

void
test_1x3v_kappa_p2(void){test_1x3v_kappa(2);}

TEST_LIST = {
  {"test_1x1v_lte_p2", test_1x1v_lte_p2},
  {"test_1x3v_lte_p2", test_1x3v_lte_p2},
  {"test_1x2v_multi_temp_max_p2", test_1x2v_multi_temp_max_p2},
  {"test_1x3v_multi_temp_max_p2", test_1x3v_multi_temp_max_p2},
  {"test_1x1v_kappa_p2", test_1x1v_kappa_p2},
  {"test_1x3v_kappa_p2", test_1x3v_kappa_p2},
  {NULL, NULL},
};
