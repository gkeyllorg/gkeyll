#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>
#include <gkyl_util.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static void
eval_flat_vierbein_test_1v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

static void
eval_flat_vierbein_gradient_test_1v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

static void
eval_flat_vierbein_test_2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
}

static void
eval_flat_vierbein_gradient_test_2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  for (int i=0; i<8; ++i) fout[i] = 0.0;
}

static void
eval_annulus_vierbein_test_2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double q_r = xn[0];

  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = q_r;
}

static void
eval_annulus_vierbein_gradient_test_2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;

  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

static void
eval_spherical_rtheta_vierbein_test_3v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1];

  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;

  fout[3] = 0.0;
  fout[4] = r;
  fout[5] = 0.0;

  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = r*sin(theta);
}

static void
eval_spherical_rtheta_vierbein_gradient_test_3v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1];

  // gradient in r
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;

  fout[3] = 0.0;
  fout[4] = 1.0;
  fout[5] = 0.0;

  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = sin(theta);

  // gradient in theta
  fout[9] = 0.0;
  fout[10] = 0.0;
  fout[11] = 0.0;

  fout[12] = 0.0;
  fout[13] = 0.0;
  fout[14] = 0.0;

  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = r*cos(theta);

  // gradient in phi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;

  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;

  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}

static void
eval_ks_rphi_vierbein_test_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];

  double rho_sq = r*r;
  double h_rr = (1.0 + 2.0*M*r/rho_sq);
  double h_pp = (rho_sq + a*a*(1.0 + 2.0*M*r/rho_sq));
  double h_rp = -a*(1.0 + 2.0*M*r/rho_sq);

  fout[0] = sqrt(h_rr);
  fout[1] = 0.0;
  fout[2] = h_rp/sqrt(h_rr);
  fout[3] = sqrt(h_pp - h_rp*h_rp/h_rr);
}

static void
eval_ks_rphi_vierbein_gradient_test_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double rho_sq = r*r;

  // gradient in r
  fout[0] = (M*(-r*r))/(pow(rho_sq, 1.5)*sqrt(rho_sq + 2.0*M*r));
  fout[1] = 0.0;
  fout[2] = -(M*a*(-r*r))/(pow(rho_sq, 1.5)*sqrt(rho_sq + 2.0*M*r));
  fout[3] = r/sqrt(rho_sq);

  // gradient in phi
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

static void
check_preset_geom(int cdim, int vdim, const double *lower, const double *upper, const int *cells,
  enum gkyl_triad_preset_geom_type preset_geom_type, evalf_t eval_vierbein,
  evalf_t eval_vierbein_gradient, void *ctx, double eps)
{
  int ndim = cdim+vdim;

  double confLower[GKYL_MAX_DIM], confUpper[GKYL_MAX_DIM];
  int confCells[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, 2);
  gkyl_cart_modal_serendip(&confBasis, cdim, 2);

  int confGhost[GKYL_MAX_DIM] = { 0 };
  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) {
    confGhost[d] = 1;
    ghost[d] = 1;
  }

  struct gkyl_range confLocal, confLocal_ext;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int num_pt_indices[3] = { 1, 6, 18 };
  int nc = confBasis.num_basis*num_pt_indices[vdim-1];
  struct gkyl_array *conf_poisson_tensor_direct = mkarr(nc, confLocal_ext.volume);
  struct gkyl_array *conf_poisson_tensor_preset = mkarr(nc, confLocal_ext.volume);

  struct gkyl_vlasov_triad_geom_inp inp_direct = {
    .use_vierbein = true,
    .use_preset_geom = false,
    .eval_vierbein = eval_vierbein,
    .eval_vierbein_gradient = eval_vierbein_gradient,
    .eval_vierbein_ctx = ctx,
    .eval_vierbein_gradient_ctx = ctx,
  };

  struct gkyl_vlasov_triad_geom_inp inp_preset = {
    .use_vierbein = true,
    .use_preset_geom = true,
    .triad_preset_geom_type = preset_geom_type,
    .eval_vierbein_ctx = ctx,
    .eval_vierbein_gradient_ctx = ctx,
  };

  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis,
    &grid, &local, basis, inp_direct, conf_poisson_tensor_direct);
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis,
    &grid, &local, basis, inp_preset, conf_poisson_tensor_preset);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *direct = gkyl_array_cfetch(conf_poisson_tensor_direct, gkyl_range_idx(&confLocal, iter.idx));
    const double *preset = gkyl_array_cfetch(conf_poisson_tensor_preset, gkyl_range_idx(&confLocal, iter.idx));

    for (int k=0; k<nc; ++k) {
      TEST_CHECK(gkyl_compare_double(preset[k], direct[k], eps));
    }
  }

  gkyl_array_release(conf_poisson_tensor_direct);
  gkyl_array_release(conf_poisson_tensor_preset);
}

static void
check_preset_vierbein_inv(int vdim, enum gkyl_triad_preset_geom_type preset_geom_type,
  const double *xn, void *ctx)
{
  evalf_t eval_vierbein = gkyl_vlasov_triad_preset_vierbein(vdim, preset_geom_type);
  evalf_t eval_vierbein_inv = gkyl_vlasov_triad_preset_vierbein_inv(vdim, preset_geom_type);

  double vierbein[GKYL_MAX_DIM*GKYL_MAX_DIM] = { 0.0 };
  double vierbein_inv[GKYL_MAX_DIM*GKYL_MAX_DIM] = { 0.0 };
  double vierbein_inv_expected[GKYL_MAX_DIM*GKYL_MAX_DIM] = { 0.0 };

  eval_vierbein(0.0, xn, vierbein, ctx);
  eval_vierbein_inv(0.0, xn, vierbein_inv, ctx);

  if (vdim == 1) {
    kernel_vierbein_inv_1v(vierbein, vierbein_inv_expected);
  }
  else if (vdim == 2) {
    kernel_vierbein_inv_2v(vierbein, vierbein_inv_expected);
  }
  else {
    kernel_vierbein_inv_3v(vierbein, vierbein_inv_expected);
  }

  for (int k=0; k<vdim*vdim; ++k) {
    TEST_CHECK(gkyl_compare_double(vierbein_inv[k], vierbein_inv_expected[k], 1e-12));
  }
}

void
test_preset_vierbein_inv()
{
  struct gkyl_triad_geom_ctx ks_ctx = {
    .mass_bh = 0.7,
    .spin_bh = 0.2,
  };

  double x_flat_1v[] = { 0.35 };
  double x_flat_2v[] = { 0.35, 0.7 };
  double x_annulus[] = { 0.35 };
  double x_spherical[] = { 1.2, 0.9 };
  double x_ks_rphi[] = { 1.2, 1.0 };

  check_preset_vierbein_inv(1, GKYL_TRIAD_FLAT, x_flat_1v, 0);
  check_preset_vierbein_inv(2, GKYL_TRIAD_FLAT, x_flat_2v, 0);
  check_preset_vierbein_inv(2, GKYL_TRIAD_ANNULUS, x_annulus, 0);
  check_preset_vierbein_inv(3, GKYL_TRIAD_SPHERICAL_RTHETA, x_spherical, 0);
  check_preset_vierbein_inv(2, GKYL_TRIAD_GR_KERR_SCHILD_RPHI, x_ks_rphi, &ks_ctx);
}

void
test_preset_1x1v_flat()
{
  double lower[] = { 0.1, -1.0 }, upper[] = { 1.0, 1.0 };
  int cells[] = { 2, 2 };

  check_preset_geom(1, 1, lower, upper, cells, GKYL_TRIAD_FLAT,
    eval_flat_vierbein_test_1v, eval_flat_vierbein_gradient_test_1v, 0, 1e-12);
}

void
test_preset_1x2v_flat()
{
  double lower[] = { 0.1, -1.0, -1.0 }, upper[] = { 1.0, 1.0, 1.0 };
  int cells[] = { 2, 2, 2 };

  check_preset_geom(1, 2, lower, upper, cells, GKYL_TRIAD_FLAT,
    eval_flat_vierbein_test_2v, eval_flat_vierbein_gradient_test_2v, 0, 1e-12);
}

void
test_preset_1x2v_annulus()
{
  double lower[] = { 0.1, -1.0, -1.0 }, upper[] = { 1.0, 1.0, 1.0 };
  int cells[] = { 2, 2, 2 };

  check_preset_geom(1, 2, lower, upper, cells, GKYL_TRIAD_ANNULUS,
    eval_annulus_vierbein_test_2v, eval_annulus_vierbein_gradient_test_2v, 0, 1e-12);
}

void
test_preset_2x2v_rphi_ks()
{
  double lower[] = { 1.2, 1.0, -1.0, -1.0 }, upper[] = { 1.201, 1.001, 1.0, 1.0 };
  int cells[] = { 2, 2, 2, 2 };
  struct gkyl_triad_geom_ctx ks_ctx = {
    .mass_bh = 0.7,
    .spin_bh = 0.2,
  };

  check_preset_geom(2, 2, lower, upper, cells, GKYL_TRIAD_GR_KERR_SCHILD_RPHI,
    eval_ks_rphi_vierbein_test_2v, eval_ks_rphi_vierbein_gradient_test_2v, &ks_ctx, 1e-8);
}

void
test_preset_2x3v_spherical_rtheta()
{
  double lower[] = { 0.5, 0.4, -1.0, -1.0, -1.0 };
  double upper[] = { 1.5, 1.2, 1.0, 1.0, 1.0 };
  int cells[] = { 2, 2, 2, 2, 2 };

  check_preset_geom(2, 3, lower, upper, cells, GKYL_TRIAD_SPHERICAL_RTHETA,
    eval_spherical_rtheta_vierbein_test_3v, eval_spherical_rtheta_vierbein_gradient_test_3v, 0, 1e-12);
}

TEST_LIST = {
  { "test_preset_vierbein_inv", test_preset_vierbein_inv },
  { "test_preset_1x1v_flat", test_preset_1x1v_flat },
  { "test_preset_1x2v_flat", test_preset_1x2v_flat },
  { "test_preset_1x2v_annulus", test_preset_1x2v_annulus },
  { "test_preset_2x2v_rphi_ks", test_preset_2x2v_rphi_ks },
  { "test_preset_2x3v_spherical_rtheta", test_preset_2x3v_spherical_rtheta },
  { NULL, NULL }
};
