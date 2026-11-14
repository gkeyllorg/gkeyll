// Tests of the direct nodal magnetic field unit vector (dg_calc_pkpm_bvar) and the
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_ten_moment.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Identity configuration-space position map (the updaters always take one).
static struct gkyl_vlasov_position_map *
identity_pos_map(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *clocal, const struct gkyl_range *clocal_ext,
  const struct gkyl_basis *cbasis
)
{
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM] = {0};
  return gkyl_vlasov_position_map_new(cgrid, clocal, clocal_ext, cbasis, inp, false);
}

// Smooth, non-trivial field with all three components varying in every direction.
static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.5 * cos(x);
  fout[1] = -0.3 * sin(y);
  fout[2] = 0.1 * cos(z);
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.2 * y;
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.1 * x;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

static void
eval_moms(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * z; // rho
  fout[1] = 1.0 + 0.3 * cos(x); // p_par
  fout[2] = 0.7 + 0.25 * cos(y) * sin(x) + 0.1 * z * z; // p_perp
  fout[3] = 0.0; // M1
}

// Fill the 2^cdim Gauss-Legendre nodes of the reference cell.
static int
gl2_nodes(int cdim, double nodes[8][3])
{
  int node_lower[3] = {0, 0, 0}, node_upper[3] = {1, 1, 1};
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, cdim, node_lower, node_upper);
  struct gkyl_range_iter node_iter;
  gkyl_range_iter_init(&node_iter, &node_range);
  while (gkyl_range_iter_next(&node_iter)) {
    long n = gkyl_range_idx(&node_range, node_iter.idx);
    for (int d = 0; d < cdim; ++d) {
      nodes[n][d] = gkyl_gauss_ordinates_2[node_iter.idx[d]];
    }
  }
  return node_range.volume;
}

static void
test_bvar(int cdim)
{
  int poly_order = 1;
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {16, 6, 4};
  int ghost[3] = {1, 1, 1};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  int nb = basis.num_basis;

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &basis);

  struct gkyl_array *em = mkarr(8 * nb, local_ext.volume);
  struct gkyl_array *bvar = mkarr(3 * nb, local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, em);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_dg_calc_pkpm_bvar *up = gkyl_dg_calc_pkpm_bvar_new(&grid, &basis, pm, false);
  gkyl_dg_calc_pkpm_bvar_advance(up, &local_ext, em, bvar);
  TEST_CHECK(gkyl_dg_calc_pkpm_bvar_num_null_b(up) == 0);

  // Legacy unit vector for a closeness check where its fallback is inactive.
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(1.0, 0.0, 0.0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);
  struct gkyl_dg_calc_em_vars *calc_legacy =
    gkyl_dg_calc_em_vars_new(&grid, &basis, &local_ext, maxwell, geom, 0.0, 0, false);
  int nsurf = nb / (poly_order + 1);
  struct gkyl_array *bvar_legacy = mkarr(9 * nb, local_ext.volume);
  struct gkyl_array *bvar_surf_legacy = mkarr(2 * cdim * 4 * nsurf, local_ext.volume);
  struct gkyl_array *cell_avg_magB2 = gkyl_array_new(GKYL_INT, 1, local_ext.volume);
  gkyl_dg_calc_em_vars_advance(calc_legacy, em, cell_avg_magB2, bvar_legacy, bvar_surf_legacy);

  double nodes[8][3];
  int nn = gl2_nodes(cdim, nodes);
  double max_legacy_diff = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *em_d = gkyl_array_cfetch(em, loc);
    const double *b_d = gkyl_array_cfetch(bvar, loc);
    const int *fallback = gkyl_array_cfetch(cell_avg_magB2, loc);
    const double *bl_d = gkyl_array_cfetch(bvar_legacy, loc);

    for (int n = 0; n < nn; ++n) {
      double B[3], b[3];
      for (int i = 0; i < 3; ++i) {
        B[i] = basis.eval_expand(nodes[n], &em_d[(3 + i) * nb]);
        b[i] = basis.eval_expand(nodes[n], &b_d[i * nb]);
      }
      double magB = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
      // b is the interpolant of B/|B| at the construction nodes: exact there.
      for (int i = 0; i < 3; ++i) {
        TEST_CHECK(fabs(b[i] - B[i] / magB) < 1e-13);
        TEST_MSG("b[%d] = %.17g, B/|B| = %.17g", i, b[i], B[i] / magB);
      }
      TEST_CHECK(fabs(b[0] * b[0] + b[1] * b[1] + b[2] * b[2] - 1.0) < 1e-13);
    }

    if (fallback[0] == 0) {
      for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < nb; ++k) {
          double diff = fabs(b_d[i * nb + k] - bl_d[i * nb + k]);
          max_legacy_diff = diff > max_legacy_diff ? diff : max_legacy_diff;
        }
      }
    }
  }
  // The legacy vector is sqrt-with-sign of a projected bb; both are second-order
  // accurate representations of the same smooth field, so they agree to the
  // projection error of the grid but not to round-off.
  TEST_CHECK(max_legacy_diff < 2.0e-2);
  TEST_MSG("max |b_new - b_legacy| = %g", max_legacy_diff);

  // Null-field policy: B = 0 in one cell gives b = 0 there and counts every node.
  int idx0[3] = {local.lower[0], cdim > 1 ? local.lower[1] : 0, cdim > 2 ? local.lower[2] : 0};
  long loc0 = gkyl_range_idx(&local_ext, idx0);
  double *em0 = gkyl_array_fetch(em, loc0);
  for (int k = 3 * nb; k < 6 * nb; ++k) {
    em0[k] = 0.0;
  }
  gkyl_dg_calc_pkpm_bvar_advance(up, &local, em, bvar);
  TEST_CHECK(gkyl_dg_calc_pkpm_bvar_num_null_b(up) == nn);
  const double *b0 = gkyl_array_cfetch(bvar, loc0);
  for (int k = 0; k < 3 * nb; ++k) {
    TEST_CHECK(b0[k] == 0.0);
  }

  gkyl_array_release(cell_avg_magB2);
  gkyl_array_release(bvar_surf_legacy);
  gkyl_array_release(bvar_legacy);
  gkyl_dg_calc_em_vars_release(calc_legacy);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(maxwell);
  gkyl_dg_calc_pkpm_bvar_release(up);
  gkyl_vlasov_position_map_release(pm);
  gkyl_array_release(bvar);
  gkyl_array_release(em);
}

// Time derivative of the field of eval_em (its magnetic part is what d_t B/|B| reads).
static void
eval_dt_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
  fout[3] = 0.3 * cos(x + y);
  fout[4] = -0.7 + 0.2 * sin(z) * cos(x);
  fout[5] = 0.5 * sin(y) + 0.1 * x * z;
}

// d_t B/|B| at the Gauss-Legendre nodes of b: (a) the nodal quotient of d_t B by |B| there;
// (b) the product the force kernels form, (I - b b) . (d_t B/|B|) with the p=1 b, is at those
// nodes the time derivative of b itself (centred difference of the b updater along d_t B).
static void
test_dtB_over_magB(int cdim)
{
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {8, 4, 3}, ghost[3] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);
  int nb = basis.num_basis;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &basis);

  struct gkyl_array *em = mkarr(8 * nb, local_ext.volume), *dt_em = mkarr(8 * nb, local_ext.volume);
  struct gkyl_array *out = mkarr(3 * nb, local_ext.volume), *bvar = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *em_p = mkarr(8 * nb, local_ext.volume), *em_m = mkarr(8 * nb, local_ext.volume);
  struct gkyl_array *b_p = mkarr(3 * nb, local_ext.volume), *b_m = mkarr(3 * nb, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&grid, &basis, 8, 8, eval_dt_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, dt_em);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_dg_calc_pkpm_bvar *up = gkyl_dg_calc_pkpm_bvar_new(&grid, &basis, pm, false);
  TEST_CHECK(gkyl_dg_calc_pkpm_bvar_dtB_over_magB(up, &local, em, dt_em, out) == 0);
  gkyl_dg_calc_pkpm_bvar_advance(up, &local, em, bvar);
  double eps = 1e-6;
  gkyl_array_set(em_p, 1.0, em);
  gkyl_array_accumulate(em_p, eps, dt_em);
  gkyl_array_set(em_m, 1.0, em);
  gkyl_array_accumulate(em_m, -eps, dt_em);
  gkyl_dg_calc_pkpm_bvar_advance(up, &local, em_p, b_p);
  gkyl_dg_calc_pkpm_bvar_advance(up, &local, em_m, b_m);

  double nodes[8][3];
  int nn = gl2_nodes(cdim, nodes);
  double err_quot = 0.0, err_dtb = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *em_d = gkyl_array_cfetch(em, loc), *dt_d = gkyl_array_cfetch(dt_em, loc);
    const double *o_d = gkyl_array_cfetch(out, loc), *b_d = gkyl_array_cfetch(bvar, loc);
    const double *bp_d = gkyl_array_cfetch(b_p, loc), *bm_d = gkyl_array_cfetch(b_m, loc);
    for (int n = 0; n < nn; ++n) {
      double B[3], dB[3], o[3], b[3], db_fd[3];
      for (int i = 0; i < 3; ++i) {
        B[i] = basis.eval_expand(nodes[n], &em_d[(3 + i) * nb]);
        dB[i] = basis.eval_expand(nodes[n], &dt_d[(3 + i) * nb]);
        o[i] = basis.eval_expand(nodes[n], &o_d[i * nb]);
        b[i] = basis.eval_expand(nodes[n], &b_d[i * nb]);
        db_fd[i] = (basis.eval_expand(nodes[n], &bp_d[i * nb]) - basis.eval_expand(nodes[n], &bm_d[i * nb])) / (2.0 * eps);
      }
      double magB = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
      double bo = b[0] * o[0] + b[1] * o[1] + b[2] * o[2];
      for (int i = 0; i < 3; ++i) {
        err_quot = fmax(err_quot, fabs(o[i] - dB[i] / magB));
        err_dtb = fmax(err_dtb, fabs(o[i] - b[i] * bo - db_fd[i]));
      }
    }
  }
  TEST_CHECK(err_quot < 1e-13);
  TEST_MSG("cdim = %d: largest error of the nodal quotient %g", cdim, err_quot);
  TEST_CHECK(err_dtb < 1e-8);
  TEST_MSG("cdim = %d: largest difference of (I - bb).(d_t B/|B|) from d_t b at the nodes %g", cdim, err_dtb);

  struct gkyl_array *arrs[] = {em, dt_em, out, bvar, em_p, em_m, b_p, b_m};
  for (int i = 0; i < 8; ++i) {
    gkyl_array_release(arrs[i]);
  }
  gkyl_dg_calc_pkpm_bvar_release(up);
  gkyl_vlasov_position_map_release(pm);
}

void test_dtB_over_magB_1x() { test_dtB_over_magB(1); }
void test_dtB_over_magB_2x() { test_dtB_over_magB(2); }
void test_dtB_over_magB_3x() { test_dtB_over_magB(3); }

void test_bvar_1x() { test_bvar(1); }
void test_bvar_2x() { test_bvar(2); }
void test_bvar_3x() { test_bvar(3); }

TEST_LIST = {
  {"test_bvar_1x", test_bvar_1x},
  {"test_bvar_2x", test_bvar_2x},
  {"test_bvar_3x", test_bvar_3x},
  {"test_dtB_over_magB_1x", test_dtB_over_magB_1x},
  {"test_dtB_over_magB_2x", test_dtB_over_magB_2x},
  {"test_dtB_over_magB_3x", test_dtB_over_magB_3x},
  {NULL, NULL},
};
