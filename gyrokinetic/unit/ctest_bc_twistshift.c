#include <gkyl_deflate_zsurf.h>
#include <gkyl_translate_dim.h>
#include <gkyl_fem_proj.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_eval_on_nodes.h>
// Test the bc_twistshift orchestrator: the twist-shift BC on its own, and the
// supersample -> twist-shift -> low-pass filter -> restrict pipeline used to
// de-alias the shift along the shear direction. The sub-cell integrals the
// shift is built from are tested separately, in ctest_twistshift_dg.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_twistshift_dg.h>

#include <math.h>

static struct gkyl_array *
mkarr(bool use_gpu, int ncomp, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
                   gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

static const double ts_lower[] = {-2.0, -1.5, -3.0};
static const double ts_upper[] = {2.0, 1.5, 3.0};
static const int ts_cells[] = {8, 8, 4};
static const int ts_ndim = 3, ts_cdim = 3, ts_bc_dir = 2;

struct ts_ctx {
  double offset; // Shift at x=0.
  double shear; // Rate at which the shift varies with x (0 = no shear).
};

static void
shift_func(double t, const double *xn, double *fout, void *ctx)
{
  struct ts_ctx *tctx = ctx;
  fout[0] = tctx->offset + tctx->shear * xn[0];
}

static void
init_donor(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = exp(-y * y / (2.0 * 0.3 * 0.3)) * (1.0 + 0.5 * cos(M_PI * x));
}

// Everything needed to build and apply the BC on the test grid.
struct ts_setup {
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_range update_r; // Local range, extended in bc_dir.
  struct gkyl_range ghost_r; // Plane the BC writes.
  struct gkyl_range skin_r; // Plane periodicity copies from.
  int ghost[GKYL_MAX_DIM];
};

// Initialize the test grid, p1 basis, and interior/ghost ranges. Select the ghost
// plane at the requested edge and the opposite skin plane that supplies its data.
static void
ts_setup_init_grid(
  struct ts_setup *s, enum gkyl_edge_loc edge, const double *lower, const double *upper,
  const int *cells
)
{
  gkyl_rect_grid_init(&s->grid, ts_ndim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->basis, ts_ndim, 1);

  for (int d = 0; d < ts_ndim; d++) {
    s->ghost[d] = 1;
  }
  gkyl_create_grid_ranges(&s->grid, s->ghost, &s->local_ext, &s->local);

  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d = 0; d < ts_ndim; d++) {
    lo[d] = s->local.lower[d];
    up[d] = s->local.upper[d];
  }
  lo[ts_bc_dir] = s->local_ext.lower[ts_bc_dir];
  up[ts_bc_dir] = s->local_ext.upper[ts_bc_dir];
  gkyl_sub_range_init(&s->update_r, &s->local_ext, lo, up);

  // The BC reads and writes the same ghost plane; periodicity fills it from
  // the skin cell at the other end.
  int slo[GKYL_MAX_DIM], sup[GKYL_MAX_DIM];
  for (int d = 0; d < ts_ndim; d++) {
    slo[d] = s->local.lower[d];
    sup[d] = s->local.upper[d];
  }
  if (edge == GKYL_LOWER_EDGE) {
    gkyl_range_shorten_from_above(&s->ghost_r, &s->update_r, ts_bc_dir, s->ghost[ts_bc_dir]);
    slo[ts_bc_dir] = sup[ts_bc_dir] = s->local.upper[ts_bc_dir];
  } else {
    gkyl_range_shorten_from_below(&s->ghost_r, &s->update_r, ts_bc_dir, s->ghost[ts_bc_dir]);
    slo[ts_bc_dir] = sup[ts_bc_dir] = s->local.lower[ts_bc_dir];
  }
  gkyl_sub_range_init(&s->skin_r, &s->local_ext, slo, sup);
}

// Project the donor and apply periodicity along bc_dir.
static struct gkyl_array *
ts_donor_new(const struct ts_setup *s)
{
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, s->basis.num_basis, s->local_ext.volume);
  gkyl_array_clear(f, 0.0);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->grid,
    .basis = &s->basis,
    .num_ret_vals = 1,
    .eval = init_donor,
    .ctx = NULL,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local, f);
  gkyl_proj_on_basis_release(proj);

  gkyl_array_copy_range_to_range(f, f, &s->ghost_r, &s->skin_r);
  return f;
}

// Apply the BC to a fresh donor field and return it. Pass shift_dg to exercise
// the discretized-shift input instead of the shift function, and a negative
// half_width for the plain twist-shift without the anti-aliasing filter.
static struct gkyl_array *
ts_run(
  const struct ts_setup *s, enum gkyl_edge_loc edge, int upsample, int half_width, double cutoff,
  struct gkyl_array *shift_dg, struct ts_ctx *tctx
)
{
  struct gkyl_array *f = ts_donor_new(s);

  struct gkyl_bc_twistshift_inp inp = {
    .bc_dir = ts_bc_dir,
    .shift_dir = 1,
    .shear_dir = 0,
    .edge = edge,
    .cdim = ts_cdim,
    .bcdir_ext_update_r = &s->update_r,
    .num_ghost = s->ghost,
    .basis = &s->basis,
    .grid = &s->grid,
    .use_gpu = false,
    .upsample_factor = upsample,
    .filter_half_width = half_width,
    .filter_cutoff_wavelength = cutoff,
  };
  if (shift_dg) {
    inp.shift_dg = shift_dg;
  } else {
    inp.shift_func = shift_func;
    inp.shift_func_ctx = tctx;
  }

  struct gkyl_bc_twistshift *up = gkyl_bc_twistshift_inew(&inp);
  gkyl_bc_twistshift_advance(up, f, f);
  gkyl_bc_twistshift_release(up);
  return f;
}

// Largest difference between two fields over a range.
static double
ts_max_diff(
  const struct gkyl_array *fa, const struct gkyl_array *fb, const struct gkyl_range *rng,
  int num_basis
)
{
  double maxd = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, rng);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(rng, iter.idx);
    const double *a = gkyl_array_cfetch(fa, linidx);
    const double *b = gkyl_array_cfetch(fb, linidx);
    for (int k = 0; k < num_basis; k++) {
      maxd = GKYL_MAX2(maxd, fabs(a[k] - b[k]));
    }
  }
  return maxd;
}

void
test_bc_twistshift_plain_matches_twistshift_dg_ho(void)
{
  // With no anti-aliasing the BC must be exactly the bare twist-shift.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 0.25};
  for (int e = 0; e < 2; e++) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init_grid(&s, edge, ts_lower, ts_upper, ts_cells);

    struct gkyl_array *f_bc = ts_run(&s, edge, 1, -1, 0.0, NULL, &tctx);

    struct gkyl_array *f_ref = ts_donor_new(&s);
    struct gkyl_twistshift_dg_inp tsinp = {
      .bc_dir = ts_bc_dir,
      .shift_dir = 1,
      .shear_dir = 0,
      .edge = edge,
      .cdim = ts_cdim,
      .bcdir_ext_update_r = &s.update_r,
      .num_ghost = s.ghost,
      .basis = &s.basis,
      .grid = &s.grid,
      .shift_func = shift_func,
      .shift_func_ctx = &tctx,
      .use_gpu = false,
    };
    struct gkyl_twistshift_dg *ts = gkyl_twistshift_dg_inew(&tsinp);
    gkyl_twistshift_dg_advance(ts, f_ref, f_ref);
    gkyl_twistshift_dg_release(ts);

    TEST_CHECK(ts_max_diff(f_bc, f_ref, &s.local_ext, s.basis.num_basis) == 0.0);
    TEST_MSG(
      "edge %d: bc_twistshift differs from twistshift_dg by %.3e", e,
      ts_max_diff(f_bc, f_ref, &s.local_ext, s.basis.num_basis)
    );

    gkyl_array_release(f_bc);
    gkyl_array_release(f_ref);
  }
}

void
test_bc_twistshift_only_ghost_plane_written_ho(void)
{
  // Nothing but the ghost plane the BC fills may change, in any configuration.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 0.25};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];
  int upsample[] = {1, 1, 2, 4};
  int half_width[] = {-1, 4, 4, 4};
  double cutoff[] = {0.0, 4.0 * dx, 2.0 * dx, 2.0 * dx};

  for (int e = 0; e < 2; e++) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init_grid(&s, edge, ts_lower, ts_upper, ts_cells);

    struct gkyl_array *f_pre = ts_donor_new(&s);

    for (int c = 0; c < 4; c++) {
      struct gkyl_array *f = ts_run(&s, edge, upsample[c], half_width[c], cutoff[c], NULL, &tctx);

      // Restore the ghost plane so any remaining difference is a stray write.
      gkyl_array_copy_range(f, f_pre, &s.ghost_r);
      double maxd = ts_max_diff(f, f_pre, &s.local_ext, s.basis.num_basis);
      TEST_CHECK(maxd == 0.0);
      TEST_MSG(
        "edge %d, config %d (upsample %d, M %d): wrote %.3e outside the ghost plane", e, c,
        upsample[c], half_width[c], maxd
      );

      gkyl_array_release(f);
    }
    gkyl_array_release(f_pre);
  }
}

void
test_bc_twistshift_identity_filter_matches_plain_ho(void)
{
  // A cutoff at the Nyquist wavelength makes the kernel the identity, so
  // filtering must leave the plain twist-shift result alone.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 0.25};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];

  struct ts_setup s;
  ts_setup_init_grid(&s, GKYL_LOWER_EDGE, ts_lower, ts_upper, ts_cells);

  struct gkyl_array *f_plain = ts_run(&s, GKYL_LOWER_EDGE, 1, -1, 0.0, NULL, &tctx);
  struct gkyl_array *f_filt = ts_run(&s, GKYL_LOWER_EDGE, 1, 4, 2.0 * dx, NULL, &tctx);

  double maxd = ts_max_diff(f_plain, f_filt, &s.ghost_r, s.basis.num_basis);
  TEST_CHECK(maxd < 1.0e-13);
  TEST_MSG("identity filter changed the result by %.3e", maxd);

  gkyl_array_release(f_plain);
  gkyl_array_release(f_filt);
}

void
test_bc_twistshift_upsample_no_shear_matches_plain_ho(void)
{
  // With a uniform shift the twist-shift creates no structure along x, so the
  // supersampled field stays an exact prolongation of the coarse one and the
  // prolong/identity-filter/restrict round trip must be lossless.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 0.0};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];

  for (int upsample = 2; upsample <= 4; upsample *= 2) {
    struct ts_setup s;
    ts_setup_init_grid(&s, GKYL_LOWER_EDGE, ts_lower, ts_upper, ts_cells);

    struct gkyl_array *f_plain = ts_run(&s, GKYL_LOWER_EDGE, 1, -1, 0.0, NULL, &tctx);
    struct gkyl_array *f_up =
      ts_run(&s, GKYL_LOWER_EDGE, upsample, 4, 2.0 * dx / upsample, NULL, &tctx);

    double maxd = ts_max_diff(f_plain, f_up, &s.ghost_r, s.basis.num_basis);
    TEST_CHECK(maxd < 1.0e-12);
    TEST_MSG("upsample %d round trip changed the result by %.3e", upsample, maxd);

    gkyl_array_release(f_plain);
    gkyl_array_release(f_up);
  }
}

void
test_bc_twistshift_upsample_six_cbc_grid_ho(void)
{
  // The CBC radial extent and 48 cells give dx_coarse/dx_fine slightly above
  // six in floating point, despite an exact integer cell-count ratio.
  double lower[] = {-0.13360368663594469, -1.5, -3.0};
  double upper[] = {0.13360368663594469, 1.5, 3.0};
  int cells[] = {48, 8, 4};
  struct ts_setup s;
  ts_setup_init_grid(&s, GKYL_LOWER_EDGE, lower, upper, cells);

  struct ts_ctx tctx = {.offset = 0.75, .shear = 0.25};
  double dx = s.grid.dx[0];
  struct gkyl_array *f = ts_run(&s, GKYL_LOWER_EDGE, 6, 4, 2.0 * dx, NULL, &tctx);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    const double *fc = gkyl_array_cfetch(f, gkyl_range_idx(&s.ghost_r, iter.idx));
    for (int k = 0; k < s.basis.num_basis; k++) {
      TEST_CHECK(isfinite(fc[k]));
    }
  }
  gkyl_array_release(f);
}

void
test_bc_twistshift_shift_dg_matches_shift_func_ho(void)
{
  // A discretized shift given on the field's own grid must give the same answer
  // as the shift function, including when the shift has to be prolongated.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 0.25};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];

  for (int upsample = 1; upsample <= 2; upsample++) {
    struct ts_setup s;
    ts_setup_init_grid(&s, GKYL_LOWER_EDGE, ts_lower, ts_upper, ts_cells);

    // The 1D DG shift on the shear cells of the update range, as the app builds it.
    struct gkyl_rect_grid xgrid;
    gkyl_rect_grid_init(&xgrid, 1, &ts_lower[0], &ts_upper[0], &ts_cells[0]);
    struct gkyl_basis xbasis;
    gkyl_cart_modal_serendip(&xbasis, 1, s.basis.poly_order);
    struct gkyl_range xrng;
    gkyl_range_init(&xrng, 1, (int[]){s.update_r.lower[0]}, (int[]){s.update_r.upper[0]});

    struct gkyl_array *shift_dg = gkyl_array_new(GKYL_DOUBLE, xbasis.num_basis, xrng.volume);
    gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&xgrid, &xbasis, 1, shift_func, &tctx);
    gkyl_eval_on_nodes_advance(ev, 0.0, &xrng, shift_dg);
    gkyl_eval_on_nodes_release(ev);

    int half_width = upsample > 1 ? 4 : -1;
    double cutoff = upsample > 1 ? 2.0 * dx : 0.0;

    struct gkyl_array *f_func =
      ts_run(&s, GKYL_LOWER_EDGE, upsample, half_width, cutoff, NULL, &tctx);
    struct gkyl_array *f_dg =
      ts_run(&s, GKYL_LOWER_EDGE, upsample, half_width, cutoff, shift_dg, &tctx);

    double maxd = ts_max_diff(f_func, f_dg, &s.ghost_r, s.basis.num_basis);
    TEST_CHECK(maxd < 1.0e-12);
    TEST_MSG("upsample %d: shift_dg differs from shift_func by %.3e", upsample, maxd);

    gkyl_array_release(f_func);
    gkyl_array_release(f_dg);
    gkyl_array_release(shift_dg);
  }
}

void
test_bc_twistshift_dealiasing_smooths_shear_direction_ho(void)
{
  // Check that the anti-aliasing filter actually smooths the field along the shear direction.
  // We test that the total variation along the x-direction of the cell averages,
  // summed over the ghost plane, is smaller for the filtered result than for the plain twist-shift.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 2.9};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];

  struct ts_setup s;
  ts_setup_init_grid(&s, GKYL_LOWER_EDGE, ts_lower, ts_upper, ts_cells);

  struct gkyl_array *f_plain = ts_run(&s, GKYL_LOWER_EDGE, 1, -1, 0.0, NULL, &tctx);
  struct gkyl_array *f_deal = ts_run(&s, GKYL_LOWER_EDGE, 4, 4, 2.0 * dx, NULL, &tctx);

  double tv[2] = {0.0, 0.0};
  struct gkyl_array *fs[] = {f_plain, f_deal};
  for (int q = 0; q < 2; q++) {
    for (int j = s.ghost_r.lower[1]; j <= s.ghost_r.upper[1]; j++) {
      for (int i = s.ghost_r.lower[0]; i < s.ghost_r.upper[0]; i++) {
        int idx0[] = {i, j, s.ghost_r.lower[2]};
        int idx1[] = {i + 1, j, s.ghost_r.lower[2]};
        const double *a = gkyl_array_cfetch(fs[q], gkyl_range_idx(&s.ghost_r, idx0));
        const double *b = gkyl_array_cfetch(fs[q], gkyl_range_idx(&s.ghost_r, idx1));
        tv[q] += fabs(b[0] - a[0]);
      }
    }
  }
  TEST_CHECK(tv[1] < tv[0]);
  TEST_MSG("total variation along x: plain %.6e, de-aliased %.6e", tv[0], tv[1]);

  gkyl_array_release(f_plain);
  gkyl_array_release(f_deal);
}

void
test_bc_twistshift_conserves_particles_ho(void)
{
  // Check that the integral is preserved.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 2.9};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];
  int upsample[] = {1, 1, 2, 4};
  int half_width[] = {-1, 4, 4, 4};
  double cutoff[] = {0.0, 4.0 * dx, 2.0 * dx, 2.0 * dx};

  for (int e = 0; e < 2; e++) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init_grid(&s, edge, ts_lower, ts_upper, ts_cells);

    // Integral over the plane the BC fills.
    struct gkyl_array_integrate *integ =
      gkyl_array_integrate_new(&s.grid, &s.basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, false);

    struct gkyl_array *f_pre = ts_donor_new(&s);
    double tot_pre;
    gkyl_array_integrate_advance(integ, f_pre, 1.0, NULL, &s.ghost_r, &s.ghost_r, &tot_pre);

    for (int c = 0; c < 4; c++) {
      struct gkyl_array *f = ts_run(&s, edge, upsample[c], half_width[c], cutoff[c], NULL, &tctx);
      double tot;
      gkyl_array_integrate_advance(integ, f, 1.0, NULL, &s.ghost_r, &s.ghost_r, &tot);

      TEST_CHECK(fabs(tot - tot_pre) < 1.0e-12 * fabs(tot_pre));
      TEST_MSG(
        "edge %d, config %d (upsample %d, M %d): %.13e -> %.13e, rel change %.3e", e, c,
        upsample[c], half_width[c], tot_pre, tot, fabs(tot - tot_pre) / fabs(tot_pre)
      );

      gkyl_array_release(f);
    }
    gkyl_array_release(f_pre);
    gkyl_array_integrate_release(integ);
  }
}

static void
init_donor_ky0(double t, const double *xn, double *fout, void *ctx)
{
  // Steep y-independent (ky = 0) radial profile.
  double x = xn[0], y = xn[1];
  double amp_ky = *(double *)ctx;
  fout[0] = 1.0 + tanh((x - 0.4) / 0.25) +
            amp_ky * cos(2.0 * M_PI * y / (ts_upper[1] - ts_lower[1])) * exp(-x * x);
}

void
test_bc_twistshift_preserves_y_average_ho(void)
{
  // The twist-shift only moves the field along y at each x, so it must leave the
  // y-average (ky = 0 part) of every x cell unchanged.
  struct ts_ctx tctx = {.offset = 0.75, .shear = 2.9};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];
  int upsample[] = {1, 4, 4};
  int half_width[] = {-1, 1, 4};
  double cutoff[] = {0.0, 2.0 * dx, 2.0 * dx};
  const char *label[] = {"no filter", "default filter (4,1,2dx)", "filter (4,4,2dx)"};
  double amp_ky[] = {0.0, 0.5};
  int ky0_coeffs[] = {0, 1, 3, 5};

  for (int e = 0; e < 2; e++) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init_grid(&s, edge, ts_lower, ts_upper, ts_cells);

    for (int d = 0; d < 2; d++) {
      for (int c = 0; c < 3; c++) {
        // Donor, with periodicity applied along bc_dir.
        struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, s.basis.num_basis, s.local_ext.volume);
        gkyl_array_clear(f, 0.0);
        gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
          .grid = &s.grid,
          .basis = &s.basis,
          .num_ret_vals = 1,
          .eval = init_donor_ky0,
          .ctx = &amp_ky[d],
        });
        gkyl_proj_on_basis_advance(proj, 0.0, &s.local, f);
        gkyl_proj_on_basis_release(proj);
        gkyl_array_copy_range_to_range(f, f, &s.ghost_r, &s.skin_r);

        struct gkyl_bc_twistshift *up = gkyl_bc_twistshift_inew(&(struct gkyl_bc_twistshift_inp){
          .bc_dir = ts_bc_dir,
          .shift_dir = 1,
          .shear_dir = 0,
          .edge = edge,
          .cdim = ts_cdim,
          .bcdir_ext_update_r = &s.update_r,
          .num_ghost = s.ghost,
          .basis = &s.basis,
          .grid = &s.grid,
          .shift_func = shift_func,
          .shift_func_ctx = &tctx,
          .use_gpu = false,
          .upsample_factor = upsample[c],
          .filter_half_width = half_width[c],
          .filter_cutoff_wavelength = cutoff[c],
        });
        gkyl_bc_twistshift_advance(up, f, f);
        gkyl_bc_twistshift_release(up);

        // Compare the y-sums of the y-independent coefficients, ghost vs skin, per x cell.
        int iz_gh = s.ghost_r.lower[ts_bc_dir], iz_sk = s.skin_r.lower[ts_bc_dir];
        double maxd = 0.0, maxref = 0.0;
        for (int ix = s.local.lower[0]; ix <= s.local.upper[0]; ix++) {
          for (int k = 0; k < 4; k++) {
            double sum_gh = 0.0, sum_sk = 0.0;
            for (int iy = s.local.lower[1]; iy <= s.local.upper[1]; iy++) {
              int idx_gh[] = {ix, iy, iz_gh}, idx_sk[] = {ix, iy, iz_sk};
              sum_gh += ((const double *)gkyl_array_cfetch(f, gkyl_range_idx(&s.local_ext, idx_gh))
              )[ky0_coeffs[k]];
              sum_sk += ((const double *)gkyl_array_cfetch(f, gkyl_range_idx(&s.local_ext, idx_sk))
              )[ky0_coeffs[k]];
            }
            maxd = GKYL_MAX2(maxd, fabs(sum_gh - sum_sk));
            maxref = GKYL_MAX2(maxref, fabs(sum_sk));
          }
        }

        TEST_CHECK(maxd < 1.0e-12 * maxref);
        TEST_MSG(
          "edge %d, ky!=0 amplitude %.1f, %s: max |<f>_y(ghost) - <f>_y(skin)| / max|<f>_y| = %.3e",
          e, amp_ky[d], label[c], maxd / maxref
        );

        gkyl_array_release(f);
      }
    }
  }
}

void
test_bc_twistshift_zero_shift_is_identity_ho(void)
{
  // Check that a zero shift leaves the ghost plane alone, for any upsampling.
  struct ts_ctx tctx = {.offset = 0.0, .shear = 0.0};
  double dx = (ts_upper[0] - ts_lower[0]) / ts_cells[0];

  for (int e = 0; e < 2; e++) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init_grid(&s, edge, ts_lower, ts_upper, ts_cells);

    struct gkyl_array *f_pre = ts_donor_new(&s);

    int upsample[] = {1, 2, 4};
    for (int c = 0; c < 3; c++) {
      int half_width = upsample[c] > 1 ? 4 : -1;
      double cutoff = upsample[c] > 1 ? 2.0 * dx / upsample[c] : 0.0;
      struct gkyl_array *f = ts_run(&s, edge, upsample[c], half_width, cutoff, NULL, &tctx);

      // The donor plane is the ghost plane periodicity just filled, so the BC
      // writing over it must reproduce what was already there.
      double maxd = ts_max_diff(f_pre, f, &s.ghost_r, s.basis.num_basis);
      TEST_CHECK(maxd < 1.0e-13);
      TEST_MSG("edge %d, upsample %d: zero shift changed the plane by %.3e", e, upsample[c], maxd);

      gkyl_array_release(f);
    }
    gkyl_array_release(f_pre);
  }
}

void
test_bc_twistshift_zero_shear_is_cell_translation_ho(void)
{
  // Check that a uniform shift is equivalent to a cell translation, for any upsampling.
  int num_cells = 3;
  double dy = (ts_upper[1] - ts_lower[1]) / ts_cells[1];
  struct ts_ctx tctx = {.offset = num_cells * dy, .shear = 0.0};

  for (int e = 0; e < 2; e++) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init_grid(&s, edge, ts_lower, ts_upper, ts_cells);

    struct gkyl_array *f_pre = ts_donor_new(&s);
    struct gkyl_array *f = ts_run(&s, edge, 1, -1, 0.0, NULL, &tctx);

    double maxd = 0.0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &s.ghost_r);
    while (gkyl_range_iter_next(&iter)) {
      int idx_do[GKYL_MAX_DIM];
      gkyl_copy_int_arr(ts_ndim, iter.idx, idx_do);
      int ny = ts_cells[1];
      idx_do[1] = ((iter.idx[1] - 1 - num_cells) % ny + ny) % ny + 1;

      const double *tar = gkyl_array_cfetch(f, gkyl_range_idx(&s.ghost_r, iter.idx));
      const double *don = gkyl_array_cfetch(f_pre, gkyl_range_idx(&s.ghost_r, idx_do));
      for (int k = 0; k < s.basis.num_basis; k++) {
        maxd = GKYL_MAX2(maxd, fabs(tar[k] - don[k]));
      }
    }
    TEST_CHECK(maxd < 1.0e-13);
    TEST_MSG(
      "edge %d: shift of %d cells is not an exact translation, off by %.3e", e, num_cells, maxd
    );

    gkyl_array_release(f);
    gkyl_array_release(f_pre);
  }
}

// Return the maximum jump between matching face quadrature points, including periodic y.
static double
ts_field_face_jump(
  const struct ts_setup *s, const struct gkyl_range *range, const struct gkyl_array *field, int dir
)
{
  double max_jump = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    int neighbor[3];
    gkyl_copy_int_arr(3, iter.idx, neighbor);
    if (iter.idx[dir] == range->upper[dir]) {
      if (dir != 1) {
        continue;
      }
      neighbor[dir] = range->lower[dir];
    } else {
      neighbor[dir]++;
    }
    const double *left = gkyl_array_cfetch(field, gkyl_range_idx(range, iter.idx));
    const double *right = gkyl_array_cfetch(field, gkyl_range_idx(range, neighbor));
    for (int qa = -1; qa <= 1; qa += 2) {
      for (int qb = -1; qb <= 1; qb += 2) {
        double point[3];
        point[dir] = 1;
        point[(dir + 1) % 3] = qa / sqrt(3.0);
        point[(dir + 2) % 3] = qb / sqrt(3.0);
        double value = s->basis.eval_expand(point, left);
        point[dir] = -1;
        max_jump = fmax(max_jump, fabs(value - s->basis.eval_expand(point, right)));
      }
    }
  }
  return max_jump;
}

// Check deflation onto the upper face and constant-in-z inflation, leaving other cells untouched.
// We use this combination in applying an FEM projection in x-y after applying TSBCs.
static void
ts_field_deflate_inflate(bool use_gpu)
{
  // Create grid and ranges.
  struct ts_setup setup;
  ts_setup_init_grid(&setup, GKYL_LOWER_EDGE, ts_lower, ts_upper, ts_cells);

  // Create 2D grid, range, basis.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, ts_lower, ts_upper, ts_cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1);
  struct gkyl_range ext, local;
  gkyl_create_grid_ranges(&grid, (int[]){1, 1}, &ext, &local);

  struct gkyl_array *host = mkarr(false, 8, setup.local_ext.volume);
  struct gkyl_array *result = mkarr(false, 8, setup.local_ext.volume);
  struct gkyl_array *input = mkarr(use_gpu, 8, host->size);
  struct gkyl_array *output = mkarr(use_gpu, 8, host->size);
  struct gkyl_array *boundary_field = mkarr(use_gpu, 4, ext.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &setup.local_ext);
  while (gkyl_range_iter_next(&iter)) {
    double *coeff = gkyl_array_fetch(host, gkyl_range_idx(&setup.local_ext, iter.idx));
    for (int mode = 0; mode < 8; mode++) {
      coeff[mode] = 0.17 * (mode + 1) + 0.03 * iter.idx[0];
    }
  }
  gkyl_array_copy(input, host);
  gkyl_array_clear(output, 17.0);

  struct gkyl_deflate_zsurf *deflate = gkyl_deflate_zsurf_new(&setup.basis, &basis, 1, use_gpu);
  struct gkyl_translate_dim *inflate =
    gkyl_translate_dim_new(2, basis, 3, setup.basis, 2, GKYL_UPPER_EDGE, use_gpu);
  gkyl_deflate_zsurf_advance(
    deflate, setup.ghost_r.lower[2], &setup.local_ext, &local, input, boundary_field, 1
  );

  gkyl_translate_dim_advance(inflate, &local, &setup.ghost_r, boundary_field, 1, output);
  gkyl_array_copy(result, output);

  double max_error = 0.0;
  gkyl_range_iter_init(&iter, &setup.local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&setup.local_ext, iter.idx);
    const double *got = gkyl_array_cfetch(result, loc);
    if (!gkyl_range_contains_idx(&setup.ghost_r, iter.idx)) {
      for (int mode = 0; mode < 8; mode++) {
        TEST_CHECK(got[mode] == 17.0);
      }
      continue;
    }
    const double *original = gkyl_array_cfetch(host, loc);
    for (int qx = -1; qx <= 1; qx += 2) {
      for (int qy = -1; qy <= 1; qy += 2) {
        double point[] = {qx / sqrt(3.0), qy / sqrt(3.0), 1};
        double expected = setup.basis.eval_expand(point, original);
        for (int edge = -1; edge <= 1; edge += 2) {
          point[2] = edge;
          max_error = fmax(max_error, fabs(expected - setup.basis.eval_expand(point, got)));
        }
      }
    }
  }
  TEST_CHECK(max_error < 1e-11);

  gkyl_deflate_zsurf_release(deflate);
  gkyl_translate_dim_release(inflate);
  gkyl_array_release(input);
  gkyl_array_release(output);
  gkyl_array_release(boundary_field);
  gkyl_array_release(host);
  gkyl_array_release(result);
}

static void
ts_field_init(double t, const double *xn, double *out, void *ctx)
{
  bool constant = *(bool *)ctx;
  out[0] = constant ? 1.0 : sin(M_PI * xn[0]) * cos(2 * M_PI * xn[1]);
}

// Check continuity and untouched cells through the shift/projection sequence, and return
// the L2 error of the FEM-projected boundary field relative to the analytic shifted field.
static double
ts_field_projection_case(int cells, bool use_gpu, bool zero_shift, bool constant, bool filtered)
{
  struct ts_setup setup;
  ts_setup_init_grid(
    &setup, GKYL_LOWER_EDGE, (double[]){0, 0, -1}, (double[]){1, 1, 1}, (int[]){cells, cells, 16}
  );
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, setup.grid.lower, setup.grid.upper, setup.grid.cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1);
  struct gkyl_range ext, local;
  gkyl_create_grid_ranges(&grid, (int[]){1, 1}, &ext, &local);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 8, setup.local_ext.volume);
  struct gkyl_array *original = gkyl_array_new(GKYL_DOUBLE, 8, host->size);
  struct gkyl_array *input = mkarr(use_gpu, 8, host->size);
  struct gkyl_array *output = mkarr(use_gpu, 8, host->size);
  struct gkyl_eval_on_nodes *eval =
    gkyl_eval_on_nodes_new(&setup.grid, &setup.basis, 1, ts_field_init, &constant);
  gkyl_eval_on_nodes_advance(eval, 0.0, &setup.local_ext, original);
  gkyl_array_copy(input, original);
  struct ts_ctx shift = {.offset = zero_shift ? 0 : 0.13, .shear = zero_shift ? 0 : 0.7};
  struct gkyl_bc_twistshift *ts = gkyl_bc_twistshift_inew(&(struct gkyl_bc_twistshift_inp){
    .bc_dir = 2,
    .shift_dir = 1,
    .shear_dir = 0,
    .edge = GKYL_LOWER_EDGE,
    .cdim = 3,
    .bcdir_ext_update_r = &setup.update_r,
    .num_ghost = setup.ghost,
    .basis = &setup.basis,
    .grid = &setup.grid,
    .shift_func = shift_func,
    .shift_func_ctx = &shift,
    .upsample_factor = filtered ? 6 : -1,
    .filter_half_width = filtered ? 4 : -1,
    .filter_cutoff_wavelength = 2 * grid.dx[0],
    .use_gpu = use_gpu,
  });
  gkyl_bc_twistshift_advance(ts, input, input);
  struct gkyl_array *boundary_field = mkarr(use_gpu, 4, ext.volume);
  struct gkyl_array *projected = mkarr(use_gpu, 4, ext.volume);
  struct gkyl_fem_proj_inp inp = {
    .grid = &grid,
    .solve_range = &local,
    .basis = &basis,
    .bcs =
      {
        .lo_type = {GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_PERIODIC},
        .up_type = {GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_PERIODIC},
        .lo_value = {constant ? 1.0 : 0.0},
        .up_value = {constant ? 1.0 : 0.0},
      },
    .use_gpu = use_gpu,
  };
  struct gkyl_fem_proj *proj = gkyl_fem_proj_inew(&inp);
  struct gkyl_deflate_zsurf *deflate = gkyl_deflate_zsurf_new(&setup.basis, &basis, 1, use_gpu);
  struct gkyl_translate_dim *inflate =
    gkyl_translate_dim_new(2, basis, 3, setup.basis, 2, GKYL_UPPER_EDGE, use_gpu);
  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(
    &setup.local, &setup.grid, &setup.basis, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, NULL, NULL, NULL,
    use_gpu
  );
  // Establish that the original operator chain has the defect on this field.
  gkyl_fem_parproj_set_rhs(parproj, input, input);
  gkyl_fem_parproj_solve(parproj, output);
  gkyl_array_copy(host, output);
  double before = ts_field_face_jump(&setup, &setup.local, host, 1);
  if (!zero_shift && !constant) {
    TEST_CHECK(before > 1e-7);
  }
  gkyl_deflate_zsurf_advance(deflate, 0, &setup.local_ext, &local, input, boundary_field, 1);
  gkyl_fem_proj_set_rhs(proj, boundary_field, NULL);
  gkyl_fem_proj_solve(proj, projected);
  gkyl_translate_dim_advance(inflate, &local, &setup.ghost_r, projected, 1, input);
  gkyl_array_copy(host, input);
  double error_sq = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &setup.local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&setup.local_ext, iter.idx);
    const double *got = gkyl_array_cfetch(host, loc), *ref = gkyl_array_cfetch(original, loc);
    if (!gkyl_range_contains_idx(&setup.ghost_r, iter.idx)) {
      for (int mode = 0; mode < 8; mode++) {
        TEST_CHECK(got[mode] == ref[mode]);
      }
      continue;
    }
    if (zero_shift || constant) {
      for (int mode = 0; mode < 8; mode++) {
        TEST_CHECK(fabs(got[mode] - ref[mode]) < 1e-11);
      }
    }
    double center[3];
    gkyl_rect_grid_cell_center(&setup.grid, iter.idx, center);
    // Three-point quadrature measures convergence to the smooth shifted function.
    const double quad[] = {-0.7745966692414834, 0, 0.7745966692414834};
    const double weights[] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    for (int qx = 0; qx < 3; qx++) {
      for (int qy = 0; qy < 3; qy++) {
        double point[] = {quad[qx], quad[qy], 1};
        double x = center[0] + 0.5 * grid.dx[0] * point[0];
        double y = center[1] + 0.5 * grid.dx[1] * point[1];
        double expected =
          constant ? 1.0 : sin(M_PI * x) * cos(2 * M_PI * (y - shift.offset - shift.shear * x));
        double err = setup.basis.eval_expand(point, got) - expected;
        error_sq += weights[qx] * weights[qy] * grid.cellVolume * err * err / 4;
      }
    }
  }
  for (int dir = 0; dir < 2; dir++) {
    TEST_CHECK(ts_field_face_jump(&setup, &setup.ghost_r, host, dir) < 1e-11);
  }
  gkyl_fem_parproj_set_rhs(parproj, input, input);
  gkyl_fem_parproj_solve(parproj, output);
  gkyl_array_copy(host, output);
  for (int dir = 0; dir < 3; dir++) {
    double jump = ts_field_face_jump(&setup, &setup.local, host, dir);
    TEST_CHECK(jump < 1e-11);
    TEST_MSG("nx=%d gpu=%d dir=%d jump=%.3e", cells, use_gpu, dir, jump);
  }
  gkyl_fem_parproj_release(parproj);
  gkyl_fem_proj_release(proj);
  gkyl_deflate_zsurf_release(deflate);
  gkyl_translate_dim_release(inflate);
  gkyl_bc_twistshift_release(ts);
  gkyl_eval_on_nodes_release(eval);
  gkyl_array_release(host);
  gkyl_array_release(original);
  gkyl_array_release(input);
  gkyl_array_release(output);
  gkyl_array_release(boundary_field);
  gkyl_array_release(projected);
  return sqrt(error_sq);
}

// Test deflation/inflation, zero shift, Preservation of a constant field, filtering, and p1
// convergence.
static void
test_ts_field_projection(bool use_gpu)
{
  ts_field_deflate_inflate(use_gpu); // Deflation/inflation
  ts_field_projection_case(8, use_gpu, true, false, false); // Zero shift.
  ts_field_projection_case(8, use_gpu, false, true, false); // Preserving a constant.
  ts_field_projection_case(8, use_gpu, false, false, true); // Filter.
  // Convergence.
  double prev = ts_field_projection_case(8, use_gpu, false, false, false);
  for (int cells = 16; cells <= 32; cells *= 2) {
    double err = ts_field_projection_case(cells, use_gpu, false, false, false);
    double order = log(prev / err) / log(2.0);
    TEST_CHECK(order > 1.8);
    TEST_MSG("nx=%d error=%.3e order=%.3f", cells, err, order);
    prev = err;
  }
}

static void
test_bc_twistshift_field_projection_ho(void)
{
  test_ts_field_projection(false);
}
#ifdef GKYL_HAVE_CUDA
static void
test_bc_twistshift_field_projection_dev(void)
{
  test_ts_field_projection(true);
}
#endif

TEST_LIST = {
  {"test_bc_twistshift_zero_shift_is_identity_ho", test_bc_twistshift_zero_shift_is_identity_ho},
  {"test_bc_twistshift_zero_shear_is_cell_translation_ho",
   test_bc_twistshift_zero_shear_is_cell_translation_ho},
  {"test_bc_twistshift_plain_matches_twistshift_dg_ho",
   test_bc_twistshift_plain_matches_twistshift_dg_ho},
  {"test_bc_twistshift_conserves_particles_ho", test_bc_twistshift_conserves_particles_ho},
  {"test_bc_twistshift_only_ghost_plane_written_ho", test_bc_twistshift_only_ghost_plane_written_ho
  },
  {"test_bc_twistshift_identity_filter_matches_plain_ho",
   test_bc_twistshift_identity_filter_matches_plain_ho},
  {"test_bc_twistshift_upsample_no_shear_matches_plain_ho",
   test_bc_twistshift_upsample_no_shear_matches_plain_ho},
  {"test_bc_twistshift_upsample_six_cbc_grid_ho", test_bc_twistshift_upsample_six_cbc_grid_ho},
  {"test_bc_twistshift_shift_dg_matches_shift_func_ho",
   test_bc_twistshift_shift_dg_matches_shift_func_ho},
  {"test_bc_twistshift_dealiasing_smooths_shear_direction_ho",
   test_bc_twistshift_dealiasing_smooths_shear_direction_ho},
  {"test_bc_twistshift_preserves_y_average_ho", test_bc_twistshift_preserves_y_average_ho},
  {"test_bc_twistshift_field_projection_ho", test_bc_twistshift_field_projection_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_bc_twistshift_field_projection_dev", test_bc_twistshift_field_projection_dev},
#endif
  {NULL, NULL}
};
