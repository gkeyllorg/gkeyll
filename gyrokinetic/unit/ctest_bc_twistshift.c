// Test the bc_twistshift orchestrator: the twist-shift BC on its own, and the
// supersample -> twist-shift -> low-pass filter -> restrict pipeline used to
// de-alias the shift along the shear direction. The sub-cell integrals the
// shift is built from are tested separately, in ctest_twistshift_dg.
#include <acutest.h>

#include <gkyl_array.h>
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

static const double ts_lower[] = {-2.0, -1.5, -3.0};
static const double ts_upper[] = { 2.0,  1.5,  3.0};
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
  fout[0] = tctx->offset + tctx->shear*xn[0];
}

static void
init_donor(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = exp(-y*y/(2.0*0.3*0.3))*(1.0 + 0.5*cos(M_PI*x));
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

static void
ts_setup_init(struct ts_setup *s, enum gkyl_edge_loc edge)
{
  gkyl_rect_grid_init(&s->grid, ts_ndim, ts_lower, ts_upper, ts_cells);
  gkyl_cart_modal_serendip(&s->basis, ts_ndim, 1);

  for (int d=0; d<ts_ndim; d++) s->ghost[d] = 1;
  gkyl_create_grid_ranges(&s->grid, s->ghost, &s->local_ext, &s->local);

  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d=0; d<ts_ndim; d++) {
    lo[d] = s->local.lower[d];
    up[d] = s->local.upper[d];
  }
  lo[ts_bc_dir] = s->local_ext.lower[ts_bc_dir];
  up[ts_bc_dir] = s->local_ext.upper[ts_bc_dir];
  gkyl_sub_range_init(&s->update_r, &s->local_ext, lo, up);

  // The BC reads and writes the same ghost plane; periodicity fills it from
  // the skin cell at the other end.
  int slo[GKYL_MAX_DIM], sup[GKYL_MAX_DIM];
  for (int d=0; d<ts_ndim; d++) {
    slo[d] = s->local.lower[d];
    sup[d] = s->local.upper[d];
  }
  if (edge == GKYL_LOWER_EDGE) {
    gkyl_range_shorten_from_above(&s->ghost_r, &s->update_r, ts_bc_dir, s->ghost[ts_bc_dir]);
    slo[ts_bc_dir] = sup[ts_bc_dir] = s->local.upper[ts_bc_dir];
  }
  else {
    gkyl_range_shorten_from_below(&s->ghost_r, &s->update_r, ts_bc_dir, s->ghost[ts_bc_dir]);
    slo[ts_bc_dir] = sup[ts_bc_dir] = s->local.lower[ts_bc_dir];
  }
  gkyl_sub_range_init(&s->skin_r, &s->local_ext, slo, sup);
}

// Project the donor and apply periodicity along bc_dir.
static struct gkyl_array*
ts_donor_new(const struct ts_setup *s)
{
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, s->basis.num_basis, s->local_ext.volume);
  gkyl_array_clear(f, 0.0);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
    .grid = &s->grid, .basis = &s->basis, .num_ret_vals = 1, .eval = init_donor, .ctx = NULL });
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local, f);
  gkyl_proj_on_basis_release(proj);

  gkyl_array_copy_range_to_range(f, f, &s->ghost_r, &s->skin_r);
  return f;
}

// The 1D DG shift on the shear cells of the update range, as the app builds it.
static struct gkyl_array*
ts_shift_dg_new(const struct ts_setup *s, struct ts_ctx *tctx)
{
  struct gkyl_rect_grid xgrid;
  gkyl_rect_grid_init(&xgrid, 1, &ts_lower[0], &ts_upper[0], &ts_cells[0]);
  struct gkyl_basis xbasis;
  gkyl_cart_modal_serendip(&xbasis, 1, s->basis.poly_order);

  struct gkyl_range xrng;
  gkyl_range_init(&xrng, 1, (int[]) {s->update_r.lower[0]}, (int[]) {s->update_r.upper[0]});

  struct gkyl_array *shift_dg = gkyl_array_new(GKYL_DOUBLE, xbasis.num_basis, xrng.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&xgrid, &xbasis, 1, shift_func, tctx);
  gkyl_eval_on_nodes_advance(ev, 0.0, &xrng, shift_dg);
  gkyl_eval_on_nodes_release(ev);
  return shift_dg;
}

// Apply the BC to a fresh donor field and return it. Pass shift_dg to exercise
// the discretized-shift input instead of the shift function.
static struct gkyl_array*
ts_run(const struct ts_setup *s, enum gkyl_edge_loc edge, int upsample, int half_width,
  double cutoff, struct gkyl_array *shift_dg, struct ts_ctx *tctx)
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
  if (shift_dg)
    inp.shift_dg = shift_dg;
  else {
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
ts_max_diff(const struct gkyl_array *fa, const struct gkyl_array *fb,
  const struct gkyl_range *rng, int num_basis)
{
  double maxd = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, rng);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(rng, iter.idx);
    const double *a = gkyl_array_cfetch(fa, linidx);
    const double *b = gkyl_array_cfetch(fb, linidx);
    for (int k=0; k<num_basis; k++)
      maxd = GKYL_MAX2(maxd, fabs(a[k]-b[k]));
  }
  return maxd;
}

// Largest difference over the whole array, ghost cells included.
static double
ts_max_diff_all(const struct ts_setup *s, const struct gkyl_array *fa,
  const struct gkyl_array *fb)
{
  double maxd = 0.0;
  for (long i=0; i<fa->size; i++) {
    const double *a = gkyl_array_cfetch(fa, i);
    const double *b = gkyl_array_cfetch(fb, i);
    for (int k=0; k<s->basis.num_basis; k++)
      maxd = GKYL_MAX2(maxd, fabs(a[k]-b[k]));
  }
  return maxd;
}

// Integral of the field over the plane the BC fills.
static double
ts_ghost_sum(const struct ts_setup *s, const struct gkyl_array *f)
{
  double tot = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->ghost_r);
  while (gkyl_range_iter_next(&iter))
    tot += ((const double *) gkyl_array_cfetch(f, gkyl_range_idx(&s->ghost_r, iter.idx)))[0];
  return tot;
}

void
test_plain_matches_twistshift_dg(void)
{
  // With no anti-aliasing the BC must be exactly the bare twist-shift.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 0.25 };
  for (int e=0; e<2; e++) {
    enum gkyl_edge_loc edge = e == 0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init(&s, edge);

    struct gkyl_array *f_bc = ts_run(&s, edge, 1, 0, 0.0, NULL, &tctx);

    struct gkyl_array *f_ref = ts_donor_new(&s);
    struct gkyl_twistshift_dg_inp tsinp = {
      .bc_dir = ts_bc_dir, .shift_dir = 1, .shear_dir = 0, .edge = edge,
      .cdim = ts_cdim, .bcdir_ext_update_r = &s.update_r, .num_ghost = s.ghost,
      .basis = &s.basis, .grid = &s.grid,
      .shift_func = shift_func, .shift_func_ctx = &tctx, .use_gpu = false,
    };
    struct gkyl_twistshift_dg *ts = gkyl_twistshift_dg_inew(&tsinp);
    gkyl_twistshift_dg_advance(ts, f_ref, f_ref);
    gkyl_twistshift_dg_release(ts);

    TEST_CHECK( ts_max_diff_all(&s, f_bc, f_ref) == 0.0 );
    TEST_MSG("edge %d: bc_twistshift differs from twistshift_dg by %.3e", e,
      ts_max_diff_all(&s, f_bc, f_ref));

    gkyl_array_release(f_bc);
    gkyl_array_release(f_ref);
  }
}

void
test_only_ghost_plane_written(void)
{
  // Nothing but the ghost plane the BC fills may change, in any configuration.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 0.25 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];
  int upsample[] = {1, 1, 2, 4};
  int half_width[] = {0, 4, 4, 4};
  double cutoff[] = {0.0, 4.0*dx, 2.0*dx, 2.0*dx};

  for (int e=0; e<2; e++) {
    enum gkyl_edge_loc edge = e == 0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init(&s, edge);

    struct gkyl_array *f_pre = ts_donor_new(&s);

    for (int c=0; c<4; c++) {
      struct gkyl_array *f = ts_run(&s, edge, upsample[c], half_width[c], cutoff[c], NULL, &tctx);

      // Restore the ghost plane so any remaining difference is a stray write.
      gkyl_array_copy_range(f, f_pre, &s.ghost_r);
      double maxd = ts_max_diff_all(&s, f, f_pre);
      TEST_CHECK( maxd == 0.0 );
      TEST_MSG("edge %d, config %d (upsample %d, M %d): wrote %.3e outside the ghost plane",
        e, c, upsample[c], half_width[c], maxd);

      gkyl_array_release(f);
    }
    gkyl_array_release(f_pre);
  }
}

void
test_identity_filter_matches_plain(void)
{
  // A cutoff at the Nyquist wavelength makes the kernel the identity, so
  // filtering must leave the plain twist-shift result alone.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 0.25 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];

  struct ts_setup s;
  ts_setup_init(&s, GKYL_LOWER_EDGE);

  struct gkyl_array *f_plain = ts_run(&s, GKYL_LOWER_EDGE, 1, 0, 0.0, NULL, &tctx);
  struct gkyl_array *f_filt = ts_run(&s, GKYL_LOWER_EDGE, 1, 4, 2.0*dx, NULL, &tctx);

  double maxd = ts_max_diff(f_plain, f_filt, &s.ghost_r, s.basis.num_basis);
  TEST_CHECK( maxd < 1.0e-13 );
  TEST_MSG("identity filter changed the result by %.3e", maxd);

  gkyl_array_release(f_plain);
  gkyl_array_release(f_filt);
}

void
test_upsample_no_shear_matches_plain(void)
{
  // With a uniform shift the twist-shift creates no structure along x, so the
  // supersampled field stays an exact refinement of the coarse one and the
  // refine/identity-filter/restrict round trip must be lossless.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 0.0 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];

  for (int upsample=2; upsample<=4; upsample*=2) {
    struct ts_setup s;
    ts_setup_init(&s, GKYL_LOWER_EDGE);

    struct gkyl_array *f_plain = ts_run(&s, GKYL_LOWER_EDGE, 1, 0, 0.0, NULL, &tctx);
    struct gkyl_array *f_up = ts_run(&s, GKYL_LOWER_EDGE, upsample, 4, 2.0*dx/upsample, NULL, &tctx);

    double maxd = ts_max_diff(f_plain, f_up, &s.ghost_r, s.basis.num_basis);
    TEST_CHECK( maxd < 1.0e-12 );
    TEST_MSG("upsample %d round trip changed the result by %.3e", upsample, maxd);

    gkyl_array_release(f_plain);
    gkyl_array_release(f_up);
  }
}

void
test_shift_dg_matches_shift_func(void)
{
  // A discretized shift given on the field's own grid must give the same answer
  // as the shift function, including when the shift has to be refined.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 0.25 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];

  for (int upsample=1; upsample<=2; upsample++) {
    struct ts_setup s;
    ts_setup_init(&s, GKYL_LOWER_EDGE);
    struct gkyl_array *shift_dg = ts_shift_dg_new(&s, &tctx);

    int half_width = upsample > 1? 4 : 0;
    double cutoff = upsample > 1? 2.0*dx : 0.0;

    struct gkyl_array *f_func = ts_run(&s, GKYL_LOWER_EDGE, upsample, half_width, cutoff, NULL, &tctx);
    struct gkyl_array *f_dg = ts_run(&s, GKYL_LOWER_EDGE, upsample, half_width, cutoff, shift_dg, &tctx);

    double maxd = ts_max_diff(f_func, f_dg, &s.ghost_r, s.basis.num_basis);
    TEST_CHECK( maxd < 1.0e-12 );
    TEST_MSG("upsample %d: shift_dg differs from shift_func by %.3e", upsample, maxd);

    gkyl_array_release(f_func);
    gkyl_array_release(f_dg);
    gkyl_array_release(shift_dg);
  }
}

void
test_dealiasing_smooths_shear_direction(void)
{
  // Check that the anti-aliasing filter actually smooths the field along the shear direction.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 2.9 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];

  struct ts_setup s;
  ts_setup_init(&s, GKYL_LOWER_EDGE);

  struct gkyl_array *f_plain = ts_run(&s, GKYL_LOWER_EDGE, 1, 0, 0.0, NULL, &tctx);
  struct gkyl_array *f_deal = ts_run(&s, GKYL_LOWER_EDGE, 4, 4, 2.0*dx, NULL, &tctx);

  // Total variation along x of the cell averages, summed over the ghost plane.
  double tv[2] = {0.0, 0.0};
  struct gkyl_array *fs[] = {f_plain, f_deal};
  for (int q=0; q<2; q++) {
    for (int j=s.ghost_r.lower[1]; j<=s.ghost_r.upper[1]; j++) {
      for (int i=s.ghost_r.lower[0]; i<s.ghost_r.upper[0]; i++) {
        int idx0[] = {i, j, s.ghost_r.lower[2]};
        int idx1[] = {i+1, j, s.ghost_r.lower[2]};
        const double *a = gkyl_array_cfetch(fs[q], gkyl_range_idx(&s.ghost_r, idx0));
        const double *b = gkyl_array_cfetch(fs[q], gkyl_range_idx(&s.ghost_r, idx1));
        tv[q] += fabs(b[0]-a[0]);
      }
    }
  }
  TEST_CHECK( tv[1] < tv[0] );
  TEST_MSG("total variation along x: plain %.6e, de-aliased %.6e", tv[0], tv[1]);

  gkyl_array_release(f_plain);
  gkyl_array_release(f_deal);
}

void
test_conserves_particles(void)
{
  // Check that the integral is preserved.
  struct ts_ctx tctx = { .offset = 0.75, .shear = 2.9 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];
  int upsample[] = {1, 1, 2, 4};
  int half_width[] = {0, 4, 4, 4};
  double cutoff[] = {0.0, 4.0*dx, 2.0*dx, 2.0*dx};

  for (int e=0; e<2; e++) {
    enum gkyl_edge_loc edge = e == 0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init(&s, edge);

    struct gkyl_array *f_pre = ts_donor_new(&s);
    double tot_pre = ts_ghost_sum(&s, f_pre);

    for (int c=0; c<4; c++) {
      struct gkyl_array *f = ts_run(&s, edge, upsample[c], half_width[c], cutoff[c], NULL, &tctx);
      double tot = ts_ghost_sum(&s, f);

      TEST_CHECK( fabs(tot-tot_pre) < 1.0e-12*fabs(tot_pre) );
      TEST_MSG("edge %d, config %d (upsample %d, M %d): %.13e -> %.13e, rel change %.3e",
        e, c, upsample[c], half_width[c], tot_pre, tot, fabs(tot-tot_pre)/fabs(tot_pre));

      gkyl_array_release(f);
    }
    gkyl_array_release(f_pre);
  }
}

void
test_zero_shift_is_identity(void)
{
  // Check that a zero shift leaves the ghost plane alone, for any upsampling.
  struct ts_ctx tctx = { .offset = 0.0, .shear = 0.0 };
  double dx = (ts_upper[0]-ts_lower[0])/ts_cells[0];

  for (int e=0; e<2; e++) {
    enum gkyl_edge_loc edge = e == 0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init(&s, edge);

    struct gkyl_array *f_pre = ts_donor_new(&s);

    int upsample[] = {1, 2, 4};
    for (int c=0; c<3; c++) {
      int half_width = upsample[c] > 1? 4 : 0;
      double cutoff = upsample[c] > 1? 2.0*dx/upsample[c] : 0.0;
      struct gkyl_array *f = ts_run(&s, edge, upsample[c], half_width, cutoff, NULL, &tctx);

      // The donor plane is the ghost plane periodicity just filled, so the BC
      // writing over it must reproduce what was already there.
      double maxd = ts_max_diff(f_pre, f, &s.ghost_r, s.basis.num_basis);
      TEST_CHECK( maxd < 1.0e-13 );
      TEST_MSG("edge %d, upsample %d: zero shift changed the plane by %.3e", e, upsample[c], maxd);

      gkyl_array_release(f);
    }
    gkyl_array_release(f_pre);
  }
}

void
test_zero_shear_is_cell_translation(void)
{
  // Check that a uniform shift is equivalent to a cell translation, for any upsampling.
  int num_cells = 3;
  double dy = (ts_upper[1]-ts_lower[1])/ts_cells[1];
  struct ts_ctx tctx = { .offset = num_cells*dy, .shear = 0.0 };

  for (int e=0; e<2; e++) {
    enum gkyl_edge_loc edge = e == 0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct ts_setup s;
    ts_setup_init(&s, edge);

    struct gkyl_array *f_pre = ts_donor_new(&s);
    struct gkyl_array *f = ts_run(&s, edge, 1, 0, 0.0, NULL, &tctx);

    double maxd = 0.0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &s.ghost_r);
    while (gkyl_range_iter_next(&iter)) {
      int idx_do[GKYL_MAX_DIM];
      gkyl_copy_int_arr(ts_ndim, iter.idx, idx_do);
      int ny = ts_cells[1];
      idx_do[1] = ((iter.idx[1]-1-num_cells) % ny + ny) % ny + 1;

      const double *tar = gkyl_array_cfetch(f, gkyl_range_idx(&s.ghost_r, iter.idx));
      const double *don = gkyl_array_cfetch(f_pre, gkyl_range_idx(&s.ghost_r, idx_do));
      for (int k=0; k<s.basis.num_basis; k++)
        maxd = GKYL_MAX2(maxd, fabs(tar[k]-don[k]));
    }
    TEST_CHECK( maxd < 1.0e-13 );
    TEST_MSG("edge %d: shift of %d cells is not an exact translation, off by %.3e",
      e, num_cells, maxd);

    gkyl_array_release(f);
    gkyl_array_release(f_pre);
  }
}

TEST_LIST = {
  { "test_zero_shift_is_identity", test_zero_shift_is_identity },
  { "test_zero_shear_is_cell_translation", test_zero_shear_is_cell_translation },
  { "test_plain_matches_twistshift_dg", test_plain_matches_twistshift_dg },
  { "test_conserves_particles", test_conserves_particles },
  { "test_only_ghost_plane_written", test_only_ghost_plane_written },
  { "test_identity_filter_matches_plain", test_identity_filter_matches_plain },
  { "test_upsample_no_shear_matches_plain", test_upsample_no_shear_matches_plain },
  { "test_shift_dg_matches_shift_func", test_shift_dg_matches_shift_func },
  { "test_dealiasing_smooths_shear_direction", test_dealiasing_smooths_shear_direction },
  { NULL, NULL },
};
