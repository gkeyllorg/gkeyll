// Tests of the PKPM configuration-space face flux (dg_pkpm_remap_flux_surf): the HLL flux of the
// two sides' distributions mapped exactly to the face frame, blended with the label-paired flux.
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_pkpm_remap_flux_kernels.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

struct test_setup {
  int cdim, pdim, nv, nc, np, nk, na;
  struct gkyl_vlasov_position_map *pos_map;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext;
  struct gkyl_array *bvar, *fluid, *u, *xi, *theta, *f, *flux, *fluid_flux, *rate_f, *rate_fluid,
    *cfl_f, *cfl_fluid, *blend;
  struct gkyl_dg_pkpm_remap_flux_surf *up;
};

// Setup with a position map in every direction (pmap = NULL: the identity map).
static void
setup_cells_map(struct test_setup *s, int cdim, int nv, double wmax, const int *nx, evalf_t pmap)
{
  int pdim = cdim + 1;
  s->cdim = cdim;
  s->pdim = pdim;
  s->nv = nv;
  double lower[4], upper[4];
  int cells[4], ccells[3], cghost[3], pghost[4];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = -1.0;
    upper[d] = 1.0;
    cells[d] = nx[d];
    ccells[d] = nx[d];
    cghost[d] = 1;
    pghost[d] = 1;
  }
  lower[cdim] = -wmax;
  upper[cdim] = wmax;
  cells[cdim] = nv;
  pghost[cdim] = 0;
  gkyl_rect_grid_init(&s->cgrid, cdim, lower, upper, ccells);
  gkyl_rect_grid_init(&s->pgrid, pdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  gkyl_create_grid_ranges(&s->pgrid, pghost, &s->plocal_ext, &s->plocal);
  struct gkyl_vlasov_position_map_inp pinp[GKYL_MAX_CDIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    pinp[d].eval_pmap = pmap;
  }
  s->pos_map =
    gkyl_vlasov_position_map_new(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis, pinp, false);

  s->nc = s->cbasis.num_basis;
  s->np = s->pbasis.num_basis;
  long cvol = s->clocal_ext.volume, pvol = s->plocal_ext.volume;
  s->bvar = mkarr(3 * s->nc, cvol);
  s->fluid = mkarr(5 * s->nc, cvol);
  s->u = mkarr(3 * s->nc, cvol);
  s->xi = mkarr(3 * s->nc, cvol);
  s->theta = mkarr(s->nc, cvol);
  s->f = mkarr(8 * s->np, pvol);
  s->rate_f = mkarr(cdim, pvol);
  s->rate_fluid = mkarr(cdim, cvol);
  s->cfl_f = mkarr(1, pvol);
  s->cfl_fluid = mkarr(1, cvol);
  s->blend = mkarr(cdim, cvol);
  gkyl_array_clear(s->f, 0.0);
  gkyl_array_clear(s->rate_f, 0.0);
  gkyl_array_clear(s->rate_fluid, 0.0);
  gkyl_array_clear(s->cfl_f, 0.0);
  gkyl_array_clear(s->cfl_fluid, 0.0);

  struct gkyl_dg_pkpm_remap_flux_surf_inp inp = {
    .phase_grid = &s->pgrid,
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .pos_map = s->pos_map,
    .num_face_nodes = 5,
    .use_gpu = false,
  };
  s->up = gkyl_dg_pkpm_remap_flux_surf_inew(&inp);
  s->nk = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(s->up);
  s->na = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(s->up);
  TEST_CHECK(s->nk == s->np / 2);
  s->flux = mkarr(2 * 8 * cdim * s->nk, pvol);
  s->fluid_flux = mkarr(5 * cdim * s->na, cvol);
  gkyl_array_clear(s->flux, 0.0);
  gkyl_array_clear(s->fluid_flux, 0.0);
}

static void
setup_cells(struct test_setup *s, int cdim, int nv, double wmax, const int *nx)
{
  setup_cells_map(s, cdim, nv, wmax, nx, NULL);
}

static void
setup(struct test_setup *s, int cdim, int nv, double wmax)
{
  const int nx[3] = {4, 3, 3};
  setup_cells(s, cdim, nv, wmax, nx);
}

static void
teardown(struct test_setup *s)
{
  gkyl_dg_pkpm_remap_flux_surf_release(s->up);
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_array_release(s->bvar);
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->u);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->f);
  gkyl_array_release(s->flux);
  gkyl_array_release(s->fluid_flux);
  gkyl_array_release(s->rate_f);
  gkyl_array_release(s->rate_fluid);
  gkyl_array_release(s->cfl_f);
  gkyl_array_release(s->cfl_fluid);
  gkyl_array_release(s->blend);
}

// Set the frame of one configuration cell to constants.
static void
set_frame(
  struct test_setup *s, long loc, double rho, const double *u, double xi, const double *b,
  double theta
)
{
  double c0 = sqrt(pow(2.0, s->cdim));
  double *fl = gkyl_array_fetch(s->fluid, loc), *ul = gkyl_array_fetch(s->u, loc);
  double *xl = gkyl_array_fetch(s->xi, loc), *bl = gkyl_array_fetch(s->bvar, loc);
  double *tl = gkyl_array_fetch(s->theta, loc);
  for (int k = 0; k < 5 * s->nc; ++k) {
    fl[k] = 0.0;
  }
  for (int k = 0; k < 3 * s->nc; ++k) {
    ul[k] = xl[k] = bl[k] = 0.0;
  }
  for (int k = 0; k < s->nc; ++k) {
    tl[k] = 0.0;
  }
  fl[0] = rho * c0;
  xl[0] = xi * c0;
  tl[0] = theta * c0;
  for (int k = 0; k < 3; ++k) {
    ul[k * s->nc] = u[k] * c0;
    bl[k * s->nc] = b[k] * c0;
  }
}

static void
set_uniform_frame(
  struct test_setup *s, double rho, const double *u, double xi, const double *b, double theta
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    set_frame(s, gkyl_range_idx(&s->clocal_ext, iter.idx), rho, u, xi, b, theta);
  }
}

static void
randomize_f(struct test_setup *s)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->plocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(s->f, gkyl_range_idx(&s->plocal_ext, iter.idx));
    for (int k = 0; k < 8 * s->np; ++k) {
      f[k] = 2.0 * rand() / RAND_MAX - 1.0;
    }
  }
}

static void
advance(struct test_setup *s)
{
  gkyl_dg_pkpm_remap_flux_surf_advance(
    s->up, &s->clocal, &s->clocal_ext, &s->plocal, &s->plocal_ext, s->bvar, s->fluid, s->u, s->xi,
    s->theta, s->f, s->flux, s->fluid_flux, s->rate_f, s->rate_fluid, s->cfl_f, s->cfl_fluid,
    s->blend
  );
}

// The orthonormal Legendre coefficients of a stored 1x flux of one block.
static void
stored_to_legendre_1x(const double *stored, double *coef)
{
  double gamma[3];
  for (int m = 0; m < 3; ++m) {
    double proj[24] = {0.0}, flux[24] = {0.0}, wpsi[1] = {1.0};
    proj[m * 8] = 1.0;
    pkpm_remap_flux_store_x_1x1v_ser_p1(wpsi, proj, flux);
    gamma[m] = flux[m];
  }
  for (int m = 0; m < 3; ++m) {
    coef[m] = stored[m] / gamma[m];
  }
}

// 1. Identity frames with b along the normal and |u_n| > w_max xi, so that no kink falls in the
// grid and the HLL flux is the flux function of the upwind side: both sides of a face receive the
// same flux, with no blending. In 1x that flux is checked against the projection of
//   (u_n + xi w) [F, G, M0_perp, N_perp] + xi [M0_n, N_n, 0, 0]
// of the upwind cell's trace on the Legendre modes of each velocity cell.
static void
test_identity_vs_label(int cdim, int dir)
{
  struct test_setup s;
  setup(&s, cdim, 8, 6.0);
  double u[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0}, xi = 1.0;
  u[dir] = 8.0;
  b[dir] = 1.0;
  set_uniform_frame(&s, 1.3, u, xi, b, 1.0);
  randomize_f(&s);
  advance(&s);

  // 4-point Gauss-Legendre rule (the flux is cubic in w) and the orthonormal Legendre modes
  static const double gn[4] = {
    -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526
  };
  static const double gw[4] = {
    0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538
  };

  double max_err = 0.0, max_val = 0.0, max_blend = 0.0, max_ref_err = 0.0;
  int nflux = 8 * s.nk;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.plocal, iter.idx);
    const double *fr = gkyl_array_cfetch(s.flux, loc);
    for (int m = 0; m < nflux; ++m) {
      double lower = fr[2 * dir * nflux + m], upper = fr[(2 * dir + 1) * nflux + m];
      max_err = fmax(max_err, fabs(lower - upper));
      max_val = fmax(max_val, fabs(lower));
    }
    long cloc = gkyl_range_idx(&s.clocal, iter.idx);
    max_blend = fmax(max_blend, ((const double *)gkyl_array_cfetch(s.blend, cloc))[dir]);

    if (cdim == 1) {
      // upwind cell of the lower face of this cell, and its trace on that face
      int up_idx[2] = {iter.idx[0] - 1, iter.idx[1]};
      const double *fu = gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, up_idx));
      double xc[2], dw = s.pgrid.dx[1], ref[8][3] = {{0.0}};
      gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
      for (int q = 0; q < 4; ++q) {
        double z[2] = {1.0, gn[q]}, X[8], J[8];
        for (int c = 0; c < 8; ++c) {
          X[c] = s.pbasis.eval_expand(z, &fu[c * s.np]);
        }
        double V = u[0] + xi * (xc[1] + 0.5 * dw * gn[q]);
        for (int c = 0; c < 8; ++c) {
          J[c] = V * X[c];
        }
        J[0] += xi * X[2];
        J[1] += xi * X[5];
        double leg[3] = {sqrt(0.5), sqrt(1.5) * gn[q], sqrt(0.625) * (3.0 * gn[q] * gn[q] - 1.0)};
        for (int c = 0; c < 8; ++c) {
          for (int m = 0; m < 3; ++m) {
            ref[c][m] += gw[q] * J[c] * leg[m];
          }
        }
      }
      for (int c = 0; c < 8; ++c) {
        double coef[3];
        stored_to_legendre_1x(&fr[c * s.nk], coef);
        for (int m = 0; m < 3; ++m) {
          max_ref_err = fmax(max_ref_err, fabs(coef[m] - ref[c][m]));
        }
      }
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-13 * max_val);
  TEST_MSG("cdim %d dir %d: max difference of the two sides %g of %g", cdim, dir, max_err, max_val);
  TEST_CHECK(max_blend == 0.0);
  TEST_MSG("blend %g", max_blend);
  if (cdim == 1) {
    TEST_CHECK(max_ref_err < 1e-13 * max_val);
    TEST_MSG("max difference from the upwind flux function %g of %g", max_ref_err, max_val);
  }

  teardown(&s);
}

void
test_identity_1x(void)
{
  test_identity_vs_label(1, 0);
}
void
test_identity_2x(void)
{
  test_identity_vs_label(2, 0);
  test_identity_vs_label(2, 1);
}
void
test_identity_3x(void)
{
  test_identity_vs_label(3, 0);
  test_identity_vs_label(3, 1);
  test_identity_vs_label(3, 2);
}

// ---- Helpers: the blocks of one phase cell from a function of w alone.
typedef void (*wfunc_t)(double w, void *ctx, double *blocks);

static void
project_cell(const struct test_setup *s, const int *idx, wfunc_t fn, void *ctx, double *f)
{
  static const double gn[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  static const double gw[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
  int cdim = s->cdim, np = s->np;
  double xc[4], basis[64];
  gkyl_rect_grid_cell_center(&s->pgrid, idx, xc);
  double dw = s->pgrid.dx[cdim];
  for (int k = 0; k < 8 * np; ++k) {
    f[k] = 0.0;
  }
  int node_lower[3] = {0, 0, 0}, node_upper[3] = {1, 1, 1};
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, cdim, node_lower, node_upper);
  struct gkyl_range_iter node_iter;
  gkyl_range_iter_init(&node_iter, &node_range);
  while (gkyl_range_iter_next(&node_iter)) {
    double z[4];
    for (int d = 0; d < cdim; ++d) {
      z[d] = gkyl_gauss_ordinates_2[node_iter.idx[d]];
    }
    for (int j = 0; j < 3; ++j) {
      z[cdim] = gn[j];
      double blocks[8];
      fn(xc[cdim] + 0.5 * dw * gn[j], ctx, blocks);
      s->pbasis.eval(z, basis);
      for (int c = 0; c < 8; ++c) {
        for (int k = 0; k < np; ++k) {
          f[c * np + k] += gw[j] * blocks[c] * basis[k];
        }
      }
    }
  }
}

struct phys_frame {
  double rho, un, xi;
};

// A compactly supported physical distribution in the frame (rho, un, xi): F = (xi/rho) g(un + xi w),
// and G with the scaling that keeps it the same physical distribution of the second moment.
static void
same_physical(double w, void *ctx, double *blocks)
{
  const struct phys_frame *fr = ctx;
  double v = fr->un + fr->xi * w, g = fabs(v) < 4.0 ? 16.0 - v * v : 0.0;
  for (int c = 0; c < 8; ++c) {
    blocks[c] = 0.0;
  }
  blocks[0] = fr->xi / fr->rho * g;
  blocks[1] = fr->xi / fr->rho * g / (fr->xi * fr->xi);
}

// 2. Every cell holds the same physical distribution in its own frame: the flux each cell receives
// through its two faces of direction dir agrees, so the face divergence vanishes.
static void
test_steady_state(int cdim, int dir, const struct phys_frame *frames)
{
  struct test_setup s;
  setup(&s, cdim, 24, 12.0);
  double b[3] = {0.0, 0.0, 0.0};
  b[dir] = 1.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    const struct phys_frame *fr = &frames[(iter.idx[dir] + 16) % 2];
    double u[3] = {0.0, 0.0, 0.0};
    u[dir] = fr->un;
    set_frame(&s, gkyl_range_idx(&s.clocal_ext, iter.idx), fr->rho, u, fr->xi, b, 1.0);
  }
  gkyl_range_iter_init(&iter, &s.plocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    struct phys_frame fr = frames[(iter.idx[dir] + 16) % 2];
    project_cell(
      &s, iter.idx, same_physical, &fr,
      gkyl_array_fetch(s.f, gkyl_range_idx(&s.plocal_ext, iter.idx))
    );
  }
  advance(&s);

  int nflux = 8 * s.nk;
  double max_err = 0.0, max_val = 0.0, max_blend = 0.0;
  int idx_up[GKYL_MAX_DIM];
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(s.pdim, iter.idx, idx_up);
    idx_up[dir] += 1;
    const double *lo = gkyl_array_cfetch(s.flux, gkyl_range_idx(&s.plocal_ext, iter.idx));
    const double *up = gkyl_array_cfetch(s.flux, gkyl_range_idx(&s.plocal_ext, idx_up));
    for (int m = 0; m < nflux; ++m) {
      double recv_lo = lo[(2 * dir + 1) * nflux + m], recv_up = up[2 * dir * nflux + m];
      max_err = fmax(max_err, fabs(recv_lo - recv_up));
      max_val = fmax(max_val, fabs(recv_lo));
    }
    max_blend = fmax(
      max_blend,
      ((const double *)gkyl_array_cfetch(s.blend, gkyl_range_idx(&s.clocal_ext, iter.idx)))[dir]
    );
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 2e-13 * max_val);
  TEST_MSG("cdim %d dir %d: %g of %g", cdim, dir, max_err, max_val);
  TEST_CHECK(max_blend == 0.0);
  TEST_MSG("blend %g", max_blend);
  teardown(&s);
}

static void
test_steady_state_all(const struct phys_frame *frames)
{
  test_steady_state(1, 0, frames);
  for (int d = 0; d < 2; ++d) {
    test_steady_state(2, d, frames);
  }
  for (int d = 0; d < 3; ++d) {
    test_steady_state(3, d, frames);
  }
}

void
test_steady_dilation(void)
{
  // xi = 1 and 2: the xi = 2 cell edges inside the xi = 1 grid coincide with xi = 1 edges in v
  struct phys_frame frames[2] = {{1.0, 0.0, 1.0}, {3.0, 0.0, 2.0}};
  test_steady_state_all(frames);
}

void
test_steady_shift(void)
{
  // a jump of u by two velocity cells: both sides' labels shift by one cell at the face
  struct phys_frame frames[2] = {{1.0, 0.0, 1.0}, {2.0, 2.0, 1.0}};
  test_steady_state_all(frames);
}

// 3. With no flux outside either grid, each side's kinetic flux has the moments of the fluid flux
// when mapped with that side's own frame.
void
test_moments(void)
{
  struct test_setup s;
  setup(&s, 1, 24, 12.0);
  const double rho[2] = {1.0, 3.0}, xi[2] = {1.0, 1.7};
  const double us[2][3] = {{0.2, 0.4, -0.3}, {1.3, -0.5, 0.6}};
  double b[3] = {1.0, 0.0, 0.0};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    int par = (iter.idx[0] + 16) % 2;
    set_frame(&s, gkyl_range_idx(&s.clocal_ext, iter.idx), rho[par], us[par], xi[par], b, 1.0);
  }
  srand(99);
  double xc[2];
  gkyl_range_iter_init(&iter, &s.plocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(s.f, gkyl_range_idx(&s.plocal_ext, iter.idx));
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    for (int k = 0; k < 8 * s.np; ++k) {
      f[k] = fabs(xc[1]) < 3.0 ? 2.0 * rand() / RAND_MAX - 1.0 : 0.0;
    }
    // the first-harmonic blocks are perpendicular to b = e_x
    for (int k = 0; k < s.np; ++k) {
      f[2 * s.np + k] = 0.0;
      f[5 * s.np + k] = 0.0;
    }
  }
  advance(&s);

  static const double gn[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  static const double gw[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
  double leg[3][3];
  for (int j = 0; j < 3; ++j) {
    leg[j][0] = 0.7071067811865475;
    leg[j][1] = 1.224744871391589 * gn[j];
    leg[j][2] = 0.7905694150420949 * (3.0 * gn[j] * gn[j] - 1.0);
  }
  double dw = s.pgrid.dx[1], max_err = 0.0, max_val = 0.0;
  for (int ix = s.clocal.lower[0]; ix <= s.clocal.upper[0]; ++ix) {
    // the face between ix - 1 and ix: sides 0 and 1
    long cloc = gkyl_range_idx(&s.clocal_ext, &ix);
    int par[2] = {(ix - 1 + 16) % 2, (ix + 16) % 2};
    const double *fluid_flux = gkyl_array_cfetch(s.fluid_flux, cloc);
    TEST_CHECK(((const double *)gkyl_array_cfetch(s.blend, cloc))[0] == 0.0);
    for (int side = 0; side < 2; ++side) {
      double mom[7] = {0.0};
      for (int k = s.plocal.lower[1]; k <= s.plocal.upper[1]; ++k) {
        int pidx[2] = {ix, k};
        gkyl_rect_grid_cell_center(&s.pgrid, pidx, xc);
        const double *flux = gkyl_array_cfetch(s.flux, gkyl_range_idx(&s.plocal_ext, pidx));
        double cc[8][3];
        for (int c = 0; c < 8; ++c) {
          stored_to_legendre_1x(&flux[side * 8 * s.nk + c * s.nk], cc[c]);
        }
        for (int j = 0; j < 3; ++j) {
          double wq = xc[1] + 0.5 * dw * gn[j], wgt = 0.5 * dw * gw[j], val[8];
          for (int c = 0; c < 8; ++c) {
            val[c] = cc[c][0] * leg[j][0] + cc[c][1] * leg[j][1] + cc[c][2] * leg[j][2];
          }
          mom[0] += wgt * val[0];
          mom[1] += wgt * wq * val[0];
          mom[2] += wgt * wq * wq * val[0];
          mom[3] += wgt * val[1];
          for (int a = 0; a < 3; ++a) {
            mom[4 + a] += wgt * val[2 + a];
          }
        }
      }
      double side_flux[5];
      pkpm_harm_fluid_moment_map(mom, rho[par[side]], us[par[side]], xi[par[side]], b, side_flux);
      for (int m = 0; m < 5; ++m) {
        max_err = fmax(max_err, fabs(side_flux[m] - fluid_flux[m]));
        max_val = fmax(max_val, fabs(fluid_flux[m]));
      }
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_val);
  TEST_MSG("moments: %g of %g", max_err, max_val);
  teardown(&s);
}

// 4. A reflecting wall (mirror ghost: u_n reversed, f(x, w) -> f(-x, -w)) has no normal mass flux.
void
test_reflecting_wall(void)
{
  struct test_setup s;
  setup(&s, 1, 16, 6.0);
  double b[3] = {1.0, 0.0, 0.0}, u_skin[3] = {0.7, 0.2, -0.1}, u_ghost[3] = {-0.7, 0.2, -0.1};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    int is_ghost = iter.idx[0] < s.clocal.lower[0];
    set_frame(
      &s, gkyl_range_idx(&s.clocal_ext, iter.idx), 1.7, is_ghost ? u_ghost : u_skin, 1.2, b, 1.0
    );
  }
  srand(7);
  randomize_f(&s);
  gkyl_range_iter_init(&iter, &s.plocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(s.f, gkyl_range_idx(&s.plocal_ext, iter.idx));
    for (int k = 2 * s.np; k < 8 * s.np; ++k) {
      f[k] = 0.0; // the coupling blocks mirror with their own signs; F and G suffice here
    }
  }
  // phase modes of the 1x1v hybrid basis: [1, x, w, x w, w^2, x w^2]
  static const double mirror[6] = {1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
  for (int k = s.plocal.lower[1]; k <= s.plocal.upper[1]; ++k) {
    int skin[2] = {s.clocal.lower[0], k}, ghost[2] = {s.clocal.lower[0] - 1, s.nv + 1 - k};
    const double *fs = gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, skin));
    double *fg = gkyl_array_fetch(s.f, gkyl_range_idx(&s.plocal_ext, ghost));
    for (int c = 0; c < 8; ++c) {
      for (int m = 0; m < s.np; ++m) {
        fg[c * s.np + m] = c < 2 ? mirror[m] * fs[c * s.np + m] : 0.0;
      }
    }
  }
  gkyl_array_clear(s.cfl_f, 0.0);
  advance(&s);
  const double *ff = gkyl_array_cfetch(s.fluid_flux, gkyl_range_idx(&s.clocal_ext, s.clocal.lower));
  TEST_CHECK(fabs(ff[0]) < 1e-13);
  TEST_MSG("wall mass flux %g", ff[0]);
  teardown(&s);
}

// 5. Disjoint velocity grids (the initial state of Noh): the flux is two non-overlapping beams.
void
test_disjoint_beams(void)
{
  struct test_setup s;
  setup(&s, 1, 8, 6.0);
  double b[3] = {1.0, 0.0, 0.0}, u_l[3] = {1000.0, 0.0, 0.0}, u_r[3] = {-1000.0, 0.0, 0.0};
  int ix = s.clocal.lower[0] + 1, il = ix - 1;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    set_frame(
      &s, gkyl_range_idx(&s.clocal_ext, iter.idx), 1.0, iter.idx[0] < ix ? u_l : u_r, 1.0, b, 1.0
    );
  }
  srand(21);
  randomize_f(&s);
  gkyl_range_iter_init(&iter, &s.plocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(s.f, gkyl_range_idx(&s.plocal_ext, iter.idx));
    for (int k = 2 * s.np; k < 8 * s.np; ++k) {
      f[k] = 0.0; // the coupling blocks do not stream along b
    }
  }
  advance(&s);
  double blend =
    ((const double *)gkyl_array_cfetch(s.blend, gkyl_range_idx(&s.clocal_ext, &ix)))[0];
  TEST_CHECK(blend > 0.0 && blend < 1.0);

  static const double gn[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  static const double gw[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
  double dw = s.pgrid.dx[1], max_err = 0.0, max_val = 0.0, xc[2];
  for (int k = s.plocal.lower[1]; k <= s.plocal.upper[1]; ++k) {
    int pl[2] = {il, k}, pr[2] = {ix, k};
    gkyl_rect_grid_cell_center(&s.pgrid, pr, xc);
    const double *f_l = gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, pl));
    const double *f_r = gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, pr));
    const double *flux = gkyl_array_cfetch(s.flux, gkyl_range_idx(&s.plocal_ext, pr));
    for (int c = 0; c < 2; ++c) {
      double coef_stored[2][3];
      for (int side = 0; side < 2; ++side) {
        double cc[3];
        stored_to_legendre_1x(&flux[side * 8 * s.nk + c * s.nk], cc);
        for (int m = 0; m < 3; ++m) {
          coef_stored[side][m] = cc[m];
        }
      }
      for (int m = 0; m < 3; ++m) {
        double remap[2] = {0.0, 0.0}, label = 0.0;
        for (int j = 0; j < 3; ++j) {
          double eta = gn[j], wq = xc[1] + 0.5 * dw * eta;
          double legv[3] = {
            0.7071067811865475, 1.224744871391589 * eta,
            0.7905694150420949 * (3.0 * eta * eta - 1.0)
          };
          // the trace of the block at the face: f = sum_k f_k b_k, x = +1 (left), -1 (right)
          double xl =
            (0.7071067811865475 * f_l[c * s.np + 0] + 1.224744871391589 * f_l[c * s.np + 1]) *
              legv[0] +
            (0.7071067811865475 * f_l[c * s.np + 2] + 1.224744871391589 * f_l[c * s.np + 3]) *
              legv[1] +
            (0.7071067811865475 * f_l[c * s.np + 4] + 1.224744871391589 * f_l[c * s.np + 5]) *
              legv[2];
          double xr =
            (0.7071067811865475 * f_r[c * s.np + 0] - 1.224744871391589 * f_r[c * s.np + 1]) *
              legv[0] +
            (0.7071067811865475 * f_r[c * s.np + 2] - 1.224744871391589 * f_r[c * s.np + 3]) *
              legv[1] +
            (0.7071067811865475 * f_r[c * s.np + 4] - 1.224744871391589 * f_r[c * s.np + 5]) *
              legv[2];
          // face frame: wbar = w + 1000 for the left grid, w - 1000 for the right
          remap[0] += gw[j] * legv[m] * (wq + 1000.0) * xl;
          remap[1] += gw[j] * legv[m] * (wq - 1000.0) * xr;
          label += gw[j] * legv[m] * wq * (wq > 0.0 ? xl : xr);
        }
        for (int side = 0; side < 2; ++side) {
          double expected = (1.0 - blend) * remap[side] + blend * label;
          max_err = fmax(max_err, fabs(coef_stored[side][m] - expected));
          max_val = fmax(max_val, fabs(expected));
        }
      }
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_val);
  TEST_MSG("beams: %g of %g, blend %g", max_err, max_val, blend);
  teardown(&s);
}

// 6. A 1x problem run along each axis of a 3x grid reproduces 1x.
struct along_ctx {
  int dir; // the axis the problem varies along
};

static void
eval_frame_em(double t, const double *xn, double *fout, void *ctx)
{
  const struct along_ctx *ac = ctx;
  double x = xn[ac->dir];
  fout[0] = 1.0 + 0.6 * sin(M_PI * x); // rho
  for (int k = 0; k < 3; ++k) {
    fout[1 + k] = 0.0;
  }
  fout[1 + ac->dir] = 0.8 * cos(M_PI * x); // u_n
  fout[4] = 0.0;
}

static void
eval_frame_xi(double t, const double *xn, double *fout, void *ctx)
{
  const struct along_ctx *ac = ctx;
  fout[0] = 1.0 + 0.4 * cos(M_PI * xn[ac->dir]);
  fout[1] = fout[2] = 0.0;
}

static void
eval_f_along(double t, const double *xn, double *fout, void *ctx)
{
  const int *dirs = ctx; // [the axis the problem varies along, the velocity direction]
  double x = xn[dirs[0]], w = xn[dirs[1]];
  fout[0] = exp(-0.5 * w * w) * (1.0 + 0.3 * sin(M_PI * x));
  fout[1] = 0.5 * exp(-0.5 * w * w);
  for (int k = 2; k < 8; ++k) {
    fout[k] = 0.0;
  }
}

struct along_ctx_f {
  int dir;
  int vdir;
};

static void
fill_along(struct test_setup *s, int dir)
{
  struct along_ctx ac = {.dir = dir};
  struct along_ctx_f af = {.dir = dir, .vdir = s->cdim};
  gkyl_proj_on_basis *proj =
    gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 3, 5, eval_frame_em, &ac);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->fluid);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 3, 3, eval_frame_xi, &ac);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->xi);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_clear(s->u, 0.0);
  gkyl_array_clear(s->bvar, 0.0);
  gkyl_array_clear(s->theta, 0.0);
  double c0 = sqrt(pow(2.0, s->cdim));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s->clocal_ext, iter.idx);
    double *u = gkyl_array_fetch(s->u, loc), *b = gkyl_array_fetch(s->bvar, loc);
    double *th = gkyl_array_fetch(s->theta, loc);
    const double *fl = gkyl_array_cfetch(s->fluid, loc);
    for (int m = 0; m < s->nc; ++m) {
      u[dir * s->nc + m] = fl[(1 + dir) * s->nc + m]; // the normal flow, as projected above
    }
    b[dir * s->nc] = c0;
    th[0] = c0;
  }
  proj = gkyl_proj_on_basis_new(&s->pgrid, &s->pbasis, 3, 8, eval_f_along, &af);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal_ext, s->f);
  gkyl_proj_on_basis_release(proj);
}

static void
test_axis(int dir)
{
  const int nx1[1] = {4}, nx3[3][3] = {{4, 1, 1}, {1, 4, 1}, {1, 1, 4}};
  struct test_setup s1, s3;
  setup_cells(&s1, 1, 8, 6.0, nx1);
  setup_cells(&s3, 3, 8, 6.0, nx3[dir]);
  fill_along(&s1, 0);
  fill_along(&s3, dir);
  advance(&s1);
  advance(&s3);
  static const int surf_map[3] = {0, 3, 7}; // the 3x surface modes (a = 0, m)
  int nf1 = 8 * s1.nk, nf3 = 8 * s3.nk, trans[3] = {1, 1, 1};
  double max_err = 0.0, max_val = 0.0;
  for (int ix = s1.clocal.lower[0]; ix <= s1.clocal.upper[0] + 1; ++ix) {
    trans[dir] = ix;
    long c1 = gkyl_range_idx(&s1.clocal_ext, &ix), c3 = gkyl_range_idx(&s3.clocal_ext, trans);
    const double *b1 = gkyl_array_cfetch(s1.blend, c1), *b3 = gkyl_array_cfetch(s3.blend, c3);
    TEST_CHECK(gkyl_compare_double(b1[0], b3[dir], 1e-13));
    TEST_MSG("blend %g vs %g", b1[0], b3[dir]);
    const double *ff1 = gkyl_array_cfetch(s1.fluid_flux, c1);
    const double *ff3 = &((const double *)gkyl_array_cfetch(s3.fluid_flux, c3))[dir * 5 * s3.na];
    for (int m = 0; m < 5; ++m) {
      // the normal momentum is component x in 1x and component dir along that axis in 3x
      int m3 = m == 1 ? 1 + dir : (m == 1 + dir ? 1 : m);
      max_err = fmax(max_err, fabs(2.0 * ff1[m] - ff3[m3 * s3.na]));
      max_val = fmax(max_val, fabs(ff1[m]));
      for (int a = 1; a < s3.na; ++a) {
        max_err = fmax(max_err, fabs(ff3[m3 * s3.na + a]));
      }
    }
    for (int k = s1.plocal.lower[1]; k <= s1.plocal.upper[1]; ++k) {
      int p1[2] = {ix, k}, p3[4] = {1, 1, 1, k};
      p3[dir] = ix;
      const double *f1 = gkyl_array_cfetch(s1.flux, gkyl_range_idx(&s1.plocal_ext, p1));
      const double *f3 = gkyl_array_cfetch(s3.flux, gkyl_range_idx(&s3.plocal_ext, p3));
      for (int side = 0; side < 2; ++side) {
        for (int c = 0; c < 8; ++c) {
          for (int m = 0; m < 3; ++m) {
            double v1 = f1[side * nf1 + c * s1.nk + m];
            double v3 = f3[(2 * dir + side) * nf3 + c * s3.nk + surf_map[m]];
            max_err = fmax(max_err, fabs(4.0 * v1 - v3));
            max_val = fmax(max_val, fabs(v1));
          }
        }
      }
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_val);
  TEST_MSG("dir %d: %g of %g", dir, max_err, max_val);
  teardown(&s1);
  teardown(&s3);
}

void
test_axis_all(void)
{
  test_axis(0);
  test_axis(1);
  test_axis(2);
}

static void
eval_maxwellian_fg(double t, const double *xn, double *fout, void *ctx)
{
  int vdir = *(int *)ctx;
  double w = xn[vdir];
  for (int c = 0; c < 8; ++c) {
    fout[c] = 0.0;
  }
  fout[0] = fout[1] = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
}

// 7. A gas at rest with uniform pressure and a density that is linear across the cells (1x): the
// thermal speed at the Gauss-Legendre nodes is sqrt(p/rho), and the face value from 1/xi^2 along
// the node line is sqrt(p/rho_face) on both sides of every face. So no mass or energy crosses a
// face and the momentum flux, the pressure, is the same on every face.
void
test_uniform_pressure(void)
{
  struct test_setup s;
  setup(&s, 1, 16, 8.0);
  const double pres = 0.7, rho_mid = 1.0, rho_slope = 0.45, node = 1.0 / sqrt(3.0);
  double u[3] = {0.0, 0.0, 0.0}, b[3] = {1.0, 0.0, 0.0}, xc[1];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.clocal_ext, iter.idx);
    gkyl_rect_grid_cell_center(&s.cgrid, iter.idx, xc);
    double rho_cell = rho_mid + rho_slope * xc[0], rho_half = 0.5 * rho_slope * s.cgrid.dx[0];
    set_frame(&s, loc, rho_cell, u, 1.0, b, 1.0);
    double *fl = gkyl_array_fetch(s.fluid, loc), *xl = gkyl_array_fetch(s.xi, loc);
    // modes [1/sqrt(2), sqrt(3/2) x]: a function with node values (lo, up) at x = -+1/sqrt(3)
    // has the coefficients ((up + lo)/sqrt(2), (up - lo)/sqrt(2))
    fl[1] = rho_half * sqrt(2.0 / 3.0);
    fl[4 * s.nc] = 1.5 * pres * sqrt(2.0);
    double xi_lo = sqrt(pres / (rho_cell - node * rho_half));
    double xi_up = sqrt(pres / (rho_cell + node * rho_half));
    xl[0] = (xi_up + xi_lo) / sqrt(2.0);
    xl[1] = (xi_up - xi_lo) / sqrt(2.0);
    xl[s.nc] = (xi_up * xi_up + xi_lo * xi_lo) / sqrt(2.0);
    xl[s.nc + 1] = (xi_up * xi_up - xi_lo * xi_lo) / sqrt(2.0);
    xl[2 * s.nc] = pres * sqrt(2.0);
  }
  int vdir = s.cdim;
  gkyl_proj_on_basis *proj =
    gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 4, 8, eval_maxwellian_fg, &vdir);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext, s.f);
  gkyl_proj_on_basis_release(proj);
  advance(&s);

  double mom_ref = 0.0, max_mass = 0.0, max_energy = 0.0, max_dev = 0.0;
  int idx[1];
  for (idx[0] = s.clocal.lower[0]; idx[0] <= s.clocal.upper[0] + 1; ++idx[0]) {
    const double *ff = gkyl_array_cfetch(s.fluid_flux, gkyl_range_idx(&s.clocal_ext, idx));
    if (idx[0] == s.clocal.lower[0]) {
      mom_ref = ff[1 * s.na];
    }
    max_mass = fmax(max_mass, fabs(ff[0]));
    max_energy = fmax(max_energy, fabs(ff[4 * s.na]));
    max_dev = fmax(max_dev, fabs(ff[1 * s.na] - mom_ref));
  }
  // the momentum flux is the pressure (the second moment of the projected Maxwellian is 1 to 1e-5)
  TEST_CHECK(fabs(mom_ref / pres - 1.0) < 1e-4);
  TEST_CHECK(max_dev < 1e-13 * fabs(mom_ref));
  TEST_CHECK(max_mass < 1e-13 * fabs(mom_ref));
  TEST_CHECK(max_energy < 1e-13 * fabs(mom_ref));
  TEST_MSG(
    "momentum flux %.15g (p = %g), deviation %g, mass flux %g, energy flux %g", mom_ref, pres,
    max_dev, max_mass, max_energy
  );
  teardown(&s);
}

// The time-step rates of configuration cell ix of a 1x setup: the kinetic rate of each of its
// velocity cells and the fluid rate.
static void
get_rates_1x(const struct test_setup *s, int ix, double *rate_kin, double *rate_fluid)
{
  for (int k = 0; k < s->nv; ++k) {
    int pidx[2] = {ix, s->plocal.lower[1] + k};
    rate_kin[k] =
      ((const double *)gkyl_array_cfetch(s->cfl_f, gkyl_range_idx(&s->plocal_ext, pidx)))[0];
  }
  rate_fluid[0] =
    ((const double *)gkyl_array_cfetch(s->cfl_fluid, gkyl_range_idx(&s->clocal_ext, &ix)))[0];
}

// Frame of a cell: density, flow and thermal speed.
struct rate_frame {
  double rho, u[3], xi;
};

// Rates of the cells of a 1x line whose frame is uniform (frame 0), or jumps from frame 0 to
// frame 1 at the upper face of cell num_lo.
static void
rates_of_line(
  int num_cells, int num_lo, const struct rate_frame *frames, double *rate_kin, double *rate_fluid
)
{
  struct test_setup s;
  const int nx[3] = {num_cells, 1, 1};
  setup_cells(&s, 1, 8, 6.0, nx);
  const double b[3] = {1.0, 0.0, 0.0};
  for (int ix = 0; ix <= num_cells + 1; ++ix) {
    const struct rate_frame *fr = ix <= num_lo ? &frames[0] : &frames[1];
    set_frame(&s, gkyl_range_idx(&s.clocal_ext, &ix), fr->rho, fr->u, fr->xi, b, 1.0);
  }
  advance(&s);
  for (int ix = 1; ix <= num_cells; ++ix) {
    get_rates_1x(&s, ix, &rate_kin[(ix - 1) * s.nv], &rate_fluid[ix - 1]);
  }
  teardown(&s);
}

// A cell takes, in each direction, the larger of the rates of its two faces: with a jump of the
// frame at one face of a line, the two cells of that face see its rate, and every other cell the
// rate of its own uniform frame.
static void
test_rate_larger_face(void)
{
  enum { NUM_CELLS = 6, NUM_LO = 3, NV = 8 };
  struct rate_frame frames[2] = {
    {.rho = 1.0, .u = {0.4, 0.0, 0.0}, .xi = 1.0}, {.rho = 0.5, .u = {-2.6, 0.0, 0.0}, .xi = 0.7}
  };
  struct rate_frame frames_lo[2] = {frames[0], frames[0]}, frames_up[2] = {frames[1], frames[1]};
  double kin_lo[NUM_CELLS * NV], fluid_lo[NUM_CELLS], kin_up[NUM_CELLS * NV], fluid_up[NUM_CELLS];
  double kin[NUM_CELLS * NV], fluid[NUM_CELLS];
  rates_of_line(NUM_CELLS, NUM_CELLS + 1, frames_lo, kin_lo, fluid_lo);
  rates_of_line(NUM_CELLS, NUM_CELLS + 1, frames_up, kin_up, fluid_up);
  rates_of_line(NUM_CELLS, NUM_LO, frames, kin, fluid);

  bool jump_is_faster = false;
  for (int ix = 1; ix <= NUM_CELLS; ++ix) {
    // Cells away from the jump: the rate of the uniform frame on their side.
    if (ix != NUM_LO && ix != NUM_LO + 1) {
      const double *kin_ref = ix < NUM_LO ? kin_lo : kin_up;
      const double *fluid_ref = ix < NUM_LO ? fluid_lo : fluid_up;
      for (int k = 0; k < NV; ++k) {
        TEST_CHECK(kin[(ix - 1) * NV + k] == kin_ref[k]);
        TEST_MSG(
          "cell %d, velocity cell %d: %g against %g", ix, k, kin[(ix - 1) * NV + k], kin_ref[k]
        );
      }
      TEST_CHECK(fluid[ix - 1] == fluid_ref[0]);
    }
  }
  // The two cells of the jump: max(uniform lower, jump) and max(jump, uniform upper), so the
  // larger of each with the other side's uniform rate is the same number.
  for (int k = 0; k < NV; ++k) {
    double below = kin[(NUM_LO - 1) * NV + k], above = kin[NUM_LO * NV + k];
    TEST_CHECK(below >= kin_lo[k] && above >= kin_up[k]);
    TEST_CHECK(fmax(below, kin_up[k]) == fmax(above, kin_lo[k]));
    TEST_MSG(
      "velocity cell %d: below %g, above %g (uniform %g, %g)", k, below, above, kin_lo[k], kin_up[k]
    );
    jump_is_faster = jump_is_faster || below > kin_lo[k];
  }
  TEST_CHECK(fluid[NUM_LO - 1] >= fluid_lo[0] && fluid[NUM_LO] >= fluid_up[0]);
  TEST_CHECK(fmax(fluid[NUM_LO - 1], fluid_up[0]) == fmax(fluid[NUM_LO], fluid_lo[0]));
  // The jump face is faster than the uniform faces of the cell below it, so the check above
  // does test that this cell sees its upper face.
  TEST_CHECK(jump_is_faster && fluid[NUM_LO - 1] > fluid_lo[0]);
}

// The rates of a cell do not depend on the range the updater is given: updating the two halves
// of the domain separately (as two ranks do) gives the rates of updating all of it at once.
static void
test_rate_decomposition(int cdim, int split_dir)
{
  struct test_setup s;
  const int nx[3] = {4, 4, 2};
  setup_cells(&s, cdim, 4, 6.0, nx);
  double b[3] = {0.6, 0.64, 0.48};
  if (cdim == 1) {
    b[0] = 1.0;
    b[1] = b[2] = 0.0;
  }
  // A different frame in every cell (ghost cells included).
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.clocal_ext, iter.idx);
    double u[3] = {sin(1.7 * loc), 0.5 * cos(0.9 * loc), 0.3 * sin(2.3 * loc + 1.0)};
    set_frame(&s, loc, 1.0 + 0.5 * sin(0.7 * loc), u, 0.8 + 0.3 * cos(1.3 * loc), b, 1.0);
  }
  advance(&s);
  struct gkyl_array *cfl_f_all = gkyl_array_clone(s.cfl_f);
  struct gkyl_array *cfl_fluid_all = gkyl_array_clone(s.cfl_fluid);

  gkyl_array_clear(s.cfl_f, 0.0);
  gkyl_array_clear(s.cfl_fluid, 0.0);
  int num_lo = nx[split_dir] / 2;
  for (int half = 0; half < 2; ++half) {
    int clo[GKYL_MAX_DIM], cup[GKYL_MAX_DIM], plo[GKYL_MAX_DIM], pup[GKYL_MAX_DIM];
    for (int d = 0; d < s.pdim; ++d) {
      plo[d] = s.plocal.lower[d];
      pup[d] = s.plocal.upper[d];
    }
    plo[split_dir] = half == 0 ? 1 : num_lo + 1;
    pup[split_dir] = half == 0 ? num_lo : nx[split_dir];
    for (int d = 0; d < cdim; ++d) {
      clo[d] = plo[d];
      cup[d] = pup[d];
    }
    struct gkyl_range crange, prange;
    gkyl_sub_range_init(&crange, &s.clocal, clo, cup);
    gkyl_sub_range_init(&prange, &s.plocal, plo, pup);
    gkyl_dg_pkpm_remap_flux_surf_advance(
      s.up, &crange, &s.clocal_ext, &prange, &s.plocal_ext, s.bvar, s.fluid, s.u, s.xi, s.theta,
      s.f, s.flux, s.fluid_flux, s.rate_f, s.rate_fluid, s.cfl_f, s.cfl_fluid, s.blend
    );
  }

  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.plocal_ext, iter.idx);
    double all = ((const double *)gkyl_array_cfetch(cfl_f_all, loc))[0];
    double halves = ((const double *)gkyl_array_cfetch(s.cfl_f, loc))[0];
    TEST_CHECK(all == halves && all > 0.0);
    TEST_MSG("kinetic rate of phase cell %ld: %.17g at once, %.17g in halves", loc, all, halves);
  }
  gkyl_range_iter_init(&iter, &s.clocal);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.clocal_ext, iter.idx);
    double all = ((const double *)gkyl_array_cfetch(cfl_fluid_all, loc))[0];
    double halves = ((const double *)gkyl_array_cfetch(s.cfl_fluid, loc))[0];
    TEST_CHECK(all == halves && all > 0.0);
    TEST_MSG("fluid rate of cell %ld: %.17g at once, %.17g in halves", loc, all, halves);
  }
  gkyl_array_release(cfl_f_all);
  gkyl_array_release(cfl_fluid_all);
  teardown(&s);
}

static void
test_rate_decomposition_all(void)
{
  test_rate_decomposition(1, 0);
  test_rate_decomposition(2, 0);
  test_rate_decomposition(2, 1);
  test_rate_decomposition(3, 2);
}

static void
pmap_quadratic(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = xn[0] + 0.3 * xn[0] * xn[0];
}

// On a position map the rate of a face takes the narrower of its two cells: with a uniform frame
// (every face has the same speeds) the rate of a cell is that of the identity map over the
// smallest Jacobian of the cell and its two neighbors (a ghost cell has that of its skin cell).
static void
test_rate_narrower_cell(void)
{
  enum { NUM_CELLS = 6, NV = 8 };
  const int nx[3] = {NUM_CELLS, 1, 1};
  const double u[3] = {0.4, 0.0, 0.0}, b[3] = {1.0, 0.0, 0.0};
  double kin_identity[NV], fluid_identity[1];
  struct test_setup s;
  setup_cells(&s, 1, NV, 6.0, nx);
  set_uniform_frame(&s, 1.0, u, 1.0, b, 1.0);
  advance(&s);
  get_rates_1x(&s, 1, kin_identity, fluid_identity);
  teardown(&s);

  setup_cells_map(&s, 1, NV, 6.0, nx, pmap_quadratic);
  set_uniform_frame(&s, 1.0, u, 1.0, b, 1.0);
  advance(&s);
  for (int ix = 1; ix <= NUM_CELLS; ++ix) {
    double jacob_min = DBL_MAX;
    for (int j = ix - 1; j <= ix + 1; ++j) {
      const double *jacob =
        gkyl_array_cfetch(s.pos_map->jacob_pos, gkyl_range_idx(&s.clocal_ext, &j));
      jacob_min = fmin(jacob_min, jacob[0]);
    }
    double kin[NV], fluid[1];
    get_rates_1x(&s, ix, kin, fluid);
    for (int k = 0; k < NV; ++k) {
      TEST_CHECK(gkyl_compare(kin[k], kin_identity[k] / jacob_min, 1e-14));
      TEST_MSG(
        "cell %d, velocity cell %d: %g against %g", ix, k, kin[k], kin_identity[k] / jacob_min
      );
    }
    TEST_CHECK(gkyl_compare(fluid[0], fluid_identity[0] / jacob_min, 1e-14));
  }
  // The map is not uniform, so the check above distinguishes the cells.
  const int first = 1, last = NUM_CELLS;
  double jacob_first =
    ((const double *)gkyl_array_cfetch(s.pos_map->jacob_pos, gkyl_range_idx(&s.clocal_ext, &first))
    )[0];
  double jacob_last =
    ((const double *)gkyl_array_cfetch(s.pos_map->jacob_pos, gkyl_range_idx(&s.clocal_ext, &last))
    )[0];
  TEST_CHECK(jacob_last > 2.0 * jacob_first);
  teardown(&s);
}

// ---- The kernels of one transverse node (gkyl_pkpm_remap_flux_kernels.h).

typedef void (*interval_kern_t)(
  const double *fc, double w_lo, double w_up, const double *side_l, const double *side_r,
  const double *cell_l, const double *cell_r, const double *Y_l, const double *Y_r, double *P_l,
  double *P_r, double *mom, double *share
);
typedef void (*label_kern_t)(
  const double *fc, double w_center, double dw, const double *X_l, const double *X_r, double *proj,
  double *mom_l, double *mom_r, double *speed
);
typedef void (*to_face_kern_t)(
  const double *psi, const double *fwd, const double *G, double *X, double *Y
);
typedef void (*face_store_kern_t)(
  const double *wpsi, const double *bck, const double *P, double *flux
);
typedef void (*store_kern_t)(const double *wpsi, const double *proj, double *flux);

static const interval_kern_t interval_kernels[3] = {
  pkpm_remap_flux_interval_x_1v_p2, pkpm_remap_flux_interval_y_1v_p2,
  pkpm_remap_flux_interval_z_1v_p2
};
static const label_kern_t label_kernels[3] = {
  pkpm_remap_flux_label_x_1v_p2, pkpm_remap_flux_label_y_1v_p2, pkpm_remap_flux_label_z_1v_p2
};
static const to_face_kern_t to_face_kernels[3] = {
  pkpm_remap_flux_to_face_1x1v_ser_p1, pkpm_remap_flux_to_face_2x1v_ser_p1,
  pkpm_remap_flux_to_face_3x1v_ser_p1
};
static const face_store_kern_t from_face_kernels[3][3] = {
  {pkpm_remap_flux_from_face_x_1x1v_ser_p1, NULL, NULL},
  {pkpm_remap_flux_from_face_x_2x1v_ser_p1, pkpm_remap_flux_from_face_y_2x1v_ser_p1, NULL},
  {pkpm_remap_flux_from_face_x_3x1v_ser_p1, pkpm_remap_flux_from_face_y_3x1v_ser_p1,
   pkpm_remap_flux_from_face_z_3x1v_ser_p1}
};
static const store_kern_t store_kernels[3][3] = {
  {pkpm_remap_flux_store_x_1x1v_ser_p1, NULL, NULL},
  {pkpm_remap_flux_store_x_2x1v_ser_p1, pkpm_remap_flux_store_y_2x1v_ser_p1, NULL},
  {pkpm_remap_flux_store_x_3x1v_ser_p1, pkpm_remap_flux_store_y_3x1v_ser_p1,
   pkpm_remap_flux_store_z_3x1v_ser_p1}
};

static double
rand_unit(void)
{
  return 2.0 * rand() / RAND_MAX - 1.0;
}

// The factors of the Legendre modes of the velocity basis, phi_m = norm_m [1, eta, 3 eta^2 - 1].
static void
vel_mode_norm(double *norm)
{
  struct gkyl_basis vel_basis;
  gkyl_cart_modal_serendip(&vel_basis, 1, 2);
  double eta = 1.0, modes[3];
  vel_basis.eval(&eta, modes);
  norm[0] = modes[0];
  norm[1] = modes[1];
  norm[2] = 0.5 * modes[2];
}

// Face coefficients of a face of direction dir whose field direction is oblique to the normal,
// with different frames on the two sides: the coupled speed c is positive and the HLL weights of
// the flux have kinks at kinks[0] < kinks[1], between which both sides contribute.
static void
oblique_face_coef(int dir, double *fc, double *kinks)
{
  double b_l[3] = {0.62, 0.55, 0.56}, b_r[3] = {0.57, 0.6, 0.56};
  pkpm_harm_face_coef(dir, 0.3, -0.1, 1.1, 0.9, b_l, b_r, 0.8, 1.2, 1, fc);
  double a0 = fc[PKPM_FC_A0], a1 = fc[PKPM_FC_A1], c = fc[PKPM_FC_C];
  TEST_CHECK(c > 0.0 && fabs(a1) > 0.0);
  kinks[0] = fmin((-c - a0) / a1, (c - a0) / a1);
  kinks[1] = fmax((-c - a0) / a1, (c - a0) / a1);
}

// 9. The interval kernel is the exact integral of the flux function against the Legendre modes
// of the cell of each side, and of its face-frame moments, on an interval between the kinks of
// the HLL weights, in each of the three regimes (the lower upwind side, both sides, the upper
// upwind side): compared with a 5-point rule of pkpm_harm_conf_flux. The shares are those of
// the 3 Gauss-Legendre points, counted as outside for a side with cell[PKPM_RF_CELL_OUT] = 1.
static void
test_interval_kernel_dir(int dir)
{
  srand(31 + dir);
  double norm[3];
  vel_mode_norm(norm);
  double fc[PKPM_FC_NUM], kinks[2];
  oblique_face_coef(dir, fc, kinks);
  double span = kinks[1] - kinks[0];
  const double ends[3][2] = {
    {kinks[0] - 1.3, kinks[0] - 0.2},
    {kinks[0] + 0.1 * span, kinks[1] - 0.15 * span},
    {kinks[1] + 0.2, kinks[1] + 1.5},
  };
  // maps of the sides (r, c), and the cell of each side (center, width)
  const double map[2][2] = {{1.2, 0.3}, {0.85, -0.2}}, vcell[2][2] = {{0.4, 1.5}, {-0.3, 1.1}};
  double side[2][PKPM_RF_SIDE_NUM], cell[2][PKPM_RF_CELL_NUM];
  for (int s = 0; s < 2; ++s) {
    side[s][PKPM_RF_SIDE_R] = map[s][0];
    side[s][PKPM_RF_SIDE_C] = map[s][1];
    side[s][PKPM_RF_SIDE_KINV] = 0.9 + 0.2 * s;
    side[s][PKPM_RF_SIDE_EG] = 1.3 - 0.4 * s;
    side[s][PKPM_RF_SIDE_CC] = 0.2 + 0.1 * s;
    for (int k = 0; k < 3; ++k) {
      side[s][PKPM_RF_SIDE_CPERP + k] = 0.3 * rand_unit();
    }
    cell[s][PKPM_RF_CELL_ETA0] = 2.0 * (map[s][1] - vcell[s][0]) / vcell[s][1];
    cell[s][PKPM_RF_CELL_ETA1] = 2.0 * map[s][0] / vcell[s][1];
  }
  cell[0][PKPM_RF_CELL_OUT] = 1.0;
  cell[1][PKPM_RF_CELL_OUT] = 0.0;

  for (int regime = 0; regime < 3; ++regime) {
    double Y[2][24], P[2][24] = {{0.0}}, mom[PKPM_RF_NUM_MOMS] = {0.0};
    double share[PKPM_RF_SHARE_NUM] = {0.0};
    for (int s = 0; s < 2; ++s) {
      for (int k = 0; k < 24; ++k) {
        Y[s][k] = rand_unit();
      }
    }
    double lo = ends[regime][0], up = ends[regime][1];
    interval_kernels[dir](
      fc, lo, up, side[0], side[1], cell[0], cell[1], Y[0], Y[1], P[0], P[1], mom, share
    );

    // reference: nq-point rule of the flux function; the shares at the 3 points of the kernel
    double P_ref[2][24] = {{0.0}}, mom_ref[PKPM_RF_NUM_MOMS] = {0.0};
    double share_ref[PKPM_RF_SHARE_NUM] = {0.0};
    for (int pass = 0; pass < 2; ++pass) {
      int nq = pass == 0 ? 5 : 3;
      for (int q = 0; q < nq; ++q) {
        double wb = 0.5 * (up + lo) + 0.5 * (up - lo) * gkyl_gauss_ordinates[nq][q];
        double om = 0.5 * (up - lo) * gkyl_gauss_weights[nq][q];
        double X[2][8], eta[2], phi_l[8], phi_r[8], sb[2], flux[8];
        for (int s = 0; s < 2; ++s) {
          eta[s] = 2.0 * (map[s][0] * wb + map[s][1] - vcell[s][0]) / vcell[s][1];
          for (int c = 0; c < 8; ++c) {
            X[s][c] =
              Y[s][c] + Y[s][8 + c] * eta[s] + Y[s][16 + c] * (3.0 * eta[s] * eta[s] - 1.0);
          }
        }
        pkpm_harm_conf_flux(dir, X[0], X[1], wb, fc, 1, phi_l, phi_r, sb);
        for (int c = 0; c < 8; ++c) {
          flux[c] = phi_l[c] + phi_r[c];
        }
        if (pass == 0) {
          for (int s = 0; s < 2; ++s) {
            double leg[3] = {norm[0], norm[1] * eta[s], norm[2] * (3.0 * eta[s] * eta[s] - 1.0)};
            for (int m = 0; m < 3; ++m) {
              for (int c = 0; c < 8; ++c) {
                P_ref[s][m * 8 + c] += om * 2.0 * map[s][0] / vcell[s][1] * leg[m] * flux[c];
              }
            }
          }
          mom_ref[0] += om * flux[0];
          mom_ref[1] += om * wb * flux[0];
          mom_ref[2] += om * wb * wb * flux[0];
          mom_ref[3] += om * flux[1];
          for (int k = 0; k < 3; ++k) {
            mom_ref[4 + k] += om * flux[2 + k];
          }
        } else {
          for (int s = 0; s < 2; ++s) {
            const double *ci = &side[s][PKPM_RF_SIDE_CPERP];
            double w_s = map[s][0] * wb + map[s][1];
            double cM = ci[0] * flux[2] + ci[1] * flux[3] + ci[2] * flux[4];
            double y_f = side[s][PKPM_RF_SIDE_KINV] * flux[0];
            double y_g = side[s][PKPM_RF_SIDE_EG] *
                         (flux[1] - cM + 0.5 * side[s][PKPM_RF_SIDE_CC] * flux[0]);
            share_ref[PKPM_RF_SHARE_MOM_TOT + s] += om * fabs(flux[0]) * fabs(w_s);
            share_ref[PKPM_RF_SHARE_EN_TOT + s] +=
              om * (fabs(y_f) * w_s * w_s + 2.0 * fabs(y_g));
          }
        }
      }
    }
    share_ref[PKPM_RF_SHARE_MOM_OUT] = share_ref[PKPM_RF_SHARE_MOM_TOT];
    share_ref[PKPM_RF_SHARE_EN_OUT] = share_ref[PKPM_RF_SHARE_EN_TOT];

    double max_err = 0.0, max_val = 0.0, max_share_err = 0.0, max_share = 0.0;
    for (int s = 0; s < 2; ++s) {
      for (int k = 0; k < 24; ++k) {
        max_err = fmax(max_err, fabs(P[s][k] - P_ref[s][k]));
        max_val = fmax(max_val, fabs(P_ref[s][k]));
      }
    }
    for (int k = 0; k < PKPM_RF_NUM_MOMS; ++k) {
      max_err = fmax(max_err, fabs(mom[k] - mom_ref[k]));
      max_val = fmax(max_val, fabs(mom_ref[k]));
    }
    for (int k = 0; k < PKPM_RF_SHARE_NUM; ++k) {
      max_share_err = fmax(max_share_err, fabs(share[k] - share_ref[k]));
      max_share = fmax(max_share, fabs(share_ref[k]));
    }
    TEST_CHECK(max_val > 0.1);
    TEST_CHECK(max_err < 1e-13 * max_val);
    TEST_MSG("dir %d regime %d: projection and moments off by %g of %g", dir, regime, max_err, max_val);
    TEST_CHECK(max_share > 0.0);
    TEST_CHECK(max_share_err < 1e-13 * max_share);
    TEST_MSG("dir %d regime %d: shares off by %g of %g", dir, regime, max_share_err, max_share);
    TEST_CHECK(share[PKPM_RF_SHARE_MOM_OUT + 1] == 0.0 && share[PKPM_RF_SHARE_EN_OUT + 1] == 0.0);
  }
}

void
test_interval_kernel(void)
{
  for (int dir = 0; dir < 3; ++dir) {
    test_interval_kernel_dir(dir);
  }
}

// 10. The label-paired kernel of a velocity cell is the interval kernel on that cell with the
// identity map on both sides (the same three points, whether or not a kink falls in the cell):
// the same projection for both sides and the same moments, the sum of the two contributions.
// Its speeds are the bounds of the flux function at the three points, and are only raised.
static void
test_label_kernel_dir(int dir)
{
  srand(57 + dir);
  double fc[PKPM_FC_NUM], kinks[2];
  oblique_face_coef(dir, fc, kinks);
  // a cell that contains the lower kink
  double w_center = kinks[0] + 0.2, dw = 1.5;
  double X[2][24];
  for (int s = 0; s < 2; ++s) {
    for (int k = 0; k < 24; ++k) {
      X[s][k] = rand_unit();
    }
  }
  double proj[24], mom_side[2][PKPM_RF_NUM_MOMS] = {{0.0}}, speed[2] = {0.0, 0.0};
  label_kernels[dir](fc, w_center, dw, X[0], X[1], proj, mom_side[0], mom_side[1], speed);

  double side[PKPM_RF_SIDE_NUM] = {0.0};
  side[PKPM_RF_SIDE_R] = 1.0;
  side[PKPM_RF_SIDE_KINV] = 1.0;
  side[PKPM_RF_SIDE_EG] = 1.0;
  double cell[PKPM_RF_CELL_NUM] = {-2.0 * w_center / dw, 2.0 / dw, 0.0};
  double P[2][24] = {{0.0}}, mom[PKPM_RF_NUM_MOMS] = {0.0}, share[PKPM_RF_SHARE_NUM] = {0.0};
  interval_kernels[dir](
    fc, w_center - 0.5 * dw, w_center + 0.5 * dw, side, side, cell, cell, X[0], X[1], P[0], P[1],
    mom, share
  );

  double max_err = 0.0, max_val = 0.0;
  for (int k = 0; k < 24; ++k) {
    max_err = fmax(max_err, fmax(fabs(proj[k] - P[0][k]), fabs(proj[k] - P[1][k])));
    max_val = fmax(max_val, fabs(proj[k]));
  }
  for (int k = 0; k < PKPM_RF_NUM_MOMS; ++k) {
    max_err = fmax(max_err, fabs(mom_side[0][k] + mom_side[1][k] - mom[k]));
    max_val = fmax(max_val, fabs(mom[k]));
  }
  TEST_CHECK(max_val > 0.1);
  TEST_CHECK(max_err < 1e-13 * max_val);
  TEST_MSG("dir %d: label-paired kernel off the interval kernel by %g of %g", dir, max_err, max_val);

  double speed_ref[2] = {0.0, 0.0};
  for (int q = 0; q < 3; ++q) {
    double wq = w_center + 0.5 * gkyl_gauss_ordinates[3][q] * dw;
    double x_l[8] = {0.0}, x_r[8] = {0.0}, phi_l[8], phi_r[8], sb[2];
    pkpm_harm_conf_flux(dir, x_l, x_r, wq, fc, 1, phi_l, phi_r, sb);
    speed_ref[0] = fmax(speed_ref[0], fmax(-sb[0], sb[1]));
    speed_ref[1] = fmax(speed_ref[1], sb[1] - sb[0]);
  }
  TEST_CHECK(speed_ref[0] > 0.0 && speed_ref[1] > 0.0);
  TEST_CHECK(fabs(speed[0] - speed_ref[0]) < 1e-14 * speed_ref[0]);
  TEST_CHECK(fabs(speed[1] - speed_ref[1]) < 1e-14 * speed_ref[1]);
  TEST_MSG("dir %d: speeds %g %g, expected %g %g", dir, speed[0], speed[1], speed_ref[0], speed_ref[1]);
  double speed_high[2] = {10.0 * speed_ref[0], 10.0 * speed_ref[1]};
  label_kernels[dir](fc, w_center, dw, X[0], X[1], proj, mom_side[0], mom_side[1], speed_high);
  TEST_CHECK(speed_high[0] == 10.0 * speed_ref[0] && speed_high[1] == 10.0 * speed_ref[1]);
}

void
test_label_kernel(void)
{
  for (int dir = 0; dir < 3; ++dir) {
    test_label_kernel_dir(dir);
  }
}

// 11. The kernels that take one velocity cell of a side to the face frame and back. The blocks
// at the node are the transverse interpolation of the modal traces, and with the factors of the
// identity map the face frame is the frame of the side. With the two halves of a face map the
// way to the face frame and back is the identity, so that the store of the returned blocks is
// the store of the blocks themselves: for all blocks when the frames differ in density, thermal
// speed and flow along the field direction, and for F, G and M0_perp when they also differ in
// the flow across it (the row of N_perp replaces a second moment by G and closes only to first
// order in that shift, so it is left out).
static void
test_face_frame_kernels_cdim(int cdim)
{
  srand(73 + cdim);
  int na = 1 << (cdim - 1), nt = 3 * na;
  double norm[3];
  vel_mode_norm(norm);
  double G[8 * 12], psi[4], wpsi[4];
  for (int k = 0; k < 8 * nt; ++k) {
    G[k] = rand_unit();
  }
  for (int a = 0; a < na; ++a) {
    psi[a] = rand_unit();
    wpsi[a] = rand_unit();
  }

  double c_zero[3] = {0.0, 0.0, 0.0}, identity[PKPM_FRAME_NUM];
  pkpm_harm_frame_factors(1.0, 1.0, c_zero, identity);
  double X[24], Y[24];
  to_face_kernels[cdim - 1](psi, identity, G, X, Y);
  double max_err = 0.0, max_val = 0.0;
  for (int m = 0; m < 3; ++m) {
    for (int c = 0; c < 8; ++c) {
      double ref = 0.0;
      for (int a = 0; a < na; ++a) {
        ref += norm[m] * psi[a] * G[c * nt + a * 3 + m];
      }
      max_err = fmax(max_err, fabs(X[m * 8 + c] - ref));
      max_val = fmax(max_val, fabs(ref));
      TEST_CHECK(Y[m * 8 + c] == X[m * 8 + c]);
    }
  }
  TEST_CHECK(max_val > 0.1);
  TEST_CHECK(max_err < 1e-14 * max_val);
  TEST_MSG("%dx: blocks at the node off by %g of %g", cdim, max_err, max_val);

  double b_bar[3] = {0.6, 0.0, 0.8}, u_s[3] = {0.3, -0.2, 0.5};
  // flow of the face frame: shifted along b_bar (all blocks), then also across it (F, G, M0_perp)
  const double u_bar[2][3] = {{0.3 + 0.4 * 0.6, -0.2, 0.5 + 0.4 * 0.8}, {0.1, 0.25, 0.35}};
  const int num_blocks[2] = {8, 5};
  for (int shift = 0; shift < 2; ++shift) {
    struct pkpm_harm_face_map map;
    pkpm_harm_face_map(1.3, u_s, 0.9, 1.1, u_bar[shift], 1.2, b_bar, &map);
    double fwd[PKPM_FRAME_NUM], bck[PKPM_FRAME_NUM];
    pkpm_harm_frame_factors(map.K, map.r, map.c_perp, fwd);
    pkpm_harm_frame_factors(map.K_inv, map.r_inv, map.c_perp_inv, bck);
    to_face_kernels[cdim - 1](psi, fwd, G, X, Y);
    double frame_diff = 0.0;
    for (int k = 0; k < 24; ++k) {
      frame_diff = fmax(frame_diff, fabs(Y[k] - X[k]));
    }
    TEST_CHECK(frame_diff > 0.01);
    for (int dir = 0; dir < cdim; ++dir) {
      double flux_back[8 * 12] = {0.0}, flux_side[8 * 12] = {0.0};
      from_face_kernels[cdim - 1][dir](wpsi, bck, Y, flux_back);
      store_kernels[cdim - 1][dir](wpsi, X, flux_side);
      max_err = 0.0, max_val = 0.0;
      for (int k = 0; k < num_blocks[shift] * nt; ++k) {
        max_err = fmax(max_err, fabs(flux_back[k] - flux_side[k]));
        max_val = fmax(max_val, fabs(flux_side[k]));
      }
      TEST_CHECK(max_val > 0.1);
      TEST_CHECK(max_err < 1e-13 * max_val);
      TEST_MSG(
        "%dx dir %d shift %d: to the face frame and back off by %g of %g", cdim, dir, shift,
        max_err, max_val
      );
    }
  }
}

void
test_face_frame_kernels(void)
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    test_face_frame_kernels_cdim(cdim);
  }
}

TEST_LIST = {
  {"test_identity_1x", test_identity_1x},
  {"test_identity_2x", test_identity_2x},
  {"test_identity_3x", test_identity_3x},
  {"test_steady_dilation", test_steady_dilation},
  {"test_steady_shift", test_steady_shift},
  {"test_moments", test_moments},
  {"test_reflecting_wall", test_reflecting_wall},
  {"test_disjoint_beams", test_disjoint_beams},
  {"test_axis_all", test_axis_all},
  {"test_uniform_pressure", test_uniform_pressure},
  {"test_rate_larger_face", test_rate_larger_face},
  {"test_rate_decomposition", test_rate_decomposition_all},
  {"test_rate_narrower_cell", test_rate_narrower_cell},
  {"test_interval_kernel", test_interval_kernel},
  {"test_label_kernel", test_label_kernel},
  {"test_face_frame_kernels", test_face_frame_kernels},
  {NULL, NULL}
};
