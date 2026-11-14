// Tests of the PKPM configuration-space face flux (dg_pkpm_remap_flux_surf), the exact remap of
// physical_jump_flux.tex blended with the label-paired flux.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
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
  struct gkyl_array *bvar, *fluid, *u, *xi, *theta, *f, *flux, *fluid_flux, *cfl_f, *cfl_fluid,
    *blend;
  struct gkyl_dg_pkpm_remap_flux_surf *up;
};

static void
setup_cells(struct test_setup *s, int cdim, int nv, double wmax, const int *nx)
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
  s->cfl_f = mkarr(1, pvol);
  s->cfl_fluid = mkarr(1, cvol);
  s->blend = mkarr(cdim, cvol);
  gkyl_array_clear(s->f, 0.0);
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
    s->theta, s->f, s->flux, s->fluid_flux, s->cfl_f, s->cfl_fluid, s->blend
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
        double leg[3] = {
          sqrt(0.5), sqrt(1.5) * gn[q], sqrt(0.625) * (3.0 * gn[q] * gn[q] - 1.0)
        };
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
  for (int corner = 0; corner < (1 << cdim); ++corner) {
    double z[4];
    for (int d = 0; d < cdim; ++d) {
      z[d] = (corner >> d) & 1 ? 0.5773502691896258 : -0.5773502691896258;
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

// 3. Every cell holds the same physical distribution in its own frame: the flux each cell receives
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
  // xi = 1 and 2: the edges of the xi = 2 grid are a subset of those of the xi = 1 grid in v
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

// 4. With no flux outside either grid, each side's kinetic flux has the moments of the fluid flux
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

// 5. A reflecting wall (mirror ghost: u_n reversed, f(x, w) -> f(-x, -w)) has no normal mass flux.
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

// 6. Disjoint velocity grids (the initial state of Noh): the flux is two non-overlapping beams.
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

// 7. A 1x problem run along each axis of a 3x grid reproduces 1x.
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

// 8. A gas at rest with uniform pressure and a density that is linear across the cells (1x): the
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
  {NULL, NULL}
};
