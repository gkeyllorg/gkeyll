// Tests of the PKPM configuration-space HLL surface flux of the eight blocks [F, G, M0_perp,
// N_perp] (dg_pkpm_conf_flux_surf, plan 17.14) against an independent nodal evaluation, and of the
// stored-flux scatter of the configuration pass (dg_pkpm_kinetic_conf) through hyper_dg:
// conservation on a periodic domain and the complete stage-A right-hand side.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// The flux both sides of every face receive, from the single-valued label-paired flux:
// [dir][mode] -> [dir][side][mode].
static struct gkyl_array *
two_sided_copy(const struct gkyl_array *flux, int cdim, int nflux)
{
  struct gkyl_array *sides = gkyl_array_new(GKYL_DOUBLE, 2 * cdim * nflux, flux->size);
  for (long i = 0; i < flux->size; ++i) {
    const double *f = gkyl_array_cfetch(flux, i);
    double *g = gkyl_array_fetch(sides, i);
    for (int d = 0; d < cdim; ++d) {
      for (int m = 0; m < nflux; ++m) {
        g[2 * d * nflux + m] = f[d * nflux + m];
        g[(2 * d + 1) * nflux + m] = f[d * nflux + m];
      }
    }
  }
  return sides;
}

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

// Periodic smooth fields on [-pi, pi]^cdim so that projection on the extended
// range fills the ghost cells with their periodic values.
static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.5 * sin(y);
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.2 * sin(x);
}

static void
eval_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  // Comparable to the thermal speed so that a_l and a_r change sign across nodes.
  fout[0] = 0.8 * sin(x) + 0.3 * cos(y);
  fout[1] = -0.7 * cos(y) * sin(z) + 0.2 * sin(x);
  fout[2] = 0.5 * sin(z) - 0.4 * cos(x) * cos(y);
}

// Thermal speed [xi, xi^2, p] and theta (positive, varying).
static void
eval_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double xi = 1.0 + 0.2 * cos(x) * sin(y) + 0.05 * sin(z);
  fout[0] = xi;
  fout[1] = xi * xi;
  fout[2] = 0.0;
}

static void
eval_theta(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.7 + 0.2 * sin(x + y) + 0.05 * cos(z);
}

static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double v = xn[cdim];
  double n = 1.0 + 0.3 * sin(x) * cos(y) + 0.1 * cos(z);
  double vt = 1.0 + 0.2 * cos(x);
  double F0 = n / sqrt(2.0 * M_PI * vt * vt) * exp(-v * v / (2.0 * vt * vt));
  fout[0] = F0 * (1.0 + 0.1 * v); // F_0 (a nonzero first moment is fine here)
  fout[1] = 0.5 * vt * vt * F0 * (1.0 - 0.1 * v);
  // generic first-harmonic blocks (not transverse: the flux is defined for any traces)
  fout[2] = F0 * (0.3 * sin(x) + 0.1 * v);
  fout[3] = F0 * (-0.2 * cos(y) + 0.15 * v * v - 0.1);
  fout[4] = F0 * (0.25 * cos(z) - 0.05 * v);
  fout[5] = 0.5 * vt * vt * F0 * (0.2 * cos(x) - 0.1 * v);
  fout[6] = 0.5 * vt * vt * F0 * (0.1 + 0.2 * sin(y) * v);
  fout[7] = 0.5 * vt * vt * F0 * (-0.3 * sin(x + z) + 0.05 * v * v);
}

// Face rule of the kinetic flux: 5 Gauss-Legendre points per transverse configuration direction,
// 3 in w.
static const double gl5_x[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
  0.9061798459386640};
static const double gl5_w[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
  0.4786286704993665, 0.2369268850561891};
static const double gl3_x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double gl3_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

// Physical flux J(X) of the eight blocks through a face with normal e_d (notes, eq.
// kinetic-face-flux), written out: J_F = V F + xi M0_d, J_G = V G + xi N_d,
// J_M0 = V M0 + xi G (e_d - b b_d), J_N = V N + xi (4 theta G - 2 theta^2 F) (e_d - b b_d).
static void
flux_ref(int d, double V, double xi, const double *b, double theta, const double *X, double *J)
{
  J[0] = V * X[0] + xi * X[2 + d];
  J[1] = V * X[1] + xi * X[5 + d];
  for (int a = 0; a < 3; ++a) {
    double pn = (a == d ? 1.0 : 0.0) - b[a] * b[d];
    J[2 + a] = V * X[2 + a] + xi * X[1] * pn;
    J[5 + a] = V * X[5 + a] + xi * (4.0 * theta * X[1] - 2.0 * theta * theta * X[0]) * pn;
  }
}

// Reference HLL flux of the eight blocks, long form (upwind branches, then the HLL average), of
// the flux function with face-averaged coefficients (plan 17.41): J(X) = a X + T(X) with
// a = {u_d} + {xi} b_d w and T with {xi}, b = {b}/|{b}|, {theta}, {q} = (q_l + q_r)/2, applied to each
// side's own state; the bounds of the flux are the averaged speeds a -+ c and zero,
// c = xi |b x e_d| sqrt((2 + sqrt 2) theta) with the averaged coefficients (plan 17.42). Returns the
// time-step speed: the largest of -S, S over the envelope of both sides' speeds V -+ c, a -+ c and 0.
static double
hll_ref(int d, const double *X_l, const double *X_r, double w, double u_l, double u_r, double xi_l,
  double xi_r, const double *b_l, const double *b_r, double th_l, double th_r, double *ghat)
{
  double ba[3], xa = 0.5 * (xi_l + xi_r), tha = 0.5 * (th_l + th_r);
  for (int a = 0; a < 3; ++a) {
    ba[a] = 0.5 * (b_l[a] + b_r[a]);
  }
  double bmag = sqrt(ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]);
  for (int a = 0; a < 3; ++a) {
    ba[a] = bmag > 0.0 ? ba[a] / bmag : 0.0;
  }
  const double *bs[3] = {b_l, b_r, ba};
  double xs[3] = {xi_l, xi_r, xa}, ths[3] = {th_l, th_r, tha}, us[3] = {u_l, u_r, 0.5 * (u_l + u_r)};
  double s_l = 0.0, s_r = 0.0, e[3] = {0.0}, env_l = 0.0, env_r = 0.0;
  e[d] = 1.0;
  for (int k = 0; k < 3; ++k) {
    double bx[3];
    for (int a = 0; a < 3; ++a) { // b x e_d
      bx[a] = bs[k][(a + 1) % 3] * e[(a + 2) % 3] - bs[k][(a + 2) % 3] * e[(a + 1) % 3];
    }
    double h = sqrt(bx[0] * bx[0] + bx[1] * bx[1] + bx[2] * bx[2]);
    double c = xs[k] * h * sqrt((2.0 + sqrt(2.0)) * ths[k]), V = us[k] + xs[k] * bs[k][d] * w;
    env_l = fmin(env_l, V - c);
    env_r = fmax(env_r, V + c);
    if (k == 2) {
      s_l = fmin(V - c, 0.0);
      s_r = fmax(V + c, 0.0);
    }
  }
  double J_l[8], J_r[8], a = us[2] + xa * ba[d] * w;
  flux_ref(d, a, xa, ba, tha, X_l, J_l);
  flux_ref(d, a, xa, ba, tha, X_r, J_r);
  for (int c = 0; c < 8; ++c) {
    if (s_l == 0.0 && s_r == 0.0) {
      ghat[c] = 0.5 * (J_l[c] + J_r[c]);
    } else if (s_l >= 0.0) {
      ghat[c] = J_l[c];
    } else if (s_r <= 0.0) {
      ghat[c] = J_r[c];
    } else {
      ghat[c] = (s_r * J_l[c] - s_l * J_r[c] + s_l * s_r * (X_r[c] - X_l[c])) / (s_r - s_l);
    }
  }
  return fmax(-env_l, env_r);
}

struct test_setup {
  struct gkyl_vlasov_position_map *pos_map; // identity configuration-space map
  int cdim, pdim;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext;
  struct gkyl_array *bvar, *u, *xi, *theta, *f, *cflrate, *flux, *face_coef;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
  struct gkyl_dg_pkpm_conf_flux_surf *calc_flux;
};

static void
setup(struct test_setup *s, int cdim)
{
  int pdim = cdim + 1;
  s->cdim = cdim;
  s->pdim = pdim;
  double lower[4] = {-M_PI, -M_PI, -M_PI, -4.0}, upper[4] = {M_PI, M_PI, M_PI, 4.0};
  int cells[4] = {6, 4, 3, 8};
  int ccells[3] = {6, 4, 3};
  int cghost[3] = {1, 1, 1}, pghost[4] = {1, 1, 1, 0};
  if (cdim < 3) {
    lower[cdim] = -4.0;
    upper[cdim] = 4.0;
    cells[cdim] = 8;
    pghost[cdim] = 0;
  }
  gkyl_rect_grid_init(&s->cgrid, cdim, lower, upper, ccells);
  gkyl_rect_grid_init(&s->pgrid, pdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  s->pos_map = identity_pos_map(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis);
  gkyl_create_grid_ranges(&s->pgrid, pghost, &s->plocal_ext, &s->plocal);

  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  struct gkyl_array *em = mkarr(8 * nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  s->xi = mkarr(3 * nc, s->clocal_ext.volume);
  s->theta = mkarr(nc, s->clocal_ext.volume);
  s->f = mkarr(8 * np, s->plocal_ext.volume);
  gkyl_array_clear(s->f, 0.0);
  s->cflrate = mkarr(1, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_xi, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 1, eval_theta, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->theta);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->pgrid, &s->pbasis, 8, 8, eval_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal_ext, s->f);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, em, s->bvar);
  gkyl_array_release(em);

  struct gkyl_dg_pkpm_conf_flux_surf_inp inp = {
    .phase_grid = &s->pgrid,
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .pos_map = s->pos_map, .use_gpu = false,
  };
  s->calc_flux = gkyl_dg_pkpm_conf_flux_surf_inew(&inp);
  int nk = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s->calc_flux);
  TEST_CHECK(nk == np / 2);
  s->flux = mkarr(8 * cdim * nk, s->plocal_ext.volume);
  s->face_coef = mkarr(gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp(s->calc_flux), s->clocal_ext.volume);
  gkyl_array_clear(s->flux, 0.0);
  gkyl_array_clear(s->cflrate, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_face_coef(s->calc_flux, &s->clocal, s->bvar, s->u, s->xi, s->theta, s->face_coef);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    s->calc_flux, &s->clocal, &s->plocal, &s->plocal_ext, s->face_coef, s->f, s->cflrate, s->flux
  );
}

static void
teardown(struct test_setup *s)
{
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_dg_pkpm_conf_flux_surf_release(s->calc_flux);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->flux);
  gkyl_array_release(s->face_coef);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->f);
  gkyl_array_release(s->u);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->bvar);
}

// Independent reference for the stored flux on the face between cells idx_l and
// idx_r in direction dir (the lower face of idx_r). Returns the maximum speed bound.
static double
reference_face(
  struct test_setup *s, int dir, const int *idx_l, const int *idx_r, const struct gkyl_basis *sbasis,
  double *flux_ref
)
{
  int cdim = s->cdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis, nk = sbasis->num_basis;
  long cl = gkyl_range_idx(&s->clocal_ext, idx_l), cr = gkyl_range_idx(&s->clocal_ext, idx_r);
  long pl = gkyl_range_idx(&s->plocal_ext, idx_l), pr = gkyl_range_idx(&s->plocal_ext, idx_r);
  const double *b_l = gkyl_array_cfetch(s->bvar, cl), *b_r = gkyl_array_cfetch(s->bvar, cr);
  const double *u_l = gkyl_array_cfetch(s->u, cl), *u_r = gkyl_array_cfetch(s->u, cr);
  const double *xi_l = gkyl_array_cfetch(s->xi, cl), *xi_r = gkyl_array_cfetch(s->xi, cr);
  const double *th_l = gkyl_array_cfetch(s->theta, cl), *th_r = gkyl_array_cfetch(s->theta, cr);
  const double *f_l = gkyl_array_cfetch(s->f, pl), *f_r = gkyl_array_cfetch(s->f, pr);
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx_r, xc);
  double wv = xc[cdim], dv = s->pgrid.dx[cdim];

  for (int k = 0; k < 8 * nk; ++k) {
    flux_ref[k] = 0.0;
  }
  double amax = 0.0;
  // Transverse configuration directions (5 nodes each), then w (3 nodes).
  int ntr = cdim - 1;
  int nnodes = 3;
  for (int d = 0; d < ntr; ++d) {
    nnodes *= 5;
  }
  for (int n = 0; n < nnodes; ++n) {
    int sub[4];
    int m = n;
    sub[ntr] = m % 3; // w fastest
    m /= 3;
    for (int d = ntr - 1; d >= 0; --d) {
      sub[d] = m % 5;
      m /= 5;
    }
    double w = 1.0, zs[4], zl[4], zr[4];
    int t = 0;
    for (int d = 0; d < cdim; ++d) {
      if (d == dir) {
        zl[d] = 1.0;
        zr[d] = -1.0;
      } else {
        zl[d] = zr[d] = zs[t] = gl5_x[sub[t]];
        w *= gl5_w[sub[t]];
        t += 1;
      }
    }
    zl[cdim] = zr[cdim] = zs[t] = gl3_x[sub[t]];
    w *= gl3_w[sub[t]];
    double v = wv + 0.5 * dv * zs[t];

    double bl[3], br[3], X_l[8], X_r[8], ghat[8];
    for (int a = 0; a < 3; ++a) {
      bl[a] = s->cbasis.eval_expand(zl, &b_l[a * nc]);
      br[a] = s->cbasis.eval_expand(zr, &b_r[a * nc]);
    }
    double xl = s->cbasis.eval_expand(zl, xi_l), xr = s->cbasis.eval_expand(zr, xi_r);
    double ul = s->cbasis.eval_expand(zl, &u_l[dir * nc]), ur = s->cbasis.eval_expand(zr, &u_r[dir * nc]);
    for (int c = 0; c < 8; ++c) {
      X_l[c] = s->pbasis.eval_expand(zl, &f_l[c * np]);
      X_r[c] = s->pbasis.eval_expand(zr, &f_r[c * np]);
    }
    amax = fmax(amax, hll_ref(dir, X_l, X_r, v, ul, ur, xl, xr, bl, br, s->cbasis.eval_expand(zl, th_l),
      s->cbasis.eval_expand(zr, th_r), ghat));
    for (int c = 0; c < 8; ++c) {
      for (int k = 0; k < nk; ++k) {
        // Surface basis function k at the surface node, divided by its value at the
        // emitter's shape-normalization point (2/5 in every variable): the stored
        // flux is the shape coefficient, and the scatter maps fold the constant back in.
        double unit[64] = {0.0};
        unit[k] = 1.0;
        double znorm[4] = {0.4, 0.4, 0.4, 0.4};
        double chi_k = sbasis->eval_expand(zs, unit) / sbasis->eval_expand(znorm, unit);
        flux_ref[c * nk + k] += w * chi_k * ghat[c];
      }
    }
  }
  return amax;
}

static void
test_conf_flux(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nk = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s.calc_flux);

  struct gkyl_basis sbasis;
  if (cdim == 1) {
    gkyl_cart_modal_serendip(&sbasis, 1, 2);
  } else {
    gkyl_cart_modal_hybrid(&sbasis, cdim - 1, 1);
  }
  TEST_CHECK(sbasis.num_basis == nk);

  double flux_ref[8 * 32];
  double max_flux = 0.0, max_err = 0.0, max_cfl_err = 0.0, max_flux_h = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    const double *flux = gkyl_array_cfetch(s.flux, pidx);
    double cfl_ref = 0.0;
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_l[4], idx_r[4];
      gkyl_copy_int_arr(s.pdim, iter.idx, idx_l);
      idx_l[dir] -= 1;
      double amax = reference_face(&s, dir, idx_l, iter.idx, &sbasis, flux_ref);
      cfl_ref += 1.5 * 2.0 / s.pgrid.dx[dir] * amax;
      for (int k = 0; k < 8 * nk; ++k) {
        double err = fabs(flux[dir * 8 * nk + k] - flux_ref[k]);
        max_err = err > max_err ? err : max_err;
        max_flux = fabs(flux_ref[k]) > max_flux ? fabs(flux_ref[k]) : max_flux;
        if (k >= 2 * nk) {
          max_flux_h = fabs(flux_ref[k]) > max_flux_h ? fabs(flux_ref[k]) : max_flux_h;
        }
      }
      // Upper ghost-owned face.
      if (iter.idx[dir] == s.plocal.upper[dir]) {
        gkyl_copy_int_arr(s.pdim, iter.idx, idx_r);
        idx_r[dir] += 1;
        long pidx_r = gkyl_range_idx(&s.plocal_ext, idx_r);
        const double *flux_g = gkyl_array_cfetch(s.flux, pidx_r);
        amax = reference_face(&s, dir, iter.idx, idx_r, &sbasis, flux_ref);
        cfl_ref += 1.5 * 2.0 / s.pgrid.dx[dir] * amax;
        for (int k = 0; k < 8 * nk; ++k) {
          double err = fabs(flux_g[dir * 8 * nk + k] - flux_ref[k]);
          max_err = err > max_err ? err : max_err;
        }
      }
    }
    const double *cfl = gkyl_array_cfetch(s.cflrate, pidx);
    double cerr = fabs(cfl[0] - cfl_ref);
    max_cfl_err = cerr > max_cfl_err ? cerr : max_cfl_err;
  }
  TEST_CHECK(max_flux > 0.0 && max_flux_h > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_flux);
  TEST_MSG("cdim=%d: max flux error %g (max flux %g, first harmonic %g)", cdim, max_err, max_flux, max_flux_h);
  TEST_CHECK(max_cfl_err < 1e-11);
  TEST_MSG("cdim=%d: max cflrate error %g", cdim, max_cfl_err);

  teardown(&s);
}

// hyper_dg over the configuration directions with the scatter surface term and no
// volume term: the increments telescope, so the cell-averaged RHS sums to zero on
// the periodic domain for both components, while individual cells are nonzero.
static void
test_scatter(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int np = s.pbasis.num_basis;

  // the surface term reads the flux each side of a face receives
  struct gkyl_array *flux_sides = two_sided_copy(s.flux, cdim, 8 * s.pbasis.num_basis / 2);

  struct gkyl_dg_pkpm_kinetic_conf_inp inp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.bvar = s.bvar, .u = s.u, .xi = s.xi, .conf_flux_surf = flux_sides},
    .pos_map = s.pos_map,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_pkpm_kinetic_conf_inew(&inp);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn, cdim, up_dirs, zero_flux, 1, false);

  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(s.cflrate, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.plocal, s.f, s.cflrate, rhs);

  // With zero gradients every block has only divergence volume terms (streaming), so the cell
  // averages change through the scattered face fluxes only, which telescope on the periodic grid.
  for (int c = 0; c < 8; ++c) {
    double sum = 0.0, mag = 0.0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &s.plocal);
    while (gkyl_range_iter_next(&iter)) {
      const double *r = gkyl_array_cfetch(rhs, gkyl_range_idx(&s.plocal_ext, iter.idx));
      sum += r[c * np];
      mag += fabs(r[c * np]);
    }
    TEST_CHECK(mag > 0.0);
    TEST_CHECK(fabs(sum) < 1e-12 * mag);
    TEST_MSG("cdim=%d: block %d sum %g, |sum| %g", cdim, c, sum, mag);
  }

  gkyl_array_release(rhs);
  gkyl_array_release(flux_sides);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}

// The full stage-A right-hand side (hyper_dg: volume kernel + scattered face fluxes)
// against the volume kernel plus an independent DG surface term built from the stored
// fluxes of every configuration direction:
//   dx1_d [ <psi_m(xi_d = -1), Fhat_l> - <psi_m(xi_d = +1), Fhat_r> ]
// with Fhat = sum_k c_k chi_k over the face variables (the other configuration
// variables and v_par), c_k = g_k chi_k(2/5, ...) (shape to orthonormal coefficient).
static void
test_scatter_ref(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int pdim = s.pdim, np = s.pbasis.num_basis, nk = np / 2;
  struct gkyl_basis sbasis;
  if (cdim == 1) {
    gkyl_cart_modal_serendip(&sbasis, 1, 2);
  } else {
    gkyl_cart_modal_hybrid(&sbasis, cdim - 1, 1);
  }
  TEST_CHECK(sbasis.num_basis == nk);

  // the surface term reads the flux each side of a face receives
  struct gkyl_array *flux_sides = two_sided_copy(s.flux, cdim, 8 * s.pbasis.num_basis / 2);

  struct gkyl_dg_pkpm_kinetic_conf_inp inp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.bvar = s.bvar, .u = s.u, .xi = s.xi, .conf_flux_surf = flux_sides},
    .pos_map = s.pos_map,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_pkpm_kinetic_conf_inew(&inp);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn, cdim, up_dirs, zero_flux, 1, false);
  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(s.cflrate, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.plocal, s.f, s.cflrate, rhs);

  // Volume kernels of the three dimensionalities.
  double max_err = 0.0, max_rhs = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    double xc[4];
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    double ref[8 * 48] = {0.0};
    const double *bv = gkyl_array_cfetch(s.bvar, cidx), *uv = gkyl_array_cfetch(s.u, cidx);
    const double *fv = gkyl_array_cfetch(s.f, pidx), *xv = gkyl_array_cfetch(s.xi, cidx);
    const double *jac = gkyl_array_cfetch(s.pos_map->jacob_pos, cidx); // identity map: ones
    if (cdim == 1) {
      pkpm_kinetic_conf_vol_1x1v_ser_p1(xc, s.pgrid.dx, jac, bv, uv, xv, fv, ref);
    } else if (cdim == 2) {
      pkpm_kinetic_conf_vol_2x1v_ser_p1(xc, s.pgrid.dx, jac, bv, uv, xv, fv, ref);
    } else {
      pkpm_kinetic_conf_vol_3x1v_ser_p1(xc, s.pgrid.dx, jac, bv, uv, xv, fv, ref);
    }
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_r[4];
      gkyl_copy_int_arr(pdim, iter.idx, idx_r);
      idx_r[dir] += 1;
      long pidx_r = gkyl_range_idx(&s.plocal_ext, idx_r);
      const double *g_l = gkyl_array_cfetch(s.flux, pidx), *g_r = gkyl_array_cfetch(s.flux, pidx_r);
      double dx1 = 2.0 / s.pgrid.dx[dir];
      int nnodes = 1;
      for (int d = 0; d < cdim; ++d) {
        nnodes *= 3;
      }
      for (int c = 0; c < 8; ++c) {
        for (int m = 0; m < np; ++m) {
          double unit_m[64] = {0.0};
          unit_m[m] = 1.0;
          double surf = 0.0;
          for (int n = 0; n < nnodes; ++n) {
            // Surface nodes: the other configuration variables then v_par, 3 each.
            int sub[4], q = n;
            for (int d = cdim - 1; d >= 0; --d) {
              sub[d] = q % 3;
              q /= 3;
            }
            double w = 1.0, zs[4], zl[4], zr[4];
            int t = 0;
            for (int d = 0; d < cdim; ++d) {
              if (d == dir) {
                zl[d] = 1.0;
                zr[d] = -1.0;
              } else {
                zl[d] = zr[d] = zs[t] = gl3_x[sub[t]];
                w *= gl3_w[sub[t]];
                t += 1;
              }
            }
            zl[cdim] = zr[cdim] = zs[t] = gl3_x[sub[t]];
            w *= gl3_w[sub[t]];
            double fl = 0.0, fr = 0.0;
            for (int k = 0; k < nk; ++k) {
              double unit_k[64] = {0.0};
              unit_k[k] = 1.0;
              double znorm[4] = {0.4, 0.4, 0.4, 0.4};
              double chi_norm = sbasis.eval_expand(znorm, unit_k), chi = sbasis.eval_expand(zs, unit_k);
              fl += g_l[dir * 8 * nk + c * nk + k] * chi_norm * chi;
              fr += g_r[dir * 8 * nk + c * nk + k] * chi_norm * chi;
            }
            // psi_m at the lower face of this cell (lower flux enters with +) and at its
            // upper face (upper flux enters with -).
            double zlo[4], zhi[4];
            for (int d = 0; d < pdim; ++d) {
              zlo[d] = zr[d];
              zhi[d] = zl[d];
            }
            surf += w * (s.pbasis.eval_expand(zlo, unit_m) * fl - s.pbasis.eval_expand(zhi, unit_m) * fr);
          }
          ref[c * np + m] += dx1 * surf;
        }
      }
    }
    const double *r = gkyl_array_cfetch(rhs, pidx);
    for (int k = 0; k < 8 * np; ++k) {
      max_err = fmax(max_err, fabs(r[k] - ref[k]));
      max_rhs = fmax(max_rhs, fabs(ref[k]));
    }
  }
  TEST_CHECK(max_rhs > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_rhs);
  TEST_MSG("cdim=%d: stage-A rhs vs reference: max error %g (max rhs %g)", cdim, max_err, max_rhs);

  gkyl_array_release(rhs);
  gkyl_array_release(flux_sides);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}


void test_scatter_ref_1x() { test_scatter_ref(1); }
void test_scatter_ref_2x() { test_scatter_ref(2); }
void test_scatter_ref_3x() { test_scatter_ref(3); }

void test_conf_flux_1x() { test_conf_flux(1); }
void test_conf_flux_2x() { test_conf_flux(2); }
void test_conf_flux_3x() { test_conf_flux(3); }
void test_scatter_1x() { test_scatter(1); }
void test_scatter_2x() { test_scatter(2); }
void test_scatter_3x() { test_scatter(3); }

TEST_LIST = {
  {"test_conf_flux_1x", test_conf_flux_1x},
  {"test_conf_flux_2x", test_conf_flux_2x},
  {"test_conf_flux_3x", test_conf_flux_3x},
  {"test_scatter_ref_1x", test_scatter_ref_1x},
  {"test_scatter_ref_2x", test_scatter_ref_2x},
  {"test_scatter_ref_3x", test_scatter_ref_3x},
  {"test_scatter_1x", test_scatter_1x},
  {"test_scatter_2x", test_scatter_2x},
  {"test_scatter_3x", test_scatter_3x},
  {NULL, NULL},
};
