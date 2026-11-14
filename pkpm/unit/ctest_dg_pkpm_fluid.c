// Tests of the momentum solve of the refactored PKPM solver: the exact modal volume kernel (with
// the first-harmonic stress and heat flux) against quadrature, and the face flux
// (the fluid flux of dg_pkpm_remap_flux_surf, the moments of the eight-block kinetic HLL flux)
// against an independent reference at interior faces, an absorbing wall and a reflecting wall.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

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

struct field_ctx {
  int cdim;
  bool uniform_rho; // for the legacy comparison
  bool b_along_x;
  bool zero_u; // for the legacy comparison
  bool inadmissible; // make p_par negative in part of the domain to exercise the HLL fallback
};

static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
  if (c->b_along_x) {
    fout[3] = 2.0;
    return;
  }
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.5 * sin(y);
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.2 * sin(x);
}

// Physical moments [rho, p_par, p_perp, M1, q_par, q_perp]
static void
eval_moms(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = c->uniform_rho ? 2.0 : 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * sin(z);
  fout[1] = 1.0 + 0.3 * cos(x) * sin(y);
  fout[2] = 0.7 + 0.25 * cos(y) * sin(x) + 0.1 * cos(z);
  fout[3] = 0.0;
  fout[4] = 0.3 * sin(x) * cos(z) + 0.1 * cos(y);
  fout[5] = -0.2 * cos(x) + 0.15 * sin(y) * sin(z);
  if (c->inadmissible && x > 1.0) {
    fout[1] = -0.2; // negative parallel pressure: P_nn <= 0 where b is along the normal
    fout[2] = 0.05;
  }
}

// Normalized moments [1, M2, theta, 0, M3, M1_G] of the same state with p = (p_par + 2 p_perp)/3
// and xi^2 = p/rho: M2 = p_par/p, theta = p_perp/p, and the third moments scaled by 1/(p xi).
static void
eval_nmoms(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6];
  eval_moms(t, xn, m, ctx);
  double p = (m[1] + 2.0 * m[2]) / 3.0, xi = sqrt(fabs(p) / m[0]);
  fout[0] = 1.0;
  fout[1] = m[1] / p;
  fout[2] = m[2] / p;
  fout[3] = 0.0;
  fout[4] = m[4] / (p * xi);
  fout[5] = m[5] / (p * xi);
}

// Thermal speed and scalar pressure [xi, xi^2, p] of the same state.
static void
eval_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6];
  eval_moms(t, xn, m, ctx);
  double p = (m[1] + 2.0 * m[2]) / 3.0;
  fout[1] = p / m[0];
  fout[0] = sqrt(fabs(fout[1]));
  fout[2] = p;
}

// Momentum M = rho u for a smooth u that stays positive along x (no stagnation switch).
static void
eval_euler(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double rho = c->uniform_rho ? 2.0 : 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * sin(z);
  fout[0] = rho * (1.5 + 0.3 * sin(x) + 0.2 * cos(y));
  fout[1] = rho * (0.4 * cos(y) * sin(z) + 0.2 * sin(x));
  fout[2] = rho * (-0.3 * sin(z) + 0.25 * cos(x) * cos(y));
  if (c->zero_u) {
    fout[0] = fout[1] = fout[2] = 0.0;
  }
}

// Fluid state [rho, rho ux, rho uy, rho uz, E] of the energy-constrained update.
static void
eval_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6], M[3];
  eval_euler(t, xn, M, ctx);
  eval_moms(t, xn, m, ctx);
  fout[0] = m[0];
  fout[1] = M[0];
  fout[2] = M[1];
  fout[3] = M[2];
  fout[4] = 0.5 * (M[0] * M[0] + M[1] * M[1] + M[2] * M[2]) / m[0] + 0.5 * (m[1] + 2.0 * m[2]);
}

// Kinetic coupling moments [p_par, p_perp, q_par, q_perp, s (3), r (3)] of the same state: s the
// parallel-perpendicular stress and r the perpendicular heat flux of the first harmonic.
static void
eval_cmoms(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double m[6];
  eval_moms(t, xn, m, ctx);
  fout[0] = m[1];
  fout[1] = m[2];
  fout[2] = m[4];
  fout[3] = m[5];
  fout[4] = 0.1 * sin(x) * cos(y);
  fout[5] = -0.15 * cos(z) + 0.05 * sin(y);
  fout[6] = 0.08 * cos(x + y);
  fout[7] = 0.2 * sin(y) * cos(z);
  fout[8] = -0.1 * cos(x) + 0.05 * sin(z);
  fout[9] = 0.12 * sin(x - y);
}


static const double gl3_x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double gl3_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
static const double gl4_x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
  0.8611363115940526};
static const double gl4_w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
  0.3478548451374538};

struct test_setup {
  struct gkyl_vlasov_position_map *pos_map; // identity configuration-space map
  int cdim;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_array *moms, *cmoms, *euler, *nmoms, *fluid, *bvar, *u, *xi;
  struct gkyl_dg_calc_pkpm_vars *vars;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
};

static void
setup(struct test_setup *s, struct field_ctx *ctx)
{
  int cdim = ctx->cdim;
  s->cdim = cdim;
  double lower[3] = {-M_PI, -M_PI, -M_PI}, upper[3] = {M_PI, M_PI, M_PI};
  int cells[3] = {8, 4, 3};
  int ghost[3] = {1, 1, 1};
  gkyl_rect_grid_init(&s->grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->basis, cdim, 1);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  s->pos_map = identity_pos_map(&s->grid, &s->local, &s->local_ext, &s->basis);
  int nc = s->basis.num_basis;

  s->moms = mkarr(6 * nc, s->local_ext.volume);
  s->cmoms = mkarr(10 * nc, s->local_ext.volume);
  s->euler = mkarr(3 * nc, s->local_ext.volume);
  s->nmoms = mkarr(6 * nc, s->local_ext.volume);
  s->fluid = mkarr(5 * nc, s->local_ext.volume);
  s->xi = mkarr(3 * nc, s->local_ext.volume);

  s->bvar = mkarr(3 * nc, s->local_ext.volume);
  s->u = mkarr(3 * nc, s->local_ext.volume);
  struct gkyl_array *em = mkarr(8 * nc, s->local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 6, eval_moms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->moms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 10, eval_cmoms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->cmoms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 3, eval_euler, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->euler);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 6, eval_nmoms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->nmoms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 5, eval_fluid, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->fluid);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 3, eval_xi, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 8, eval_em, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, em);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->grid, &s->basis, s->pos_map, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->local_ext, em, s->bvar);
  gkyl_array_release(em);

  s->vars = gkyl_dg_calc_pkpm_vars_new(&s->grid, &s->basis, s->pos_map, false);
  // the flow M/rho by the nodal quotient (the weak division at p=1); xi is set above
  struct gkyl_array *xi_scratch = mkarr(3 * nc, s->local_ext.volume);
  gkyl_dg_calc_pkpm_vars_prim(s->vars, &s->local_ext, s->fluid, s->u, xi_scratch);
  gkyl_array_release(xi_scratch);
}

static void
teardown(struct test_setup *s)
{
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_dg_calc_pkpm_vars_release(s->vars);
  gkyl_array_release(s->u);
  gkyl_array_release(s->bvar);
  gkyl_array_release(s->euler);
  gkyl_array_release(s->moms);
  gkyl_array_release(s->cmoms);
  gkyl_array_release(s->nmoms);
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->xi);
}

static struct gkyl_dg_eqn *
make_eqn(struct test_setup *s, const struct gkyl_array *flux)
{
  struct gkyl_dg_pkpm_fluid_inp inp = {
    .conf_basis = &s->basis,
    .conf_range = &s->local,
    .auxfields = {.vlasov_pkpm_moms = s->cmoms, .bvar = s->bvar, .u = s->u, .fluid_flux_surf = flux},
    .pos_map = s->pos_map, .use_gpu = false,
  };
  return gkyl_dg_pkpm_fluid_inew(&inp);
}

// Pressure tensor P = p_perp I + (p_par - p_perp) b b + s b + b s, heat flux
// q = b (q_par + 2 q_perp)/2 + r, density and energy at reference point z of cell cidx from the p=1
// factors.
static void
pressure_at(
  const struct test_setup *s, long cidx, const double *z, double P[3][3], double q[3], double *rho,
  double *E
)
{
  int nc = s->basis.num_basis;
  const double *m = gkyl_array_cfetch(s->moms, cidx), *b = gkyl_array_cfetch(s->bvar, cidx);
  const double *fl = gkyl_array_cfetch(s->fluid, cidx), *cm = gkyl_array_cfetch(s->cmoms, cidx);
  double bq[3], sq[3], rq[3];
  for (int i = 0; i < 3; ++i) {
    bq[i] = s->basis.eval_expand(z, &b[i * nc]);
    sq[i] = s->basis.eval_expand(z, &cm[(4 + i) * nc]);
    rq[i] = s->basis.eval_expand(z, &cm[(7 + i) * nc]);
  }
  *rho = s->basis.eval_expand(z, &fl[0]);
  *E = s->basis.eval_expand(z, &fl[4 * nc]);
  double ppar = s->basis.eval_expand(z, &m[nc]), pperp = s->basis.eval_expand(z, &m[2 * nc]);
  double h = s->basis.eval_expand(z, &m[4 * nc]) + 2.0 * s->basis.eval_expand(z, &m[5 * nc]);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      P[i][j] = (ppar - pperp) * bq[i] * bq[j] + (i == j ? pperp : 0.0) + sq[i] * bq[j] + bq[i] * sq[j];
    }
    q[i] = 0.5 * bq[i] * h + rq[i];
  }
}

static void
test_vol(int cdim)
{
  struct field_ctx ctx = {.cdim = cdim, .uniform_rho = false, .b_along_x = false};
  struct test_setup s;
  setup(&s, &ctx);
  int nc = s.basis.num_basis;
  struct gkyl_dg_eqn *eqn = make_eqn(&s, NULL);

  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *M = gkyl_array_cfetch(s.fluid, cidx), *u = gkyl_array_cfetch(s.u, cidx);
    double xc[3], out[64] = {0.0}, ref[64] = {0.0};
    gkyl_rect_grid_cell_center(&s.grid, iter.idx, xc);
    eqn->vol_term(eqn, xc, s.grid.dx, iter.idx, M, out);

    int nnodes = 1;
    for (int d = 0; d < cdim; ++d) {
      nnodes *= 4;
    }
    for (int n = 0; n < nnodes; ++n) {
      double z[3], w = 1.0;
      int m = n;
      for (int d = cdim - 1; d >= 0; --d) {
        z[d] = gl4_x[m % 4];
        w *= gl4_w[m % 4];
        m /= 4;
      }
      double P[3][3], q[3], rho, E;
      pressure_at(&s, cidx, z, P, q, &rho, &E);
      double Mq[3], uq[3];
      for (int i = 0; i < 3; ++i) {
        Mq[i] = s.basis.eval_expand(z, &M[(1 + i) * nc]);
        uq[i] = s.basis.eval_expand(z, &u[i * nc]);
      }
      for (int k = 0; k < nc; ++k) {
        double unit[8] = {0.0};
        unit[k] = 1.0;
        for (int d = 0; d < cdim; ++d) {
          double dphi = 2.0 / s.grid.dx[d] * s.basis.eval_grad_expand(d, z, unit);
          double Pu = 0.0;
          for (int i = 0; i < 3; ++i) {
            ref[(1 + i) * nc + k] += w * dphi * (Mq[i] * uq[d] + P[i][d]);
            Pu += P[i][d] * uq[i];
          }
          ref[0 * nc + k] += w * dphi * Mq[d];
          ref[4 * nc + k] += w * dphi * (uq[d] * E + Pu + q[d]);
        }
      }
    }
    for (int k = 0; k < 5 * nc; ++k) {
      double err = fabs(out[k] - ref[k]);
      max_err = err > max_err ? err : max_err;
      max_ref = fabs(ref[k]) > max_ref ? fabs(ref[k]) : max_ref;
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_ref);
  TEST_MSG("cdim=%d: max volume error %g (max reference %g)", cdim, max_err, max_ref);

  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}

// Uniform state of the face-flux tests.
struct absorb_ctx {
  int cdim;
  double rho, xi, u[3], b[3];
};

static void
eval_absorb_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  double uu = c->u[0] * c->u[0] + c->u[1] * c->u[1] + c->u[2] * c->u[2];
  fout[0] = c->rho;
  for (int i = 0; i < 3; ++i) {
    fout[1 + i] = c->rho * c->u[i];
  }
  fout[4] = 0.5 * c->rho * uu + 1.5 * c->rho * c->xi * c->xi;
}

static void
eval_absorb_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  for (int i = 0; i < 3; ++i) {
    fout[i] = c->u[i];
  }
}

static void
eval_absorb_b(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  for (int i = 0; i < 3; ++i) {
    fout[i] = c->b[i];
  }
}

static void
eval_absorb_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  fout[0] = c->xi;
  fout[1] = c->xi * c->xi;
  fout[2] = c->rho * c->xi * c->xi;
}

static void
eval_absorb_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  double w = xn[c->cdim];
  fout[0] = fout[1] = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
}

static void
eval_one(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

// One side of a face for the reference flux: the uniform traces of a cell and whether it holds the
// Maxwellian distribution F = G = M(w), M0_perp = N_perp = 0 (or is an empty absorbing ghost).
struct ref_side {
  double rho, xi, theta, u[3], b[3];
  bool full;
};

// Reference fluid flux through a face with normal e_d between two uniform sides (notes, eqs.
// kinetic-face-flux, vector-HLL and fluid-flux-map), written out independently of the kernels:
// per w the side contributions of the eight-block HLL Phi^l = (S_R J_l - S_L S_R X_l)/ds,
// Phi^r = (-S_L J_r + S_L S_R X_r)/ds of the flux function with face-averaged coefficients
// (plan 17.41), J(X) = a X + T(X), a = {u_d} + {xi} b_d w, T with {xi}, b = {b}/|{b}|, {theta}, applied to
// each side's state, the bounds S_L = min(a - c, 0), S_R = max(a + c, 0) with the averaged
// c = xi |b x e_d| sqrt((2 + sqrt 2) theta) (plan 17.42), then the moment map of each contribution
// with its own frame, integrated over w in [-wmax, wmax] by Simpson's rule.
static void
ref_kinetic_flux(int d, const struct ref_side *sl, const struct ref_side *sr, double wmax, double *F)
{
  const struct ref_side *sd[2] = {sl, sr};
  double mom[2][7] = {{0.0}}; // per side: int Phi_F, int w Phi_F, int w^2 Phi_F, int Phi_G, int Phi_M0
  int n = 48000;
  double hw = 2.0 * wmax / n;
  for (int i = 0; i <= n; ++i) {
    double w = -wmax + i * hw, wt = (i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    double maxw = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
    // face averages of the coefficients (index 2) next to those of the two sides
    struct ref_side avg = {.xi = 0.5 * (sl->xi + sr->xi), .theta = 0.5 * (sl->theta + sr->theta)};
    for (int a = 0; a < 3; ++a) {
      avg.u[a] = 0.5 * (sl->u[a] + sr->u[a]);
      avg.b[a] = 0.5 * (sl->b[a] + sr->b[a]);
    }
    double bmag = sqrt(avg.b[0] * avg.b[0] + avg.b[1] * avg.b[1] + avg.b[2] * avg.b[2]);
    for (int a = 0; a < 3; ++a) {
      avg.b[a] = bmag > 0.0 ? avg.b[a] / bmag : 0.0;
    }
    double X[2][8], J[2][8];
    double V = avg.u[d] + avg.xi * w * avg.b[d];
    double bx[3], e[3] = {0.0, 0.0, 0.0};
    e[d] = 1.0;
    for (int a = 0; a < 3; ++a) {
      bx[a] = avg.b[(a + 1) % 3] * e[(a + 2) % 3] - avg.b[(a + 2) % 3] * e[(a + 1) % 3];
    }
    double c = avg.xi * sqrt(bx[0] * bx[0] + bx[1] * bx[1] + bx[2] * bx[2]) * sqrt((2.0 + sqrt(2.0)) * avg.theta);
    double s_l = fmin(V - c, 0.0), s_r = fmax(V + c, 0.0);
    for (int k = 0; k < 2; ++k) {
      for (int m = 0; m < 8; ++m) {
        X[k][m] = 0.0;
      }
      if (sd[k]->full) {
        X[k][0] = X[k][1] = maxw;
      }
      J[k][0] = V * X[k][0] + avg.xi * X[k][2 + d];
      J[k][1] = V * X[k][1] + avg.xi * X[k][5 + d];
      for (int a = 0; a < 3; ++a) {
        double pn = (a == d ? 1.0 : 0.0) - avg.b[a] * avg.b[d];
        J[k][2 + a] = V * X[k][2 + a] + avg.xi * X[k][1] * pn;
        J[k][5 + a] = V * X[k][5 + a] +
          avg.xi * (4.0 * avg.theta * X[k][1] - 2.0 * avg.theta * avg.theta * X[k][0]) * pn;
      }
    }
    double ds = s_r - s_l;
    if (!(ds > 0.0)) {
      continue;
    }
    for (int k = 0; k < 2; ++k) {
      double phi[8];
      for (int m = 0; m < 8; ++m) {
        phi[m] = k == 0 ? (s_r * J[0][m] - s_l * s_r * X[0][m]) / ds : (-s_l * J[1][m] + s_l * s_r * X[1][m]) / ds;
      }
      mom[k][0] += wt * phi[0];
      mom[k][1] += wt * w * phi[0];
      mom[k][2] += wt * w * w * phi[0];
      mom[k][3] += wt * phi[1];
      for (int a = 0; a < 3; ++a) {
        mom[k][4 + a] += wt * phi[2 + a];
      }
    }
  }
  for (int m = 0; m < 5; ++m) {
    F[m] = 0.0;
  }
  for (int k = 0; k < 2; ++k) {
    const struct ref_side *q = sd[k];
    double I[7];
    for (int m = 0; m < 7; ++m) {
      I[m] = mom[k][m] * hw / 3.0;
    }
    double bu = 0.0, uu = 0.0, uIM = 0.0;
    for (int a = 0; a < 3; ++a) {
      bu += q->b[a] * q->u[a];
      uu += q->u[a] * q->u[a];
      uIM += q->u[a] * I[4 + a];
    }
    F[0] += q->rho * I[0];
    for (int a = 0; a < 3; ++a) {
      F[1 + a] += q->rho * (q->u[a] * I[0] + q->xi * q->b[a] * I[1] + q->xi * I[4 + a]);
    }
    F[4] += q->rho * (0.5 * uu * I[0] + q->xi * bu * I[1] + q->xi * uIM + q->xi * q->xi * (0.5 * I[2] + I[3]));
  }
}

// Fluid fluxes as moments of the kinetic flux of the eight blocks for a uniform isotropic state
// (Maxwellian F = G, theta = 1, no first harmonic): (a) interior faces, where the HLL flux reduces
// to the exact flux function, so that the fluid flux is the ideal one (isotropic pressure, its
// work), and (b) an absorbing wall (ghost F = 0, theta copied), both against the independent
// reference above (the two sides share one frame there, so pairing them at equal w is exact);
// (c) a reflecting wall (mirrored ghost), where the mass, energy and tangential momentum fluxes
// cancel exactly and the normal momentum flux has a closed form. The flow and the field angle put the kinks of the HLL
// weights (V = 0, V = +-c) on velocity cell edges, so the kernels' Gauss-Legendre integration of the
// projected Maxwellian is exact and the reference differs from them by the projection error only.
static void
test_kinetic_flux(int cdim)
{
  double wmax = 8.0;
  int nv = 96; // dv = 1/6
  struct absorb_ctx ctx = {.cdim = cdim, .rho = 2.0, .xi = 1.5};
  // c/(xi b_x) = tan(angle) sqrt(2 + sqrt 2) = 2/3 in 2x (a whole number of velocity cells)
  double angle = cdim > 1 ? atan(2.0 / (3.0 * sqrt(2.0 + sqrt(2.0)))) : 0.0;
  ctx.b[0] = cos(angle);
  ctx.b[1] = sin(angle);
  ctx.b[2] = 0.0;
  // V = 0 at w0 = 1/3 (a velocity cell edge): u_n = -xi b_n w0
  ctx.u[0] = -ctx.xi * ctx.b[0] / 3.0;
  ctx.u[1] = 0.1;
  ctx.u[2] = -0.2;

  double lower[4] = {-1.0, -1.0, -1.0, -wmax}, upper[4] = {1.0, 1.0, 1.0, wmax};
  int cells[4] = {4, 3, 2, nv}, ghost[4] = {1, 1, 1, 0};
  lower[cdim] = -wmax;
  upper[cdim] = wmax;
  cells[cdim] = nv;
  ghost[cdim] = 0;
  struct gkyl_rect_grid cgrid, pgrid;
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&pgrid, cdim + 1, lower, upper, cells);
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext;
  gkyl_create_grid_ranges(&cgrid, ghost, &clocal_ext, &clocal);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&cgrid, &clocal, &clocal_ext, &cbasis);
  gkyl_create_grid_ranges(&pgrid, ghost, &plocal_ext, &plocal);
  int nc = cbasis.num_basis;

  struct gkyl_array *fluid = mkarr(5 * nc, clocal_ext.volume), *u = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_array *bvar = mkarr(3 * nc, clocal_ext.volume), *xi = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_array *theta = mkarr(nc, clocal_ext.volume);
  struct gkyl_array *f = mkarr(8 * pbasis.num_basis, plocal_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 5, eval_absorb_fluid, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, fluid);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 3, eval_absorb_u, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 3, eval_absorb_b, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, bvar);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 3, eval_absorb_xi, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 1, eval_one, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, theta);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_clear(f, 0.0);
  proj = gkyl_proj_on_basis_new(&pgrid, &pbasis, 8, 2, eval_absorb_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &plocal_ext, f);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_dg_pkpm_remap_flux_surf_inp inp = {
    .phase_grid = &pgrid,
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .pos_map = pm,
    .num_face_nodes = 3,
    .use_gpu = false,
  };
  struct gkyl_dg_pkpm_remap_flux_surf *up = gkyl_dg_pkpm_remap_flux_surf_inew(&inp);
  int ns = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(up);
  TEST_CHECK(ns == nc / 2);
  struct gkyl_array *flux = mkarr(5 * cdim * ns, clocal_ext.volume);
  // the kinetic flux, the kinetic rate and the blend weight are computed alongside (not checked)
  int nk = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(up);
  struct gkyl_array *kin_flux = mkarr(2 * 8 * cdim * nk, plocal_ext.volume);
  struct gkyl_array *kin_cfl = mkarr(1, plocal_ext.volume);
  struct gkyl_array *blend = mkarr(cdim, clocal_ext.volume);
  struct gkyl_array *cfl = mkarr(1, clocal_ext.volume);
  gkyl_array_clear(kin_flux, 0.0);
  gkyl_array_clear(kin_cfl, 0.0);
  gkyl_array_clear(blend, 0.0);
  gkyl_array_clear(cfl, 0.0);
  // modal coefficient of a constant on the surface basis: sqrt(2)^(cdim-1)
  double c0 = pow(sqrt(2.0), cdim - 1);
  double uu = 0.0;
  for (int i = 0; i < 3; ++i) {
    uu += ctx.u[i] * ctx.u[i];
  }
  double scale = ctx.rho * (uu + ctx.xi * ctx.xi) * (fabs(ctx.u[0]) + ctx.xi);
  struct ref_side in = {.rho = ctx.rho, .xi = ctx.xi, .theta = 1.0, .full = true};
  for (int i = 0; i < 3; ++i) {
    in.u[i] = ctx.u[i];
    in.b[i] = ctx.b[i];
  }

  // (a) uniform interior faces in every direction: the exact flux, P = rho xi^2 I
  gkyl_array_clear(flux, 0.0);
  gkyl_dg_pkpm_remap_flux_surf_advance(
    up, &clocal, &clocal_ext, &plocal, &plocal_ext, bvar, fluid, u, xi, theta, f, kin_flux, flux,
    kin_cfl, cfl, blend
  );
  {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &clocal);
    while (gkyl_range_iter_next(&iter)) {
      const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&clocal_ext, iter.idx));
      for (int dir = 0; dir < cdim; ++dir) {
        double un = ctx.u[dir], ref[5], refk[5];
        ref[0] = ctx.rho * un;
        for (int i = 0; i < 3; ++i) {
          ref[1 + i] = ctx.rho * ctx.u[i] * un + (i == dir ? ctx.rho * ctx.xi * ctx.xi : 0.0);
        }
        ref[4] = un * (0.5 * ctx.rho * uu + 2.5 * ctx.rho * ctx.xi * ctx.xi);
        ref_kinetic_flux(dir, &in, &in, wmax, refk);
        for (int c = 0; c < 5; ++c) {
          double got = fl[dir * 5 * ns + c * ns] / c0;
          TEST_CHECK(gkyl_compare(got, ref[c], 2e-4 * scale) && gkyl_compare(refk[c], ref[c], 2e-4 * scale));
          TEST_MSG("cdim=%d dir=%d comp=%d: got %g expected %g (reference %g)", cdim, dir, c, got, ref[c], refk[c]);
          for (int k = 1; k < ns; ++k) {
            TEST_CHECK(fabs(fl[dir * 5 * ns + c * ns + k]) < 1e-12 * scale);
          }
        }
      }
    }
  }

  // (b) absorbing wall in x: zero distribution in the x ghost cells; fluid, u, b, xi and theta
  // continue the uniform interior state (as the app's mirrored fluid and theta do)
  int dir = 0;
  for (int e = 0; e < 2; ++e) {
    struct gkyl_range skin, ghost_r;
    gkyl_skin_ghost_ranges(&skin, &ghost_r, dir, e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE, &plocal_ext, ghost);
    gkyl_array_clear_range(f, 0.0, &ghost_r);
  }
  gkyl_array_clear(flux, 0.0);
  gkyl_dg_pkpm_remap_flux_surf_advance(
    up, &clocal, &clocal_ext, &plocal, &plocal_ext, bvar, fluid, u, xi, theta, f, kin_flux, flux,
    kin_cfl, cfl, blend
  );
  struct ref_side empty = in;
  empty.full = false;
  for (int e = 0; e < 2; ++e) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct gkyl_range skin, ghost_r;
    gkyl_skin_ghost_ranges(&skin, &ghost_r, dir, edge, &clocal_ext, ghost);
    double ref[5];
    if (e == 0) {
      ref_kinetic_flux(dir, &empty, &in, wmax, ref);
    } else {
      ref_kinetic_flux(dir, &in, &empty, wmax, ref);
    }
    // the flux of the lower wall is stored in the skin cell, that of the upper wall in
    // the ghost cell beyond it
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &skin);
    while (gkyl_range_iter_next(&iter)) {
      int idx[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx);
      if (e == 1) {
        idx[dir] += 1;
      }
      const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&clocal_ext, idx));
      for (int c = 0; c < 5; ++c) {
        double got = fl[dir * 5 * ns + c * ns] / c0;
        TEST_CHECK(gkyl_compare(got, ref[c], 1e-4 * scale));
        TEST_MSG("cdim=%d absorb edge=%d comp=%d: got %g expected %g", cdim, e, c, got, ref[c]);
        for (int k = 1; k < ns; ++k) {
          TEST_CHECK(fabs(fl[dir * 5 * ns + c * ns + k]) < 1e-12 * scale);
        }
      }
    }
    TEST_CHECK(e == 0 ? ref[0] < 0.0 : ref[0] > 0.0);
  }

  // (c) reflecting wall in x: mirrored ghost (F(-w), -u_n, -M_n; the Maxwellian is even). In 1x the
  // field along the normal is line-tied (an external field carried into the ghost unchanged); in 2x
  // and 3x the oblique field has the PEC image (-b_n) of a fluctuating normal field.
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 1, eval_one, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, theta);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&pgrid, &pbasis, 8, 2, eval_absorb_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &plocal_ext, f);
  gkyl_proj_on_basis_release(proj);
  for (int e = 0; e < 2; ++e) {
    struct gkyl_range skin, ghost_r;
    gkyl_skin_ghost_ranges(&skin, &ghost_r, dir, e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE, &clocal_ext, ghost);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost_r);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&clocal_ext, iter.idx);
      double *ug = gkyl_array_fetch(u, loc), *bg = gkyl_array_fetch(bvar, loc), *flg = gkyl_array_fetch(fluid, loc);
      for (int k = 0; k < nc; ++k) {
        ug[dir * nc + k] *= -1.0;
        flg[(1 + dir) * nc + k] *= -1.0;
        if (cdim > 1) {
          bg[dir * nc + k] *= -1.0;
        }
      }
    }
  }
  gkyl_array_clear(flux, 0.0);
  gkyl_dg_pkpm_remap_flux_surf_advance(
    up, &clocal, &clocal_ext, &plocal, &plocal_ext, bvar, fluid, u, xi, theta, f, kin_flux, flux,
    kin_cfl, cfl, blend
  );
  // Normal momentum flux of the wall. The face frame is at rest along the normal, and the two
  // sides are paired at equal velocity in it.
  // 1x (b along the normal): the HLL weights reduce to upwinding, so the flux is that of the gas
  // and its image, twice the momentum flux of the particles moving towards the wall.
  // 2x, 3x (the face-averaged b is tangential): nothing streams through the face; the drifts
  // -/+ u_n across b enter the face frame as first harmonics -/+ (u_n/xi) F and as the
  // perpendicular energy u_n^2/2, so the flux is the pressure rho (xi^2 + u_n^2/2) and the HLL
  // dissipation rho c u_n of the jump, c = xi sqrt(2 + sqrt 2).
  double un = ctx.u[dir], ref_mn[2];
  if (cdim == 1) {
    double toward_lower = (un * un + ctx.xi * ctx.xi) * 0.5 * erfc(un / (sqrt(2.0) * ctx.xi)) -
      ctx.xi * un * exp(-0.5 * un * un / (ctx.xi * ctx.xi)) / sqrt(2.0 * M_PI);
    ref_mn[0] = 2.0 * ctx.rho * toward_lower;
    ref_mn[1] = 2.0 * ctx.rho * (un * un + ctx.xi * ctx.xi - toward_lower);
  } else {
    double pres = ctx.rho * (ctx.xi * ctx.xi + 0.5 * un * un);
    double diss = ctx.rho * ctx.xi * sqrt(2.0 + sqrt(2.0)) * un;
    ref_mn[0] = pres - diss;
    ref_mn[1] = pres + diss;
  }
  for (int e = 0; e < 2; ++e) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct gkyl_range skin, ghost_r;
    gkyl_skin_ghost_ranges(&skin, &ghost_r, dir, edge, &clocal_ext, ghost);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &skin);
    while (gkyl_range_iter_next(&iter)) {
      int idx[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx);
      if (e == 1) {
        idx[dir] += 1;
      }
      const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&clocal_ext, idx));
      double got_rho = fl[dir * 5 * ns] / c0, got_mn = fl[dir * 5 * ns + (1 + dir) * ns] / c0, got_e = fl[dir * 5 * ns + 4 * ns] / c0;
      TEST_CHECK(fabs(got_rho) < 1e-12 * scale);
      TEST_CHECK(fabs(got_e) < 1e-12 * scale);
      TEST_MSG("cdim=%d reflect edge=%d: mass flux %g energy flux %g", cdim, e, got_rho, got_e);
      TEST_CHECK(gkyl_compare(got_mn, ref_mn[e], 1e-4 * scale));
      TEST_MSG("cdim=%d reflect edge=%d: normal momentum flux %g expected %g", cdim, e, got_mn, ref_mn[e]);
      for (int i = 0; i < 3; ++i) {
        if (i != dir) {
          TEST_CHECK(fabs(fl[dir * 5 * ns + (1 + i) * ns]) < 1e-12 * scale);
          TEST_MSG("cdim=%d reflect edge=%d: tangential momentum flux %d: %g", cdim, e, i, fl[dir * 5 * ns + (1 + i) * ns]);
        }
      }
    }
  }

  gkyl_array_release(cfl);
  gkyl_array_release(flux);
  gkyl_array_release(kin_flux);
  gkyl_array_release(kin_cfl);
  gkyl_array_release(blend);
  gkyl_dg_pkpm_remap_flux_surf_release(up);
  gkyl_vlasov_position_map_release(pm);
  gkyl_array_release(f);
  gkyl_array_release(theta);
  gkyl_array_release(xi);
  gkyl_array_release(bvar);
  gkyl_array_release(u);
  gkyl_array_release(fluid);
}

void test_kinetic_flux_1x() { test_kinetic_flux(1); }
void test_kinetic_flux_2x() { test_kinetic_flux(2); }
void test_vol_1x() { test_vol(1); }
void test_vol_2x() { test_vol(2); }
void test_vol_3x() { test_vol(3); }

TEST_LIST = {
  {"test_vol_1x", test_vol_1x},
  {"test_vol_2x", test_vol_2x},
  {"test_vol_3x", test_vol_3x},
  {"test_kinetic_flux_1x", test_kinetic_flux_1x},
  {"test_kinetic_flux_2x", test_kinetic_flux_2x},
  {NULL, NULL},
};
