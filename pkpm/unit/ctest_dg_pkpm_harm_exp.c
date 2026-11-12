// Tests of the gyrophase-mixing integrating factor of the first-harmonic blocks (dg_pkpm_harm_exp),
// driven through the three stages exactly as the SSP-RK3 driver of the app does, with the gyration
// explicit in the stage right-hand side: the nodal transforms are exact inverses, a step with the
// mixing active is e^(-nu h) R_RK3(h Omega J) on the plane perpendicular to b (exact relaxation, the
// SSP-RK3 gyration), states inside the bound are SSP-RK3 bit for bit, and the mixing relaxes a column
// to the bound with int M0_perp dw kept zero.
#include <math.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_stage_weights.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_exp_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// One phase cell of 8 blocks (the first-harmonic blocks are 2-7) and its configuration cell.
struct hx_cell {
  int cdim, np, nc;
  struct gkyl_basis cbasis;
  struct gkyl_rect_grid pgrid;
  struct gkyl_range crange, prange;
  struct gkyl_dg_pkpm_harm_exp *up;
  struct gkyl_array *em, *coef, *nu, *act;
  double om[3]; // (q/m) B (qbym = 1), uniform
};

static struct hx_cell
hx_cell_new(int cdim, double mix, double om0)
{
  struct hx_cell c = {.cdim = cdim};
  gkyl_cart_modal_serendip(&c.cbasis, cdim, 1);
  c.nc = c.cbasis.num_basis;
  c.np = 3 * c.nc;
  int lo[GKYL_MAX_DIM], hi[GKYL_MAX_DIM];
  for (int d = 0; d < cdim + 1; ++d) {
    lo[d] = hi[d] = 1;
  }
  gkyl_range_init(&c.crange, cdim, lo, hi);
  gkyl_range_init(&c.prange, cdim + 1, lo, hi);
  double plo[GKYL_MAX_DIM] = {0.0}, pup[GKYL_MAX_DIM] = {1.0, 1.0, 1.0, 1.0};
  int pcells[GKYL_MAX_DIM] = {1, 1, 1, 1};
  gkyl_rect_grid_init(&c.pgrid, cdim + 1, plo, pup, pcells);
  c.up = gkyl_dg_pkpm_harm_exp_new(&(struct gkyl_dg_pkpm_harm_exp_inp) {
    .phase_grid = &c.pgrid, .conf_basis = &c.cbasis, .qbym = 1.0, .mix = mix, .use_gpu = false});
  c.em = gkyl_array_new(GKYL_DOUBLE, 8 * c.nc, 1);
  c.coef = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_coef(c.up), 1);
  c.nu = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_rate(c.up), 1);
  c.act = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  double n = sqrt(0.29);
  c.om[0] = om0 * 0.3 / n;
  c.om[1] = om0 * -0.2 / n;
  c.om[2] = om0 * 0.4 / n;
  // a spatially uniform field (only the constant mode of the configuration basis)
  double *em = gkyl_array_fetch(c.em, 0);
  for (int k = 0; k < 8 * c.nc; ++k) {
    em[k] = 0.0;
  }
  for (int a = 0; a < 3; ++a) {
    em[(3 + a) * c.nc] = c.om[a] * pow(sqrt(2.0), cdim);
  }
  return c;
}

static void
hx_cell_release(struct hx_cell *c)
{
  gkyl_dg_pkpm_harm_exp_release(c->up);
  gkyl_array_release(c->em);
  gkyl_array_release(c->coef);
  gkyl_array_release(c->nu);
  gkyl_array_release(c->act);
}

// rhs = U x omega (the gyration) in the first-harmonic blocks, F and G get zero. A uniform omega
// makes the Galerkin product with U the cross product of the coefficient vectors.
static void
hx_rhs(const struct hx_cell *c, const double *f, double *rhs)
{
  const double *om = c->om;
  for (int k = 0; k < 8 * c->np; ++k) {
    rhs[k] = 0.0;
  }
  for (int v = 0; v < 2; ++v) {
    for (int k = 0; k < c->np; ++k) {
      double u[3];
      for (int a = 0; a < 3; ++a) {
        u[a] = f[(2 + 3 * v + a) * c->np + k];
      }
      rhs[(2 + 3 * v) * c->np + k] = u[1] * om[2] - u[2] * om[1];
      rhs[(3 + 3 * v) * c->np + k] = u[2] * om[0] - u[0] * om[2];
      rhs[(4 + 3 * v) * c->np + k] = u[0] * om[1] - u[1] * om[0];
    }
  }
}

// One step of the app's SSP-RK3 driver (Shu-Osher form, forward Euler stages) with the mixing
// integrating factor written into the first-harmonic blocks of each stage output.
static void
hx_step(struct hx_cell *c, double h, bool use_mix, double *f)
{
  int n = 8 * c->np;
  double fin[8 * 24], rhs[8 * 24];
  struct gkyl_array *fa = gkyl_array_new(GKYL_DOUBLE, n, 1), *oa = gkyl_array_new(GKYL_DOUBLE, n, 1);
  double *f0 = gkyl_array_fetch(fa, 0), *out = gkyl_array_fetch(oa, 0);
  gkyl_dg_pkpm_harm_exp_coef(c->up, &c->crange, c->em, c->coef);
  gkyl_array_clear(c->act, 0.0);
  for (int k = 0; k < n; ++k) {
    f0[k] = fin[k] = f[k];
  }
  for (int s = 1; s <= 3; ++s) {
    hx_rhs(c, fin, rhs);
    for (int k = 0; k < n; ++k) {
      out[k] = fin[k] + h * rhs[k];
    }
    if (use_mix) {
      gkyl_dg_pkpm_harm_exp_advance(c->up, s, h, &c->crange, &c->prange, c->coef, fa, 0, c->nu, c->act, oa);
    }
    for (int k = 0; k < n; ++k) {
      fin[k] = s == 1 ? out[k] : (s == 2 ? 0.75 * f[k] + 0.25 * out[k] : f[k] / 3.0 + 2.0 / 3.0 * out[k]);
    }
  }
  for (int k = 0; k < n; ++k) {
    f[k] = fin[k];
  }
  gkyl_array_release(fa);
  gkyl_array_release(oa);
}

// A state from nodal first-harmonic vectors perpendicular to omega with |M0_perp| = amp (1 + 0.1 p)
// at phase node p (r_M = sqrt2 |M0_perp| with F = G = 1) and |N_perp| half of that (r_N < r_M), and
// F = G = 1 at the nodes.
static void
hx_init_nodal(const struct hx_cell *c, double amp, double *f)
{
  double mag = sqrt(c->om[0] * c->om[0] + c->om[1] * c->om[1] + c->om[2] * c->om[2]);
  double oh[3] = {c->om[0] / mag, c->om[1] / mag, c->om[2] / mag};
  // e1 perpendicular to oh, e2 = oh x e1
  double e1[3] = {oh[1], -oh[0], 0.0}, n1 = sqrt(oh[0] * oh[0] + oh[1] * oh[1]);
  for (int a = 0; a < 3; ++a) {
    e1[a] /= n1;
  }
  double e2[3] = {oh[1] * e1[2] - oh[2] * e1[1], oh[2] * e1[0] - oh[0] * e1[2], oh[0] * e1[1] - oh[1] * e1[0]};
  for (int k = 0; k < 8 * c->np; ++k) {
    f[k] = 0.0;
  }
  // F = G = 1: the constant mode of the orthonormal basis is 2^(-(cdim+1)/2)
  f[0] = f[c->np] = pow(sqrt(2.0), c->cdim + 1);
  double Xn[24 * 6];
  for (int p = 0; p < c->np; ++p) {
    double am = amp * (1.0 + 0.1 * p) / sqrt(2.0);
    for (int v = 0; v < 2; ++v) {
      double ph = 0.7 * p + 1.3 * v, s = v == 0 ? am : 0.5 * am;
      for (int a = 0; a < 3; ++a) {
        Xn[6 * p + 3 * v + a] = s * (cos(ph) * e1[a] + sin(ph) * e2[a]);
      }
    }
  }
  pkpm_hx_to_modal(c->cdim, Xn, f);
}

static void
test_transforms(int cdim)
{
  struct hx_cell c = hx_cell_new(cdim, 0.0, 1.0);
  double f[8 * 24], g[8 * 24], fn[24 * 6];
  for (int k = 0; k < 8 * c.np; ++k) {
    f[k] = sin(0.7 * k + 0.1);
    g[k] = -1.0;
  }
  pkpm_hx_to_nodes(cdim, f, fn);
  pkpm_hx_to_modal(cdim, fn, g);
  for (int k = 0; k < 2 * c.np; ++k) {
    TEST_CHECK(g[k] == -1.0); // F and G untouched
  }
  double err = 0.0;
  for (int k = 2 * c.np; k < 8 * c.np; ++k) {
    err = fmax(err, fabs(g[k] - f[k]));
  }
  TEST_CHECK(err < 1e-14);
  TEST_MSG("cdim %d: round trip error %g", cdim, err);
  hx_cell_release(&c);
}

static void test_transforms_1x(void) { test_transforms(1); }
static void test_transforms_2x(void) { test_transforms(2); }
static void test_transforms_3x(void) { test_transforms(3); }

// The mixing active at every node (r = 1.5 (1 + 0.1 p) at phase node p), Omega h = 0.5: one step is,
// node by node, e^(-nu h) R_RK3(h Omega J) on the plane perpendicular to b (the damping and the
// rotation commute there, so the integrating factor is exact), with nu the rate of the start of the
// step; uncapped (mix 2) and at the cap ln(r)/h (mix 100, the step ends on the bound up to R_RK3).
static void
test_damp_step(int cdim)
{
  double h = 0.1, mixes[2] = {2.0, 100.0};
  for (int im = 0; im < 2; ++im) {
    struct hx_cell c = hx_cell_new(cdim, mixes[im], 5.0);
    double f[8 * 24], f0[8 * 24];
    hx_init_nodal(&c, 1.5, f);
    for (int k = 0; k < 8 * c.np; ++k) {
      f0[k] = f[k];
    }
    hx_step(&c, h, true, f);
    double Mn[24 * 6], Xn[24 * 6], fgn[24 * 2];
    pkpm_hx_to_nodes(cdim, f0, Mn);
    pkpm_hx_fg_to_nodes(cdim, f0, fgn);
    pkpm_hx_to_nodes(cdim, f, Xn);
    double oh[3];
    double mag = pkpm_hx_axis(c.om, oh), y = mag * h, err = 0.0, vmax = 0.0;
    int ncap = 0;
    for (int p = 0; p < c.np; ++p) {
      double nu = pkpm_hx_rate(c.up->mix, mag, h, fgn[2 * p], fgn[2 * p + 1], &Mn[6 * p], &Mn[6 * p + 3]);
      TEST_CHECK(nu > 0.0);
      const double *m0 = &Mn[6 * p];
      double r = sqrt(2.0) * sqrt(m0[0] * m0[0] + m0[1] * m0[1] + m0[2] * m0[2]) / sqrt(fgn[2 * p] * fgn[2 * p + 1]);
      ncap += fabs(nu - log(r) / h) < 1e-12 * nu;
      double d = exp(-nu * h);
      for (int v = 0; v < 2; ++v) {
        const double *m = &Mn[6 * p + 3 * v];
        double mxo[3] = {m[1] * oh[2] - m[2] * oh[1], m[2] * oh[0] - m[0] * oh[2], m[0] * oh[1] - m[1] * oh[0]};
        for (int a = 0; a < 3; ++a) {
          double ref = d * ((1.0 - 0.5 * y * y) * m[a] + (y - y * y * y / 6.0) * mxo[a]);
          err = fmax(err, fabs(Xn[6 * p + 3 * v + a] - ref));
          vmax = fmax(vmax, fabs(m[a]));
        }
      }
    }
    TEST_CHECK(err < 1e-13 * vmax);
    TEST_MSG("cdim %d, mix %g: max nodal error %g (max |M| %g)", cdim, mixes[im], err, vmax);
    TEST_CHECK(im == 0 ? ncap == 0 : ncap == c.np); // uncapped at mix 2, capped everywhere at 100
    bool fg_same = true;
    for (int k = 0; k < 2 * c.np; ++k) {
      fg_same = fg_same && f[k] == f0[k];
    }
    TEST_CHECK(fg_same);
    hx_cell_release(&c);
  }
}

static void test_damp_step_1x(void) { test_damp_step(1); }
static void test_damp_step_2x(void) { test_damp_step(2); }
static void test_damp_step_3x(void) { test_damp_step(3); }

// Inside the bound (r <= 0.74) with the mixing on, and with it off, the step is SSP-RK3 bit for bit.
static void
test_inactive(void)
{
  double mixes[2] = {10.0, -1.0};
  for (int im = 0; im < 2; ++im) {
    struct hx_cell c = hx_cell_new(2, mixes[im], 5.0);
    double f[8 * 24], g[8 * 24];
    hx_init_nodal(&c, 0.35, f);
    for (int k = 0; k < 8 * c.np; ++k) {
      g[k] = f[k];
    }
    for (int s = 0; s < 5; ++s) {
      hx_step(&c, 0.1, true, f);
      hx_step(&c, 0.1, false, g);
    }
    bool same = true;
    for (int k = 0; k < 8 * c.np; ++k) {
      same = same && f[k] == g[k];
    }
    TEST_CHECK(same);
    TEST_MSG("mix %g", mixes[im]);
    hx_cell_release(&c);
  }
}

// The mixing rate of one node: zero inside the bound, mix |Omega| (1 - 1/r) beyond it, the full rate
// where F or G <= 0, capped at ln(r)/h.
static void
test_rate(void)
{
  double m0[3] = {0.0, 0.0, 0.0}, n[3] = {0.0, 0.0, 0.0};
  // r_M = sqrt2 |M0|/sqrt(F G) = 3 with F = G = 1
  m0[0] = 3.0 / sqrt(2.0);
  TEST_CHECK(fabs(pkpm_hx_rate(1.0, 2.0, 1e-3, 1.0, 1.0, m0, n) - 2.0 * (1.0 - 1.0 / 3.0)) < 1e-14);
  // capped: ln(3)/h < 2 (1 - 1/3) for h = 1
  TEST_CHECK(fabs(pkpm_hx_rate(1.0, 2.0, 1.0, 1.0, 1.0, m0, n) - log(3.0)) < 1e-14);
  // r_N = |N| sqrt(F)/G^(3/2) = 2 with F = 1, G = 4: |N| = 16
  m0[0] = 0.0;
  n[1] = 16.0;
  TEST_CHECK(fabs(pkpm_hx_rate(1.0, 1.0, 1e-3, 1.0, 4.0, m0, n) - 0.5) < 1e-14);
  // realizable: no damping
  n[1] = 7.0;
  TEST_CHECK(pkpm_hx_rate(1.0, 1.0, 1e-3, 1.0, 4.0, m0, n) == 0.0);
  // F or G <= 0: the full rate; no gyration: none
  TEST_CHECK(pkpm_hx_rate(2.0, 3.0, 1e-3, -1e-9, 4.0, m0, n) == 6.0);
  TEST_CHECK(pkpm_hx_rate(2.0, 0.0, 1e-3, -1e-9, 4.0, m0, n) == 0.0);
}

// A column of w cells (1x1v, one configuration cell) driven by the gyration alone, with the mixing.
struct hx_col {
  int nv, np;
  struct gkyl_basis cbasis;
  struct gkyl_rect_grid pgrid;
  struct gkyl_range crange, prange;
  struct gkyl_dg_pkpm_harm_exp *up;
  struct gkyl_array *em, *coef, *nu, *act, *col;
  struct gkyl_array *f, *f1, *fnew, *rhs;
};

static struct hx_col
hx_col_new(int nv, double mix)
{
  struct hx_col c = {.nv = nv};
  gkyl_cart_modal_serendip(&c.cbasis, 1, 1);
  c.np = 3 * c.cbasis.num_basis;
  int lo[2] = {1, 1}, chi[1] = {1}, phi[2] = {1, nv};
  gkyl_range_init(&c.crange, 1, lo, chi);
  gkyl_range_init(&c.prange, 2, lo, phi);
  double plo[2] = {0.0, -0.5 * nv}, pup[2] = {1.0, 0.5 * nv};
  int pcells[2] = {1, nv};
  gkyl_rect_grid_init(&c.pgrid, 2, plo, pup, pcells);
  c.up = gkyl_dg_pkpm_harm_exp_new(&(struct gkyl_dg_pkpm_harm_exp_inp) {
    .phase_grid = &c.pgrid, .conf_basis = &c.cbasis, .qbym = 1.0, .mix = mix, .use_gpu = false});
  c.em = gkyl_array_new(GKYL_DOUBLE, 8 * c.cbasis.num_basis, 1);
  c.coef = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_coef(c.up), 1);
  c.col = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_col(c.up), 1);
  c.act = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  c.nu = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_rate(c.up), nv);
  c.f = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, nv);
  c.f1 = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, nv);
  c.fnew = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, nv);
  c.rhs = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, nv);
  // B along z with (q/m)|B| = 1
  double *em = gkyl_array_fetch(c.em, 0);
  for (int k = 0; k < 8 * c.cbasis.num_basis; ++k) em[k] = 0.0;
  em[5 * c.cbasis.num_basis] = sqrt(2.0);
  return c;
}

static void
hx_col_release(struct hx_col *c)
{
  gkyl_dg_pkpm_harm_exp_release(c->up);
  struct gkyl_array *arrs[] = {c->em, c->coef, c->col, c->act, c->nu, c->f, c->f1, c->fnew, c->rhs};
  for (int i = 0; i < 9; ++i) gkyl_array_release(arrs[i]);
}

// rhs = U x (q/m) B (the gyration) in the first-harmonic blocks.
static void
hx_col_rhs(const struct hx_col *c, const struct gkyl_array *f, struct gkyl_array *rhs)
{
  for (int iv = 0; iv < c->nv; ++iv) {
    const double *u = gkyl_array_cfetch(f, iv);
    double *r = gkyl_array_fetch(rhs, iv);
    for (int k = 0; k < 8 * c->np; ++k) r[k] = 0.0;
    for (int v = 0; v < 2; ++v) {
      for (int k = 0; k < c->np; ++k) {
        // om = (0, 0, 1): U x om = (U_y, -U_x, 0)
        r[(2 + 3 * v) * c->np + k] = u[(3 + 3 * v) * c->np + k];
        r[(3 + 3 * v) * c->np + k] = -u[(2 + 3 * v) * c->np + k];
      }
    }
  }
}

// One forward Euler stage of the driver with the mixing integrating factor and the restoration.
static void
hx_col_stage(struct hx_col *c, int stage, double h, struct gkyl_array *fin, struct gkyl_array *fout)
{
  if (stage == 1) {
    gkyl_dg_pkpm_harm_exp_coef(c->up, &c->crange, c->em, c->coef);
    gkyl_array_clear(c->act, 0.0);
  }
  hx_col_rhs(c, fin, c->rhs);
  gkyl_array_set(fout, 1.0, fin);
  gkyl_array_accumulate(fout, h, c->rhs);
  gkyl_dg_pkpm_harm_exp_advance(c->up, stage, h, &c->crange, &c->prange, c->coef, c->f, 0, c->nu, c->act, fout);
  gkyl_dg_pkpm_harm_exp_restore(c->up, &c->crange, &c->prange, c->act, c->col, fout);
}

static void
hx_col_step(struct hx_col *c, double h)
{
  hx_col_stage(c, 1, h, c->f, c->f1);
  hx_col_stage(c, 2, h, c->f1, c->fnew);
  gkyl_array_set(c->f1, 0.75, c->f);
  gkyl_array_accumulate(c->f1, 0.25, c->fnew);
  hx_col_stage(c, 3, h, c->f1, c->fnew);
  gkyl_array_scale(c->f, 1.0 / 3.0);
  gkyl_array_accumulate(c->f, 2.0 / 3.0, c->fnew);
}

// F = G = a Gaussian of the cell centers (cell-constant), M0_perp along x: 3 F/sqrt2 in cell k0 minus
// the F-shaped compensation that makes int M0_perp dw = 0 (r = 2.61 in cell k0 of 8, k0 = 2).
static void
hx_col_init(struct hx_col *c, int k0)
{
  double c0 = 2.0; // the constant mode of the orthonormal 1x1v basis is 1/2
  double sF = 0.0, sM = 0.0;
  for (int iv = 0; iv < c->nv; ++iv) {
    double *f = gkyl_array_fetch(c->f, iv);
    for (int k = 0; k < 8 * c->np; ++k) f[k] = 0.0;
    double w = -0.5 * c->nv + iv + 0.5, F = exp(-0.5 * w * w);
    f[0] = c0 * F;
    f[c->np] = c0 * F;
    if (iv == k0) {
      f[2 * c->np] = c0 * 3.0 * F / sqrt(2.0);
    }
    sF += F;
    sM += iv == k0 ? 3.0 * F / sqrt(2.0) : 0.0;
  }
  for (int iv = 0; iv < c->nv; ++iv) {
    double *f = gkyl_array_fetch(c->f, iv);
    f[2 * c->np] -= sM / sF * f[0];
  }
}

// max r over the nodes, and int M0_x dw relative to its largest cell value
static void
hx_col_diag(const struct hx_col *c, double *rmax, double *mom)
{
  *rmax = 0.0;
  double s = 0.0, smax = 0.0;
  for (int iv = 0; iv < c->nv; ++iv) {
    const double *f = gkyl_array_cfetch(c->f, iv);
    double fn[24 * 6], fgn[24 * 2];
    pkpm_hx_to_nodes(1, f, fn);
    pkpm_hx_fg_to_nodes(1, f, fgn);
    for (int p = 0; p < 6; ++p) {
      double F = fgn[2 * p], G = fgn[2 * p + 1];
      double am = sqrt(fn[6 * p] * fn[6 * p] + fn[6 * p + 1] * fn[6 * p + 1] + fn[6 * p + 2] * fn[6 * p + 2]);
      *rmax = fmax(*rmax, sqrt(2.0) * am / sqrt(F * G));
    }
    s += f[2 * c->np];
    smax = fmax(smax, fabs(f[2 * c->np]));
  }
  *mom = fabs(s) / smax;
}

// Mixing at mix |Omega| h = 0.01: the excess decays to the bound (dr/dt = -(r - 1) unforced, here
// with the gyration and the restoration), int M0_perp dw stays zero, F and G are untouched.
static void
test_mix_relax(void)
{
  struct hx_col c = hx_col_new(8, 1.0);
  hx_col_init(&c, 2);
  struct gkyl_array *f0 = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, c.nv);
  gkyl_array_copy(f0, c.f);
  double rmax, mom, mom_max = 0.0;
  hx_col_diag(&c, &rmax, &mom);
  TEST_CHECK(rmax > 2.6 && rmax < 2.62);
  for (int s = 0; s < 1000; ++s) {
    hx_col_step(&c, 0.01);
    hx_col_diag(&c, &rmax, &mom);
    mom_max = fmax(mom_max, mom);
  }
  TEST_CHECK(rmax < 1.0 + 1e-3);
  TEST_MSG("max r after 10/Omega: %g", rmax);
  TEST_CHECK(mom_max < 1e-13);
  TEST_MSG("max |int M0 dw|/max|M0| over the run: %g", mom_max);
  bool fg_same = true;
  for (int iv = 0; iv < c.nv; ++iv) {
    const double *a = gkyl_array_cfetch(c.f, iv), *b = gkyl_array_cfetch(f0, iv);
    for (int k = 0; k < 2 * c.np; ++k) fg_same = fg_same && a[k] == b[k];
  }
  TEST_CHECK(fg_same);
  gkyl_array_release(f0);
  hx_col_release(&c);
}

// Stiff mixing (mix |Omega| h = 50): every step takes the excess to the bound (the cap ln(r)/h); what
// is left is the share of the restoration's frame shift that returns to the damped cell, which
// shrinks geometrically (2.61 -> 1.21 -> 1.027 -> 1.0035).
static void
test_mix_stiff(void)
{
  struct hx_col c = hx_col_new(8, 100.0);
  hx_col_init(&c, 2);
  double rmax, mom, r[3];
  for (int s = 0; s < 3; ++s) {
    hx_col_step(&c, 0.5);
    hx_col_diag(&c, &r[s], &mom);
    TEST_CHECK(mom < 1e-13);
    TEST_MSG("step %d: |int M0 dw|/max|M0|: %g", s + 1, mom);
  }
  TEST_CHECK(r[0] < 1.25 && r[1] < 1.05 && r[2] < 1.01);
  TEST_MSG("max r after 1, 2, 3 steps: %g %g %g", r[0], r[1], r[2]);
  hx_col_release(&c);
}

// Realizable everywhere: with mixing on the update is that without it, bit for bit.
static void
test_mix_inactive(void)
{
  struct hx_col a = hx_col_new(8, 1.0), b = hx_col_new(8, 0.0);
  hx_col_init(&a, 2);
  hx_col_init(&b, 2);
  for (int iv = 0; iv < 8; ++iv) {
    double *fa = gkyl_array_fetch(a.f, iv), *fb = gkyl_array_fetch(b.f, iv);
    for (int k = 2 * a.np; k < 8 * a.np; ++k) {
      fa[k] *= 0.2;
      fb[k] *= 0.2;
    }
  }
  for (int s = 0; s < 50; ++s) {
    hx_col_step(&a, 0.05);
    hx_col_step(&b, 0.05);
  }
  bool same = true;
  for (int iv = 0; iv < 8; ++iv) {
    const double *fa = gkyl_array_cfetch(a.f, iv), *fb = gkyl_array_cfetch(b.f, iv);
    for (int k = 0; k < 8 * a.np; ++k) same = same && fa[k] == fb[k];
  }
  TEST_CHECK(same);
  hx_col_release(&a);
  hx_col_release(&b);
}

// Stage weights: when the driver combines the state M at the start of the step and the stage value
// with w_M M + w_Y X, w_Y = c_Y lam / (c_M + c_Y lam) at each configuration node, the stage value
// built with lam makes that combination the integrating-factor stage with the same weights,
// w_M e^(.) M + w_Y e^(.) Y. e^(.) M and e^(.) Y are taken from two calls with the SSP-RK3 weights
// (Y = 0 gives c_M / c_Y (e^(.) M - M)). The combination is also formed with
// gkyl_pkpm_stage_weights, which must use the same configuration nodes.
static void
test_stage_weights(int cdim)
{
  double h = 0.1;
  struct hx_cell c = hx_cell_new(cdim, 2.0, 5.0);
  int n = 8 * c.np;
  struct gkyl_basis pbasis;
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  struct gkyl_pkpm_stage_weights *sw = gkyl_pkpm_stage_weights_new(&c.cbasis, &pbasis, false);
  struct gkyl_array *fa = gkyl_array_new(GKYL_DOUBLE, n, 1), *ya = gkyl_array_new(GKYL_DOUBLE, n, 1);
  struct gkyl_array *xa = gkyl_array_new(GKYL_DOUBLE, n, 1), *za = gkyl_array_new(GKYL_DOUBLE, n, 1);
  struct gkyl_array *la = gkyl_array_new(GKYL_DOUBLE, n, 1), *ca = gkyl_array_new(GKYL_DOUBLE, n, 1);
  struct gkyl_array *lam = gkyl_array_new(GKYL_DOUBLE, c.nc, 1), *lam_c = gkyl_array_new(GKYL_DOUBLE, c.nc, 1);
  double *m = gkyl_array_fetch(fa, 0), *y = gkyl_array_fetch(ya, 0), rhs[8 * 24];
  hx_init_nodal(&c, 1.5, m);
  hx_rhs(&c, m, rhs);
  // a temperature factor that differs at every configuration node (heating and cooling)
  double *lam_d = gkyl_array_fetch(lam, 0);
  for (int q = 0; q < c.nc; ++q) {
    lam_d[q] = q % 2 ? 0.4 + 0.1 * q : 3.0 + 2.0 * q;
  }
  // the rates of the step (stage 1)
  gkyl_dg_pkpm_harm_exp_coef(c.up, &c.crange, c.em, c.coef);
  gkyl_array_clear(c.act, 0.0);
  gkyl_array_copy(xa, fa);
  gkyl_dg_pkpm_harm_exp_advance(c.up, 1, h, &c.crange, &c.prange, c.coef, fa, 0, c.nu, c.act, xa);

  for (int stage = 2; stage <= 3; ++stage) {
    double coef_fe = stage == 2 ? 0.25 : 2.0 / 3.0, coef_old = 1.0 - coef_fe;
    for (int k = 0; k < n; ++k) {
      y[k] = m[k] + (0.3 + 0.2 * stage) * h * rhs[k];
    }
    // SSP-RK3 weights: X = c_M/c_Y (e M - M) + e Y, and with Y = 0 the first term alone
    gkyl_array_copy(xa, ya);
    gkyl_dg_pkpm_harm_exp_advance(c.up, stage, h, &c.crange, &c.prange, c.coef, fa, 0, c.nu, c.act, xa);
    gkyl_array_clear(za, 0.0);
    gkyl_dg_pkpm_harm_exp_advance(c.up, stage, h, &c.crange, &c.prange, c.coef, fa, 0, c.nu, c.act, za);
    // the stage value for the weighted combination
    gkyl_array_copy(la, ya);
    gkyl_dg_pkpm_harm_exp_advance(c.up, stage, h, &c.crange, &c.prange, c.coef, fa, lam, c.nu, c.act, la);
    // the weighted combination of the driver from the plain one
    gkyl_array_set(ca, coef_old, fa);
    gkyl_array_accumulate(ca, coef_fe, la);
    gkyl_array_copy(lam_c, lam);
    gkyl_pkpm_stage_weights_advance(sw, &c.crange, &c.prange, coef_fe, fa, la, lam_c, ca);

    double Mn[24 * 6], Xn[24 * 6], Zn[24 * 6], Ln[24 * 6], Cn[24 * 6], err = 0.0, err_comb = 0.0, vmax = 0.0;
    pkpm_hx_to_nodes(cdim, m, Mn);
    pkpm_hx_to_nodes(cdim, gkyl_array_cfetch(xa, 0), Xn);
    pkpm_hx_to_nodes(cdim, gkyl_array_cfetch(za, 0), Zn);
    pkpm_hx_to_nodes(cdim, gkyl_array_cfetch(la, 0), Ln);
    pkpm_hx_to_nodes(cdim, gkyl_array_cfetch(ca, 0), Cn);
    for (int p = 0; p < c.np; ++p) {
      double lam_q = lam_d[pkpm_hx_conf_nodal_idx(cdim, p / 3)];
      double w_fe = coef_fe * lam_q / (coef_old + coef_fe * lam_q), w_old = 1.0 - w_fe;
      for (int k = 0; k < 6; ++k) {
        double em = Mn[6 * p + k] + Zn[6 * p + k] * coef_fe / coef_old; // e M
        double ey = Xn[6 * p + k] - Zn[6 * p + k]; // e Y
        double ref = w_old * em + w_fe * ey;
        err = fmax(err, fabs(w_old * Mn[6 * p + k] + w_fe * Ln[6 * p + k] - ref));
        err_comb = fmax(err_comb, fabs(Cn[6 * p + k] - ref));
        vmax = fmax(vmax, fabs(Mn[6 * p + k]));
      }
    }
    TEST_CHECK(err < 1e-13 * vmax);
    TEST_MSG("cdim %d, stage %d: stage value error %g (max |M| %g)", cdim, stage, err, vmax);
    TEST_CHECK(err_comb < 1e-13 * vmax);
    TEST_MSG("cdim %d, stage %d: combination error %g (max |M| %g)", cdim, stage, err_comb, vmax);
    const double *lam_c_d = gkyl_array_cfetch(lam_c, 0);
    for (int q = 0; q < c.nc; ++q) {
      TEST_CHECK(gkyl_compare_double(lam_c_d[q], coef_old + coef_fe * lam_d[q], 1e-15));
    }
  }
  gkyl_pkpm_stage_weights_release(sw);
  gkyl_array_release(fa);
  gkyl_array_release(ya);
  gkyl_array_release(xa);
  gkyl_array_release(za);
  gkyl_array_release(la);
  gkyl_array_release(ca);
  gkyl_array_release(lam);
  gkyl_array_release(lam_c);
  hx_cell_release(&c);
}

static void test_stage_weights_1x(void) { test_stage_weights(1); }
static void test_stage_weights_2x(void) { test_stage_weights(2); }
static void test_stage_weights_3x(void) { test_stage_weights(3); }

TEST_LIST = {
  {"test_transforms_1x", test_transforms_1x},
  {"test_transforms_2x", test_transforms_2x},
  {"test_transforms_3x", test_transforms_3x},
  {"test_damp_step_1x", test_damp_step_1x},
  {"test_damp_step_2x", test_damp_step_2x},
  {"test_damp_step_3x", test_damp_step_3x},
  {"test_inactive", test_inactive},
  {"test_rate", test_rate},
  {"test_mix_relax", test_mix_relax},
  {"test_mix_stiff", test_mix_stiff},
  {"test_mix_inactive", test_mix_inactive},
  {"test_stage_weights_1x", test_stage_weights_1x},
  {"test_stage_weights_2x", test_stage_weights_2x},
  {"test_stage_weights_3x", test_stage_weights_3x},
  {NULL, NULL},
};
