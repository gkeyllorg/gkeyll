// Tests of the first-harmonic coefficients of the normalized PKPM kinetic equations at a point
// (gkyl_pkpm_harm_coeffs.h): every coefficient and term against an independent construction with
// explicit 3x3 matrices (I_perp = I - b b^T, U, B_1 = xi grad b), the bounds of the spectrum of the
// coupling part of the w flux (sign changes of det(T - lambda I): all six roots real and inside),
// the pairing that keeps d_t b out of the energy row, and the w flux from node values (the stored
// coefficients evaluated at a node) against the pointwise definitions.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_pkpm_harm_coeffs.h>

#include "ctest_pkpm_harm_ref.h"

// Deterministic pseudo-random numbers in [-1, 1].
static double
rnd(unsigned *state)
{
  *state = *state * 1664525u + 1013904223u;
  return 2.0 * ((*state >> 8) / 16777216.0) - 1.0;
}

struct point {
  int cdim;
  double b[3], u[3], xi, theta, dtb[3], gu[9], gb[9], gxi[3], gln[3], qbym, B[3];
};

static void
random_point(struct point *p, int cdim, unsigned *st, bool unit_b)
{
  p->cdim = cdim;
  double bb = 0.0;
  for (int a = 0; a < 3; ++a) {
    p->b[a] = rnd(st);
    bb += p->b[a] * p->b[a];
    p->u[a] = rnd(st);
    p->dtb[a] = rnd(st);
    p->B[a] = rnd(st);
    p->gxi[a] = rnd(st);
    p->gln[a] = rnd(st);
  }
  if (unit_b) {
    for (int a = 0; a < 3; ++a) {
      p->b[a] /= sqrt(bb);
    }
  }
  for (int k = 0; k < 9; ++k) {
    p->gu[k] = rnd(st);
    p->gb[k] = rnd(st);
  }
  p->xi = 1.0 + 0.5 * rnd(st);
  p->theta = 1.0 + 0.5 * rnd(st);
  p->qbym = 2.0 * rnd(st);
}

static void
set_coeffs(struct pkpm_harm_coeffs *c, const struct point *p)
{
  pkpm_harm_coeffs_set(
    c, p->cdim, p->b, p->u, p->xi, p->theta, p->dtb, p->gu, p->gb, p->gxi, p->gln, p->qbym, p->B
  );
}

// out = A v for a 3x3 matrix.
static void
matvec(double A[3][3], const double *v, double *out)
{
  for (int a = 0; a < 3; ++a) {
    out[a] = A[a][0] * v[0] + A[a][1] * v[1] + A[a][2] * v[2];
  }
}

// Independent construction of the terms from the definitions of the notes with explicit matrices.
struct ref_terms {
  double flux[8], src[8], T[8];
};

static void
reference_terms(const struct point *p, double w, const double *Ux, struct ref_terms *r)
{
  int cdim = p->cdim;
  double P[3][3], U[3][3], UT[3][3], Db[3][3], B1[3][3], B1s[3][3], Us[3][3];
  double g[3] = {0.0}, ell[3] = {0.0};
  for (int a = 0; a < 3; ++a) {
    for (int d = 0; d < 3; ++d) {
      P[a][d] = (a == d ? 1.0 : 0.0) - p->b[a] * p->b[d];
      U[a][d] = d < cdim ? p->gu[d * 3 + a] : 0.0;
      Db[a][d] = d < cdim ? p->gb[d * 3 + a] : 0.0;
      B1[a][d] = p->xi * Db[a][d];
    }
  }
  for (int a = 0; a < 3; ++a) {
    for (int d = 0; d < 3; ++d) {
      UT[a][d] = U[d][a];
      B1s[a][d] = B1[a][d] + B1[d][a];
      Us[a][d] = U[a][d] + U[d][a];
    }
  }
  for (int d = 0; d < cdim; ++d) {
    g[d] = p->gxi[d];
    ell[d] = p->xi * p->gln[d];
  }
  double trU = 0.0, trDb = 0.0, bUb = 0.0, bg = 0.0, bl = 0.0, Ub[3], UTb[3], Dbu[3], Dbb[3];
  matvec(U, p->b, Ub);
  matvec(UT, p->b, UTb);
  matvec(Db, p->u, Dbu);
  matvec(Db, p->b, Dbb);
  for (int a = 0; a < 3; ++a) {
    trU += U[a][a];
    trDb += Db[a][a];
    bUb += p->b[a] * Ub[a];
    bg += p->b[a] * g[a];
    bl += p->b[a] * ell[a];
  }
  double beta = p->xi * trDb;
  // C = D_t b + xi w (b.grad) b, D_t b = I_perp d_t B/|B| + (u.grad) b
  double v[3], C[3], t[3], fp[3], Pd[3];
  matvec(P, p->dtb, Pd);
  for (int a = 0; a < 3; ++a) {
    C[a] = Pd[a] + Dbu[a] + p->xi * w * Dbb[a];
  }
  for (int a = 0; a < 3; ++a) {
    v[a] = C[a] - UTb[a] - w * g[a];
  }
  matvec(P, v, t);
  for (int a = 0; a < 3; ++a) {
    v[a] = Ub[a] + C[a];
  }
  matvec(P, v, fp);
  for (int a = 0; a < 3; ++a) {
    fp[a] *= -w;
  }
  double PB1s[3][3], PUs[3][3], PU[3][3], PB1[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int d = 0; d < 3; ++d) {
      PB1s[a][d] = PU[a][d] = PUs[a][d] = PB1[a][d] = 0.0;
      for (int k = 0; k < 3; ++k) {
        PB1s[a][d] += P[a][k] * B1s[k][d];
        PUs[a][d] += P[a][k] * Us[k][d];
        PU[a][d] += P[a][k] * U[k][d];
        PB1[a][d] += P[a][k] * B1[k][d];
      }
    }
  }
  double th = p->theta, F = Ux[0], G = Ux[1];
  const double *M = &Ux[2], *N = &Ux[5];
  double H = 4.0 * th * G - 2.0 * th * th * F, O[3];
  for (int a = 0; a < 3; ++a) {
    O[a] = 6.0 * th * N[a] - 6.0 * th * th * M[a];
  }
  // R_B1 = (1/2)(beta I + P (B_1 + B_1^T)), R_E1 = (1/2)(eta I + P (U + U^T) + w P (B_1 + B_1^T) + 2 w bg I)
  double eta = trU - bUb + w * beta + 2.0 * w * bg;
  double RN[3], RO[3], PB1sN[3], PB1sO[3], PUsN[3], PUM[3], PUN[3], PB1M[3], PB1N[3];
  matvec(PB1s, N, PB1sN);
  matvec(PB1s, O, PB1sO);
  matvec(PUs, N, PUsN);
  matvec(PU, M, PUM);
  matvec(PU, N, PUN);
  matvec(PB1, M, PB1M);
  matvec(PB1, N, PB1N);
  double qB[3] = {p->qbym * p->B[0], p->qbym * p->B[1], p->qbym * p->B[2]};
  double MxB[3] = {M[1] * qB[2] - M[2] * qB[1], M[2] * qB[0] - M[0] * qB[2], M[0] * qB[1] - M[1] * qB[0]};
  double NxB[3] = {N[1] * qB[2] - N[2] * qB[1], N[2] * qB[0] - N[0] * qB[2], N[0] * qB[1] - N[1] * qB[0]};
  double plg[3], pl3g[3];
  for (int a = 0; a < 3; ++a) {
    RN[a] = 0.5 * (beta * N[a] + PB1sN[a]);
    RO[a] = 0.5 * (beta * O[a] + PB1sO[a]);
    v[a] = ell[a] + g[a];
  }
  matvec(P, v, plg);
  for (int a = 0; a < 3; ++a) {
    v[a] = ell[a] + 3.0 * g[a];
  }
  matvec(P, v, pl3g);
  double tM = 0.0, tN = 0.0, fM = 0.0, lM = 0.0, l2gN = 0.0, CM = 0.0, CN = 0.0;
  for (int a = 0; a < 3; ++a) {
    tM += t[a] * M[a];
    tN += t[a] * N[a];
    fM += fp[a] * M[a];
    lM += ell[a] * M[a];
    l2gN += (ell[a] + 2.0 * g[a]) * N[a];
    CM += C[a] * M[a];
    CN += C[a] * N[a];
  }
  double ak = -w * bUb - w * w * bg, gam = w * bg, sig = -w * bl;
  r->T[0] = beta * G + tM;
  r->T[1] = beta * H + tN;
  r->flux[0] = ak * F + r->T[0];
  r->flux[1] = ak * G + r->T[1];
  r->src[0] = sig * F - lM;
  r->src[1] = (bUb - trU - w * (beta + 2.0 * bg + bl)) * G + fM - l2gN;
  for (int a = 0; a < 3; ++a) {
    r->T[2 + a] = G * t[a] + RN[a];
    r->T[5 + a] = H * t[a] + RO[a];
    r->flux[2 + a] = ak * M[a] + r->T[2 + a];
    r->flux[5 + a] = ak * N[a] + r->T[5 + a];
    double LM = -PUM[a] - w * PB1M[a] - gam * M[a] + MxB[a];
    double LN = -PUN[a] - w * PB1N[a] - gam * N[a] + NxB[a];
    double RE = 0.5 * (eta * N[a] + PUsN[a] + w * PB1sN[a] + 2.0 * gam * N[a]);
    r->src[2 + a] = fp[a] * F + sig * M[a] + LM - G * plg[a] - p->b[a] * (CM + beta * G);
    r->src[5 + a] = 2.0 * G * fp[a] + sig * N[a] + LN - RE - H * pl3g[a] - p->b[a] * (CN + beta * H);
  }
}

static void
test_terms(int cdim)
{
  unsigned st = 12345u + cdim;
  double max_err = 0.0, max_val = 0.0;
  for (int it = 0; it < 200; ++it) {
    struct point p;
    random_point(&p, cdim, &st, it % 2 == 0);
    struct pkpm_harm_coeffs c;
    set_coeffs(&c, &p);
    double w = 3.0 * rnd(&st), U[8];
    for (int k = 0; k < 8; ++k) {
      U[k] = rnd(&st);
    }
    struct ref_terms r;
    reference_terms(&p, w, U, &r);
    double T[8], flux[8], src[8];
    pkpm_harm_coupling_flux(&c, w, U, T);
    pkpm_harm_flux_known(&c, w, U, flux);
    pkpm_harm_source(&c, w, U, src);
    for (int k = 0; k < 8; ++k) {
      max_err = fmax(max_err, fabs(T[k] - r.T[k]));
      max_err = fmax(max_err, fabs(flux[k] - r.flux[k]));
      max_err = fmax(max_err, fabs(src[k] - r.src[k]));
      max_val = fmax(max_val, fmax(fabs(r.T[k]), fmax(fabs(r.flux[k]), fabs(r.src[k]))));
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-13 * max_val);
  TEST_MSG("cdim=%d: max term error %g (max value %g)", cdim, max_err, max_val);
}

// Determinant of an n x n matrix (Gaussian elimination with partial pivoting; A is overwritten).
static double
det_n(int n, double A[6][6])
{
  double det = 1.0;
  for (int k = 0; k < n; ++k) {
    int piv = k;
    for (int i = k + 1; i < n; ++i) {
      if (fabs(A[i][k]) > fabs(A[piv][k])) {
        piv = i;
      }
    }
    if (A[piv][k] == 0.0) {
      return 0.0;
    }
    if (piv != k) {
      for (int j = 0; j < n; ++j) {
        double tmp = A[k][j];
        A[k][j] = A[piv][j];
        A[piv][j] = tmp;
      }
      det = -det;
    }
    det *= A[k][k];
    for (int i = k + 1; i < n; ++i) {
      double fac = A[i][k] / A[k][k];
      for (int j = k; j < n; ++j) {
        A[i][j] -= fac * A[k][j];
      }
    }
  }
  return det;
}

// The coupling T on the six independent components [F, G, M0 . e_1, M0 . e_2, N . e_1, N . e_2]
// (e_1, e_2 an orthonormal basis of the plane perpendicular to the unit b): its spectrum is real
// and lies in the bounds. Counted by the sign changes of det(T - lambda I) on a fine grid:
// exactly six inside [lo, hi], none outside.
static void
test_speed_bounds(int cdim)
{
  unsigned st = 777u + cdim;
  int num_cases = 0, num_ok = 0;
  double max_slack = 0.0;
  for (int it = 0; it < 40; ++it) {
    struct point p;
    random_point(&p, cdim, &st, true);
    struct pkpm_harm_coeffs c;
    set_coeffs(&c, &p);
    double w = 3.0 * rnd(&st);
    double e1[3], e2[3];
    // e_1 from the coordinate axis least aligned with b, e_2 = b x e_1
    int amin = 0;
    for (int a = 1; a < 3; ++a) {
      if (fabs(p.b[a]) < fabs(p.b[amin])) {
        amin = a;
      }
    }
    double ax[3] = {0.0, 0.0, 0.0};
    ax[amin] = 1.0;
    pkpm_harm_perp(p.b, ax, e1);
    double n1 = sqrt(pkpm_harm_dot(e1, e1));
    for (int a = 0; a < 3; ++a) {
      e1[a] /= n1;
    }
    e2[0] = p.b[1] * e1[2] - p.b[2] * e1[1];
    e2[1] = p.b[2] * e1[0] - p.b[0] * e1[2];
    e2[2] = p.b[0] * e1[1] - p.b[1] * e1[0];
    double T[6][6];
    for (int k = 0; k < 6; ++k) {
      double U[8] = {0.0}, out[8];
      if (k < 2) {
        U[k] = 1.0;
      } else {
        const double *e = (k % 2 == 0) ? e1 : e2;
        int off = k < 4 ? 2 : 5;
        for (int a = 0; a < 3; ++a) {
          U[off + a] = e[a];
        }
      }
      pkpm_harm_coupling_flux(&c, w, U, out);
      T[0][k] = out[0];
      T[1][k] = out[1];
      T[2][k] = pkpm_harm_dot(&out[2], e1);
      T[3][k] = pkpm_harm_dot(&out[2], e2);
      T[4][k] = pkpm_harm_dot(&out[5], e1);
      T[5][k] = pkpm_harm_dot(&out[5], e2);
    }
    double lo, hi;
    pkpm_harm_speed_bounds(&c, w, &lo, &hi);
    double span = hi - lo, a0 = lo - span, a1 = hi + span;
    int ngrid = 200000, inside = 0, outside = 0;
    double prev = 0.0, lam_min = 1e300, lam_max = -1e300;
    for (int i = 0; i <= ngrid; ++i) {
      double lam = a0 + (a1 - a0) * i / ngrid, A[6][6];
      for (int r = 0; r < 6; ++r) {
        for (int s = 0; s < 6; ++s) {
          A[r][s] = T[r][s] - (r == s ? lam : 0.0);
        }
      }
      double d = det_n(6, A);
      if (i > 0 && d * prev < 0.0) {
        double lam_c = lam - 0.5 * (a1 - a0) / ngrid;
        if (lam_c >= lo && lam_c <= hi) {
          inside += 1;
        } else {
          outside += 1;
        }
        lam_min = fmin(lam_min, lam_c);
        lam_max = fmax(lam_max, lam_c);
      }
      prev = d;
    }
    num_cases += 1;
    if (inside == 6 && outside == 0) {
      num_ok += 1;
      max_slack = fmax(max_slack, fmax(lam_min - lo, hi - lam_max) / span);
    } else {
      TEST_MSG("cdim=%d case %d: %d roots inside [%g, %g], %d outside", cdim, it, inside, lo, hi, outside);
    }
  }
  TEST_CHECK(num_ok == num_cases);
  TEST_MSG("cdim=%d: %d of %d cases with six real roots inside the bounds (max relative slack %g)", cdim,
    num_ok, num_cases, max_slack);
}

// d_t b enters the energy row only through the paired terms 2 w t.M0 (moment of the F flux) and
// 2 f_perp.M0 (G source), which cancel: the combination is independent of d_t B/|B|.
static void
test_energy_pairing(int cdim)
{
  unsigned st = 4242u + cdim;
  double max_diff = 0.0, max_val = 0.0;
  for (int it = 0; it < 100; ++it) {
    struct point p;
    random_point(&p, cdim, &st, it % 2 == 0);
    double w = 3.0 * rnd(&st), U[8];
    for (int k = 0; k < 8; ++k) {
      U[k] = rnd(&st);
    }
    double e[2];
    for (int s = 0; s < 2; ++s) {
      if (s == 1) {
        for (int a = 0; a < 3; ++a) {
          p.dtb[a] += 5.0 * rnd(&st);
        }
      }
      struct pkpm_harm_coeffs c;
      set_coeffs(&c, &p);
      double flux[8], src[8];
      pkpm_harm_flux_known(&c, w, U, flux);
      pkpm_harm_source(&c, w, U, src);
      e[s] = 2.0 * w * flux[0] + 2.0 * src[1];
      max_val = fmax(max_val, fabs(2.0 * w * flux[0]));
    }
    max_diff = fmax(max_diff, fabs(e[1] - e[0]));
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_diff < 1e-13 * max_val);
  TEST_MSG("cdim=%d: energy-row dependence on d_t B/|B| %g (scale %g)", cdim, max_diff, max_val);
}

// The w flux from node values (pkpm_harm_node_*, what the w-face and closure kernels use): with the
// node values formed from the pointwise coefficients, alpha, the coupling T and the speed bounds
// equal the pointwise definitions above.
static void
test_node_functions(int cdim)
{
  unsigned st = 4242u + cdim;
  double max_err_T = 0.0, max_T = 0.0, max_err_b = 0.0, max_b = 0.0, max_err_a = 0.0;
  for (int it = 0; it < 40; ++it) {
    struct point p;
    random_point(&p, cdim, &st, true);
    struct pkpm_harm_coeffs c;
    set_coeffs(&c, &p);
    double th = c.theta, a0 = rnd(&st), ch = rnd(&st), w = 3.0 * rnd(&st);
    double nv[PKPM_HS_NODE_NV] = {0.0};
    nv[PKPM_HS_KW] = -c.K;
    nv[PKPM_HS_BG] = -c.bg;
    for (int a = 0; a < 3; ++a) {
      double tj[2] = {c.t0[a], c.t1[a]};
      for (int j = 0; j < 2; ++j) {
        nv[PKPM_HS_TF(a, j)] = tj[j];
        nv[PKPM_HS_TG(a, j)] = tj[j];
        nv[PKPM_HS_TM(a, j)] = tj[j];
        nv[PKPM_HS_TNG(a, j)] = 4.0 * th * tj[j];
        nv[PKPM_HS_TNF(a, j)] = -2.0 * th * th * tj[j];
      }
    }
    for (int col = 0; col < 3; ++col) {
      double e[3] = {0.0, 0.0, 0.0}, Re[3];
      e[col] = 1.0;
      pkpm_harm_rb1(&c, e, Re);
      for (int a = 0; a < 3; ++a) {
        nv[PKPM_HS_RMN(a, col)] = Re[a];
        nv[PKPM_HS_RNN(a, col)] = 6.0 * th * Re[a];
        nv[PKPM_HS_RNM(a, col)] = -6.0 * th * th * Re[a];
      }
    }
    nv[PKPM_HS_BF] = c.beta;
    nv[PKPM_HS_BGG] = 4.0 * th * c.beta;
    nv[PKPM_HS_BGF] = -2.0 * th * th * c.beta;
    nv[PKPM_HS_NODE_THETA] = th;
    for (int a = 0; a < 3; ++a) {
      nv[PKPM_HS_NODE_B + a] = p.b[a];
    }
    nv[PKPM_HS_NODE_A0] = a0;
    nv[PKPM_HS_NODE_CH] = ch;
    pkpm_harm_node_spec(nv);

    double U[8], T[8], Tref[8], lo, hi, lo_ref, hi_ref;
    for (int k = 0; k < 8; ++k) {
      U[k] = rnd(&st);
    }
    pkpm_harm_node_coupling(nv, w, U, T);
    pkpm_harm_coupling_flux(&c, w, U, Tref);
    for (int k = 0; k < 8; ++k) {
      max_err_T = fmax(max_err_T, fabs(T[k] - Tref[k]));
      max_T = fmax(max_T, fabs(Tref[k]));
    }
    pkpm_harm_node_bounds(nv, w, &lo, &hi);
    pkpm_harm_speed_bounds(&c, w, &lo_ref, &hi_ref);
    max_err_b = fmax(max_err_b, fmax(fabs(lo - lo_ref), fabs(hi - hi_ref)));
    max_b = fmax(max_b, fmax(fabs(lo_ref), fabs(hi_ref)));
    max_err_a = fmax(max_err_a, fabs(pkpm_harm_node_alpha(nv, w) - (a0 + (ch - c.K) * w - c.bg * w * w)));
  }
  TEST_CHECK(max_err_T < 1e-13 * max_T);
  TEST_MSG("cdim=%d: coupling from node values, max error %g (max %g)", cdim, max_err_T, max_T);
  TEST_CHECK(max_err_b < 1e-12 * max_b);
  TEST_MSG("cdim=%d: speed bounds from node values, max error %g (max %g)", cdim, max_err_b, max_b);
  TEST_CHECK(max_err_a < 1e-13);
  TEST_MSG("cdim=%d: alpha from node values, max error %g", cdim, max_err_a);
}

void test_node_functions_1x() { test_node_functions(1); }
void test_node_functions_2x() { test_node_functions(2); }
void test_node_functions_3x() { test_node_functions(3); }
void test_terms_1x() { test_terms(1); }
void test_terms_2x() { test_terms(2); }
void test_terms_3x() { test_terms(3); }
void test_speed_bounds_1x() { test_speed_bounds(1); }
void test_speed_bounds_2x() { test_speed_bounds(2); }
void test_speed_bounds_3x() { test_speed_bounds(3); }
void test_energy_pairing_1x() { test_energy_pairing(1); }
void test_energy_pairing_3x() { test_energy_pairing(3); }

TEST_LIST = {
  {"test_terms_1x", test_terms_1x},
  {"test_terms_2x", test_terms_2x},
  {"test_terms_3x", test_terms_3x},
  {"test_speed_bounds_1x", test_speed_bounds_1x},
  {"test_speed_bounds_2x", test_speed_bounds_2x},
  {"test_speed_bounds_3x", test_speed_bounds_3x},
  {"test_energy_pairing_1x", test_energy_pairing_1x},
  {"test_energy_pairing_3x", test_energy_pairing_3x},
  {"test_node_functions_1x", test_node_functions_1x},
  {"test_node_functions_2x", test_node_functions_2x},
  {"test_node_functions_3x", test_node_functions_3x},
  {NULL, NULL},
};
