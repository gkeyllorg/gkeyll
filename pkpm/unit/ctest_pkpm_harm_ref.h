#pragma once

// Reference definitions for the unit tests of the first harmonic of the normalized PKPM kinetic
// equations at one point of configuration space: C, t, f_perp, R_B1, L_perp, R_E1, the coupling
// part of the w flux, its speed bounds and the known sources, written out from the equations.
// The kernels form the same quantities at their own nodes from the stored operator
// (gkylcas/pkpm/pkpm/pkpm-harm-ops.mac); the tests compare the two.
//
// The distributions are the eight blocks [F, G, M0_perp (x, y, z), N_perp (x, y, z)]. The inputs
// are the values at the point of p=1 fields; every product is formed here, at the point, with
// I_perp v = v - b (b.v) for the interpolated b. The coefficients are the known parts: the closure
// parts (c_0 in sigma_0, alpha_0 = b.a, c in Gamma_w, a_perp = I_perp a) are added by the kernels
// of the velocity-space pass.

#include <math.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_util.h>

struct pkpm_harm_coeffs {
  double b[3]; // Field direction.
  double theta; // T_perp/T, the Laguerre scale of the H and O closures.
  double qbB[3]; // (q/m) B: the gyration Omega v x b = (q/m) v x B.
  double U[3][3]; // U[a][d] = d_d u_a (zero for d >= cdim).
  double B1[3][3]; // B_1 = xi grad b: B1[a][d] = xi d_d b_a (zero for d >= cdim).
  double g[3]; // grad xi (zero beyond cdim).
  double ell[3]; // ell_xi = xi grad ln n (zero beyond cdim).
  double K; // b.U b.
  double divu; // div u (physical compression; the closure c_0 replaces it in sigma_0).
  double beta; // xi div b, the mirror rate (as in the F and G equations).
  double bg; // b.grad xi: Gamma_w = -c + w bg.
  double bl; // b.ell_xi: sigma_0 = c_0 - w bl.
  double C0[3],
    C1[3]; // C = C0 + w C1 = D_t b + xi w (b.grad) b, D_t b = I_perp d_t B/|B| + (u.grad) b.
  double t0[3], t1[3]; // t = t0 + w t1 = I_perp (C - U^T b - w g).
  double fp1[3], fp2[3]; // f_perp (known part) = w fp1 + w^2 fp2 = -w I_perp (U b + C).
  double plg[3], pl3g[3]; // I_perp (ell_xi + g), I_perp (ell_xi + 3 g).
};

/**
 * Form the known first-harmonic coefficients at a point.
 *
 * @param c Output coefficients.
 * @param cdim Configuration-space dimensions (the gradients have cdim directions).
 * @param b Field direction b at the point.
 * @param u Flow velocity at the point.
 * @param xi Thermal speed at the point.
 * @param theta T_perp/T at the point.
 * @param dtb d_t B/|B| at the point.
 * @param gu Gradient of u at the point, gu[d*3 + j] = d_d u_j.
 * @param gb Weak gradient of b at the point, gb[d*3 + j] = d_d b_j.
 * @param gxi Gradient of xi at the point, gxi[d] = d_d xi.
 * @param gln grad ln n at the point (cdim values), or 0 when the caller needs no ell_xi.
 * @param qbym Charge to mass ratio of the species.
 * @param B Magnetic field at the point, or 0 when the caller needs no gyration.
 */
GKYL_CU_DH static inline void
pkpm_harm_coeffs_set(
  struct pkpm_harm_coeffs *c, int cdim, const double *b, const double *u, double xi, double theta,
  const double *dtb, const double *gu, const double *gb, const double *gxi, const double *gln,
  double qbym, const double *B
)
{
  for (int a = 0; a < 3; ++a) {
    c->b[a] = b[a];
    c->qbB[a] = B ? qbym * B[a] : 0.0;
    c->g[a] = 0.0;
    c->ell[a] = 0.0;
    for (int d = 0; d < 3; ++d) {
      c->U[a][d] = 0.0;
      c->B1[a][d] = 0.0;
    }
  }
  c->theta = theta;
  for (int d = 0; d < cdim; ++d) {
    c->g[d] = gxi[d];
    c->ell[d] = gln ? xi * gln[d] : 0.0;
    for (int a = 0; a < 3; ++a) {
      c->U[a][d] = gu[d * 3 + a];
      c->B1[a][d] = xi * gb[d * 3 + a];
    }
  }

  double K = 0.0, divu = 0.0, divb = 0.0, bg = 0.0, bl = 0.0;
  double Ub[3] = {0.0}, UTb[3] = {0.0}, ugb[3] = {0.0}, bgb[3] = {0.0};
  for (int d = 0; d < cdim; ++d) {
    divu += c->U[d][d];
    divb += gb[d * 3 + d];
    bg += b[d] * c->g[d];
    bl += b[d] * c->ell[d];
    for (int a = 0; a < 3; ++a) {
      K += b[d] * b[a] * c->U[a][d];
      Ub[a] += c->U[a][d] * b[d];
      UTb[d] += c->U[a][d] * b[a];
      ugb[a] += u[d] * gb[d * 3 + a]; // (u.grad) b
      bgb[a] += xi * b[d] * gb[d * 3 + a]; // xi (b.grad) b
    }
  }
  c->K = K;
  c->divu = divu;
  c->beta = xi * divb;
  c->bg = bg;
  c->bl = bl;

  // C as defined: I_perp only inside d_t b = I_perp d_t B/|B|; the
  // convective derivatives (u.grad) b and xi (b.grad) b are perpendicular in the continuum and get
  // no projection of their own (t and f_perp carry theirs).
  double v[3];
  pkpm_harm_perp(b, dtb, c->C0);
  for (int a = 0; a < 3; ++a) {
    c->C0[a] += ugb[a];
    c->C1[a] = bgb[a];
  }
  for (int a = 0; a < 3; ++a) {
    v[a] = c->C0[a] - UTb[a];
  }
  pkpm_harm_perp(b, v, c->t0);
  for (int a = 0; a < 3; ++a) {
    v[a] = c->C1[a] - c->g[a];
  }
  pkpm_harm_perp(b, v, c->t1);
  for (int a = 0; a < 3; ++a) {
    v[a] = -(Ub[a] + c->C0[a]);
  }
  pkpm_harm_perp(b, v, c->fp1);
  for (int a = 0; a < 3; ++a) {
    v[a] = -c->C1[a];
  }
  pkpm_harm_perp(b, v, c->fp2);
  for (int a = 0; a < 3; ++a) {
    v[a] = c->ell[a] + c->g[a];
  }
  pkpm_harm_perp(b, v, c->plg);
  for (int a = 0; a < 3; ++a) {
    v[a] = c->ell[a] + 3.0 * c->g[a];
  }
  pkpm_harm_perp(b, v, c->pl3g);
}

// t at w.
GKYL_CU_DH static inline void
pkpm_harm_t(const struct pkpm_harm_coeffs *c, double w, double *t)
{
  for (int a = 0; a < 3; ++a) {
    t[a] = c->t0[a] + w * c->t1[a];
  }
}

// out = I_perp (B_1 + B_1^T) v.
GKYL_CU_DH static inline void
pkpm_harm_b1sym(const struct pkpm_harm_coeffs *c, const double *v, double *out)
{
  double s[3];
  for (int a = 0; a < 3; ++a) {
    s[a] = 0.0;
    for (int d = 0; d < 3; ++d) {
      s[a] += c->B1[a][d] * v[d] + c->B1[d][a] * v[d];
    }
  }
  pkpm_harm_perp(c->b, s, out);
}

// out = R_B1 v = (1/2) [beta v + I_perp (B_1 + B_1^T) v].
GKYL_CU_DH static inline void
pkpm_harm_rb1(const struct pkpm_harm_coeffs *c, const double *v, double *out)
{
  double s[3];
  pkpm_harm_b1sym(c, v, s);
  for (int a = 0; a < 3; ++a) {
    out[a] = 0.5 * (c->beta * v[a] + s[a]);
  }
}

// out = L_perp v (known part) = -I_perp U v - w I_perp B_1 v - w bg v + (q/m) v x B.
GKYL_CU_DH static inline void
pkpm_harm_lperp(const struct pkpm_harm_coeffs *c, double w, const double *v, double *out)
{
  double Uv[3], Bv[3];
  for (int a = 0; a < 3; ++a) {
    Uv[a] = 0.0;
    Bv[a] = 0.0;
    for (int d = 0; d < 3; ++d) {
      Uv[a] += c->U[a][d] * v[d];
      Bv[a] += c->B1[a][d] * v[d];
    }
  }
  pkpm_harm_perp(c->b, Uv, Uv);
  pkpm_harm_perp(c->b, Bv, Bv);
  double gam = w * c->bg;
  const double *qb = c->qbB;
  double vxB[3] = {
    v[1] * qb[2] - v[2] * qb[1], v[2] * qb[0] - v[0] * qb[2], v[0] * qb[1] - v[1] * qb[0]
  };
  for (int a = 0; a < 3; ++a) {
    out[a] = -Uv[a] - w * Bv[a] - gam * v[a] + vxB[a];
  }
}

// out = R_E1 v (known part) = (1/2) [eta v + I_perp (U + U^T) v + w I_perp (B_1 + B_1^T) v + 2 w bg v],
// eta = div u - K + w beta + 2 w bg.
GKYL_CU_DH static inline void
pkpm_harm_re1(const struct pkpm_harm_coeffs *c, double w, const double *v, double *out)
{
  double Us[3], Bs[3];
  for (int a = 0; a < 3; ++a) {
    Us[a] = 0.0;
    for (int d = 0; d < 3; ++d) {
      Us[a] += c->U[a][d] * v[d] + c->U[d][a] * v[d];
    }
  }
  pkpm_harm_perp(c->b, Us, Us);
  pkpm_harm_b1sym(c, v, Bs);
  double gam = w * c->bg;
  double eta = c->divu - c->K + w * c->beta + 2.0 * gam;
  for (int a = 0; a < 3; ++a) {
    out[a] = 0.5 * (eta * v[a] + Us[a] + w * Bs[a] + 2.0 * gam * v[a]);
  }
}

/**
 * Coupling part T(U) of the w flux, alpha U + T(U), at w:
 *   T_F = beta G + t.M0,  T_G = beta H + t.N,  T_M0 = G t + R_B1 N,  T_N = H t + R_B1 O,
 * with the closures H = 4 theta G - 2 theta^2 F, O = 6 theta N - 6 theta^2 M0.
 */
GKYL_CU_DH static inline void
pkpm_harm_coupling_flux(const struct pkpm_harm_coeffs *c, double w, const double *U, double *out)
{
  double th = c->theta, t[3], O[3], RN[3], RO[3];
  pkpm_harm_t(c, w, t);
  double F = U[0], G = U[1];
  const double *M = &U[2], *N = &U[5];
  double H = 4.0 * th * G - 2.0 * th * th * F;
  for (int a = 0; a < 3; ++a) {
    O[a] = 6.0 * th * N[a] - 6.0 * th * th * M[a];
  }
  pkpm_harm_rb1(c, N, RN);
  pkpm_harm_rb1(c, O, RO);
  out[0] = c->beta * G + pkpm_harm_dot(t, M);
  out[1] = c->beta * H + pkpm_harm_dot(t, N);
  for (int a = 0; a < 3; ++a) {
    out[2 + a] = G * t[a] + RN[a];
    out[5 + a] = H * t[a] + RO[a];
  }
}

/**
 * Bounds of the spectrum of the coupling T at w (real):
 *   [min(cand) - |t| sqrt((2 + sqrt 2) theta), max(cand) + |t| sqrt((2 + sqrt 2) theta)],
 *   cand = {beta theta (2 +- sqrt 2), r theta (3 +- sqrt 3)}, r the eigenvalues of R_B1 on the
 * plane perpendicular to b. Exact for t = 0 and for B_1 = 0.
 */
GKYL_CU_DH static inline void
pkpm_harm_speed_bounds(const struct pkpm_harm_coeffs *c, double w, double *lo, double *hi)
{
  double th = fmax(c->theta, 0.0), t[3];
  pkpm_harm_t(c, w, t);
  double tn = sqrt(pkpm_harm_dot(t, t));

  // S = B_1 + B_1^T restricted to the plane perpendicular to the unit direction of b.
  double S[3][3], bh[3], bb = pkpm_harm_dot(c->b, c->b);
  double ib = bb > 0.0 ? 1.0 / sqrt(bb) : 0.0;
  for (int a = 0; a < 3; ++a) {
    bh[a] = c->b[a] * ib;
    for (int d = 0; d < 3; ++d) {
      S[a][d] = c->B1[a][d] + c->B1[d][a];
    }
  }
  double Sb[3], trS = 0.0, trS2 = 0.0;
  for (int a = 0; a < 3; ++a) {
    trS += S[a][a];
    Sb[a] = 0.0;
    for (int d = 0; d < 3; ++d) {
      Sb[a] += S[a][d] * bh[d];
      trS2 += S[a][d] * S[a][d];
    }
  }
  double bSb = pkpm_harm_dot(bh, Sb);
  double trE = trS - bSb, trE2 = trS2 - 2.0 * pkpm_harm_dot(Sb, Sb) + bSb * bSb;
  double disc = sqrt(fmax(0.5 * trE2 - 0.25 * trE * trE, 0.0));
  double r_lo = 0.5 * (c->beta + 0.5 * trE - disc), r_hi = 0.5 * (c->beta + 0.5 * trE + disc);

  const double s2m = 0.5857864376269049, s2p = 3.414213562373095; // 2 -+ sqrt 2
  const double s3m = 1.2679491924311228, s3p = 4.732050807568877; // 3 -+ sqrt 3
  double cand[6] = {c->beta * th * s2m, c->beta * th * s2p, r_lo * th * s3m,
                    r_lo * th * s3p,    r_hi * th * s3m,    r_hi * th * s3p};
  double cmin = cand[0], cmax = cand[0];
  for (int k = 1; k < 6; ++k) {
    cmin = fmin(cmin, cand[k]);
    cmax = fmax(cmax, cand[k]);
  }
  double ts = tn * sqrt(s2p * th);
  *lo = cmin - ts;
  *hi = cmax + ts;
}

/**
 * Known part of the w flux of the eight blocks at w (the w-flux entries of the stored operator,
 * pkpm-harm-ops.mac): alpha_k U + T(U) for every block, alpha_k = -w K - w^2 b.grad xi.
 */
GKYL_CU_DH static inline void
pkpm_harm_flux_known(const struct pkpm_harm_coeffs *c, double w, const double *U, double *out)
{
  double T[8];
  pkpm_harm_coupling_flux(c, w, U, T);
  double ak = -w * c->K - w * w * c->bg;
  for (int k = 0; k < 8; ++k) {
    out[k] = ak * U[k] + T[k];
  }
}

/**
 * Known part of the sources of the eight blocks at w (the source entries of the stored operator,
 * pkpm-harm-ops.mac; the gyration included):
 *   F:  sigma_0 F - ell_xi.M0
 *   G:  (K - div u - w (beta + 2 b.grad xi + b.ell_xi)) G + f_perp.M0 - (ell_xi + 2 g).N
 *   M0: f_perp F + sigma_0 M0 + L_perp M0 - G I_perp (ell_xi + g) - b (C.M0 + beta G)
 *   N:  2 G f_perp + sigma_0 N + L_perp N - R_E1 N - H I_perp (ell_xi + 3 g) - b (C.N + beta H)
 * with sigma_0 = -w b.ell_xi, f_perp = -w I_perp (U b + C) and the known L_perp, R_E1.
 */
GKYL_CU_DH static inline void
pkpm_harm_source(const struct pkpm_harm_coeffs *c, double w, const double *U, double *out)
{
  double th = c->theta;
  double F = U[0], G = U[1];
  const double *M = &U[2], *N = &U[5];
  double H = 4.0 * th * G - 2.0 * th * th * F;
  double sig = -w * c->bl;
  double fp[3], C[3], l2g[3], LM[3], LN[3], RN[3];
  for (int a = 0; a < 3; ++a) {
    fp[a] = w * c->fp1[a] + w * w * c->fp2[a];
    C[a] = c->C0[a] + w * c->C1[a];
    l2g[a] = c->ell[a] + 2.0 * c->g[a];
  }
  pkpm_harm_lperp(c, w, M, LM);
  pkpm_harm_lperp(c, w, N, LN);
  pkpm_harm_re1(c, w, N, RN);
  double CM = pkpm_harm_dot(C, M), CN = pkpm_harm_dot(C, N);
  out[0] = sig * F - pkpm_harm_dot(c->ell, M);
  out[1] = (c->K - c->divu - w * (c->beta + 2.0 * c->bg + c->bl)) * G + pkpm_harm_dot(fp, M) -
           pkpm_harm_dot(l2g, N);
  for (int a = 0; a < 3; ++a) {
    out[2 + a] = fp[a] * F + sig * M[a] + LM[a] - G * c->plg[a] - c->b[a] * (CM + c->beta * G);
    out[5 + a] =
      2.0 * G * fp[a] + sig * N[a] + LN[a] - RN[a] - H * c->pl3g[a] - c->b[a] * (CN + c->beta * H);
  }
}
