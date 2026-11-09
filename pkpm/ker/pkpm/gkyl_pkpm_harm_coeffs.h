#pragma once

// First-harmonic coefficients of the normalized PKPM kinetic equations at one point of
// configuration space: the single definition of C, t, f_perp, R_B1, L_perp, R_E1, the coupling
// part of the w flux and its speed bounds, and the configuration-space face flux, shared by the
// code that needs them at its own evaluation points (w surface-flux kernel pkpm_vel_flux_vpar_*,
// closure kernel pkpm_vars_closure_*, the rate estimate of gkyl_dg_pkpm_harm_vol, configuration
// face-flux kernels pkpm_conf_flux_* and pkpm_fluid_flux_accum_*). Equations normalized-F/G/M/N
// of First_Fourier_Harmonic_Notes/normalized_first_harmonic.tex; plan sections 17.11 and 17.14.
//
// The distributions are the eight blocks [F, G, M0_perp (x, y, z), N_perp (x, y, z)]. The inputs
// are the values at the point of p=1 fields; every product is formed here, at the point, with
// I_perp v = v - b (b.v) for the interpolated b. The coefficients are the known parts: the closure
// parts (c_0 in sigma_0, alpha_0 = b.a, c in Gamma_w, a_perp = I_perp a) are added by the stage-B
// kernels.

#include <math.h>
#include <gkyl_pkpm_harm_slots.h>
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
  double C0[3], C1[3]; // C = C0 + w C1 = D_t b + xi w (b.grad) b, D_t b = I_perp d_t B/|B| + (u.grad) b.
  double t0[3], t1[3]; // t = t0 + w t1 = I_perp (C - U^T b - w g).
  double fp1[3], fp2[3]; // f_perp (known part) = w fp1 + w^2 fp2 = -w I_perp (U b + C).
  double plg[3], pl3g[3]; // I_perp (ell_xi + g), I_perp (ell_xi + 3 g).
};

// out = I_perp v (out may alias v).
GKYL_CU_DH static inline void
pkpm_harm_perp(const double *b, const double *v, double *out)
{
  double bv = b[0] * v[0] + b[1] * v[1] + b[2] * v[2];
  for (int i = 0; i < 3; ++i) {
    out[i] = v[i] - b[i] * bv;
  }
}

GKYL_CU_DH static inline double
pkpm_harm_dot(const double *a, const double *b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

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

  // C as defined (normalized_first_harmonic.tex): I_perp only inside d_t b = I_perp d_t B/|B|; the
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
  double vxB[3] = {v[1] * qb[2] - v[2] * qb[1], v[2] * qb[0] - v[0] * qb[2], v[0] * qb[1] - v[1] * qb[0]};
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
 * Bounds of the spectrum of the coupling T at w (real; scratch study of plan 17.11):
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
  double cand[6] = {c->beta * th * s2m, c->beta * th * s2p, r_lo * th * s3m, r_lo * th * s3p,
    r_hi * th * s3m, r_hi * th * s3p};
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
    out[2 + a] =
      fp[a] * F + sig * M[a] + LM[a] - G * c->plg[a] - c->b[a] * (CM + c->beta * G);
    out[5 + a] = 2.0 * G * fp[a] + sig * N[a] + LN[a] - RN[a] - H * c->pl3g[a] -
      c->b[a] * (CN + c->beta * H);
  }
}

/**
 * Rate estimate of the known first-harmonic sources at w (time-step bound): the largest
 * block-row sum of the coefficient magnitudes (Gershgorin, vector couplings by their Euclidean
 * norms, operators by their Frobenius norms). The gyration is excluded: it is the exponential part
 * of the time integrator (plan 17.6).
 */
GKYL_CU_DH static inline double
pkpm_harm_source_rate(const struct pkpm_harm_coeffs *c, double w)
{
  double th = fabs(c->theta), nU = 0.0, nB = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int d = 0; d < 3; ++d) {
      nU += c->U[a][d] * c->U[a][d];
      nB += c->B1[a][d] * c->B1[a][d];
    }
  }
  nU = sqrt(nU);
  nB = sqrt(nB);
  double fp[3], C[3], l2g[3];
  for (int a = 0; a < 3; ++a) {
    fp[a] = w * c->fp1[a] + w * w * c->fp2[a];
    C[a] = c->C0[a] + w * c->C1[a];
    l2g[a] = c->ell[a] + 2.0 * c->g[a];
  }
  double nfp = sqrt(pkpm_harm_dot(fp, fp)), nC = sqrt(pkpm_harm_dot(C, C));
  double gam = fabs(w * c->bg), sig = fabs(w * c->bl), aw = fabs(w);
  double eta = fabs(c->divu - c->K + w * c->beta + 2.0 * w * c->bg);
  double lp = nU + aw * nB + gam; // |L_perp| without the gyration
  double re = 0.5 * eta + nU + aw * nB + gam; // |R_E1|
  double hc = 4.0 * th + 2.0 * th * th; // |H| per unit max(|F|, |G|)
  double rF = sqrt(pkpm_harm_dot(c->ell, c->ell));
  double rG = nfp + sqrt(pkpm_harm_dot(l2g, l2g));
  double rM = nfp + sig + lp + sqrt(pkpm_harm_dot(c->plg, c->plg)) + nC + fabs(c->beta);
  double rN = 2.0 * nfp + sig + lp + re + nC +
    hc * (sqrt(pkpm_harm_dot(c->pl3g, c->pl3g)) + fabs(c->beta));
  return fmax(fmax(rF, rG), fmax(rM, rN));
}

// Face coefficients of the configuration-space flux at one face node (plan 17.42), stored with
// stride nn (the number of transverse face nodes) as fc[PKPM_FC_* nn]:
#define PKPM_FC_A0 0 // {u_d}: a = A0 + A1 w is the speed of the flux function
#define PKPM_FC_A1 1 // {xi} b_d
#define PKPM_FC_C 2 // {xi} |b x e_d| sqrt((2 + sqrt 2) max({theta}, 0)): its coupled speed
#define PKPM_FC_XI 3 // {xi}: T_F = XI M0_d, T_G = XI N_d
#define PKPM_FC_P 4 // P_k = {xi} (delta_kd - b_k b_d), k = 0, 1, 2 (b = {b}/|{b}|): T_M0 = G P, T_N = H P
#define PKPM_FC_H1 7 // 4 {theta}: H = H1 G - H2 F
#define PKPM_FC_H2 8 // 2 {theta}^2
#define PKPM_FC_UL 9 // u_d, xi b_d and c of the l side (lower cell)
#define PKPM_FC_SL 10
#define PKPM_FC_CL 11
#define PKPM_FC_UR 12 // u_d, xi b_d and c of the r side (upper cell)
#define PKPM_FC_SR 13
#define PKPM_FC_CR 14
#define PKPM_FC_NUM 15

/**
 * Face coefficients of the flux of the eight blocks through a configuration face with normal e_d
 * at one face node (plans 17.14, 17.41, 17.42; discrete_energy_streaming.tex), from the traces of
 * u_d, xi, b and theta of the cells below (l) and above (r) the face. The flux function is
 * evaluated with the face averages {q} = (q_l + q_r)/2 of its coefficients, the field direction
 * renormalized, b = {b}/|{b}| (zero where {b} = 0, as in d_t B/|B|):
 *   J(X) = a X + T(X),  a = {u_d} + {xi} b_d w,
 *   T_F = {xi} M0_d,  T_G = {xi} N_d,  T_M0 = {xi} G (e_d - b b_d),  T_N = {xi} H (e_d - b b_d),
 * H = 4 {theta} G - 2 {theta}^2 F; a coefficient jump then enters the face only through the
 * average, which the central-trace weak gradients of the coefficients (the frame sources) pair
 * with. The coupled speeds are vn -+ c with c = xi |b x e_d| sqrt((2 + sqrt 2) theta) (the fastest
 * of vn +- xi h sqrt((2 +- sqrt 2) theta)); the per-side speeds are kept for the time step.
 * Written as fc[q*nn], q = PKPM_FC_*.
 */
GKYL_CU_DH static inline void
pkpm_harm_face_coef(
  int d, double un_l, double un_r, double xi_l, double xi_r, const double *b_l, const double *b_r,
  double theta_l, double theta_r, int nn, double *fc
)
{
  const double s2p = 3.414213562373095; // 2 + sqrt 2
  int d1 = (d + 1) % 3, d2 = (d + 2) % 3;
  double xa = 0.5 * (xi_l + xi_r), tha = 0.5 * (theta_l + theta_r), ba[3];
  for (int k = 0; k < 3; ++k) {
    ba[k] = 0.5 * (b_l[k] + b_r[k]);
  }
  double bmag = sqrt(ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]);
  double ibmag = bmag > 0.0 ? 1.0 / bmag : 0.0;
  for (int k = 0; k < 3; ++k) {
    ba[k] *= ibmag;
  }
  fc[PKPM_FC_A0 * nn] = 0.5 * (un_l + un_r);
  fc[PKPM_FC_A1 * nn] = xa * ba[d];
  fc[PKPM_FC_C * nn] = xa * sqrt(ba[d1] * ba[d1] + ba[d2] * ba[d2]) * sqrt(s2p * fmax(tha, 0.0));
  fc[PKPM_FC_XI * nn] = xa;
  for (int k = 0; k < 3; ++k) {
    fc[(PKPM_FC_P + k) * nn] = xa * ((k == d ? 1.0 : 0.0) - ba[k] * ba[d]);
  }
  fc[PKPM_FC_H1 * nn] = 4.0 * tha;
  fc[PKPM_FC_H2 * nn] = 2.0 * tha * tha;
  fc[PKPM_FC_UL * nn] = un_l;
  fc[PKPM_FC_SL * nn] = xi_l * b_l[d];
  fc[PKPM_FC_CL * nn] = xi_l * sqrt(b_l[d1] * b_l[d1] + b_l[d2] * b_l[d2]) * sqrt(s2p * fmax(theta_l, 0.0));
  fc[PKPM_FC_UR * nn] = un_r;
  fc[PKPM_FC_SR * nn] = xi_r * b_r[d];
  fc[PKPM_FC_CR * nn] = xi_r * sqrt(b_r[d1] * b_r[d1] + b_r[d2] * b_r[d2]) * sqrt(s2p * fmax(theta_r, 0.0));
}

/**
 * Flux of the eight blocks through a configuration face with normal e_d at one face node and one
 * w, from the face coefficients fc of the node (pkpm_harm_face_coef, stride nn), as the two side
 * contributions phi_l + phi_r; each side keeps its own trace X_s. The HLL flux of J(X) = a X + T(X),
 *   phi_l = (S_R J(X_l) - S_L S_R X_l)/ds,   phi_r = (S_L S_R X_r - S_L J(X_r))/ds,
 * ds = S_R - S_L, is the product of averages a {X} + {T} plus a nonnegative dissipation. Its bounds
 * are the wave speeds a and a -+ c of the frozen-coefficient flux, S_L = min(a - c, 0) and
 * S_R = max(a + c, 0): the streaming part is a {X} - nu [X]/2 with nu = max(c, |a|) >= |a|, the
 * dissipation the linear energy estimate needs (discrete_energy_streaming.tex, lemma 1). The time
 * step sees the fastest signal at the node
 * instead, the envelope of the speeds of both sides and of the face flux,
 *   min(vn_l - c_l, vn_r - c_r, a - c, 0),   max(vn_l + c_l, vn_r + c_r, a + c, 0),
 * vn_s = u_d,s + xi_s b_d,s w, returned in sb for the rate (plan 17.42). If ds = 0 (then a = 0 and
 * c = 0) the flux is the central (T_l + T_r)/2.
 */
GKYL_CU_DH static inline void
pkpm_harm_conf_flux(
  int d, const double *X_l, const double *X_r, double w, const double *fc, int nn, double *phi_l,
  double *phi_r, double *sb
)
{
  double a = fc[PKPM_FC_A0 * nn] + fc[PKPM_FC_A1 * nn] * w, c = fc[PKPM_FC_C * nn];
  double vn_l = fc[PKPM_FC_UL * nn] + fc[PKPM_FC_SL * nn] * w, c_l = fc[PKPM_FC_CL * nn];
  double vn_r = fc[PKPM_FC_UR * nn] + fc[PKPM_FC_SR * nn] * w, c_r = fc[PKPM_FC_CR * nn];
  double s_l = fmin(a - c, 0.0), s_r = fmax(a + c, 0.0);
  sb[0] = fmin(fmin(vn_l - c_l, vn_r - c_r), s_l);
  sb[1] = fmax(fmax(vn_l + c_l, vn_r + c_r), s_r);

  // the coupling T of both sides
  double xa = fc[PKPM_FC_XI * nn], h1 = fc[PKPM_FC_H1 * nn], h2 = fc[PKPM_FC_H2 * nn];
  double T_l[8], T_r[8];
  T_l[0] = xa * X_l[2 + d];
  T_l[1] = xa * X_l[5 + d];
  T_r[0] = xa * X_r[2 + d];
  T_r[1] = xa * X_r[5 + d];
  double H_l = h1 * X_l[1] - h2 * X_l[0], H_r = h1 * X_r[1] - h2 * X_r[0];
  for (int k = 0; k < 3; ++k) {
    double pk = fc[(PKPM_FC_P + k) * nn];
    T_l[2 + k] = X_l[1] * pk;
    T_l[5 + k] = H_l * pk;
    T_r[2 + k] = X_r[1] * pk;
    T_r[5 + k] = H_r * pk;
  }
  double ds = s_r - s_l;
  if (!(ds > 0.0)) {
    for (int k = 0; k < 8; ++k) {
      phi_l[k] = 0.5 * T_l[k];
      phi_r[k] = 0.5 * T_r[k];
    }
    return;
  }
  double w_l = s_r / ds, w_r = -s_l / ds, diss = s_l * s_r / ds;
  double A_l = w_l * a - diss, A_r = w_r * a + diss;
  for (int k = 0; k < 8; ++k) {
    phi_l[k] = A_l * X_l[k] + w_l * T_l[k];
    phi_r[k] = A_r * X_r[k] + w_r * T_r[k];
  }
}

// Speed bounds of the flux of pkpm_harm_conf_flux at one face node and w, for the time step.
GKYL_CU_DH static inline void
pkpm_harm_conf_bounds(double w, const double *fc, int nn, double *sb)
{
  double a = fc[PKPM_FC_A0 * nn] + fc[PKPM_FC_A1 * nn] * w, c = fc[PKPM_FC_C * nn];
  double vn_l = fc[PKPM_FC_UL * nn] + fc[PKPM_FC_SL * nn] * w, c_l = fc[PKPM_FC_CL * nn];
  double vn_r = fc[PKPM_FC_UR * nn] + fc[PKPM_FC_SR * nn] * w, c_r = fc[PKPM_FC_CR * nn];
  sb[0] = fmin(fmin(vn_l - c_l, vn_r - c_r), fmin(a - c, 0.0));
  sb[1] = fmax(fmax(vn_l + c_l, vn_r + c_r), fmax(a + c, 0.0));
}

// The map of one side of a face to the face frame and back (gkyl_dg_pkpm_remap_flux_surf).
struct pkpm_harm_face_map {
  double r, K, c_par, c_perp[3]; // to the face frame: w_s = r wbar + c_par
  double r_inv, K_inv, c_perp_inv[3]; // from the face frame to this side
};

/**
 * Map of one side of a face to the face frame with averaged rho, u and xi and field direction
 * bbar (physical_jump_flux.tex, section "The exact remap"): r = xibar/xi_s,
 * c = (ubar - u_s)/xi_s split along bbar, K = rho_s xibar/(rhobar xi_s). Needs rho_s, xi_s > 0.
 */
GKYL_CU_DH static inline void
pkpm_harm_face_map(
  double rho_s, const double *u_s, double xi_s, double rho_bar, const double *u_bar, double xi_bar,
  const double *b_bar, struct pkpm_harm_face_map *map
)
{
  double c_fwd[3], c_bck[3];
  for (int k = 0; k < 3; ++k) {
    c_fwd[k] = (u_bar[k] - u_s[k]) / xi_s;
    c_bck[k] = (u_s[k] - u_bar[k]) / xi_bar;
  }
  map->r = xi_bar / xi_s;
  map->K = rho_s * xi_bar / (rho_bar * xi_s);
  map->r_inv = 1.0 / map->r;
  map->K_inv = 1.0 / map->K;
  map->c_par = c_fwd[0] * b_bar[0] + c_fwd[1] * b_bar[1] + c_fwd[2] * b_bar[2];
  double c_par_bck = c_bck[0] * b_bar[0] + c_bck[1] * b_bar[1] + c_bck[2] * b_bar[2];
  for (int k = 0; k < 3; ++k) {
    map->c_perp[k] = c_fwd[k] - map->c_par * b_bar[k];
    map->c_perp_inv[k] = c_bck[k] - c_par_bck * b_bar[k];
  }
}

/**
 * Exact change of frame of the eight blocks X = [F, G, M0_perp, N_perp] at the source label to
 * the blocks Y at the target label (physical_jump_flux.tex, eq. frame-change). The factors do not
 * depend on w, and the map is linear, so it also applies to Legendre coefficients in w.
 */
GKYL_CU_DH static inline void
pkpm_harm_frame_change(const double *X, double K, double r, const double *c, double *Y)
{
  const double F = X[0], G = X[1], *M = &X[2], *N = &X[5];
  double cM = c[0] * M[0] + c[1] * M[1] + c[2] * M[2];
  double cc = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
  Y[0] = K * F;
  Y[1] = K / (r * r) * (G - cM + 0.5 * cc * F);
  for (int k = 0; k < 3; ++k) {
    Y[2 + k] = K / r * (M[k] - c[k] * F);
    Y[5 + k] = K / (r * r * r) * (N[k] - 2.0 * G * c[k] + c[k] * cM + 0.5 * cc * M[k] - 0.5 * cc * F * c[k]);
  }
}

/**
 * Fluid flux [F_rho, F_x, F_y, F_z, F_E] from the moments mom = [I_0, I_1, I_2, I_G, I_M (3)] of
 * the flux of the blocks in a frame with density rho, flow u, thermal speed xi and field
 * direction b: the moment map of the normalized distribution to the lab-frame fluid flux.
 */
GKYL_CU_DH static inline void
pkpm_harm_fluid_moment_map(
  const double *mom, double rho, const double *u, double xi, const double *b, double *flux
)
{
  double I0 = mom[0], I1 = mom[1], I2 = mom[2], J0 = mom[3];
  const double *IM = &mom[4];
  double bu = b[0] * u[0] + b[1] * u[1] + b[2] * u[2], uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  double uIM = u[0] * IM[0] + u[1] * IM[1] + u[2] * IM[2];
  flux[0] = rho * I0;
  for (int k = 0; k < 3; ++k) {
    flux[1 + k] = rho * (u[k] * I0 + xi * b[k] * I1) + rho * xi * IM[k];
  }
  flux[4] = 0.5 * rho * (uu * I0 + 2.0 * xi * bu * I1 + xi * xi * I2) + rho * xi * xi * J0 + rho * xi * uIM;
}

// ---- The w flux at a node (gkyl_pkpm_harm_slots.h): nv holds the raw w-flux coefficients at the
// node, theta, b, alpha_0 and c there (PKPM_HS_NODE_*). These are the
// single definition used by the w-face kernels and the closure kernels (plan 17.15).

// alpha = alpha_0 + (c - K) w - (b.grad xi) w^2 at a node.
GKYL_CU_DH static inline double
pkpm_harm_node_alpha(const double *nv, double w)
{
  return nv[PKPM_HS_NODE_A0] + (nv[PKPM_HS_NODE_CH] + nv[PKPM_HS_KW]) * w + nv[PKPM_HS_BG] * w * w;
}

// The known part of alpha, -K w - (b.grad xi) w^2, at a node.
GKYL_CU_DH static inline double
pkpm_harm_node_alpha_k(const double *nv, double w)
{
  return nv[PKPM_HS_KW] * w + nv[PKPM_HS_BG] * w * w;
}

// Coupling part T(U) of the w flux at a node: F <- t.M0 + beta G, G <- t.N + beta H,
// M0 <- G t + R_B1 N, N <- H t + R_B1 O.
GKYL_CU_DH static inline void
pkpm_harm_node_coupling(const double *nv, double w, const double *U, double *T)
{
  double F = U[0], G = U[1];
  const double *M = &U[2], *N = &U[5];
  double tF = nv[PKPM_HS_BF] * G, tG = nv[PKPM_HS_BGG] * G + nv[PKPM_HS_BGF] * F;
  for (int a = 0; a < 3; ++a) {
    tF += (nv[PKPM_HS_TF(a, 0)] + w * nv[PKPM_HS_TF(a, 1)]) * M[a];
    tG += (nv[PKPM_HS_TG(a, 0)] + w * nv[PKPM_HS_TG(a, 1)]) * N[a];
  }
  T[0] = tF;
  T[1] = tG;
  for (int a = 0; a < 3; ++a) {
    double tM = (nv[PKPM_HS_TM(a, 0)] + w * nv[PKPM_HS_TM(a, 1)]) * G;
    double tN = (nv[PKPM_HS_TNG(a, 0)] + w * nv[PKPM_HS_TNG(a, 1)]) * G +
      (nv[PKPM_HS_TNF(a, 0)] + w * nv[PKPM_HS_TNF(a, 1)]) * F;
    for (int c = 0; c < 3; ++c) {
      tM += nv[PKPM_HS_RMN(a, c)] * N[c];
      tN += nv[PKPM_HS_RNN(a, c)] * N[c] + nv[PKPM_HS_RNM(a, c)] * M[c];
    }
    T[2 + a] = tM;
    T[5 + a] = tN;
  }
}

/**
 * The w-independent part of the bounds of the spectrum of the coupling at a node (as
 * pkpm_harm_speed_bounds, plan 17.11): the candidates cand = {beta theta (2 +- sqrt 2),
 * r theta (3 +- sqrt 3)}, r the eigenvalues of R_B1 on the plane perpendicular to b, give
 * nv[PKPM_HS_NODE_SLO] = min(cand) and nv[PKPM_HS_NODE_SHI] = max(cand), and
 * nv[PKPM_HS_NODE_ST] = sqrt((2 + sqrt 2) theta). Called once per node after the w-flux values,
 * theta and b are set.
 */
GKYL_CU_DH static inline void
pkpm_harm_node_spec(double *nv)
{
  double th = fmax(nv[PKPM_HS_NODE_THETA], 0.0), beta = nv[PKPM_HS_BF];
  // M = P R_B1 P on the plane perpendicular to the unit direction of b (its third eigenvalue is
  // 0), P = I - bh bh: tr M = tr(R P) and tr M^2 = tr((R P)^2) (P is a projector), R P = R - (R bh) bh.
  const double *b = &nv[PKPM_HS_NODE_B];
  double bb = pkpm_harm_dot(b, b), ib = bb > 0.0 ? 1.0 / sqrt(bb) : 0.0, bh[3], Rb[3], RP[3][3];
  for (int a = 0; a < 3; ++a) {
    bh[a] = b[a] * ib;
  }
  for (int a = 0; a < 3; ++a) {
    Rb[a] = nv[PKPM_HS_RMN(a, 0)] * bh[0] + nv[PKPM_HS_RMN(a, 1)] * bh[1] + nv[PKPM_HS_RMN(a, 2)] * bh[2];
  }
  double trM = 0.0, trM2 = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int c = 0; c < 3; ++c) {
      RP[a][c] = nv[PKPM_HS_RMN(a, c)] - Rb[a] * bh[c];
    }
    trM += RP[a][a];
  }
  for (int a = 0; a < 3; ++a) {
    for (int c = 0; c < 3; ++c) {
      trM2 += RP[a][c] * RP[c][a];
    }
  }
  double disc = sqrt(fmax(2.0 * trM2 - trM * trM, 0.0));
  double r_lo = 0.5 * (trM - disc), r_hi = 0.5 * (trM + disc);
  const double s2m = 0.5857864376269049, s2p = 3.414213562373095; // 2 -+ sqrt 2
  const double s3m = 1.2679491924311228, s3p = 4.732050807568877; // 3 -+ sqrt 3
  double cand[6] = {beta * th * s2m, beta * th * s2p, r_lo * th * s3m, r_lo * th * s3p, r_hi * th * s3m,
    r_hi * th * s3p};
  double cmin = cand[0], cmax = cand[0];
  for (int k = 1; k < 6; ++k) {
    cmin = fmin(cmin, cand[k]);
    cmax = fmax(cmax, cand[k]);
  }
  nv[PKPM_HS_NODE_SLO] = cmin;
  nv[PKPM_HS_NODE_SHI] = cmax;
  nv[PKPM_HS_NODE_ST] = sqrt(s2p * th);
}

/**
 * Bounds of the spectrum of the coupling at a node at w (pkpm_harm_node_spec must have been called):
 *   [min(cand) - |t| sqrt((2 + sqrt 2) theta), max(cand) + |t| sqrt((2 + sqrt 2) theta)],
 * t = t0 + w t1.
 */
GKYL_CU_DH static inline void
pkpm_harm_node_bounds(const double *nv, double w, double *lo, double *hi)
{
  double t[3];
  for (int a = 0; a < 3; ++a) {
    t[a] = nv[PKPM_HS_TF(a, 0)] + w * nv[PKPM_HS_TF(a, 1)];
  }
  double ts = sqrt(pkpm_harm_dot(t, t)) * nv[PKPM_HS_NODE_ST];
  *lo = nv[PKPM_HS_NODE_SLO] - ts;
  *hi = nv[PKPM_HS_NODE_SHI] + ts;
}

/**
 * The part of a bound over |w| <= W of the w-flux speeds at a node that does not involve the
 * closures (pkpm_harm_node_spec must have been called). The speeds are alpha(w) + s, s in the
 * spectrum bounds of the coupling, alpha(w) = alpha_0 + (c - K) w - (b.grad xi) w^2 and t = t0 + w t1;
 * with absolute coefficients evaluated at W no |w| <= W exceeds
 *   |alpha_0| + |c - K| W + |b.grad xi| W^2 + max(|SLO|, |SHI|) + (|t0| + |t1| W) sqrt((2 + sqrt 2) theta),
 * of which this returns everything after the first two terms (the column time-step rate of the
 * velocity-space terms, plan 17.27).
 */
GKYL_CU_DH static inline double
pkpm_harm_node_speed_rest(const double *nv, double W)
{
  double t0[3], t1[3];
  for (int a = 0; a < 3; ++a) {
    t0[a] = nv[PKPM_HS_TF(a, 0)];
    t1[a] = nv[PKPM_HS_TF(a, 1)];
  }
  return fabs(nv[PKPM_HS_BG]) * W * W + fmax(fabs(nv[PKPM_HS_NODE_SLO]), fabs(nv[PKPM_HS_NODE_SHI]))
    + (sqrt(pkpm_harm_dot(t0, t0)) + W * sqrt(pkpm_harm_dot(t1, t1))) * nv[PKPM_HS_NODE_ST];
}

/**
 * HLL flux of the eight blocks at a w face, at a node: the flux alpha U + T(U) with the speeds
 * S_L, S_R = alpha + the bounds of T's spectrum, exact upwinding when S_L >= 0 (left state) or
 * S_R <= 0 (right state), the HLL average otherwise. T is linear, so the HLL average is
 * (alpha V + T(V) + S_L S_R (U_r - U_l))/(S_R - S_L) with V = S_R U_l - S_L U_r: one coupling per
 * node in every branch. Returns max(|S_L|, |S_R|).
 */
GKYL_CU_DH static inline double
pkpm_harm_node_wflux(const double *nv, double w, const double *U_l, const double *U_r, double *Fhat)
{
  double alpha = pkpm_harm_node_alpha(nv, w), lo, hi;
  pkpm_harm_node_bounds(nv, w, &lo, &hi);
  double sl = alpha + lo, sr = alpha + hi;
  bool hll = sl < 0.0 && sr > 0.0;
  double V[8], T[8];
  const double *X = sl >= 0.0 ? U_l : U_r;
  if (hll) {
    for (int c = 0; c < 8; ++c) {
      V[c] = sr * U_l[c] - sl * U_r[c];
    }
    X = V;
  }
  pkpm_harm_node_coupling(nv, w, X, T);
  if (hll) {
    double ids = 1.0 / (sr - sl), slsr = sl * sr;
    for (int c = 0; c < 8; ++c) {
      Fhat[c] = (alpha * V[c] + T[c] + slsr * (U_r[c] - U_l[c])) * ids;
    }
  } else {
    for (int c = 0; c < 8; ++c) {
      Fhat[c] = alpha * X[c] + T[c];
    }
  }
  return fmax(fabs(sl), fabs(sr));
}
