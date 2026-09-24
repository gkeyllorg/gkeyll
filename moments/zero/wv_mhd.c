#include <math.h>
#include <float.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_mhd.h>
#include <gkyl_wv_mhd.h>

// Make indexing cleaner and clearer
#define DN (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)
#define BX (5)
#define BY (6)
#define BZ (7)
#define PSI_GLM (8)

#define sq(x) ((x) * (x))

static inline void
cons_to_riem_8(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<8; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons_8(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<8; ++i)
    qout[i] = win[i];
}

static inline void
cons_to_riem_9(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<9; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons_9(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<9; ++i)
    qout[i] = win[i];
}

static inline void
rot_to_local_rect(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  // Mass density is a scalar
  qlocal[0] = qglobal[0];
  // Rotate momentum to local coordinates
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
  // Total energy is a scalar
  qlocal[4] = qglobal[4];
  // Rotate B to local coordinates
  qlocal[5] = qglobal[5]*norm[0] + qglobal[6]*norm[1] + qglobal[7]*norm[2];
  qlocal[6] = qglobal[5]*tau1[0] + qglobal[6]*tau1[1] + qglobal[7]*tau1[2];
  qlocal[7] = qglobal[5]*tau2[0] + qglobal[6]*tau2[1] + qglobal[7]*tau2[2];
}

static inline void
rot_to_global_rect(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  // Mass density is a scalar
  qglobal[0] = qlocal[0];
  // Rotate momentum back to global coordinates
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
  // Total energy is a scalar
  qglobal[4] = qlocal[4];
  // Rotate B back to global coordinates
  qglobal[5] = qlocal[5]*norm[0] + qlocal[6]*tau1[0] + qlocal[7]*tau2[0];
  qglobal[6] = qlocal[5]*norm[1] + qlocal[6]*tau1[1] + qlocal[7]*tau2[1];
  qglobal[7] = qlocal[5]*norm[2] + qlocal[6]*tau1[2] + qlocal[7]*tau2[2];
}

static inline void
rot_to_local_rect_glm(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  rot_to_local_rect(eqn, tau1, tau2, norm, qglobal, qlocal);
  qlocal[8] = qglobal[8];
}

static inline void
rot_to_global_rect_glm(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  rot_to_global_rect(eqn, tau1, tau2, norm, qlocal, qglobal);
  qglobal[8] = qlocal[8];
}

struct wv_mhd {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_mhd_div_constraint divergence_constraint; // divB correction
  double glm_ch; // factor to use in GLM scheme
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
};

static void
mhd_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_mhd *mhd = container_of(base, struct wv_mhd, eqn);
  gkyl_free(mhd);
}

// Computing waves and waves speeds from Roe linearization.
// Following Cargo & Gallice 1997 section 4.2.
// FIXME: is CG97 linearization consistent with the flux difference when the jump
// in Bx is nonzero?
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *dQ, const double *ql, const double *qr, double *waves, double *ev)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gamma = mhd->gas_gamma;

  //////////////////////////////////////////////////////////////////////////////
  // STEP 1: COMPUTE PRIMITIVE VARIABLES                                      //
  //////////////////////////////////////////////////////////////////////////////
  double ul = ql[MX] / ql[DN], ur = qr[MX] / qr[DN];
  double vl = ql[MY] / ql[DN], vr = qr[MY] / qr[DN];
  double wl = ql[MZ] / ql[DN], wr = qr[MZ] / qr[DN];
  double pl = gkyl_mhd_pressure(gamma, ql);
  double pr = gkyl_mhd_pressure(gamma, qr);
  double pbl = 0.5 * (sq(ql[BX]) + sq(ql[BY]) + sq(ql[BZ]));
  double pbr = 0.5 * (sq(qr[BX]) + sq(qr[BY]) + sq(qr[BZ]));
  // total enthalpy (density) in CG97 eq. 2.2
  double Hl = (ql[ER] + pl + pbl) / ql[DN];
  double Hr = (qr[ER] + pr + pbr) / qr[DN];

  //////////////////////////////////////////////////////////////////////////////
  // STEP 2: COMPUTE ROE AVERAGES OF PRIMITIVE VARIABLES                      //
  //////////////////////////////////////////////////////////////////////////////
  double srrhol = sqrt(ql[DN]);
  double srrhor = sqrt(qr[DN]);
  double sl = srrhol / (srrhol + srrhor);
  double sr = srrhor / (srrhol + srrhor);

  double rho = srrhol * srrhor;
  double u = sl * ul + sr * ur;
  double v = sl * vl + sr * vr;
  double w = sl * wl + sr * wr;
  double H = sl * Hl + sr * Hr;  // total enthalpy
  double Bx = sr * ql[BX] + sl * qr[BX];
  double By = sr * ql[BY] + sl * qr[BY];
  double Bz = sr * ql[BZ] + sl * qr[BZ];

  //////////////////////////////////////////////////////////////////////////////
  // STEP 3: COMPUTE CHARACTERASTIC WAVE SPEEDS AND OTHER USEFUL QUANTITIES   //
  //////////////////////////////////////////////////////////////////////////////
  // CG97 eq. 4.12, including jump in Bx seems to give correct jump in pressure
  // according to the equation bewteen CG97 eq. 4.15 and 4.16; X may also be
  // computed from CG97 eq. 4.15
  double X = (sq(dQ[BX]) + sq(dQ[BY]) + sq(dQ[BZ])) / (2*sq(srrhol+srrhor));

  // CG97 eq 4.17, wave speeds
  double ca2 = Bx*Bx/rho; // for alfven speed due to normal B field
  double b2 = (Bx*Bx+By*By+Bz*Bz) / rho; // for alfven speed due to full B field
  double v2 = u*u+v*v+w*w;
  double Hgas = H - b2;  // enthalpy of the gas
  double a2 = (2-gamma)*X + (gamma-1)*(Hgas-0.5*v2);  // for sound speed

  double astar2 = a2 + b2;
  double cf2 = (astar2 + sqrt(sq(astar2)-4*a2*ca2)) / 2;  // fast wave speed
  double cs2 = (astar2 - sqrt(sq(astar2)-4*a2*ca2)) / 2;  // slow wave speed

  double a = sqrt(a2);  // sound speed
  double ca = sqrt(ca2);  // alfven speed due to normal B field
  double cf = sqrt(cf2);  // fast magnetosonic speed
  double cs = sqrt(cs2);  // slow magnetosonic speed

  // CG97 eq 4.21, S, beta, and alpha
  double S = Bx >=0? 1: -1;

  double Bt = sqrt(By*By+Bz*Bz);
  double betay, betaz;
  if (Bt>0) { // TODO compare wtih a tiny number
    betay = By / Bt;
    betaz = Bz / Bt;
  } else {
    betay = betaz = 1/sqrt(2);
  }

  double alphaf2 = (a2 - cs2) / (cf2 - cs2);
  double alphas2 = (cf2 - a2) / (cf2 - cs2);
  double alphaf = sqrt(alphaf2);
  double alphas = sqrt(alphas2);

  //////////////////////////////////////////////////////////////////////////////
  // STEP 4: COMPUTE WAVE SPEEDS AND FLUCTUATIONS DUE TO ROE LINEARIZATION    //
  //                                                                          //
  // For each eigensolution of the Roe linearizing matrix, compute:           //
  // 0. Eigenvalue, which is the wave speed. The results are stored in ev     //
  //     array and sorted by their values, that is,                           //
  //    u-c_fast<u-c_alfven<u-c_slow<u<u+c_slow<u+c_alfven<u+c_fast           //
  // 1. Projection coefficient of jumps in conservative state on to the       //
  //    left-eigenvector according to CG97 eq. 4.20.                          //
  // 2. The wave (fluctuations) due to this mode.                             //
  //////////////////////////////////////////////////////////////////////////////
  // CG97 eq. 4.20 gave projection coefficients in primitive varialbles.
  // Alternatively, one may compute the full left-eigenvectors and dot it with
  // dQ. The left-eigenvector of the Roe matrix can be found in Stone et al.
  // (2008), Appendix B2.
  double drho = dQ[DN];
  double du = ur - ul;
  double dv = vr - vl;
  double dw = wr - wl;
  double dBy = dQ[BY];
  double dBz = dQ[BZ];
  double dp = pr - pl; // consistent the CG97 eq. between 4.15 and 4.16

  const int meqns = eqn->num_equations;
  const int mwaves = eqn->num_waves;
  for (int i=0; i<meqns*mwaves; ++i) waves[i] = 0.0;
  double *wv;
  double eta[mwaves];

  /////////////////////////////
  // Fast magnetosonic waves //
  /////////////////////////////

  // For projection coefficients in eq. 4.20 for magnetosonic waves
  double t1 = alphaf * (X*drho + dp);
  double t2 = alphas * cs * rho * S * (betay*dv + betaz*dw);
  double t3 = alphaf * cf * rho * du;
  double t4 = alphas * sqrt(rho) * a * (betay*dBy + betaz*dBz);
  // For right eigenvectors for magnetosonic waves in eq. 4.19
  double t5 = alphas * cs * S * (v*betay + w*betaz);
  double t6 = alphas * a * Bt / sqrt(rho);

  ev[0] = u - cf;
  // Coefficients of projection onto the left-eigenvector
  eta[0] = (t1 + t2 - t3 + t4) / 2; // eq. 4.20
  eta[0] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[0*meqns];
  wv[DN] = eta[0]*alphaf;
  wv[MX] = eta[0]*(alphaf*(u - cf));
  wv[MY] = eta[0]*(alphaf*v + alphas*cs*betay*S);
  wv[MZ] = eta[0]*(alphaf*w + alphas*cs*betaz*S);
  wv[ER] = eta[0]*(alphaf*(Hgas-u*cf)+t5+t6); // -t6 in CG97 typo? +t6 in SG08
  wv[BY] = eta[0]*alphas*a*betay/sqrt(rho);
  wv[BZ] = eta[0]*alphas*a*betaz/sqrt(rho);

  ev[6] = u + cf;
  // Coefficients of projection onto the left-eigenvector
  eta[6] = (t1 - t2 + t3 + t4) / 2; // eq. 4.20
  eta[6] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[6*meqns];
  wv[DN] = eta[6]*alphaf;
  wv[MX] = eta[6]*(alphaf*(u + cf));
  wv[MY] = eta[6]*(alphaf*v - alphas*cs*betay*S);
  wv[MZ] = eta[6]*(alphaf*w - alphas*cs*betaz*S);
  wv[ER] = eta[6]*(alphaf*(Hgas+u*cf)-t5+t6); // -t6 in CG97 typo? +t6 in SG08
  wv[BY] = eta[6]*alphas*a*betay/sqrt(rho);
  wv[BZ] = eta[6]*alphas*a*betaz/sqrt(rho);

  //////////////////
  // Alfven waves //
  //////////////////

  // For projection coefficients in eq. 4.20 for alfven waves
  double t7 = betay*dw - betaz*dv;
  double t8 = (S/sqrt(rho)) * (betay*dBz-betaz*dBy);

  ev[1] = u - ca;
  // Coefficients of projection onto the left-eigenvector
  eta[1] = (t7 + t8) / 2; // eg. 4.20
  eta[1] *= rho; // merge rho in right-eigenvector in 4.18 into coefficents
  wv = &waves[1*meqns];
  // Fluctuations due to this mode
  wv[MY] = -eta[1]*betaz;
  wv[MZ] = eta[1]*betay;
  wv[ER] = -eta[1]*(v*betaz - w*betay);
  wv[BY] = -eta[1]*S*betaz/sqrt(rho);
  wv[BZ] = eta[1]*S*betay/sqrt(rho);

  ev[5] = u + ca;
  // Coefficients of projection onto the left-eigenvector
  eta[5] = (-t7 + t8) / 2; // eq. 4.20
  eta[5] *= rho; // merge rho in right-eigenvector in 4.18 into coefficents
  wv = &waves[5*meqns];
  // Fluctuations due to this mode
  wv[MY] = eta[5]*betaz;
  wv[MZ] = -eta[5]*betay;
  wv[ER] = +eta[5]*(v*betaz - w*betay);
  wv[BY] = -eta[5]*S*betaz/sqrt(rho);
  wv[BZ] = eta[5]*S*betay/sqrt(rho);

  /////////////////////////////
  // Slow magnetosonic waves //
  /////////////////////////////

  // For projection coefficients in eq. 4.20 for magnetosonic waves
  // Compared to the fast wave coefficients, alphas/alphaf & cs/cf are switched
  t1 = alphas * (X*drho + dp);
  t2 = alphaf * cf * rho * S * (betay*dv + betaz*dw);
  t3 = alphas * cs * rho * du;
  t4 = alphaf * sqrt(rho) * a * (betay*dBy + betaz*dBz);
  // for right eigenvectors for magnetosonic waves in eq. 4.19
  t5 = alphaf * cf * S * (v*betay + w*betaz);
  t6 = alphaf * a * Bt / sqrt(rho);

  ev[2] = u - cs;
  // Coefficients of projection onto the left-eigenvector
  eta[2] = (t1 - t2 - t3 - t4) / 2; // eq. 4.20
  eta[2] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[2*meqns];
  wv[DN] = eta[2]*alphas;
  wv[MX] = eta[2]*(alphas*(u - cs));
  wv[MY] = eta[2]*(alphas*v - alphaf*cf*betay*S);
  wv[MZ] = eta[2]*(alphas*w - alphaf*cf*betaz*S);
  wv[ER] = eta[2]*(alphas*(Hgas-u*cs) - t5 - t6);
  wv[BY] = -eta[2]*alphaf*a*betay/sqrt(rho);
  wv[BZ] = -eta[2]*alphaf*a*betaz/sqrt(rho);

  ev[4] = u + cs;
  // Coefficients of projection onto the left-eigenvector
  eta[4] = (t1 + t2 + t3 - t4) / 2; // eq. 4.20
  eta[4] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[4*meqns];
  wv[DN] = eta[4]*alphas;
  wv[MX] = eta[4]*(alphas*(u + cs));
  wv[MY] = eta[4]*(alphas*v + alphaf*cf*betay*S);
  wv[MZ] = eta[4]*(alphas*w + alphaf*cf*betaz*S);
  wv[ER] = eta[4]*(alphas*(Hgas+u*cs) + t5 - t6);
  wv[BY] = -eta[4]*alphaf*a*betay/sqrt(rho);
  wv[BZ] = -eta[4]*alphaf*a*betaz/sqrt(rho);

  //////////////////
  // Enropy wave  //
  //////////////////
  ev[3] = u;
  // Coefficients of projection onto the left-eigenvector
  eta[3] = (a2 - X) * drho - dp; // eq. 4.20
  eta[3] /= a2; // merge a2 in right-eigenvector in 4.18 into coefficents
  // Fluctuations due to this mode
  wv = &waves[3*meqns];
  wv[DN] = eta[3];
  wv[MX] = eta[3]*u;
  wv[MY] = eta[3]*v;
  wv[MZ] = eta[3]*w;
  wv[ER] = eta[3]*(v2/2 + X*(gamma-2)/(gamma-1));

  double max_speed =  fabs(u) + cf;

  ////////////////////
  // div(B) wave(s) //
  ////////////////////

  // For the eight-wave scheme, advect the jump in Bx at the speed u.
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_EIGHT_WAVES) {
    ev[7] = u;
    wv = &waves[7*meqns];
    wv[BX] = dQ[BX];
  }

  // For the GLM Bx and psi waves, solve the linear Riemann problem.
  // TODO create and use a separate glm RP solver
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_GLM)
  {
    double ch = mhd->glm_ch;

    // L = 0.5*(-ch, 1), R = (-1/ch, 1)
    ev[7] = -ch;
    eta[7] = 0.5 * (-dQ[BX]*ch+dQ[PSI_GLM]);
    wv = &waves[7*meqns];
    wv[BX] = -eta[7]/ch;
    wv[PSI_GLM] = eta[7];

    // L = 0.5*(+ch, 1), R = (+1/ch, 1)
    ev[8] = ch;
    eta[8] = 0.5 * (dQ[BX]*ch+dQ[PSI_GLM]);
    wv = &waves[8*meqns];
    wv[BX] = eta[8]/ch;
    wv[PSI_GLM] = eta[8];

    max_speed = max_speed > ch ? max_speed : ch;
  }

  return max_speed;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s, double *amdq,
  double *apdq)
{
  int meqn = eqn->num_equations;
  for (int i=0; i<meqn; ++i) {
    amdq[i] = fmin(0.0, s[0]) * waves[i];
    apdq[i] = fmax(0.0, s[0]) * waves[i];
    for (int mw=1; mw<eqn->num_waves; ++mw) {
      amdq[i] += fmin(0.0, s[mw]) * waves[meqn*mw+i];
      apdq[i] += fmax(0.0, s[mw]) * waves[meqn*mw+i];
    }
  }
}

// Computing waves and waves speeds from Lax fluxes
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *dQ, const double *ql, const double *qr, double *waves, double *ev)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  int meqn = eqn->num_equations;
  double gas_gamma = mhd->gas_gamma;

#if 0
  double sl = gkyl_mhd_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_mhd_max_abs_speed(gas_gamma, qr);
  ev[0] = 0.5*(sl+sr);
#else
  double amax = gkyl_mhd_max_abs_speed_roe(gas_gamma, ql, qr);
  ev[0] = -amax;
  ev[1] = amax;
#endif

  double fl[10], fr[10];
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_GLM) {
    double ch = mhd->glm_ch;
    gkyl_glm_mhd_flux(gas_gamma, ch, ql, fl);
    gkyl_glm_mhd_flux(gas_gamma, ch, qr, fr);
  } else {
    gkyl_mhd_flux(gas_gamma, ql, fl);
    gkyl_mhd_flux(gas_gamma, qr, fr);
  }

  double *w0 = &waves[0], *w1 = &waves[meqn];
  for (int i=0; i<meqn; ++i) {
    w0[i] = 0.5*((qr[i]-ql[i]) - (fr[i]-fl[i])/amax);
    w1[i] = 0.5*((qr[i]-ql[i]) + (fr[i]-fl[i])/amax);
  }

  return ev[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = eqn->num_equations;
  const double *w0 = &waves[0], *w1 = &waves[meqn];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<meqn; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

// HLLD, Miyoshi & Kusano (2005), JCP, 208(1), 315-344
static double
wave_hlld(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double g = mhd->gas_gamma;
  int meqn = eqn->num_equations;

  // notations based on Miyoshi & Kusano: r:rho, u:ux, p:pressure, l:left,
  // r:right, s:star, m:middle/intermidiate
  double rl = ql[DN];
  double ul = ql[MX] / rl;
  double vl = ql[MY] / rl;
  double wl = ql[MZ] / rl;
  double pl = gkyl_mhd_pressure(g, ql);
  double ptl = pl + 0.5 * (sq(ql[BX]) + sq(ql[BY]) + sq(ql[BZ]));

  double rr = qr[DN];
  double ur = qr[MX] / rr;
  double vr = qr[MY] / rr;
  double wr = qr[MZ] / rr;
  double pr = gkyl_mhd_pressure(g, qr);
  double ptr = pr + 0.5 * (sq(qr[BX]) + sq(qr[BY]) + sq(qr[BZ]));

  // STEP 1. compute min and max wave speeds
  double cf_l = gkyl_mhd_fast_speed(g, ql);
  double smin_l = ul - cf_l;
  double smax_l = ul + cf_l;
  double cf_r = gkyl_mhd_fast_speed(g, qr);
  double smin_r = ur - cf_r;
  double smax_r = ur + cf_r;
#if 0
  // Estimation eq. 12 by Davis SIAM J. Sci. Statist. Comput., 9 (1988), p. 445
  double sl = smin_l < smin_r ? smin_l : smin_r;
  double sr = smax_l > smax_r ? smax_l : smax_r;
#else
  // Estimation eq. 13 by Einfeldt et al. J. Comput. Phys., 92 (1991), p. 273
  double buf[4];
  gkyl_mhd_eigen_speeds_roe(g, ql, qr, buf);
  double u_roe = buf[0], cf_roe = buf[3];
  double smin_roe = u_roe - cf_roe;
  double smax_roe = u_roe + cf_roe;
  double sl = smin_l < smin_roe ? smin_l : smin_roe;
  double sr = smax_r > smax_roe ? smax_r : smax_roe;
#endif

  // FIXME Miyoshi & Kusano did not specify Bx
  double Bx = (sr*qr[BX] - sl*ql[BX]) / (sr - sl);
  double sign = Bx > 0? 1 : -1;

  // STEP 2. compute intermediate wave speeds
  // middle wave speed,eq. 39
  double tmp = 1 / (rr*(sr-ur) - rl*(sl-ul));
  double sm = (rr*ur*(sr-ur) - rl*ul*(sl-ul) - ptr+ptl) * tmp;
  // eq. 41, p^*_T
  double pt = ((sr-ur)*rr*ptl-(sl-ul)*rl*ptr+rl*rr*(sr-ur)*(sl-ul)*(ur-ul))*tmp;

  // sl: outer left, ssl: inner left, ssr: inner right, sr: outer right
  double rsl = rl * (sl-ul) / (sl-sm); // eq. 43
  double rsr = rr * (sr-ur) / (sr-sm); // eq. 43
  double sqrtl = sqrt(rsl);
  double sqrtr = sqrt(rsr);
  double ssl = sm - Bx*sign / sqrtl; // eq. 51
  double ssr = sm + Bx*sign / sqrtr; // eq. 51

  // STEP 3. compute intermediate states
  // outer left, inner left, inner right, outer right; s: star, ss: two star
  double qsl[8], qssl[8], qssr[8], qsr[8];
  double tmp1, tmp2, tmp3; // convenience temporary variables

  // left and right outer intermediate states
  tmp = 1 / (rl*(sl-ul)*(sl-sm)-Bx*Bx);
  tmp1 = Bx * (sm-ul) * tmp;
  tmp2 = (rl*sq(sl-ul)-sq(Bx)) * tmp;
  double usl = sm; // eq. 39
  double vsl = vl - ql[BY] * tmp1; // eq. 44
  double wsl = wl - ql[BZ] * tmp1; // eq. 46
  qsl[DN] = rsl;
  qsl[MX] = qsl[DN] * usl;
  qsl[MY] = qsl[DN] * vsl;
  qsl[MZ] = qsl[DN] * wsl;
  qsl[BX] = Bx; // FIXME
  qsl[BY] = ql[BY] * tmp2; // eq. 45
  qsl[BZ] = ql[BZ] * tmp2; // eq. 47
  tmp3 = ul*ql[BX]+vl*ql[BY]+wl*ql[BZ] - (usl*qsl[BX]+vsl*qsl[BY]+wsl*qsl[BZ]);
  qsl[ER] = ((sl-ul)*ql[ER] -ptl*ul + pt*sm + Bx*tmp3) / (sl-sm); // eq. 48

  tmp = 1 / (rr*(sr-ur)*(sr-sm)-Bx*Bx);
  tmp1 = Bx * (sm-ur) * tmp;
  tmp2 = (rr*sq(sr-ur)-sq(Bx)) * tmp;
  double usr = sm; // eq. 39
  double vsr = vr - qr[BY] * tmp1; // eq. 44
  double wsr = wr - qr[BZ] * tmp1; // eq. 46
  qsr[DN] = rsr ;
  qsr[MX] = qsr[DN] * usr;
  qsr[MY] = qsr[DN] * vsr;
  qsr[MZ] = qsr[DN] * wsr;
  qsr[BX] = Bx; // FIXME
  qsr[BY] = qr[BY] * tmp2; // eq. 45
  qsr[BZ] = qr[BZ] * tmp2; // eq. 47
  tmp3 = ur*qr[BX]+vr*qr[BY]+wr*qr[BZ] - (usr*qsr[BX]+vsr*qsr[BY]+wsr*qsr[BZ]);
  qsr[ER] = ((sr-ur)*qr[ER] - ptr*ur + pt*sm + Bx*tmp3) / (sr-sm); // eq. 48

  // left and right inner intermediate states
  tmp = 1 / (sqrtl + sqrtr);
  double uss = sm; // eq. 39
  // eq. 59, 60
  double vss = (sqrtl*vsl + sqrtr*vsr + (qsr[BY] - qsl[BY]) * sign) * tmp;
  double wss = (sqrtl*wsl + sqrtr*wsr + (qsr[BZ] - qsl[BZ]) * sign) * tmp;
  // eq. 61, 62
  tmp1 = sqrtl*sqrtr*sign;
  double Byss = (sqrtl*qsr[BY] + sqrtr*qsl[BY] + (vsr-vsl)*tmp1) * tmp;
  double Bzss = (sqrtl*qsr[BZ] + sqrtr*qsl[BZ] + (wsr-wsl)*tmp1) * tmp;
  tmp2 = uss*Bx + vss*Byss + wss*Bzss;

  qssl[DN] = qsl[DN]; // eq. 49
  qssl[MX] = qssl[DN] * uss; // eq. 39
  qssl[MY] = qssl[DN] * vss; // eq. 55
  qssl[MZ] = qssl[DN] * wss; // eq. 56
  qssl[BX] = Bx; // FIXME
  qssl[BY] = Byss; // eq. 56
  qssl[BZ] = Bzss; // eq. 56
  tmp3 = usl*qsl[BX]+vsl*qsl[BY]+wsl*qsl[BZ] - tmp2;
  qssl[ER] = qsl[ER] - sqrtl*tmp3*sign; // eq. 63

  qssr[DN] = qsr[DN]; // eq. 49
  qssr[MX] = qssr[DN] * uss; // eq. 39
  qssr[MY] = qssr[DN] * vss; // eq. 55
  qssr[MZ] = qssr[DN] * wss; // eq. 56
  qssr[BX] = Bx; // FIXME
  qssr[BY] = Byss; // eq. 56
  qssr[BZ] = Bzss; // eq. 56
  tmp3 = usr*qsr[BX]+vsr*qsr[BY]+wsr*qsr[BZ] - tmp2;
  qssr[ER] = qsr[ER] + sqrtr*tmp3*sign; // eq. 63

  // STEP 4. collect all waves and wave speeds
  speeds[0] = sl;
  speeds[1] = ssl;
  speeds[2] = sm;
  speeds[3] = ssr;
  speeds[4] = sr;

  double *wv;

  wv = waves;
  for (int i=0; i<8; ++i)  wv[i] = qsl[i] - ql[i];

  wv += meqn;
  for (int i=0; i<8; ++i)  wv[i] = qssl[i] - qsl[i];

  wv += meqn;
  for (int i=0; i<8; ++i)  wv[i] = qssr[i] - qssl[i];

  wv += meqn;
  for (int i=0; i<8; ++i)  wv[i] = qsr[i] - qssr[i];

  wv += meqn;
  for (int i=0; i<8; ++i)  wv[i] = qr[i] - qsr[i];

  double max_speed = sr;

  // For the eight-wave scheme, advect the jump in Bx at the speed u.
  // XXX is this correct?
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_EIGHT_WAVES) {
    speeds[5] = sm;
    wv += meqn;
    for (int i=0; i<8; ++i) wv[i] = 0.0;
    wv[BX] = dQ[BX];
  }

  // For the GLM Bx and psi waves, solve the linear Riemann problem.
  // XXX is this correct? TODO create and use a separate glm RP solver
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_GLM)
  {
    for (int w=0; w<5; ++w)
    {
      waves[w*meqn + PSI_GLM] = 0.0;
    }

    double ch = mhd->glm_ch;

    // L = 0.5*(-ch, 1), R = (-1/ch, 1)
    speeds[5] = -ch;
    wv += meqn;
    for (int i=0; i<8; ++i) wv[i] = 0.0;
    double eta = 0.5 * (-dQ[BX]*ch+dQ[PSI_GLM]);
    wv[BX] = -eta/ch;
    wv[PSI_GLM] = eta;

    // L = 0.5*(+ch, 1), R = (+1/ch, 1)
    speeds[6] = ch;
    wv += meqn;
    for (int i=0; i<8; ++i) wv[i] = 0.0;
    eta = 0.5 * (dQ[BX]*ch+dQ[PSI_GLM]);
    wv[BX] = eta/ch;
    wv[PSI_GLM] = eta;

    max_speed = max_speed > ch ? max_speed : ch;
  }

  return max_speed;
}

static void
qfluct_hlld(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
            const double *waves, const double *s, double *amdq, double *apdq)
{
  int meqn = eqn->num_equations;
  int mwave = eqn->num_waves;
  for (int i=0; i<meqn; ++i) {
    amdq[i] = fmin(0.0, s[0]) * waves[i];
    apdq[i] = fmax(0.0, s[0]) * waves[i];
    for (int mw=1; mw<mwave; ++mw) {
      amdq[i] += fmin(0.0, s[mw]) * waves[meqn*mw+i];
      apdq[i] += fmax(0.0, s[mw]) * waves[meqn*mw+i];
    }
  }
}

static double
wave_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_roe(eqn, delta, ql, qr, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_hlld_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_hlld(eqn, delta, ql, qr, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct_hlld_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_hlld(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir, double *waves, double *s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);

  if (q[0] < 0.0)
    return false;

  double pr = gkyl_mhd_pressure(mhd->gas_gamma, q);
  if (pr < 0.0)
    return false;

  return true;
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return gkyl_mhd_max_abs_speed(mhd->gas_gamma, q);
}

static inline void
mhd_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 8; i++) {
    sout[i] = 0.0;
  }
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_new(const struct gkyl_wv_mhd_inp *inp)
{
  struct wv_mhd *mhd = gkyl_malloc(sizeof(struct wv_mhd));

  mhd->eqn.type = GKYL_EQN_MHD;
  mhd->gas_gamma = inp->gas_gamma;

  switch (inp->rp_type) {
    case WV_MHD_RP_ROE:
      mhd->eqn.num_equations = 8;
      mhd->eqn.num_waves = 7;  
      mhd->eqn.waves_func = wave_roe_l;
      mhd->eqn.qfluct_func = qfluct_roe_l;
      break;

    case WV_MHD_RP_HLLD:
      mhd->eqn.num_equations = 8;
      mhd->eqn.num_waves = 5;  
      mhd->eqn.waves_func = wave_hlld_l;
      mhd->eqn.qfluct_func = qfluct_hlld_l;
      break;
      
    case WV_MHD_RP_LAX:
      mhd->eqn.num_equations = 8;
      mhd->eqn.num_waves = 2;
      mhd->eqn.waves_func = wave_lax_l;
      mhd->eqn.qfluct_func = qfluct_lax_l;
      break;      
  }

  mhd->eqn.check_inv_func = check_inv;
  mhd->eqn.max_speed_func = max_speed;
  mhd->eqn.rotate_to_local_func = rot_to_local_rect;
  mhd->eqn.rotate_to_global_func = rot_to_global_rect;

  mhd->eqn.cons_to_riem = cons_to_riem_8;
  mhd->eqn.riem_to_cons = riem_to_cons_8;

  mhd->divergence_constraint = inp->divergence_constraint;
  switch (inp->divergence_constraint) {
    case GKYL_MHD_DIVB_NONE:
      break;

    case GKYL_MHD_DIVB_EIGHT_WAVES:
      if (inp->rp_type != WV_MHD_RP_LAX)
        mhd->eqn.num_waves += 1;
      break;

    case GKYL_MHD_DIVB_GLM:
      mhd->eqn.num_equations += 1;
      if (inp->rp_type != WV_MHD_RP_LAX)
        mhd->eqn.num_waves += 2;
      mhd->eqn.cons_to_riem = cons_to_riem_9;
      mhd->eqn.riem_to_cons = riem_to_cons_9;
      mhd->eqn.rotate_to_local_func = rot_to_local_rect_glm;
      mhd->eqn.rotate_to_global_func = rot_to_global_rect_glm;
      mhd->glm_ch = inp->glm_ch;
      mhd->glm_alpha = inp->glm_alpha;
      break;
  }

  mhd->eqn.num_diag = mhd->eqn.num_equations;
  // probably want to change this to store magnetic, internal and KE 
  mhd->eqn.cons_to_diag = gkyl_default_cons_to_diag;

  mhd->eqn.source_func = mhd_source;

  mhd->eqn.ref_count = gkyl_ref_count_init(mhd_free);

  mhd->eqn.embed_geo = NULL;

  return &mhd->eqn;
}

////////////////////////////////
// member getters and setters //
////////////////////////////////

double
gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->gas_gamma;
}

double
gkyl_wv_mhd_divergence_constraint(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->divergence_constraint;
}

double
gkyl_wv_mhd_glm_ch(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->glm_ch;
}

double
gkyl_wv_mhd_glm_alpha(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->glm_alpha;
}

void
gkyl_wv_mhd_set_glm_ch(struct gkyl_wv_eqn* eqn, const double glm_ch)
{
  struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  mhd->glm_ch = glm_ch;
}

