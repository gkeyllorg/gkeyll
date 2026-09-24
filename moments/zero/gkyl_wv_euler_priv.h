#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_euler {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
};

/**
 * Free Euler eqn object.
 *
 * @param ref Reference counter for Euler eqn
 */
void gkyl_euler_free(const struct gkyl_ref_count *ref);

/**
 * Computes the scalar pressure given the conserved variables.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 */
GKYL_CU_DH
static inline double
gkyl_euler_pressure(double gas_gamma, const double q[5])
{
  return (gas_gamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
}

/**
 * Compute primitive variables given conserved variables.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @param v Primitive variables (output)
 */
GKYL_CU_DH
static inline void
gkyl_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  v[0] = q[0];
  v[1] = q[1]/q[0];
  v[2] = q[2]/q[0];
  v[3] = q[3]/q[0];
  v[4] = gkyl_euler_pressure(gas_gamma, q);
}

/**
 * Compute maximum absolute speed.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
GKYL_CU_DH
static inline double
gkyl_euler_max_abs_speed(double gas_gamma, const double q[5])
{
  double v[5] = {0.0};
  gkyl_euler_prim_vars(gas_gamma, q, v);
  double u2 = sqrt(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  double pr = v[4];
  return fabs(u2) + sqrt(gas_gamma*pr/q[0]);
}

/**
 * Compute flux. Assumes rotation to local coordinate system.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
GKYL_CU_DH
static void
gkyl_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[1]/q[0];

  flux[0] = q[1]; // rho*u
  flux[1] = q[1]*u + pr; // rho*u*u + pr
  flux[2] = q[2]*u; // rho*v*u
  flux[3] = q[3]*u; // rho*w*u
  flux[4] = (q[4]+pr)*u; // (E+p)*u
}

GKYL_CU_DH
static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    wout[i] = qin[i];
}
GKYL_CU_DH
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    qout[i] = win[i];
}

// Euler perfectly reflecting wall
GKYL_CU_DH
static void
euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density and pressure
  ghost[0] = skin[0];
  ghost[4] = skin[4];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = skin[2];
  ghost[3] = skin[3];
}

// Euler no-slip wall
GKYL_CU_DH
static void
euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density and pressure
  ghost[0] = skin[0];
  ghost[4] = skin[4];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
}

GKYL_CU_DH
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
  qlocal[4] = qglobal[4];
}

GKYL_CU_DH
static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
  qglobal[4] = qlocal[4];
}

GKYL_CU_DH
static void
wave_embed_absorb(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = 0.0;
  qphi[1] = 0.0;
  qphi[2] = 0.0;
  qphi[3] = 0.0;
  qphi[4] = 0.0;
}

GKYL_CU_DH
static void
wave_embed_reflect(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = q[0];
  qphi[1] = -q[1];
  qphi[2] = q[2];
  qphi[3] = q[3];
  qphi[4] = q[4];
}

GKYL_CU_DH
static double
wave_embedded(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, const double phil,
  const double phir, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;

  double qphi[5] = {0.0};
  double deltaphi[5] = {0.0};

  double amax = 0.0;

  double fl[5], fr[5] = {0.0};
  if ((phil < 0.0) && (phir > 0.0)) {
    double rhor = qr[0];
    double ur = qr[1]/qr[0];
    double pr = gkyl_euler_pressure(gas_gamma, qr);
    double sr = fabs(ur) + sqrt(gas_gamma*pr/rhor);
    amax = sr;
   
    eqn->embed_geo->embed_func(qr, qphi, deltaphi, eqn->embed_geo->ctx);

    gkyl_euler_flux(gas_gamma, qphi, fl);
    gkyl_euler_flux(gas_gamma, qr, fr);

    double *w0 = &waves[0], *w1 = &waves[5];
    for (int i=0; i<5; ++i) {
      w1[i] = 0.5*((qr[i]-qphi[i]) + (fr[i]-fl[i])/amax);
    }
  }
  
  if ((phir < 0.0) && (phil > 0.0)) {
    double rhol = ql[0];
    double ul = ql[1]/ql[0];
    double pl = gkyl_euler_pressure(gas_gamma, ql);
    double sl = fabs(ul) + sqrt(gas_gamma*pl/rhol);
    amax = sl;
    
    eqn->embed_geo->embed_func(ql, qphi, deltaphi, eqn->embed_geo->ctx);

    gkyl_euler_flux(gas_gamma, ql, fl);
    gkyl_euler_flux(gas_gamma, qphi, fr);
    
    double *w0 = &waves[0], *w1 = &waves[5];
    for (int i=0; i<5; ++i) {
      w0[i] = 0.5*((qphi[i]-ql[i]) - (fr[i]-fl[i])/amax);
    }
  }

  s[0] = -amax;
  s[1] = amax;
  
  return s[1];
}

// Waves and speeds using Lax fluxes
GKYL_CU_DH
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;

  double rhol = ql[0], rhor = qr[0];
  double ul = ql[1]/ql[0], ur = qr[1]/qr[0];
  double pl = gkyl_euler_pressure(gas_gamma, ql), pr = gkyl_euler_pressure(gas_gamma, qr);
  double sl = fabs(ul) + sqrt(gas_gamma*pl/rhol), sr = fabs(ur) + sqrt(gas_gamma*pr/rhor);
  double amax = fmax(sl, sr);

  double fl[5], fr[5];
  gkyl_euler_flux(gas_gamma, ql, fl);
  gkyl_euler_flux(gas_gamma, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[5];
  for (int i=0; i<5; ++i) {
    w0[i] = 0.5*((qr[i]-ql[i]) - (fr[i]-fl[i])/amax);
    w1[i] = 0.5*((qr[i]-ql[i]) + (fr[i]-fl[i])/amax);
  }

  s[0] = -amax;
  s[1] = amax;
  
  return s[1];
}

GKYL_CU_DH
static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

GKYL_CU_DH
static double
wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double phil,
  double phir, double *waves, double *s)
{
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);
}

GKYL_CU_DH
static void
qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, double phil,
  double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// project column vector delta onto the right eigenvectors of the flux Jacobian
// evaluated at avg = {u, v, w, enth}
GKYL_CU_DH
static double
proj_onto_euler_eigvect(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *avg, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;
  double g1 = gas_gamma - 1;

  double u = avg[0], v = avg[1], w = avg[2], enth = avg[3];

  // See http://ammar-hakim.org/sj/euler-eigensystem.html for notation
  // and meaning of these terms
  double q2 = u*u+v*v+w*w;
  double aa2 = g1*(enth-0.5*q2);
  double a = sqrt(aa2);
  double g1a2 = g1/aa2, euv = enth-q2;

  // Compute projections of jump
  double a4 = g1a2*(euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  double a2 = delta[2] - v*delta[0];
  double a3 = delta[3] - w*delta[0];
  double a5 = 0.5*(delta[1] + (a-u)*delta[0] - a*a4)/a;
  double a1 = delta[0] - a4 - a5;

  double *wv;
  // Wave 1: eigenvalue is u-c
  wv = &waves[0];
  wv[0] = a1;
  wv[1] = a1*(u-a);
  wv[2] = a1*v;
  wv[3] = a1*w;
  wv[4] = a1*(enth-u*a);
  s[0] = u-a;

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a4;
  wv[1] = a4*u;
  wv[2] = a4*v + a2;
  wv[3] = a4*w + a3;
  wv[4] = a4*0.5*q2 + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wv = &waves[10];
  wv[0] = a5;
  wv[1] = a5*(u+a);
  wv[2] = a5*v;
  wv[3] = a5*w;
  wv[4] = a5*(enth+u*a);
  s[2] = u+a;
  
  return fabs(u)+a;
}

GKYL_CU_DH
inline static void
roe_avg(const struct gkyl_wv_eqn *eqn, const double *ql,
    const double *qr, double *avg)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;

  double rhol = ql[0], rhor = qr[0];
  double pl = gkyl_euler_pressure(gas_gamma, ql), pr = gkyl_euler_pressure(gas_gamma, qr);

  // Roe averages: see Roe's original 1981 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[1]*ravgl1 + qr[1]*ravgr1)*ravg2;
  double v = (ql[2]*ravgl1 + qr[2]*ravgr1)*ravg2;
  double w = (ql[3]*ravgl1 + qr[3]*ravgr1)*ravg2;
  double enth = ((ql[4]+pl)*ravgl1 + (qr[4]+pr)*ravgr1)*ravg2;

  avg[0] = u; avg[1] = v; avg[2] = w; avg[3] = enth;
}

// Waves and speeds using Roe averaging
GKYL_CU_DH
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  double avg[4];
  roe_avg(eqn, ql, qr, avg);
  return proj_onto_euler_eigvect(eqn, delta, avg, waves, s);
}

GKYL_CU_DH
static void
qfluct_roe(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5], *w2 = &waves[10];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
}

GKYL_CU_DH
static double
wave_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double phil,
  double phir, double *waves, double *s)
{
  // clear waves and wave speeds
  int mwaves = (type == GKYL_WV_HIGH_ORDER_FLUX) ? eqn->num_waves : 2;
  int meqn = eqn->num_equations;
  for (int i=0; i<mwaves; ++i) {
    double *w = &waves[i*meqn];
    for (int j=0; j<meqn; ++j)
      w[j] = 0.0;
    s[i] = 0.0;
  }
  // Check if inside an embedded surface
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else {
    if (type == GKYL_WV_HIGH_ORDER_FLUX)
      return wave_roe(eqn, delta, ql, qr, waves, s);
    else
      return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // can't happen
}

GKYL_CU_DH
static void
qfluct_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, double phil, double phir, const double *waves,
  const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX && (phil > 0.0) && (phir > 0.0))
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// HLL
GKYL_CU_DH
static void
states_hll_common(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, double state[8])
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double g = euler->gas_gamma;

  // r->rho, u->u_x, p->pressure, c->sound speed; l->left, r->right, m->middle
  double rl = ql[0];
  double ul = ql[1] / rl;
  double pl = gkyl_euler_pressure(g, ql);

  double rr = qr[0];
  double ur = qr[1] / rr;
  double pr = gkyl_euler_pressure(g, qr);

  // STEP 1. compute min and max wave speeds; here, the pressure-based speed
  // estimates Toro (10.59) are used
  //   first, estimate middle pressure; using PVRS Toro (10.61) but Toro Fig 9.4
  //   presents an adaptive method to switch between PVRS and TRRS (10.63) and
  //   TSRS (10.65); for ideal gases, max(0, p_PVRS) might also be good
  double cl = sqrt(g*pl/rl);
  double cr = sqrt(g*pr/rr);
  double ra = 0.5 * (rl + rr);
  double ca = 0.5 * (cl + cr);
  double pm = 0.5 * (pl + pr) - 0.5 * (ur - ul) * ra * ca;
  //   second, compue the q coefficients in Toro (10.60)
  double coeffl = pm <= pl ? 1 : sqrt(1 + 0.5 * (1 + 1/g) * (pm/pl - 1));
  double coeffr = pm <= pr ? 1 : sqrt(1 + 0.5 * (1 + 1/g) * (pm/pr - 1));
  //   finally, compute Toro (10.59)
  double sl = ul - cl * coeffl;
  double sr = ur + cr * coeffr;

  state[0] = rl;
  state[1] = ul;
  state[2] = pl;
  state[3] = rr;
  state[4] = ur;
  state[5] = pr;
  state[6] = sl;
  state[7] = sr;
}

// HLL
GKYL_CU_DH
static void
states_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  double *speeds, double *qm)
{
  // STEP 1. compute min and max wave speeds
  double state[8];
  states_hll_common(eqn, ql, qr, state);
  double ul = state[1], pl = state[2];
  double ur = state[4], pr = state[5];
  double sl = state[6], sr = state[7];

  // STEP 2. compute left and right fluxes
  double fl[5];
  fl[0] = ql[0] * ul;
  fl[1] = ql[1] * ul + pl;
  fl[2] = ql[2] * ul;
  fl[3] = ql[3] * ul;
  fl[4] = ql[4] * ul + pl * ul;

  double fr[5];
  fr[0] = qr[0] * ur;
  fr[1] = qr[1] * ur + pr;
  fr[2] = qr[2] * ur;
  fr[3] = qr[3] * ur;
  fr[4] = qr[4] * ur + pr * ur;

  // STEP 3. compute standard HLL intermediate states
  qm[0] = (sr*qr[0]-sl*ql[0]+fl[0]-fr[0]) / (sr-sl);
  qm[1] = (sr*qr[1]-sl*ql[1]+fl[1]-fr[1]) / (sr-sl);
  qm[2] = (sr*qr[2]-sl*ql[2]+fl[2]-fr[2]) / (sr-sl);
  qm[3] = (sr*qr[3]-sl*ql[3]+fl[3]-fr[3]) / (sr-sl);
  qm[4] = (sr*qr[4]-sl*ql[4]+fl[4]-fr[4]) / (sr-sl);

  // STEP 4. collect all speeds
  speeds[0] = sl;
  speeds[1] = sr;
}

GKYL_CU_DH  
static double
wave_hll(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  double qm[5];
  states_hll(eqn, ql, qr, speeds, qm);

  double *wv;

  wv = waves;
  for (int i=0; i<5; ++i)  wv[i] = qm[i] - ql[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qr[i] - qm[i];

  return fmax(fabs(speeds[0]), fabs(speeds[1]));
}

GKYL_CU_DH
static void
qfluct_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

GKYL_CU_DH
static double
wave_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double phil,
  double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_hll(eqn, delta, ql, qr, waves, s);
  else // FIXME perhaps not needed
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

GKYL_CU_DH
static void
qfluct_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, double phil,
  double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// HLLC
GKYL_CU_DH
static void
states_hllc(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  double *speeds, double *qml, double *qmr)
{
  // STEP 1. compute min and max wave speeds
  double state[8];
  states_hll_common(eqn, ql, qr, state);
  double rl = state[0], ul = state[1], pl = state[2];
  double rr = state[3], ur = state[4], pr = state[5];
  double sl = state[6], sr = state[7];

  // STEP 2. compute middle wave speed, Toro (10.37)
  double sm = (pr-pl+rl*ul*(sl-ul)-rr*ur*(sr-ur)) / (rl*(sl-ul)-rr*(sr-ur));

  // STEP 3. compute left and right intermediate states, Toro (10.39)
  qml[0] = rl * (sl-ul) / (sl-sm);
  qml[1] = qml[0] * sm;
  qml[2] = qml[0] * ql[2] / ql[0];
  qml[3] = qml[0] * ql[3] / ql[0];
  qml[4] = qml[0] * (ql[4]/rl + (sm-ul) * (sm + pl / rl / (sl-ul)));

  qmr[0] = rr * (sr-ur) / (sr-sm);
  qmr[1] = qmr[0] * sm;
  qmr[2] = qmr[0] * qr[2] / qr[0];
  qmr[3] = qmr[0] * qr[3] / qr[0];
  qmr[4] = qmr[0] * (qr[4]/rr + (sm-ur) * (sm + pr / rr / (sr-ur)));

  // STEP 4. collect all speeds
  speeds[0] = sl;
  speeds[1] = sm;
  speeds[2] = sr;
}

GKYL_CU_DH  
static double
wave_hllc(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  double qml[5], qmr[5];
  states_hllc(eqn, ql, qr, speeds, qml, qmr);

  double *wv;

  wv = waves;
  for (int i=0; i<5; ++i)  wv[i] = qml[i] - ql[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qmr[i] - qml[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qr[i] - qmr[i];

  return fmax(fabs(speeds[0]), fabs(speeds[2]));
}

GKYL_CU_DH
static void
qfluct_hllc(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5], *w2 = &waves[10];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
}

GKYL_CU_DH
static double
wave_hllc_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir, double *waves, double *s)
{
  // clear waves and wave speeds
  for (int i=0; i<3; ++i) { // mwaves
    double *w = &waves[5*i];
    for (int j=0; j<5; ++j) // meqn
      w[j] = 0.0;
    s[i] = 0.0;
  }
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else {
    if (type == GKYL_WV_HIGH_ORDER_FLUX)
      return wave_hllc(eqn, delta, ql, qr, waves, s);
    else
      return wave_lax(eqn, delta, ql, qr, waves, s);
  }
  return 0.0; // can't happen
}

GKYL_CU_DH
static void
qfluct_hllc_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX && (phil > 0.0) && (phir > 0.0))
    return qfluct_hllc(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

GKYL_CU_DH
static void
qfluct_hllc_direct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq)
{
  double s[3], qml[5], qmr[5];
  states_hllc(eqn, ql, qr, s, qml, qmr);

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<5; ++i) {
    double w0 = qml[i] - ql[i];
    double w1 = qmr[i] - qml[i];
    double w2 = qr[i] - qmr[i];
    amdq[i] = s0m*w0 + s1m*w1 + s2m*w2;
    apdq[i] = s0p*w0 + s1p*w1 + s2p*w2;
  }
}

GKYL_CU_DH
static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);

  double fr[5], fl[5];
  gkyl_euler_flux(euler->gas_gamma, ql, fl);
  gkyl_euler_flux(euler->gas_gamma, qr, fr);

  for (int m=0; m<5; ++m) flux_jump[m] = fr[m]-fl[m];

  double amaxl = gkyl_euler_max_abs_speed(euler->gas_gamma, ql);
  double amaxr = gkyl_euler_max_abs_speed(euler->gas_gamma, qr);  

  return fmax(amaxl, amaxr);
}

GKYL_CU_DH
static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  
  if (q[0] < 0.0)
    return false;

  double pr = gkyl_euler_pressure(euler->gas_gamma, q);
  if (pr < 0.0)
    return false;

  return true;
}

GKYL_CU_DH
static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return gkyl_euler_max_abs_speed(euler->gas_gamma, q);
}

GKYL_CU_DH
static inline void
euler_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // density and moment as copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  double ke = 0.5*(qin[1]*qin[1] + qin[2]*qin[2] + qin[3]*qin[3])/qin[0];
  diag[4] = ke; 
  diag[5] = qin[4]-ke;
}

GKYL_CU_DH
static inline void
euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 5; i++) {
    sout[i] = 0.0;
  }
}
