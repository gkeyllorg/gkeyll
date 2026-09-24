#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_sr_euler.h>
#include <gkyl_wv_sr_euler.h>

struct wv_sr_euler {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
};

static void
sr_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_sr_euler *sr_euler = container_of(base, struct wv_sr_euler, eqn);
  gkyl_free(sr_euler);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    qout[i] = win[i];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  // Mass density and energy are scalars
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1];
  // Rotate momentum to local coordinates
  qlocal[2] = qglobal[2]*norm[0] + qglobal[3]*norm[1] + qglobal[4]*norm[2];
  qlocal[3] = qglobal[2]*tau1[0] + qglobal[3]*tau1[1] + qglobal[4]*tau1[2];
  qlocal[4] = qglobal[2]*tau2[0] + qglobal[3]*tau2[1] + qglobal[4]*tau2[2];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  // Mass density and energy are scalars
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1];
  // Rotate momentum back to global coordinates
  qglobal[2] = qlocal[2]*norm[0] + qlocal[3]*tau1[0] + qlocal[4]*tau2[0];
  qglobal[3] = qlocal[2]*norm[1] + qlocal[3]*tau1[1] + qlocal[4]*tau2[1];
  qglobal[4] = qlocal[2]*norm[2] + qlocal[3]*tau1[2] + qlocal[4]*tau2[2];
  
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  const double phil, const double phir, double *waves, double *s)
{
  const struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);
  double vl[5], vr[5];
  double gas_gamma = sr_euler->gas_gamma;
  double g1 = gas_gamma - 1.;
  double gFrac = gas_gamma/g1;

  // Get prim variables rho, p, u, v, w. 
  gkyl_sr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_sr_euler_prim_vars(gas_gamma, qr, vr);
  double gammal = 1. / sqrt(1. - (vl[2]*vl[2] + vl[3]*vl[3] + vl[4]*vl[4]));
  double pl = vl[1];
  double gammar = 1. / sqrt(1. - (vr[2]*vr[2] + vr[3]*vr[3] + vr[4]*vr[4]));
  double pr = vr[1];
		    
  //Equation numbers in all of the following follows Eulderink AASS 110, 587 (1995)
  //Roe Averages
  //double Kl = sqrt(rhol + gFrac*pl),  Kr = sqrt(rhor + gFrac*pr); //sqrt(10.3)
   
  double Kl = sqrt(ql[1] + pl) / gammal,  Kr = sqrt(qr[1] + pr) / gammar;  //sqrt(10.3)
  double ravgK = 1./(Kl + Kr);
  double v0 = (Kl*gammal + Kr*gammar)*ravgK; //10.7
  double v1 = (Kl*gammal*vl[2] + Kr*gammar*vr[2])*ravgK;
  double v2 = (Kl*gammal*vl[3] + Kr*gammar*vr[3])*ravgK;  
  double v3 = (Kl*gammal*vl[4] + Kr*gammar*vr[4])*ravgK;
  double v4 = (pl/Kl + pr/Kr)*ravgK;
  double cm = 1. - gFrac*v4, cp = 1. + gFrac*v4;
  
  double vava = -v0*v0 + v1*v1 + v2*v2 + v3*v3; //v_alpha v^alpha

  double s2 = 0.5*gas_gamma*v4*(1-vava) - 0.5*g1*(1+vava); //10.13
  double e = v0*v0 - v1*v1; //10.14
  double y = sqrt((1-gas_gamma*v4)*e + s2); //10.14

  
  // Compute projections of jump, Eq 10.16
  double k = v0*delta[1] - v1*delta[2];
  double vada = -v0*delta[1] + v1*delta[2] + v2*delta[3] + v3*delta[4]; // v_alpha Delta^alpha
  double a1 = -(s2*k + sqrt(s2)*y*(v0*delta[2] - v1*delta[1]) + g1*e*(delta[0]+cp*vada))/(2.*e*s2);
  double a2 = -(s2*k - sqrt(s2)*y*(v0*delta[2] - v1*delta[1]) + g1*e*(delta[0]+cp*vada))/(2.*e*s2);
  double a3 = (2.*s2*k + g1*e*(delta[0]+cp*vada))/(e*s2);
  double a4 = delta[3] - k*v2 / e;
  double a5 = delta[4] - k*v3 / e;

  double *wv;
  // Wave 1: eigenvalue is lambda- 
  wv = &waves[0];
  wv[0] = a1*cm;
  wv[1] = a1*(v0-sqrt(s2)*v1/y);
  wv[2] = a1*(v1-sqrt(s2)*v0/y);
  wv[3] = a1*v2;
  wv[4] = a1*v3;
  s[0] = ((1.-gas_gamma*v4)*v0*v1 - sqrt(s2)*y)/((1.-gas_gamma*v4)*v0*v0+s2);//10.12

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a3*(cm+s2/g1) - a4*cp*v2 - a5*cp*v3;
  wv[1] = a3*v0;
  wv[2] = a3*v1;
  wv[3] = a3*v2 + a4;
  wv[4] = a3*v3 + a5;
  s[1] = v1 / v0;

  // Wave 3: eigenvalue is lambda+
  wv = &waves[10];
  wv[0] = a2*cm;
  wv[1] = a2*(v0+sqrt(s2)*v1/y);
  wv[2] = a2*(v1+sqrt(s2)*v0/y);
  wv[3] = a2*v2;
  wv[4] = a2*v3;;
  s[2] = ((1.-gas_gamma*v4)*v0*v1 + sqrt(s2)*y)/((1.-gas_gamma*v4)*v0*v0+s2);//10.12
  
  return ((1.-gas_gamma*v4)*v0*fabs(v1) + sqrt(s2)*y)/((1.-gas_gamma*v4)*v0*v0+s2);
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
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

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return true; // TODO
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);
  return gkyl_sr_euler_max_abs_speed(sr_euler->gas_gamma, q);
}

static inline void
sr_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 5; i++) {
    sout[i] = 0.0;
  }
}

struct gkyl_wv_eqn*
gkyl_wv_sr_euler_new(double gas_gamma)
{
  struct wv_sr_euler *sr_euler = gkyl_malloc(sizeof(struct wv_sr_euler));

  sr_euler->eqn.type = GKYL_EQN_SR_EULER;
  sr_euler->eqn.num_equations = 5;
  sr_euler->eqn.num_waves = 3;
  sr_euler->eqn.num_diag = 5;
  
  sr_euler->gas_gamma = gas_gamma;
  sr_euler->eqn.waves_func = wave_roe;
  sr_euler->eqn.qfluct_func = qfluct_roe;

  sr_euler->eqn.check_inv_func = check_inv;
  sr_euler->eqn.max_speed_func = max_speed;
  sr_euler->eqn.rotate_to_local_func = rot_to_local;
  sr_euler->eqn.rotate_to_global_func = rot_to_global;

  sr_euler->eqn.cons_to_riem = cons_to_riem;
  sr_euler->eqn.riem_to_cons = riem_to_cons;

  sr_euler->eqn.cons_to_diag = gkyl_default_cons_to_diag;

  sr_euler->eqn.source_func = sr_euler_source;

  sr_euler->eqn.ref_count = gkyl_ref_count_init(sr_euler_free);

  return &sr_euler->eqn;
}

double
gkyl_wv_sr_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);
  return sr_euler->gas_gamma;
}
