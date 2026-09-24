#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_coldfluid.h>

#define RHOU 1
#define RHOV 2
#define RHOW 3

#define SQ(x) ((x)*(x))

struct wv_coldfluid {
  struct gkyl_wv_eqn eqn; // base object
};

static void
coldfluid_flux(const double q[4], double flux[4])
{
  double u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u; // rho*u*u
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
}

static inline void
coldfluid_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // density and moment as copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  double ke = 0.5*(qin[1]*qin[1] + qin[2]*qin[2] + qin[3]*qin[3])/qin[0];
  diag[4] = ke;
}

static void
coldfluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_coldfluid *coldfluid = container_of(base, struct wv_coldfluid, eqn);
  gkyl_free(coldfluid);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    qout[i] = win[i];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  const double phil, const double phir, double *waves, double *s)
{
  double f[4];
  double ur = qr[RHOU]/qr[0], ul = ql[RHOU]/ql[0];

  double *wv = 0;

  if ((ul < 0) && (0 < ur)) { // vacuum intermediate state will be formed
    coldfluid_flux(ql, f);
    wv = &waves[0];
    for(int m=0; m<4; ++m) wv[m] = -f[m];
    s[0] = ul;

    coldfluid_flux(qr, f);
    wv = &waves[4];
    for(int m=0; m<4; ++m) wv[m] = f[m];
    s[1] = ur;
  }
  else {
    // no vacuum state
    double rl = ql[0];
    double rr = qr[0];
    // compute Roe averaged speed
    double uav = (sqrt(rl)*ul + sqrt(rr)*ur)/(sqrt(rl)+sqrt(rr));
            
    if(uav<0) {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = 0.0;
    }
    else {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = 0;

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];
    }
    s[0] = uav;
    s[1] = uav;
  }

  return fmax(fabs(s[0]), fabs(s[1]));
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0)
        amdq[m] += s[mw]*wv[m];
      else
        apdq[m] += s[mw]*wv[m];
    }
  }
}

static void
ffluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0) {
        amdq[m] += wv[m];
      }
      else if (s[mw] > 0.0) {
        apdq[m] += wv[m];
      }
      else {
        amdq[m] += 0.5*wv[m];
        apdq[m] += 0.5*wv[m];
      }
    }
  }
}

static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  double fr[4], fl[4];
  coldfluid_flux(ql, fl);
  coldfluid_flux(qr, fr);

  for (int m=0; m<4; ++m) flux_jump[m] = fr[m]-fl[m];

  double amaxl = ql[RHOU]/ql[0];
  double amaxr = qr[RHOU]/qr[0];

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return q[0] > 0.0;
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_coldfluid *coldfluid = container_of(eqn, struct wv_coldfluid, eqn);
  return fabs(q[RHOU]/q[0]);
}

static inline void
coldfluid_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 4; i++) {
    sout[i] = 0.0;
  }
}

struct gkyl_wv_eqn*
gkyl_wv_coldfluid_new(void)
{
  struct wv_coldfluid *coldfluid = gkyl_malloc(sizeof(struct wv_coldfluid));

  coldfluid->eqn.type = GKYL_EQN_COLDFLUID;
  coldfluid->eqn.num_equations = 4;
  coldfluid->eqn.num_waves = 2;
  coldfluid->eqn.num_diag = 5; // KE is final component
  
  coldfluid->eqn.waves_func = wave_roe;
  coldfluid->eqn.qfluct_func = qfluct_roe;
  coldfluid->eqn.ffluct_func = ffluct_roe;

  coldfluid->eqn.flux_jump = flux_jump;
  coldfluid->eqn.check_inv_func = check_inv;
  coldfluid->eqn.max_speed_func = max_speed;
  coldfluid->eqn.rotate_to_local_func = rot_to_local;
  coldfluid->eqn.rotate_to_global_func = rot_to_global;

  coldfluid->eqn.cons_to_riem = cons_to_riem;
  coldfluid->eqn.riem_to_cons = riem_to_cons;

  coldfluid->eqn.cons_to_diag = coldfluid_cons_to_diag;

  coldfluid->eqn.source_func = coldfluid_source;

  coldfluid->eqn.ref_count = gkyl_ref_count_init(coldfluid_free);

  coldfluid->eqn.embed_geo = NULL;

  return &coldfluid->eqn;
}
