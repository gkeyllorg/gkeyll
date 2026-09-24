#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler_mixture.h>
#include <gkyl_wv_euler_mixture_priv.h>

void
gkyl_euler_mixture_prim_vars(int num_species, double* gas_gamma_s, const double* q, double* v)
{
  double rho_total = q[0];
  double momx_total = q[1];
  double momy_total = q[2];
  double momz_total = q[3];
  double E_total = q[4];

  double *vol_frac_cons_s = gkyl_malloc(sizeof(double[num_species - 1]));
  for (int i = 0; i < num_species - 1; i++) {
    vol_frac_cons_s[i] = q[5 + i];
  }

  double *rho_cons_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_cons_s[i] = q[4 + num_species + i];
  }

  double vx_total = momx_total / rho_total;
  double vy_total = momy_total / rho_total;
  double vz_total = momz_total / rho_total;

  double *vol_frac_s = gkyl_malloc(sizeof(double[num_species]));
  double vol_frac_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    vol_frac_s[i] = vol_frac_cons_s[i] / rho_total;
    vol_frac_total += vol_frac_s[i];
  }
  vol_frac_s[num_species - 1] = 1.0 - vol_frac_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = rho_cons_s[i] / vol_frac_s[i];
  }

  double *p_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    p_s[i] = (gas_gamma_s[i] - 1.0) * (E_total - (0.5 * rho_total * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))));
  }

  double p_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    p_total += vol_frac_s[i] * p_s[i];
  }

  v[0] = rho_total;
  v[1] = vx_total;
  v[2] = vy_total;
  v[3] = vz_total;
  v[4] = p_total;
  for (int i = 0; i < num_species - 1; i++) {
    v[5 + i] = vol_frac_s[i];
  }
  for (int i = 0; i < num_species; i++) {
    v[4 + num_species + i] = rho_s[i];
  }

  gkyl_free(vol_frac_cons_s);
  gkyl_free(rho_cons_s);
  gkyl_free(vol_frac_s);
  gkyl_free(rho_s);
  gkyl_free(p_s);
}

static inline double
gkyl_euler_mixture_max_abs_speed(int num_species, double* gas_gamma_s, const double* q)
{
  double *v = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, q, v);

  double vx_total = v[1];
  double vy_total = v[2];
  double vz_total = v[3];
  double p_total = v[4];

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = v[4 + num_species + i];
  }

  double v_mag = sqrt((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total));

  double max_abs_speed = 0.0;
  for (int i = 0; i < num_species; i++) {
    if (fabs(v_mag) + sqrt(gas_gamma_s[i] * (p_total / rho_s[i])) > max_abs_speed) {
      max_abs_speed = fabs(v_mag) + sqrt(gas_gamma_s[i] * (p_total / rho_s[i]));
    }
  }

  gkyl_free(v);
  gkyl_free(rho_s);

  return max_abs_speed;
}

void
gkyl_euler_mixture_flux(int num_species, double* gas_gamma_s, const double* q, double* flux)
{
  double *v = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, q, v);

  double rho_total = v[0];
  double vx_total = v[1];
  double vy_total = v[2];
  double vz_total = v[3];
  double p_total = v[4];
  double E_total = q[4];

  double *vol_frac_s = gkyl_malloc(sizeof(double[num_species]));
  double vol_frac_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    vol_frac_s[i] = v[5 + i];
    vol_frac_total += vol_frac_s[i];
  }
  vol_frac_s[num_species - 1] = 1.0 - vol_frac_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = v[4 + num_species + i];
  }

  flux[0] = rho_total * vx_total;
  flux[1] = (rho_total * (vx_total * vx_total)) + p_total;
  flux[2] = rho_total * (vx_total * vy_total);
  flux[3] = rho_total * (vx_total * vz_total);
  flux[4] = (E_total * vx_total) + (vx_total * p_total);
  
  for (int i = 0; i < num_species - 1; i++) {
    flux[5 + i] = rho_total * (vx_total * vol_frac_s[i]);
  }
  for (int i = 0; i < num_species; i++) {
    flux[4 + num_species + i] = vol_frac_s[i] * (vx_total * rho_s[i]);
  }

  gkyl_free(v);
  gkyl_free(vol_frac_s);
  gkyl_free(rho_s);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 4 + (2 * num_species); i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 4 + (2 * num_species); i++) {
    qout[i] = win[i];
  }
}

static void
euler_mixture_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
euler_mixture_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    if (i > 0 && i < 4) {
      ghost[i] = -skin[i];
    }
    else {
      ghost[i] = skin[i];
    }
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    qlocal[i] = qglobal[i];
  }

  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    qglobal[i] = qlocal[i];
  }

  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double* gas_gamma_s = euler_mixture->gas_gamma_s;

  double sl = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, ql);
  double sr = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  double *fr = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[4 + (2 * num_species)];
  for (int i = 0; i < 4 + (2 * num_species); i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  gkyl_free(fl);
  gkyl_free(fr);

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  const double *w0 = &waves[0], *w1 = &waves[4 + (2 * num_species)];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  double rho_l = ql[0];
  double rho_r = qr[0];

  double *vl = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  double *vr = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, ql, vl);
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, qr, vr);
  double p_l = vl[4];
  double p_r = vr[4];

  double sqrt_rho_l = sqrt(rho_l);
  double sqrt_rho_r = sqrt(rho_r);
  double roe_avg_l = 1.0 / sqrt_rho_l;
  double roe_avg_r = 1.0 / sqrt_rho_r;
  double roe_avg_sq = 1.0 / (sqrt_rho_l + sqrt_rho_r);

  double vx = ((ql[1] * roe_avg_l) + (qr[1] * roe_avg_r)) * roe_avg_sq;
  double vy = ((ql[2] * roe_avg_l) + (qr[2] * roe_avg_r)) * roe_avg_sq;
  double vz = ((ql[3] * roe_avg_l) + (qr[3] * roe_avg_r)) * roe_avg_sq;
  double enth = (((ql[4] + p_l) * roe_avg_l) + ((qr[4] + p_r) * roe_avg_r)) * roe_avg_sq;

  double *vol_frac_s = gkyl_malloc(sizeof(double[num_species]));
  double vol_frac_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    vol_frac_s[i] = ((ql[5 + i] * roe_avg_l) + (qr[5 + i] * roe_avg_r)) * roe_avg_sq;
    vol_frac_total += vol_frac_s[i];
  }
  vol_frac_s[num_species - 1] = 1.0 - vol_frac_total;

  double gas_gamma_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    gas_gamma_total += vol_frac_s[i] * gas_gamma_s[i];
  }

  double vel_sq = ((vx * vx) + (vy * vy) + (vz * vz));
  double a_sq = (gas_gamma_total - 1.0) * (enth - (0.5 * vel_sq));
  double a = sqrt(a_sq);
  double gamma1_over_a_sq = (gas_gamma_total - 1.0) / a_sq;
  double internal_enth = enth - vel_sq;

  double a4 = gamma1_over_a_sq * ((internal_enth * delta[0]) + (vx * delta[1]) + (vy * delta[2]) + (vz * delta[3]) - delta[4]);
  double a2 = delta[2] - (vy * delta[0]);
  double a3 = delta[3] - (vz * delta[0]);
  double a5 = 0.5 * (delta[1] + (((a - vx) * delta[0]) - (a * a4))) / a;
  double a1 = delta[0] - a4 - a5;

  double *a_vol = gkyl_malloc(sizeof(double[num_species - 1]));
  for (int i = 0; i < num_species; i++) {
    a_vol[i] = delta[5 + i] - (vol_frac_s[i] * delta[0]);
  }

  double *a_rho1 = gkyl_malloc(sizeof(double[num_species]));
  double *a_rho2 = gkyl_malloc(sizeof(double[num_species]));
  double *a_rho3 = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    a_rho2[i] = gamma1_over_a_sq *  -((internal_enth * delta[4 + num_species + i]) + (vx * delta[1]) + (vy * delta[2]) + (vz * delta[3]) - delta[4]);
    a_rho3[i] = 0.5 * (delta[1] + (((a - vx) * delta[4 + num_species + i]) - (a * a_rho2[i]))) / a;
    a_rho1[i] = delta[4 + num_species + i] - a_rho2[i] - a_rho3[i];
  }

  double *wv;
  wv = &waves[0];
  wv[0] = a1;
  wv[1] = a1 * (vx - a);
  wv[2] = a1 * vy;
  wv[3] = a1 * vz;
  wv[4] = a1 * (enth - (vx * a));
  for (int i = 0; i < num_species - 1; i++) {
    wv[5 + i] = a1 * vol_frac_s[i];
  }
  for (int i = 0; i < num_species; i++) {
    wv[4 + num_species + i] = a_rho1[i];
  }
  s[0] = vx - a;

  wv = &waves[4 + (2 * num_species)];
  wv[0] = a4;
  wv[1] = a4 * vx;
  wv[2] = (a4 * vy) + a2;
  wv[3] = (a4 * vz) + a3;
  wv[4] = (0.5 * a4 * vel_sq) + (a2 * vy) + (a3 * vz);
  for (int i = 0; i < num_species - 1; i++) {
    wv[5 + i] = (a4 * vol_frac_s[i]) + a_vol[i];
  }
  for (int i = 0; i < num_species; i++) {
    wv[4 + num_species + i] = a_rho2[i];
  }
  s[1] = vx;

  wv = &waves[8 + (4 * num_species)];
  wv[0] = a5;
  wv[1] = a5 * (vx + a);
  wv[2] = a5 * vy;
  wv[3] = a5 * vz;
  wv[4] = a5 * (enth + (vx * a));
  for (int i = 0; i < num_species - 1; i++) {
    wv[5 + i] = a5 * vol_frac_s[i];
  }
  for (int i = 0; i < num_species; i++) {
    wv[4 + num_species + i] = a_rho3[i];
  }
  s[2] = vx + a;

  gkyl_free(vl);
  gkyl_free(vr);
  gkyl_free(vol_frac_s);
  gkyl_free(a_vol);
  gkyl_free(a_rho1);
  gkyl_free(a_rho2);
  gkyl_free(a_rho3);

  return fabs(vx) + a;
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  
  const double *w0 = &waves[0], *w1 = &waves[4 + (2 * num_species)], *w2 = &waves[8 + (4 * num_species)];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]);
  }
}

static double
wave_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_roe(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  double *fr = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  double *fl = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, qr, fr);

  for (int m = 0; m < 4 + (2 * num_species); m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, ql);
  double amaxr = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, qr);
  
  gkyl_free(fr);
  gkyl_free(fl);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  double *v = gkyl_malloc(sizeof(double[4 + (2 * num_species)]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, q, v);

  for (int i = 0; i < num_species; i++) {
    if (v[4 + num_species + i] < 0.0) {
      gkyl_free(v);
      return false;
    }
  }
  
  double vol_frac_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    vol_frac_total += v[5 + i];
  }
  if (vol_frac_total > 1.0) {
    gkyl_free(v);
    return false;
  }

  if (v[0] < 0.0 || v[4] < 0.0) {
    gkyl_free(v);
    return false;
  }
  else {
    gkyl_free(v);
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double* gas_gamma_s = euler_mixture->gas_gamma_s;

  return gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, q);
}

static inline void
euler_mixture_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    diag[i] = qin[i];
  }
}

static inline void
euler_mixture_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + (2 * num_species); i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_euler_mixture_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_euler_mixture *euler_mixture = container_of(base->on_dev, struct wv_euler_mixture, eqn);
    gkyl_cu_free(euler_mixture);
  }

  struct wv_euler_mixture *euler_mixture = container_of(base, struct wv_euler_mixture, eqn);
  gkyl_free(euler_mixture);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_mixture_new(int num_species, double* gas_gamma_s, bool use_gpu)
{
  return gkyl_wv_euler_mixture_inew(&(struct gkyl_wv_euler_mixture_inp) {
      .num_species = num_species,
      .gas_gamma_s = gas_gamma_s,
      .rp_type = WV_EULER_MIXTURE_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_euler_mixture_inew(const struct gkyl_wv_euler_mixture_inp* inp)
{
  struct wv_euler_mixture *euler_mixture = gkyl_malloc(sizeof(struct wv_euler_mixture));

  euler_mixture->eqn.type = GKYL_EQN_EULER_MIXTURE;
  euler_mixture->eqn.num_equations = 4 + (2 * inp->num_species);
  euler_mixture->eqn.num_diag = 4 + (2 * inp->num_species);

  euler_mixture->num_species = inp->num_species;
  euler_mixture->gas_gamma_s = inp->gas_gamma_s;

  if (inp->rp_type == WV_EULER_MIXTURE_RP_LAX) {
    euler_mixture->eqn.num_waves = 2;
    euler_mixture->eqn.waves_func = wave_lax_l;
    euler_mixture->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_EULER_MIXTURE_RP_ROE) {
    euler_mixture->eqn.num_waves = 3;
    euler_mixture->eqn.waves_func = wave_roe_l;
    euler_mixture->eqn.qfluct_func = qfluct_roe_l;
  }

  euler_mixture->eqn.flux_jump = flux_jump;
  euler_mixture->eqn.check_inv_func = check_inv;
  euler_mixture->eqn.max_speed_func = max_speed;
  euler_mixture->eqn.rotate_to_local_func = rot_to_local;
  euler_mixture->eqn.rotate_to_global_func = rot_to_global;
  
  euler_mixture->eqn.wall_bc_func = euler_mixture_wall;
  euler_mixture->eqn.no_slip_bc_func = euler_mixture_no_slip;

  euler_mixture->eqn.cons_to_riem = cons_to_riem;
  euler_mixture->eqn.riem_to_cons = riem_to_cons;

  euler_mixture->eqn.cons_to_diag = euler_mixture_cons_to_diag;

  euler_mixture->eqn.source_func = euler_mixture_source;

  euler_mixture->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler_mixture->eqn.flags);
  euler_mixture->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_mixture_free);
  euler_mixture->eqn.on_dev = &euler_mixture->eqn; // On the CPU, the equation object points ot itself.

  euler_mixture->eqn.embed_geo = NULL;

  return &euler_mixture->eqn;
}

int
gkyl_wv_euler_mixture_num_species(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  return num_species;
}

double*
gkyl_wv_euler_mixture_gas_gamma_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  return gas_gamma_s;
}
