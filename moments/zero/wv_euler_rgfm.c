#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler_rgfm.h>
#include <gkyl_wv_euler_rgfm_priv.h>

void
gkyl_euler_rgfm_prim_vars(int num_species, double* gas_gamma_s, const double* q, double* v)
{
  double rho_total = q[0];
  double momx_total = q[1];
  double momy_total = q[2];
  double momz_total = q[3];
  double E_total = q[4];
  double reinit_param = q[5 + (2 * num_species)];

  double *level_set_cons_s = gkyl_malloc(sizeof(double[num_species - 1]));
  for (int i = 0; i < num_species - 1; i++) {
    level_set_cons_s[i] = q[5 + i];
  }

  double *rho_cons_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_cons_s[i] = q[4 + num_species + i];
  }

  double vx_total = momx_total / rho_total;
  double vy_total = momy_total / rho_total;
  double vz_total = momz_total / rho_total;

  double *level_set_s = gkyl_malloc(sizeof(double[num_species]));
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_s[i] = level_set_cons_s[i] / rho_total;
    level_set_total += level_set_s[i];
  }
  level_set_s[num_species - 1] = 1.0 - level_set_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = rho_cons_s[i] / level_set_s[i];
  }

  double *p_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    p_s[i] = (gas_gamma_s[i] - 1.0) * (E_total - (0.5 * rho_total * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))));
  }

  double p_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    p_total += level_set_s[i] * p_s[i];
  }

  v[0] = rho_total;
  v[1] = vx_total;
  v[2] = vy_total;
  v[3] = vz_total;
  v[4] = p_total;
  for (int i = 0; i < num_species - 1; i++) {
    v[5 + i] = level_set_s[i];
  }
  for (int i = 0; i < num_species; i++) {
    v[4 + num_species + i] = rho_s[i];
  }
  v[4 + (2 * num_species)] = reinit_param;

  gkyl_free(level_set_cons_s);
  gkyl_free(rho_cons_s);
  gkyl_free(level_set_s);
  gkyl_free(rho_s);
  gkyl_free(p_s);
}

static inline double
gkyl_euler_rgfm_max_abs_speed(int num_species, double* gas_gamma_s, const double* q)
{
  double *v = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  gkyl_euler_rgfm_prim_vars(num_species, gas_gamma_s, q, v);

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
gkyl_euler_rgfm_flux(int num_species, double* gas_gamma_s, const double* q, double* flux)
{
  double *v = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  gkyl_euler_rgfm_prim_vars(num_species, gas_gamma_s, q, v);

  double rho_total = v[0];
  double vx_total = v[1];
  double vy_total = v[2];
  double vz_total = v[3];
  double p_total = v[4];
  double E_total = q[4];
  double reinit_param = v[4 + (2 * num_species)];

  double *level_set_s = gkyl_malloc(sizeof(double[num_species]));
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_s[i] = v[5 + i];
    level_set_total += level_set_s[i];
  }
  level_set_s[num_species - 1] = 1.0 - level_set_total;

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
    flux[5 + i] = rho_total * (vx_total * level_set_s[i]);
  }
  for (int i = 0; i < num_species; i++) {
    flux[4 + num_species + i] = level_set_s[i] * (vx_total * rho_s[i]);
  }

  flux[4 + (2 * num_species)] = 0.0;

  gkyl_free(v);
  gkyl_free(level_set_s);
  gkyl_free(rho_s);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 5 + (2 * num_species); i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 5 + (2 * num_species); i++) {
    qout[i] = win[i];
  }
}

static void
euler_rgfm_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  for (int i = 0; i < 5 + (2 * num_species); i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
euler_rgfm_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  for (int i = 0; i < 5 + (2 * num_species); i++) {
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
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  for (int i = 0; i < 5 + (2 * num_species); i++) {
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
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  for (int i = 0; i < 5 + (2 * num_species); i++) {
    qglobal[i] = qlocal[i];
  }

  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;
  double* gas_gamma_s = euler_rgfm->gas_gamma_s;
  int reinit_freq = euler_rgfm->reinit_freq;

  double sl = gkyl_euler_rgfm_max_abs_speed(num_species, gas_gamma_s, ql);
  double sr = gkyl_euler_rgfm_max_abs_speed(num_species, gas_gamma_s, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  double *fr = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  gkyl_euler_rgfm_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_euler_rgfm_flux(num_species, gas_gamma_s, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[5 + (2 * num_species)];
  for (int i = 0; i < 5 + (2 * num_species); i++) {
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
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  const double *w0 = &waves[0], *w1 = &waves[5 + (2 * num_species)];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 5 + (2 * num_species); i++) {
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
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;
  double *gas_gamma_s = euler_rgfm->gas_gamma_s;
  int reinit_freq = euler_rgfm->reinit_freq;

  double *fr = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  double *fl = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  gkyl_euler_rgfm_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_euler_rgfm_flux(num_species, gas_gamma_s, qr, fr);

  for (int m = 0; m < 5 + (2 * num_species); m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_euler_rgfm_max_abs_speed(num_species, gas_gamma_s, ql);
  double amaxr = gkyl_euler_rgfm_max_abs_speed(num_species, gas_gamma_s, qr);
  
  gkyl_free(fr);
  gkyl_free(fl);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;
  double *gas_gamma_s = euler_rgfm->gas_gamma_s;

  double *v = gkyl_malloc(sizeof(double[5 + (2 * num_species)]));
  gkyl_euler_rgfm_prim_vars(num_species, gas_gamma_s, q, v);

  for (int i = 0; i < num_species; i++) {
    if (v[5 + num_species + i] < 0.0) {
      gkyl_free(v);
      return false;
    }
  }
  
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_total += v[5 + i];
  }
  if (level_set_total > 1.0) {
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
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;
  double* gas_gamma_s = euler_rgfm->gas_gamma_s;

  return gkyl_euler_rgfm_max_abs_speed(num_species, gas_gamma_s, q);
}

static inline void
euler_rgfm_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  for (int i = 0; i < 5 + (2 * num_species); i++) {
    diag[i] = qin[i];
  }
}

static inline void
euler_rgfm_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  for (int i = 0; i < 5 + (2 * num_species); i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_euler_rgfm_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_euler_rgfm *euler_rgfm = container_of(base->on_dev, struct wv_euler_rgfm, eqn);
    gkyl_cu_free(euler_rgfm);
  }

  struct wv_euler_rgfm *euler_rgfm = container_of(base, struct wv_euler_rgfm, eqn);
  gkyl_free(euler_rgfm);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_rgfm_new(int num_species, double* gas_gamma_s, int reinit_freq, bool use_gpu)
{
  return gkyl_wv_euler_rgfm_inew(&(struct gkyl_wv_euler_rgfm_inp) {
      .num_species = num_species,
      .gas_gamma_s = gas_gamma_s,
      .reinit_freq = reinit_freq,
      .rp_type = WV_EULER_RGFM_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_euler_rgfm_inew(const struct gkyl_wv_euler_rgfm_inp* inp)
{
  struct wv_euler_rgfm *euler_rgfm = gkyl_malloc(sizeof(struct wv_euler_rgfm));

  euler_rgfm->eqn.type = GKYL_EQN_EULER_RGFM;
  euler_rgfm->eqn.num_equations = 5 + (2 * inp->num_species);
  euler_rgfm->eqn.num_diag = 5 + (2 * inp->num_species);

  euler_rgfm->num_species = inp->num_species;
  euler_rgfm->gas_gamma_s = inp->gas_gamma_s;
  euler_rgfm->reinit_freq = inp->reinit_freq;

  if (inp->rp_type == WV_EULER_RGFM_RP_LAX) {
    euler_rgfm->eqn.num_waves = 2;
    euler_rgfm->eqn.waves_func = wave_lax_l;
    euler_rgfm->eqn.qfluct_func = qfluct_lax_l;
  }

  euler_rgfm->eqn.flux_jump = flux_jump;
  euler_rgfm->eqn.check_inv_func = check_inv;
  euler_rgfm->eqn.max_speed_func = max_speed;
  euler_rgfm->eqn.rotate_to_local_func = rot_to_local;
  euler_rgfm->eqn.rotate_to_global_func = rot_to_global;
  
  euler_rgfm->eqn.wall_bc_func = euler_rgfm_wall;
  euler_rgfm->eqn.no_slip_bc_func = euler_rgfm_no_slip;

  euler_rgfm->eqn.cons_to_riem = cons_to_riem;
  euler_rgfm->eqn.riem_to_cons = riem_to_cons;

  euler_rgfm->eqn.cons_to_diag = euler_rgfm_cons_to_diag;

  euler_rgfm->eqn.source_func = euler_rgfm_source;

  euler_rgfm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler_rgfm->eqn.flags);
  euler_rgfm->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_rgfm_free);
  euler_rgfm->eqn.on_dev = &euler_rgfm->eqn; // On the CPU, the equation object points ot itself.

  euler_rgfm->eqn.embed_geo = NULL;

  return &euler_rgfm->eqn;
}

int
gkyl_wv_euler_rgfm_num_species(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;

  return num_species;
}

double*
gkyl_wv_euler_rgfm_gas_gamma_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  double *gas_gamma_s = euler_rgfm->gas_gamma_s;

  return gas_gamma_s;
}

int
gkyl_wv_euler_rgfm_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int reinit_freq = euler_rgfm->reinit_freq;

  return reinit_freq;
}
