#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_reactive_euler.h>
#include <gkyl_wv_reactive_euler_priv.h>

void
gkyl_reactive_euler_prim_vars(double gas_gamma, double energy_of_formation, const double q[6], double v[6])
{
  double rho = q[0];
  double momx = q[1];
  double momy = q[2];
  double momz = q[3];
  double Etot = q[4];
  double reaction_density = q[5];

  double specific_internal_energy = (Etot / rho) - (0.5 * ((momx * momx) + (momy * momy) + (momz * momz)) / (rho * rho)) -
    (energy_of_formation * ((reaction_density / rho) - 1.0));

  v[0] = rho;
  v[1] = momx / rho;
  v[2] = momy / rho;
  v[3] = momz / rho;
  v[4] = specific_internal_energy * (gas_gamma - 1.0) * rho;
  v[5] = reaction_density / rho;
}

static inline double
gkyl_reactive_euler_max_abs_speed(double gas_gamma, double energy_of_formation, const double q[6])
{
  double v[6] = { 0.0 };
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, q, v);
  
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double v_mag = sqrt((vx * vx) + (vy * vy) + (vz * vz));

  return fabs(v_mag) + sqrt(gas_gamma * (p / rho));
}

void
gkyl_reactive_euler_flux(double gas_gamma, double energy_of_formation, const double q[6], double flux[6])
{
  double v[6] = { 0.0 };
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, q, v);

  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];
  double reaction_progress = v[5];

  double Etot = q[4];

  flux[0] = rho * vx;
  flux[1] = (rho * (vx * vx)) + p;
  flux[2] = rho * (vx * vy);
  flux[3] = rho * (vx * vz);
  flux[4] = (Etot + p) * vx;
  flux[5] = rho * (vx * reaction_progress);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 6; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 6; i++) {
    qout[i] = win[i];
  }
}

static void
reactive_euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 6; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
reactive_euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];
  ghost[5] = skin[5];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];
  qlocal[5] = qglobal[5];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];
  qglobal[5] = qlocal[5];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;
  double energy_of_formation = reactive_euler->energy_of_formation;

  double sl = gkyl_reactive_euler_max_abs_speed(gas_gamma, energy_of_formation, ql);
  double sr = gkyl_reactive_euler_max_abs_speed(gas_gamma, energy_of_formation, qr);
  double amax = fmax(sl, sr);

  double fl[6], fr[6];
  gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, ql, fl);
  gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[6];
  for (int i = 0; i < 6; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[6];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 6; i++) {
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
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;
  double energy_of_formation = reactive_euler->energy_of_formation;

  double rho_l = ql[0];
  double rho_r = qr[0];

  double vl[6] = { 0.0 };
  double vr[6] = { 0.0 };
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, ql, vl);
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, qr, vr);
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
  double reaction_progress = ((ql[5] * roe_avg_l) + (qr[5] * roe_avg_r)) * roe_avg_sq;

  double vel_sq = ((vx * vx) + (vy * vy) + (vz * vz));
  double a_sq = (gas_gamma - 1.0) * (enth - (0.5 * vel_sq) - (energy_of_formation * (reaction_progress - 1.0)));
  double a = sqrt(a_sq);
  double gamma1_over_a_sq = (gas_gamma - 1.0) / a_sq;
  double internal_enth = enth - vel_sq;

  double a4 = gamma1_over_a_sq * ((internal_enth * delta[0]) + (vx * delta[1]) + (vy * delta[2]) + (vz * delta[3]) - delta[4]);
  double a2 = delta[2] - (vy * delta[0]);
  double a3 = delta[3] - (vz * delta[0]);
  double a5 = 0.5 * (delta[1] + (((a - vx) * delta[0]) - (a * a4))) / a;
  double a1 = delta[0] - a4 - a5;
  double a6 = delta[5] - (reaction_progress * delta[0]);

  double *wv;
  wv = &waves[0];
  wv[0] = a1;
  wv[1] = a1 * (vx - a);
  wv[2] = a1 * vy;
  wv[3] = a1 * vz;
  wv[4] = a1 * (enth - (vx * a));
  wv[5] = a1 * reaction_progress;
  s[0] = vx - a;

  wv = &waves[6];
  wv[0] = a4;
  wv[1] = a4 * vx;
  wv[2] = (a4 * vy) + a2;
  wv[3] = (a4 * vz) + a3;
  wv[4] = (0.5 * a4 * vel_sq) + (a2 * vy) + (a3 * vz);
  wv[5] = (a4 * reaction_progress) + a6;
  s[1] = vx;

  wv = &waves[12];
  wv[0] = a5;
  wv[1] = a5 * (vx + a);
  wv[2] = a5 * vy;
  wv[3] = a5 * vz;
  wv[4] = a5 * (enth + (vx * a));
  wv[5] = a5 * reaction_progress;
  s[2] = vx + a;

  return fabs(vx) + a;
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[6], *w2 = &waves[12];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 6; i++) {
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
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;
  double energy_of_formation = reactive_euler->energy_of_formation;

  double fr[6], fl[6];
  gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, ql, fl);
  gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, qr, fr);

  for (int m = 0; m < 6; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_reactive_euler_max_abs_speed(gas_gamma, energy_of_formation, ql);
  double amaxr = gkyl_reactive_euler_max_abs_speed(gas_gamma, energy_of_formation, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;
  double energy_of_formation = reactive_euler->energy_of_formation;

  double v[6] = { 0.0 };
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, q, v);

  if (v[0] < 0.0 || v[4] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;
  double energy_of_formation = reactive_euler->energy_of_formation;

  return gkyl_reactive_euler_max_abs_speed(gas_gamma, energy_of_formation, q);
}

static inline void
reactive_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 5; i++) {
    diag[i] = qin[i];
  }
}

static inline void
reactive_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;
  double specific_heat_capacity = reactive_euler->specific_heat_capacity;
  double energy_of_formation = reactive_euler->energy_of_formation;
  double ignition_temperature = reactive_euler->ignition_temperature;
  double reaction_rate = reactive_euler->reaction_rate;

  double v[6] = { 0.0 };
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, qin, v);

  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double reaction_progress = v[5];

  double specific_internal_energy = (qin[4] / rho) - (0.5 * ((vx * vx) + (vy * vy) + (vz * vz))) -
    (energy_of_formation * (reaction_progress - 1.0));
  double temperature = specific_internal_energy / specific_heat_capacity;

  for (int i = 0; i < 5; i++) {
    sout[i] = 0.0;
  }

  if (temperature > ignition_temperature) {
    sout[5] = -(rho * reaction_progress * reaction_rate);
  }
  else {
    sout[5] = 0.0;
  }
}

void
gkyl_reactive_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_reactive_euler *reactive_euler = container_of(base->on_dev, struct wv_reactive_euler, eqn);
    gkyl_cu_free(reactive_euler);
  }

  struct wv_reactive_euler *reactive_euler = container_of(base, struct wv_reactive_euler, eqn);
  gkyl_free(reactive_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_reactive_euler_new(double gas_gamma, double specific_heat_capacity, double energy_of_formation, double ignition_temperature,
  double reaction_rate, bool use_gpu)
{
  return gkyl_wv_reactive_euler_inew(&(struct gkyl_wv_reactive_euler_inp) {
      .gas_gamma = gas_gamma,
      .specific_heat_capacity = specific_heat_capacity,
      .energy_of_formation = energy_of_formation,
      .ignition_temperature = ignition_temperature,
      .reaction_rate = reaction_rate,
      .rp_type = WV_REACTIVE_EULER_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_reactive_euler_inew(const struct gkyl_wv_reactive_euler_inp* inp)
{
  struct wv_reactive_euler *reactive_euler = gkyl_malloc(sizeof(struct wv_reactive_euler));

  reactive_euler->eqn.type = GKYL_EQN_REACTIVE_EULER;
  reactive_euler->eqn.num_equations = 6;
  reactive_euler->eqn.num_diag = 5;

  reactive_euler->gas_gamma = inp->gas_gamma;
  reactive_euler->specific_heat_capacity = inp->specific_heat_capacity;
  reactive_euler->energy_of_formation = inp->energy_of_formation;
  reactive_euler->ignition_temperature = inp->ignition_temperature;
  reactive_euler->reaction_rate = inp->reaction_rate;

  if (inp->rp_type == WV_REACTIVE_EULER_RP_LAX) {
    reactive_euler->eqn.num_waves = 2;
    reactive_euler->eqn.waves_func = wave_lax_l;
    reactive_euler->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_REACTIVE_EULER_RP_ROE) {
    reactive_euler->eqn.num_waves = 3;
    reactive_euler->eqn.waves_func = wave_roe_l;
    reactive_euler->eqn.qfluct_func = qfluct_roe_l;
  }

  reactive_euler->eqn.flux_jump = flux_jump;
  reactive_euler->eqn.check_inv_func = check_inv;
  reactive_euler->eqn.max_speed_func = max_speed;
  reactive_euler->eqn.rotate_to_local_func = rot_to_local;
  reactive_euler->eqn.rotate_to_global_func = rot_to_global;

  reactive_euler->eqn.wall_bc_func = reactive_euler_wall;
  reactive_euler->eqn.no_slip_bc_func = reactive_euler_no_slip;

  reactive_euler->eqn.cons_to_riem = cons_to_riem;
  reactive_euler->eqn.riem_to_cons = riem_to_cons;

  reactive_euler->eqn.cons_to_diag = reactive_euler_cons_to_diag;

  reactive_euler->eqn.source_func = reactive_euler_source;

  reactive_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(reactive_euler->eqn.flags);
  reactive_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_reactive_euler_free);
  reactive_euler->eqn.on_dev = &reactive_euler->eqn; // On the CPU, the equation object points to itself.

  reactive_euler->eqn.embed_geo = NULL;

  return &reactive_euler->eqn;
}

double
gkyl_wv_reactive_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double gas_gamma = reactive_euler->gas_gamma;

  return gas_gamma;
}

double
gkyl_wv_reactive_euler_specific_heat_capacity(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double specific_heat_capacity = reactive_euler->specific_heat_capacity;

  return specific_heat_capacity;
}

double
gkyl_wv_reactive_euler_energy_of_formation(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double energy_of_formation = reactive_euler->energy_of_formation;

  return energy_of_formation;
}

double
gkyl_wv_reactive_euler_ignition_temperature(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double ignition_temperature = reactive_euler->ignition_temperature;

  return ignition_temperature;
}

double
gkyl_wv_reactive_euler_reaction_rate(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_reactive_euler *reactive_euler = container_of(eqn, struct wv_reactive_euler, eqn);
  double reaction_rate = reactive_euler->reaction_rate;

  return reaction_rate;
}
