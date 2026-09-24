#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_medium.h>
#include <gkyl_wv_gr_medium_priv.h>

void
gkyl_gr_medium_prim_vars(double gas_gamma, const double q[15], double v[15])
{
  double exp_2a = q[0];

  double a_dt = q[1], a_dx = q[2];
  double b_dt = q[3], b_dx = q[4];
  double c_dt = q[5], c_dx = q[6];
  
  double a_dt_dx = q[7], a_dx_dx = q[8];
  double b_dt_dx = q[9], b_dx_dx = q[10];
  double c_dt_dx = q[11], c_dx_dx = q[12];

  double Etot = q[13];
  double mom = q[14];

  double rho = (1.0 / (gas_gamma - 1.0)) * ((-0.5 * (2.0 - gas_gamma) * Etot) + sqrt((0.25 * (2.0 - gas_gamma) * (2.0 - gas_gamma) * Etot * Etot) +
    ((gas_gamma - 1.0) * ((Etot * Etot) - (mom * mom)))));

  double vel = 0.0;
  if (fabs(mom) > pow(10.0, -8.0)) {
    vel = ((gas_gamma * rho) / (2.0 * mom)) * (sqrt(1.0 + ((4 * mom * mom) / ((gas_gamma * gas_gamma) * (rho * rho)))) - 1.0);
  }

  v[0] = exp_2a;

  v[1] = a_dt; v[2] = a_dx;
  v[3] = b_dt; v[4] = b_dx;
  v[5] = c_dt; v[6] = c_dx;

  v[7] = a_dt_dx; v[8] = a_dx_dx;
  v[9] = b_dt_dx; v[10] = b_dx_dx;
  v[11] = c_dt_dx; v[12] = c_dx_dx;

  v[13] = rho;
  v[14] = vel;
}

static inline double
gkyl_gr_medium_max_abs_speed(double gas_gamma, const double q[15])
{
  double v[15] = { 0.0 };
  gkyl_gr_medium_prim_vars(gas_gamma, q, v);

  double vel = v[14];
  double c_s = gas_gamma - 1.0;

  //return fabs(vel) + c_s;
  return 1.0; // Return speed of light.
}

void
gkyl_gr_medium_flux(double gas_gamma, double kappa, const double q[15], double flux[15])
{
  double v[15] = { 0.0 };
  gkyl_gr_medium_prim_vars(gas_gamma, q, v);
  double exp_2a = v[0];
  
  double a_dt = v[1], a_dx = v[2];
  double b_dt = v[3], b_dx = v[4];
  double c_dt = v[5], c_dx = v[6];

  double a_dt_dx = v[7], a_dx_dx = v[8];
  double b_dt_dx = v[9], b_dx_dx = v[10];
  double c_dt_dx = v[11], c_dx_dx = v[12];

  double rho = v[13];
  double vel = v[14];
  double p = (gas_gamma - 1.0) * rho;

  double W = 1.0 / sqrt(1.0 - (vel * vel));
  if (vel * vel > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(pow(10.0, -8.0));
  }

  double Etot = ((rho + p) * (W * W)) - p;
  double mom = (rho + p) * vel * (W * W);

  flux[0] = 0.0;

  flux[1] = 0.0; flux[2] = 0.0;
  flux[3] = 0.0; flux[4] = 0.0;
  flux[5] = 0.0; flux[6] = 0.0;

  flux[7] = -a_dx_dx + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p)));
  flux[8] = -a_dt_dx;
  flux[9] = -b_dx_dx - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p)));
  flux[10] = -b_dt_dx;
  flux[11] = -c_dx_dx;
  flux[12] = -c_dt_dx;

  flux[13] = mom;
  flux[14] = (mom * vel) + p;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 15; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 15; i++) {
    qout[i] = win[i];
  }
}

static void
gr_medium_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 14; i++) {
    ghost[i] = skin[i];
  }
  
  ghost[14] = -skin[14];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  for (int i = 0; i < 15; i++) {
    qlocal[i] = qglobal[i];
  }
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  for (int i = 0; i < 15; i++) {
    qglobal[i] = qlocal[i];
  }
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;
  double kappa = gr_medium->kappa;

  double sl = gkyl_gr_medium_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_medium_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[15], fr[15];
  gkyl_gr_medium_flux(gas_gamma, kappa, ql, fl);
  gkyl_gr_medium_flux(gas_gamma, kappa, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[15];
  for (int i = 0; i < 15; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[15];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 15; i++) {
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
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;
  double kappa = gr_medium->kappa;

  double fr[15], fl[15];
  gkyl_gr_medium_flux(gas_gamma, kappa, ql, fl);
  gkyl_gr_medium_flux(gas_gamma, kappa, qr, fr);

  for (int m = 0; m < 15; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_gr_medium_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_gr_medium_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;

  double v[15] = { 0.0 };
  gkyl_gr_medium_prim_vars(gas_gamma, q, v);

  if (v[13] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;

  return gkyl_gr_medium_max_abs_speed(gas_gamma, q);
}

static inline void
gr_medium_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 15; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_medium_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;
  double kappa = gr_medium->kappa;

  double v[15] = { 0.0 };
  gkyl_gr_medium_prim_vars(gas_gamma, qin, v);
  double exp_2a = v[0];

  double a_dt = v[1], a_dx = v[2];
  double b_dt = v[3], b_dx = v[4];
  double c_dt = v[5], c_dx = v[6];

  double a_dt_dx = v[7], a_dx_dx = v[8];
  double b_dt_dx = v[9], b_dx_dx = v[10];
  double c_dt_dx = v[11], c_dx_dx = v[12];

  double rho = v[13];
  double vel = v[14];
  double p = (gas_gamma - 1.0) * rho;

  double W = 1.0 / sqrt(1.0 - (vel * vel));
  if (vel * vel > 1.0 - pow(1.0, -8.0)) {
    W = 1.0 / sqrt(pow(1.0, -8.0));
  }

  double Etot = ((rho + p) * (W * W)) - p;
  double mom = (rho + p) * vel * (W * W);

  sout[0] = 2.0 * a_dt * exp_2a;

  sout[1] = a_dx_dx + (b_dt * b_dt) - (b_dx * b_dx) - (c_dt * c_dt) + (c_dx * c_dx) - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p)));
  sout[2] = a_dt_dx;
  sout[3] = b_dx_dx - (2.0 * (b_dt * b_dt)) + (2.0 * (b_dx * b_dx)) + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p)));
  sout[4] = b_dt_dx;
  sout[5] = c_dx_dx - (2.0 * ((b_dt * c_dt) - (b_dx * c_dx)));
  sout[6] = c_dt_dx;

  sout[7] = (2.0 * (b_dt * b_dt_dx)) - (2.0 * (b_dx * b_dx_dx)) - (2.0 * (c_dt * c_dt_dx)) + (2.0 * (c_dx * c_dx_dx));
  sout[8] = 0.0;
  sout[9] = -(4.0 * (b_dt * b_dt_dx)) + (4.0 * (b_dx * b_dx_dx));
  sout[10] = 0.0;
  sout[11] = -2.0 * ((b_dt * c_dt_dx) - (b_dx * c_dx_dx) + (b_dt_dx * c_dt) - (b_dx_dx * c_dx));
  sout[12] = 0.0;

  sout[13] = (-Etot * (a_dt + (2.0 * b_dt))) - (2.0 * mom * (a_dx + b_dx)) - (((mom * vel) + p) * a_dt) - (2.0 * p * b_dt);
  sout[14] = (-Etot * a_dx) - (2.0 * mom * (a_dt + b_dt)) - (((mom * vel) + p) * (a_dx + (2.0 * b_dx))) + (2.0 * p * b_dx);
}

void
gkyl_gr_medium_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_medium *gr_medium = container_of(base->on_dev, struct wv_gr_medium, eqn);
    gkyl_cu_free(gr_medium);
  }

  struct wv_gr_medium *gr_medium = container_of(base, struct wv_gr_medium, eqn);
  gkyl_free(gr_medium);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_medium_new(double gas_gamma, double kappa, bool use_gpu)
{
  return gkyl_wv_gr_medium_inew(&(struct gkyl_wv_gr_medium_inp) {
      .gas_gamma = gas_gamma,
      .kappa = kappa,
      .rp_type = WV_GR_MEDIUM_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_medium_inew(const struct gkyl_wv_gr_medium_inp* inp)
{
  struct wv_gr_medium *gr_medium = gkyl_malloc(sizeof(struct wv_gr_medium));

  gr_medium->eqn.type = GKYL_EQN_GR_MEDIUM;
  gr_medium->eqn.num_equations = 15;
  gr_medium->eqn.num_diag = 15;

  gr_medium->gas_gamma = inp->gas_gamma;
  gr_medium->kappa = inp->kappa;

  if (inp->rp_type == WV_GR_MEDIUM_RP_LAX) {
    gr_medium->eqn.num_waves = 2;
    gr_medium->eqn.waves_func = wave_lax_l;
    gr_medium->eqn.qfluct_func = qfluct_lax_l;
  }

  gr_medium->eqn.flux_jump = flux_jump;
  gr_medium->eqn.check_inv_func = check_inv;
  gr_medium->eqn.max_speed_func = max_speed;
  gr_medium->eqn.rotate_to_local_func = rot_to_local;
  gr_medium->eqn.rotate_to_global_func = rot_to_global;

  gr_medium->eqn.wall_bc_func = gr_medium_wall;
  
  gr_medium->eqn.cons_to_riem = cons_to_riem;
  gr_medium->eqn.riem_to_cons = riem_to_cons;

  gr_medium->eqn.cons_to_diag = gr_medium_cons_to_diag;

  gr_medium->eqn.source_func = gr_medium_source;

  gr_medium->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_medium->eqn.flags);
  gr_medium->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_medium_free);
  gr_medium->eqn.on_dev = &gr_medium->eqn; // On the CPU, the equation object points to itself.

  gr_medium->eqn.embed_geo = NULL;

  return &gr_medium->eqn;
}

double
gkyl_wv_gr_medium_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;

  return gas_gamma;
}

double
gkyl_wv_gr_medium_kappa(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double kappa = gr_medium->kappa;

  return kappa;
}
