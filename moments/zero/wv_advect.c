// ** Fully formally-verified solvers for the linear advection equation **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_linear_advection_lax_hyperbolicity.rkt **
// ** Proof of CFL stability: ../proofs/finite_volume/proof_linear_advection_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function: ../proofs/finite_volume/proof_linear_advection_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_linear_advection_roe_hyperbolicity.rkt **
// ** Proof of flux conservation (jump continuity): ../proofs/finite_volume/proof_linear_advection_roe_flux_conservation.rkt **

#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_advect.h>
#include <gkyl_wv_advect_priv.h>

static inline double
gkyl_advect_max_abs_speed(double a, const double* q)
{
  return fabs(a);
}

void
gkyl_advect_flux(double a, const double* q, double* flux)
{
  flux[0] = (a * q[0]);
}

void
gkyl_advect_flux_deriv(double a, const double* q, double* flux_deriv)
{
  flux_deriv[0] = a;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  wout[0] = qin[0];
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  qout[0] = win[0];
}

static void
advect_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  ghost[0] = skin[0];
}

static void
advect_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  ghost[0] = skin[0];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);
  double a = advect->a; // Advection speed.

  double sl = gkyl_advect_max_abs_speed(a, ql);
  double sr = gkyl_advect_max_abs_speed(a, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double));
  double *fr = gkyl_malloc(sizeof(double));
  gkyl_advect_flux(a, ql, fl);
  gkyl_advect_flux(a, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[1];
  w0[0] = 0.5 * ((qr[0] - ql[0]) - (fr[0] - fl[0]) / amax);
  w1[0] = 0.5 * ((qr[0] - ql[0]) + (fr[0] - fl[0]) / amax);

  s[0] = -amax;
  s[1] = amax;

  gkyl_free(fl);
  gkyl_free(fr);

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[1];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  amdq[0] = (s0m * w0[0]) + (s1m * w1[0]);
  apdq[0] = (s0p * w0[0]) + (s1p * w1[0]);
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
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);
  double a = advect->a; // Additional simulation parameter.

  double *fl_deriv = gkyl_malloc(sizeof(double));
  double *fr_deriv = gkyl_malloc(sizeof(double));
  gkyl_advect_flux_deriv(a, ql, fl_deriv);
  gkyl_advect_flux_deriv(a, qr, fr_deriv);

  double a_roe = 0.5 * (fl_deriv[0] + fr_deriv[0]);

  double *w0 = &waves[0];
  w0[0] = delta[0];

  s[0] = a_roe;

  gkyl_free(fl_deriv);
  gkyl_free(fr_deriv);

  return s[0];
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0];

  if (s[0] < 0.0) {
    amdq[0] = s[0] * w0[0];
    apdq[0] = 0.0;
  }
  else {
    amdq[0] = 0.0;
    apdq[0] = s[0] * w0[0];
  }
}

static double
wave(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_roe(eqn, delta, ql, qr, phil, phir, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_roe(eqn, ql, qr, phil, phir, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);
  double a = advect->a; // Advection speed.

  double *fr = gkyl_malloc(sizeof(double));
  double *fl = gkyl_malloc(sizeof(double));
  gkyl_advect_flux(a, ql, fl);
  gkyl_advect_flux(a, qr, fr);

  flux_jump[0] = fr[0] - fl[0];

  double amaxl = gkyl_advect_max_abs_speed(a, ql);
  double amaxr = gkyl_advect_max_abs_speed(a, qr);
  
  gkyl_free(fr);
  gkyl_free(fl);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  return true; // All states are assumed to be valid.
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);
  double a = advect->a; // Advection speed.

  return gkyl_advect_max_abs_speed(a, q);
}

static inline void
advect_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  diag[0] = qin[0];
}

static inline void
advect_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  sout[0] = 0.0;
}

void
gkyl_advect_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_advect *advect = container_of(base->on_dev, struct wv_advect, eqn);
    gkyl_cu_free(advect);
  }

  struct wv_advect *advect = container_of(base, struct wv_advect, eqn);
  gkyl_free(advect);
}

struct gkyl_wv_eqn*
gkyl_wv_advect_new(double a, bool use_gpu)
{
  return gkyl_wv_advect_inew(&(struct gkyl_wv_advect_inp) {
      .a = a,
      .rp_type = WV_ADVECT_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_advect_inew(const struct gkyl_wv_advect_inp* inp)
{
  struct wv_advect *advect = gkyl_malloc(sizeof(struct wv_advect));

  advect->eqn.type = GKYL_EQN_ADVECTION;
  advect->eqn.num_equations = 1;
  advect->eqn.num_diag = 1;

  advect->a = inp->a;

  if (inp->rp_type == WV_ADVECT_RP_ROE) {
    advect->eqn.num_waves = 1;
    advect->eqn.waves_func = wave;
    advect->eqn.qfluct_func = qfluct;
  }
  else if (inp->rp_type == WV_ADVECT_RP_LAX) {
    advect->eqn.num_waves = 2;
    advect->eqn.waves_func = wave_lax_l;
    advect->eqn.qfluct_func = qfluct_lax_l;
  }

  advect->eqn.flux_jump = flux_jump;
  advect->eqn.check_inv_func = check_inv;
  advect->eqn.max_speed_func = max_speed;
  advect->eqn.rotate_to_local_func = rot_to_local;
  advect->eqn.rotate_to_global_func = rot_to_global;
  
  advect->eqn.wall_bc_func = advect_wall;
  advect->eqn.no_slip_bc_func = advect_no_slip;

  advect->eqn.cons_to_riem = cons_to_riem;
  advect->eqn.riem_to_cons = riem_to_cons;

  advect->eqn.cons_to_diag = advect_cons_to_diag;

  advect->eqn.source_func = advect_source;

  advect->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(advect->eqn.flags);
  advect->eqn.ref_count = gkyl_ref_count_init(gkyl_advect_free);
  advect->eqn.on_dev = &advect->eqn; // On the CPU, the equation object points to itself.

  advect->eqn.embed_geo = NULL;

  return &advect->eqn;
}
