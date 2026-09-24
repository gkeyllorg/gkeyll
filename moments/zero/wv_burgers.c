// ** Fully formally-verified solvers for the inviscid Burgers' equation **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_inviscid_burgers_lax_hyperbolicity.rkt **
// ** Proof of CFL stability: ../proofs/finite_volume/proof_inviscid_burgers_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function: ../proofs/finite_volume/proof_inviscid_burgers_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_inviscid_burgers_roe_hyperbolicity.rkt **
// ** Proof of flux conservation (jump continuity): ../proofs/finite_volume/proof_inviscid_burgers_roe_flux_conservation.rkt **

#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_burgers.h>
#include <gkyl_wv_burgers_priv.h>

static inline double
gkyl_burgers_max_abs_speed(const double* q)
{
  return fabs(q[0]);
}

void
gkyl_burgers_flux(const double* q, double* flux)
{
  flux[0] = (0.5 * q[0] * q[0]);
}

void
gkyl_burgers_flux_deriv(const double* q, double* flux_deriv)
{
  flux_deriv[0] = q[0];
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
burgers_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  ghost[0] = skin[0];
}

static void
burgers_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
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
  double sl = gkyl_burgers_max_abs_speed(ql);
  double sr = gkyl_burgers_max_abs_speed(qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double));
  double *fr = gkyl_malloc(sizeof(double));
  gkyl_burgers_flux(ql, fl);
  gkyl_burgers_flux(qr, fr);

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
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_burgers *burgers = container_of(eqn, struct wv_burgers, eqn);
  
  double *fl_deriv = gkyl_malloc(sizeof(double));
  double *fr_deriv = gkyl_malloc(sizeof(double));
  gkyl_burgers_flux_deriv(ql, fl_deriv);
  gkyl_burgers_flux_deriv(qr, fr_deriv);

  double a_roe = 0.5 * (fl_deriv[0] + fr_deriv[0]);

  double *w0 = &waves[0];
  w0[0] = delta[0];

  s[0] = a_roe;

  gkyl_free(fl_deriv);
  gkyl_free(fr_deriv);

  return s[0];
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
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
  double *fr = gkyl_malloc(sizeof(double));
  double *fl = gkyl_malloc(sizeof(double));
  gkyl_burgers_flux(ql, fl);
  gkyl_burgers_flux(qr, fr);

  flux_jump[0] = fr[0] - fl[0];

  double amaxl = gkyl_burgers_max_abs_speed(ql);
  double amaxr = gkyl_burgers_max_abs_speed(qr);
  
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
  return gkyl_burgers_max_abs_speed(q);
}

static inline void
burgers_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  diag[0] = qin[0];
}

static inline void
burgers_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  sout[0] = 0.0;
}

void
gkyl_burgers_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_burgers *burgers = container_of(base->on_dev, struct wv_burgers, eqn);
    gkyl_cu_free(burgers);
  }

  struct wv_burgers *burgers = container_of(base, struct wv_burgers, eqn);
  gkyl_free(burgers);
}

struct gkyl_wv_eqn*
gkyl_wv_burgers_new(bool use_gpu)
{
  return gkyl_wv_burgers_inew(&(struct gkyl_wv_burgers_inp) {
      .rp_type = WV_BURGERS_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_burgers_inew(const struct gkyl_wv_burgers_inp* inp)
{
  struct wv_burgers *burgers = gkyl_malloc(sizeof(struct wv_burgers));

  burgers->eqn.type = GKYL_EQN_BURGERS;
  burgers->eqn.num_equations = 1;
  burgers->eqn.num_diag = 1;

  if (inp->rp_type == WV_BURGERS_RP_ROE) {
    burgers->eqn.num_waves = 1;
    burgers->eqn.waves_func = wave_roe_l;
    burgers->eqn.qfluct_func = qfluct_roe_l;
  }
  else if (inp->rp_type == WV_BURGERS_RP_LAX) {
    burgers->eqn.num_waves = 2;
    burgers->eqn.waves_func = wave_lax_l;
    burgers->eqn.qfluct_func = qfluct_lax_l;
  }

  burgers->eqn.flux_jump = flux_jump;
  burgers->eqn.check_inv_func = check_inv;
  burgers->eqn.max_speed_func = max_speed;
  burgers->eqn.rotate_to_local_func = rot_to_local;
  burgers->eqn.rotate_to_global_func = rot_to_global;
  
  burgers->eqn.wall_bc_func = burgers_wall;
  burgers->eqn.no_slip_bc_func = burgers_no_slip;

  burgers->eqn.cons_to_riem = cons_to_riem;
  burgers->eqn.riem_to_cons = riem_to_cons;

  burgers->eqn.cons_to_diag = burgers_cons_to_diag;

  burgers->eqn.source_func = burgers_source;

  burgers->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(burgers->eqn.flags);
  burgers->eqn.ref_count = gkyl_ref_count_init(gkyl_burgers_free);
  burgers->eqn.on_dev = &burgers->eqn; // On the CPU, the equation object points to itself.

  burgers->eqn.embed_geo = NULL;

  return &burgers->eqn;
}
