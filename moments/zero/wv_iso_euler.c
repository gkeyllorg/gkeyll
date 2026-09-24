// ** Partially formally-verified solvers for the isothermal Euler equations **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation (rho and mom_x components): ../proofs/finite_volume/proof_isothermal_euler_mom_x_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (rho and mom_x components): ../proofs/finite_volume/proof_isothermal_euler_mom_x_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_strict_hyperbolicity.rkt **
// ** Proof of CFL stability (rho and mom_x components): ../proofs/finite_volume/proof_isothermal_euler_mom_x_lax_cfl_stability.rkt **
// ** Proof of CFL stability (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (rho and mom_x components): NOT PROVEN **
// ** Proof of local Lipschitz continuity of discrete flux function (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation (rho and mom_x components): NOT PROVEN **
// ** Proof of hyperbolicity preservation (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (rho and mom_x components): NOT PROVEN **
// ** Proof of strict hyperbolicity preservation (mom_y and mom_z components): NOT PROVEN **
// ** Proof of flux conservation (jump continuity, rho and mom_x components): NOT PROVEN **
// ** Proof of flux conservation (jump continuity, mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_flux_conservation.rkt **

#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_iso_euler_priv.h>

static inline double
gkyl_iso_euler_max_abs_speed(double vt, const double* q)
{
  return fmax(fabs(((q[1] / q[0]) - vt)), fabs(((q[1] / q[0]) + vt)));
}

void
gkyl_iso_euler_flux(double vt, const double* q, double* flux)
{
  flux[0] = q[1];
  flux[1] = (((q[1] * q[1]) / q[0]) + (q[0] * vt * vt));
  flux[2] = (q[2] * (q[1] / q[0]));
  flux[3] = (q[3] * (q[1] / q[0]));
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 4; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 4; i++) {
    qout[i] = win[i];
  }
}

static void
iso_euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  // Copy density.
  ghost[0] = skin[0];

  // Zero normal for the momentum.
  ghost[1] = -skin[1];
  ghost[2] = skin[2];
  ghost[3] = skin[3];
}

static void
iso_euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  // Copy density.
  ghost[0] = skin[0];

  // Zero normal and tangent for the momentum.
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  // Density is a scalar (so remains unchanged).
  qlocal[0] = qglobal[0];

  // Rotate momentum vector to local coordinates.
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  // Density is a scalar (so remains unchanged).
  qglobal[0] = qlocal[0];

  // Rotate momentum vector to global coordinates.
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt; // Thermal velocity.

  double sl = gkyl_iso_euler_max_abs_speed(vt, ql);
  double sr = gkyl_iso_euler_max_abs_speed(vt, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double) * 4);
  double *fr = gkyl_malloc(sizeof(double) * 4);
  gkyl_iso_euler_flux(vt, ql, fl);
  gkyl_iso_euler_flux(vt, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[4];
  for (int i = 0; i < 4; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[4];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 4; i++) {
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
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt; // Thermal velocity.

  double a0 = (delta[0] * (vt + ((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])))) / vt / 2.0) - (delta[1] / vt / 2.0);
  double a1 = delta[2] - (delta[0] * ((ql[2] * (1.0 / sqrt(ql[0]))) + (qr[2] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0]))));
  double a2 = delta[3] - (delta[0] * ((ql[3] * (1.0 / sqrt(ql[0]))) + (qr[3] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0]))));
  double a3 = (delta[0] * (vt - ((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])))) / vt / 2.0) + (delta[1] / vt / 2.0);

  double *w0 = &waves[0 * 4], *w1 = &waves[1 * 4], *w2 = &waves[2 * 4];
  for (int i = 0; i < 4; i++) {
    w0[i] = 0.0; w1[i] = 0.0; w2[i] = 0.0;
  }

  w0[0] = a0;
  w0[1] = a0 * (((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0]))) - vt);
  w0[2] = a0 * ((ql[2] * (1.0 / sqrt(ql[0]))) + (qr[2] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])));
  w0[3] = a0 * ((ql[3] * (1.0 / sqrt(ql[0]))) + (qr[3] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])));
  s[0] = ((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0]))) - vt;

  w1[0] = 0.0;
  w1[1] = 0.0;
  w1[2] = a1;
  w1[3] = a2;
  s[1] = ((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])));

  w2[0] = a3;
  w2[1] = a3 * (((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0]))) + vt);
  w2[2] = a3 * ((ql[2] * (1.0 / sqrt(ql[0]))) + (qr[2] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])));
  w2[3] = a3 * ((ql[3] * (1.0 / sqrt(ql[0]))) + (qr[3] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])));
  s[2] = ((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0]))) + vt;

  return fabs(((ql[1] * (1.0 / sqrt(ql[0]))) + (qr[1] * (1.0 / sqrt(qr[0])))) * (1.0 / (sqrt(ql[0]) + sqrt(qr[0])))) + vt;
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0 * 4], *w1 = &waves[1 * 4], *w2 = &waves[2 * 4];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 4; i++) {
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
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt; // Thermal velocity.

  double *fr = gkyl_malloc(sizeof(double) * 4);
  double *fl = gkyl_malloc(sizeof(double) * 4);
  gkyl_iso_euler_flux(vt, ql, fl);
  gkyl_iso_euler_flux(vt, qr, fr);

  for (int i = 0; i < 4; i++) {
    flux_jump[i] = fr[i] - fl[i];
  }

  double amaxl = gkyl_iso_euler_max_abs_speed(vt, ql);
  double amaxr = gkyl_iso_euler_max_abs_speed(vt, qr);
  
  gkyl_free(fr);
  gkyl_free(fl);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  return q[0] > 0.0; // Density must be positive.
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt; // Thermal velocity.

  return gkyl_iso_euler_max_abs_speed(vt, q);
}

static inline void
iso_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 4; i++) {
    diag[i] = qin[i];
  }
}

static inline void
iso_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 4; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_iso_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_iso_euler *iso_euler = container_of(base->on_dev, struct wv_iso_euler, eqn);
    gkyl_cu_free(iso_euler);
  }

  struct wv_iso_euler *iso_euler = container_of(base, struct wv_iso_euler, eqn);
  gkyl_free(iso_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_new(double vt, bool use_gpu)
{
  return gkyl_wv_iso_euler_inew(&(struct gkyl_wv_iso_euler_inp) {
      .vt = vt,
      .rp_type = WV_ISO_EULER_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_inew(const struct gkyl_wv_iso_euler_inp* inp)
{
  struct wv_iso_euler *iso_euler = gkyl_malloc(sizeof(struct wv_iso_euler));

  iso_euler->eqn.type = GKYL_EQN_ISO_EULER;
  iso_euler->eqn.num_equations = 4;
  iso_euler->eqn.num_diag = 4;

  iso_euler->vt = inp->vt;

  if (inp->rp_type == WV_ISO_EULER_RP_LAX) {
    iso_euler->eqn.num_waves = 2;
    iso_euler->eqn.waves_func = wave_lax_l;
    iso_euler->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_ISO_EULER_RP_ROE) {
    iso_euler->eqn.num_waves = 3;
    iso_euler->eqn.waves_func = wave_roe_l;
    iso_euler->eqn.qfluct_func = qfluct_roe_l;
  }

  iso_euler->eqn.flux_jump = flux_jump;
  iso_euler->eqn.check_inv_func = check_inv;
  iso_euler->eqn.max_speed_func = max_speed;
  iso_euler->eqn.rotate_to_local_func = rot_to_local;
  iso_euler->eqn.rotate_to_global_func = rot_to_global;
  
  iso_euler->eqn.wall_bc_func = iso_euler_wall;
  iso_euler->eqn.no_slip_bc_func = iso_euler_no_slip;

  iso_euler->eqn.cons_to_riem = cons_to_riem;
  iso_euler->eqn.riem_to_cons = riem_to_cons;

  iso_euler->eqn.cons_to_diag = iso_euler_cons_to_diag;

  iso_euler->eqn.source_func = iso_euler_source;

  iso_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(iso_euler->eqn.flags);
  iso_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_iso_euler_free);
  iso_euler->eqn.on_dev = &iso_euler->eqn; // On the CPU, the equation object points to itself.

  iso_euler->eqn.embed_geo = NULL;

  return &iso_euler->eqn;
}

double
gkyl_wv_iso_euler_vt(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt;

  return vt;
}
