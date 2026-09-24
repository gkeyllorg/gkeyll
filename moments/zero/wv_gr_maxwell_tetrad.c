#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_maxwell_tetrad.h>
#include <gkyl_wv_gr_maxwell_tetrad_priv.h>

void
gkyl_gr_maxwell_tetrad_flux(double light_speed, double e_fact, double b_fact, const double q[26], double flux[26])
{
  double Ex = q[0], Ey = q[1], Ez = q[2];
  double Bx = q[3], By = q[4], Bz = q[5];

  double phi = q[6];
  double psi = q[7];

  bool in_excision_region = false;
  if (q[21] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    flux[0] = e_fact * (light_speed * light_speed) * phi;
    flux[1] = (light_speed * light_speed) * Bz;
    flux[2] = -(light_speed * light_speed) * By;
    flux[3] = b_fact * psi;
    flux[4] = -Ez;
    flux[5] = Ey;
    flux[6] = e_fact * Ex;
    flux[7] = b_fact * (light_speed * light_speed) * Bx;

    for (int i = 8; i < 26; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 26; i++) {
      flux[i] = 0.0;
    }
  }
}

void
gkyl_gr_maxwell_tetrad_flux_correction(double light_speed, double e_fact, double b_fact, const double q[26], const double flux_sr[26], double flux_gr[26])
{
  // The flux transformation is _almost_ purely geometrical, but requires knowledge of Ex and Bx for hyperbolic divergence cleaning.
  double Ex = q[0];
  double Bx = q[3];

  double lapse = q[8];
  double shift_x = q[9];
  double shift_y = q[10];
  double shift_z = q[11];

  bool in_excision_region = false;
  if (q[21] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    flux_gr[0] = flux_sr[0];
    flux_gr[1] = (lapse * flux_sr[1]) - (light_speed * light_speed) * ((shift_x * flux_sr[5]) - (shift_y * Ex));
    flux_gr[2] = (lapse * flux_sr[2]) + (light_speed * light_speed) * ((shift_x * flux_sr[4]) + (shift_z * Ex));
    flux_gr[3] = flux_sr[3];
    flux_gr[4] = (lapse * flux_sr[4]) + ((shift_x * (flux_sr[2] / (light_speed * light_speed))) + (shift_y * Bx));
    flux_gr[5] = (lapse * flux_sr[5]) - ((shift_x * (flux_sr[1] / (light_speed * light_speed))) - (shift_z * Bx));
    flux_gr[6] = flux_sr[6];
    flux_gr[7] = flux_sr[7];

    for (int i = 8; i < 26; i++) {
      flux_gr[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 26; i++) {
      flux_gr[i] = 0.0;
    }
  }
}

static inline double
gkyl_gr_maxwell_tetrad_max_abs_speed(double light_speed, const double q[26])
{
  bool in_excision_region = false;
  if (q[21] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double lapse = q[8];

    double spatial_metric[3][3];
    spatial_metric[0][0] = q[12]; spatial_metric[0][1] = q[13]; spatial_metric[0][2] = q[14];
    spatial_metric[1][0] = q[15]; spatial_metric[1][1] = q[16]; spatial_metric[1][2] = q[17];
    spatial_metric[2][0] = q[18]; spatial_metric[2][1] = q[19]; spatial_metric[2][2] = q[20];

    double spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    return light_speed * sqrt(spatial_metric_det) * lapse;
  }
  else {
    return pow(10.0, -8.0);
  }
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 26; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 26; i++) {
    qout[i] = win[i];
  }
}

static void
gr_maxwell_tetrad_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  // Zero tangent for the electric field.
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];

  // Zero normal for the magnetic field.
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
  ghost[5] = skin[5];

  // Correction potentials.
  ghost[6] = -skin[6];
  ghost[7] = skin[7];

  for (int i = 8; i < 26; i++) {
    ghost[i] = skin[i];
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
 // Rotate electric field vector to local coordinates.
  qlocal[0] = (qglobal[0] * norm[0]) + (qglobal[1] * norm[1]) + (qglobal[2] * norm[2]);
  qlocal[1] = (qglobal[0] * tau1[0]) + (qglobal[1] * tau1[1]) + (qglobal[2] * tau1[2]);
  qlocal[2] = (qglobal[0] * tau2[0]) + (qglobal[1] * tau2[1]) + (qglobal[2] * tau2[2]);

  // Rotate magnetic field vector to local coordinates.
  qlocal[3] = (qglobal[3] * norm[0]) + (qglobal[4] * norm[1]) + (qglobal[5] * norm[2]);
  qlocal[4] = (qglobal[3] * tau1[0]) + (qglobal[4] * tau1[1]) + (qglobal[5] * tau1[2]);
  qlocal[5] = (qglobal[3] * tau2[0]) + (qglobal[4] * tau2[1]) + (qglobal[5] * tau2[2]);

  // Correction potentials are scalars (so remain unchanged).
  qlocal[6] = qglobal[6];
  qlocal[7] = qglobal[7];

  // Lapse function is a scalar (so remains unchanged).
  qlocal[8] = qglobal[8];

  // Rotate shift vector to local coordinates.
  qlocal[9] = (qglobal[9] * norm[0]) + (qglobal[10] * norm[1]) + (qglobal[11] * norm[2]);
  qlocal[10] = (qglobal[9] * tau1[0]) + (qglobal[10] * tau1[1]) + (qglobal[11] * tau1[2]);
  qlocal[11] = (qglobal[9] * tau2[0]) + (qglobal[10] * tau2[1]) + (qglobal[11] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[12] * norm[0]) + (qglobal[13] * norm[1]) + (qglobal[14] * norm[2]);
  r1[1] = (qglobal[12] * tau1[0]) + (qglobal[13] * tau1[1]) + (qglobal[14] * tau1[2]);
  r1[2] = (qglobal[12] * tau2[0]) + (qglobal[13] * tau2[1]) + (qglobal[14] * tau2[2]);

  r2[0] = (qglobal[15] * norm[0]) + (qglobal[16] * norm[1]) + (qglobal[17] * norm[2]);
  r2[1] = (qglobal[15] * tau1[0]) + (qglobal[16] * tau1[1]) + (qglobal[17] * tau1[2]);
  r2[2] = (qglobal[15] * tau2[0]) + (qglobal[16] * tau2[1]) + (qglobal[17] * tau2[2]);

  r3[0] = (qglobal[18] * norm[0]) + (qglobal[19] * norm[1]) + (qglobal[20] * norm[2]);
  r3[1] = (qglobal[18] * tau1[0]) + (qglobal[19] * tau1[1]) + (qglobal[20] * tau1[2]);
  r3[2] = (qglobal[18] * tau2[0]) + (qglobal[19] * tau2[1]) + (qglobal[20] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * norm[1]) + (r3[0] * norm[2]);
  v1[1] = (r1[0] * tau1[0]) + (r2[0] * tau1[1]) + (r3[0] * tau1[2]);
  v1[2] = (r1[0] * tau2[0]) + (r2[0] * tau2[1]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * norm[1]) + (r3[1] * norm[2]);
  v2[1] = (r1[1] * tau1[0]) + (r2[1] * tau1[1]) + (r3[1] * tau1[2]);
  v2[2] = (r1[1] * tau2[0]) + (r2[1] * tau2[1]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * norm[1]) + (r3[2] * norm[2]);
  v3[1] = (r1[2] * tau1[0]) + (r2[2] * tau1[1]) + (r3[2] * tau1[2]);
  v3[2] = (r1[2] * tau2[0]) + (r2[2] * tau2[1]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to local coordinates.
  qlocal[12] = v1[0]; qlocal[13] = v1[1]; qlocal[14] = v1[2];
  qlocal[15] = v2[0]; qlocal[16] = v2[1]; qlocal[17] = v2[2];
  qlocal[18] = v3[0]; qlocal[19] = v3[1]; qlocal[20] = v3[2];

  // Excision parameter is a scalar (so remains unchanged).
  qlocal[21] = qglobal[21];

  qlocal[22] = qglobal[22];
  qlocal[23] = (qglobal[23] * norm[0]) + (qglobal[24] * norm[1]) + (qglobal[25] * norm[2]);
  qlocal[24] = (qglobal[23] * tau1[0]) + (qglobal[24] * tau1[1]) + (qglobal[25] * tau1[2]);
  qlocal[25] = (qglobal[23] * tau2[0]) + (qglobal[24] * tau2[1]) + (qglobal[25] * tau2[2]);
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  // Rotate electric field vector to global coordinates.
  qglobal[0] = (qlocal[0] * norm[0]) + (qlocal[1] * tau1[0]) + (qlocal[2] * tau2[0]);
  qglobal[1] = (qlocal[0] * norm[1]) + (qlocal[1] * tau1[1]) + (qlocal[2] * tau2[1]);
  qglobal[2] = (qlocal[0] * norm[2]) + (qlocal[1] * tau1[2]) + (qlocal[2] * tau2[2]);

  // Rotate magnetic field vector to global coordinates.
  qglobal[3] = (qlocal[3] * norm[0]) + (qlocal[4] * tau1[0]) + (qlocal[5] * tau2[0]);
  qglobal[4] = (qlocal[3] * norm[1]) + (qlocal[4] * tau1[1]) + (qlocal[5] * tau2[1]);
  qglobal[5] = (qlocal[3] * norm[2]) + (qlocal[4] * tau1[2]) + (qlocal[5] * tau2[2]);

  // Correction potentials are scalars (so remain unchanged).
  qglobal[6] = qlocal[6];
  qglobal[7] = qlocal[7];

  // Lapse function is a scalar (so remains unchanged).
  qglobal[8] = qlocal[8];

  // Rotate shift vector to global coordinates.
  qglobal[9] = (qlocal[9] * norm[0]) + (qlocal[10] * tau1[0]) + (qlocal[11] * tau2[0]);
  qglobal[10] = (qlocal[9] * norm[1]) + (qlocal[10] * tau1[1]) + (qlocal[11] * tau2[1]);
  qglobal[11] = (qlocal[9] * norm[2]) + (qlocal[10] * tau1[2]) + (qlocal[11] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[12] * norm[0]) + (qlocal[13] * tau1[0]) + (qlocal[14] * tau2[0]);
  r1[1] = (qlocal[12] * norm[1]) + (qlocal[13] * tau1[1]) + (qlocal[14] * tau2[1]);
  r1[2] = (qlocal[12] * norm[2]) + (qlocal[13] * tau1[2]) + (qlocal[14] * tau2[2]);

  r2[0] = (qlocal[15] * norm[0]) + (qlocal[16] * tau1[0]) + (qlocal[17] * tau2[0]);
  r2[1] = (qlocal[15] * norm[1]) + (qlocal[16] * tau1[1]) + (qlocal[17] * tau2[1]);
  r2[2] = (qlocal[15] * norm[2]) + (qlocal[16] * tau1[2]) + (qlocal[17] * tau2[2]);

  r3[0] = (qlocal[18] * norm[0]) + (qlocal[19] * tau1[0]) + (qlocal[20] * tau2[0]);
  r3[1] = (qlocal[18] * norm[1]) + (qlocal[19] * tau1[1]) + (qlocal[20] * tau2[1]);
  r3[2] = (qlocal[18] * norm[2]) + (qlocal[19] * tau1[2]) + (qlocal[20] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * tau1[0]) + (r3[0] * tau2[0]);
  v1[1] = (r1[0] * norm[1]) + (r2[0] * tau1[1]) + (r3[0] * tau2[1]);
  v1[2] = (r1[0] * norm[2]) + (r2[0] * tau1[2]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * tau1[0]) + (r3[1] * tau2[0]);
  v2[1] = (r1[1] * norm[1]) + (r2[1] * tau1[1]) + (r3[1] * tau2[1]);
  v2[2] = (r1[1] * norm[2]) + (r2[1] * tau1[2]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * tau1[0]) + (r3[2] * tau2[0]);
  v3[1] = (r1[2] * norm[1]) + (r2[2] * tau1[1]) + (r3[2] * tau2[1]);
  v3[2] = (r1[2] * norm[2]) + (r2[2] * tau1[2]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to global coordinates.
  qglobal[12] = v1[0]; qglobal[13] = v1[1]; qglobal[14] = v1[2];
  qglobal[15] = v2[0]; qglobal[16] = v2[1]; qglobal[17] = v2[2];
  qglobal[18] = v3[0]; qglobal[19] = v3[1]; qglobal[20] = v3[2];

  // Excision parameter is a scalar (so remains unchanged).
  qglobal[21] = qlocal[21];

  qglobal[22] = qlocal[22];
  qglobal[23] = (qlocal[23] * norm[0]) + (qlocal[24] * tau1[0]) + (qlocal[25] * tau2[0]);
  qglobal[24] = (qlocal[23] * norm[1]) + (qlocal[24] * tau1[1]) + (qlocal[25] * tau2[1]);
  qglobal[25] = (qlocal[23] * norm[2]) + (qlocal[24] * tau1[2]) + (qlocal[25] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  double light_speed = gr_maxwell_tetrad->light_speed;
  double e_fact = gr_maxwell_tetrad->e_fact;
  double b_fact = gr_maxwell_tetrad->b_fact;

  double sl = gkyl_gr_maxwell_tetrad_max_abs_speed(light_speed, ql);
  double sr = gkyl_gr_maxwell_tetrad_max_abs_speed(light_speed, qr);
  double amax = fmax(sl, sr);

  double fl_sr[26], fr_sr[26];
  gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, ql, fl_sr);
  gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, qr, fr_sr);

  double fl_gr[26], fr_gr[26];
  gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, qr, fr_sr, fr_gr);

  bool in_excision_region_l = false;
  if (ql[21] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[21] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[26];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 26; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr_gr[i] - fl_gr[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr_gr[i] - fl_gr[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 26; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[26];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 26; i++) {
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
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  double light_speed = gr_maxwell_tetrad->light_speed;
  double e_fact = gr_maxwell_tetrad->e_fact;
  double b_fact = gr_maxwell_tetrad->b_fact;

  double a1 = 0.5 * (delta[3] - ((1.0 / light_speed) * delta[7]));
  double a2 = 0.5 * (delta[3] + ((1.0 / light_speed) * delta[7]));
  double a3 = 0.5 * (delta[0] - (light_speed * delta[6]));
  double a4 = 0.5 * (delta[0] + (light_speed * delta[6]));
  double a5 = 0.5 * (delta[1] - (light_speed * delta[5]));
  double a6 = 0.5 * ((light_speed * delta[4]) + delta[2]);
  double a7 = 0.5 * ((light_speed * delta[5]) + delta[1]);
  double a8 = 0.5 * (delta[2] - (light_speed * delta[4]));

  for (int i = 0; i < 26 * 6; i++) {
    waves[i] = 0.0;
  }

  double *wv ;
  wv = &waves[0 * 26];
  wv[3] = a1;
  wv[7] = -a1 * light_speed;
  s[0] = -light_speed * b_fact;

  wv = &waves[1 * 26];
  wv[3] = a2;
  wv[7] = a2 * light_speed;
  s[1] = light_speed * b_fact;

  wv = &waves[2 * 26];
  wv[0] = a3;
  wv[6] = -a3 * (1.0 / light_speed);
  s[2] = -light_speed * e_fact;

  wv = &waves[3 * 26];
  wv[0] = a4;
  wv[6] = a4 * (1.0 / light_speed);
  s[3] = light_speed * e_fact;

  wv = &waves[4 * 26];
  wv[1] = a5;
  wv[2] = a6;
  wv[4] = a6 * (1.0 / light_speed);
  wv[5] = -a5 * (1.0 / light_speed);
  s[4] = -light_speed;

  wv = &waves[5 * 26];
  wv[1] = a7;
  wv[2] = a8;
  wv[4] = -a8 * (1.0 / light_speed);
  wv[5] = a7 * (1.0 / light_speed);
  s[5] = light_speed;

  return light_speed;
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0 * 26], *w1 = &waves[1 * 26], *w2 = &waves[2 * 26];
  const double *w3 = &waves[3 * 26], *w4 = &waves[4 * 26], *w5 = &waves[5 * 26];

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);

  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  double s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  for (int i = 0; i < 26; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]) + (s3m * w3[i]) + (s4m * w4[i]) + (s5m * w5[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]) + (s3p * w3[i]) + (s4p * w4[i]) + (s5p * w5[i]);
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
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  const double light_speed = gr_maxwell_tetrad->light_speed;
  const double e_fact = gr_maxwell_tetrad->e_fact;
  const double b_fact = gr_maxwell_tetrad->b_fact;

  double fr_sr[26], fl_sr[26];
  gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, ql, fl_sr);
  gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, qr, fr_sr);

  double fr_gr[26], fl_gr[26];
  gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, qr, fr_sr, fr_gr);

  bool in_excision_region_l = false;
  if (ql[21] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[21] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 26; m++) {
      flux_jump[m] = fr_gr[m] - fl_gr[m];
    }
  }
  else {
    for (int m = 0; m < 26; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_maxwell_tetrad_max_abs_speed(light_speed, ql);
  double amaxr = gkyl_gr_maxwell_tetrad_max_abs_speed(light_speed, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  // No invalid states for general relativistic Maxwell in the tetrad basis.
  return true;
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  double light_speed = gr_maxwell_tetrad->light_speed;

  return gkyl_gr_maxwell_tetrad_max_abs_speed(light_speed, q);
}

static inline void
gr_maxwell_tetrad_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 8; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_maxwell_tetrad_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 26; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_gr_maxwell_tetrad_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(base->on_dev, struct wv_gr_maxwell_tetrad, eqn);
    gkyl_cu_free(gr_maxwell_tetrad);
  }

  struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(base, struct wv_gr_maxwell_tetrad, eqn);
  gkyl_free(gr_maxwell_tetrad);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_tetrad_new(double light_speed, double e_fact, double b_fact, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_maxwell_tetrad_inew(&(struct gkyl_wv_gr_maxwell_tetrad_inp) {
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_MAXWELL_TETRAD_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_tetrad_inew(const struct gkyl_wv_gr_maxwell_tetrad_inp* inp)
{
  struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = gkyl_malloc(sizeof(struct wv_gr_maxwell_tetrad));

  gr_maxwell_tetrad->eqn.type = GKYL_EQN_GR_MAXWELL_TETRAD;
  gr_maxwell_tetrad->eqn.num_equations = 26;
  gr_maxwell_tetrad->eqn.num_diag = 8;

  gr_maxwell_tetrad->light_speed = inp->light_speed;
  gr_maxwell_tetrad->e_fact = inp->e_fact;
  gr_maxwell_tetrad->b_fact = inp->b_fact;
  gr_maxwell_tetrad->spacetime_gauge = inp->spacetime_gauge;
  gr_maxwell_tetrad->reinit_freq = inp->reinit_freq;
  gr_maxwell_tetrad->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_MAXWELL_TETRAD_RP_LAX) {
    gr_maxwell_tetrad->eqn.num_waves = 2;
    gr_maxwell_tetrad->eqn.waves_func = wave_lax_l;
    gr_maxwell_tetrad->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_MAXWELL_TETRAD_RP_ROE) {
    gr_maxwell_tetrad->eqn.num_waves = 6;
    gr_maxwell_tetrad->eqn.waves_func = wave_roe_l;
    gr_maxwell_tetrad->eqn.qfluct_func = qfluct_roe_l;
  }

  gr_maxwell_tetrad->eqn.flux_jump = flux_jump;
  gr_maxwell_tetrad->eqn.check_inv_func = check_inv;
  gr_maxwell_tetrad->eqn.max_speed_func = max_speed;
  gr_maxwell_tetrad->eqn.rotate_to_local_func = rot_to_local;
  gr_maxwell_tetrad->eqn.rotate_to_global_func = rot_to_global;

  gr_maxwell_tetrad->eqn.wall_bc_func = gr_maxwell_tetrad_wall;

  gr_maxwell_tetrad->eqn.cons_to_riem = cons_to_riem;
  gr_maxwell_tetrad->eqn.riem_to_cons = riem_to_cons;

  gr_maxwell_tetrad->eqn.cons_to_diag = gr_maxwell_tetrad_cons_to_diag;

  gr_maxwell_tetrad->eqn.source_func = gr_maxwell_tetrad_source;

  gr_maxwell_tetrad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_maxwell_tetrad->eqn.flags);
  gr_maxwell_tetrad->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_maxwell_tetrad_free);
  gr_maxwell_tetrad->eqn.on_dev = &gr_maxwell_tetrad->eqn; // On the CPU, the equation object points to itself.

  gr_maxwell_tetrad->eqn.embed_geo = NULL;

  return &gr_maxwell_tetrad->eqn;
}

double
gkyl_wv_gr_maxwell_tetrad_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  double light_speed = gr_maxwell_tetrad->light_speed;

  return light_speed;
}

double
gkyl_wv_gr_maxwell_tetrad_e_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  double e_fact = gr_maxwell_tetrad->e_fact;

  return e_fact;
}

double
gkyl_wv_gr_maxwell_tetrad_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  double b_fact = gr_maxwell_tetrad->b_fact;

  return b_fact;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_maxwell_tetrad_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_maxwell_tetrad->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_maxwell_tetrad_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  int reinit_freq = gr_maxwell_tetrad->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_maxwell_tetrad_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_maxwell_tetrad->spacetime;

  return spacetime;
}
