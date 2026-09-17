#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_mat.h>
#include <gkyl_moment_em_coupling_priv.h>
#include <gkyl_sources_explicit_priv.h>
#include <gkyl_sources_implicit_priv.h>

void
pressure_tensor_rotate(double q_over_m, double dt, const double* em, const double* ext_em, double p_tensor_old[6], double p_tensor_rhs[6],
  double p_tensor_new[6])
{
  double Bx = em[3] + ext_em[3];
  double By = em[4] + ext_em[4];
  double Bz = em[5] + ext_em[5];

  double dt1 = 0.5 * dt;
  double dt1_sq = dt1 * dt1;
  double dt1_cu = dt1_sq * dt1;
  double dt1_qu = dt1_cu * dt1;

  double Bx_sq = Bx * Bx;
  double Bx_cu = Bx_sq * Bx;
  double Bx_qu = Bx_cu * Bx;

  double By_sq = By * By;
  double By_cu = By_sq * By;
  double By_qu = By_cu * By;

  double Bz_sq = Bz * Bz;
  double Bz_cu = Bz_sq * Bz;
  double Bz_qu = Bz_cu * Bz;

  double q_over_m_sq = q_over_m * q_over_m;
  double q_over_m_cu = q_over_m_sq * q_over_m;
  double q_over_m_qu = q_over_m_cu * q_over_m;

  double denom = 1.0 + (5.0 * (Bx_sq + By_sq + Bz_sq) * dt1_sq * q_over_m_sq) + (4.0 * (Bx_sq + By_sq + Bz_sq) * (Bx_sq + By_sq + Bz_sq) * dt1_qu * q_over_m_qu);

  p_tensor_new[0] = 2.0 * (p_tensor_rhs[0] + 2.0 * dt1 * (Bz * p_tensor_rhs[1] - By * p_tensor_rhs[2]) * q_over_m + dt1_sq * (5.0 * Bx_sq * p_tensor_rhs[0] +
    2.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + Bz_sq * (3.0 * p_tensor_rhs[0] + 2.0 * p_tensor_rhs[3]) - 4.0 * By * Bz * p_tensor_rhs[4] +
    By_sq * (3.0 * p_tensor_rhs[0] + 2.0 * p_tensor_rhs[5])) * q_over_m_sq + 2.0 * dt1_cu * (4.0 * Bx_sq * (Bz * p_tensor_rhs[1] - By * p_tensor_rhs[2]) -
    (By_sq + Bz_sq) * (-(Bz * p_tensor_rhs[1]) + By * p_tensor_rhs[2]) - 3.0 * Bx * (By_sq * p_tensor_rhs[4] - Bz_sq * p_tensor_rhs[4] + By * Bz *
    (-p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * dt1_qu * (2.0 * Bx_qu * p_tensor_rhs[0] + 4.0 * Bx_cu * (By * p_tensor_rhs[1] +
    Bz * p_tensor_rhs[2]) - 2.0 * Bx * (By_sq + Bz_sq) * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + (By_sq + Bz_sq) * (Bz_sq * (p_tensor_rhs[0] +
    p_tensor_rhs[3]) - 2.0 * By * Bz * p_tensor_rhs[4] + By_sq * (p_tensor_rhs[0] + p_tensor_rhs[5])) + Bx_sq * (4.0 * By * Bz * p_tensor_rhs[4] + By_sq * (3.0 *
    p_tensor_rhs[3] + p_tensor_rhs[5]) + Bz_sq * (p_tensor_rhs[3] + 3.0 * p_tensor_rhs[5]))) * q_over_m_qu) / denom - p_tensor_old[0];

  p_tensor_new[1] = 2.0 * (p_tensor_rhs[1] + dt1 * (Bx * p_tensor_rhs[2] + Bz * (-p_tensor_rhs[0] + p_tensor_rhs[3]) - By * p_tensor_rhs[4]) * q_over_m + dt1_sq *
    (4.0 * Bx_sq * p_tensor_rhs[1] + 4.0 * By_sq * p_tensor_rhs[1] + Bz_sq * p_tensor_rhs[1] + 3.0 * By * Bz * p_tensor_rhs[2] + Bx * (3.0 * Bz * p_tensor_rhs[4] +
    By * (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]))) * q_over_m_sq + dt1_cu * (4.0 * Bx_cu * p_tensor_rhs[2] - 2.0 * Bx * (By_sq + Bz_sq) * 
    p_tensor_rhs[2] + Bz_cu * (-p_tensor_rhs[0] + p_tensor_rhs[3]) - 4.0 * By_cu * p_tensor_rhs[4] + 2.0 * By * Bz_sq * p_tensor_rhs[4] - By_sq * Bz *
    (p_tensor_rhs[0] - 4.0 * p_tensor_rhs[3] + 3.0 * p_tensor_rhs[5]) + Bx_sq * (2.0 * By * p_tensor_rhs[4] + Bz * (-4.0 * p_tensor_rhs[0] + p_tensor_rhs[3] +
    3.0 * p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * Bx * By * dt1_qu * (6.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 6.0 * By * Bz * p_tensor_rhs[4] -
    Bz_sq * (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]) + Bx_sq * (2.0 * p_tensor_rhs[0] - p_tensor_rhs[3] - p_tensor_rhs[5]) - By_sq *
    (p_tensor_rhs[0] - 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m_qu) / denom - p_tensor_old[1];

  p_tensor_new[2] = 2.0 * (p_tensor_rhs[2] + dt1 * (-(Bx * p_tensor_rhs[1]) + Bz * p_tensor_rhs[4] + By * (p_tensor_rhs[0] - p_tensor_rhs[5])) * q_over_m + dt1_sq *
    (3.0 * By * Bz * p_tensor_rhs[1] + 4.0 * Bx_sq * p_tensor_rhs[2] + By_sq * p_tensor_rhs[2] + 4.0 * Bz_sq * p_tensor_rhs[2] + Bx * (3.0 * By * p_tensor_rhs[4] +
    Bz * (p_tensor_rhs[0] - 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_sq + dt1_cu * (-4.0 * Bx_cu * p_tensor_rhs[1] + 2.0 * Bx * (By_sq + Bz_sq) *
    p_tensor_rhs[1] - 2.0 * By_sq * Bz * p_tensor_rhs[4] + 4.0 * Bz_cu * p_tensor_rhs[4] + By * Bz_sq * (p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3] - 4.0 *
    p_tensor_rhs[5]) + By_cu * (p_tensor_rhs[0] - p_tensor_rhs[5]) - Bx_sq * (2.0 * Bz * p_tensor_rhs[4] + By * (-4.0 * p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3] +
    p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * Bx * Bz * dt1_qu * (6.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 6.0 * By * Bz * p_tensor_rhs[4] -
    Bz_sq * (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]) + Bx_sq * (2.0 * p_tensor_rhs[0] - p_tensor_rhs[3] - p_tensor_rhs[5]) - By_sq * (p_tensor_rhs[0] -
    2.0 * p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m_qu) / denom - p_tensor_old[2];
  
  p_tensor_new[3] = 2.0 * (p_tensor_rhs[3] + (-2.0 * Bz * dt1 * p_tensor_rhs[1] + 2.0 * Bx * dt1 * p_tensor_rhs[4]) * q_over_m + dt1_sq * (2.0 * Bx * By *
    p_tensor_rhs[1] + 5.0 * By_sq * p_tensor_rhs[3] + Bz_sq * (2.0 * p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3]) + Bz * (-4.0 * Bx * p_tensor_rhs[2] + 2.0 * By *
    p_tensor_rhs[4]) + Bx_sq * (3.0 * p_tensor_rhs[3] + 2.0 * p_tensor_rhs[5])) * q_over_m_sq + 2.0 * dt1_cu * (Bx_sq * (-(Bz * p_tensor_rhs[1]) + 3.0 * By *
    p_tensor_rhs[2]) - Bz * (4.0 * By_sq * p_tensor_rhs[1] + Bz_sq * p_tensor_rhs[1] + 3.0 * By * Bz * p_tensor_rhs[2]) + Bx_cu * p_tensor_rhs[4] + Bx *
    (4.0 * By_sq * p_tensor_rhs[4] + Bz_sq * p_tensor_rhs[4] + 3.0 * By * Bz * (-p_tensor_rhs[0] + p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * dt1_qu * (-2.0 *
    Bx_cu * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 2.0 * Bx * (2.0 * By_sq - Bz_sq) * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 2.0 * By_qu *
    p_tensor_rhs[3] + Bz_qu * (p_tensor_rhs[0] + p_tensor_rhs[3]) + 4.0 * By_cu * Bz * p_tensor_rhs[4] - 2.0 * By * Bz_cu * p_tensor_rhs[4] + Bx_qu * (p_tensor_rhs[3] +
    p_tensor_rhs[5]) + By_sq * Bz_sq * (p_tensor_rhs[0] + 3.0 * p_tensor_rhs[5]) + Bx_sq * (-2.0 * By * Bz * p_tensor_rhs[4] + By_sq * (3.0 * p_tensor_rhs[0] +
    p_tensor_rhs[5]) + Bz_sq * (p_tensor_rhs[0] + 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_qu) / denom - p_tensor_old[3];

  p_tensor_new[4] = 2.0 * (p_tensor_rhs[4] + dt1 * (By * p_tensor_rhs[1] - Bz * p_tensor_rhs[2] + Bx * (-p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m + dt1_sq *
    (3.0 * Bx * Bz * p_tensor_rhs[1] + Bx_sq * p_tensor_rhs[4] + 4.0 * By_sq * p_tensor_rhs[4] + 4.0 * Bz_sq * p_tensor_rhs[4] + By * (3.0 * Bx * p_tensor_rhs[2] +
    Bz * (-2.0 * p_tensor_rhs[0] + p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_sq + dt1_cu * (4.0 * By_cu * p_tensor_rhs[1] - 2.0 * By * Bz_sq * p_tensor_rhs[1] +
    2.0 * By_sq * Bz * p_tensor_rhs[2] - 4.0 * Bz_cu * p_tensor_rhs[2] + Bx_sq * (-2.0 * By * p_tensor_rhs[1] + 2.0 * Bz * p_tensor_rhs[2]) + Bx_cu * (-p_tensor_rhs[3] +
    p_tensor_rhs[5]) + Bx * (-(Bz_sq * (3.0 * p_tensor_rhs[0] + p_tensor_rhs[3] - 4.0 * p_tensor_rhs[5])) + By_sq * (3.0 * p_tensor_rhs[0] - 4.0 * p_tensor_rhs[3] +
    p_tensor_rhs[5]))) * q_over_m_cu - 2.0 * By * Bz * dt1_qu * (-6.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) - 6.0 * By * Bz * p_tensor_rhs[4] + Bz_sq *
    (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]) + By_sq * (p_tensor_rhs[0] - 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5]) + Bx_sq * (-2.0 * p_tensor_rhs[0] +
    p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m_qu) / denom - p_tensor_old[4];

  p_tensor_new[5] = 2.0 * (p_tensor_rhs[5] + 2.0 * dt1 * (By * p_tensor_rhs[2] - Bx * p_tensor_rhs[4]) * q_over_m + dt1_sq * (2.0 * Bx * Bz * p_tensor_rhs[2] + By *
    (-4.0 * Bx * p_tensor_rhs[1] + 2 * Bz * p_tensor_rhs[4]) + 5.0 * Bz_sq * p_tensor_rhs[5] + By_sq * (2.0 * p_tensor_rhs[0] + 3.0 * p_tensor_rhs[5]) + Bx_sq * 
    (2.0 * p_tensor_rhs[3] + 3.0 * p_tensor_rhs[5])) * q_over_m_sq - 2.0 * dt1_cu * (Bx_sq * (3.0 * Bz * p_tensor_rhs[1] - By * p_tensor_rhs[2]) - By * (3.0 * By * Bz *
    p_tensor_rhs[1] + By_sq * p_tensor_rhs[2] + 4.0 * Bz_sq * p_tensor_rhs[2]) + Bx_cu * p_tensor_rhs[4] + Bx * (3.0 * By * Bz * (-p_tensor_rhs[0] + p_tensor_rhs[3]) +
    By_sq * p_tensor_rhs[4] + 4.0 * Bz_sq * p_tensor_rhs[4])) * q_over_m_cu + 2.0 * dt1_qu * (-2.0 * Bx_cu * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) - 2.0 * Bx *
    (By_sq - 2.0 * Bz_sq) * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + By_sq * Bz_sq * (p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3]) - 2.0 * By_cu * Bz *
    p_tensor_rhs[4] + 4.0 * By * Bz_cu * p_tensor_rhs[4] + 2.0 * Bz_qu * p_tensor_rhs[5] + By_qu * (p_tensor_rhs[0] + p_tensor_rhs[5]) + Bx_qu * (p_tensor_rhs[3] +
    p_tensor_rhs[5]) + Bx_sq * (Bz_sq * (3.0 * p_tensor_rhs[0] + p_tensor_rhs[3]) - 2.0 * By * Bz * p_tensor_rhs[4] + By_sq * (p_tensor_rhs[0] + p_tensor_rhs[3] + 2.0 *
    p_tensor_rhs[5]))) * q_over_m_qu) / denom - p_tensor_old[5];
}

void
implicit_em_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, 
  double fluid_rhs_s[GKYL_MAX_SPECIES][4], double* fluid_s[GKYL_MAX_SPECIES],
  const double *app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* ext_em)
{
  int nfluids = mom_em->nfluids;
  double epsilon0 = mom_em->epsilon0;

  double Bx = em[3] + ext_em[3];
  double By = em[4] + ext_em[4];
  double Bz = em[5] + ext_em[5];
  double B_mag = sqrt((Bx * Bx) + (By * By) + (Bz * Bz));
  
  double bx = 0.0, by = 0.0, bz = 0.0;
  if (B_mag > 0.0) {
    bx = Bx / B_mag;
    by = By / B_mag;
    bz = Bz / B_mag;
  }

  double q_over_m[GKYL_MAX_SPECIES];
  double wc_dt[GKYL_MAX_SPECIES];
  double wp_dt_sq[GKYL_MAX_SPECIES];
  double J_old[GKYL_MAX_SPECIES][3];
  double J[GKYL_MAX_SPECIES][3];

  double scale_fact_E = 1.0;
  double scale_fact_curr = 1.0;
  double t_ramp_E = mom_em->t_ramp_E;
  double t_ramp_curr = mom_em->t_ramp_curr;

  double w0_sq = 0.0;
  double gam_sq = 0.0;
  double delta = 0.0;
  double Kx = 0.0, Ky = 0.0, Kz = 0.0;

  for (int i = 0; i < nfluids; i++) {
    double q = mom_em->param[i].charge;
    double m = mom_em->param[i].mass;
    q_over_m[i] = q / m;

    const double *app_accel = app_accel_s[i];

    double rho = fluid_rhs_s[i][0];
    double mom_x = fluid_rhs_s[i][1], mom_y = fluid_rhs_s[i][2], mom_z = fluid_rhs_s[i][3];

    J_old[i][0] = mom_x * q_over_m[i];
    J_old[i][1] = mom_y * q_over_m[i];
    J_old[i][2] = mom_z * q_over_m[i];

    if (mom_em->ramp_app_E) {
      scale_fact_E = fmin(1.0, t_curr / t_ramp_E);
    }

    J[i][0] = J_old[i][0] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_em[0] * scale_fact_E) + app_accel[0]));
    J[i][1] = J_old[i][1] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_em[1] * scale_fact_E) + app_accel[1]));
    J[i][2] = J_old[i][2] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_em[2] * scale_fact_E) + app_accel[2]));

    wc_dt[i] = q_over_m[i] * B_mag * dt;
    wp_dt_sq[i] = (rho * (q_over_m[i] * q_over_m[i]) * (dt * dt)) / epsilon0;

    double denom = 1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0);
    w0_sq += wp_dt_sq[i] / denom;
    gam_sq += (wp_dt_sq[i] * (wc_dt[i] * wc_dt[i])) / denom;
    delta += (wp_dt_sq[i] * wc_dt[i]) / denom;

    Kx -= (dt / denom) * (J[i][0] + (((wc_dt[i] * wc_dt[i]) / 4.0) * bx * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((by * J[i][2]) - (bz * J[i][1]))));
    Ky -= (dt / denom) * (J[i][1] + (((wc_dt[i] * wc_dt[i]) / 4.0) * by * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((bz * J[i][0]) - (bx * J[i][2]))));
    Kz -= (dt / denom) * (J[i][2] + (((wc_dt[i] * wc_dt[i]) / 4.0) * bz * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((bx * J[i][1]) - (by * J[i][0]))));
  }

  double Delta_sq = (delta * delta) / (1.0 + (w0_sq / 4.0));

  double Fx_old = em[0] * epsilon0;
  double Fy_old = em[1] * epsilon0;
  double Fz_old = em[2] * epsilon0;

  if (mom_em->ramp_app_curr) {
    scale_fact_curr = fmin(1.0, t_curr / t_ramp_curr);
  }

  double Fx = Fx_old - (0.5 * dt * app_current[0] * scale_fact_curr);
  double Fy = Fy_old - (0.5 * dt * app_current[1] * scale_fact_curr);
  double Fz = Fz_old - (0.5 * dt * app_current[2] * scale_fact_curr);

  double Fx_K = Fx + (0.5 * Kx);
  double Fy_K = Fy + (0.5 * Ky);
  double Fz_K = Fz + (0.5 * Kz);

  double Fx_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fx_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * bx * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((by * Fz_K) - (bz * Fy_K))));
  double Fy_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fy_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * by * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((bz * Fx_K) - (bx * Fz_K))));
  double Fz_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fz_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * bz * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((bx * Fy_K) - (by * Fx_K))));

  if (mom_em->static_field) {
    em[0] = Fx_old / epsilon0;
    em[1] = Fy_old / epsilon0;
    em[2] = Fz_old / epsilon0;
  }  
  else {
    em[0] = ((2.0 * Fx_bar) - Fx_old) / epsilon0;
    em[1] = ((2.0 * Fy_bar) - Fy_old) / epsilon0;
    em[2] = ((2.0 * Fz_bar) - Fz_old) / epsilon0;
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double Jx_star = J[i][0] + (Fx_bar * ((wp_dt_sq[i] / dt) / 2.0));
    double Jy_star = J[i][1] + (Fy_bar * ((wp_dt_sq[i] / dt) / 2.0));
    double Jz_star = J[i][2] + (Fz_bar * ((wp_dt_sq[i] / dt) / 2.0));

    double Jx_new = ((2.0 * (Jx_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * bx * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((by * Jz_star) - (bz * Jy_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][0];
    double Jy_new = ((2.0 * (Jy_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * by * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((bz * Jx_star) - (bx * Jz_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][1];
    double Jz_new = ((2.0 * (Jz_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * bz * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((bx * Jy_star) - (by * Jx_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][2];
    
    f[1] = Jx_new / q_over_m[i];
    f[2] = Jy_new / q_over_m[i];
    f[3] = Jz_new / q_over_m[i];
  }
}

void
implicit_neut_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt,
  double fluid_rhs_s[GKYL_MAX_SPECIES][4], double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], const double* p_rhs_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double *app_accel = app_accel_s[i];

    double rho = fluid_rhs_s[i][0];

    f[1] += dt * rho * app_accel[0];
    f[2] += dt * rho * app_accel[1];
    f[3] += dt * rho * app_accel[2];

    // Pressure-tensor source for neutral ten-moment fluids. The charged path
    // folds both the heat-flux closure increment and the k0 integrating-factor
    // relaxation into pressure_tensor_rotate, which vanishes at q/m = 0; here we
    // apply them explicitly: a forward-Euler increment by the closure rhs,
    // followed by relaxation toward the scalar pressure at rate nu = v_th * k0.
    // With no closure p_rhs is zero and (typically) k0 = 0, so this is a no-op.
    if (mom_em->param[i].type == GKYL_EQN_TEN_MOMENT) {
      const double *p_rhs = p_rhs_s[i];
      double k0 = mom_em->param[i].k0;

      double rho_f = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];

      // Thermal pressure tensor (subtract the bulk-flow contribution).
      double p_tensor[6];
      p_tensor[0] = f[4] - ((mom_x * mom_x) / rho_f);
      p_tensor[1] = f[5] - ((mom_x * mom_y) / rho_f);
      p_tensor[2] = f[6] - ((mom_x * mom_z) / rho_f);
      p_tensor[3] = f[7] - ((mom_y * mom_y) / rho_f);
      p_tensor[4] = f[8] - ((mom_y * mom_z) / rho_f);
      p_tensor[5] = f[9] - ((mom_z * mom_z) / rho_f);

      double p_scalar = (1.0 / 3.0) * (p_tensor[0] + p_tensor[3] + p_tensor[5]);
      double v_th = sqrt(p_scalar / rho_f);
      double exp_nu = exp(v_th * k0 * dt);

      // Forward-Euler increment from the closure heat-flux divergence.
      for (int n = 0; n < 6; n++) {
        p_tensor[n] += dt * p_rhs[4 + n];
      }

      // Integrating-factor relaxation toward isotropy (toward p_scalar).
      p_tensor[0] = ((p_tensor[0] - p_scalar) / exp_nu) + p_scalar;
      p_tensor[1] = p_tensor[1] / exp_nu;
      p_tensor[2] = p_tensor[2] / exp_nu;
      p_tensor[3] = ((p_tensor[3] - p_scalar) / exp_nu) + p_scalar;
      p_tensor[4] = p_tensor[4] / exp_nu;
      p_tensor[5] = ((p_tensor[5] - p_scalar) / exp_nu) + p_scalar;

      // Write back the conserved pressure-tensor components.
      f[4] = p_tensor[0] + ((mom_x * mom_x) / rho_f);
      f[5] = p_tensor[1] + ((mom_x * mom_y) / rho_f);
      f[6] = p_tensor[2] + ((mom_x * mom_z) / rho_f);
      f[7] = p_tensor[3] + ((mom_y * mom_y) / rho_f);
      f[8] = p_tensor[4] + ((mom_y * mom_z) / rho_f);
      f[9] = p_tensor[5] + ((mom_z * mom_z) / rho_f);
    }
  }
}

void
implicit_collision_source_update(const gkyl_moment_em_coupling* mom_em, double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
  for (int i  = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      nu_base[i][j] = (mom_em->nu_base)[i][j];
    }
  }

  double nu[nfluids * nfluids];
  for (int i = 0; i < nfluids; i++) {
    double *nu_i = nu + (nfluids * i);

    for (int j = 0; j < nfluids; j++) {
      double rho = fluid_s[j][0];

      nu_i[j] = nu_base[i][j] * rho;
    }
  }

  double lhs[nfluids][nfluids];
  double rhs[nfluids][3];
  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      lhs[i][j] = 0.0;
    }

    rhs[i][0] = 0.0; rhs[i][1] = 0.0; rhs[i][2] = 0.0;
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];

    rhs[i][0] = mom_x / rho;
    rhs[i][1] = mom_y / rho;
    rhs[i][2] = mom_z / rho;

    lhs[i][i] = 1.0;
    double* nu_i = nu + (nfluids * i);

    for (int j = 0; j < nfluids; j++) {
      if (i == j) {
        lhs[i][i] += 0.5 * dt * nu_i[i];
      }
      else {
        double dt_nu_ij = 0.5 * dt * nu_i[j];
        lhs[i][i] += dt_nu_ij;
        lhs[i][j] -= dt_nu_ij;
      }
    }
  }

  struct gkyl_mat *lhs_mat = gkyl_mat_new(nfluids, nfluids, 0.0);
  struct gkyl_mat *rhs_mat = gkyl_mat_new(nfluids, 3, 0.0);

  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      gkyl_mat_set(lhs_mat, i, j, lhs[i][j]);
    }

    gkyl_mat_set(rhs_mat, i, 0, rhs[i][0]);
    gkyl_mat_set(rhs_mat, i, 1, rhs[i][1]);
    gkyl_mat_set(rhs_mat, i, 2, rhs[i][2]);
  }

  gkyl_mem_buff sol_buff = gkyl_mem_buff_new(sizeof(long[nfluids]));
  bool status = gkyl_mat_linsolve_lu(lhs_mat, rhs_mat, gkyl_mem_buff_data(sol_buff));

  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < 3; j++) {
      rhs[i][j] = gkyl_mat_get(rhs_mat, i, j);
    }
  }

  double rhs_T[nfluids][1];
  for (int i = 0; i < nfluids; i++) {
    rhs_T[i][0] = 0.0;

    for (int j = 0; j < nfluids; j++) {
      lhs[i][j] = 0.0;
    }
  }

  double T[nfluids];
  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    double m = mom_em->param[i].mass;

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];
    double E = f[4];
    double internal_energy = E - (0.5 * ((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho);

    T[i] = (internal_energy / rho) * m;
    rhs_T[i][0] = T[i];
    lhs[i][i] = 1.0;

    double *nu_i = nu + (nfluids * i);

    for (int j = 0; j < nfluids; j++) {
      if (i != j) {
        double m_j = mom_em->param[j].mass;
        
        double du_sq = ((rhs[i][0] - rhs[j][0]) * (rhs[i][0] - rhs[j][0])) + ((rhs[i][1] - rhs[j][1]) * (rhs[i][1] - rhs[j][1])) +
          ((rhs[i][2] - rhs[j][2]) * (rhs[i][2] - rhs[j][2]));
        double coeff_ij = (dt * nu_i[j] * m) / (m + m_j);

        rhs_T[i][0] += 0.5 * coeff_ij * m_j * du_sq;
        lhs[i][i] += coeff_ij;
        lhs[i][j] -= coeff_ij;
      }
    }
  }

  struct gkyl_mat *lhs_T_mat = gkyl_mat_new(nfluids, nfluids, 0.0);
  struct gkyl_mat *rhs_T_mat = gkyl_mat_new(nfluids, 1, 0.0);

  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      gkyl_mat_set(lhs_T_mat, i, j, lhs[i][j]);
    }

    gkyl_mat_set(rhs_T_mat, i, 0, rhs_T[i][0]);
  }

  status = gkyl_mat_linsolve_lu(lhs_T_mat, rhs_T_mat, gkyl_mem_buff_data(sol_buff));

  for (int i = 0; i < nfluids; i++) {
    rhs_T[i][0] = gkyl_mat_get(rhs_T_mat, i, 0);
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    double m = mom_em->param[i].mass;

    double rho = f[0];

    f[4] = ((2.0 * rhs_T[i][0] - T[i]) * rho) / m;
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];
    double E = f[4];

    f[1] = 2.0 * rho * rhs[i][0] - mom_x;
    f[2] = 2.0 * rho * rhs[i][1] - mom_y;
    f[3] = 2.0 * rho * rhs[i][2] - mom_z;

    mom_x = f[1], mom_y = f[2], mom_z = f[3];
    f[4] = E + (0.5 * ((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho);
  }

  gkyl_mem_buff_release(sol_buff);
  gkyl_mat_release(lhs_mat);
  gkyl_mat_release(rhs_mat);
  gkyl_mat_release(lhs_T_mat);
  gkyl_mat_release(rhs_T_mat);
}

void
implicit_frictional_source_update_half(const gkyl_moment_em_coupling* mom_em, const double Z, const double T_elc, const double Lambda_ee,
  double t_curr, const double dt, double* f_elc_old, double* f_ion_old, double* f_elc_new, double* f_ion_new,
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em_old, double* em_new, const double* app_current, const double* ext_em)
{
  int nfluids = mom_em->nfluids;
  double pi = M_PI;

  if (nfluids == 2) {
    double mass_elc = mom_em->param[0].mass;
    double mass_ion = mom_em->param[1].mass;
    double charge_elc = mom_em->param[0].charge;
    double charge_ion = mom_em->param[1].charge;
    double epsilon0 = mom_em->epsilon0;

    const double *app_accel_elc = app_accel_s[0];
    const double *app_accel_ion = app_accel_s[1];

    double Ex = em_old[0] + ext_em[0];
    double Ey = em_old[1] + ext_em[1];
    double Ez = em_old[2] + ext_em[2];

    double rho_elc = f_elc_old[0];
    double rho_ion = f_ion_old[0];

    double u_elc = f_elc_old[1], v_elc = f_elc_old[2], w_elc = f_elc_old[3];
    double u_ion = f_ion_old[1], v_ion = f_ion_old[2], w_ion = f_ion_old[3];

    double E_elc = f_elc_old[4];
    double E_ion = f_ion_old[4];

    double n_elc = rho_elc / mass_elc;
    double n_ion = rho_ion / mass_ion;

    double sigma_elc = -charge_elc * n_elc;
    double sigma_ion = charge_ion * n_ion;
    double s_elc = sigma_elc / rho_elc;
    double s_ion = sigma_ion / rho_ion;

    double tau_ei = (1.0 / Z) * ((3.0 * sqrt(mass_elc) * ((4.0 * pi * epsilon0) * (4.0 * pi * epsilon0)) * pow(T_elc, 3.0 / 2.0)) /
      (4.0 * sqrt(2.0 * pi) * n_elc * exp(4.0) * log(Lambda_ee)));
    double alpha_par = 1.0 - (pow(Z, 2.0 / 3.0) / ((1.46 * pow(Z, 2.0 / 3.0)) - (0.33 * pow (Z, 1.0 / 3.0)) + 0.888));

    double A_ee = 1.0 + ((0.5 * dt) * (1.0 / tau_ei) * alpha_par);
    double A_ei = -(0.5 * dt) * (1.0 / tau_ei) * alpha_par;
    double C_e = -(0.5 * dt) * s_elc;

    double A_ie = -(0.5 * dt) * (1.0 / tau_ei) * (rho_elc / rho_ion) * alpha_par;
    double A_ii = 1.0 + ((0.5 * dt) * (rho_elc / rho_ion) * (1.0 / tau_ei) * alpha_par);
    double C_i = -(0.5 * dt) * s_ion;

    double D_e = (dt / (2.0 * epsilon0)) * sigma_elc;
    double D_i = (dt / (2.0 * epsilon0)) * sigma_ion;

    double det = (A_ei * (-A_ie + (C_i * D_e))) + (C_e * ((-A_ii * D_e) + (A_ie * D_i))) + (A_ee * (A_ii - (C_i * D_i)));
    double mat_11 = (A_ii - (C_i * D_i)) / det;
    double mat_12 = (A_ei - (C_e * D_i)) / det;
    double mat_13 = (-(A_ii * C_e) + (A_ei * C_i)) / det;
    double mat_21 = (A_ie - (C_i * D_e)) / det;
    double mat_22 = (A_ee - (C_e * D_e)) / det;
    double mat_23 = ((A_ie * C_e) - (A_ee * C_i)) / det;
    double mat_31 = (-(A_ii * D_e) + (A_ie * D_i)) / det;
    double mat_32 = ((A_ei * D_e) - (A_ee * D_i)) / det;
    double mat_33 = (-(A_ei * A_ie) + (A_ee * A_ii)) / det;

    f_elc_new[1] = (mat_11 * (u_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[0]))) +
      (mat_12 * (u_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[0]))) +
      (mat_13 * (Ex - ((dt / (2.0 * epsilon0)) * app_current[0])));
    f_elc_new[2] = (mat_11 * (v_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[1]))) +
      (mat_12 * (v_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[1]))) +
      (mat_13 * (Ey - ((dt / (2.0 * epsilon0)) * app_current[1])));
    f_elc_new[3] = (mat_11 * (w_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[2]))) +
      (mat_12 * (w_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[2]))) +
      (mat_13 * (Ez - ((dt / (2.0 * epsilon0)) * app_current[2])));

    f_ion_new[1] = (mat_21 * (u_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[0]))) +
      (mat_22 * (u_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[0]))) +
      (mat_23 * (Ex - ((dt / (2.0 * epsilon0)) * app_current[0])));
    f_ion_new[2] = (mat_21 * (v_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[1]))) +
      (mat_22 * (v_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[1]))) +
      (mat_23 * (Ey - ((dt / (2.0 * epsilon0)) * app_current[1])));
    f_ion_new[3] = (mat_21 * (w_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[2]))) +
      (mat_22 * (w_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[2]))) +
      (mat_23 * (Ez - ((dt / (2.0 * epsilon0)) * app_current[2])));
    
    em_new[0] = (mat_31 * (u_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[0]))) +
      (mat_32 * (u_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[0]))) +
      (mat_33 * (Ex - ((dt / (2.0 * epsilon0)) * app_current[0])));
    em_new[1] = (mat_31 * (v_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[1]))) +
      (mat_32 * (v_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[1]))) +
      (mat_33 * (Ey - ((dt / (2.0 * epsilon0)) * app_current[1])));
    em_new[2] = (mat_31 * (w_elc + (0.5 * dt * (1.0 / rho_elc) * app_accel_elc[2]))) +
      (mat_32 * (w_ion + (0.5 * dt * (1.0 / rho_ion) * app_accel_ion[2]))) +
      (mat_33 * (Ez - ((dt / (2.0 * epsilon0)) * app_current[2])));
    
    em_new[0] -= ext_em[0];
    em_new[1] -= ext_em[1];
    em_new[2] -= ext_em[2];

    for (int i = 3; i < 8; i++) {
      em_new[i] = em_old[i];
    }

    f_elc_new[4] = f_elc_old[4];
    f_ion_new[4] = f_ion_old[4];

    f_elc_new[0] = f_elc_old[0];
    f_ion_new[0] = f_ion_old[0];
  }
}

void
implicit_frictional_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* ext_em)
{
  int nfluids = mom_em->nfluids;

  if (nfluids == 2) {
    double *f_elc = fluid_s[0];
    double *f_ion = fluid_s[1];

    double Z = mom_em->friction_Z;
    double T_elc = mom_em->friction_T_elc;
    double Lambda_ee = mom_em->friction_Lambda_ee;

    double f_elc_stage1[5], f_elc_new[5];
    double f_ion_stage1[5], f_ion_new[5];
    double f_elc_old[5], f_ion_old[5];
    double em_stage1[8], em_new[8], em_old[8];

    for (int i = 0; i < 5; i++) {
      f_elc_old[i] = f_elc[i];
      f_ion_old[i] = f_ion[i];
    }
    
    for (int i = 0; i < 8; i++) {
      em_old[i] = em[i];
    }

    implicit_frictional_source_update_half(mom_em, Z, T_elc, Lambda_ee, t_curr, dt, f_elc_old, f_ion_old, f_elc_stage1, f_ion_stage1,
      app_accel_s, em_old, em_stage1, app_current, ext_em);
    implicit_frictional_source_update_half(mom_em, Z, T_elc, Lambda_ee, t_curr + (0.5 * dt), dt, f_elc_stage1, f_ion_stage1, f_elc_new, f_ion_new,
      app_accel_s, em_stage1, em_new, app_current, ext_em);

    for (int i = 0; i < 5; i++) {
      f_elc[i] = f_elc_new[i];
      f_ion[i] = f_ion_new[i];
    }

    for (int i = 0; i < 8; i++) {
      em[i] = em_new[i];
    }
  }
}

void
implicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], const double* p_rhs_s[GKYL_MAX_SPECIES], double* em, const double* app_current,
  const double* ext_em, const double* nT_sources_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double ke_old[GKYL_MAX_SPECIES];
  double energy_old[GKYL_MAX_SPECIES];
  double fluid_rhs[GKYL_MAX_SPECIES][4];
  double p_tensor_old[6], p_tensor_rhs[6], p_tensor_new[GKYL_MAX_SPECIES][6];

  for (int i = 0; i < nfluids; i++) {
    const double *f = fluid_s[i];
    const double *p_rhs = p_rhs_s[i];

    double q = mom_em->param[i].charge;
    double m = mom_em->param[i].mass;
    double k0 = mom_em->param[i].k0;

    // Setup RHS of implicit solve with fluid variables at known time-step
    // includes potential contributions from transport terms/density & momentum sources
    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];
    double rho_rhs = p_rhs[0];
    double mom_x_rhs = p_rhs[1], mom_y_rhs = p_rhs[2], mom_z_rhs = p_rhs[3];

    fluid_rhs[i][0] = rho + (0.5 * dt * rho_rhs);
    fluid_rhs[i][1] = mom_x + (0.5 * dt * mom_x_rhs);
    fluid_rhs[i][2] = mom_y + (0.5 * dt * mom_y_rhs);
    fluid_rhs[i][3] = mom_z + (0.5 * dt * mom_z_rhs);

    if (mom_em->param[i].type == GKYL_EQN_EULER) {
      double energy = f[4];

      // Include potential contributions from transport terms to energy
      double energy_rhs = p_rhs[4];

      // kinetic energy at known time (including potential transport terms)
      ke_old[i] = 0.5 * (((fluid_rhs[i][1] * fluid_rhs[i][1]) 
        + (fluid_rhs[i][2] * fluid_rhs[i][2]) 
        + (fluid_rhs[i][3] * fluid_rhs[i][3])) / fluid_rhs[i][0]);

      // total energy at known time (including potential transport terms)
      energy_old[i] = energy + (0.5 * dt * energy_rhs);
    }
    else if (mom_em->param[i].type == GKYL_EQN_TEN_MOMENT) {
      double q_over_m = q / m;

      double p11 = f[4], p12 = f[5], p13 = f[6];
      double p22 = f[7], p23 = f[8], p33 = f[9];

      p_tensor_old[0] = p11 - ((mom_x * mom_x) / rho);
      p_tensor_old[1] = p12 - ((mom_x * mom_y) / rho);
      p_tensor_old[2] = p13 - ((mom_x * mom_z) / rho);
      p_tensor_old[3] = p22 - ((mom_y * mom_y) / rho);
      p_tensor_old[4] = p23 - ((mom_y * mom_z) / rho);
      p_tensor_old[5] = p33 - ((mom_z * mom_z) / rho);

      double p11_rhs = p_rhs[4], p12_rhs = p_rhs[5], p13_rhs = p_rhs[6];
      double p22_rhs = p_rhs[7], p23_rhs = p_rhs[8], p33_rhs = p_rhs[9];

      p_tensor_rhs[0] = p_tensor_old[0] + (0.5 * dt * p11_rhs);
      p_tensor_rhs[1] = p_tensor_old[1] + (0.5 * dt * p12_rhs);
      p_tensor_rhs[2] = p_tensor_old[2] + (0.5 * dt * p13_rhs);
      p_tensor_rhs[3] = p_tensor_old[3] + (0.5 * dt * p22_rhs);
      p_tensor_rhs[4] = p_tensor_old[4] + (0.5 * dt * p23_rhs);
      p_tensor_rhs[5] = p_tensor_old[5] + (0.5 * dt * p33_rhs);

      double p = (1.0 / 3.0) * (p_tensor_old[0] + p_tensor_old[3] + p_tensor_old[5]);
      double v_th = sqrt(p / rho);
      double nu = v_th * k0;
      double exp_nu = exp(nu * dt);

      if (mom_em->is_charged_species) {
        pressure_tensor_rotate(q_over_m, dt, em, ext_em, p_tensor_old, p_tensor_rhs, p_tensor_new[i]);
      }

      p_tensor_new[i][0] = ((p_tensor_new[i][0] - p) / exp_nu) + p;
      p_tensor_new[i][1] = p_tensor_new[i][1] / exp_nu;
      p_tensor_new[i][2] = p_tensor_new[i][2] / exp_nu;
      p_tensor_new[i][3] = ((p_tensor_new[i][3] - p) / exp_nu) + p;
      p_tensor_new[i][4] = p_tensor_new[i][4] / exp_nu;
      p_tensor_new[i][5] = ((p_tensor_new[i][5] - p) / exp_nu) + p;
    }
  }

  if (mom_em->is_charged_species) {
    implicit_em_source_update(mom_em, t_curr, dt, fluid_rhs, fluid_s, app_accel_s, em, app_current, ext_em);
  }
  else {
    implicit_neut_source_update(mom_em, t_curr, dt, fluid_rhs, fluid_s, app_accel_s, p_rhs_s);
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    if (mom_em->param[i].type == GKYL_EQN_EULER) {
      double rho = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];
      double energy = f[4];
      
      // Energy at new time is new kinetic energy plus potential contribution from
      // transport terms. We use a time-centered approach even though the update
      // from the transport terms is a simple forward Euler.
      f[4] = (0.5 * ((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho) 
              + 2.0*energy_old[i] - energy - ke_old[i];
    }
    // As I do not understand how the source terms for gradient-based closure interact with the source terms for the expanding-box
    // model, I'm disabling the former whenever the latter are present, pro tem. This should be updated! -JG 07/25/24
    else if (mom_em->param[i].type == GKYL_EQN_TEN_MOMENT && mom_em->param[i].charge != 0.0) {
      double rho = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];

      f[4] = ((mom_x * mom_x) / rho) + p_tensor_new[i][0];
      f[5] = ((mom_x * mom_y) / rho) + p_tensor_new[i][1];
      f[6] = ((mom_x * mom_z) / rho) + p_tensor_new[i][2];
      f[7] = ((mom_y * mom_y) / rho) + p_tensor_new[i][3];
      f[8] = ((mom_y * mom_z) / rho) + p_tensor_new[i][4];
      f[9] = ((mom_z * mom_z) / rho) + p_tensor_new[i][5];
    }
  }

  if (mom_em->has_collision) {
    implicit_collision_source_update(mom_em, dt, fluid_s);
  }

  // These terms are currently handled by their own specialized forcing solver(s). To be revisited...
  if (mom_em->has_nT_sources) {
    explicit_nT_source_update(mom_em, dt, fluid_s, nT_sources_s);
  }
  if (mom_em->has_frictional_sources) {
    if (mom_em->use_explicit_friction) {
      explicit_frictional_source_update(mom_em, t_curr, dt, fluid_s);
    }
    else {
      implicit_frictional_source_update(mom_em, t_curr, dt, fluid_s, app_accel_s, em, app_current, ext_em);
    }
  }
  if (mom_em->has_volume_sources) {
    explicit_volume_source_update(mom_em, t_curr, dt, fluid_s, em, ext_em);
  }
  if (mom_em->has_reactive_sources) {
    explicit_reactive_source_update(mom_em, t_curr, dt, fluid_s);
  }
  if (mom_em->has_einstein_medium_sources) {
    explicit_medium_source_update(mom_em, t_curr, dt, fluid_s);
  }
  if (mom_em->has_gr_ultra_rel_sources) {
    explicit_gr_ultra_rel_source_update(mom_em, t_curr, dt, fluid_s);
  }
  if (mom_em->has_gr_euler_sources) {
    explicit_gr_euler_source_update(mom_em, t_curr, dt, fluid_s);
  }
  if (mom_em->has_gr_twofluid_sources) {
    // Subcycling capabilities: ratio set to 1 by default.
    for (int i = 0; i < 1; i++) {
      explicit_gr_twofluid_source_update(mom_em, t_curr, 1.0 * dt, fluid_s);
    }
  }
  if (mom_em->has_vacuum_einstein_sources) {
    explicit_vacuum_einstein_source_update(mom_em, t_curr, dt, fluid_s);
  }
  if (mom_em->has_vacuum_einstein_conformal_sources) {
    explicit_vacuum_einstein_conformal_source_update(mom_em, t_curr, dt, fluid_s);
  }
  if (mom_em->has_gr_mhd_sources) {
    explicit_gr_mhd_source_update(mom_em, t_curr, dt, fluid_s);
  }
}
