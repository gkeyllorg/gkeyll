#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_mhd_tetrad.h>
#include <gkyl_wv_gr_mhd_tetrad_priv.h>

void
gkyl_gr_mhd_tetrad_flux(double gas_gamma, double light_speed, double b_fact, const double q[75], double flux[75])
{
  double v[75] = { 0.0 };
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];
  double psi = v[8];

  bool in_excision_region = false;
  if (v[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double W = 1.0 / (sqrt(1.0 - ((vx * vx) + (vy * vy) + (vz * vz))));
    if ((vx * vx) + (vy * vy) + (vz * vz) > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double b0 = W * ((mag_x * vx) + (mag_y * vy) + (mag_z * vz));
    double b1 = (mag_x / W) + (b0 * vx);
    double b2 = (mag_y / W) + (b0 * vy);
    double b3 = (mag_z / W) + (b0 * vz);
    double b_sq = ((mag_x * mag_x) + (mag_y * mag_y) + (mag_z * mag_z) + (b0 * b0)) / (W * W);

    double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
    double p_star = p + (0.5 * b_sq);

    double D = rho * W;
    double Sx = (rho * h_star * (W * W) * vx) - (b0 * b1);
    double Sy = (rho * h_star * (W * W) * vy) - (b0 * b2);
    double Sz = (rho * h_star * (W * W) * vz) - (b0 * b3);
    double Etot = (rho * h_star * (W * W)) - p_star - (b0 * b0) - (rho * W);

    flux[0] = D * vx;
    flux[1] = (Sx * vx) + p_star - ((b1 * mag_x) / W);
    flux[2] = (Sy * vx) - ((b2 * mag_x) / W);
    flux[3] = (Sz * vx) - ((b3 * mag_x) / W);
    flux[4] = (Etot * vx) + (p_star * vx) - ((b0 * mag_x) / W);

    flux[5] = (vx * mag_x) - (vx * mag_x) + (b_fact * psi);
    flux[6] = (vx * mag_y) - (vy * mag_x);
    flux[7] = (vx * mag_z) - (vz * mag_x);
    flux[8] = b_fact * (light_speed * light_speed) * mag_x;

    for (int i = 9; i < 75; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 75; i++) {
      flux[i] = 0.0;
    }
  }
}

void
gkyl_gr_mhd_tetrad_flux_correction(double gas_gamma, double light_speed, double b_fact, const double q[75], const double flux_sr[75], double flux_gr[75])
{
  double v[75] = { 0.0 };
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];
  double psi = v[8];

  double lapse = v[9];
  double shift_x = v[10];
  double shift_y = v[11];
  double shift_z = v[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[13]; spatial_metric[0][1] = v[14]; spatial_metric[0][2] = v[15];
  spatial_metric[1][0] = v[16]; spatial_metric[1][1] = v[17]; spatial_metric[1][2] = v[18];
  spatial_metric[2][0] = v[19]; spatial_metric[2][1] = v[20]; spatial_metric[2][2] = v[21];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  bool in_excision_region = false;
  if (v[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W_flat = 1.0 / (sqrt(1.0 - ((vx * vx) + (vy * vy) + (vz * vz))));
    if ((vx * vx) + (vy * vy) + (vz * vz) > 1.0 - pow(10.0, -8.0)) {
      W_flat = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double W_curved = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W_curved = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double b0_flat = W_flat * ((mag_x * vx) + (mag_y * vy) + (mag_z * vz));
    double b1_flat = (mag_x / W_flat) + (b0_flat * vx);
    double b2_flat = (mag_y / W_flat) + (b0_flat * vy);
    double b3_flat = (mag_z / W_flat) + (b0_flat * vz);
    double b_sq_flat = ((mag_x * mag_x) + (mag_y * mag_y) + (mag_z * mag_z) + (b0_flat * b0_flat)) / (W_flat * W_flat);

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
    for (int i = 0; i < 3; i++) {
      cov_vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }
    
    double b0_curved = 0.0;
    for (int i = 0; i < 3; i++) {
      b0_curved += W_curved * mag[i] * (cov_vel[i] / lapse);
    }

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W_curved / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W_curved * vel[i]) - (shift[i] * (W_curved / lapse));
    }

    double b_curved[3];
    for (int i = 0; i < 3; i++) {
      b_curved[i] = (mag[i] + (lapse * b0_curved * spacetime_vel[i + 1])) / W_curved;
    }

    double b_sq_curved = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq_curved += (mag[i] * cov_mag[i]) / (W_curved * W_curved);
    }
    b_sq_curved += ((lapse * lapse) * (b0_curved * b0_curved)) / (W_curved * W_curved);

    double cov_b[3];
    for (int i = 0; i < 3; i++) {
      cov_b[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_b[i] += spatial_metric[i][j] * b_curved[j];
      }
    }

    double h_star_flat = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_flat / rho);
    double p_star_flat = p + (0.5 * b_sq_flat);

    double h_star_curved = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_curved / rho);
    double p_star_curved = p + (0.5 * b_sq_curved);

    if (fabs(vx) < pow(10.0, -12.0)) {
      if (vx > 0.0) {
        vx = pow(10.0, -12.0);
      }
      else {
        vx = -pow(10.0, -12.0);
      }
    }
    if (fabs(vy) < pow(10.0, -12.0)) {
      if (vy > 0.0) {
        vy = pow(10.0, -12.0);
      }
      else {
        vy = -pow(10.0, -12.0);
      }
    }
    if (fabs(vz) < pow(10.0, -12.0)) {
      if (vz > 0.0) {
        vz = pow(10.0, -12.0);
      }
      else {
        vz = -pow(10.0, -12.0);
      }
    }

    flux_gr[0] = (lapse * sqrt(spatial_det)) * ((flux_sr[0] * (vx - (shift_x / lapse)) * W_curved) / (vx * W_flat));

    flux_gr[1] = (lapse * sqrt(spatial_det)) * (((((((flux_sr[1] - p_star_flat + ((b1_flat * mag_x) / W_flat)) / vx) + (b0_flat * b1_flat)) *
      ((h_star_curved * (W_curved * W_curved) * cov_vel[0]) / (h_star_flat * (W_flat * W_flat) * vx))) - (lapse * b0_curved * cov_b[0])) *
      (vx - (shift_x / lapse))) + p_star_curved - ((cov_b[0] * mag_x) / W_curved));
    flux_gr[2] = (lapse * sqrt(spatial_det)) * (((((((flux_sr[2] + ((b2_flat * mag_x) / W_flat)) / vx) + (b0_flat * b2_flat)) *
      ((h_star_curved * (W_curved * W_curved) * cov_vel[1]) / (h_star_flat * (W_flat * W_flat) * vy))) - (lapse * b0_curved * cov_b[1])) *
      (vx - (shift_x / lapse))) - ((cov_b[1] * mag_x) / W_curved));
    flux_gr[3] = (lapse * sqrt(spatial_det)) * (((((((flux_sr[3] + ((b3_flat * mag_x) / W_flat)) / vx) + (b0_flat * b3_flat)) *
      ((h_star_curved * (W_curved * W_curved) * cov_vel[2]) / (h_star_flat * (W_flat * W_flat) * vz))) - (lapse * b0_curved * cov_b[2])) *
      (vx - (shift_x / lapse))) - ((cov_b[2] * mag_x) / W_curved));
    
    flux_gr[4] = (lapse * sqrt(spatial_det)) * (((((((flux_sr[4] - (p_star_flat * vx) + ((b0_flat * mag_x) / W_flat)) / vx) + p_star_flat +
      (b0_flat * b0_flat) + (rho * W_flat)) * ((h_star_curved * (W_curved * W_curved)) / (h_star_flat * (W_flat * W_flat)))) - p_star_curved -
      ((lapse * lapse) * (b0_curved * b0_curved)) - (rho * W_curved)) * (vx - (shift_x / lapse))) + (p_star_curved * vx) -
      (lapse * ((b0_curved * mag_x) / W_curved)));

    flux_gr[5] = (lapse * sqrt(spatial_det)) * flux_sr[5];
    flux_gr[6] = (lapse * sqrt(spatial_det)) * (flux_sr[6] - ((shift_x / lapse) * mag_y) + ((shift_y / lapse) * mag_x));
    flux_gr[7] = (lapse * sqrt(spatial_det)) * (flux_sr[7] - ((shift_x / lapse) * mag_z) + ((shift_z / lapse) * mag_x));
    flux_gr[8] = (lapse * sqrt(spatial_det)) * flux_sr[8];

    for (int i = 9; i < 75; i++) {
      flux_gr[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 75; i++) {
      flux_gr[i] = 0.0;
    }
  }
}

void
gkyl_gr_mhd_tetrad_prim_vars(double gas_gamma, const double q[75], double v[75])
{
  double lapse = q[9];
  double shift_x = q[10];
  double shift_y = q[11];
  double shift_z = q[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = q[13]; spatial_metric[0][1] = q[14]; spatial_metric[0][2] = q[15];
  spatial_metric[1][0] = q[16]; spatial_metric[1][1] = q[17]; spatial_metric[1][2] = q[18];
  spatial_metric[2][0] = q[19]; spatial_metric[2][1] = q[20]; spatial_metric[2][2] = q[21];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_tetrad_inv_spatial_metric(q, &inv_spatial_metric);
  
  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[22]; extrinsic_curvature[0][1] = q[23]; extrinsic_curvature[0][2] = q[24];
  extrinsic_curvature[1][0] = q[25]; extrinsic_curvature[1][1] = q[26]; extrinsic_curvature[1][2] = q[27];
  extrinsic_curvature[2][0] = q[28]; extrinsic_curvature[2][1] = q[29]; extrinsic_curvature[2][2] = q[30];

  double lapse_der[3];
  lapse_der[0] = q[32];
  lapse_der[1] = q[33];
  lapse_der[2] = q[34];

  double shift_der[3][3];
  shift_der[0][0] = q[35]; shift_der[0][1] = q[36]; shift_der[0][2] = q[37];
  shift_der[1][0] = q[38]; shift_der[1][1] = q[39]; shift_der[1][2] = q[40];
  shift_der[2][0] = q[41]; shift_der[2][1] = q[42]; shift_der[2][2] = q[43];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[44]; spatial_metric_der[0][0][1] = q[45]; spatial_metric_der[0][0][2] = q[46];
  spatial_metric_der[0][1][0] = q[47]; spatial_metric_der[0][1][1] = q[48]; spatial_metric_der[0][1][2] = q[49];
  spatial_metric_der[0][2][0] = q[50]; spatial_metric_der[0][2][1] = q[51]; spatial_metric_der[0][2][2] = q[52];

  spatial_metric_der[1][0][0] = q[53]; spatial_metric_der[1][0][1] = q[54]; spatial_metric_der[1][0][2] = q[55];
  spatial_metric_der[1][1][0] = q[56]; spatial_metric_der[1][1][1] = q[57]; spatial_metric_der[1][1][2] = q[58];
  spatial_metric_der[1][2][0] = q[59]; spatial_metric_der[1][2][1] = q[60]; spatial_metric_der[1][2][2] = q[61];

  spatial_metric_der[0][0][0] = q[62]; spatial_metric_der[0][0][1] = q[63]; spatial_metric_der[0][0][2] = q[64];
  spatial_metric_der[0][1][0] = q[65]; spatial_metric_der[0][1][1] = q[66]; spatial_metric_der[0][1][2] = q[67];
  spatial_metric_der[0][2][0] = q[68]; spatial_metric_der[0][2][1] = q[69]; spatial_metric_der[0][2][2] = q[70];

  double evol_param = q[71];
  double x = q[72];
  double y = q[73];
  double z = q[74];

  bool in_excision_region = false;
  if (q[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det);
    double momy = q[2] / sqrt(spatial_det);
    double momz = q[3] / sqrt(spatial_det);
    double Etot = q[4] / sqrt(spatial_det);

    double mag_x = q[5] / sqrt(spatial_det);
    double mag_y = q[6] / sqrt(spatial_det);
    double mag_z = q[7] / sqrt(spatial_det);
    double psi = q[8] / sqrt(spatial_det);

    double cov_mom[3];
    cov_mom[0] = momx; cov_mom[1] = momy; cov_mom[2] = momz;

    double mom[3];
    for (int i = 0; i < 3; i++) {
      mom[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        mom[i] += inv_spatial_metric[i][j] * cov_mom[j];
      }
    }

    double M_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      M_sq += (mom[i] * cov_mom[i]);
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }
    
    double mag_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      mag_sq += mag[i] * cov_mag[i];
    }

    double tau_star = 0.0;
    for (int i = 0; i < 3; i++) {
      tau_star += (mag[i] * cov_mom[i]);
    }

    double p_guess = 0.0, p_new = 0.0, rho = 0.0, z = 0.0;
    int iter = 0;

    while (iter < 100) {
      double a = Etot + D + p_guess + (0.5 * mag_sq);
      double d = 0.5 * ((M_sq * mag_sq) - (tau_star * tau_star));

      double phi = acos((1.0 / a) * sqrt((27.0 * d) / (4.0 * a)));
      double epsilon1 = ((1.0 / 3.0) * a) - ((2.0 / 3.0) * a * cos(((2.0 / 3.0) * phi) + ((2.0 / 3.0) * M_PI)));
      z = epsilon1 - mag_sq;

      double v_sq = ((M_sq * (z * z)) + ((tau_star * tau_star) * (mag_sq + (2.0 * z)))) / ((z * z) * (mag_sq + z) * (mag_sq + z));
      double W = 1.0 / sqrt(1.0 - v_sq);
      rho = D / W;
      double h = z / (W * W * rho);

      p_new = rho * (h - 1.0) * ((gas_gamma - 1.0) / gas_gamma);

      if (fabs(p_guess - p_new) < pow(10.0, -15.0)) {
        iter = 100;
      }
      else {
        iter += 1;
        p_guess = p_new;
      }
    }

    if (p_new < pow(10.0, -8.0)) {
      p_new = pow(10.0, -8.0);
    }
    if (rho < pow(10.0, -8.0)) {
      rho = pow(10.0, -8.0);
    }

    double vel[3];
    for (int i = 0; i < 3; i++) {
      vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        vel[i] += (inv_spatial_metric[i][j] * cov_mom[j]) / (z + mag_sq);
        vel[i] += ((mag[j] * cov_mom[j]) * mag[i]) / (z * (z + mag_sq));
      }
    }

    v[0] = rho;
    v[1] = vel[0];
    v[2] = vel[1];
    v[3] = vel[2];
    v[4] = p_new;

    v[5] = mag_x;
    v[6] = mag_y;
    v[7] = mag_z;
    v[8] = psi;

    v[9] = lapse;
    v[10] = shift_x;
    v[11] = shift_y;
    v[12] = shift_z;

    v[13] = spatial_metric[0][0]; v[14] = spatial_metric[0][1]; v[15] = spatial_metric[0][2];
    v[16] = spatial_metric[1][0]; v[17] = spatial_metric[1][1]; v[18] = spatial_metric[1][2];
    v[19] = spatial_metric[2][0]; v[20] = spatial_metric[2][1]; v[21] = spatial_metric[2][2];

    v[22] = extrinsic_curvature[0][0]; v[23] = extrinsic_curvature[0][1]; v[24] = extrinsic_curvature[0][2];
    v[25] = extrinsic_curvature[1][0]; v[26] = extrinsic_curvature[1][1]; v[27] = extrinsic_curvature[1][2];
    v[28] = extrinsic_curvature[2][0]; v[29] = extrinsic_curvature[2][1]; v[30] = extrinsic_curvature[2][2];

    v[31] = 1.0;

    v[32] = lapse_der[0];
    v[33] = lapse_der[1];
    v[34] = lapse_der[2];

    v[35] = shift_der[0][0]; v[36] = shift_der[0][1]; v[37] = shift_der[0][2];
    v[38] = shift_der[1][0]; v[39] = shift_der[1][1]; v[40] = shift_der[1][2];
    v[41] = shift_der[2][0]; v[42] = shift_der[2][1]; v[43] = shift_der[2][2];

    v[44] = spatial_metric_der[0][0][0]; v[45] = spatial_metric_der[0][0][1]; v[46] = spatial_metric_der[0][0][2];
    v[47] = spatial_metric_der[0][1][0]; v[48] = spatial_metric_der[0][1][1]; v[49] = spatial_metric_der[0][1][2];
    v[50] = spatial_metric_der[0][2][0]; v[51] = spatial_metric_der[0][2][1]; v[52] = spatial_metric_der[0][2][2];

    v[53] = spatial_metric_der[1][0][0]; v[54] = spatial_metric_der[1][0][1]; v[55] = spatial_metric_der[1][0][2];
    v[56] = spatial_metric_der[1][1][0]; v[57] = spatial_metric_der[1][1][1]; v[58] = spatial_metric_der[1][1][2];
    v[59] = spatial_metric_der[1][2][0]; v[60] = spatial_metric_der[1][2][1]; v[61] = spatial_metric_der[1][2][2];

    v[62] = spatial_metric_der[2][0][0]; v[63] = spatial_metric_der[2][0][1]; v[64] = spatial_metric_der[2][0][2];
    v[65] = spatial_metric_der[2][1][0]; v[66] = spatial_metric_der[2][1][1]; v[67] = spatial_metric_der[2][1][2];
    v[68] = spatial_metric_der[2][2][0]; v[69] = spatial_metric_der[2][2][1]; v[70] = spatial_metric_der[2][2][2];

    v[71] = evol_param;
    v[72] = x;
    v[73] = y;
    v[74] = z;
  }
  else {
    for (int i = 0; i < 75; i++) {
      v[i] = 0.0;
    }
    
    v[31] = -1.0;
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

void 
gkyl_gr_mhd_tetrad_inv_spatial_metric(const double q[75], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[13]; spatial_metric[0][1] = q[14]; spatial_metric[0][2] = q[15];
  spatial_metric[1][0] = q[16]; spatial_metric[1][1] = q[17]; spatial_metric[1][2] = q[18];
  spatial_metric[2][0] = q[19]; spatial_metric[2][1] = q[20]; spatial_metric[2][2] = q[21];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_spatial_metric)[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }
}

void
gkyl_gr_mhd_tetrad_stress_energy_tensor(double gas_gamma, const double q[75], double ***stress_energy)
{
  double v[75] = { 0.0 };
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];
  double psi = v[8];

  double lapse = v[9];
  double shift_x = v[10];
  double shift_y = v[11];
  double shift_z = v[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[13]; spatial_metric[0][1] = v[14]; spatial_metric[0][2] = v[15];
  spatial_metric[1][0] = v[16]; spatial_metric[1][1] = v[17]; spatial_metric[1][2] = v[18];
  spatial_metric[2][0] = v[19]; spatial_metric[2][1] = v[20]; spatial_metric[2][2] = v[21];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_tetrad_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
    for (int i = 0; i < 3; i++) {
      cov_vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }
    
    double b0 = 0.0;
    for (int i = 0; i < 3; i++) {
      b0 += W * mag[i] * (cov_vel[i] / lapse);
    }

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double cov_b[3];
    for (int i = 0; i < 3; i++) {
      cov_b[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_b[i] += spatial_metric[i][j] * b[j];
      }
    }

    double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
    double p_star = p + (0.5 * b_sq);

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double spacetime_b[4];
    spacetime_b[0] = b0;
    for (int i = 0; i < 3; i++) {
      spacetime_b[i + 1] = b[i];
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = (rho * h_star * spacetime_vel[i] * spacetime_vel[j]) + (p_star * inv_spacetime_metric[i][j]) - (spacetime_b[i] * spacetime_b[j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

static inline double
gkyl_gr_mhd_tetrad_max_abs_speed(double gas_gamma, const double q[75])
{
  double v[75] = { 0.0 };
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];
  double psi = v[8];

  double lapse = v[9];
  double shift_x = v[10];
  double shift_y = v[11];
  double shift_z = v[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[13]; spatial_metric[0][1] = v[14]; spatial_metric[0][2] = v[15];
  spatial_metric[1][0] = v[16]; spatial_metric[1][1] = v[17]; spatial_metric[1][2] = v[18];
  spatial_metric[2][0] = v[19]; spatial_metric[2][1] = v[20]; spatial_metric[2][2] = v[21];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_tetrad_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
    for (int i = 0; i < 3; i++) {
      cov_vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }
    
    double b0 = 0.0;
    for (int i = 0; i < 3; i++) {
      b0 += W * mag[i] * (cov_vel[i] / lapse);
    }

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
    double C = (rho * h) + b_sq;

    double num = (gas_gamma * p) / rho;
    double den = 1.0 + ((p / rho) * (gas_gamma) / (gas_gamma - 1.0));
    double c_s = sqrt(num / den);

    double entropy_eigs[3];
    double fast_alfven_eigs[3];
    double slow_alfven_eigs[3];
    double fast_magnetosonic_eigs[3];
    double slow_magnetosonic_eigs[3];

    for (int i = 0; i < 3; i++) {
      entropy_eigs[i] = (lapse * vel[i] - shift[i]);

      fast_alfven_eigs[i] = (b[i] + (sqrt(C) * spacetime_vel[i + 1])) / (b0 + (sqrt(C) * spacetime_vel[0]));
      slow_alfven_eigs[i] = (b[i] - (sqrt(C) * spacetime_vel[i + 1])) / (b0 - (sqrt(C) * spacetime_vel[0]));

      fast_magnetosonic_eigs[i] = sqrt((fast_alfven_eigs[i] * fast_alfven_eigs[i]) + ((c_s * c_s) * (1.0 - (fast_alfven_eigs[i] * fast_alfven_eigs[i]))));
      slow_magnetosonic_eigs[i] = sqrt((slow_alfven_eigs[i] * slow_alfven_eigs[i]) + ((c_s * c_s) * (1.0 - (slow_alfven_eigs[i] * slow_alfven_eigs[i]))));
    }

    double max_eig = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(entropy_eigs[i]) > max_eig) {
        max_eig = fabs(entropy_eigs[i]);
      }
      if (fabs(fast_alfven_eigs[i]) > max_eig) {
        max_eig = fabs(fast_alfven_eigs[i]);
      }
      if (fabs(slow_alfven_eigs[i]) > max_eig) {
        max_eig = fabs(slow_alfven_eigs[i]);
      }
      if (fabs(fast_magnetosonic_eigs[i]) > max_eig) {
        max_eig = fabs(fast_magnetosonic_eigs[i]);
      }
      if (fabs(slow_magnetosonic_eigs[i]) > max_eig) {
        max_eig = fabs(slow_magnetosonic_eigs[i]);
      }
    }

    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_spatial_metric[i]);
    }
    gkyl_free(inv_spatial_metric);

    return fabs(v_sq) + max_eig;
  }
  else {
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_spatial_metric[i]);
    }
    gkyl_free(inv_spatial_metric);

    return pow(10.0, -8.0);
  }
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 75; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 75; i++) {
    qout[i] = win[i];
  }
}

static void
gr_mhd_tetrad_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 75; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_mhd_tetrad_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];

  for (int i = 5; i < 75; i++) {
    ghost[i] = skin[i];
  }
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

  qlocal[5] = (qglobal[5] * norm[0]) + (qglobal[6] * norm[1]) + (qglobal[7] * norm[2]);
  qlocal[6] = (qglobal[5] * tau1[0]) + (qglobal[6] * tau1[1]) + (qglobal[7] * tau1[2]);
  qlocal[7] = (qglobal[5] * tau2[0]) + (qglobal[6] * tau2[1]) + (qglobal[7] * tau2[2]);
  qlocal[8] = qglobal[8];

  qlocal[9] = qglobal[9];
  qlocal[10] = (qglobal[10] * norm[0]) + (qglobal[11] * norm[1]) + (qglobal[12] * norm[2]);
  qlocal[11] = (qglobal[10] * tau1[0]) + (qglobal[11] * tau1[1]) + (qglobal[12] * tau1[2]);
  qlocal[12] = (qglobal[10] * tau2[0]) + (qglobal[11] * tau2[1]) + (qglobal[12] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[13] * norm[0]) + (qglobal[14] * norm[1]) + (qglobal[15] * norm[2]);
  r1[1] = (qglobal[13] * tau1[0]) + (qglobal[14] * tau1[1]) + (qglobal[15] * tau1[2]);
  r1[2] = (qglobal[13] * tau2[0]) + (qglobal[14] * tau2[1]) + (qglobal[15] * tau2[2]);

  r2[0] = (qglobal[16] * norm[0]) + (qglobal[17] * norm[1]) + (qglobal[18] * norm[2]);
  r2[1] = (qglobal[16] * tau1[0]) + (qglobal[17] * tau1[1]) + (qglobal[18] * tau1[2]);
  r2[2] = (qglobal[16] * tau2[0]) + (qglobal[17] * tau2[1]) + (qglobal[18] * tau2[2]);

  r3[0] = (qglobal[19] * norm[0]) + (qglobal[20] * norm[1]) + (qglobal[21] * norm[2]);
  r3[1] = (qglobal[19] * tau1[0]) + (qglobal[20] * tau1[1]) + (qglobal[21] * tau1[2]);
  r3[2] = (qglobal[19] * tau2[0]) + (qglobal[20] * tau2[1]) + (qglobal[21] * tau2[2]);

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

  // Rotate spatial metric tensor to local coordinate frame.
  qlocal[13] = v1[0]; qlocal[14] = v1[1]; qlocal[15] = v1[2];
  qlocal[16] = v2[0]; qlocal[17] = v2[1]; qlocal[18] = v2[2];
  qlocal[19] = v3[0]; qlocal[20] = v3[1]; qlocal[21] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[22] * norm[0]) + (qglobal[23] * norm[1]) + (qglobal[24] * norm[2]);
  extr_r1[1] = (qglobal[22] * tau1[0]) + (qglobal[23] * tau1[1]) + (qglobal[24] * tau1[2]);
  extr_r1[2] = (qglobal[22] * tau2[0]) + (qglobal[23] * tau2[1]) + (qglobal[24] * tau2[2]);

  extr_r2[0] = (qglobal[25] * norm[0]) + (qglobal[26] * norm[1]) + (qglobal[27] * norm[2]);
  extr_r2[1] = (qglobal[25] * tau1[0]) + (qglobal[26] * tau1[1]) + (qglobal[27] * tau1[2]);
  extr_r2[2] = (qglobal[25] * tau2[0]) + (qglobal[26] * tau2[1]) + (qglobal[27] * tau2[2]);

  extr_r3[0] = (qglobal[28] * norm[0]) + (qglobal[29] * norm[1]) + (qglobal[30] * norm[2]);
  extr_r3[1] = (qglobal[28] * tau1[0]) + (qglobal[29] * tau1[1]) + (qglobal[30] * tau1[2]);
  extr_r3[2] = (qglobal[28] * tau2[0]) + (qglobal[29] * tau2[1]) + (qglobal[30] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * norm[1]) + (extr_r3[0] * norm[2]);
  inv_v1[1] = (extr_r1[0] * tau1[0]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau1[2]);
  inv_v1[2] = (extr_r1[0] * tau2[0]) + (extr_r2[0] * tau2[1]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * norm[1]) + (extr_r3[1] * norm[2]);
  inv_v2[1] = (extr_r1[1] * tau1[0]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau1[2]);
  inv_v2[2] = (extr_r1[1] * tau2[0]) + (extr_r2[1] * tau2[1]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * norm[1]) + (extr_r3[2] * norm[2]);
  inv_v3[1] = (extr_r1[2] * tau1[0]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau1[2]);
  inv_v3[2] = (extr_r1[2] * tau2[0]) + (extr_r2[2] * tau2[1]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor to local coordinate frame.
  qlocal[22] = inv_v1[0]; qlocal[23] = inv_v1[1]; qlocal[24] = inv_v1[2];
  qlocal[25] = inv_v2[0]; qlocal[26] = inv_v2[1]; qlocal[27] = inv_v2[2];
  qlocal[28] = inv_v3[0]; qlocal[29] = inv_v3[1]; qlocal[30] = inv_v3[2];

  qlocal[31] = qglobal[31];

  qlocal[32] = (qglobal[32] * norm[0]) + (qglobal[33] * norm[1]) + (qglobal[34] * norm[2]);
  qlocal[33] = (qglobal[32] * tau1[0]) + (qglobal[33] * tau1[1]) + (qglobal[34] * tau1[2]);
  qlocal[34] = (qglobal[32] * tau2[0]) + (qglobal[33] * tau2[1]) + (qglobal[34] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qglobal[35] * norm[0]) + (qglobal[36] * norm[1]) + (qglobal[37] * norm[2]);
  shiftder_r1[1] = (qglobal[35] * tau1[0]) + (qglobal[36] * tau1[1]) + (qglobal[37] * tau1[2]);
  shiftder_r1[2] = (qglobal[35] * tau2[0]) + (qglobal[36] * tau2[1]) + (qglobal[37] * tau2[2]);

  shiftder_r2[0] = (qglobal[38] * norm[0]) + (qglobal[39] * norm[1]) + (qglobal[40] * norm[2]);
  shiftder_r2[1] = (qglobal[38] * tau1[0]) + (qglobal[39] * tau1[1]) + (qglobal[40] * tau1[2]);
  shiftder_r2[2] = (qglobal[38] * tau2[0]) + (qglobal[39] * tau2[1]) + (qglobal[40] * tau2[2]);

  shiftder_r3[0] = (qglobal[41] * norm[0]) + (qglobal[42] * norm[1]) + (qglobal[43] * norm[2]);
  shiftder_r3[1] = (qglobal[41] * tau1[0]) + (qglobal[42] * tau1[1]) + (qglobal[43] * tau1[2]);
  shiftder_r3[2] = (qglobal[41] * tau2[0]) + (qglobal[42] * tau2[1]) + (qglobal[43] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * norm[1]) + (shiftder_r3[0] * norm[2]);
  shiftder_v1[1] = (shiftder_r1[0] * tau1[0]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau1[2]);
  shiftder_v1[2] = (shiftder_r1[0] * tau2[0]) + (shiftder_r2[0] * tau2[1]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * norm[1]) + (shiftder_r3[1] * norm[2]);
  shiftder_v2[1] = (shiftder_r1[1] * tau1[0]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau1[2]);
  shiftder_v2[2] = (shiftder_r1[1] * tau2[0]) + (shiftder_r2[1] * tau2[1]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * norm[1]) + (shiftder_r3[2] * norm[2]);
  shiftder_v3[1] = (shiftder_r1[2] * tau1[0]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau1[2]);
  shiftder_v3[2] = (shiftder_r1[2] * tau2[0]) + (shiftder_r2[2] * tau2[1]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative to local coordinate frame.
  qlocal[35] = shiftder_v1[0]; qlocal[36] = shiftder_v1[1]; qlocal[37] = shiftder_v1[2];
  qlocal[38] = shiftder_v2[0]; qlocal[39] = shiftder_v2[1]; qlocal[40] = shiftder_v2[2];
  qlocal[41] = shiftder_v3[0]; qlocal[42] = shiftder_v3[1]; qlocal[43] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qglobal[44] * norm[0]) + (qglobal[45] * norm[1]) + (qglobal[46] * norm[2]);
  r11[1] = (qglobal[44] * tau1[0]) + (qglobal[45] * tau1[1]) + (qglobal[46] * tau1[2]);
  r11[2] = (qglobal[44] * tau2[0]) + (qglobal[45] * tau2[1]) + (qglobal[46] * tau2[2]);

  r12[0] = (qglobal[47] * norm[0]) + (qglobal[48] * norm[1]) + (qglobal[49] * norm[2]);
  r12[1] = (qglobal[47] * tau1[0]) + (qglobal[48] * tau1[1]) + (qglobal[49] * tau1[2]);
  r12[2] = (qglobal[47] * tau2[0]) + (qglobal[48] * tau2[1]) + (qglobal[49] * tau2[2]);

  r13[0] = (qglobal[50] * norm[0]) + (qglobal[51] * norm[1]) + (qglobal[52] * norm[2]);
  r13[1] = (qglobal[50] * tau1[0]) + (qglobal[51] * tau1[1]) + (qglobal[52] * tau1[2]);
  r13[2] = (qglobal[50] * tau2[0]) + (qglobal[51] * tau2[1]) + (qglobal[52] * tau2[2]);

  r21[0] = (qglobal[53] * norm[0]) + (qglobal[54] * norm[1]) + (qglobal[55] * norm[2]);
  r21[1] = (qglobal[53] * tau1[0]) + (qglobal[54] * tau1[1]) + (qglobal[55] * tau1[2]);
  r21[2] = (qglobal[53] * tau2[0]) + (qglobal[54] * tau2[1]) + (qglobal[55] * tau2[2]);

  r22[0] = (qglobal[56] * norm[0]) + (qglobal[57] * norm[1]) + (qglobal[58] * norm[2]);
  r22[1] = (qglobal[56] * tau1[0]) + (qglobal[57] * tau1[1]) + (qglobal[58] * tau1[2]);
  r22[2] = (qglobal[56] * tau2[0]) + (qglobal[57] * tau2[1]) + (qglobal[58] * tau2[2]);

  r23[0] = (qglobal[59] * norm[0]) + (qglobal[60] * norm[1]) + (qglobal[61] * norm[2]);
  r23[1] = (qglobal[59] * tau1[0]) + (qglobal[60] * tau1[1]) + (qglobal[61] * tau1[2]);
  r23[2] = (qglobal[59] * tau2[0]) + (qglobal[60] * tau2[1]) + (qglobal[61] * tau2[2]);

  r31[0] = (qglobal[62] * norm[0]) + (qglobal[63] * norm[1]) + (qglobal[64] * norm[2]);
  r31[1] = (qglobal[62] * tau1[0]) + (qglobal[63] * tau1[1]) + (qglobal[64] * tau1[2]);
  r31[2] = (qglobal[62] * tau2[0]) + (qglobal[63] * tau2[1]) + (qglobal[64] * tau2[2]);

  r32[0] = (qglobal[65] * norm[0]) + (qglobal[66] * norm[1]) + (qglobal[67] * norm[2]);
  r32[1] = (qglobal[65] * tau1[0]) + (qglobal[66] * tau1[1]) + (qglobal[67] * tau1[2]);
  r32[2] = (qglobal[65] * tau2[0]) + (qglobal[66] * tau2[1]) + (qglobal[67] * tau2[2]);

  r33[0] = (qglobal[68] * norm[0]) + (qglobal[69] * norm[1]) + (qglobal[70] * norm[2]);
  r33[1] = (qglobal[68] * tau1[0]) + (qglobal[69] * tau1[1]) + (qglobal[70] * tau1[2]);
  r33[2] = (qglobal[68] * tau2[0]) + (qglobal[69] * tau2[1]) + (qglobal[70] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * norm[1]) + (r13[0] * norm[2]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * norm[1]) + (r13[1] * norm[2]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * norm[1]) + (r13[2] * norm[2]);

  s12[0] = (r11[0] * tau1[0]) + (r12[0] * tau1[1]) + (r13[0] * tau1[2]);
  s12[1] = (r11[1] * tau1[0]) + (r12[1] * tau1[1]) + (r13[1] * tau1[2]);
  s12[2] = (r11[2] * tau1[0]) + (r12[2] * tau1[1]) + (r13[2] * tau1[2]);

  s13[0] = (r11[0] * tau2[0]) + (r12[0] * tau2[1]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * tau2[0]) + (r12[1] * tau2[1]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * tau2[0]) + (r12[2] * tau2[1]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * norm[1]) + (r23[0] * norm[2]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * norm[1]) + (r23[1] * norm[2]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * norm[1]) + (r23[2] * norm[2]);

  s22[0] = (r21[0] * tau1[0]) + (r22[0] * tau1[1]) + (r23[0] * tau1[2]);
  s22[1] = (r21[1] * tau1[0]) + (r22[1] * tau1[1]) + (r23[1] * tau1[2]);
  s22[2] = (r21[2] * tau1[0]) + (r22[2] * tau1[1]) + (r23[2] * tau1[2]);

  s23[0] = (r21[0] * tau2[0]) + (r22[0] * tau2[1]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * tau2[0]) + (r22[1] * tau2[1]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * tau2[0]) + (r22[2] * tau2[1]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * norm[1]) + (r33[0] * norm[2]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * norm[1]) + (r33[1] * norm[2]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * norm[1]) + (r33[2] * norm[2]);

  s32[0] = (r31[0] * tau1[0]) + (r32[0] * tau1[1]) + (r33[0] * tau1[2]);
  s32[1] = (r31[1] * tau1[0]) + (r32[1] * tau1[1]) + (r33[1] * tau1[2]);
  s32[2] = (r31[2] * tau1[0]) + (r32[2] * tau1[1]) + (r33[2] * tau1[2]);

  s33[0] = (r31[0] * tau2[0]) + (r32[0] * tau2[1]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * tau2[0]) + (r32[1] * tau2[1]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * tau2[0]) + (r32[2] * tau2[1]) + (r33[2] * tau2[2]);
  
  // Rotate spatial metric tensor derivative to local coordinate frame.
  qlocal[44] = (s11[0] * norm[0]) + (s21[0] * norm[1]) + (s31[0] * norm[2]);
  qlocal[45] = (s11[1] * norm[0]) + (s21[1] * norm[1]) + (s31[1] * norm[2]);
  qlocal[46] = (s11[2] * norm[0]) + (s21[2] * norm[1]) + (s31[2] * norm[2]);

  qlocal[47] = (s12[0] * norm[0]) + (s22[0] * norm[1]) + (s32[0] * norm[2]);
  qlocal[48] = (s12[1] * norm[0]) + (s22[1] * norm[1]) + (s32[1] * norm[2]);
  qlocal[49] = (s12[2] * norm[0]) + (s22[2] * norm[1]) + (s32[2] * norm[2]);

  qlocal[50] = (s13[0] * norm[0]) + (s23[0] * norm[1]) + (s33[0] * norm[2]);
  qlocal[51] = (s13[1] * norm[0]) + (s23[1] * norm[1]) + (s33[1] * norm[2]);
  qlocal[52] = (s13[2] * norm[0]) + (s23[2] * norm[1]) + (s33[2] * norm[2]);

  qlocal[53] = (s11[0] * tau1[0]) + (s21[0] * tau1[1]) + (s31[0] * tau1[2]);
  qlocal[54] = (s11[1] * tau1[0]) + (s21[1] * tau1[1]) + (s31[1] * tau1[2]);
  qlocal[55] = (s11[2] * tau1[0]) + (s21[2] * tau1[1]) + (s31[2] * tau1[2]);

  qlocal[56] = (s12[0] * tau1[0]) + (s22[0] * tau1[1]) + (s32[0] * tau1[2]);
  qlocal[57] = (s12[1] * tau1[0]) + (s22[1] * tau1[1]) + (s32[1] * tau1[2]);
  qlocal[58] = (s12[2] * tau1[0]) + (s22[2] * tau1[1]) + (s32[2] * tau1[2]);

  qlocal[59] = (s13[0] * tau1[0]) + (s23[0] * tau1[1]) + (s33[0] * tau1[2]);
  qlocal[60] = (s13[1] * tau1[0]) + (s23[1] * tau1[1]) + (s33[1] * tau1[2]);
  qlocal[61] = (s13[2] * tau1[0]) + (s23[2] * tau1[1]) + (s33[2] * tau1[2]);

  qlocal[62] = (s11[0] * tau2[0]) + (s21[0] * tau2[1]) + (s31[0] * tau2[2]);
  qlocal[63] = (s11[1] * tau2[0]) + (s21[1] * tau2[1]) + (s31[1] * tau2[2]);
  qlocal[64] = (s11[2] * tau2[0]) + (s21[2] * tau2[1]) + (s31[2] * tau2[2]);

  qlocal[65] = (s12[0] * tau2[0]) + (s22[0] * tau2[1]) + (s32[0] * tau2[2]);
  qlocal[66] = (s12[1] * tau2[0]) + (s22[1] * tau2[1]) + (s32[1] * tau2[2]);
  qlocal[67] = (s12[2] * tau2[0]) + (s22[2] * tau2[1]) + (s32[2] * tau2[2]);

  qlocal[68] = (s13[0] * tau2[0]) + (s23[0] * tau2[1]) + (s33[0] * tau2[2]);
  qlocal[69] = (s13[1] * tau2[0]) + (s23[1] * tau2[1]) + (s33[1] * tau2[2]);
  qlocal[70] = (s13[2] * tau2[0]) + (s23[2] * tau2[1]) + (s33[2] * tau2[2]);

  qlocal[71] = qglobal[71];
  qlocal[72] = (qglobal[72] * norm[0]) + (qglobal[73] * norm[1]) + (qglobal[74] * norm[2]);
  qlocal[73] = (qglobal[72] * tau1[0]) + (qglobal[73] * tau1[1]) + (qglobal[74] * tau1[2]);
  qlocal[74] = (qglobal[72] * tau2[0]) + (qglobal[73] * tau2[1]) + (qglobal[74] * tau2[2]);
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

  qglobal[5] = (qlocal[5] * norm[0]) + (qlocal[6] * tau1[0]) + (qlocal[7] * tau2[0]);
  qglobal[6] = (qlocal[5] * norm[1]) + (qlocal[6] * tau1[1]) + (qlocal[7] * tau2[1]);
  qglobal[7] = (qlocal[5] * norm[2]) + (qlocal[6] * tau1[2]) + (qlocal[7] * tau2[2]);
  qglobal[8] = qlocal[8];

  qglobal[9] = qlocal[9];
  qglobal[10] = (qlocal[10] * norm[0]) + (qlocal[11] * tau1[0]) + (qlocal[12] * tau2[0]);
  qglobal[11] = (qlocal[10] * norm[1]) + (qlocal[11] * tau1[1]) + (qlocal[12] * tau2[1]);
  qglobal[12] = (qlocal[10] * norm[2]) + (qlocal[11] * tau1[2]) + (qlocal[12] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[13] * norm[0]) + (qlocal[14] * tau1[0]) + (qlocal[15] * tau2[0]);
  r1[1] = (qlocal[13] * norm[1]) + (qlocal[14] * tau1[1]) + (qlocal[15] * tau2[1]);
  r1[2] = (qlocal[13] * norm[2]) + (qlocal[14] * tau1[2]) + (qlocal[15] * tau2[2]);

  r2[0] = (qlocal[16] * norm[0]) + (qlocal[17] * tau1[0]) + (qlocal[18] * tau2[0]);
  r2[1] = (qlocal[16] * norm[1]) + (qlocal[17] * tau1[1]) + (qlocal[18] * tau2[1]);
  r2[2] = (qlocal[16] * norm[2]) + (qlocal[17] * tau1[2]) + (qlocal[18] * tau2[2]);

  r3[0] = (qlocal[19] * norm[0]) + (qlocal[20] * tau1[0]) + (qlocal[21] * tau2[0]);
  r3[1] = (qlocal[19] * norm[1]) + (qlocal[20] * tau1[1]) + (qlocal[21] * tau2[1]);
  r3[2] = (qlocal[19] * norm[2]) + (qlocal[20] * tau1[2]) + (qlocal[21] * tau2[2]);

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

  // Rotate spatial metric tensor back to global coordinate frame.
  qglobal[13] = v1[0]; qglobal[14] = v1[1]; qglobal[15] = v1[2];
  qglobal[16] = v2[0]; qglobal[17] = v2[1]; qglobal[18] = v2[2];
  qglobal[19] = v3[0]; qglobal[20] = v3[1]; qglobal[21] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[22] * norm[0]) + (qlocal[23] * tau1[0]) + (qlocal[24] * tau2[0]);
  extr_r1[1] = (qlocal[22] * norm[1]) + (qlocal[23] * tau1[1]) + (qlocal[24] * tau2[1]);
  extr_r1[2] = (qlocal[22] * norm[2]) + (qlocal[23] * tau1[2]) + (qlocal[24] * tau2[2]);

  extr_r2[0] = (qlocal[25] * norm[0]) + (qlocal[26] * tau1[0]) + (qlocal[27] * tau2[0]);
  extr_r2[1] = (qlocal[25] * norm[1]) + (qlocal[26] * tau1[1]) + (qlocal[27] * tau2[1]);
  extr_r2[2] = (qlocal[25] * norm[2]) + (qlocal[26] * tau1[2]) + (qlocal[27] * tau2[2]);

  extr_r3[0] = (qlocal[28] * norm[0]) + (qlocal[29] * tau1[0]) + (qlocal[30] * tau2[0]);
  extr_r3[1] = (qlocal[28] * norm[1]) + (qlocal[29] * tau1[1]) + (qlocal[30] * tau2[1]);
  extr_r3[2] = (qlocal[28] * norm[2]) + (qlocal[29] * tau1[2]) + (qlocal[30] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * tau1[0]) + (extr_r3[0] * tau2[0]);
  inv_v1[1] = (extr_r1[0] * norm[1]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau2[1]);
  inv_v1[2] = (extr_r1[0] * norm[2]) + (extr_r2[0] * tau1[2]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * tau1[0]) + (extr_r3[1] * tau2[0]);
  inv_v2[1] = (extr_r1[1] * norm[1]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau2[1]);
  inv_v2[2] = (extr_r1[1] * norm[2]) + (extr_r2[1] * tau1[2]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * tau1[0]) + (extr_r3[2] * tau2[0]);
  inv_v3[1] = (extr_r1[2] * norm[1]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau2[1]);
  inv_v3[2] = (extr_r1[2] * norm[2]) + (extr_r2[2] * tau1[2]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor back to global coordinate frame.
  qglobal[22] = inv_v1[0]; qglobal[23] = inv_v1[1]; qglobal[24] = inv_v1[2];
  qglobal[25] = inv_v2[0]; qglobal[26] = inv_v2[1]; qglobal[27] = inv_v2[2];
  qglobal[28] = inv_v3[0]; qglobal[29] = inv_v3[1]; qglobal[30] = inv_v3[2];

  qglobal[31] = qlocal[31];

  qglobal[32] = (qlocal[32] * norm[0]) + (qlocal[33] * tau1[0]) + (qlocal[34] * tau2[0]);
  qglobal[33] = (qlocal[32] * norm[1]) + (qlocal[33] * tau1[1]) + (qlocal[34] * tau2[1]);
  qglobal[34] = (qlocal[32] * norm[2]) + (qlocal[33] * tau1[2]) + (qlocal[34] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qlocal[35] * norm[0]) + (qlocal[36] * tau1[0]) + (qlocal[37] * tau2[0]);
  shiftder_r1[1] = (qlocal[35] * norm[1]) + (qlocal[36] * tau1[1]) + (qlocal[37] * tau2[1]);
  shiftder_r1[2] = (qlocal[35] * norm[2]) + (qlocal[36] * tau1[2]) + (qlocal[37] * tau2[2]);

  shiftder_r2[0] = (qlocal[38] * norm[0]) + (qlocal[39] * tau1[0]) + (qlocal[40] * tau2[0]);
  shiftder_r2[1] = (qlocal[38] * norm[1]) + (qlocal[39] * tau1[1]) + (qlocal[40] * tau2[1]);
  shiftder_r2[2] = (qlocal[38] * norm[2]) + (qlocal[39] * tau1[2]) + (qlocal[40] * tau2[2]);

  shiftder_r3[0] = (qlocal[41] * norm[0]) + (qlocal[42] * tau1[0]) + (qlocal[43] * tau2[0]);
  shiftder_r3[1] = (qlocal[41] * norm[1]) + (qlocal[42] * tau1[1]) + (qlocal[43] * tau2[1]);
  shiftder_r3[2] = (qlocal[41] * norm[2]) + (qlocal[42] * tau1[2]) + (qlocal[43] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * tau1[0]) + (shiftder_r3[0] * tau2[0]);
  shiftder_v1[1] = (shiftder_r1[0] * norm[1]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau2[1]);
  shiftder_v1[2] = (shiftder_r1[0] * norm[2]) + (shiftder_r2[0] * tau1[2]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * tau1[0]) + (shiftder_r3[1] * tau2[0]);
  shiftder_v2[1] = (shiftder_r1[1] * norm[1]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau2[1]);
  shiftder_v2[2] = (shiftder_r1[1] * norm[2]) + (shiftder_r2[1] * tau1[2]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * tau1[0]) + (shiftder_r3[2] * tau2[0]);
  shiftder_v3[1] = (shiftder_r1[2] * norm[1]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau2[1]);
  shiftder_v3[2] = (shiftder_r1[2] * norm[2]) + (shiftder_r2[2] * tau1[2]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative back to global coordinate frame.
  qglobal[35] = shiftder_v1[0]; qglobal[36] = shiftder_v1[1]; qglobal[37] = shiftder_v1[2];
  qglobal[38] = shiftder_v2[0]; qglobal[39] = shiftder_v2[1]; qglobal[40] = shiftder_v2[2];
  qglobal[41] = shiftder_v3[0]; qglobal[42] = shiftder_v3[1]; qglobal[43] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qlocal[44] * norm[0]) + (qlocal[45] * tau1[0]) + (qlocal[46] * tau2[0]);
  r11[1] = (qlocal[44] * norm[1]) + (qlocal[45] * tau1[1]) + (qlocal[46] * tau2[1]);
  r11[2] = (qlocal[44] * norm[2]) + (qlocal[45] * tau1[2]) + (qlocal[46] * tau2[2]);

  r12[0] = (qlocal[47] * norm[0]) + (qlocal[48] * tau1[0]) + (qlocal[49] * tau2[0]);
  r12[1] = (qlocal[47] * norm[1]) + (qlocal[48] * tau1[1]) + (qlocal[49] * tau2[1]);
  r12[2] = (qlocal[47] * norm[2]) + (qlocal[48] * tau1[2]) + (qlocal[49] * tau2[2]);

  r13[0] = (qlocal[50] * norm[0]) + (qlocal[51] * tau1[0]) + (qlocal[52] * tau2[0]);
  r13[1] = (qlocal[50] * norm[1]) + (qlocal[51] * tau1[1]) + (qlocal[52] * tau2[1]);
  r13[2] = (qlocal[50] * norm[2]) + (qlocal[51] * tau1[2]) + (qlocal[52] * tau2[2]);

  r21[0] = (qlocal[53] * norm[0]) + (qlocal[54] * tau1[0]) + (qlocal[55] * tau2[0]);
  r21[1] = (qlocal[53] * norm[1]) + (qlocal[54] * tau1[1]) + (qlocal[55] * tau2[1]);
  r21[2] = (qlocal[53] * norm[2]) + (qlocal[54] * tau1[2]) + (qlocal[55] * tau2[2]);

  r22[0] = (qlocal[56] * norm[0]) + (qlocal[57] * tau1[0]) + (qlocal[58] * tau2[0]);
  r22[1] = (qlocal[56] * norm[1]) + (qlocal[57] * tau1[1]) + (qlocal[58] * tau2[1]);
  r22[2] = (qlocal[56] * norm[2]) + (qlocal[57] * tau1[2]) + (qlocal[58] * tau2[2]);

  r23[0] = (qlocal[59] * norm[0]) + (qlocal[60] * tau1[0]) + (qlocal[61] * tau2[0]);
  r23[1] = (qlocal[59] * norm[1]) + (qlocal[60] * tau1[1]) + (qlocal[61] * tau2[1]);
  r23[2] = (qlocal[59] * norm[2]) + (qlocal[60] * tau1[2]) + (qlocal[61] * tau2[2]);

  r31[0] = (qlocal[62] * norm[0]) + (qlocal[63] * tau1[0]) + (qlocal[64] * tau2[0]);
  r31[1] = (qlocal[62] * norm[1]) + (qlocal[63] * tau1[1]) + (qlocal[64] * tau2[1]);
  r31[2] = (qlocal[62] * norm[2]) + (qlocal[63] * tau1[2]) + (qlocal[64] * tau2[2]);

  r32[0] = (qlocal[65] * norm[0]) + (qlocal[66] * tau1[0]) + (qlocal[67] * tau2[0]);
  r32[1] = (qlocal[65] * norm[1]) + (qlocal[66] * tau1[1]) + (qlocal[67] * tau2[1]);
  r32[2] = (qlocal[65] * norm[2]) + (qlocal[66] * tau1[2]) + (qlocal[67] * tau2[2]);

  r33[0] = (qlocal[68] * norm[0]) + (qlocal[69] * tau1[0]) + (qlocal[70] * tau2[0]);
  r33[1] = (qlocal[68] * norm[1]) + (qlocal[69] * tau1[1]) + (qlocal[70] * tau2[1]);
  r33[2] = (qlocal[68] * norm[2]) + (qlocal[69] * tau1[2]) + (qlocal[70] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * tau1[0]) + (r13[0] * tau2[0]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * tau1[0]) + (r13[1] * tau2[0]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * tau1[0]) + (r13[2] * tau2[0]);

  s12[0] = (r11[0] * norm[1]) + (r12[0] * tau1[1]) + (r13[0] * tau2[1]);
  s12[1] = (r11[1] * norm[1]) + (r12[1] * tau1[1]) + (r13[1] * tau2[1]);
  s12[2] = (r11[2] * norm[1]) + (r12[2] * tau1[1]) + (r13[2] * tau2[1]);

  s13[0] = (r11[0] * norm[2]) + (r12[0] * tau1[2]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * norm[2]) + (r12[1] * tau1[2]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * norm[2]) + (r12[2] * tau1[2]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * tau1[0]) + (r23[0] * tau2[0]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * tau1[0]) + (r23[1] * tau2[0]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * tau1[0]) + (r23[2] * tau2[0]);

  s22[0] = (r21[0] * norm[1]) + (r22[0] * tau1[1]) + (r23[0] * tau2[1]);
  s22[1] = (r21[1] * norm[1]) + (r22[1] * tau1[1]) + (r23[1] * tau2[1]);
  s22[2] = (r21[2] * norm[1]) + (r22[2] * tau1[1]) + (r23[2] * tau2[1]);

  s23[0] = (r21[0] * norm[2]) + (r22[0] * tau1[2]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * norm[2]) + (r22[1] * tau1[2]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * norm[2]) + (r22[2] * tau1[2]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * tau1[0]) + (r33[0] * tau2[0]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * tau1[0]) + (r33[1] * tau2[0]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * tau1[0]) + (r33[2] * tau2[0]);

  s32[0] = (r31[0] * norm[1]) + (r32[0] * tau1[1]) + (r33[0] * tau2[1]);
  s32[1] = (r31[1] * norm[1]) + (r32[1] * tau1[1]) + (r33[1] * tau2[1]);
  s32[2] = (r31[2] * norm[1]) + (r32[2] * tau1[1]) + (r33[2] * tau2[1]);

  s33[0] = (r31[0] * norm[2]) + (r32[0] * tau1[2]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * norm[2]) + (r32[1] * tau1[2]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * norm[2]) + (r32[2] * tau1[2]) + (r33[2] * tau2[2]);

  // Rotate spatial metric tensor derivative back to global coordinate frame.
  qglobal[44] = (s11[0] * norm[0]) + (s21[0] * tau1[0]) + (s31[0] * tau2[0]);
  qglobal[45] = (s11[1] * norm[0]) + (s21[1] * tau1[0]) + (s31[1] * tau2[0]);
  qglobal[46] = (s11[2] * norm[0]) + (s21[2] * tau1[0]) + (s31[2] * tau2[0]);

  qglobal[47] = (s12[0] * norm[0]) + (s22[0] * tau1[0]) + (s32[0] * tau2[0]);
  qglobal[48] = (s12[1] * norm[0]) + (s22[1] * tau1[0]) + (s32[1] * tau2[0]);
  qglobal[49] = (s12[2] * norm[0]) + (s22[2] * tau1[0]) + (s32[2] * tau2[0]);

  qglobal[50] = (s13[0] * norm[0]) + (s23[0] * tau1[0]) + (s33[0] * tau2[0]);
  qglobal[51] = (s13[1] * norm[0]) + (s23[1] * tau1[0]) + (s33[1] * tau2[0]);
  qglobal[52] = (s13[2] * norm[0]) + (s23[2] * tau1[0]) + (s33[2] * tau2[0]);

  qglobal[53] = (s11[0] * norm[1]) + (s21[0] * tau1[1]) + (s31[0] * tau2[1]);
  qglobal[54] = (s11[1] * norm[1]) + (s21[1] * tau1[1]) + (s31[1] * tau2[1]);
  qglobal[55] = (s11[2] * norm[1]) + (s21[2] * tau1[1]) + (s31[2] * tau2[1]);

  qglobal[56] = (s12[0] * norm[1]) + (s22[0] * tau1[1]) + (s32[0] * tau2[1]);
  qglobal[57] = (s12[1] * norm[1]) + (s22[1] * tau1[1]) + (s32[1] * tau2[1]);
  qglobal[58] = (s12[2] * norm[1]) + (s22[2] * tau1[1]) + (s32[2] * tau2[1]);

  qglobal[59] = (s13[0] * norm[1]) + (s23[0] * tau1[1]) + (s33[0] * tau2[1]);
  qglobal[60] = (s13[1] * norm[1]) + (s23[1] * tau1[1]) + (s33[1] * tau2[1]);
  qglobal[61] = (s13[2] * norm[1]) + (s23[2] * tau1[1]) + (s33[2] * tau2[1]);

  qglobal[62] = (s11[0] * norm[2]) + (s21[0] * tau1[2]) + (s31[0] * tau2[2]);
  qglobal[63] = (s11[1] * norm[2]) + (s21[1] * tau1[2]) + (s31[1] * tau2[2]);
  qglobal[64] = (s11[2] * norm[2]) + (s21[2] * tau1[2]) + (s31[2] * tau2[2]);

  qglobal[65] = (s12[0] * norm[2]) + (s22[0] * tau1[2]) + (s32[0] * tau2[2]);
  qglobal[66] = (s12[1] * norm[2]) + (s22[1] * tau1[2]) + (s32[1] * tau2[2]);
  qglobal[67] = (s12[2] * norm[2]) + (s22[2] * tau1[2]) + (s32[2] * tau2[2]);

  qglobal[68] = (s13[0] * norm[2]) + (s23[0] * tau1[2]) + (s33[0] * tau2[2]);
  qglobal[69] = (s13[1] * norm[2]) + (s23[1] * tau1[2]) + (s33[1] * tau2[2]);
  qglobal[70] = (s13[2] * norm[2]) + (s23[2] * tau1[2]) + (s33[2] * tau2[2]);

  qglobal[71] = qlocal[71];
  qglobal[72] = (qlocal[72] * norm[0]) + (qlocal[73] * tau1[0]) + (qlocal[74] * tau2[0]);
  qglobal[73] = (qlocal[72] * norm[1]) + (qlocal[73] * tau1[1]) + (qlocal[74] * tau2[1]);
  qglobal[74] = (qlocal[72] * norm[2]) + (qlocal[73] * tau1[2]) + (qlocal[74] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;
  double light_speed = gr_mhd_tetrad->light_speed;
  double b_fact = gr_mhd_tetrad->b_fact;

  double sl = gkyl_gr_mhd_tetrad_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_mhd_tetrad_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl_sr[75], fr_sr[75];
  gkyl_gr_mhd_tetrad_flux(gas_gamma, light_speed, b_fact, ql, fl_sr);
  gkyl_gr_mhd_tetrad_flux(gas_gamma, light_speed, b_fact, qr, fr_sr);

  double fl_gr[75], fr_gr[75];
  gkyl_gr_mhd_tetrad_flux_correction(gas_gamma, light_speed, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_mhd_tetrad_flux_correction(gas_gamma, light_speed, b_fact, qr, fr_sr, fr_gr);

  bool in_excision_region_l = false;
  if (ql[31] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[31] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[75];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 75; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr_gr[i] - fl_gr[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr_gr[i] - fl_gr[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 75; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[75];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 75; i++) {
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
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;
  double light_speed = gr_mhd_tetrad->light_speed;
  double b_fact = gr_mhd_tetrad->b_fact;

  double vl[75], vr[75];
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];

  double mag_x_l = vl[5];
  double mag_y_l = vl[6];
  double mag_z_l = vl[7];
  double psi_l = vl[8];

  double lapse_l = vl[9];
  double shift_x_l = vl[10];
  double shift_y_l = vl[11];
  double shift_z_l = vl[12];

  double spatial_metric_l[3][3];
  spatial_metric_l[0][0] = vl[13]; spatial_metric_l[0][1] = vl[14]; spatial_metric_l[0][2] = vl[15];
  spatial_metric_l[1][0] = vl[16]; spatial_metric_l[1][1] = vl[17]; spatial_metric_l[1][2] = vl[18];
  spatial_metric_l[2][0] = vl[19]; spatial_metric_l[2][1] = vl[20]; spatial_metric_l[2][2] = vl[21];

  double **inv_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_tetrad_inv_spatial_metric(ql, &inv_spatial_metric_l);

  bool in_excision_region_l = false;
  if (vl[31] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  double vel_l[3];
  double v_sq_l = 0.0;
  vel_l[0] = vx_l; vel_l[1] = vy_l; vel_l[2] = vz_l;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
    }
  }

  double W_l = 1.0 / sqrt(1.0 - v_sq_l);
  if (v_sq_l > 1.0 - pow(10.0, -8.0)) {
    W_l = 1.0 / sqrt(pow(10.0, -8.0));
  }

  double mag_l[3];
  mag_l[0] = mag_x_l; mag_l[1] = mag_y_l; mag_l[2] = mag_z_l;

  double cov_mag_l[3];
  for (int i = 0; i < 3; i++) {
    cov_mag_l[i] = 0.0;

    for (int j = 0; j < 3; j++) {
      cov_mag_l[i] += spatial_metric_l[i][j] * mag_l[j];
    }
  }

  double cov_vel_l[3];
  for (int i = 0; i < 3; i++) {
    cov_vel_l[i] = 0.0;

    for (int j = 0; j < 3; j++) {
      cov_vel_l[i] += spatial_metric_l[i][j] * vel_l[j];
    }
  }

  double b0_l = 0.0;
  for (int i = 0; i < 3; i++) {
    b0_l += W_l * mag_l[i] * (cov_vel_l[i] / lapse_l);
  }

  double shift_l[3];
  shift_l[0] = shift_x_l; shift_l[1] = shift_y_l; shift_l[2] = shift_z_l;

  double spacetime_vel_l[4];
  spacetime_vel_l[0] = W_l / lapse_l;
  for (int i = 0; i < 3; i++) {
    spacetime_vel_l[i + 1] = (W_l * vel_l[i]) - (shift_l[i] * (W_l / lapse_l));
  }

  double b_l[3];
  for (int i = 0; i < 3; i++) {
    b_l[i] = (mag_l[i] + (lapse_l * b0_l * spacetime_vel_l[i + 1])) / W_l;
  }

  double b_sq_l = 0.0;
  for (int i = 0; i < 3; i++) {
    b_sq_l += (mag_l[i] * cov_mag_l[i]) / (W_l * W_l);
  }
  b_sq_l += ((lapse_l * lapse_l) * (b0_l * b0_l)) / (W_l * W_l);

  double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
  double C_l = (rho_l * h_l) + b_sq_l;

  double num_l = (gas_gamma * p_l) / rho_l;
  double den_l = 1.0 + ((p_l / rho_l) * (gas_gamma) / (gas_gamma - 1.0));
  double c_sl = sqrt(num_l / den_l);

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];

  double mag_x_r = vr[5];
  double mag_y_r = vr[6];
  double mag_z_r = vr[7];
  double psi_r = vr[8];

  double lapse_r = vr[9];
  double shift_x_r = vr[10];
  double shift_y_r = vr[11];
  double shift_z_r = vr[12];

  double spatial_metric_r[3][3];
  spatial_metric_r[0][0] = vr[13]; spatial_metric_r[0][1] = vr[14]; spatial_metric_r[0][2] = vr[15];
  spatial_metric_r[1][0] = vr[16]; spatial_metric_r[1][1] = vr[17]; spatial_metric_r[1][2] = vr[18];
  spatial_metric_r[2][0] = vr[19]; spatial_metric_r[2][1] = vr[20]; spatial_metric_r[2][2] = vr[21];

  double **inv_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_tetrad_inv_spatial_metric(ql, &inv_spatial_metric_r);

  bool in_excision_region_r = false;
  if (vr[31] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double vel_r[3];
  double v_sq_r = 0.0;
  vel_r[0] = vx_r; vel_r[1] = vy_r; vel_r[2] = vz_r;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
    }
  }

  double W_r = 1.0 / sqrt(1.0 - v_sq_r);
  if (v_sq_r > 1.0 - pow(10.0, -8.0)) {
    W_r = 1.0 / sqrt(pow(10.0, -8.0));
  }

  double mag_r[3];
  mag_r[0] = mag_x_r; mag_r[1] = mag_y_r; mag_r[2] = mag_z_r;

  double cov_mag_r[3];
  for (int i = 0; i < 3; i++) {
    cov_mag_r[i] = 0.0;

    for (int j = 0; j < 3; j++) {
      cov_mag_r[i] += spatial_metric_r[i][j] * mag_r[j];
    }
  }

  double cov_vel_r[3];
  for (int i = 0; i < 3; i++) {
    cov_vel_r[i] = 0.0;

    for (int j = 0; j < 3; j++) {
      cov_vel_r[i] += spatial_metric_r[i][j] * vel_r[j];
    }
  }

  double b0_r = 0.0;
  for (int i = 0; i < 3; i++) {
    b0_r += W_r * mag_r[i] * (cov_vel_r[i] / lapse_r);
  }

  double shift_r[3];
  shift_r[0] = shift_x_r; shift_r[1] = shift_y_r; shift_r[2] = shift_z_r;

  double spacetime_vel_r[4];
  spacetime_vel_r[0] = W_r / lapse_r;
  for (int i = 0; i < 3; i++) {
    spacetime_vel_r[i + 1] = (W_r * vel_r[i]) - (shift_r[i] * (W_r / lapse_r));
  }

  double b_r[3];
  for (int i = 0; i < 3; i++) {
    b_r[i] = (mag_r[i] + (lapse_r * b0_r * spacetime_vel_r[i + 1])) / W_r;
  }

  double b_sq_r = 0.0;
  for (int i = 0; i < 3; i++) {
    b_sq_r += (mag_r[i] * cov_mag_r[i]) / (W_r * W_r);
  }
  b_sq_r += ((lapse_r * lapse_r) * (b0_r * b0_r)) / (W_r * W_r);

  double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));
  double C_r = (rho_r * h_r) + b_sq_r;

  double num_r = (gas_gamma * p_r) / rho_r;
  double den_r = 1.0 + ((p_r / rho_r) * (gas_gamma) / (gas_gamma - 1.0));
  double c_sr = sqrt(num_r / den_r);

  double fast_alfven_eig_l = (b_l[0] + (sqrt(C_l) * spacetime_vel_l[1])) / (b0_l + (sqrt(C_l) * spacetime_vel_l[0]));
  double slow_alfven_eig_l = (b_l[0] - (sqrt(C_l) * spacetime_vel_l[1])) / (b0_l - (sqrt(C_l) * spacetime_vel_l[0]));
  
  double fast_magnetosonic_eig_l = sqrt((fast_alfven_eig_l * fast_alfven_eig_l) + ((c_sl * c_sl) * (1.0 - (fast_alfven_eig_l * fast_alfven_eig_l))));
  double slow_magnetosonic_eig_l = sqrt((slow_alfven_eig_l * slow_alfven_eig_l) + ((c_sl * c_sl) * (1.0 - (slow_alfven_eig_l * slow_alfven_eig_l))));

  double max_eig_l = fmax(fast_magnetosonic_eig_l, slow_magnetosonic_eig_l);

  double fast_alfven_eig_r = (b_r[0] + (sqrt(C_r) * spacetime_vel_r[1])) / (b0_r + (sqrt(C_r) * spacetime_vel_r[0]));
  double slow_alfven_eig_r = (b_r[0] - (sqrt(C_r) * spacetime_vel_r[1])) / (b0_r - (sqrt(C_r) * spacetime_vel_r[0]));
  
  double fast_magnetosonic_eig_r = sqrt((fast_alfven_eig_r * fast_alfven_eig_r) + ((c_sr * c_sr) * (1.0 - (fast_alfven_eig_r * fast_alfven_eig_r))));
  double slow_magnetosonic_eig_r = sqrt((slow_alfven_eig_r * slow_alfven_eig_r) + ((c_sr * c_sr) * (1.0 - (slow_alfven_eig_r * slow_alfven_eig_r))));

  double max_eig_r = fmax(fast_magnetosonic_eig_r, slow_magnetosonic_eig_r);

  double vx_avg = 0.5 * (vx_l + vx_r);
  double max_eig_avg = 0.5 * (max_eig_l + max_eig_r);

  double sl = (vx_avg - max_eig_avg) / (1.0 - (vx_avg * max_eig_avg));
  double sr = (vx_avg + max_eig_avg) / (1.0 + (vx_avg * max_eig_avg));

  double fl_sr[75], fr_sr[75];
  gkyl_gr_mhd_tetrad_flux(gas_gamma, light_speed, b_fact, ql, fl_sr);
  gkyl_gr_mhd_tetrad_flux(gas_gamma, light_speed, b_fact, qr, fr_sr);

  double fl_gr[75], fr_gr[75];
  gkyl_gr_mhd_tetrad_flux_correction(gas_gamma, light_speed, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_mhd_tetrad_flux_correction(gas_gamma, light_speed, b_fact, qr, fr_sr, fr_gr);

  double qm[75];
  for (int i = 0; i < 75; i++) {
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl_gr[i] - fr_gr[i])) / (sr - sl);
  }

  double *w0 = &waves[0], *w1 = &waves[75];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 75; i++) {
      w0[i] = qm[i] - ql[i];
      w1[i] = qr[i] - qm[i];
    }
  }
  else {
    for (int i = 0; i < 75; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = sl;
  s[1] = sr;

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric_l[i]);
    gkyl_free(inv_spatial_metric_r[i]);
  }
  gkyl_free(inv_spatial_metric_l);
  gkyl_free(inv_spatial_metric_r);

  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[75];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 75; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hll(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;
  double light_speed = gr_mhd_tetrad->light_speed;
  double b_fact = gr_mhd_tetrad->b_fact;

  double fr_sr[75], fl_sr[75];
  gkyl_gr_mhd_tetrad_flux(gas_gamma, light_speed, b_fact, ql, fl_sr);
  gkyl_gr_mhd_tetrad_flux(gas_gamma, light_speed, b_fact, qr, fr_sr);

  double fr_gr[75], fl_gr[75];
  gkyl_gr_mhd_tetrad_flux_correction(gas_gamma, light_speed, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_mhd_tetrad_flux_correction(gas_gamma, light_speed, b_fact, qr, fr_sr, fr_gr);

  bool in_excision_region_l = false;
  if (ql[31] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[31] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 75; m++) {
      flux_jump[m] = fr_gr[m] - fl_gr[m];
    }
  }
  else {
    for (int m = 0; m < 75; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_mhd_tetrad_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_gr_mhd_tetrad_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;

  double v[75] = { 0.0 };
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, q, v);

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
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;

  return gkyl_gr_mhd_tetrad_max_abs_speed(gas_gamma, q);
}

static inline void
gr_mhd_tetrad_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 8; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_mhd_tetrad_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;

  double v[75] = { 0.0 };
  gkyl_gr_mhd_tetrad_prim_vars(gas_gamma, qin, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];
  double psi = v[8];

  double lapse = v[9];
  double shift_x = v[10];
  double shift_y = v[11];
  double shift_z = v[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[13]; spatial_metric[0][1] = v[14]; spatial_metric[0][2] = v[15];
  spatial_metric[1][0] = v[16]; spatial_metric[1][1] = v[17]; spatial_metric[1][2] = v[18];
  spatial_metric[2][0] = v[19]; spatial_metric[2][1] = v[20]; spatial_metric[2][2] = v[21];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = v[22]; extrinsic_curvature[0][1] = v[23]; extrinsic_curvature[0][2] = v[24];
  extrinsic_curvature[1][0] = v[25]; extrinsic_curvature[1][1] = v[26]; extrinsic_curvature[1][2] = v[27];
  extrinsic_curvature[2][0] = v[28]; extrinsic_curvature[2][1] = v[29]; extrinsic_curvature[2][2] = v[30];

  double lapse_der[3];
  lapse_der[0] = v[32];
  lapse_der[1] = v[33];
  lapse_der[2] = v[34];

  double shift_der[3][3];
  shift_der[0][0] = v[35]; shift_der[0][1] = v[36]; shift_der[0][2] = v[37];
  shift_der[1][0] = v[38]; shift_der[1][1] = v[39]; shift_der[1][2] = v[40];
  shift_der[2][0] = v[41]; shift_der[2][1] = v[42]; shift_der[2][2] = v[43];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = v[44]; spatial_metric_der[0][0][1] = v[45]; spatial_metric_der[0][0][2] = v[46];
  spatial_metric_der[0][1][0] = v[47]; spatial_metric_der[0][1][1] = v[48]; spatial_metric_der[0][1][2] = v[49];
  spatial_metric_der[0][2][0] = v[50]; spatial_metric_der[0][2][1] = v[51]; spatial_metric_der[0][2][2] = v[52];

  spatial_metric_der[1][0][0] = v[53]; spatial_metric_der[1][0][1] = v[54]; spatial_metric_der[1][0][2] = v[55];
  spatial_metric_der[1][1][0] = v[56]; spatial_metric_der[1][1][1] = v[57]; spatial_metric_der[1][1][2] = v[58];
  spatial_metric_der[1][2][0] = v[59]; spatial_metric_der[1][2][1] = v[60]; spatial_metric_der[1][2][2] = v[61];

  spatial_metric_der[0][0][0] = v[62]; spatial_metric_der[0][0][1] = v[63]; spatial_metric_der[0][0][2] = v[64];
  spatial_metric_der[0][1][0] = v[65]; spatial_metric_der[0][1][1] = v[66]; spatial_metric_der[0][1][2] = v[67];
  spatial_metric_der[0][2][0] = v[68]; spatial_metric_der[0][2][1] = v[69]; spatial_metric_der[0][2][2] = v[70];

  double **stress_energy = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    stress_energy[i] = gkyl_malloc(sizeof(double[4]));
  }

  gkyl_gr_mhd_tetrad_stress_energy_tensor(gas_gamma, qin, &stress_energy);

  bool in_excision_region = false;
  if (v[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
    for (int i = 0; i < 3; i++) {
      cov_vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }
    
    double b0 = 0.0;
    for (int i = 0; i < 3; i++) {
      b0 += W * mag[i] * (cov_vel[i] / lapse);
    }

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double cov_b[3];
    for (int i = 0; i < 3; i++) {
      cov_b[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_b[i] += spatial_metric[i][j] * b[j];
      }
    }

    double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);

    double mom[3];
    for (int i = 0; i < 3; i++) {
      mom[i] = (rho * h_star * (W * W) * cov_vel[i]) - (lapse * b0 * cov_b[i]);
    }

    // Energy density source.
    sout[4] = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sout[4] += (lapse * sqrt(spatial_det)) * (stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        sout[4] += (lapse * sqrt(spatial_det)) * (2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        sout[4] += (lapse * sqrt(spatial_det)) * (stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      sout[4] -= (lapse * sqrt(spatial_det)) * (stress_energy[0][0] * shift[i] * lapse_der[i]);
      sout[4] -= (lapse * sqrt(spatial_det)) * (stress_energy[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      sout[1 + j] = -(lapse * sqrt(spatial_det)) * (stress_energy[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          sout[1 + j] += (lapse * sqrt(spatial_det)) * (0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          sout[1 + j] += (lapse * sqrt(spatial_det)) * (0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        sout[1 + j] += (lapse * sqrt(spatial_det)) * ((mom[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          sout[1 + j] += (lapse * sqrt(spatial_det)) * (stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 75; i++) {
      sout[i] = 0.0;
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(stress_energy[i]);
  }
  gkyl_free(stress_energy);
}

void
gkyl_gr_mhd_tetrad_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(base->on_dev, struct wv_gr_mhd_tetrad, eqn);
    gkyl_cu_free(gr_mhd_tetrad);
  }

  struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(base, struct wv_gr_mhd_tetrad, eqn);
  gkyl_free(gr_mhd_tetrad);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_mhd_tetrad_new(double gas_gamma, double light_speed, double b_fact, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq,
  struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_mhd_tetrad_inew(&(struct gkyl_wv_gr_mhd_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .light_speed = light_speed,
      .b_fact = b_fact,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_MHD_TETRAD_RP_HLL,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_mhd_tetrad_inew(const struct gkyl_wv_gr_mhd_tetrad_inp* inp)
{
  struct wv_gr_mhd_tetrad *gr_mhd_tetrad = gkyl_malloc(sizeof(struct wv_gr_mhd_tetrad));

  gr_mhd_tetrad->eqn.type = GKYL_EQN_GR_MHD_TETRAD;
  gr_mhd_tetrad->eqn.num_equations = 75;
  gr_mhd_tetrad->eqn.num_diag = 8;

  gr_mhd_tetrad->gas_gamma = inp->gas_gamma;
  gr_mhd_tetrad->light_speed = inp->light_speed;
  gr_mhd_tetrad->b_fact = inp->b_fact;
  gr_mhd_tetrad->spacetime_gauge = inp->spacetime_gauge;
  gr_mhd_tetrad->reinit_freq = inp->reinit_freq;
  gr_mhd_tetrad->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_MHD_TETRAD_RP_LAX) {
    gr_mhd_tetrad->eqn.num_waves = 2;
    gr_mhd_tetrad->eqn.waves_func = wave_lax_l;
    gr_mhd_tetrad->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_MHD_TETRAD_RP_HLL) {
    gr_mhd_tetrad->eqn.num_waves = 2;
    gr_mhd_tetrad->eqn.waves_func = wave_hll_l;
    gr_mhd_tetrad->eqn.qfluct_func = qfluct_hll_l;
  }

  gr_mhd_tetrad->eqn.flux_jump = flux_jump;
  gr_mhd_tetrad->eqn.check_inv_func = check_inv;
  gr_mhd_tetrad->eqn.max_speed_func = max_speed;
  gr_mhd_tetrad->eqn.rotate_to_local_func = rot_to_local;
  gr_mhd_tetrad->eqn.rotate_to_global_func = rot_to_global;

  gr_mhd_tetrad->eqn.wall_bc_func = gr_mhd_tetrad_wall;
  gr_mhd_tetrad->eqn.no_slip_bc_func = gr_mhd_tetrad_no_slip;

  gr_mhd_tetrad->eqn.cons_to_riem = cons_to_riem;
  gr_mhd_tetrad->eqn.riem_to_cons = riem_to_cons;

  gr_mhd_tetrad->eqn.cons_to_diag = gr_mhd_tetrad_cons_to_diag;

  gr_mhd_tetrad->eqn.source_func = gr_mhd_tetrad_source;

  gr_mhd_tetrad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_mhd_tetrad->eqn.flags);
  gr_mhd_tetrad->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_mhd_tetrad_free);
  gr_mhd_tetrad->eqn.on_dev = &gr_mhd_tetrad->eqn; // On the CPU, the equation object points to itself.

  gr_mhd_tetrad->eqn.embed_geo = NULL;

  return &gr_mhd_tetrad->eqn;
}

double
gkyl_wv_gr_mhd_tetrad_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double gas_gamma = gr_mhd_tetrad->gas_gamma;

  return gas_gamma;
}

double
gkyl_wv_gr_mhd_tetrad_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double light_speed = gr_mhd_tetrad->light_speed;

  return light_speed;
}

double
gkyl_wv_gr_mhd_tetrad_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  double b_fact = gr_mhd_tetrad->b_fact;

  return b_fact;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_mhd_tetrad_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_mhd_tetrad->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_mhd_tetrad_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  int reinit_freq = gr_mhd_tetrad->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_mhd_tetrad_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_mhd_tetrad->spacetime;

  return spacetime;
}
