#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>

void
gkyl_gr_euler_flux(double gas_gamma, const double q[71], double flux[71])
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
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

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    flux[0] = (lapse * sqrt(spatial_det)) * (rho * W * (vx - (shift_x / lapse)));
    flux[1] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vx * (vx - (shift_x / lapse))) + p);
    flux[2] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vy * (vx - (shift_x / lapse))));
    flux[3] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vz * (vx - (shift_x / lapse))));
    flux[4] = (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vx - (shift_x / lapse)) + (p * vx));

    for (int i = 5; i < 71; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      flux[i] = 0.0;
    }
  }
}

void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[71], double v[71])
{
  double lapse = q[5];
  double shift_x = q[6];
  double shift_y = q[7];
  double shift_z = q[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = q[9]; spatial_metric[0][1] = q[10]; spatial_metric[0][2] = q[11];
  spatial_metric[1][0] = q[12]; spatial_metric[1][1] = q[13]; spatial_metric[1][2] = q[14];
  spatial_metric[2][0] = q[15]; spatial_metric[2][1] = q[16]; spatial_metric[2][2] = q[17];
  
  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[18]; extrinsic_curvature[0][1] = q[19]; extrinsic_curvature[0][2] = q[20];
  extrinsic_curvature[1][0] = q[21]; extrinsic_curvature[1][1] = q[22]; extrinsic_curvature[1][2] = q[23];
  extrinsic_curvature[2][0] = q[24]; extrinsic_curvature[2][1] = q[25]; extrinsic_curvature[2][2] = q[26];

  double lapse_der[3];
  lapse_der[0] = q[28];
  lapse_der[1] = q[29];
  lapse_der[2] = q[30];

  double shift_der[3][3];
  shift_der[0][0] = q[31]; shift_der[0][1] = q[32]; shift_der[0][2] = q[33];
  shift_der[1][0] = q[34]; shift_der[1][1] = q[35]; shift_der[1][2] = q[36];
  shift_der[2][0] = q[37]; shift_der[2][1] = q[38]; shift_der[2][2] = q[39];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[40]; spatial_metric_der[0][0][1] = q[41]; spatial_metric_der[0][0][2] = q[42];
  spatial_metric_der[0][1][0] = q[43]; spatial_metric_der[0][1][1] = q[44]; spatial_metric_der[0][1][2] = q[45];
  spatial_metric_der[0][2][0] = q[46]; spatial_metric_der[0][2][1] = q[47]; spatial_metric_der[0][2][2] = q[48];

  spatial_metric_der[1][0][0] = q[49]; spatial_metric_der[1][0][1] = q[50]; spatial_metric_der[1][0][2] = q[51];
  spatial_metric_der[1][1][0] = q[52]; spatial_metric_der[1][1][1] = q[53]; spatial_metric_der[1][1][2] = q[54];
  spatial_metric_der[1][2][0] = q[55]; spatial_metric_der[1][2][1] = q[56]; spatial_metric_der[1][2][2] = q[57];

  spatial_metric_der[0][0][0] = q[58]; spatial_metric_der[0][0][1] = q[59]; spatial_metric_der[0][0][2] = q[60];
  spatial_metric_der[0][1][0] = q[61]; spatial_metric_der[0][1][1] = q[62]; spatial_metric_der[0][1][2] = q[63];
  spatial_metric_der[0][2][0] = q[64]; spatial_metric_der[0][2][1] = q[65]; spatial_metric_der[0][2][2] = q[66];

  double evol_param = q[67];
  double x = q[68];
  double y = q[69];
  double z = q[70];

  bool in_excision_region = false;
  if (q[27] < pow(10.0, -8.0)) {
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

    double C = D / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
    double C0 = (D + Etot) / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
    if (((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)) < pow(10.0, -8.0)) {
      C = D / sqrt(pow(10.0, -8.0));
      C0 = (D + Etot) / sqrt(pow(10.0, -8.0));
    }

    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0 * C0) - 1.0) + 1.0 - (C * C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0 * C0) - 1.0;
    double eta = 2.0 * C *((gas_gamma - 1.0) / gas_gamma);

    double guess = 1.0;
    int iter = 0;

    while (iter < 100) {
      double poly = (alpha4 * (guess * guess * guess) * (guess - eta)) + (alpha2 * (guess * guess)) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess) + (4.0 * alpha4 * (guess * guess * guess)) - (3.0 * eta * alpha4 * (guess * guess));

      double guess_new = guess - (poly / poly_der);

      if (fabs(guess - guess_new) < pow(10.0, -8.0)) {
        iter = 100;
      }
      else {
        iter += 1;
        guess = guess_new;
      }
    }

    double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0) / gas_gamma) * ((1.0 - (C * guess)) / ((C0 * C0) * (guess * guess))))));
    double h = 1.0 / (C * guess);

    v[0] = D / W; 
    v[1] = momx / (v[0] * h * (W * W));
    v[2] = momy / (v[0] * h * (W * W));
    v[3] = momz / (v[0] * h * (W * W));
    v[4] = (v[0] * h * (W * W)) - D - Etot;

    if (v[0] < pow(10.0, -8.0)) {
      v[0] = pow(10.0, -8.0);
    }
    if (v[4] < pow(10.0, -8.0)) {
      v[4] = pow(10.0, -8.0);
    }

    v[5] = lapse;
    v[6] = shift_x;
    v[7] = shift_y;
    v[8] = shift_z;

    v[9] = spatial_metric[0][0]; v[10] = spatial_metric[0][1]; v[11] = spatial_metric[0][2];
    v[12] = spatial_metric[1][0]; v[13] = spatial_metric[1][1]; v[14] = spatial_metric[1][2];
    v[15] = spatial_metric[2][0]; v[16] = spatial_metric[2][1]; v[17] = spatial_metric[2][2];

    v[18] = extrinsic_curvature[0][0]; v[19] = extrinsic_curvature[0][1]; v[20] = extrinsic_curvature[0][2];
    v[21] = extrinsic_curvature[1][0]; v[22] = extrinsic_curvature[1][1]; v[23] = extrinsic_curvature[1][2];
    v[24] = extrinsic_curvature[2][0]; v[25] = extrinsic_curvature[2][1]; v[26] = extrinsic_curvature[2][2];

    v[27] = 1.0;

    v[28] = lapse_der[0];
    v[29] = lapse_der[1];
    v[30] = lapse_der[2];

    v[31] = shift_der[0][0]; v[32] = shift_der[0][1]; v[33] = shift_der[0][2];
    v[34] = shift_der[1][0]; v[35] = shift_der[1][1]; v[36] = shift_der[1][2];
    v[37] = shift_der[2][0]; v[38] = shift_der[2][1]; v[39] = shift_der[2][2];

    v[40] = spatial_metric_der[0][0][0]; v[41] = spatial_metric_der[0][0][1]; v[42] = spatial_metric_der[0][0][2];
    v[43] = spatial_metric_der[0][1][0]; v[44] = spatial_metric_der[0][1][1]; v[45] = spatial_metric_der[0][1][2];
    v[46] = spatial_metric_der[0][2][0]; v[47] = spatial_metric_der[0][2][1]; v[48] = spatial_metric_der[0][2][2];

    v[49] = spatial_metric_der[1][0][0]; v[50] = spatial_metric_der[1][0][1]; v[51] = spatial_metric_der[1][0][2];
    v[52] = spatial_metric_der[1][1][0]; v[53] = spatial_metric_der[1][1][1]; v[54] = spatial_metric_der[1][1][2];
    v[55] = spatial_metric_der[1][2][0]; v[56] = spatial_metric_der[1][2][1]; v[57] = spatial_metric_der[1][2][2];

    v[58] = spatial_metric_der[2][0][0]; v[59] = spatial_metric_der[2][0][1]; v[60] = spatial_metric_der[2][0][2];
    v[61] = spatial_metric_der[2][1][0]; v[62] = spatial_metric_der[2][1][1]; v[63] = spatial_metric_der[2][1][2];
    v[64] = spatial_metric_der[2][2][0]; v[65] = spatial_metric_der[2][2][1]; v[66] = spatial_metric_der[2][2][2];

    v[67] = evol_param;
    v[68] = x;
    v[69] = y;
    v[70] = z;
  }
  else {
    for (int i = 0; i < 71; i++) {
      v[i] = 0.0;
    }
    
    v[27] = -1.0;
  }
}

void 
gkyl_gr_euler_inv_spatial_metric(const double q[71], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[9]; spatial_metric[0][1] = q[10]; spatial_metric[0][2] = q[11];
  spatial_metric[1][0] = q[12]; spatial_metric[1][1] = q[13]; spatial_metric[1][2] = q[14];
  spatial_metric[2][0] = q[15]; spatial_metric[2][1] = q[16]; spatial_metric[2][2] = q[17];

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
gkyl_gr_euler_stress_energy_tensor(double gas_gamma, const double q[71], double ***stress_energy)
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];
  double shift_y = v[7];
  double shift_z = v[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
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

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    spacetime_vel[1] = (W * vx) - (shift_x * (W / lapse));
    spacetime_vel[2] = (W * vy) - (shift_y * (W / lapse));
    spacetime_vel[3] = (W * vz) - (shift_z * (W / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

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

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = (rho * h * spacetime_vel[i] * spacetime_vel[j]) + (p * inv_spacetime_metric[i][j]);
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
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[71])
{
  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];
  double shift_y = v[7];
  double shift_z = v[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(q, &inv_spatial_metric);

  double num = (gas_gamma * p) / rho;
  double den = 1.0 + ((p / rho) * (gas_gamma) / (gas_gamma - 1.0));
  double c_s = sqrt(num / den);

  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  bool curved_spacetime = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
      else {
        if (fabs(spatial_metric[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
    }
  }
  if (fabs(lapse - 1.0) > pow(10.0, -8.0) || fabs(shift_x) > pow(10.0, -8.0) || fabs(shift_y) > pow(10.0, -8.0) ||
    fabs(shift_z) > pow(10.0, -8.0)) {
    curved_spacetime = true;
  }

  if (!in_excision_region) {
    if (curved_spacetime) {
      double vel[3];
      double v_sq = 0.0;
      vel[0] = vx; vel[1] = vy; vel[2] = vz;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double shift[3];
      shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

      double material_eigs[3];
      double fast_acoustic_eigs[3];
      double slow_acoustic_eigs[3];

      for (int i = 0; i < 3; i++) {
        material_eigs[i] = (lapse * vel[i]) - shift[i];

        fast_acoustic_eigs[i] = (lapse / (1.0 - (v_sq * (c_s * c_s)))) * ((vel[i] * (1.0 - (c_s * c_s))) +
          (c_s * sqrt((1.0 - v_sq) * (inv_spatial_metric[i][i] * (1.0 - (v_sq * (c_s * c_s))) - (vel[i] * vel[i]) * (1.0 - (c_s * c_s)))))) - shift[i];
        
        slow_acoustic_eigs[i] = (lapse / (1.0 - (v_sq * (c_s * c_s)))) * ((vel[i] * (1.0 - (c_s * c_s))) -
          (c_s * sqrt((1.0 - v_sq) * (inv_spatial_metric[i][i] * (1.0 - (v_sq * (c_s * c_s))) - (vel[i] * vel[i]) * (1.0 - (c_s * c_s)))))) - shift[i];
      }

      double max_eig = 0.0;
      for (int i = 0; i < 3; i++) {
        if (fabs(material_eigs[i]) > max_eig) {
          max_eig = fabs(material_eigs[i]);
        }
        if (fabs(fast_acoustic_eigs[i]) > max_eig) {
          max_eig = fabs(fast_acoustic_eigs[i]);
        }
        if (fabs(slow_acoustic_eigs[i]) > max_eig) {
          max_eig = fabs(slow_acoustic_eigs[i]);
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(inv_spatial_metric);

      return fabs(v_sq) + max_eig;
    }
    else {
      double v_sq = sqrt((vx * vx) + (vy * vy) + (vz * vz));

      for (int i = 0; i < 3; i++) {
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(inv_spatial_metric);

      return fabs(v_sq) + c_s;
    }
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
  for (int i = 0; i < 71; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 71; i++) {
    qout[i] = win[i];
  }
}

static void
gr_euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 71; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];

  for (int i = 5; i < 71; i++) {
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

  qlocal[5] = qglobal[5];
  qlocal[6] = (qglobal[6] * norm[0]) + (qglobal[7] * norm[1]) + (qglobal[8] * norm[2]);
  qlocal[7] = (qglobal[6] * tau1[0]) + (qglobal[7] * tau1[1]) + (qglobal[8] * tau1[2]);
  qlocal[8] = (qglobal[6] * tau2[0]) + (qglobal[7] * tau2[1]) + (qglobal[8] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[9] * norm[0]) + (qglobal[10] * norm[1]) + (qglobal[11] * norm[2]);
  r1[1] = (qglobal[9] * tau1[0]) + (qglobal[10] * tau1[1]) + (qglobal[11] * tau1[2]);
  r1[2] = (qglobal[9] * tau2[0]) + (qglobal[10] * tau2[1]) + (qglobal[11] * tau2[2]);

  r2[0] = (qglobal[12] * norm[0]) + (qglobal[13] * norm[1]) + (qglobal[14] * norm[2]);
  r2[1] = (qglobal[12] * tau1[0]) + (qglobal[13] * tau1[1]) + (qglobal[14] * tau1[2]);
  r2[2] = (qglobal[12] * tau2[0]) + (qglobal[13] * tau2[1]) + (qglobal[14] * tau2[2]);

  r3[0] = (qglobal[15] * norm[0]) + (qglobal[16] * norm[1]) + (qglobal[17] * norm[2]);
  r3[1] = (qglobal[15] * tau1[0]) + (qglobal[16] * tau1[1]) + (qglobal[17] * tau1[2]);
  r3[2] = (qglobal[15] * tau2[0]) + (qglobal[16] * tau2[1]) + (qglobal[17] * tau2[2]);

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
  qlocal[9] = v1[0]; qlocal[10] = v1[1]; qlocal[11] = v1[2];
  qlocal[12] = v2[0]; qlocal[13] = v2[1]; qlocal[14] = v2[2];
  qlocal[15] = v3[0]; qlocal[16] = v3[1]; qlocal[17] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[18] * norm[0]) + (qglobal[19] * norm[1]) + (qglobal[20] * norm[2]);
  extr_r1[1] = (qglobal[18] * tau1[0]) + (qglobal[19] * tau1[1]) + (qglobal[20] * tau1[2]);
  extr_r1[2] = (qglobal[18] * tau2[0]) + (qglobal[19] * tau2[1]) + (qglobal[20] * tau2[2]);

  extr_r2[0] = (qglobal[21] * norm[0]) + (qglobal[22] * norm[1]) + (qglobal[23] * norm[2]);
  extr_r2[1] = (qglobal[21] * tau1[0]) + (qglobal[22] * tau1[1]) + (qglobal[23] * tau1[2]);
  extr_r2[2] = (qglobal[21] * tau2[0]) + (qglobal[22] * tau2[1]) + (qglobal[23] * tau2[2]);

  extr_r3[0] = (qglobal[24] * norm[0]) + (qglobal[25] * norm[1]) + (qglobal[26] * norm[2]);
  extr_r3[1] = (qglobal[24] * tau1[0]) + (qglobal[25] * tau1[1]) + (qglobal[26] * tau1[2]);
  extr_r3[2] = (qglobal[24] * tau2[0]) + (qglobal[25] * tau2[1]) + (qglobal[26] * tau2[2]);

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
  qlocal[18] = inv_v1[0]; qlocal[19] = inv_v1[1]; qlocal[20] = inv_v1[2];
  qlocal[21] = inv_v2[0]; qlocal[22] = inv_v2[1]; qlocal[23] = inv_v2[2];
  qlocal[24] = inv_v3[0]; qlocal[25] = inv_v3[1]; qlocal[26] = inv_v3[2];

  qlocal[27] = qglobal[27];

  qlocal[28] = (qglobal[28] * norm[0]) + (qglobal[29] * norm[1]) + (qglobal[30] * norm[2]);
  qlocal[29] = (qglobal[28] * tau1[0]) + (qglobal[29] * tau1[1]) + (qglobal[30] * tau1[2]);
  qlocal[30] = (qglobal[28] * tau2[0]) + (qglobal[29] * tau2[1]) + (qglobal[30] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qglobal[31] * norm[0]) + (qglobal[32] * norm[1]) + (qglobal[33] * norm[2]);
  shiftder_r1[1] = (qglobal[31] * tau1[0]) + (qglobal[32] * tau1[1]) + (qglobal[33] * tau1[2]);
  shiftder_r1[2] = (qglobal[31] * tau2[0]) + (qglobal[32] * tau2[1]) + (qglobal[33] * tau2[2]);

  shiftder_r2[0] = (qglobal[34] * norm[0]) + (qglobal[35] * norm[1]) + (qglobal[36] * norm[2]);
  shiftder_r2[1] = (qglobal[34] * tau1[0]) + (qglobal[35] * tau1[1]) + (qglobal[36] * tau1[2]);
  shiftder_r2[2] = (qglobal[34] * tau2[0]) + (qglobal[35] * tau2[1]) + (qglobal[36] * tau2[2]);

  shiftder_r3[0] = (qglobal[37] * norm[0]) + (qglobal[38] * norm[1]) + (qglobal[39] * norm[2]);
  shiftder_r3[1] = (qglobal[37] * tau1[0]) + (qglobal[38] * tau1[1]) + (qglobal[39] * tau1[2]);
  shiftder_r3[2] = (qglobal[37] * tau2[0]) + (qglobal[38] * tau2[1]) + (qglobal[39] * tau2[2]);

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
  qlocal[31] = shiftder_v1[0]; qlocal[32] = shiftder_v1[1]; qlocal[33] = shiftder_v1[2];
  qlocal[34] = shiftder_v2[0]; qlocal[35] = shiftder_v2[1]; qlocal[36] = shiftder_v2[2];
  qlocal[37] = shiftder_v3[0]; qlocal[38] = shiftder_v3[1]; qlocal[39] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qglobal[40] * norm[0]) + (qglobal[41] * norm[1]) + (qglobal[42] * norm[2]);
  r11[1] = (qglobal[40] * tau1[0]) + (qglobal[41] * tau1[1]) + (qglobal[42] * tau1[2]);
  r11[2] = (qglobal[40] * tau2[0]) + (qglobal[41] * tau2[1]) + (qglobal[42] * tau2[2]);

  r12[0] = (qglobal[43] * norm[0]) + (qglobal[44] * norm[1]) + (qglobal[45] * norm[2]);
  r12[1] = (qglobal[43] * tau1[0]) + (qglobal[44] * tau1[1]) + (qglobal[45] * tau1[2]);
  r12[2] = (qglobal[43] * tau2[0]) + (qglobal[44] * tau2[1]) + (qglobal[45] * tau2[2]);

  r13[0] = (qglobal[46] * norm[0]) + (qglobal[47] * norm[1]) + (qglobal[48] * norm[2]);
  r13[1] = (qglobal[46] * tau1[0]) + (qglobal[47] * tau1[1]) + (qglobal[48] * tau1[2]);
  r13[2] = (qglobal[46] * tau2[0]) + (qglobal[47] * tau2[1]) + (qglobal[48] * tau2[2]);

  r21[0] = (qglobal[49] * norm[0]) + (qglobal[50] * norm[1]) + (qglobal[51] * norm[2]);
  r21[1] = (qglobal[49] * tau1[0]) + (qglobal[50] * tau1[1]) + (qglobal[51] * tau1[2]);
  r21[2] = (qglobal[49] * tau2[0]) + (qglobal[50] * tau2[1]) + (qglobal[51] * tau2[2]);

  r22[0] = (qglobal[52] * norm[0]) + (qglobal[53] * norm[1]) + (qglobal[54] * norm[2]);
  r22[1] = (qglobal[52] * tau1[0]) + (qglobal[53] * tau1[1]) + (qglobal[54] * tau1[2]);
  r22[2] = (qglobal[52] * tau2[0]) + (qglobal[53] * tau2[1]) + (qglobal[54] * tau2[2]);

  r23[0] = (qglobal[55] * norm[0]) + (qglobal[56] * norm[1]) + (qglobal[57] * norm[2]);
  r23[1] = (qglobal[55] * tau1[0]) + (qglobal[56] * tau1[1]) + (qglobal[57] * tau1[2]);
  r23[2] = (qglobal[55] * tau2[0]) + (qglobal[56] * tau2[1]) + (qglobal[57] * tau2[2]);

  r31[0] = (qglobal[58] * norm[0]) + (qglobal[59] * norm[1]) + (qglobal[60] * norm[2]);
  r31[1] = (qglobal[58] * tau1[0]) + (qglobal[59] * tau1[1]) + (qglobal[60] * tau1[2]);
  r31[2] = (qglobal[58] * tau2[0]) + (qglobal[59] * tau2[1]) + (qglobal[60] * tau2[2]);

  r32[0] = (qglobal[61] * norm[0]) + (qglobal[62] * norm[1]) + (qglobal[63] * norm[2]);
  r32[1] = (qglobal[61] * tau1[0]) + (qglobal[62] * tau1[1]) + (qglobal[63] * tau1[2]);
  r32[2] = (qglobal[61] * tau2[0]) + (qglobal[62] * tau2[1]) + (qglobal[63] * tau2[2]);

  r33[0] = (qglobal[64] * norm[0]) + (qglobal[65] * norm[1]) + (qglobal[66] * norm[2]);
  r33[1] = (qglobal[64] * tau1[0]) + (qglobal[65] * tau1[1]) + (qglobal[66] * tau1[2]);
  r33[2] = (qglobal[64] * tau2[0]) + (qglobal[65] * tau2[1]) + (qglobal[66] * tau2[2]);

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
  qlocal[40] = (s11[0] * norm[0]) + (s21[0] * norm[1]) + (s31[0] * norm[2]);
  qlocal[41] = (s11[1] * norm[0]) + (s21[1] * norm[1]) + (s31[1] * norm[2]);
  qlocal[42] = (s11[2] * norm[0]) + (s21[2] * norm[1]) + (s31[2] * norm[2]);

  qlocal[43] = (s12[0] * norm[0]) + (s22[0] * norm[1]) + (s32[0] * norm[2]);
  qlocal[44] = (s12[1] * norm[0]) + (s22[1] * norm[1]) + (s32[1] * norm[2]);
  qlocal[45] = (s12[2] * norm[0]) + (s22[2] * norm[1]) + (s32[2] * norm[2]);

  qlocal[46] = (s13[0] * norm[0]) + (s23[0] * norm[1]) + (s33[0] * norm[2]);
  qlocal[47] = (s13[1] * norm[0]) + (s23[1] * norm[1]) + (s33[1] * norm[2]);
  qlocal[48] = (s13[2] * norm[0]) + (s23[2] * norm[1]) + (s33[2] * norm[2]);

  qlocal[49] = (s11[0] * tau1[0]) + (s21[0] * tau1[1]) + (s31[0] * tau1[2]);
  qlocal[50] = (s11[1] * tau1[0]) + (s21[1] * tau1[1]) + (s31[1] * tau1[2]);
  qlocal[51] = (s11[2] * tau1[0]) + (s21[2] * tau1[1]) + (s31[2] * tau1[2]);

  qlocal[52] = (s12[0] * tau1[0]) + (s22[0] * tau1[1]) + (s32[0] * tau1[2]);
  qlocal[53] = (s12[1] * tau1[0]) + (s22[1] * tau1[1]) + (s32[1] * tau1[2]);
  qlocal[54] = (s12[2] * tau1[0]) + (s22[2] * tau1[1]) + (s32[2] * tau1[2]);

  qlocal[55] = (s13[0] * tau1[0]) + (s23[0] * tau1[1]) + (s33[0] * tau1[2]);
  qlocal[56] = (s13[1] * tau1[0]) + (s23[1] * tau1[1]) + (s33[1] * tau1[2]);
  qlocal[57] = (s13[2] * tau1[0]) + (s23[2] * tau1[1]) + (s33[2] * tau1[2]);

  qlocal[58] = (s11[0] * tau2[0]) + (s21[0] * tau2[1]) + (s31[0] * tau2[2]);
  qlocal[59] = (s11[1] * tau2[0]) + (s21[1] * tau2[1]) + (s31[1] * tau2[2]);
  qlocal[60] = (s11[2] * tau2[0]) + (s21[2] * tau2[1]) + (s31[2] * tau2[2]);

  qlocal[61] = (s12[0] * tau2[0]) + (s22[0] * tau2[1]) + (s32[0] * tau2[2]);
  qlocal[62] = (s12[1] * tau2[0]) + (s22[1] * tau2[1]) + (s32[1] * tau2[2]);
  qlocal[63] = (s12[2] * tau2[0]) + (s22[2] * tau2[1]) + (s32[2] * tau2[2]);

  qlocal[64] = (s13[0] * tau2[0]) + (s23[0] * tau2[1]) + (s33[0] * tau2[2]);
  qlocal[65] = (s13[1] * tau2[0]) + (s23[1] * tau2[1]) + (s33[1] * tau2[2]);
  qlocal[66] = (s13[2] * tau2[0]) + (s23[2] * tau2[1]) + (s33[2] * tau2[2]);

  qlocal[67] = qglobal[67];
  qlocal[68] = (qglobal[68] * norm[0]) + (qglobal[69] * norm[1]) + (qglobal[70] * norm[2]);
  qlocal[69] = (qglobal[68] * tau1[0]) + (qglobal[69] * tau1[1]) + (qglobal[70] * tau1[2]);
  qlocal[70] = (qglobal[68] * tau2[0]) + (qglobal[69] * tau2[1]) + (qglobal[70] * tau2[2]);
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
  qglobal[6] = (qlocal[6] * norm[0]) + (qlocal[7] * tau1[0]) + (qlocal[8] * tau2[0]);
  qglobal[7] = (qlocal[6] * norm[1]) + (qlocal[7] * tau1[1]) + (qlocal[8] * tau2[1]);
  qglobal[8] = (qlocal[6] * norm[2]) + (qlocal[7] * tau1[2]) + (qlocal[8] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[9] * norm[0]) + (qlocal[10] * tau1[0]) + (qlocal[11] * tau2[0]);
  r1[1] = (qlocal[9] * norm[1]) + (qlocal[10] * tau1[1]) + (qlocal[11] * tau2[1]);
  r1[2] = (qlocal[9] * norm[2]) + (qlocal[10] * tau1[2]) + (qlocal[11] * tau2[2]);

  r2[0] = (qlocal[12] * norm[0]) + (qlocal[13] * tau1[0]) + (qlocal[14] * tau2[0]);
  r2[1] = (qlocal[12] * norm[1]) + (qlocal[13] * tau1[1]) + (qlocal[14] * tau2[1]);
  r2[2] = (qlocal[12] * norm[2]) + (qlocal[13] * tau1[2]) + (qlocal[14] * tau2[2]);

  r3[0] = (qlocal[15] * norm[0]) + (qlocal[16] * tau1[0]) + (qlocal[17] * tau2[0]);
  r3[1] = (qlocal[15] * norm[1]) + (qlocal[16] * tau1[1]) + (qlocal[17] * tau2[1]);
  r3[2] = (qlocal[15] * norm[2]) + (qlocal[16] * tau1[2]) + (qlocal[17] * tau2[2]);

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
  qglobal[9] = v1[0]; qglobal[10] = v1[1]; qglobal[11] = v1[2];
  qglobal[12] = v2[0]; qglobal[13] = v2[1]; qglobal[14] = v2[2];
  qglobal[15] = v3[0]; qglobal[16] = v3[1]; qglobal[17] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[18] * norm[0]) + (qlocal[19] * tau1[0]) + (qlocal[20] * tau2[0]);
  extr_r1[1] = (qlocal[18] * norm[1]) + (qlocal[19] * tau1[1]) + (qlocal[20] * tau2[1]);
  extr_r1[2] = (qlocal[18] * norm[2]) + (qlocal[19] * tau1[2]) + (qlocal[20] * tau2[2]);

  extr_r2[0] = (qlocal[21] * norm[0]) + (qlocal[22] * tau1[0]) + (qlocal[23] * tau2[0]);
  extr_r2[1] = (qlocal[21] * norm[1]) + (qlocal[22] * tau1[1]) + (qlocal[23] * tau2[1]);
  extr_r2[2] = (qlocal[21] * norm[2]) + (qlocal[22] * tau1[2]) + (qlocal[23] * tau2[2]);

  extr_r3[0] = (qlocal[24] * norm[0]) + (qlocal[25] * tau1[0]) + (qlocal[26] * tau2[0]);
  extr_r3[1] = (qlocal[24] * norm[1]) + (qlocal[25] * tau1[1]) + (qlocal[26] * tau2[1]);
  extr_r3[2] = (qlocal[24] * norm[2]) + (qlocal[25] * tau1[2]) + (qlocal[26] * tau2[2]);

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
  qglobal[18] = inv_v1[0]; qglobal[19] = inv_v1[1]; qglobal[20] = inv_v1[2];
  qglobal[21] = inv_v2[0]; qglobal[22] = inv_v2[1]; qglobal[23] = inv_v2[2];
  qglobal[24] = inv_v3[0]; qglobal[25] = inv_v3[1]; qglobal[26] = inv_v3[2];

  qglobal[27] = qlocal[27];

  qglobal[28] = (qlocal[28] * norm[0]) + (qlocal[29] * tau1[0]) + (qlocal[30] * tau2[0]);
  qglobal[29] = (qlocal[28] * norm[1]) + (qlocal[29] * tau1[1]) + (qlocal[30] * tau2[1]);
  qglobal[30] = (qlocal[28] * norm[2]) + (qlocal[29] * tau1[2]) + (qlocal[30] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qlocal[31] * norm[0]) + (qlocal[32] * tau1[0]) + (qlocal[33] * tau2[0]);
  shiftder_r1[1] = (qlocal[31] * norm[1]) + (qlocal[32] * tau1[1]) + (qlocal[33] * tau2[1]);
  shiftder_r1[2] = (qlocal[31] * norm[2]) + (qlocal[32] * tau1[2]) + (qlocal[33] * tau2[2]);

  shiftder_r2[0] = (qlocal[34] * norm[0]) + (qlocal[35] * tau1[0]) + (qlocal[36] * tau2[0]);
  shiftder_r2[1] = (qlocal[34] * norm[1]) + (qlocal[35] * tau1[1]) + (qlocal[36] * tau2[1]);
  shiftder_r2[2] = (qlocal[34] * norm[2]) + (qlocal[35] * tau1[2]) + (qlocal[36] * tau2[2]);

  shiftder_r3[0] = (qlocal[37] * norm[0]) + (qlocal[38] * tau1[0]) + (qlocal[39] * tau2[0]);
  shiftder_r3[1] = (qlocal[37] * norm[1]) + (qlocal[38] * tau1[1]) + (qlocal[39] * tau2[1]);
  shiftder_r3[2] = (qlocal[37] * norm[2]) + (qlocal[38] * tau1[2]) + (qlocal[39] * tau2[2]);

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
  qglobal[31] = shiftder_v1[0]; qglobal[32] = shiftder_v1[1]; qglobal[33] = shiftder_v1[2];
  qglobal[34] = shiftder_v2[0]; qglobal[35] = shiftder_v2[1]; qglobal[36] = shiftder_v2[2];
  qglobal[37] = shiftder_v3[0]; qglobal[38] = shiftder_v3[1]; qglobal[39] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qlocal[40] * norm[0]) + (qlocal[41] * tau1[0]) + (qlocal[42] * tau2[0]);
  r11[1] = (qlocal[40] * norm[1]) + (qlocal[41] * tau1[1]) + (qlocal[42] * tau2[1]);
  r11[2] = (qlocal[40] * norm[2]) + (qlocal[41] * tau1[2]) + (qlocal[42] * tau2[2]);

  r12[0] = (qlocal[43] * norm[0]) + (qlocal[44] * tau1[0]) + (qlocal[45] * tau2[0]);
  r12[1] = (qlocal[43] * norm[1]) + (qlocal[44] * tau1[1]) + (qlocal[45] * tau2[1]);
  r12[2] = (qlocal[43] * norm[2]) + (qlocal[44] * tau1[2]) + (qlocal[45] * tau2[2]);

  r13[0] = (qlocal[46] * norm[0]) + (qlocal[47] * tau1[0]) + (qlocal[48] * tau2[0]);
  r13[1] = (qlocal[46] * norm[1]) + (qlocal[47] * tau1[1]) + (qlocal[48] * tau2[1]);
  r13[2] = (qlocal[46] * norm[2]) + (qlocal[47] * tau1[2]) + (qlocal[48] * tau2[2]);

  r21[0] = (qlocal[49] * norm[0]) + (qlocal[50] * tau1[0]) + (qlocal[51] * tau2[0]);
  r21[1] = (qlocal[49] * norm[1]) + (qlocal[50] * tau1[1]) + (qlocal[51] * tau2[1]);
  r21[2] = (qlocal[49] * norm[2]) + (qlocal[50] * tau1[2]) + (qlocal[51] * tau2[2]);

  r22[0] = (qlocal[52] * norm[0]) + (qlocal[53] * tau1[0]) + (qlocal[54] * tau2[0]);
  r22[1] = (qlocal[52] * norm[1]) + (qlocal[53] * tau1[1]) + (qlocal[54] * tau2[1]);
  r22[2] = (qlocal[52] * norm[2]) + (qlocal[53] * tau1[2]) + (qlocal[54] * tau2[2]);

  r23[0] = (qlocal[55] * norm[0]) + (qlocal[56] * tau1[0]) + (qlocal[57] * tau2[0]);
  r23[1] = (qlocal[55] * norm[1]) + (qlocal[56] * tau1[1]) + (qlocal[57] * tau2[1]);
  r23[2] = (qlocal[55] * norm[2]) + (qlocal[56] * tau1[2]) + (qlocal[57] * tau2[2]);

  r31[0] = (qlocal[58] * norm[0]) + (qlocal[59] * tau1[0]) + (qlocal[60] * tau2[0]);
  r31[1] = (qlocal[58] * norm[1]) + (qlocal[59] * tau1[1]) + (qlocal[60] * tau2[1]);
  r31[2] = (qlocal[58] * norm[2]) + (qlocal[59] * tau1[2]) + (qlocal[60] * tau2[2]);

  r32[0] = (qlocal[61] * norm[0]) + (qlocal[62] * tau1[0]) + (qlocal[63] * tau2[0]);
  r32[1] = (qlocal[61] * norm[1]) + (qlocal[62] * tau1[1]) + (qlocal[63] * tau2[1]);
  r32[2] = (qlocal[61] * norm[2]) + (qlocal[62] * tau1[2]) + (qlocal[63] * tau2[2]);

  r33[0] = (qlocal[64] * norm[0]) + (qlocal[65] * tau1[0]) + (qlocal[66] * tau2[0]);
  r33[1] = (qlocal[64] * norm[1]) + (qlocal[65] * tau1[1]) + (qlocal[66] * tau2[1]);
  r33[2] = (qlocal[64] * norm[2]) + (qlocal[65] * tau1[2]) + (qlocal[66] * tau2[2]);

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
  qglobal[40] = (s11[0] * norm[0]) + (s21[0] * tau1[0]) + (s31[0] * tau2[0]);
  qglobal[41] = (s11[1] * norm[0]) + (s21[1] * tau1[0]) + (s31[1] * tau2[0]);
  qglobal[42] = (s11[2] * norm[0]) + (s21[2] * tau1[0]) + (s31[2] * tau2[0]);

  qglobal[43] = (s12[0] * norm[0]) + (s22[0] * tau1[0]) + (s32[0] * tau2[0]);
  qglobal[44] = (s12[1] * norm[0]) + (s22[1] * tau1[0]) + (s32[1] * tau2[0]);
  qglobal[45] = (s12[2] * norm[0]) + (s22[2] * tau1[0]) + (s32[2] * tau2[0]);

  qglobal[46] = (s13[0] * norm[0]) + (s23[0] * tau1[0]) + (s33[0] * tau2[0]);
  qglobal[47] = (s13[1] * norm[0]) + (s23[1] * tau1[0]) + (s33[1] * tau2[0]);
  qglobal[48] = (s13[2] * norm[0]) + (s23[2] * tau1[0]) + (s33[2] * tau2[0]);

  qglobal[49] = (s11[0] * norm[1]) + (s21[0] * tau1[1]) + (s31[0] * tau2[1]);
  qglobal[50] = (s11[1] * norm[1]) + (s21[1] * tau1[1]) + (s31[1] * tau2[1]);
  qglobal[51] = (s11[2] * norm[1]) + (s21[2] * tau1[1]) + (s31[2] * tau2[1]);

  qglobal[52] = (s12[0] * norm[1]) + (s22[0] * tau1[1]) + (s32[0] * tau2[1]);
  qglobal[53] = (s12[1] * norm[1]) + (s22[1] * tau1[1]) + (s32[1] * tau2[1]);
  qglobal[54] = (s12[2] * norm[1]) + (s22[2] * tau1[1]) + (s32[2] * tau2[1]);

  qglobal[55] = (s13[0] * norm[1]) + (s23[0] * tau1[1]) + (s33[0] * tau2[1]);
  qglobal[56] = (s13[1] * norm[1]) + (s23[1] * tau1[1]) + (s33[1] * tau2[1]);
  qglobal[57] = (s13[2] * norm[1]) + (s23[2] * tau1[1]) + (s33[2] * tau2[1]);

  qglobal[58] = (s11[0] * norm[2]) + (s21[0] * tau1[2]) + (s31[0] * tau2[2]);
  qglobal[59] = (s11[1] * norm[2]) + (s21[1] * tau1[2]) + (s31[1] * tau2[2]);
  qglobal[60] = (s11[2] * norm[2]) + (s21[2] * tau1[2]) + (s31[2] * tau2[2]);

  qglobal[61] = (s12[0] * norm[2]) + (s22[0] * tau1[2]) + (s32[0] * tau2[2]);
  qglobal[62] = (s12[1] * norm[2]) + (s22[1] * tau1[2]) + (s32[1] * tau2[2]);
  qglobal[63] = (s12[2] * norm[2]) + (s22[2] * tau1[2]) + (s32[2] * tau2[2]);

  qglobal[64] = (s13[0] * norm[2]) + (s23[0] * tau1[2]) + (s33[0] * tau2[2]);
  qglobal[65] = (s13[1] * norm[2]) + (s23[1] * tau1[2]) + (s33[1] * tau2[2]);
  qglobal[66] = (s13[2] * norm[2]) + (s23[2] * tau1[2]) + (s33[2] * tau2[2]);

  qglobal[67] = qlocal[67];
  qglobal[68] = (qlocal[68] * norm[0]) + (qlocal[69] * tau1[0]) + (qlocal[70] * tau2[0]);
  qglobal[69] = (qlocal[68] * norm[1]) + (qlocal[69] * tau1[1]) + (qlocal[70] * tau2[1]);
  qglobal[70] = (qlocal[68] * norm[2]) + (qlocal[69] * tau1[2]) + (qlocal[70] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double sl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[71], fr[71];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[27] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[27] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[71];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 71; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[71];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 71; i++) {
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
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double vl[71], vr[71];
  double gas_gamma = gr_euler->gas_gamma;
  
  gkyl_gr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_euler_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];

  double Etot_l = ql[4];
  double Etot_r = qr[4];

  double W_l = 1.0 / sqrt(1.0 - ((vx_l * vx_l) + (vy_l * vy_l) + (vz_l * vz_l)));
  double W_r = 1.0 / sqrt(1.0 - ((vx_r * vx_r) + (vy_r * vy_r) + (vz_r * vz_r)));

  double K_l = sqrt(Etot_l + p_l) / W_l;
  double K_r = sqrt(Etot_r + p_r) / W_r;
  double K_avg = 1.0 / (K_l + K_r);

  double v0 = ((K_l * W_l) + (K_r * W_r)) * K_avg;
  double v1 = ((K_l * W_l * vx_l) + (K_r * W_r * vx_r)) * K_avg;
  double v2 = ((K_l * W_l * vy_l) + (K_r * W_r * vy_r)) * K_avg;
  double v3 = ((K_l * W_l * vz_l) + (K_r * W_r * vz_r)) * K_avg;
  double v4 = ((p_l / K_l) + (p_r / K_r)) * K_avg;

  double c_minus = 1.0 - ((gas_gamma / (gas_gamma - 1.0)) * v4);
  double c_plus = 1.0 + ((gas_gamma / (gas_gamma - 1.0)) * v4);

  double v_alpha_sq = -(v0 * v0) + (v1 * v1) + (v2 * v2) + (v3 * v3);
  double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_alpha_sq)) - (0.5 * (gas_gamma - 1.0) * (1.0 + v_alpha_sq));
  double energy = (v0 * v0) - (v1 * v1);
  double y = sqrt(((1.0 - (gas_gamma * v4)) * energy) + s_sq);

  double k = (v0 * delta[4]) - (v1 * delta[1]);
  double v_delta = (-v0 * delta[4]) + (v1 * delta[1]) + (v2 * delta[2]) + (v3 * delta[3]);
  double a1 = -((s_sq * k) + (sqrt(s_sq) * y * ((v0 * delta[1]) - (v1 * delta[4])) + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta))))) / (2.0 * energy * s_sq);
  double a2 = -((s_sq * k) - (sqrt(s_sq) * y * ((v0 * delta[1]) - (v1 * delta[4])) + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta))))) / (2.0 * energy * s_sq);
  double a3 = ((2.0 * s_sq * k) + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta)))) / (energy * s_sq);
  double a4 = delta[2] - ((k * v2) / energy);
  double a5 = delta[3] - ((k * v3) / energy);

  for (int i = 0; i < 71 * 3; i++) {
    waves[i] = 0;
  }

  double *wv;
  wv = &waves[0 * 71];
  wv[0] = a1 * c_minus;
  wv[1] = a1 * (v1 - ((sqrt(s_sq) * v0) / y));
  wv[2] = a1 * v2;
  wv[3] = a1 * v3;
  wv[4] = a1 * (v0 - ((sqrt(s_sq) * v1) / y));
  s[0] = (((1.0 - (gas_gamma * v4)) * v0 * v1) - (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  wv = &waves[1 * 71];
  wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0)))) - (a4 * c_plus * v2) - (a5 * c_plus * v3);
  wv[1] = a3 * v1;
  wv[2] = (a3 * v2) + a4;
  wv[3] = (a3 * v3) + a5;
  wv[4] = a3 * v0;
  s[1] = v1 / v0;

  wv = &waves[2 * 71];
  wv[0] = a2 * c_minus;
  wv[1] = a2 * (v1 + ((sqrt(s_sq) * v0) / y));
  wv[2] = a2 * v2;
  wv[3] = a2 * v3;
  wv[4] = a2 * (v0 + ((sqrt(s_sq) * v1) / y));
  s[2] = (((1.0 - (gas_gamma * v4)) * v0 * v1) + (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  return (((1.0 - (gas_gamma * v4)) * v0 * fabs(v1)) + (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double* w0 = &waves[0 * 71];
  const double* w1 = &waves[1 * 71];
  const double* w2 = &waves[2 * 71];

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]);
  }
  for (int i = 5; i < 71; i++) {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
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
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double vl[71], vr[71];
  gkyl_gr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_euler_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];

  double lapse_l = vl[5];
  double shift_xl = vl[6];
  double shift_yl = vl[7];
  double shift_zl = vl[8];

  double spatial_metric_l[3][3];
  spatial_metric_l[0][0] = vl[9]; spatial_metric_l[0][1] = vl[10]; spatial_metric_l[0][2] = vl[11];
  spatial_metric_l[1][0] = vl[12]; spatial_metric_l[1][1] = vl[13]; spatial_metric_l[1][2] = vl[14];
  spatial_metric_l[2][0] = vl[15]; spatial_metric_l[2][1] = vl[16]; spatial_metric_l[2][2] = vl[17];

  double **inv_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(ql, &inv_spatial_metric_l);

  bool in_excision_region_l = false;
  if (vl[27] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool curved_spacetime_l = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric_l[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime_l = true;
        }
      }
      else {
        if (fabs(spatial_metric_l[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime_l = true;
        }
      }
    }
  }
  if (fabs(lapse_l - 1.0) > pow(10.0, -8.0) || fabs(shift_xl) > pow(10.0, -8.0) || fabs(shift_yl) > pow(10.0, -8.0) ||
    fabs(shift_zl) > pow(10.0, -8.0)) {
    curved_spacetime_l = true;
  }

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];

  double lapse_r = vr[5];
  double shift_xr = vr[6];
  double shift_yr = vr[7];
  double shift_zr = vr[8];

  double spatial_metric_r[3][3];
  spatial_metric_r[0][0] = vr[9]; spatial_metric_r[0][1] = vr[10]; spatial_metric_r[0][2] = vr[11];
  spatial_metric_r[1][0] = vr[12]; spatial_metric_r[1][1] = vr[13]; spatial_metric_r[1][2] = vr[14];
  spatial_metric_r[2][0] = vr[15]; spatial_metric_r[2][1] = vr[16]; spatial_metric_r[2][2] = vr[17];

  double **inv_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_euler_inv_spatial_metric(qr, &inv_spatial_metric_r);

  bool in_excision_region_r = false;
  if (vr[27] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  bool curved_spacetime_r = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric_r[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime_r = true;
        }
      }
      else {
        if (fabs(spatial_metric_l[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime_r = true;
        }
      }
    }
  }
  if (fabs(lapse_r - 1.0) > pow(10.0, -8.0) || fabs(shift_xr) > pow(10.0, -8.0) || fabs(shift_yr) > pow(10.0, -8.0) ||
    fabs(shift_zr) > pow(10.0, -8.0)) {
    curved_spacetime_r = true;
  }

  double num_l = (gas_gamma * p_l) / rho_l;
  double den_l = 1.0 + ((p_l / rho_l) * (gas_gamma) / (gas_gamma - 1.0));
  double c_sl = sqrt(num_l / den_l);

  double num_r = (gas_gamma * p_r) / rho_r;
  double den_r = 1.0 + ((p_r / rho_r) * (gas_gamma) / (gas_gamma - 1.0));
  double c_sr = sqrt(num_r / den_r);

  double vx_avg = 0.5 * (vx_l + vx_r);
  double cs_avg = 0.5 * (c_sl + c_sr);

  double sl, sr;

  if (curved_spacetime_l || curved_spacetime_r) {
    double vel_l[3];
    double v_sq_l = 0.0;
    vel_l[0] = vx_l; vel_l[1] = vy_l; vel_l[2] = vz_l;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
      }
    }

    double shift_l[3];
    shift_l[0] = shift_xl; shift_l[1] = shift_yl; shift_l[2] = shift_zl;

    double material_eigs_l[3];
    double fast_acoustic_eigs_l[3];
    double slow_acoustic_eigs_l[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_l[i] = (lapse_l * vel_l[i]) - shift_l[i];

      fast_acoustic_eigs_l[i] = (lapse_l / (1.0 - (v_sq_l * (c_sl * c_sl)))) * ((vel_l[i] * (1.0 - (c_sl * c_sl))) +
        (c_sl * sqrt((1.0 - v_sq_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_l * (c_sl * c_sl))) - (vel_l[i] * vel_l[i]) * (1.0 - (c_sl * c_sl)))))) - shift_l[i];
      
      slow_acoustic_eigs_l[i] = (lapse_l / (1.0 - (v_sq_l * (c_sl * c_sl)))) * ((vel_l[i] * (1.0 - (c_sl * c_sl))) -
        (c_sl * sqrt((1.0 - v_sq_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_l * (c_sl * c_sl))) - (vel_l[i] * vel_l[i]) * (1.0 - (c_sl * c_sl)))))) - shift_l[i];
    }

    double max_eig_l = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_l[i]) > max_eig_l) {
        max_eig_l = fabs(material_eigs_l[i]);
      }
      if (fabs(fast_acoustic_eigs_l[i]) > max_eig_l) {
        max_eig_l = fabs(fast_acoustic_eigs_l[i]);
      }
      if (fabs(slow_acoustic_eigs_l[i]) > max_eig_l) {
        max_eig_l = fabs(slow_acoustic_eigs_l[i]);
      }
    }

    double vel_r[3];
    double v_sq_r = 0.0;
    vel_r[0] = vx_r; vel_r[1] = vy_r; vel_r[2] = vz_r;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
      }
    }

    double shift_r[3];
    shift_r[0] = shift_xr; shift_r[1] = shift_yr; shift_r[2] = shift_zr;

    double material_eigs_r[3];
    double fast_acoustic_eigs_r[3];
    double slow_acoustic_eigs_r[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_r[i] = (lapse_r * vel_r[i]) - shift_r[i];

      fast_acoustic_eigs_r[i] = (lapse_r / (1.0 - (v_sq_r * (c_sl * c_sl)))) * ((vel_r[i] * (1.0 - (c_sl * c_sl))) +
        (c_sl * sqrt((1.0 - v_sq_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_r * (c_sl * c_sl))) - (vel_r[i] * vel_r[i]) * (1.0 - (c_sl * c_sl)))))) - shift_r[i];
      
      slow_acoustic_eigs_r[i] = (lapse_r / (1.0 - (v_sq_r * (c_sl * c_sl)))) * ((vel_r[i] * (1.0 - (c_sl * c_sl))) -
        (c_sl * sqrt((1.0 - v_sq_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_r * (c_sl * c_sl))) - (vel_r[i] * vel_r[i]) * (1.0 - (c_sl * c_sl)))))) - shift_r[i];
    }

    double max_eig_r = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_r[i]) > max_eig_r) {
        max_eig_r = fabs(material_eigs_r[i]);
      }
      if (fabs(fast_acoustic_eigs_r[i]) > max_eig_r) {
        max_eig_r = fabs(fast_acoustic_eigs_r[i]);
      }
      if (fabs(slow_acoustic_eigs_r[i]) > max_eig_r) {
        max_eig_r = fabs(slow_acoustic_eigs_r[i]);
      }
    }

    double max_eig_avg = 0.5 * (max_eig_l + max_eig_r);
    
    sl = (vx_avg - max_eig_avg) / (1.0 - (vx_avg * max_eig_avg));
    sr = (vx_avg + max_eig_avg) / (1.0 + (vx_avg * max_eig_avg));
  }
  else {
    sl = (vx_avg - cs_avg) / (1.0 - (vx_avg * cs_avg));
    sr = (vx_avg + cs_avg) / (1.0 + (vx_avg * cs_avg));
  }

  double fl[71], fr[71];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  double qm[71];
  for (int i = 0; i < 71; i++) {
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl[i] - fr[i])) / (sr - sl);
  }

  double *w0 = &waves[0], *w1 = &waves[71];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 71; i++) {
      w0[i] = qm[i] - ql[i];
      w1[i] = qr[i] - qm[i];
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[71];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 71; i++) {
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
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double fr[71], fl[71];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[27] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[27] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 71; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
  }
  else {
    for (int m = 0; m < 71; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);

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
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gkyl_gr_euler_max_abs_speed(gas_gamma, q);
}

static inline void
gr_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 5; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double v[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, qin, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double lapse = v[5];
  double shift_x = v[6];
  double shift_y = v[7];
  double shift_z = v[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[9]; spatial_metric[0][1] = v[10]; spatial_metric[0][2] = v[11];
  spatial_metric[1][0] = v[12]; spatial_metric[1][1] = v[13]; spatial_metric[1][2] = v[14];
  spatial_metric[2][0] = v[15]; spatial_metric[2][1] = v[16]; spatial_metric[2][2] = v[17];

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = v[18]; extrinsic_curvature[0][1] = v[19]; extrinsic_curvature[0][2] = v[20];
  extrinsic_curvature[1][0] = v[21]; extrinsic_curvature[1][1] = v[22]; extrinsic_curvature[1][2] = v[23];
  extrinsic_curvature[2][0] = v[24]; extrinsic_curvature[2][1] = v[25]; extrinsic_curvature[2][2] = v[26];

  double lapse_der[3];
  lapse_der[0] = v[28];
  lapse_der[1] = v[29];
  lapse_der[2] = v[30];

  double shift_der[3][3];
  shift_der[0][0] = v[31]; shift_der[0][1] = v[32]; shift_der[0][2] = v[33];
  shift_der[1][0] = v[34]; shift_der[1][1] = v[35]; shift_der[1][2] = v[36];
  shift_der[2][0] = v[37]; shift_der[2][1] = v[38]; shift_der[2][2] = v[39];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = v[40]; spatial_metric_der[0][0][1] = v[41]; spatial_metric_der[0][0][2] = v[42];
  spatial_metric_der[0][1][0] = v[43]; spatial_metric_der[0][1][1] = v[44]; spatial_metric_der[0][1][2] = v[45];
  spatial_metric_der[0][2][0] = v[46]; spatial_metric_der[0][2][1] = v[47]; spatial_metric_der[0][2][2] = v[48];

  spatial_metric_der[1][0][0] = v[49]; spatial_metric_der[1][0][1] = v[50]; spatial_metric_der[1][0][2] = v[51];
  spatial_metric_der[1][1][0] = v[52]; spatial_metric_der[1][1][1] = v[53]; spatial_metric_der[1][1][2] = v[54];
  spatial_metric_der[1][2][0] = v[55]; spatial_metric_der[1][2][1] = v[56]; spatial_metric_der[1][2][2] = v[57];

  spatial_metric_der[0][0][0] = v[58]; spatial_metric_der[0][0][1] = v[59]; spatial_metric_der[0][0][2] = v[60];
  spatial_metric_der[0][1][0] = v[61]; spatial_metric_der[0][1][1] = v[62]; spatial_metric_der[0][1][2] = v[63];
  spatial_metric_der[0][2][0] = v[64]; spatial_metric_der[0][2][1] = v[65]; spatial_metric_der[0][2][2] = v[66];

  double **stress_energy = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    stress_energy[i] = gkyl_malloc(sizeof(double[4]));
  }

  gkyl_gr_euler_stress_energy_tensor(gas_gamma, qin, &stress_energy);

  bool in_excision_region = false;
  if (v[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

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
      W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
    
    double mom[3];
    mom[0] = (rho * h) * (W * W) * vx;
    mom[1] = (rho * h) * (W * W) * vy;
    mom[2] = (rho * h) * (W * W) * vz;

    // Energy density source.
    sout[4] = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sout[4] += stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j];
        sout[4] += 2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j];
        sout[4] += stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j];
      }

      sout[4] += stress_energy[0][0] * shift[i] * lapse_der[i];
      sout[4] -= stress_energy[0][i + 1] * lapse_der[i];
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      sout[1 + j] = -stress_energy[0][0] * lapse * lapse_der[j];

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          sout[1 + j] += 0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l];
          sout[1 + j] += 0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l];
        }

        sout[1 + j] += (mom[k] / lapse) * shift_der[j][k];

        for (int i = 0; i < 3; i++) {
          sout[1 + j] += stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k];
        }
      }
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      sout[i] = 0.0;
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(stress_energy[i]);
  }
  gkyl_free(stress_energy);
}

void
gkyl_gr_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_euler *gr_euler = container_of(base->on_dev, struct wv_gr_euler, eqn);
    gkyl_cu_free(gr_euler);
  }

  struct wv_gr_euler *gr_euler = container_of(base, struct wv_gr_euler, eqn);
  gkyl_free(gr_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_new(double gas_gamma, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_EULER_RP_HLL,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_inew(const struct gkyl_wv_gr_euler_inp* inp)
{
  struct wv_gr_euler *gr_euler = gkyl_malloc(sizeof(struct wv_gr_euler));

  gr_euler->eqn.type = GKYL_EQN_GR_EULER;
  gr_euler->eqn.num_equations = 71;
  gr_euler->eqn.num_diag = 5;

  gr_euler->gas_gamma = inp->gas_gamma;
  gr_euler->spacetime_gauge = inp->spacetime_gauge;
  gr_euler->reinit_freq = inp->reinit_freq;
  gr_euler->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_EULER_RP_LAX) {
    gr_euler->eqn.num_waves = 2;
    gr_euler->eqn.waves_func = wave_lax_l;
    gr_euler->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_EULER_RP_ROE) {
    gr_euler->eqn.num_waves = 3;
    gr_euler->eqn.waves_func = wave_roe_l;
    gr_euler->eqn.qfluct_func = qfluct_roe_l;
  }
  else if (inp->rp_type == WV_GR_EULER_RP_HLL) {
    gr_euler->eqn.num_waves = 2;
    gr_euler->eqn.waves_func = wave_hll_l;
    gr_euler->eqn.qfluct_func = qfluct_hll_l;
  }

  gr_euler->eqn.flux_jump = flux_jump;
  gr_euler->eqn.check_inv_func = check_inv;
  gr_euler->eqn.max_speed_func = max_speed;
  gr_euler->eqn.rotate_to_local_func = rot_to_local;
  gr_euler->eqn.rotate_to_global_func = rot_to_global;

  gr_euler->eqn.wall_bc_func = gr_euler_wall;
  gr_euler->eqn.no_slip_bc_func = gr_euler_no_slip;

  gr_euler->eqn.cons_to_riem = cons_to_riem;
  gr_euler->eqn.riem_to_cons = riem_to_cons;

  gr_euler->eqn.cons_to_diag = gr_euler_cons_to_diag;

  gr_euler->eqn.source_func = gr_euler_source;

  gr_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_euler->eqn.flags);
  gr_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_free);
  gr_euler->eqn.on_dev = &gr_euler->eqn; // On the CPU, the equation object points to itself.

  gr_euler->eqn.embed_geo = NULL;

  return &gr_euler->eqn;
}

double
gkyl_wv_gr_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gas_gamma;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_euler_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_euler->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_euler_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  int reinit_freq = gr_euler->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_euler_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_euler->spacetime;

  return spacetime;
}
