#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_twofluid_tetrad.h>
#include <gkyl_wv_gr_twofluid_tetrad_priv.h>

void
gkyl_gr_twofluid_tetrad_flux(double gas_gamma_elc, double gas_gamma_ion, double light_speed, double e_fact, double b_fact, const double q[84], double flux[84])
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_elc = v[0];
  double vx_elc = v[1];
  double vy_elc = v[2];
  double vz_elc = v[3];
  double p_elc = v[4];

  double rho_ion = v[5];
  double vx_ion = v[6];
  double vy_ion = v[7];
  double vz_ion = v[8];
  double p_ion = v[9];

  double Ex = v[10], Ey = v[11], Ez = v[12];
  double Bx = v[13], By = v[14], Bz = v[15];

  double phi = v[16];
  double psi = v[17];
  
  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double W_elc = 1.0 / (sqrt(1.0 - ((vx_elc * vx_elc) + (vy_elc * vy_elc) + (vz_elc * vz_elc))));
    if ((vx_elc * vx_elc) + (vy_elc * vy_elc) + (vz_elc * vz_elc) > 1.0 - pow(10.0, -8.0)) {
      W_elc = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_elc = 1.0 + ((p_elc / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));

    flux[0] = rho_elc * W_elc * vx_elc;
    flux[1] = (rho_elc * h_elc * (W_elc * W_elc) * (vx_elc * vx_elc)) + p_elc;
    flux[2] = rho_elc * h_elc * (W_elc * W_elc) * (vy_elc * vx_elc);
    flux[3] = rho_elc * h_elc * (W_elc * W_elc) * (vz_elc * vx_elc);
    flux[4] = ((rho_elc * h_elc * (W_elc * W_elc)) -  (rho_elc * W_elc)) * vx_elc;

    double W_ion = 1.0 / (sqrt(1.0 - ((vx_ion * vx_ion) + (vy_ion * vy_ion) + (vz_ion * vz_ion))));
    if ((vx_ion * vx_ion) + (vy_ion * vy_ion) + (vz_ion * vz_ion) > 1.0 - pow(10.0, -8.0)) {
      W_ion = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_ion = 1.0 + ((p_ion / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

    flux[5] = rho_ion * W_ion * vx_ion;
    flux[6] = (rho_ion * h_ion * (W_ion * W_ion) * (vx_ion * vx_ion)) + p_ion;
    flux[7] = rho_ion * h_ion * (W_ion * W_ion) * (vy_ion * vx_ion);
    flux[8] = rho_ion * h_ion * (W_ion * W_ion) * (vz_ion * vx_ion);
    flux[9] = ((rho_ion * h_ion * (W_ion * W_ion)) -  (rho_ion * W_ion)) * vx_ion;

    flux[10] = e_fact * (light_speed * light_speed) * phi;
    flux[11] = (light_speed * light_speed) * Bz;
    flux[12] = -(light_speed * light_speed) * By;
    flux[13] = b_fact * psi;
    flux[14] = -Ez;
    flux[15] = Ey;
    flux[16] = e_fact * Ex;
    flux[17] = b_fact * (light_speed * light_speed) * Bx;

    for (int i = 18; i < 84; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 84; i++) {
      flux[i] = 0.0;
    }
  }
}

void
gkyl_gr_twofluid_tetrad_flux_correction(double gas_gamma_elc, double gas_gamma_ion, double light_speed, double e_fact, double b_fact, const double q[84], const double flux_sr[84], double flux_gr[84])
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_elc = v[0];
  double vx_elc = v[1];
  double vy_elc = v[2];
  double vz_elc = v[3];
  double p_elc = v[4];

  double rho_ion = v[5];
  double vx_ion = v[6];
  double vy_ion = v[7];
  double vz_ion = v[8];
  double p_ion = v[9];

  double Ex = v[10];
  double Bx = v[13];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel_elc[3];
    double v_sq_elc = 0.0;
    vel_elc[0] = vx_elc; vel_elc[1] = vy_elc; vel_elc[2] = vz_elc;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_elc += spatial_metric[i][j] * vel_elc[i] * vel_elc[j];
      }
    }

    double W_flat_elc = 1.0 / (sqrt(1.0 - ((vx_elc * vx_elc) + (vy_elc * vy_elc) + (vz_elc * vz_elc))));
    if ((vx_elc * vx_elc) + (vy_elc * vy_elc) + (vz_elc * vz_elc) > 1.0 - pow(10.0, -8.0)) {
      W_flat_elc = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double W_curved_elc = 1.0 / (sqrt(1.0 - v_sq_elc));
    if (v_sq_elc > 1.0 - pow(10.0, -8.0)) {
      W_curved_elc = 1.0 / sqrt(pow(10.0, -8.0));
    }

    if (fabs(vx_elc) < pow(10.0, -8.0)) {
      if (vx_elc > 0.0) {
        vx_elc = pow(10.0, -8.0);
      }
      else {
        vx_elc = -pow(10.0, -8.0);
      }
    }

    flux_gr[0] = (lapse * sqrt(spatial_det)) * ((flux_sr[0] * (vx_elc - (shift_x / lapse)) * W_curved_elc) / (vx_elc * W_flat_elc));
    flux_gr[1] = (lapse * sqrt(spatial_det)) * ((((flux_sr[1] - p_elc) * (vx_elc - (shift_x / lapse)) * (W_curved_elc * W_curved_elc)) / (vx_elc * (W_flat_elc * W_flat_elc))) + p_elc);
    flux_gr[2] = (lapse * sqrt(spatial_det)) * ((flux_sr[2] * (vx_elc - (shift_x / lapse)) * (W_curved_elc * W_curved_elc)) / (vx_elc * (W_flat_elc * W_flat_elc)));
    flux_gr[3] = (lapse * sqrt(spatial_det)) * ((flux_sr[3] * (vx_elc - (shift_x / lapse)) * (W_curved_elc * W_curved_elc)) / (vx_elc * (W_flat_elc * W_flat_elc)));
    flux_gr[4] = (lapse * sqrt(spatial_det)) * (((((flux_sr[4] + (rho_elc * vx_elc * W_flat_elc)) * (W_curved_elc * W_curved_elc)) / (vx_elc * (W_flat_elc * W_flat_elc))) - p_elc -
      (rho_elc * W_curved_elc)) * (vx_elc - (shift_x / lapse)) + (p_elc * vx_elc));
    
    double vel_ion[3];
    double v_sq_ion = 0.0;
    vel_ion[0] = vx_ion; vel_ion[1] = vy_ion; vel_ion[2] = vz_ion;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_ion += spatial_metric[i][j] * vel_ion[i] * vel_ion[j];
      }
    }

    double W_flat_ion = 1.0 / (sqrt(1.0 - ((vx_ion * vx_ion) + (vy_ion * vy_ion) + (vz_ion * vz_ion))));
    if ((vx_ion * vx_ion) + (vy_ion * vy_ion) + (vz_ion * vz_ion) > 1.0 - pow(10.0, -8.0)) {
      W_flat_ion = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double W_curved_ion = 1.0 / (sqrt(1.0 - v_sq_ion));
    if (v_sq_ion > 1.0 - pow(10.0, -8.0)) {
      W_curved_ion = 1.0 / sqrt(pow(10.0, -8.0));
    }

    if (fabs(vx_ion) < pow(10.0, -8.0)) {
      if (vx_ion > 0.0) {
        vx_ion = pow(10.0, -8.0);
      }
      else {
        vx_ion = -pow(10.0, -8.0);
      }
    }

    flux_gr[5] = (lapse * sqrt(spatial_det)) * ((flux_sr[5] * (vx_ion - (shift_x / lapse)) * W_curved_ion) / (vx_ion * W_flat_ion));
    flux_gr[6] = (lapse * sqrt(spatial_det)) * ((((flux_sr[6] - p_ion) * (vx_ion - (shift_x / lapse)) * (W_curved_ion * W_curved_ion)) / (vx_ion * (W_flat_ion * W_flat_ion))) + p_ion);
    flux_gr[7] = (lapse * sqrt(spatial_det)) * ((flux_sr[7] * (vx_ion - (shift_x / lapse)) * (W_curved_ion * W_curved_ion)) / (vx_ion * (W_flat_ion * W_flat_ion)));
    flux_gr[8] = (lapse * sqrt(spatial_det)) * ((flux_sr[8] * (vx_ion - (shift_x / lapse)) * (W_curved_ion * W_curved_ion)) / (vx_ion * (W_flat_ion * W_flat_ion)));
    flux_gr[9] = (lapse * sqrt(spatial_det)) * (((((flux_sr[9] + (rho_ion * vx_ion * W_flat_ion)) * (W_curved_ion * W_curved_ion)) / (vx_ion * (W_flat_ion * W_flat_ion))) - p_ion -
      (rho_ion * W_curved_ion)) * (vx_ion - (shift_x / lapse)) + (p_ion * vx_ion));

    flux_gr[10] = flux_sr[10];
    flux_gr[11] = (lapse * flux_sr[11]) - (light_speed * light_speed) * ((shift_x * flux_sr[15]) - (shift_y * Ex));
    flux_gr[12] = (lapse * flux_sr[12]) + (light_speed * light_speed) * ((shift_x * flux_sr[14]) + (shift_z * Ex));
    flux_gr[13] = flux_sr[13];
    flux_gr[14] = (lapse * flux_sr[14]) + ((shift_x * (flux_sr[12] / (light_speed * light_speed))) + (shift_y * Bx));
    flux_gr[15] = (lapse * flux_sr[15]) - ((shift_x * (flux_sr[11] / (light_speed * light_speed))) - (shift_z * Bx));
    flux_gr[16] = flux_sr[16];
    flux_gr[17] = flux_sr[17];

    for (int i = 18; i < 84; i++) {
      flux_gr[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 84; i++) {
      flux_gr[i] = 0.0;
    }
  }
}

void
gkyl_gr_twofluid_tetrad_prim_vars(double gas_gamma_elc, double gas_gamma_ion, const double q[84], double v[84])
{
  double Dx = q[10], Dy = q[11], Dz = q[12];
  double Bx = q[13], By = q[14], Bz = q[15];

  double phi = q[16];
  double psi = q[17];

  double lapse = q[18];
  double shift_x = q[19];
  double shift_y = q[20];
  double shift_z = q[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = q[22]; spatial_metric[0][1] = q[23]; spatial_metric[0][2] = q[24];
  spatial_metric[1][0] = q[25]; spatial_metric[1][1] = q[26]; spatial_metric[1][2] = q[27];
  spatial_metric[2][0] = q[28]; spatial_metric[2][1] = q[29]; spatial_metric[2][2] = q[30];
  
  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[31]; extrinsic_curvature[0][1] = q[32]; extrinsic_curvature[0][2] = q[33];
  extrinsic_curvature[1][0] = q[34]; extrinsic_curvature[1][1] = q[35]; extrinsic_curvature[1][2] = q[36];
  extrinsic_curvature[2][0] = q[37]; extrinsic_curvature[2][1] = q[38]; extrinsic_curvature[2][2] = q[39];

  double lapse_der[3];
  lapse_der[0] = q[41];
  lapse_der[1] = q[42];
  lapse_der[2] = q[43];

  double shift_der[3][3];
  shift_der[0][0] = q[44]; shift_der[0][1] = q[45]; shift_der[0][2] = q[46];
  shift_der[1][0] = q[47]; shift_der[1][1] = q[48]; shift_der[1][2] = q[49];
  shift_der[2][0] = q[50]; shift_der[2][1] = q[51]; shift_der[2][2] = q[52];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[53]; spatial_metric_der[0][0][1] = q[54]; spatial_metric_der[0][0][2] = q[55];
  spatial_metric_der[0][1][0] = q[56]; spatial_metric_der[0][1][1] = q[57]; spatial_metric_der[0][1][2] = q[58];
  spatial_metric_der[0][2][0] = q[59]; spatial_metric_der[0][2][1] = q[60]; spatial_metric_der[0][2][2] = q[61];

  spatial_metric_der[1][0][0] = q[62]; spatial_metric_der[1][0][1] = q[63]; spatial_metric_der[1][0][2] = q[64];
  spatial_metric_der[1][1][0] = q[65]; spatial_metric_der[1][1][1] = q[66]; spatial_metric_der[1][1][2] = q[67];
  spatial_metric_der[1][2][0] = q[68]; spatial_metric_der[1][2][1] = q[69]; spatial_metric_der[1][2][2] = q[70];

  spatial_metric_der[0][0][0] = q[71]; spatial_metric_der[0][0][1] = q[72]; spatial_metric_der[0][0][2] = q[73];
  spatial_metric_der[0][1][0] = q[74]; spatial_metric_der[0][1][1] = q[75]; spatial_metric_der[0][1][2] = q[76];
  spatial_metric_der[0][2][0] = q[77]; spatial_metric_der[0][2][1] = q[78]; spatial_metric_der[0][2][2] = q[79];

  double evol_param = q[80];
  double x = q[81];
  double y = q[82];
  double z = q[83];

  bool in_excision_region = false;
  if (q[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D_elc = q[0] / sqrt(spatial_det);
    double momx_elc = q[1] / sqrt(spatial_det);
    double momy_elc = q[2] / sqrt(spatial_det);
    double momz_elc = q[3] / sqrt(spatial_det);
    double Etot_elc = q[4] / sqrt(spatial_det);

    double C_elc = D_elc / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    double C0_elc = (D_elc + Etot_elc) / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    if (((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)) < pow(10.0, -8.0)) {
      C_elc = D_elc / sqrt(pow(10.0, -8.0));
      C0_elc = (D_elc + Etot_elc) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_elc = -1.0 / (gas_gamma_elc * gas_gamma_elc);
    double alpha1_elc = -2.0 * C_elc * ((gas_gamma_elc - 1.0) / (gas_gamma_elc * gas_gamma_elc));
    double alpha2_elc = ((gas_gamma_elc - 2.0) / gas_gamma_elc) * ((C0_elc * C0_elc) - 1.0) + 1.0 - (C_elc * C_elc) * ((gas_gamma_elc - 1.0) / gas_gamma_elc) *
      ((gas_gamma_elc - 1.0) / gas_gamma_elc);
    double alpha4_elc = (C0_elc * C0_elc) - 1.0;
    double eta_elc = 2.0 * C_elc *((gas_gamma_elc - 1.0) / gas_gamma_elc);

    double guess_elc = 1.0;
    int iter_elc = 0;

    while (iter_elc < 100) {
      double poly_elc = (alpha4_elc * (guess_elc * guess_elc * guess_elc) * (guess_elc - eta_elc)) + (alpha2_elc * (guess_elc * guess_elc)) +
        (alpha1_elc * guess_elc) + alpha0_elc;
      double poly_der_elc = alpha1_elc + (2.0 * alpha2_elc * guess_elc) + (4.0 * alpha4_elc * (guess_elc * guess_elc * guess_elc)) -
        (3.0 * eta_elc * alpha4_elc * (guess_elc * guess_elc));

      double guess_new_elc = guess_elc - (poly_elc / poly_der_elc);

      if (fabs(guess_elc - guess_new_elc) < pow(10.0, -8.0)) {
        iter_elc = 100;
      }
      else {
        iter_elc += 1;
        guess_elc = guess_new_elc;
      }
    }

    double W_elc = 0.5 * C0_elc * guess_elc * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((1.0 - (C_elc * guess_elc)) /
      ((C0_elc * C0_elc) * (guess_elc * guess_elc))))));
    double h_elc = 1.0 / (C_elc * guess_elc);

    v[0] = D_elc / W_elc; 
    v[1] = momx_elc / (v[0] * h_elc * (W_elc * W_elc));
    v[2] = momy_elc / (v[0] * h_elc * (W_elc * W_elc));
    v[3] = momz_elc / (v[0] * h_elc * (W_elc * W_elc));
    v[4] = (v[0] * h_elc * (W_elc * W_elc)) - D_elc - Etot_elc;

    if (v[0] < pow(10.0, -8.0)) {
      v[0] = pow(10.0, -8.0);
    }
    if (v[4] < pow(10.0, -8.0)) {
      v[4] = pow(10.0, -8.0);
    }

    double D_ion = q[5] / sqrt(spatial_det);
    double momx_ion = q[6] / sqrt(spatial_det);
    double momy_ion = q[7] / sqrt(spatial_det);
    double momz_ion = q[8] / sqrt(spatial_det);
    double Etot_ion = q[9] / sqrt(spatial_det);

    double C_ion = D_ion / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    double C0_ion = (D_ion + Etot_ion) / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    if (((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)) < pow(10.0, -8.0)) {
      C_ion = D_ion / sqrt(pow(10.0, -8.0));
      C0_ion = (D_ion + Etot_ion) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_ion = -1.0 / (gas_gamma_ion * gas_gamma_ion);
    double alpha1_ion = -2.0 * C_ion * ((gas_gamma_ion - 1.0) / (gas_gamma_ion * gas_gamma_ion));
    double alpha2_ion = ((gas_gamma_ion - 2.0) / gas_gamma_ion) * ((C0_ion * C0_ion) - 1.0) + 1.0 - (C_ion * C_ion) * ((gas_gamma_ion - 1.0) / gas_gamma_ion) *
      ((gas_gamma_ion - 1.0) / gas_gamma_ion);
    double alpha4_ion = (C0_ion * C0_ion) - 1.0;
    double eta_ion = 2.0 * C_ion *((gas_gamma_ion - 1.0) / gas_gamma_ion);

    double guess_ion = 1.0;
    int iter_ion = 0;

    while (iter_ion < 100) {
      double poly_ion = (alpha4_ion * (guess_ion * guess_ion * guess_ion) * (guess_ion - eta_ion)) + (alpha2_ion * (guess_ion * guess_ion)) +
        (alpha1_ion * guess_ion) + alpha0_ion;
      double poly_der_ion = alpha1_ion + (2.0 * alpha2_ion * guess_ion) + (4.0 * alpha4_ion * (guess_ion * guess_ion * guess_ion)) -
        (3.0 * eta_ion * alpha4_ion * (guess_ion * guess_ion));

      double guess_new_ion = guess_ion - (poly_ion / poly_der_ion);

      if (fabs(guess_ion - guess_new_ion) < pow(10.0, -8.0)) {
        iter_ion = 100;
      }
      else {
        iter_ion += 1;
        guess_ion = guess_new_ion;
      }
    }

    double W_ion = 0.5 * C0_ion * guess_ion * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((1.0 - (C_ion * guess_ion)) /
      ((C0_ion * C0_ion) * (guess_ion * guess_ion))))));
    double h_ion = 1.0 / (C_ion * guess_ion);

    v[5] = D_ion / W_ion; 
    v[6] = momx_ion / (v[5] * h_ion * (W_ion * W_ion));
    v[7] = momy_ion / (v[5] * h_ion * (W_ion * W_ion));
    v[8] = momz_ion / (v[5] * h_ion * (W_ion * W_ion));
    v[9] = (v[5] * h_ion * (W_ion * W_ion)) - D_ion - Etot_ion;

    if (v[5] < pow(10.0, -8.0)) {
      v[5] = pow(10.0, -8.0);
    }
    if (v[9] < pow(10.0, -8.0)) {
      v[9] = pow(10.0, -8.0);
    }

    v[10] = Dx; v[11] = Dy; v[12] = Dz;
    v[13] = Bx; v[14] = By; v[15] = Bz;

    v[16] = phi;
    v[17] = psi;

    v[18] = lapse;
    v[19] = shift_x;
    v[20] = shift_y;
    v[21] = shift_z;

    v[22] = spatial_metric[0][0]; v[23] = spatial_metric[0][1]; v[24] = spatial_metric[0][2];
    v[25] = spatial_metric[1][0]; v[26] = spatial_metric[1][1]; v[27] = spatial_metric[1][2];
    v[28] = spatial_metric[2][0]; v[29] = spatial_metric[2][1]; v[30] = spatial_metric[2][2];

    v[31] = extrinsic_curvature[0][0]; v[32] = extrinsic_curvature[0][1]; v[33] = extrinsic_curvature[0][2];
    v[34] = extrinsic_curvature[1][0]; v[35] = extrinsic_curvature[1][1]; v[36] = extrinsic_curvature[1][2];
    v[37] = extrinsic_curvature[2][0]; v[38] = extrinsic_curvature[2][1]; v[39] = extrinsic_curvature[2][2];

    v[40] = 1.0;

    v[41] = lapse_der[0];
    v[42] = lapse_der[1];
    v[43] = lapse_der[2];

    v[44] = shift_der[0][0]; v[45] = shift_der[0][1]; v[46] = shift_der[0][2];
    v[47] = shift_der[1][0]; v[48] = shift_der[1][1]; v[49] = shift_der[1][2];
    v[50] = shift_der[2][0]; v[51] = shift_der[2][1]; v[52] = shift_der[2][2];

    v[53] = spatial_metric_der[0][0][0]; v[54] = spatial_metric_der[0][0][1]; v[55] = spatial_metric_der[0][0][2];
    v[56] = spatial_metric_der[0][1][0]; v[57] = spatial_metric_der[0][1][1]; v[58] = spatial_metric_der[0][1][2];
    v[59] = spatial_metric_der[0][2][0]; v[60] = spatial_metric_der[0][2][1]; v[61] = spatial_metric_der[0][2][2];

    v[62] = spatial_metric_der[1][0][0]; v[63] = spatial_metric_der[1][0][1]; v[64] = spatial_metric_der[1][0][2];
    v[65] = spatial_metric_der[1][1][0]; v[66] = spatial_metric_der[1][1][1]; v[67] = spatial_metric_der[1][1][2];
    v[68] = spatial_metric_der[1][2][0]; v[69] = spatial_metric_der[1][2][1]; v[70] = spatial_metric_der[1][2][2];

    v[71] = spatial_metric_der[2][0][0]; v[72] = spatial_metric_der[2][0][1]; v[73] = spatial_metric_der[2][0][2];
    v[74] = spatial_metric_der[2][1][0]; v[75] = spatial_metric_der[2][1][1]; v[76] = spatial_metric_der[2][1][2];
    v[77] = spatial_metric_der[2][2][0]; v[78] = spatial_metric_der[2][2][1]; v[79] = spatial_metric_der[2][2][2];

    v[80] = evol_param;
    v[81] = x;
    v[82] = y;
    v[83] = z;
  }
  else {
    for (int i = 0; i < 84; i++) {
      v[i] = 0.0;
    }
    
    v[40] = -1.0;
  }
}

void 
gkyl_gr_twofluid_tetrad_inv_spatial_metric(const double q[84], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[22]; spatial_metric[0][1] = q[23]; spatial_metric[0][2] = q[24];
  spatial_metric[1][0] = q[25]; spatial_metric[1][1] = q[26]; spatial_metric[1][2] = q[27];
  spatial_metric[2][0] = q[28]; spatial_metric[2][1] = q[29]; spatial_metric[2][2] = q[30];

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
gkyl_gr_twofluid_tetrad_stress_energy_tensor_elc(double gas_gamma_elc, double gas_gamma_ion, const double q[84], double ***stress_energy_elc)
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_elc = v[0];
  double vx_elc = v[1];
  double vy_elc = v[2];
  double vz_elc = v[3];
  double p_elc = v[4];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_tetrad_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel_elc[3];
    double v_sq_elc = 0.0;
    vel_elc[0] = vx_elc; vel_elc[1] = vy_elc; vel_elc[2] = vz_elc;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_elc += spatial_metric[i][j] * vel_elc[i] * vel_elc[j];
      }
    }

    double W_elc = 1.0 / sqrt(1.0 - v_sq_elc);
    if (v_sq_elc > 1.0 - pow(10.0, -8.0)) {
      W_elc = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_elc = 1.0 + ((p_elc / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));

    double spacetime_vel_elc[4];
    spacetime_vel_elc[0] = W_elc / lapse;
    spacetime_vel_elc[1] = (W_elc * vx_elc) - (shift_x * (W_elc / lapse));
    spacetime_vel_elc[2] = (W_elc * vy_elc) - (shift_y * (W_elc / lapse));
    spacetime_vel_elc[3] = (W_elc * vz_elc) - (shift_z * (W_elc / lapse));

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
        (*stress_energy_elc)[i][j] = (rho_elc * h_elc * spacetime_vel_elc[i] * spacetime_vel_elc[j]) + (p_elc * inv_spacetime_metric[i][j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy_elc)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

void
gkyl_gr_twofluid_tetrad_stress_energy_tensor_ion(double gas_gamma_elc, double gas_gamma_ion, const double q[84], double ***stress_energy_ion)
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_ion = v[5];
  double vx_ion = v[6];
  double vy_ion = v[7];
  double vz_ion = v[8];
  double p_ion = v[9];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_tetrad_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel_ion[3];
    double v_sq_ion = 0.0;
    vel_ion[0] = vx_ion; vel_ion[1] = vy_ion; vel_ion[2] = vz_ion;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_ion += spatial_metric[i][j] * vel_ion[i] * vel_ion[j];
      }
    }

    double W_ion = 1.0 / sqrt(1.0 - v_sq_ion);
    if (v_sq_ion > 1.0 - pow(10.0, -8.0)) {
      W_ion = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_ion = 1.0 + ((p_ion / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

    double spacetime_vel_ion[4];
    spacetime_vel_ion[0] = W_ion / lapse;
    spacetime_vel_ion[1] = (W_ion * vx_ion) - (shift_x * (W_ion / lapse));
    spacetime_vel_ion[2] = (W_ion * vy_ion) - (shift_y * (W_ion / lapse));
    spacetime_vel_ion[3] = (W_ion * vz_ion) - (shift_z * (W_ion / lapse));

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
        (*stress_energy_ion)[i][j] = (rho_ion * h_ion * spacetime_vel_ion[i] * spacetime_vel_ion[j]) + (p_ion * inv_spacetime_metric[i][j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy_ion)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

static inline double
gkyl_gr_twofluid_tetrad_max_abs_speed(double gas_gamma_elc, double gas_gamma_ion, double light_speed, const double q[84])
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_elc = v[0];
  double vx_elc = v[1];
  double vy_elc = v[2];
  double vz_elc = v[3];
  double p_elc = v[4];

  double rho_ion = v[5];
  double vx_ion = v[6];
  double vy_ion = v[7];
  double vz_ion = v[8];
  double p_ion = v[9];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_tetrad_inv_spatial_metric(q, &inv_spatial_metric);

  double num_elc = (gas_gamma_elc * p_elc) / rho_elc;
  double den_elc = 1.0 + ((p_elc / rho_elc) * (gas_gamma_elc) / (gas_gamma_elc - 1.0));
  double c_s_elc = sqrt(num_elc / den_elc);

  double num_ion = (gas_gamma_ion * p_ion) / rho_ion;
  double den_ion = 1.0 + ((p_ion / rho_ion) * (gas_gamma_ion) / (gas_gamma_ion - 1.0));
  double c_s_ion = sqrt(num_ion / den_ion);

  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
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
      double vel_elc[3];
      double v_sq_elc = 0.0;
      vel_elc[0] = vx_elc; vel_elc[1] = vy_elc; vel_elc[2] = vz_elc;
      
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_elc += spatial_metric[i][j] * vel_elc[i] * vel_elc[j];
        }
      }

      double shift[3];
      shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

      double material_eigs_elc[3];
      double fast_acoustic_eigs_elc[3];
      double slow_acoustic_eigs_elc[3];

      for (int i = 0; i < 3; i++) {
        material_eigs_elc[i] = (lapse * vel_elc[i]) - shift[i];

        fast_acoustic_eigs_elc[i] = (lapse / (1.0 - (v_sq_elc * (c_s_elc * c_s_elc)))) * ((vel_elc[i] * (1.0 - (c_s_elc * c_s_elc))) +
          (c_s_elc * sqrt((1.0 - v_sq_elc) * (inv_spatial_metric[i][i] * (1.0 - (v_sq_elc * (c_s_elc * c_s_elc))) -
          (vel_elc[i] * vel_elc[i]) * (1.0 - (c_s_elc * c_s_elc)))))) - shift[i];
        
        slow_acoustic_eigs_elc[i] = (lapse / (1.0 - (v_sq_elc * (c_s_elc * c_s_elc)))) * ((vel_elc[i] * (1.0 - (c_s_elc * c_s_elc))) -
          (c_s_elc * sqrt((1.0 - v_sq_elc) * (inv_spatial_metric[i][i] * (1.0 - (v_sq_elc * (c_s_elc * c_s_elc))) -
          (vel_elc[i] * vel_elc[i]) * (1.0 - (c_s_elc * c_s_elc)))))) - shift[i];
      }

      double vel_ion[3];
      double v_sq_ion = 0.0;
      vel_ion[0] = vx_ion; vel_ion[1] = vy_ion; vel_ion[2] = vz_ion;
      
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_ion += spatial_metric[i][j] * vel_ion[i] * vel_ion[j];
        }
      }

      double material_eigs_ion[3];
      double fast_acoustic_eigs_ion[3];
      double slow_acoustic_eigs_ion[3];

      for (int i = 0; i < 3; i++) {
        material_eigs_ion[i] = (lapse * vel_ion[i]) - shift[i];

        fast_acoustic_eigs_ion[i] = (lapse / (1.0 - (v_sq_ion * (c_s_ion * c_s_ion)))) * ((vel_ion[i] * (1.0 - (c_s_ion * c_s_ion))) +
          (c_s_ion * sqrt((1.0 - v_sq_ion) * (inv_spatial_metric[i][i] * (1.0 - (v_sq_ion * (c_s_ion * c_s_ion))) -
          (vel_ion[i] * vel_ion[i]) * (1.0 - (c_s_ion * c_s_ion)))))) - shift[i];
        
        slow_acoustic_eigs_ion[i] = (lapse / (1.0 - (v_sq_ion * (c_s_ion * c_s_ion)))) * ((vel_ion[i] * (1.0 - (c_s_ion * c_s_ion))) -
          (c_s_ion * sqrt((1.0 - v_sq_ion) * (inv_spatial_metric[i][i] * (1.0 - (v_sq_ion * (c_s_ion * c_s_ion))) -
          (vel_ion[i] * vel_ion[i]) * (1.0 - (c_s_ion * c_s_ion)))))) - shift[i];
      }

      double max_eig = 0.0;
      for (int i = 0; i < 3; i++) {
        if (fabs(material_eigs_elc[i]) > max_eig) {
          max_eig = fabs(material_eigs_elc[i]);
        }
        if (fabs(fast_acoustic_eigs_elc[i]) > max_eig) {
          max_eig = fabs(fast_acoustic_eigs_elc[i]);
        }
        if (fabs(slow_acoustic_eigs_elc[i]) > max_eig) {
          max_eig = fabs(slow_acoustic_eigs_elc[i]);
        }

        if (fabs(material_eigs_ion[i]) > max_eig) {
          max_eig = fabs(material_eigs_ion[i]);
        }
        if (fabs(fast_acoustic_eigs_ion[i]) > max_eig) {
          max_eig = fabs(fast_acoustic_eigs_ion[i]);
        }
        if (fabs(slow_acoustic_eigs_ion[i]) > max_eig) {
          max_eig = fabs(slow_acoustic_eigs_ion[i]);
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(inv_spatial_metric);

      double v_sq = 0.0;
      if (v_sq_elc > v_sq_ion) {
        v_sq = v_sq_elc;
      }
      else {
        v_sq = v_sq_ion;
      }

      if (fabs(v_sq) + max_eig > light_speed * sqrt(spatial_metric_det) * lapse) {
        return fabs(v_sq) + max_eig;
      }
      else {
        return light_speed * sqrt(spatial_metric_det) * lapse;
      }
    }
    else {
      double v_sq_elc = sqrt((vx_elc * vx_elc) + (vy_elc * vy_elc) + (vz_elc * vz_elc));
      double v_sq_ion = sqrt((vx_ion * vx_ion) + (vy_ion * vy_ion) + (vz_ion * vz_ion));

      double v_sq = 0.0;
      if (v_sq_elc > v_sq_ion) {
        v_sq = v_sq_elc;
      }
      else {
        v_sq = v_sq_ion;
      }

      double c_s = 0.0;
      if (c_s_elc > c_s_ion) {
        c_s = c_s_elc;
      }
      else {
        c_s = c_s_ion;
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(inv_spatial_metric);

      if (fabs(v_sq) + c_s > light_speed * sqrt(spatial_metric_det) * lapse) {
        return fabs(v_sq) + c_s;
      }
      else {
        return light_speed * sqrt(spatial_metric_det) * lapse;
      }
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
  for (int i = 0; i < 84; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 84; i++) {
    qout[i] = win[i];
  }
}

static void
gr_twofluid_tetrad_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 84; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
  ghost[6] = -ghost[6];
}

static void
gr_twofluid_tetrad_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];

  for (int i = 6; i < 9; i++) {
    ghost[i] = -skin[i];
  }

  ghost[5] = skin[5];
  ghost[9] = skin[9];

  for (int i = 10; i < 84; i++) {
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
  qlocal[9] = qglobal[9];

  qlocal[10] = (qglobal[10] * norm[0]) + (qglobal[11] * norm[1]) + (qglobal[12] * norm[2]);
  qlocal[11] = (qglobal[10] * tau1[0]) + (qglobal[11] * tau1[1]) + (qglobal[12] * tau1[2]);
  qlocal[12] = (qglobal[10] * tau2[0]) + (qglobal[11] * tau2[1]) + (qglobal[12] * tau2[2]);

  qlocal[13] = (qglobal[13] * norm[0]) + (qglobal[14] * norm[1]) + (qglobal[15] * norm[2]);
  qlocal[14] = (qglobal[13] * tau1[0]) + (qglobal[14] * tau1[1]) + (qglobal[15] * tau1[2]);
  qlocal[15] = (qglobal[13] * tau2[0]) + (qglobal[14] * tau2[1]) + (qglobal[15] * tau2[2]);

  qlocal[16] = qglobal[16];
  qlocal[17] = qglobal[17];

  qlocal[18] = qglobal[18];
  qlocal[19] = (qglobal[19] * norm[0]) + (qglobal[20] * norm[1]) + (qglobal[21] * norm[2]);
  qlocal[20] = (qglobal[19] * tau1[0]) + (qglobal[20] * tau1[1]) + (qglobal[21] * tau1[2]);
  qlocal[21] = (qglobal[19] * tau2[0]) + (qglobal[20] * tau2[1]) + (qglobal[21] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[22] * norm[0]) + (qglobal[23] * norm[1]) + (qglobal[24] * norm[2]);
  r1[1] = (qglobal[22] * tau1[0]) + (qglobal[23] * tau1[1]) + (qglobal[24] * tau1[2]);
  r1[2] = (qglobal[22] * tau2[0]) + (qglobal[23] * tau2[1]) + (qglobal[24] * tau2[2]);

  r2[0] = (qglobal[25] * norm[0]) + (qglobal[26] * norm[1]) + (qglobal[27] * norm[2]);
  r2[1] = (qglobal[25] * tau1[0]) + (qglobal[26] * tau1[1]) + (qglobal[27] * tau1[2]);
  r2[2] = (qglobal[25] * tau2[0]) + (qglobal[26] * tau2[1]) + (qglobal[27] * tau2[2]);

  r3[0] = (qglobal[28] * norm[0]) + (qglobal[29] * norm[1]) + (qglobal[30] * norm[2]);
  r3[1] = (qglobal[28] * tau1[0]) + (qglobal[29] * tau1[1]) + (qglobal[30] * tau1[2]);
  r3[2] = (qglobal[28] * tau2[0]) + (qglobal[29] * tau2[1]) + (qglobal[30] * tau2[2]);

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
  qlocal[22] = v1[0]; qlocal[23] = v1[1]; qlocal[24] = v1[2];
  qlocal[25] = v2[0]; qlocal[26] = v2[1]; qlocal[27] = v2[2];
  qlocal[28] = v3[0]; qlocal[29] = v3[1]; qlocal[30] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[31] * norm[0]) + (qglobal[32] * norm[1]) + (qglobal[33] * norm[2]);
  extr_r1[1] = (qglobal[31] * tau1[0]) + (qglobal[32] * tau1[1]) + (qglobal[33] * tau1[2]);
  extr_r1[2] = (qglobal[31] * tau2[0]) + (qglobal[32] * tau2[1]) + (qglobal[33] * tau2[2]);

  extr_r2[0] = (qglobal[34] * norm[0]) + (qglobal[35] * norm[1]) + (qglobal[36] * norm[2]);
  extr_r2[1] = (qglobal[34] * tau1[0]) + (qglobal[35] * tau1[1]) + (qglobal[36] * tau1[2]);
  extr_r2[2] = (qglobal[34] * tau2[0]) + (qglobal[35] * tau2[1]) + (qglobal[36] * tau2[2]);

  extr_r3[0] = (qglobal[37] * norm[0]) + (qglobal[38] * norm[1]) + (qglobal[39] * norm[2]);
  extr_r3[1] = (qglobal[37] * tau1[0]) + (qglobal[38] * tau1[1]) + (qglobal[39] * tau1[2]);
  extr_r3[2] = (qglobal[37] * tau2[0]) + (qglobal[38] * tau2[1]) + (qglobal[39] * tau2[2]);

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
  qlocal[31] = inv_v1[0]; qlocal[32] = inv_v1[1]; qlocal[33] = inv_v1[2];
  qlocal[34] = inv_v2[0]; qlocal[35] = inv_v2[1]; qlocal[36] = inv_v2[2];
  qlocal[37] = inv_v3[0]; qlocal[38] = inv_v3[1]; qlocal[39] = inv_v3[2];

  qlocal[40] = qglobal[40];

  qlocal[41] = (qglobal[41] * norm[0]) + (qglobal[42] * norm[1]) + (qglobal[43] * norm[2]);
  qlocal[42] = (qglobal[41] * tau1[0]) + (qglobal[42] * tau1[1]) + (qglobal[43] * tau1[2]);
  qlocal[43] = (qglobal[41] * tau2[0]) + (qglobal[42] * tau2[1]) + (qglobal[43] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qglobal[44] * norm[0]) + (qglobal[45] * norm[1]) + (qglobal[46] * norm[2]);
  shiftder_r1[1] = (qglobal[44] * tau1[0]) + (qglobal[45] * tau1[1]) + (qglobal[46] * tau1[2]);
  shiftder_r1[2] = (qglobal[44] * tau2[0]) + (qglobal[45] * tau2[1]) + (qglobal[46] * tau2[2]);

  shiftder_r2[0] = (qglobal[47] * norm[0]) + (qglobal[48] * norm[1]) + (qglobal[49] * norm[2]);
  shiftder_r2[1] = (qglobal[47] * tau1[0]) + (qglobal[48] * tau1[1]) + (qglobal[49] * tau1[2]);
  shiftder_r2[2] = (qglobal[47] * tau2[0]) + (qglobal[48] * tau2[1]) + (qglobal[49] * tau2[2]);

  shiftder_r3[0] = (qglobal[50] * norm[0]) + (qglobal[51] * norm[1]) + (qglobal[52] * norm[2]);
  shiftder_r3[1] = (qglobal[50] * tau1[0]) + (qglobal[51] * tau1[1]) + (qglobal[52] * tau1[2]);
  shiftder_r3[2] = (qglobal[50] * tau2[0]) + (qglobal[51] * tau2[1]) + (qglobal[52] * tau2[2]);

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
  qlocal[44] = shiftder_v1[0]; qlocal[45] = shiftder_v1[1]; qlocal[46] = shiftder_v1[2];
  qlocal[47] = shiftder_v2[0]; qlocal[48] = shiftder_v2[1]; qlocal[49] = shiftder_v2[2];
  qlocal[50] = shiftder_v3[0]; qlocal[51] = shiftder_v3[1]; qlocal[52] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qglobal[53] * norm[0]) + (qglobal[54] * norm[1]) + (qglobal[55] * norm[2]);
  r11[1] = (qglobal[53] * tau1[0]) + (qglobal[54] * tau1[1]) + (qglobal[55] * tau1[2]);
  r11[2] = (qglobal[53] * tau2[0]) + (qglobal[54] * tau2[1]) + (qglobal[55] * tau2[2]);

  r12[0] = (qglobal[56] * norm[0]) + (qglobal[57] * norm[1]) + (qglobal[58] * norm[2]);
  r12[1] = (qglobal[56] * tau1[0]) + (qglobal[57] * tau1[1]) + (qglobal[58] * tau1[2]);
  r12[2] = (qglobal[56] * tau2[0]) + (qglobal[57] * tau2[1]) + (qglobal[58] * tau2[2]);

  r13[0] = (qglobal[59] * norm[0]) + (qglobal[60] * norm[1]) + (qglobal[61] * norm[2]);
  r13[1] = (qglobal[59] * tau1[0]) + (qglobal[60] * tau1[1]) + (qglobal[61] * tau1[2]);
  r13[2] = (qglobal[59] * tau2[0]) + (qglobal[60] * tau2[1]) + (qglobal[61] * tau2[2]);

  r21[0] = (qglobal[62] * norm[0]) + (qglobal[63] * norm[1]) + (qglobal[64] * norm[2]);
  r21[1] = (qglobal[62] * tau1[0]) + (qglobal[63] * tau1[1]) + (qglobal[64] * tau1[2]);
  r21[2] = (qglobal[62] * tau2[0]) + (qglobal[63] * tau2[1]) + (qglobal[64] * tau2[2]);

  r22[0] = (qglobal[65] * norm[0]) + (qglobal[66] * norm[1]) + (qglobal[67] * norm[2]);
  r22[1] = (qglobal[65] * tau1[0]) + (qglobal[66] * tau1[1]) + (qglobal[67] * tau1[2]);
  r22[2] = (qglobal[65] * tau2[0]) + (qglobal[66] * tau2[1]) + (qglobal[67] * tau2[2]);

  r23[0] = (qglobal[68] * norm[0]) + (qglobal[69] * norm[1]) + (qglobal[70] * norm[2]);
  r23[1] = (qglobal[68] * tau1[0]) + (qglobal[69] * tau1[1]) + (qglobal[70] * tau1[2]);
  r23[2] = (qglobal[68] * tau2[0]) + (qglobal[69] * tau2[1]) + (qglobal[70] * tau2[2]);

  r31[0] = (qglobal[71] * norm[0]) + (qglobal[72] * norm[1]) + (qglobal[73] * norm[2]);
  r31[1] = (qglobal[71] * tau1[0]) + (qglobal[72] * tau1[1]) + (qglobal[73] * tau1[2]);
  r31[2] = (qglobal[71] * tau2[0]) + (qglobal[72] * tau2[1]) + (qglobal[73] * tau2[2]);

  r32[0] = (qglobal[74] * norm[0]) + (qglobal[75] * norm[1]) + (qglobal[76] * norm[2]);
  r32[1] = (qglobal[74] * tau1[0]) + (qglobal[75] * tau1[1]) + (qglobal[76] * tau1[2]);
  r32[2] = (qglobal[74] * tau2[0]) + (qglobal[75] * tau2[1]) + (qglobal[76] * tau2[2]);

  r33[0] = (qglobal[77] * norm[0]) + (qglobal[78] * norm[1]) + (qglobal[79] * norm[2]);
  r33[1] = (qglobal[77] * tau1[0]) + (qglobal[78] * tau1[1]) + (qglobal[79] * tau1[2]);
  r33[2] = (qglobal[77] * tau2[0]) + (qglobal[78] * tau2[1]) + (qglobal[79] * tau2[2]);

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
  qlocal[53] = (s11[0] * norm[0]) + (s21[0] * norm[1]) + (s31[0] * norm[2]);
  qlocal[54] = (s11[1] * norm[0]) + (s21[1] * norm[1]) + (s31[1] * norm[2]);
  qlocal[55] = (s11[2] * norm[0]) + (s21[2] * norm[1]) + (s31[2] * norm[2]);

  qlocal[56] = (s12[0] * norm[0]) + (s22[0] * norm[1]) + (s32[0] * norm[2]);
  qlocal[57] = (s12[1] * norm[0]) + (s22[1] * norm[1]) + (s32[1] * norm[2]);
  qlocal[58] = (s12[2] * norm[0]) + (s22[2] * norm[1]) + (s32[2] * norm[2]);

  qlocal[59] = (s13[0] * norm[0]) + (s23[0] * norm[1]) + (s33[0] * norm[2]);
  qlocal[60] = (s13[1] * norm[0]) + (s23[1] * norm[1]) + (s33[1] * norm[2]);
  qlocal[61] = (s13[2] * norm[0]) + (s23[2] * norm[1]) + (s33[2] * norm[2]);

  qlocal[62] = (s11[0] * tau1[0]) + (s21[0] * tau1[1]) + (s31[0] * tau1[2]);
  qlocal[63] = (s11[1] * tau1[0]) + (s21[1] * tau1[1]) + (s31[1] * tau1[2]);
  qlocal[64] = (s11[2] * tau1[0]) + (s21[2] * tau1[1]) + (s31[2] * tau1[2]);

  qlocal[65] = (s12[0] * tau1[0]) + (s22[0] * tau1[1]) + (s32[0] * tau1[2]);
  qlocal[66] = (s12[1] * tau1[0]) + (s22[1] * tau1[1]) + (s32[1] * tau1[2]);
  qlocal[67] = (s12[2] * tau1[0]) + (s22[2] * tau1[1]) + (s32[2] * tau1[2]);

  qlocal[68] = (s13[0] * tau1[0]) + (s23[0] * tau1[1]) + (s33[0] * tau1[2]);
  qlocal[69] = (s13[1] * tau1[0]) + (s23[1] * tau1[1]) + (s33[1] * tau1[2]);
  qlocal[70] = (s13[2] * tau1[0]) + (s23[2] * tau1[1]) + (s33[2] * tau1[2]);

  qlocal[71] = (s11[0] * tau2[0]) + (s21[0] * tau2[1]) + (s31[0] * tau2[2]);
  qlocal[72] = (s11[1] * tau2[0]) + (s21[1] * tau2[1]) + (s31[1] * tau2[2]);
  qlocal[73] = (s11[2] * tau2[0]) + (s21[2] * tau2[1]) + (s31[2] * tau2[2]);

  qlocal[74] = (s12[0] * tau2[0]) + (s22[0] * tau2[1]) + (s32[0] * tau2[2]);
  qlocal[75] = (s12[1] * tau2[0]) + (s22[1] * tau2[1]) + (s32[1] * tau2[2]);
  qlocal[76] = (s12[2] * tau2[0]) + (s22[2] * tau2[1]) + (s32[2] * tau2[2]);

  qlocal[77] = (s13[0] * tau2[0]) + (s23[0] * tau2[1]) + (s33[0] * tau2[2]);
  qlocal[78] = (s13[1] * tau2[0]) + (s23[1] * tau2[1]) + (s33[1] * tau2[2]);
  qlocal[79] = (s13[2] * tau2[0]) + (s23[2] * tau2[1]) + (s33[2] * tau2[2]);

  qlocal[80] = qglobal[80];
  qlocal[81] = (qglobal[81] * norm[0]) + (qglobal[82] * norm[1]) + (qglobal[83] * norm[2]);
  qlocal[82] = (qglobal[81] * tau1[0]) + (qglobal[82] * tau1[1]) + (qglobal[83] * tau1[2]);
  qlocal[83] = (qglobal[81] * tau2[0]) + (qglobal[82] * tau2[1]) + (qglobal[83] * tau2[2]);
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
  qglobal[9] = qlocal[9];

  qglobal[10] = (qlocal[10] * norm[0]) + (qlocal[11] * tau1[0]) + (qlocal[12] * tau2[0]);
  qglobal[11] = (qlocal[10] * norm[1]) + (qlocal[11] * tau1[1]) + (qlocal[12] * tau2[1]);
  qglobal[12] = (qlocal[10] * norm[2]) + (qlocal[11] * tau1[2]) + (qlocal[12] * tau2[2]);

  qglobal[13] = (qlocal[13] * norm[0]) + (qlocal[14] * tau1[0]) + (qlocal[15] * tau2[0]);
  qglobal[14] = (qlocal[13] * norm[1]) + (qlocal[14] * tau1[1]) + (qlocal[15] * tau2[1]);
  qglobal[15] = (qlocal[13] * norm[2]) + (qlocal[14] * tau1[2]) + (qlocal[15] * tau2[2]);

  qglobal[16] = qlocal[16];
  qglobal[17] = qlocal[17];

  qglobal[18] = qlocal[18];
  qglobal[19] = (qlocal[19] * norm[0]) + (qlocal[20] * tau1[0]) + (qlocal[21] * tau2[0]);
  qglobal[20] = (qlocal[19] * norm[1]) + (qlocal[20] * tau1[1]) + (qlocal[21] * tau2[1]);
  qglobal[21] = (qlocal[19] * norm[2]) + (qlocal[20] * tau1[2]) + (qlocal[21] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[22] * norm[0]) + (qlocal[23] * tau1[0]) + (qlocal[24] * tau2[0]);
  r1[1] = (qlocal[22] * norm[1]) + (qlocal[23] * tau1[1]) + (qlocal[24] * tau2[1]);
  r1[2] = (qlocal[22] * norm[2]) + (qlocal[23] * tau1[2]) + (qlocal[24] * tau2[2]);

  r2[0] = (qlocal[25] * norm[0]) + (qlocal[26] * tau1[0]) + (qlocal[27] * tau2[0]);
  r2[1] = (qlocal[25] * norm[1]) + (qlocal[26] * tau1[1]) + (qlocal[27] * tau2[1]);
  r2[2] = (qlocal[25] * norm[2]) + (qlocal[26] * tau1[2]) + (qlocal[27] * tau2[2]);

  r3[0] = (qlocal[28] * norm[0]) + (qlocal[29] * tau1[0]) + (qlocal[30] * tau2[0]);
  r3[1] = (qlocal[28] * norm[1]) + (qlocal[29] * tau1[1]) + (qlocal[30] * tau2[1]);
  r3[2] = (qlocal[28] * norm[2]) + (qlocal[29] * tau1[2]) + (qlocal[30] * tau2[2]);

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
  qglobal[22] = v1[0]; qglobal[23] = v1[1]; qglobal[24] = v1[2];
  qglobal[25] = v2[0]; qglobal[26] = v2[1]; qglobal[27] = v2[2];
  qglobal[28] = v3[0]; qglobal[29] = v3[1]; qglobal[30] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[31] * norm[0]) + (qlocal[32] * tau1[0]) + (qlocal[33] * tau2[0]);
  extr_r1[1] = (qlocal[31] * norm[1]) + (qlocal[32] * tau1[1]) + (qlocal[33] * tau2[1]);
  extr_r1[2] = (qlocal[31] * norm[2]) + (qlocal[32] * tau1[2]) + (qlocal[33] * tau2[2]);

  extr_r2[0] = (qlocal[34] * norm[0]) + (qlocal[35] * tau1[0]) + (qlocal[36] * tau2[0]);
  extr_r2[1] = (qlocal[34] * norm[1]) + (qlocal[35] * tau1[1]) + (qlocal[36] * tau2[1]);
  extr_r2[2] = (qlocal[34] * norm[2]) + (qlocal[35] * tau1[2]) + (qlocal[36] * tau2[2]);

  extr_r3[0] = (qlocal[37] * norm[0]) + (qlocal[38] * tau1[0]) + (qlocal[39] * tau2[0]);
  extr_r3[1] = (qlocal[37] * norm[1]) + (qlocal[38] * tau1[1]) + (qlocal[39] * tau2[1]);
  extr_r3[2] = (qlocal[37] * norm[2]) + (qlocal[38] * tau1[2]) + (qlocal[39] * tau2[2]);

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
  qglobal[31] = inv_v1[0]; qglobal[32] = inv_v1[1]; qglobal[33] = inv_v1[2];
  qglobal[34] = inv_v2[0]; qglobal[35] = inv_v2[1]; qglobal[36] = inv_v2[2];
  qglobal[37] = inv_v3[0]; qglobal[38] = inv_v3[1]; qglobal[39] = inv_v3[2];

  qglobal[40] = qlocal[40];

  qglobal[41] = (qlocal[41] * norm[0]) + (qlocal[42] * tau1[0]) + (qlocal[43] * tau2[0]);
  qglobal[42] = (qlocal[41] * norm[1]) + (qlocal[42] * tau1[1]) + (qlocal[43] * tau2[1]);
  qglobal[43] = (qlocal[41] * norm[2]) + (qlocal[42] * tau1[2]) + (qlocal[43] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qlocal[44] * norm[0]) + (qlocal[45] * tau1[0]) + (qlocal[46] * tau2[0]);
  shiftder_r1[1] = (qlocal[44] * norm[1]) + (qlocal[45] * tau1[1]) + (qlocal[46] * tau2[1]);
  shiftder_r1[2] = (qlocal[44] * norm[2]) + (qlocal[45] * tau1[2]) + (qlocal[46] * tau2[2]);

  shiftder_r2[0] = (qlocal[47] * norm[0]) + (qlocal[48] * tau1[0]) + (qlocal[49] * tau2[0]);
  shiftder_r2[1] = (qlocal[47] * norm[1]) + (qlocal[48] * tau1[1]) + (qlocal[49] * tau2[1]);
  shiftder_r2[2] = (qlocal[47] * norm[2]) + (qlocal[48] * tau1[2]) + (qlocal[49] * tau2[2]);

  shiftder_r3[0] = (qlocal[50] * norm[0]) + (qlocal[51] * tau1[0]) + (qlocal[52] * tau2[0]);
  shiftder_r3[1] = (qlocal[50] * norm[1]) + (qlocal[51] * tau1[1]) + (qlocal[52] * tau2[1]);
  shiftder_r3[2] = (qlocal[50] * norm[2]) + (qlocal[51] * tau1[2]) + (qlocal[52] * tau2[2]);

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
  qglobal[44] = shiftder_v1[0]; qglobal[45] = shiftder_v1[1]; qglobal[46] = shiftder_v1[2];
  qglobal[47] = shiftder_v2[0]; qglobal[48] = shiftder_v2[1]; qglobal[49] = shiftder_v2[2];
  qglobal[50] = shiftder_v3[0]; qglobal[51] = shiftder_v3[1]; qglobal[52] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qlocal[53] * norm[0]) + (qlocal[54] * tau1[0]) + (qlocal[55] * tau2[0]);
  r11[1] = (qlocal[53] * norm[1]) + (qlocal[54] * tau1[1]) + (qlocal[55] * tau2[1]);
  r11[2] = (qlocal[53] * norm[2]) + (qlocal[54] * tau1[2]) + (qlocal[55] * tau2[2]);

  r12[0] = (qlocal[56] * norm[0]) + (qlocal[57] * tau1[0]) + (qlocal[58] * tau2[0]);
  r12[1] = (qlocal[56] * norm[1]) + (qlocal[57] * tau1[1]) + (qlocal[58] * tau2[1]);
  r12[2] = (qlocal[56] * norm[2]) + (qlocal[57] * tau1[2]) + (qlocal[58] * tau2[2]);

  r13[0] = (qlocal[59] * norm[0]) + (qlocal[60] * tau1[0]) + (qlocal[61] * tau2[0]);
  r13[1] = (qlocal[59] * norm[1]) + (qlocal[60] * tau1[1]) + (qlocal[61] * tau2[1]);
  r13[2] = (qlocal[59] * norm[2]) + (qlocal[60] * tau1[2]) + (qlocal[61] * tau2[2]);

  r21[0] = (qlocal[62] * norm[0]) + (qlocal[63] * tau1[0]) + (qlocal[64] * tau2[0]);
  r21[1] = (qlocal[62] * norm[1]) + (qlocal[63] * tau1[1]) + (qlocal[64] * tau2[1]);
  r21[2] = (qlocal[62] * norm[2]) + (qlocal[63] * tau1[2]) + (qlocal[64] * tau2[2]);

  r22[0] = (qlocal[65] * norm[0]) + (qlocal[66] * tau1[0]) + (qlocal[67] * tau2[0]);
  r22[1] = (qlocal[65] * norm[1]) + (qlocal[66] * tau1[1]) + (qlocal[67] * tau2[1]);
  r22[2] = (qlocal[65] * norm[2]) + (qlocal[66] * tau1[2]) + (qlocal[67] * tau2[2]);

  r23[0] = (qlocal[68] * norm[0]) + (qlocal[69] * tau1[0]) + (qlocal[70] * tau2[0]);
  r23[1] = (qlocal[68] * norm[1]) + (qlocal[69] * tau1[1]) + (qlocal[70] * tau2[1]);
  r23[2] = (qlocal[68] * norm[2]) + (qlocal[69] * tau1[2]) + (qlocal[70] * tau2[2]);

  r31[0] = (qlocal[71] * norm[0]) + (qlocal[72] * tau1[0]) + (qlocal[73] * tau2[0]);
  r31[1] = (qlocal[71] * norm[1]) + (qlocal[72] * tau1[1]) + (qlocal[73] * tau2[1]);
  r31[2] = (qlocal[71] * norm[2]) + (qlocal[72] * tau1[2]) + (qlocal[73] * tau2[2]);

  r32[0] = (qlocal[74] * norm[0]) + (qlocal[75] * tau1[0]) + (qlocal[76] * tau2[0]);
  r32[1] = (qlocal[74] * norm[1]) + (qlocal[75] * tau1[1]) + (qlocal[76] * tau2[1]);
  r32[2] = (qlocal[74] * norm[2]) + (qlocal[75] * tau1[2]) + (qlocal[76] * tau2[2]);

  r33[0] = (qlocal[77] * norm[0]) + (qlocal[78] * tau1[0]) + (qlocal[79] * tau2[0]);
  r33[1] = (qlocal[77] * norm[1]) + (qlocal[78] * tau1[1]) + (qlocal[79] * tau2[1]);
  r33[2] = (qlocal[77] * norm[2]) + (qlocal[78] * tau1[2]) + (qlocal[79] * tau2[2]);

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
  qglobal[53] = (s11[0] * norm[0]) + (s21[0] * tau1[0]) + (s31[0] * tau2[0]);
  qglobal[54] = (s11[1] * norm[0]) + (s21[1] * tau1[0]) + (s31[1] * tau2[0]);
  qglobal[55] = (s11[2] * norm[0]) + (s21[2] * tau1[0]) + (s31[2] * tau2[0]);

  qglobal[56] = (s12[0] * norm[0]) + (s22[0] * tau1[0]) + (s32[0] * tau2[0]);
  qglobal[57] = (s12[1] * norm[0]) + (s22[1] * tau1[0]) + (s32[1] * tau2[0]);
  qglobal[58] = (s12[2] * norm[0]) + (s22[2] * tau1[0]) + (s32[2] * tau2[0]);

  qglobal[59] = (s13[0] * norm[0]) + (s23[0] * tau1[0]) + (s33[0] * tau2[0]);
  qglobal[60] = (s13[1] * norm[0]) + (s23[1] * tau1[0]) + (s33[1] * tau2[0]);
  qglobal[61] = (s13[2] * norm[0]) + (s23[2] * tau1[0]) + (s33[2] * tau2[0]);

  qglobal[62] = (s11[0] * norm[1]) + (s21[0] * tau1[1]) + (s31[0] * tau2[1]);
  qglobal[63] = (s11[1] * norm[1]) + (s21[1] * tau1[1]) + (s31[1] * tau2[1]);
  qglobal[64] = (s11[2] * norm[1]) + (s21[2] * tau1[1]) + (s31[2] * tau2[1]);

  qglobal[65] = (s12[0] * norm[1]) + (s22[0] * tau1[1]) + (s32[0] * tau2[1]);
  qglobal[66] = (s12[1] * norm[1]) + (s22[1] * tau1[1]) + (s32[1] * tau2[1]);
  qglobal[67] = (s12[2] * norm[1]) + (s22[2] * tau1[1]) + (s32[2] * tau2[1]);

  qglobal[68] = (s13[0] * norm[1]) + (s23[0] * tau1[1]) + (s33[0] * tau2[1]);
  qglobal[69] = (s13[1] * norm[1]) + (s23[1] * tau1[1]) + (s33[1] * tau2[1]);
  qglobal[70] = (s13[2] * norm[1]) + (s23[2] * tau1[1]) + (s33[2] * tau2[1]);

  qglobal[71] = (s11[0] * norm[2]) + (s21[0] * tau1[2]) + (s31[0] * tau2[2]);
  qglobal[72] = (s11[1] * norm[2]) + (s21[1] * tau1[2]) + (s31[1] * tau2[2]);
  qglobal[73] = (s11[2] * norm[2]) + (s21[2] * tau1[2]) + (s31[2] * tau2[2]);

  qglobal[74] = (s12[0] * norm[2]) + (s22[0] * tau1[2]) + (s32[0] * tau2[2]);
  qglobal[75] = (s12[1] * norm[2]) + (s22[1] * tau1[2]) + (s32[1] * tau2[2]);
  qglobal[76] = (s12[2] * norm[2]) + (s22[2] * tau1[2]) + (s32[2] * tau2[2]);

  qglobal[77] = (s13[0] * norm[2]) + (s23[0] * tau1[2]) + (s33[0] * tau2[2]);
  qglobal[78] = (s13[1] * norm[2]) + (s23[1] * tau1[2]) + (s33[1] * tau2[2]);
  qglobal[79] = (s13[2] * norm[2]) + (s23[2] * tau1[2]) + (s33[2] * tau2[2]);

  qglobal[80] = qlocal[80];
  qglobal[81] = (qlocal[81] * norm[0]) + (qlocal[82] * tau1[0]) + (qlocal[83] * tau2[0]);
  qglobal[82] = (qlocal[81] * norm[1]) + (qlocal[82] * tau1[1]) + (qlocal[83] * tau2[1]);
  qglobal[83] = (qlocal[81] * norm[2]) + (qlocal[82] * tau1[2]) + (qlocal[83] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;
  double light_speed = gr_twofluid_tetrad->light_speed;
  double e_fact = gr_twofluid_tetrad->e_fact;
  double b_fact = gr_twofluid_tetrad->b_fact;

  double sl = gkyl_gr_twofluid_tetrad_max_abs_speed(gas_gamma_elc, gas_gamma_ion, light_speed, ql);
  double sr = gkyl_gr_twofluid_tetrad_max_abs_speed(gas_gamma_elc, gas_gamma_ion, light_speed, qr);
  double amax = fmax(sl, sr);

  double fl_sr[84], fr_sr[84];
  gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql, fl_sr);
  gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr, fr_sr);

  double fl_gr[84], fr_gr[84];
  gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr, fr_sr, fr_gr);

  bool in_excision_region_l = false;
  if (ql[40] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[40] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[84];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 84; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr_gr[i] - fl_gr[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr_gr[i] - fl_gr[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 84; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[84];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 84; i++) {
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
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;
  double light_speed = gr_twofluid_tetrad->light_speed;
  double e_fact = gr_twofluid_tetrad->e_fact;
  double b_fact = gr_twofluid_tetrad->b_fact;

  double vl[84], vr[84];
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, ql, vl);
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, qr, vr);

  double rho_elc_l = vl[0];
  double vx_elc_l = vl[1];
  double vy_elc_l = vl[2];
  double vz_elc_l = vl[3];
  double p_elc_l = vl[4];

  double rho_ion_l = vl[5];
  double vx_ion_l = vl[6];
  double vy_ion_l = vl[7];
  double vz_ion_l = vl[8];
  double p_ion_l = vl[9];

  double lapse_l = vl[18];
  double shift_x_l = vl[19];
  double shift_y_l = vl[20];
  double shift_z_l = vl[21];

  double spatial_metric_l[3][3];
  spatial_metric_l[0][0] = vl[22]; spatial_metric_l[0][1] = vl[23]; spatial_metric_l[0][2] = vl[24];
  spatial_metric_l[1][0] = vl[25]; spatial_metric_l[1][1] = vl[26]; spatial_metric_l[1][2] = vl[27];
  spatial_metric_l[2][0] = vl[28]; spatial_metric_l[2][1] = vl[29]; spatial_metric_l[2][2] = vl[30];

  double spatial_metric_det_l = (spatial_metric_l[0][0] * ((spatial_metric_l[1][1] * spatial_metric_l[2][2]) - (spatial_metric_l[2][1] * spatial_metric_l[1][2]))) -
    (spatial_metric_l[0][1] * ((spatial_metric_l[1][0] * spatial_metric_l[2][2]) - (spatial_metric_l[1][2] * spatial_metric_l[2][0]))) +
    (spatial_metric_l[0][2] * ((spatial_metric_l[1][0] * spatial_metric_l[2][1]) - (spatial_metric_l[1][1] * spatial_metric_l[2][0])));

  double **inv_spatial_metric_l= gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_tetrad_inv_spatial_metric(ql, &inv_spatial_metric_l);

  bool in_excision_region_l = false;
  if (vl[40] < pow(10.0, -8.0)) {
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
  if (fabs(lapse_l - 1.0) > pow(10.0, -8.0) || fabs(shift_x_l) > pow(10.0, -8.0) || fabs(shift_y_l) > pow(10.0, -8.0) ||
    fabs(shift_z_l) > pow(10.0, -8.0)) {
    curved_spacetime_l = true;
  }

  double rho_elc_r = vr[0];
  double vx_elc_r = vr[1];
  double vy_elc_r = vr[2];
  double vz_elc_r = vr[3];
  double p_elc_r = vr[4];

  double rho_ion_r = vr[5];
  double vx_ion_r = vr[6];
  double vy_ion_r = vr[7];
  double vz_ion_r = vr[8];
  double p_ion_r = vr[9];

  double lapse_r = vr[18];
  double shift_x_r = vr[19];
  double shift_y_r = vr[20];
  double shift_z_r = vr[21];

  double spatial_metric_r[3][3];
  spatial_metric_r[0][0] = vr[22]; spatial_metric_r[0][1] = vr[23]; spatial_metric_r[0][2] = vr[24];
  spatial_metric_r[1][0] = vr[25]; spatial_metric_r[1][1] = vr[26]; spatial_metric_r[1][2] = vr[27];
  spatial_metric_r[2][0] = vr[28]; spatial_metric_r[2][1] = vr[29]; spatial_metric_r[2][2] = vr[30];

  double spatial_metric_det_r = (spatial_metric_r[0][0] * ((spatial_metric_r[1][1] * spatial_metric_r[2][2]) - (spatial_metric_r[2][1] * spatial_metric_r[1][2]))) -
    (spatial_metric_r[0][1] * ((spatial_metric_r[1][0] * spatial_metric_r[2][2]) - (spatial_metric_r[1][2] * spatial_metric_r[2][0]))) +
    (spatial_metric_r[0][2] * ((spatial_metric_r[1][0] * spatial_metric_r[2][1]) - (spatial_metric_r[1][1] * spatial_metric_r[2][0])));

  double **inv_spatial_metric_r= gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_tetrad_inv_spatial_metric(ql, &inv_spatial_metric_r);

  bool in_excision_region_r = false;
  if (vr[40] < pow(10.0, -8.0)) {
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
        if (fabs(spatial_metric_r[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime_r = true;
        }
      }
    }
  }
  if (fabs(lapse_r - 1.0) > pow(10.0, -8.0) || fabs(shift_x_r) > pow(10.0, -8.0) || fabs(shift_y_r) > pow(10.0, -8.0) ||
    fabs(shift_z_r) > pow(10.0, -8.0)) {
    curved_spacetime_r = true;
  }

  double num_elc_l = (gas_gamma_elc * p_elc_l) / rho_elc_l;
  double den_elc_l = 1.0 + ((p_elc_l / rho_elc_l) * (gas_gamma_elc) / (gas_gamma_elc - 1.0));
  double c_s_elc_l = sqrt(num_elc_l / den_elc_l);

  double num_ion_l = (gas_gamma_ion * p_ion_l) / rho_ion_l;
  double den_ion_l = 1.0 + ((p_ion_l / rho_ion_l) * (gas_gamma_ion) / (gas_gamma_ion - 1.0));
  double c_s_ion_l = sqrt(num_ion_l / den_ion_l);

  double num_elc_r = (gas_gamma_elc * p_elc_r) / rho_elc_r;
  double den_elc_r = 1.0 + ((p_elc_r / rho_elc_r) * (gas_gamma_elc) / (gas_gamma_elc - 1.0));
  double c_s_elc_r = sqrt(num_elc_r / den_elc_r);

  double num_ion_r = (gas_gamma_ion * p_ion_r) / rho_ion_r;
  double den_ion_r = 1.0 + ((p_ion_r / rho_ion_r) * (gas_gamma_ion) / (gas_gamma_ion - 1.0));
  double c_s_ion_r = sqrt(num_ion_r / den_ion_r);

  double vx_avg_elc = 0.5 * (vx_elc_l + vx_elc_r);
  double cs_avg_elc = 0.5 * (c_s_elc_l + c_s_elc_r);

  double vx_avg_ion = 0.5 * (vx_ion_l + vx_ion_r);
  double cs_avg_ion = 0.5 * (c_s_ion_l + c_s_ion_r);

  double sl_elc, sr_elc;
  double sl_ion, sr_ion;
  double sl_em, sr_em;

  if (curved_spacetime_l || curved_spacetime_r) {
    double vel_elc_l[3];
    double v_sq_elc_l = 0.0;
    vel_elc_l[0] = vx_elc_l; vel_elc_l[1] = vy_elc_l; vel_elc_l[2] = vz_elc_l;
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_elc_l += spatial_metric_l[i][j] * vel_elc_l[i] * vel_elc_l[j];
      }
    }

    double shift_l[3];
    shift_l[0] = shift_x_l; shift_l[1] = shift_y_l; shift_l[2] = shift_z_l;

    double material_eigs_elc_l[3];
    double fast_acoustic_eigs_elc_l[3];
    double slow_acoustic_eigs_elc_l[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_elc_l[i] = (lapse_l * vel_elc_l[i]) - shift_l[i];

      fast_acoustic_eigs_elc_l[i] = (lapse_l / (1.0 - (v_sq_elc_l * (c_s_elc_l * c_s_elc_l)))) * ((vel_elc_l[i] * (1.0 - (c_s_elc_l * c_s_elc_l))) +
        (c_s_elc_l * sqrt((1.0 - v_sq_elc_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_elc_l * (c_s_elc_l * c_s_elc_l))) -
        (vel_elc_l[i] * vel_elc_l[i]) * (1.0 - (c_s_elc_l * c_s_elc_l)))))) - shift_l[i];
      
      slow_acoustic_eigs_elc_l[i] = (lapse_l / (1.0 - (v_sq_elc_l * (c_s_elc_l * c_s_elc_l)))) * ((vel_elc_l[i] * (1.0 - (c_s_elc_l * c_s_elc_l))) -
        (c_s_elc_l * sqrt((1.0 - v_sq_elc_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_elc_l * (c_s_elc_l * c_s_elc_l))) -
        (vel_elc_l[i] * vel_elc_l[i]) * (1.0 - (c_s_elc_l * c_s_elc_l)))))) - shift_l[i];
    }

    double max_eig_elc_l = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_elc_l[i]) > max_eig_elc_l) {
        max_eig_elc_l = fabs(material_eigs_elc_l[i]);
      }
      if (fabs(fast_acoustic_eigs_elc_l[i]) > max_eig_elc_l) {
        max_eig_elc_l = fabs(fast_acoustic_eigs_elc_l[i]);
      }
      if (fabs(slow_acoustic_eigs_elc_l[i]) > max_eig_elc_l) {
        max_eig_elc_l = fabs(slow_acoustic_eigs_elc_l[i]);
      }
    }

    double vel_elc_r[3];
    double v_sq_elc_r = 0.0;
    vel_elc_r[0] = vx_elc_r; vel_elc_r[1] = vy_elc_r; vel_elc_r[2] = vz_elc_r;
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_elc_r += spatial_metric_r[i][j] * vel_elc_r[i] * vel_elc_r[j];
      }
    }

    double shift_r[3];
    shift_r[0] = shift_x_r; shift_r[1] = shift_y_r; shift_r[2] = shift_z_r;

    double material_eigs_elc_r[3];
    double fast_acoustic_eigs_elc_r[3];
    double slow_acoustic_eigs_elc_r[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_elc_r[i] = (lapse_r * vel_elc_r[i]) - shift_r[i];

      fast_acoustic_eigs_elc_r[i] = (lapse_r / (1.0 - (v_sq_elc_r * (c_s_elc_r * c_s_elc_r)))) * ((vel_elc_r[i] * (1.0 - (c_s_elc_r * c_s_elc_r))) +
        (c_s_elc_r * sqrt((1.0 - v_sq_elc_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_elc_r * (c_s_elc_r * c_s_elc_r))) -
        (vel_elc_r[i] * vel_elc_r[i]) * (1.0 - (c_s_elc_r * c_s_elc_r)))))) - shift_r[i];
      
      slow_acoustic_eigs_elc_r[i] = (lapse_r / (1.0 - (v_sq_elc_r * (c_s_elc_r * c_s_elc_r)))) * ((vel_elc_r[i] * (1.0 - (c_s_elc_r * c_s_elc_r))) -
        (c_s_elc_r * sqrt((1.0 - v_sq_elc_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_elc_r * (c_s_elc_r * c_s_elc_r))) -
        (vel_elc_r[i] * vel_elc_r[i]) * (1.0 - (c_s_elc_r * c_s_elc_r)))))) - shift_r[i];
    }

    double max_eig_elc_r = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_elc_r[i]) > max_eig_elc_r) {
        max_eig_elc_r = fabs(material_eigs_elc_r[i]);
      }
      if (fabs(fast_acoustic_eigs_elc_r[i]) > max_eig_elc_r) {
        max_eig_elc_r = fabs(fast_acoustic_eigs_elc_r[i]);
      }
      if (fabs(slow_acoustic_eigs_elc_r[i]) > max_eig_elc_r) {
        max_eig_elc_r = fabs(slow_acoustic_eigs_elc_r[i]);
      }
    }

    double max_eig_avg_elc = 0.5 * (max_eig_elc_l + max_eig_elc_r);

    sl_elc = (vx_avg_elc - max_eig_avg_elc) / (1.0 - (vx_avg_elc * max_eig_avg_elc));
    sr_elc = (vx_avg_elc + max_eig_avg_elc) / (1.0 + (vx_avg_elc * max_eig_avg_elc));

    double vel_ion_l[3];
    double v_sq_ion_l = 0.0;
    vel_ion_l[0] = vx_ion_l; vel_ion_l[1] = vy_ion_l; vel_ion_l[2] = vz_ion_l;
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_ion_l += spatial_metric_l[i][j] * vel_ion_l[i] * vel_ion_l[j];
      }
    }

    double material_eigs_ion_l[3];
    double fast_acoustic_eigs_ion_l[3];
    double slow_acoustic_eigs_ion_l[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_ion_l[i] = (lapse_l * vel_ion_l[i]) - shift_l[i];

      fast_acoustic_eigs_ion_l[i] = (lapse_l / (1.0 - (v_sq_ion_l * (c_s_ion_l * c_s_ion_l)))) * ((vel_ion_l[i] * (1.0 - (c_s_ion_l * c_s_ion_l))) +
        (c_s_ion_l * sqrt((1.0 - v_sq_ion_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_ion_l * (c_s_ion_l * c_s_ion_l))) -
        (vel_ion_l[i] * vel_ion_l[i]) * (1.0 - (c_s_ion_l * c_s_ion_l)))))) - shift_l[i];
      
      slow_acoustic_eigs_ion_l[i] = (lapse_l / (1.0 - (v_sq_ion_l * (c_s_ion_l * c_s_ion_l)))) * ((vel_ion_l[i] * (1.0 - (c_s_ion_l * c_s_ion_l))) -
        (c_s_ion_l * sqrt((1.0 - v_sq_ion_l) * (inv_spatial_metric_l[i][i] * (1.0 - (v_sq_ion_l * (c_s_ion_l * c_s_ion_l))) -
        (vel_ion_l[i] * vel_ion_l[i]) * (1.0 - (c_s_ion_l * c_s_ion_l)))))) - shift_l[i];
    }

    double max_eig_ion_l = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_ion_l[i]) > max_eig_ion_l) {
        max_eig_ion_l = fabs(material_eigs_ion_l[i]);
      }
      if (fabs(fast_acoustic_eigs_ion_l[i]) > max_eig_ion_l) {
        max_eig_ion_l = fabs(fast_acoustic_eigs_ion_l[i]);
      }
      if (fabs(slow_acoustic_eigs_ion_l[i]) > max_eig_ion_l) {
        max_eig_ion_l = fabs(slow_acoustic_eigs_ion_l[i]);
      }
    }

    double vel_ion_r[3];
    double v_sq_ion_r = 0.0;
    vel_ion_r[0] = vx_ion_r; vel_ion_r[1] = vy_ion_r; vel_ion_r[2] = vz_ion_r;
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_ion_r += spatial_metric_r[i][j] * vel_ion_r[i] * vel_ion_r[j];
      }
    }

    double material_eigs_ion_r[3];
    double fast_acoustic_eigs_ion_r[3];
    double slow_acoustic_eigs_ion_r[3];

    for (int i = 0; i < 3; i++) {
      material_eigs_ion_r[i] = (lapse_r * vel_ion_r[i]) - shift_r[i];

      fast_acoustic_eigs_ion_r[i] = (lapse_r / (1.0 - (v_sq_ion_r * (c_s_ion_r * c_s_ion_r)))) * ((vel_ion_r[i] * (1.0 - (c_s_ion_r * c_s_ion_r))) +
        (c_s_ion_r * sqrt((1.0 - v_sq_ion_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_ion_r * (c_s_ion_r * c_s_ion_r))) -
        (vel_ion_r[i] * vel_ion_r[i]) * (1.0 - (c_s_ion_r * c_s_ion_r)))))) - shift_r[i];
      
      slow_acoustic_eigs_ion_r[i] = (lapse_r / (1.0 - (v_sq_ion_r * (c_s_ion_r * c_s_ion_r)))) * ((vel_ion_r[i] * (1.0 - (c_s_ion_r * c_s_ion_r))) -
        (c_s_ion_r * sqrt((1.0 - v_sq_ion_r) * (inv_spatial_metric_r[i][i] * (1.0 - (v_sq_ion_r * (c_s_ion_r * c_s_ion_r))) -
        (vel_ion_r[i] * vel_ion_r[i]) * (1.0 - (c_s_ion_r * c_s_ion_r)))))) - shift_r[i];
    }

    double max_eig_ion_r = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(material_eigs_ion_r[i]) > max_eig_ion_r) {
        max_eig_ion_r = fabs(material_eigs_ion_r[i]);
      }
      if (fabs(fast_acoustic_eigs_ion_r[i]) > max_eig_ion_r) {
        max_eig_ion_r = fabs(fast_acoustic_eigs_ion_r[i]);
      }
      if (fabs(slow_acoustic_eigs_ion_r[i]) > max_eig_ion_r) {
        max_eig_ion_r = fabs(slow_acoustic_eigs_ion_r[i]);
      }
    }

    double max_eig_avg_ion = 0.5 * (max_eig_ion_l + max_eig_ion_r);

    sl_ion = (vx_avg_ion - max_eig_avg_ion) / (1.0 - (vx_avg_ion * max_eig_avg_ion));
    sr_ion = (vx_avg_ion + max_eig_avg_ion) / (1.0 + (vx_avg_ion * max_eig_avg_ion));

    sl_em = -light_speed * sqrt(spatial_metric_det_l) * lapse_l;
    sr_em = light_speed * sqrt(spatial_metric_det_r) * lapse_r;
  }
  else {
    sl_elc = (vx_avg_elc - cs_avg_elc) / (1.0 - (vx_avg_elc * cs_avg_elc));
    sr_elc = (vx_avg_elc + cs_avg_elc) / (1.0 + (vx_avg_elc * cs_avg_elc));

    sl_ion = (vx_avg_ion - cs_avg_ion) / (1.0 - (vx_avg_ion * cs_avg_ion));
    sr_ion = (vx_avg_ion + cs_avg_ion) / (1.0 + (vx_avg_ion * cs_avg_ion));

    sl_em = -light_speed * sqrt(spatial_metric_det_l) * lapse_l;
    sr_em = light_speed * sqrt(spatial_metric_det_r) * lapse_r;
  }

  double fl_sr[84], fr_sr[84];
  gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql, fl_sr);
  gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr, fr_sr);

  double fl_gr[84], fr_gr[84];
  gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr, fr_sr, fr_gr);

  double qm[84];
  for (int i = 0; i < 5; i++) {
    qm[i] = ((sr_elc * qr[i]) - (sl_elc * ql[i]) + (fl_gr[i] - fr_gr[i])) / (sr_elc - sl_elc);
  }
  for (int i = 5; i < 10; i++) {
    qm[i] = ((sr_ion * qr[i]) - (sl_ion * ql[i]) + (fl_gr[i] - fr_gr[i])) / (sr_ion - sl_ion);
  }
  for (int i = 10; i < 84; i++) {
    qm[i] = ((sr_em * qr[i]) - (sl_em * ql[i]) + (fl_gr[i] - fr_gr[i])) / (sr_em - sl_em);
  }

  double *w0 = &waves[0 * 84], *w1 = &waves[1 * 84], *w2 = &waves[2 * 84], *w3 = &waves[3 * 84], *w4 = &waves[4 * 84], *w5 = &waves[5 * 84];

  for (int i = 0; i < 84; i++) {
    w0[i] = 0.0;
    w1[i] = 0.0;
    w2[i] = 0.0;
    w3[i] = 0.0;
    w4[i] = 0.0;
    w5[i] = 0.0;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 5; i++) {
      w0[i] = qm[i] - ql[i];
      w1[i] = qr[i] - qm[i];
    }

    for (int i = 5; i < 10; i++) {
      w2[i] = qm[i] - ql[i];
      w3[i] = qr[i] - qm[i];
    }

    for (int i = 10; i < 84; i++) {
      w4[i] = qm[i] - ql[i];
      w5[i] = qr[i] - qm[i];
    }
  }

  s[0] = sl_elc;
  s[1] = sr_elc;
  s[2] = sl_ion;
  s[3] = sr_ion;
  s[4] = sl_em;
  s[5] = sr_em;

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric_l[i]);
    gkyl_free(inv_spatial_metric_r[i]);
  }
  gkyl_free(inv_spatial_metric_l);
  gkyl_free(inv_spatial_metric_r);

  return fmax(fmax(fmax(fabs(sl_elc), fabs(sr_elc)), fmax(fabs(sl_ion), fabs(sr_ion))), fmax(fabs(sl_em), fabs(sr_em)));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0 * 84], *w1 = &waves[1 * 84], *w2 = &waves[2 * 84], *w3 = &waves[3 * 84], *w4 = &waves[4 * 84], *w5 = &waves[5 * 84];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]), s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]), s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  for (int i = 0; i < 84; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]) + (s3m * w3[i]) + (s4m * w4[i]) + (s5m * w5[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]) + (s3p * w3[i]) + (s4p * w4[i]) + (s5p * w5[i]);
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
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;
  double light_speed = gr_twofluid_tetrad->light_speed;
  double e_fact = gr_twofluid_tetrad->e_fact;
  double b_fact = gr_twofluid_tetrad->b_fact;

  double fr_sr[84], fl_sr[84];
  gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql, fl_sr);
  gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr, fr_sr);

  double fr_gr[84], fl_gr[84];
  gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql, fl_sr, fl_gr);
  gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr, fr_sr, fr_gr);

  bool in_excision_region_l = false;
  if (ql[40] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[40] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 84; m++) {
      flux_jump[m] = fr_gr[m] - fl_gr[m];
    }
  }
  else {
    for (int m = 0; m < 84; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_twofluid_tetrad_max_abs_speed(gas_gamma_elc, gas_gamma_ion, light_speed, ql);
  double amaxr = gkyl_gr_twofluid_tetrad_max_abs_speed(gas_gamma_elc, gas_gamma_ion, light_speed, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;

  double v[84] = { 0.0 };
  gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);

  if (v[0] < 0.0 || v[4] < 0.0 || v[5] < 0.0 || v[9] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;
  double light_speed = gr_twofluid_tetrad->light_speed;

  return gkyl_gr_twofluid_tetrad_max_abs_speed(gas_gamma_elc, gas_gamma_ion, light_speed, q);
}

static inline void
gr_twofluid_tetrad_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 5; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_twofluid_tetrad_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 84; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_gr_twofluid_tetrad_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(base->on_dev, struct wv_gr_twofluid_tetrad, eqn);
    gkyl_cu_free(gr_twofluid_tetrad);
  }

  struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(base, struct wv_gr_twofluid_tetrad, eqn);
  gkyl_free(gr_twofluid_tetrad);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_twofluid_tetrad_new(double mass_elc, double mass_ion, double charge_elc, double charge_ion, double gas_gamma_elc, double gas_gamma_ion,
  double light_speed, double e_fact, double b_fact, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime,
  bool use_gpu)
{
  return gkyl_wv_gr_twofluid_tetrad_inew(&(struct gkyl_wv_gr_twofluid_tetrad_inp) {
      .mass_elc = mass_elc,
      .mass_ion = mass_ion,
      .charge_elc = charge_elc,
      .charge_ion = charge_ion,
      .gas_gamma_elc = gas_gamma_elc,
      .gas_gamma_ion = gas_gamma_ion,
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_TWOFLUID_TETRAD_RP_HLL,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_twofluid_tetrad_inew(const struct gkyl_wv_gr_twofluid_tetrad_inp* inp)
{
  struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = gkyl_malloc(sizeof(struct wv_gr_twofluid_tetrad));

  gr_twofluid_tetrad->eqn.type = GKYL_EQN_GR_TWOFLUID_TETRAD;
  gr_twofluid_tetrad->eqn.num_equations = 84;
  gr_twofluid_tetrad->eqn.num_diag = 5;

  gr_twofluid_tetrad->mass_elc = inp->mass_elc;
  gr_twofluid_tetrad->mass_ion = inp->mass_ion;
  gr_twofluid_tetrad->charge_elc = inp->charge_elc;
  gr_twofluid_tetrad->charge_ion = inp->charge_ion;
  gr_twofluid_tetrad->gas_gamma_elc = inp->gas_gamma_elc;
  gr_twofluid_tetrad->gas_gamma_ion = inp->gas_gamma_ion;
  gr_twofluid_tetrad->light_speed = inp->light_speed;
  gr_twofluid_tetrad->e_fact = inp->e_fact;
  gr_twofluid_tetrad->b_fact = inp->b_fact;

  gr_twofluid_tetrad->spacetime_gauge = inp->spacetime_gauge;
  gr_twofluid_tetrad->reinit_freq = inp->reinit_freq;
  gr_twofluid_tetrad->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_TWOFLUID_TETRAD_RP_LAX) {
    gr_twofluid_tetrad->eqn.num_waves = 2;
    gr_twofluid_tetrad->eqn.waves_func = wave_lax_l;
    gr_twofluid_tetrad->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_TWOFLUID_TETRAD_RP_HLL) {
    gr_twofluid_tetrad->eqn.num_waves = 6;
    gr_twofluid_tetrad->eqn.waves_func = wave_hll_l;
    gr_twofluid_tetrad->eqn.qfluct_func = qfluct_hll_l;
  }

  gr_twofluid_tetrad->eqn.flux_jump = flux_jump;
  gr_twofluid_tetrad->eqn.check_inv_func = check_inv;
  gr_twofluid_tetrad->eqn.max_speed_func = max_speed;
  gr_twofluid_tetrad->eqn.rotate_to_local_func = rot_to_local;
  gr_twofluid_tetrad->eqn.rotate_to_global_func = rot_to_global;

  gr_twofluid_tetrad->eqn.wall_bc_func = gr_twofluid_tetrad_wall;
  gr_twofluid_tetrad->eqn.no_slip_bc_func = gr_twofluid_tetrad_no_slip;

  gr_twofluid_tetrad->eqn.cons_to_riem = cons_to_riem;
  gr_twofluid_tetrad->eqn.riem_to_cons = riem_to_cons;

  gr_twofluid_tetrad->eqn.cons_to_diag = gr_twofluid_tetrad_cons_to_diag;

  gr_twofluid_tetrad->eqn.source_func = gr_twofluid_tetrad_source;

  gr_twofluid_tetrad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_twofluid_tetrad->eqn.flags);
  gr_twofluid_tetrad->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_twofluid_tetrad_free);
  gr_twofluid_tetrad->eqn.on_dev = &gr_twofluid_tetrad->eqn; // On the CPU, the equation object points to itself.

  return &gr_twofluid_tetrad->eqn;
}

double
gkyl_wv_gr_twofluid_tetrad_mass_elc(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double mass_elc = gr_twofluid_tetrad->mass_elc;

  return mass_elc;
}

double
gkyl_wv_gr_twofluid_tetrad_mass_ion(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double mass_ion = gr_twofluid_tetrad->mass_ion;

  return mass_ion;
}

double
gkyl_wv_gr_twofluid_tetrad_charge_elc(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double charge_elc = gr_twofluid_tetrad->charge_elc;

  return charge_elc;
}

double
gkyl_wv_gr_twofluid_tetrad_charge_ion(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double charge_ion = gr_twofluid_tetrad->charge_ion;

  return charge_ion;
}

double
gkyl_wv_gr_twofluid_tetrad_gas_gamma_elc(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_elc = gr_twofluid_tetrad->gas_gamma_elc;

  return gas_gamma_elc;
}

double
gkyl_wv_gr_twofluid_tetrad_gas_gamma_ion(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double gas_gamma_ion = gr_twofluid_tetrad->gas_gamma_ion;

  return gas_gamma_ion;
}

double
gkyl_wv_gr_twofluid_tetrad_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double light_speed = gr_twofluid_tetrad->light_speed;

  return light_speed;
}

double
gkyl_wv_gr_twofluid_tetrad_e_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double e_fact = gr_twofluid_tetrad->e_fact;

  return e_fact;
}

double
gkyl_wv_gr_twofluid_tetrad_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  double b_fact = gr_twofluid_tetrad->b_fact;

  return b_fact;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_twofluid_tetrad_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_twofluid_tetrad->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_twofluid_tetrad_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  int reinit_freq = gr_twofluid_tetrad->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_twofluid_tetrad_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_twofluid_tetrad->spacetime;

  return spacetime;
}
