#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_twofluid_tetrad.h>
#include <gkyl_wv_gr_twofluid_tetrad_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_twofluid_tetrad_basic_minkowski_ho()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_twofluid_tetrad = gkyl_wv_gr_twofluid_tetrad_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion,
    light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_twofluid_tetrad->num_equations == 84 );
  TEST_CHECK( gr_twofluid_tetrad->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double rho_elc = rho_ion * mass_elc / mass_ion;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double he = 1.0 + ((p / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi = 1.0 + ((p / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double q[84];
      q[0] = sqrt(spatial_det) * rho_elc * W;
      q[1] = sqrt(spatial_det) * rho_elc * he * (W * W) * u;
      q[2] = sqrt(spatial_det) * rho_elc * he * (W * W) * v;
      q[3] = sqrt(spatial_det) * rho_elc * he * (W * W) * w;
      q[4] = sqrt(spatial_det) * ((rho_elc * he * (W * W)) - p - (rho_elc * W));

      q[5] = sqrt(spatial_det) * rho_ion * W;
      q[6] = sqrt(spatial_det) * rho_ion * hi * (W * W) * u;
      q[7] = sqrt(spatial_det) * rho_ion * hi * (W * W) * v;
      q[8] = sqrt(spatial_det) * rho_ion * hi * (W * W) * w;
      q[9] = sqrt(spatial_det) * ((rho_ion * hi * (W * W)) - p - (rho_ion * W));

      q[10] = Dx; q[11] = Dy; q[12] = Dz;
      q[13] = Bx; q[14] = By; q[15] = Bz;
      q[16] = phi; q[17] = psi;

      q[18] = lapse;
      q[19] = shift[0]; q[20] = shift[1]; q[21] = shift[2];

      q[22] = spatial_metric[0][0]; q[23] = spatial_metric[0][1]; q[24] = spatial_metric[0][2];
      q[25] = spatial_metric[1][0]; q[26] = spatial_metric[1][1]; q[27] = spatial_metric[1][2];
      q[28] = spatial_metric[2][0]; q[29] = spatial_metric[2][1]; q[30] = spatial_metric[2][2];

      q[31] = extrinsic_curvature[0][0]; q[32] = extrinsic_curvature[0][1]; q[33] = extrinsic_curvature[0][2];
      q[34] = extrinsic_curvature[1][0]; q[35] = extrinsic_curvature[1][1]; q[36] = extrinsic_curvature[1][2];
      q[37] = extrinsic_curvature[2][0]; q[38] = extrinsic_curvature[2][1]; q[39] = extrinsic_curvature[2][2];

      q[40] = 1.0;

      q[41] = lapse_der[0]; q[42] = lapse_der[1]; q[43] = lapse_der[2];
      q[44] = shift_der[0][0]; q[45] = shift_der[0][1]; q[46] = shift_der[0][2];
      q[47] = shift_der[1][0]; q[48] = shift_der[1][1]; q[49] = shift_der[1][2];
      q[50] = shift_der[2][0]; q[51] = shift_der[2][1]; q[52] = shift_der[2][2];

      q[53] = spatial_metric_der[0][0][0]; q[54] = spatial_metric_der[0][0][1]; q[55] = spatial_metric_der[0][0][2];
      q[56] = spatial_metric_der[0][1][0]; q[57] = spatial_metric_der[0][1][1]; q[58] = spatial_metric_der[0][1][2];
      q[59] = spatial_metric_der[0][2][0]; q[60] = spatial_metric_der[0][2][1]; q[61] = spatial_metric_der[0][2][2];

      q[62] = spatial_metric_der[1][0][0]; q[63] = spatial_metric_der[1][0][1]; q[64] = spatial_metric_der[1][0][2];
      q[65] = spatial_metric_der[1][1][0]; q[66] = spatial_metric_der[1][1][1]; q[67] = spatial_metric_der[1][1][2];
      q[68] = spatial_metric_der[1][2][0]; q[69] = spatial_metric_der[1][2][1]; q[70] = spatial_metric_der[1][2][2];

      q[71] = spatial_metric_der[2][0][0]; q[72] = spatial_metric_der[2][0][1]; q[73] = spatial_metric_der[2][0][2];
      q[74] = spatial_metric_der[2][1][0]; q[75] = spatial_metric_der[2][1][1]; q[76] = spatial_metric_der[2][1][2];
      q[77] = spatial_metric_der[2][2][0]; q[78] = spatial_metric_der[2][2][1]; q[79] = spatial_metric_der[2][2][2];

      q[80] = 0.0;
      q[81] = x; q[82] = y; q[83] = 0.0;

      double prims[84];
      gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, prims);
      
      TEST_CHECK( gkyl_compare(prims[0], rho_elc, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[1], u, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[2], v, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[3], w, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[4], p, 1e-12) );

      TEST_CHECK( gkyl_compare(prims[5], rho_ion, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[6], u, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[7], v, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[8], w, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[9], p, 1e-12) );

      double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
      double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
      double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

      double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
      double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
      double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

      double fluxes[3][18] = {
        { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
          (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
          e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
          -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
        { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
          (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
          -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
          -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
        { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
          (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
          (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
          e_fact * Dz, b_fact * (light_speed * light_speed) * Bz },
      };

      double norm[3][3] = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
      };

      double tau1[3][3] = {
        { 0.0, 1.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
      };

      double tau2[3][3] = {
        { 0.0, 0.0, 1.0 },
        { 0.0, 0.0, -1.0 },
        { 0.0, 1.0, 0.0 },
      };

      double q_local[84], flux_local_sr[84], flux_local_gr[84], flux[84];
      for (int d = 0; d < 3; d++) {
        gr_twofluid_tetrad->rotate_to_local_func(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local_sr);
        gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
        gr_twofluid_tetrad->rotate_to_global_func(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

        for (int i = 0; i < 18; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-8) );
        }
      }

      double q_l[84], q_g[84];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 84; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
        }

        double w1[84], q1[84];
        gr_twofluid_tetrad->cons_to_riem(gr_twofluid_tetrad, q_local, q_local, w1);
        gr_twofluid_tetrad->riem_to_cons(gr_twofluid_tetrad, q_local, w1, q1);

        for (int i = 0; i < 84; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_twofluid_tetrad_basic_schwarzschild_ho()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_twofluid_tetrad = gkyl_wv_gr_twofluid_tetrad_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion,
    light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_twofluid_tetrad->num_equations == 84 );
  TEST_CHECK( gr_twofluid_tetrad->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double rho_elc = rho_ion * mass_elc / mass_ion;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double he = 1.0 + ((p / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi = 1.0 + ((p / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      if (!in_excision_region) {
        double q[84];
        q[0] = sqrt(spatial_det) * rho_elc * W;
        q[1] = sqrt(spatial_det) * rho_elc * he * (W * W) * u;
        q[2] = sqrt(spatial_det) * rho_elc * he * (W * W) * v;
        q[3] = sqrt(spatial_det) * rho_elc * he * (W * W) * w;
        q[4] = sqrt(spatial_det) * ((rho_elc * he * (W * W)) - p - (rho_elc * W));

        q[5] = sqrt(spatial_det) * rho_ion * W;
        q[6] = sqrt(spatial_det) * rho_ion * hi * (W * W) * u;
        q[7] = sqrt(spatial_det) * rho_ion * hi * (W * W) * v;
        q[8] = sqrt(spatial_det) * rho_ion * hi * (W * W) * w;
        q[9] = sqrt(spatial_det) * ((rho_ion * hi * (W * W)) - p - (rho_ion * W));

        q[10] = Dx; q[11] = Dy; q[12] = Dz;
        q[13] = Bx; q[14] = By; q[15] = Bz;
        q[16] = phi; q[17] = psi;

        q[18] = lapse;
        q[19] = shift[0]; q[20] = shift[1]; q[21] = shift[2];

        q[22] = spatial_metric[0][0]; q[23] = spatial_metric[0][1]; q[24] = spatial_metric[0][2];
        q[25] = spatial_metric[1][0]; q[26] = spatial_metric[1][1]; q[27] = spatial_metric[1][2];
        q[28] = spatial_metric[2][0]; q[29] = spatial_metric[2][1]; q[30] = spatial_metric[2][2];

        q[31] = extrinsic_curvature[0][0]; q[32] = extrinsic_curvature[0][1]; q[33] = extrinsic_curvature[0][2];
        q[34] = extrinsic_curvature[1][0]; q[35] = extrinsic_curvature[1][1]; q[36] = extrinsic_curvature[1][2];
        q[37] = extrinsic_curvature[2][0]; q[38] = extrinsic_curvature[2][1]; q[39] = extrinsic_curvature[2][2];

        q[40] = 1.0;

        q[41] = lapse_der[0]; q[42] = lapse_der[1]; q[43] = lapse_der[2];
        q[44] = shift_der[0][0]; q[45] = shift_der[0][1]; q[46] = shift_der[0][2];
        q[47] = shift_der[1][0]; q[48] = shift_der[1][1]; q[49] = shift_der[1][2];
        q[50] = shift_der[2][0]; q[51] = shift_der[2][1]; q[52] = shift_der[2][2];

        q[53] = spatial_metric_der[0][0][0]; q[54] = spatial_metric_der[0][0][1]; q[55] = spatial_metric_der[0][0][2];
        q[56] = spatial_metric_der[0][1][0]; q[57] = spatial_metric_der[0][1][1]; q[58] = spatial_metric_der[0][1][2];
        q[59] = spatial_metric_der[0][2][0]; q[60] = spatial_metric_der[0][2][1]; q[61] = spatial_metric_der[0][2][2];

        q[62] = spatial_metric_der[1][0][0]; q[63] = spatial_metric_der[1][0][1]; q[64] = spatial_metric_der[1][0][2];
        q[65] = spatial_metric_der[1][1][0]; q[66] = spatial_metric_der[1][1][1]; q[67] = spatial_metric_der[1][1][2];
        q[68] = spatial_metric_der[1][2][0]; q[69] = spatial_metric_der[1][2][1]; q[70] = spatial_metric_der[1][2][2];

        q[71] = spatial_metric_der[2][0][0]; q[72] = spatial_metric_der[2][0][1]; q[73] = spatial_metric_der[2][0][2];
        q[74] = spatial_metric_der[2][1][0]; q[75] = spatial_metric_der[2][1][1]; q[76] = spatial_metric_der[2][1][2];
        q[77] = spatial_metric_der[2][2][0]; q[78] = spatial_metric_der[2][2][1]; q[79] = spatial_metric_der[2][2][2];

        q[80] = 0.0;
        q[81] = x; q[82] = y; q[83] = 0.0;

        double prims[84];
        gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho_elc, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        TEST_CHECK( gkyl_compare(prims[5], rho_ion, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[6], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[7], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[8], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[9], p, 1e-1) );

        double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
        double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
        double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

        double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
        double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
        double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

        double fluxes[3][18] = {
          { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
            (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
            e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
            -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
          { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
            (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
            -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
            -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
          { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
            (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
            (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
            e_fact * Dz, b_fact * (light_speed * light_speed) * Bz },
        };

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

      double q_local[84], flux_local_sr[84], flux_local_gr[84], flux[84];
      for (int d = 0; d < 3; d++) {
        gr_twofluid_tetrad->rotate_to_local_func(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local_sr);
        gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
        gr_twofluid_tetrad->rotate_to_global_func(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 18; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[84], q_g[84];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 84; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[84], q1[84];
          gr_twofluid_tetrad->cons_to_riem(gr_twofluid_tetrad, q_local, q_local, w1);
          gr_twofluid_tetrad->riem_to_cons(gr_twofluid_tetrad, q_local, w1, q1);

          for (int i = 0; i < 84; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_twofluid_tetrad_basic_kerr_ho()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_twofluid_tetrad = gkyl_wv_gr_twofluid_tetrad_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion,
    light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_twofluid_tetrad->num_equations == 84 );
  TEST_CHECK( gr_twofluid_tetrad->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double rho_elc = rho_ion * mass_elc / mass_ion;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double he = 1.0 + ((p / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi = 1.0 + ((p / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      if (!in_excision_region) {
        double q[84];
        q[0] = sqrt(spatial_det) * rho_elc * W;
        q[1] = sqrt(spatial_det) * rho_elc * he * (W * W) * u;
        q[2] = sqrt(spatial_det) * rho_elc * he * (W * W) * v;
        q[3] = sqrt(spatial_det) * rho_elc * he * (W * W) * w;
        q[4] = sqrt(spatial_det) * ((rho_elc * he * (W * W)) - p - (rho_elc * W));

        q[5] = sqrt(spatial_det) * rho_ion * W;
        q[6] = sqrt(spatial_det) * rho_ion * hi * (W * W) * u;
        q[7] = sqrt(spatial_det) * rho_ion * hi * (W * W) * v;
        q[8] = sqrt(spatial_det) * rho_ion * hi * (W * W) * w;
        q[9] = sqrt(spatial_det) * ((rho_ion * hi * (W * W)) - p - (rho_ion * W));

        q[10] = Dx; q[11] = Dy; q[12] = Dz;
        q[13] = Bx; q[14] = By; q[15] = Bz;
        q[16] = phi; q[17] = psi;

        q[18] = lapse;
        q[19] = shift[0]; q[20] = shift[1]; q[21] = shift[2];

        q[22] = spatial_metric[0][0]; q[23] = spatial_metric[0][1]; q[24] = spatial_metric[0][2];
        q[25] = spatial_metric[1][0]; q[26] = spatial_metric[1][1]; q[27] = spatial_metric[1][2];
        q[28] = spatial_metric[2][0]; q[29] = spatial_metric[2][1]; q[30] = spatial_metric[2][2];

        q[31] = extrinsic_curvature[0][0]; q[32] = extrinsic_curvature[0][1]; q[33] = extrinsic_curvature[0][2];
        q[34] = extrinsic_curvature[1][0]; q[35] = extrinsic_curvature[1][1]; q[36] = extrinsic_curvature[1][2];
        q[37] = extrinsic_curvature[2][0]; q[38] = extrinsic_curvature[2][1]; q[39] = extrinsic_curvature[2][2];

        q[40] = 1.0;

        q[41] = lapse_der[0]; q[42] = lapse_der[1]; q[43] = lapse_der[2];
        q[44] = shift_der[0][0]; q[45] = shift_der[0][1]; q[46] = shift_der[0][2];
        q[47] = shift_der[1][0]; q[48] = shift_der[1][1]; q[49] = shift_der[1][2];
        q[50] = shift_der[2][0]; q[51] = shift_der[2][1]; q[52] = shift_der[2][2];

        q[53] = spatial_metric_der[0][0][0]; q[54] = spatial_metric_der[0][0][1]; q[55] = spatial_metric_der[0][0][2];
        q[56] = spatial_metric_der[0][1][0]; q[57] = spatial_metric_der[0][1][1]; q[58] = spatial_metric_der[0][1][2];
        q[59] = spatial_metric_der[0][2][0]; q[60] = spatial_metric_der[0][2][1]; q[61] = spatial_metric_der[0][2][2];

        q[62] = spatial_metric_der[1][0][0]; q[63] = spatial_metric_der[1][0][1]; q[64] = spatial_metric_der[1][0][2];
        q[65] = spatial_metric_der[1][1][0]; q[66] = spatial_metric_der[1][1][1]; q[67] = spatial_metric_der[1][1][2];
        q[68] = spatial_metric_der[1][2][0]; q[69] = spatial_metric_der[1][2][1]; q[70] = spatial_metric_der[1][2][2];

        q[71] = spatial_metric_der[2][0][0]; q[72] = spatial_metric_der[2][0][1]; q[73] = spatial_metric_der[2][0][2];
        q[74] = spatial_metric_der[2][1][0]; q[75] = spatial_metric_der[2][1][1]; q[76] = spatial_metric_der[2][1][2];
        q[77] = spatial_metric_der[2][2][0]; q[78] = spatial_metric_der[2][2][1]; q[79] = spatial_metric_der[2][2][2];

        q[80] = 0.0;
        q[81] = x; q[82] = y; q[83] = 0.0;

        double prims[84];
        gkyl_gr_twofluid_tetrad_prim_vars(gas_gamma_elc, gas_gamma_ion, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho_elc, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        TEST_CHECK( gkyl_compare(prims[5], rho_ion, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[6], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[7], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[8], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[9], p, 1e-1) );

        double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
        double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
        double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

        double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
        double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
        double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

        double fluxes[3][18] = {
          { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
            (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
            e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
            -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
          { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
            (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
            -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
            -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
          { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
            (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
            (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
            e_fact * Dz, b_fact * (light_speed * light_speed) * Bz },
        };

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

      double q_local[84], flux_local_sr[84], flux_local_gr[84], flux[84];
      for (int d = 0; d < 3; d++) {
        gr_twofluid_tetrad->rotate_to_local_func(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local_sr);
        gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
        gr_twofluid_tetrad->rotate_to_global_func(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 18; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[84], q_g[84];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 84; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }
          
          double w1[84], q1[84];
          gr_twofluid_tetrad->cons_to_riem(gr_twofluid_tetrad, q_local, q_local, w1);
          gr_twofluid_tetrad->riem_to_cons(gr_twofluid_tetrad, q_local, w1, q1);

          for (int i = 0; i < 84; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_twofluid_tetrad_waves_minkowski_ho()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_twofluid_tetrad = gkyl_wv_gr_twofluid_tetrad_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion,
    light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_elc_l = rho_ion_l * mass_elc / mass_ion;

      double rho_ion_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;
      double rho_elc_r = rho_ion_r * mass_elc / mass_ion;

      double Dx_l = 0.1, Dy_l = 0.2, Dz_l = 0.3;
      double Bx_l = 0.4, By_l = 0.5, Bz_l = 0.6;
      double phi_l = 0.0, psi_l = 0.0;

      double Dx_r = 1.0, Dy_r = 0.9, Dz_r = 0.8;
      double Bx_r = 0.7, By_r = 0.6, Bz_r = 0.5;
      double phi_r = 0.0, psi_r = 0.0;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);
      double he_l = 1.0 + ((p_l / rho_elc_l) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi_l = 1.0 + ((p_l / rho_ion_l) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));
      double he_r = 1.0 + ((p_r / rho_elc_r) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi_r = 1.0 + ((p_r / rho_ion_r) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      double ql[84], qr[84];
      ql[0] = sqrt(spatial_det_l) * rho_elc_l * W_l;
      ql[1] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * u_l;
      ql[2] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * v_l;
      ql[3] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * w_l;
      ql[4] = sqrt(spatial_det_l) * ((rho_elc_l * he_l * (W_l * W_l)) - p_l - (rho_elc_l * W_l));

      ql[5] = sqrt(spatial_det_l) * rho_ion_l * W_l;
      ql[6] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * u_l;
      ql[7] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * v_l;
      ql[8] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * w_l;
      ql[9] = sqrt(spatial_det_l) * ((rho_ion_l * hi_l * (W_l * W_l)) - p_l - (rho_ion_l * W_l));

      ql[10] = Dx_l; ql[11] = Dy_l; ql[12] = Dz_l;
      ql[13] = Bx_l; ql[14] = By_l; ql[15] = Bz_l;
      ql[16] = phi_l; ql[17] = psi_l;

      ql[18] = lapse_l;
      ql[19] = shift_l[0]; ql[20] = shift_l[1]; ql[21] = shift_l[2];

      ql[22] = spatial_metric_l[0][0]; ql[23] = spatial_metric_l[0][1]; ql[24] = spatial_metric_l[0][2];
      ql[25] = spatial_metric_l[1][0]; ql[26] = spatial_metric_l[1][1]; ql[27] = spatial_metric_l[1][2];
      ql[28] = spatial_metric_l[2][0]; ql[29] = spatial_metric_l[2][1]; ql[30] = spatial_metric_l[2][2];

      ql[31] = extrinsic_curvature_l[0][0]; ql[32] = extrinsic_curvature_l[0][1]; ql[33] = extrinsic_curvature_l[0][2];
      ql[34] = extrinsic_curvature_l[1][0]; ql[35] = extrinsic_curvature_l[1][1]; ql[36] = extrinsic_curvature_l[1][2];
      ql[37] = extrinsic_curvature_l[2][0]; ql[38] = extrinsic_curvature_l[2][1]; ql[39] = extrinsic_curvature_l[2][2];

      ql[40] = 1.0;

      ql[41] = lapse_der_l[0]; ql[42] = lapse_der_l[1]; ql[43] = lapse_der_l[2];
      ql[44] = shift_der_l[0][0]; ql[45] = shift_der_l[0][1]; ql[46] = shift_der_l[0][2];
      ql[47] = shift_der_l[1][0]; ql[48] = shift_der_l[1][1]; ql[49] = shift_der_l[1][2];
      ql[50] = shift_der_l[2][0]; ql[51] = shift_der_l[2][1]; ql[52] = shift_der_l[2][2];

      ql[53] = spatial_metric_der_l[0][0][0]; ql[54] = spatial_metric_der_l[0][0][1]; ql[55] = spatial_metric_der_l[0][0][2];
      ql[56] = spatial_metric_der_l[0][1][0]; ql[57] = spatial_metric_der_l[0][1][1]; ql[58] = spatial_metric_der_l[0][1][2];
      ql[59] = spatial_metric_der_l[0][2][0]; ql[60] = spatial_metric_der_l[0][2][1]; ql[61] = spatial_metric_der_l[0][2][2];

      ql[62] = spatial_metric_der_l[1][0][0]; ql[63] = spatial_metric_der_l[1][0][1]; ql[64] = spatial_metric_der_l[1][0][2];
      ql[65] = spatial_metric_der_l[1][1][0]; ql[66] = spatial_metric_der_l[1][1][1]; ql[67] = spatial_metric_der_l[1][1][2];
      ql[68] = spatial_metric_der_l[1][2][0]; ql[69] = spatial_metric_der_l[1][2][1]; ql[70] = spatial_metric_der_l[1][2][2];

      ql[71] = spatial_metric_der_l[2][0][0]; ql[72] = spatial_metric_der_l[2][0][1]; ql[73] = spatial_metric_der_l[2][0][2];
      ql[74] = spatial_metric_der_l[2][1][0]; ql[75] = spatial_metric_der_l[2][1][1]; ql[76] = spatial_metric_der_l[2][1][2];
      ql[77] = spatial_metric_der_l[2][2][0]; ql[78] = spatial_metric_der_l[2][2][1]; ql[79] = spatial_metric_der_l[2][2][2];

      ql[80] = 0.0;
      ql[81] = x - 0.5; ql[82] = y; ql[83] = 0.0;

      qr[0] = sqrt(spatial_det_r) * rho_elc_r * W_r;
      qr[1] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * u_r;
      qr[2] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * v_r;
      qr[3] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * w_r;
      qr[4] = sqrt(spatial_det_r) * ((rho_elc_r * he_r * (W_r * W_r)) - p_r - (rho_elc_r * W_r));

      qr[5] = sqrt(spatial_det_r) * rho_ion_r * W_r;
      qr[6] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * u_r;
      qr[7] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * v_r;
      qr[8] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * w_r;
      qr[9] = sqrt(spatial_det_r) * ((rho_ion_r * hi_r * (W_r * W_r)) - p_r - (rho_ion_r * W_r));

      qr[10] = Dx_r; qr[11] = Dy_r; qr[12] = Dz_r;
      qr[13] = Bx_r; qr[14] = By_r; qr[15] = Bz_r;
      qr[16] = phi_r; qr[17] = psi_r;

      qr[18] = lapse_r;
      qr[19] = shift_r[0]; qr[20] = shift_r[1]; qr[21] = shift_r[2];

      qr[22] = spatial_metric_r[0][0]; qr[23] = spatial_metric_r[0][1]; qr[24] = spatial_metric_r[0][2];
      qr[25] = spatial_metric_r[1][0]; qr[26] = spatial_metric_r[1][1]; qr[27] = spatial_metric_r[1][2];
      qr[28] = spatial_metric_r[2][0]; qr[29] = spatial_metric_r[2][1]; qr[30] = spatial_metric_r[2][2];

      qr[31] = extrinsic_curvature_r[0][0]; qr[32] = extrinsic_curvature_r[0][1]; qr[33] = extrinsic_curvature_r[0][2];
      qr[34] = extrinsic_curvature_r[1][0]; qr[35] = extrinsic_curvature_r[1][1]; qr[36] = extrinsic_curvature_r[1][2];
      qr[37] = extrinsic_curvature_r[2][0]; qr[38] = extrinsic_curvature_r[2][1]; qr[39] = extrinsic_curvature_r[2][2];

      qr[40] = 1.0;

      qr[41] = lapse_der_r[0]; qr[42] = lapse_der_r[1]; qr[43] = lapse_der_r[2];
      qr[44] = shift_der_r[0][0]; qr[45] = shift_der_r[0][1]; qr[46] = shift_der_r[0][2];
      qr[47] = shift_der_r[1][0]; qr[48] = shift_der_r[1][1]; qr[49] = shift_der_r[1][2];
      qr[50] = shift_der_r[2][0]; qr[51] = shift_der_r[2][1]; qr[52] = shift_der_r[2][2];

      qr[53] = spatial_metric_der_r[0][0][0]; qr[54] = spatial_metric_der_r[0][0][1]; qr[55] = spatial_metric_der_r[0][0][2];
      qr[56] = spatial_metric_der_r[0][1][0]; qr[57] = spatial_metric_der_r[0][1][1]; qr[58] = spatial_metric_der_r[0][1][2];
      qr[59] = spatial_metric_der_r[0][2][0]; qr[60] = spatial_metric_der_r[0][2][1]; qr[61] = spatial_metric_der_r[0][2][2];

      qr[62] = spatial_metric_der_r[1][0][0]; qr[63] = spatial_metric_der_r[1][0][1]; qr[64] = spatial_metric_der_r[1][0][2];
      qr[65] = spatial_metric_der_r[1][1][0]; qr[66] = spatial_metric_der_r[1][1][1]; qr[67] = spatial_metric_der_r[1][1][2];
      qr[68] = spatial_metric_der_r[1][2][0]; qr[69] = spatial_metric_der_r[1][2][1]; qr[70] = spatial_metric_der_r[1][2][2];

      qr[71] = spatial_metric_der_r[2][0][0]; qr[72] = spatial_metric_der_r[2][0][1]; qr[73] = spatial_metric_der_r[2][0][2];
      qr[74] = spatial_metric_der_r[2][1][0]; qr[75] = spatial_metric_der_r[2][1][1]; qr[76] = spatial_metric_der_r[2][1][2];
      qr[77] = spatial_metric_der_r[2][2][0]; qr[78] = spatial_metric_der_r[2][2][1]; qr[79] = spatial_metric_der_r[2][2][2];

      qr[80] = 0.0;
      qr[81] = x + 0.5; qr[82] = y; qr[83] = 0.0;

      double norm[3][3] = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
      };

      double tau1[3][3] = {
        { 0.0, 1.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
      };

      double tau2[3][3] = {
        { 0.0, 0.0, 1.0 },
        { 0.0, 0.0, -1.0 },
        { 0.0, 1.0, 0.0 },
      };

      for (int d = 0; d < 3; d++) {
        double speeds[3], waves[3 * 84], waves_local[3 * 84];

        double ql_local[84], qr_local[84];
        gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
        gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

        double delta[84];
        for (int i = 0; i < 84; i++) {
          delta[i] = qr_local[i] - ql_local[i];
        }

        gkyl_wv_eqn_waves(gr_twofluid_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

        double apdq_local[84], amdq_local[84];
        gkyl_wv_eqn_qfluct(gr_twofluid_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

        for (int i = 0; i < 3; i++) {
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 84], &waves[i * 84]);
        }

        double apdq[84], amdq[84];
        gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
        gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

        double fl_local_sr[84], fr_local_sr[84];
        gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql_local, fl_local_sr);
        gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr_local, fr_local_sr);

        double fl_local_gr[84], fr_local_gr[84];
        gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql_local, fl_local_sr, fl_local_gr);
        gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr_local, fr_local_sr, fr_local_gr);

        double fl[84], fr[84];
        gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
        gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

        for (int i = 0; i < 84; i++) {
          TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_twofluid_tetrad_waves_schwarzschild_ho()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_twofluid_tetrad = gkyl_wv_gr_twofluid_tetrad_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion,
    light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_elc_l = rho_ion_l * mass_elc / mass_ion;

      double rho_ion_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;
      double rho_elc_r = rho_ion_r * mass_elc / mass_ion;

      double Dx_l = 0.1, Dy_l = 0.2, Dz_l = 0.3;
      double Bx_l = 0.4, By_l = 0.5, Bz_l = 0.6;
      double phi_l = 0.0, psi_l = 0.0;

      double Dx_r = 1.0, Dy_r = 0.9, Dz_r = 0.8;
      double Bx_r = 0.7, By_r = 0.6, Bz_r = 0.5;
      double phi_r = 0.0, psi_r = 0.0;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);
      double he_l = 1.0 + ((p_l / rho_elc_l) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi_l = 1.0 + ((p_l / rho_ion_l) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));
      double he_r = 1.0 + ((p_r / rho_elc_r) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi_r = 1.0 + ((p_r / rho_ion_r) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[84], qr[84];
        ql[0] = sqrt(spatial_det_l) * rho_elc_l * W_l;
        ql[1] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * u_l;
        ql[2] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * v_l;
        ql[3] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * w_l;
        ql[4] = sqrt(spatial_det_l) * ((rho_elc_l * he_l * (W_l * W_l)) - p_l - (rho_elc_l * W_l));

        ql[5] = sqrt(spatial_det_l) * rho_ion_l * W_l;
        ql[6] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * u_l;
        ql[7] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * v_l;
        ql[8] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * w_l;
        ql[9] = sqrt(spatial_det_l) * ((rho_ion_l * hi_l * (W_l * W_l)) - p_l - (rho_ion_l * W_l));

        ql[10] = Dx_l; ql[11] = Dy_l; ql[12] = Dz_l;
        ql[13] = Bx_l; ql[14] = By_l; ql[15] = Bz_l;
        ql[16] = phi_l; ql[17] = psi_l;

        ql[18] = lapse_l;
        ql[19] = shift_l[0]; ql[20] = shift_l[1]; ql[21] = shift_l[2];

        ql[22] = spatial_metric_l[0][0]; ql[23] = spatial_metric_l[0][1]; ql[24] = spatial_metric_l[0][2];
        ql[25] = spatial_metric_l[1][0]; ql[26] = spatial_metric_l[1][1]; ql[27] = spatial_metric_l[1][2];
        ql[28] = spatial_metric_l[2][0]; ql[29] = spatial_metric_l[2][1]; ql[30] = spatial_metric_l[2][2];

        ql[31] = extrinsic_curvature_l[0][0]; ql[32] = extrinsic_curvature_l[0][1]; ql[33] = extrinsic_curvature_l[0][2];
        ql[34] = extrinsic_curvature_l[1][0]; ql[35] = extrinsic_curvature_l[1][1]; ql[36] = extrinsic_curvature_l[1][2];
        ql[37] = extrinsic_curvature_l[2][0]; ql[38] = extrinsic_curvature_l[2][1]; ql[39] = extrinsic_curvature_l[2][2];

        ql[40] = 1.0;

        ql[41] = lapse_der_l[0]; ql[42] = lapse_der_l[1]; ql[43] = lapse_der_l[2];
        ql[44] = shift_der_l[0][0]; ql[45] = shift_der_l[0][1]; ql[46] = shift_der_l[0][2];
        ql[47] = shift_der_l[1][0]; ql[48] = shift_der_l[1][1]; ql[49] = shift_der_l[1][2];
        ql[50] = shift_der_l[2][0]; ql[51] = shift_der_l[2][1]; ql[52] = shift_der_l[2][2];

        ql[53] = spatial_metric_der_l[0][0][0]; ql[54] = spatial_metric_der_l[0][0][1]; ql[55] = spatial_metric_der_l[0][0][2];
        ql[56] = spatial_metric_der_l[0][1][0]; ql[57] = spatial_metric_der_l[0][1][1]; ql[58] = spatial_metric_der_l[0][1][2];
        ql[59] = spatial_metric_der_l[0][2][0]; ql[60] = spatial_metric_der_l[0][2][1]; ql[61] = spatial_metric_der_l[0][2][2];

        ql[62] = spatial_metric_der_l[1][0][0]; ql[63] = spatial_metric_der_l[1][0][1]; ql[64] = spatial_metric_der_l[1][0][2];
        ql[65] = spatial_metric_der_l[1][1][0]; ql[66] = spatial_metric_der_l[1][1][1]; ql[67] = spatial_metric_der_l[1][1][2];
        ql[68] = spatial_metric_der_l[1][2][0]; ql[69] = spatial_metric_der_l[1][2][1]; ql[70] = spatial_metric_der_l[1][2][2];

        ql[71] = spatial_metric_der_l[2][0][0]; ql[72] = spatial_metric_der_l[2][0][1]; ql[73] = spatial_metric_der_l[2][0][2];
        ql[74] = spatial_metric_der_l[2][1][0]; ql[75] = spatial_metric_der_l[2][1][1]; ql[76] = spatial_metric_der_l[2][1][2];
        ql[77] = spatial_metric_der_l[2][2][0]; ql[78] = spatial_metric_der_l[2][2][1]; ql[79] = spatial_metric_der_l[2][2][2];

        ql[80] = 0.0;
        ql[81] = x - 0.5; ql[82] = y; ql[83] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_elc_r * W_r;
        qr[1] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * u_r;
        qr[2] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * v_r;
        qr[3] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * w_r;
        qr[4] = sqrt(spatial_det_r) * ((rho_elc_r * he_r * (W_r * W_r)) - p_r - (rho_elc_r * W_r));

        qr[5] = sqrt(spatial_det_r) * rho_ion_r * W_r;
        qr[6] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * u_r;
        qr[7] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * v_r;
        qr[8] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * w_r;
        qr[9] = sqrt(spatial_det_r) * ((rho_ion_r * hi_r * (W_r * W_r)) - p_r - (rho_ion_r * W_r));

        qr[10] = Dx_r; qr[11] = Dy_r; qr[12] = Dz_r;
        qr[13] = Bx_r; qr[14] = By_r; qr[15] = Bz_r;
        qr[16] = phi_r; qr[17] = psi_r;

        qr[18] = lapse_r;
        qr[19] = shift_r[0]; qr[20] = shift_r[1]; qr[21] = shift_r[2];

        qr[22] = spatial_metric_r[0][0]; qr[23] = spatial_metric_r[0][1]; qr[24] = spatial_metric_r[0][2];
        qr[25] = spatial_metric_r[1][0]; qr[26] = spatial_metric_r[1][1]; qr[27] = spatial_metric_r[1][2];
        qr[28] = spatial_metric_r[2][0]; qr[29] = spatial_metric_r[2][1]; qr[30] = spatial_metric_r[2][2];

        qr[31] = extrinsic_curvature_r[0][0]; qr[32] = extrinsic_curvature_r[0][1]; qr[33] = extrinsic_curvature_r[0][2];
        qr[34] = extrinsic_curvature_r[1][0]; qr[35] = extrinsic_curvature_r[1][1]; qr[36] = extrinsic_curvature_r[1][2];
        qr[37] = extrinsic_curvature_r[2][0]; qr[38] = extrinsic_curvature_r[2][1]; qr[39] = extrinsic_curvature_r[2][2];

        qr[40] = 1.0;

        qr[41] = lapse_der_r[0]; qr[42] = lapse_der_r[1]; qr[43] = lapse_der_r[2];
        qr[44] = shift_der_r[0][0]; qr[45] = shift_der_r[0][1]; qr[46] = shift_der_r[0][2];
        qr[47] = shift_der_r[1][0]; qr[48] = shift_der_r[1][1]; qr[49] = shift_der_r[1][2];
        qr[50] = shift_der_r[2][0]; qr[51] = shift_der_r[2][1]; qr[52] = shift_der_r[2][2];

        qr[53] = spatial_metric_der_r[0][0][0]; qr[54] = spatial_metric_der_r[0][0][1]; qr[55] = spatial_metric_der_r[0][0][2];
        qr[56] = spatial_metric_der_r[0][1][0]; qr[57] = spatial_metric_der_r[0][1][1]; qr[58] = spatial_metric_der_r[0][1][2];
        qr[59] = spatial_metric_der_r[0][2][0]; qr[60] = spatial_metric_der_r[0][2][1]; qr[61] = spatial_metric_der_r[0][2][2];

        qr[62] = spatial_metric_der_r[1][0][0]; qr[63] = spatial_metric_der_r[1][0][1]; qr[64] = spatial_metric_der_r[1][0][2];
        qr[65] = spatial_metric_der_r[1][1][0]; qr[66] = spatial_metric_der_r[1][1][1]; qr[67] = spatial_metric_der_r[1][1][2];
        qr[68] = spatial_metric_der_r[1][2][0]; qr[69] = spatial_metric_der_r[1][2][1]; qr[70] = spatial_metric_der_r[1][2][2];

        qr[71] = spatial_metric_der_r[2][0][0]; qr[72] = spatial_metric_der_r[2][0][1]; qr[73] = spatial_metric_der_r[2][0][2];
        qr[74] = spatial_metric_der_r[2][1][0]; qr[75] = spatial_metric_der_r[2][1][1]; qr[76] = spatial_metric_der_r[2][1][2];
        qr[77] = spatial_metric_der_r[2][2][0]; qr[78] = spatial_metric_der_r[2][2][1]; qr[79] = spatial_metric_der_r[2][2][2];

        qr[80] = 0.0;
        qr[81] = x + 0.5; qr[82] = y; qr[83] = 0.0;

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

        for (int d = 0; d < 3; d++) {
          double speeds[3], waves[3 * 84], waves_local[3 * 84];

          double ql_local[84], qr_local[84];
          gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[84];
          for (int i = 0; i < 84; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_twofluid_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[84], amdq_local[84];
          gkyl_wv_eqn_qfluct(gr_twofluid_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 3; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 84], &waves[i * 84]);
          }

          double apdq[84], amdq[84];
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local_sr[84], fr_local_sr[84];
          gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql_local, fl_local_sr);
          gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr_local, fr_local_sr);

          double fl_local_gr[84], fr_local_gr[84];
          gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql_local, fl_local_sr, fl_local_gr);
          gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr_local, fr_local_sr, fr_local_gr);

          double fl[84], fr[84];
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

          for (int i = 0; i < 84; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_twofluid_tetrad_waves_kerr_ho()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_twofluid_tetrad = gkyl_wv_gr_twofluid_tetrad_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion,
    light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_elc_l = rho_ion_l * mass_elc / mass_ion;

      double rho_ion_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;
      double rho_elc_r = rho_ion_r * mass_elc / mass_ion;

      double Dx_l = 0.1, Dy_l = 0.2, Dz_l = 0.3;
      double Bx_l = 0.4, By_l = 0.5, Bz_l = 0.6;
      double phi_l = 0.0, psi_l = 0.0;

      double Dx_r = 1.0, Dy_r = 0.9, Dz_r = 0.8;
      double Bx_r = 0.7, By_r = 0.6, Bz_r = 0.5;
      double phi_r = 0.0, psi_r = 0.0;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);
      double he_l = 1.0 + ((p_l / rho_elc_l) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi_l = 1.0 + ((p_l / rho_ion_l) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));
      double he_r = 1.0 + ((p_r / rho_elc_r) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi_r = 1.0 + ((p_r / rho_ion_r) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[84], qr[84];
        ql[0] = sqrt(spatial_det_l) * rho_elc_l * W_l;
        ql[1] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * u_l;
        ql[2] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * v_l;
        ql[3] = sqrt(spatial_det_l) * rho_elc_l * he_l * (W_l * W_l) * w_l;
        ql[4] = sqrt(spatial_det_l) * ((rho_elc_l * he_l * (W_l * W_l)) - p_l - (rho_elc_l * W_l));

        ql[5] = sqrt(spatial_det_l) * rho_ion_l * W_l;
        ql[6] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * u_l;
        ql[7] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * v_l;
        ql[8] = sqrt(spatial_det_l) * rho_ion_l * hi_l * (W_l * W_l) * w_l;
        ql[9] = sqrt(spatial_det_l) * ((rho_ion_l * hi_l * (W_l * W_l)) - p_l - (rho_ion_l * W_l));

        ql[10] = Dx_l; ql[11] = Dy_l; ql[12] = Dz_l;
        ql[13] = Bx_l; ql[14] = By_l; ql[15] = Bz_l;
        ql[16] = phi_l; ql[17] = psi_l;

        ql[18] = lapse_l;
        ql[19] = shift_l[0]; ql[20] = shift_l[1]; ql[21] = shift_l[2];

        ql[22] = spatial_metric_l[0][0]; ql[23] = spatial_metric_l[0][1]; ql[24] = spatial_metric_l[0][2];
        ql[25] = spatial_metric_l[1][0]; ql[26] = spatial_metric_l[1][1]; ql[27] = spatial_metric_l[1][2];
        ql[28] = spatial_metric_l[2][0]; ql[29] = spatial_metric_l[2][1]; ql[30] = spatial_metric_l[2][2];

        ql[31] = extrinsic_curvature_l[0][0]; ql[32] = extrinsic_curvature_l[0][1]; ql[33] = extrinsic_curvature_l[0][2];
        ql[34] = extrinsic_curvature_l[1][0]; ql[35] = extrinsic_curvature_l[1][1]; ql[36] = extrinsic_curvature_l[1][2];
        ql[37] = extrinsic_curvature_l[2][0]; ql[38] = extrinsic_curvature_l[2][1]; ql[39] = extrinsic_curvature_l[2][2];

        ql[40] = 1.0;

        ql[41] = lapse_der_l[0]; ql[42] = lapse_der_l[1]; ql[43] = lapse_der_l[2];
        ql[44] = shift_der_l[0][0]; ql[45] = shift_der_l[0][1]; ql[46] = shift_der_l[0][2];
        ql[47] = shift_der_l[1][0]; ql[48] = shift_der_l[1][1]; ql[49] = shift_der_l[1][2];
        ql[50] = shift_der_l[2][0]; ql[51] = shift_der_l[2][1]; ql[52] = shift_der_l[2][2];

        ql[53] = spatial_metric_der_l[0][0][0]; ql[54] = spatial_metric_der_l[0][0][1]; ql[55] = spatial_metric_der_l[0][0][2];
        ql[56] = spatial_metric_der_l[0][1][0]; ql[57] = spatial_metric_der_l[0][1][1]; ql[58] = spatial_metric_der_l[0][1][2];
        ql[59] = spatial_metric_der_l[0][2][0]; ql[60] = spatial_metric_der_l[0][2][1]; ql[61] = spatial_metric_der_l[0][2][2];

        ql[62] = spatial_metric_der_l[1][0][0]; ql[63] = spatial_metric_der_l[1][0][1]; ql[64] = spatial_metric_der_l[1][0][2];
        ql[65] = spatial_metric_der_l[1][1][0]; ql[66] = spatial_metric_der_l[1][1][1]; ql[67] = spatial_metric_der_l[1][1][2];
        ql[68] = spatial_metric_der_l[1][2][0]; ql[69] = spatial_metric_der_l[1][2][1]; ql[70] = spatial_metric_der_l[1][2][2];

        ql[71] = spatial_metric_der_l[2][0][0]; ql[72] = spatial_metric_der_l[2][0][1]; ql[73] = spatial_metric_der_l[2][0][2];
        ql[74] = spatial_metric_der_l[2][1][0]; ql[75] = spatial_metric_der_l[2][1][1]; ql[76] = spatial_metric_der_l[2][1][2];
        ql[77] = spatial_metric_der_l[2][2][0]; ql[78] = spatial_metric_der_l[2][2][1]; ql[79] = spatial_metric_der_l[2][2][2];

        ql[80] = 0.0;
        ql[81] = x - 0.5; ql[82] = y; ql[83] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_elc_r * W_r;
        qr[1] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * u_r;
        qr[2] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * v_r;
        qr[3] = sqrt(spatial_det_r) * rho_elc_r * he_r * (W_r * W_r) * w_r;
        qr[4] = sqrt(spatial_det_r) * ((rho_elc_r * he_r * (W_r * W_r)) - p_r - (rho_elc_r * W_r));

        qr[5] = sqrt(spatial_det_r) * rho_ion_r * W_r;
        qr[6] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * u_r;
        qr[7] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * v_r;
        qr[8] = sqrt(spatial_det_r) * rho_ion_r * hi_r * (W_r * W_r) * w_r;
        qr[9] = sqrt(spatial_det_r) * ((rho_ion_r * hi_r * (W_r * W_r)) - p_r - (rho_ion_r * W_r));

        qr[10] = Dx_r; qr[11] = Dy_r; qr[12] = Dz_r;
        qr[13] = Bx_r; qr[14] = By_r; qr[15] = Bz_r;
        qr[16] = phi_r; qr[17] = psi_r;

        qr[18] = lapse_r;
        qr[19] = shift_r[0]; qr[20] = shift_r[1]; qr[21] = shift_r[2];

        qr[22] = spatial_metric_r[0][0]; qr[23] = spatial_metric_r[0][1]; qr[24] = spatial_metric_r[0][2];
        qr[25] = spatial_metric_r[1][0]; qr[26] = spatial_metric_r[1][1]; qr[27] = spatial_metric_r[1][2];
        qr[28] = spatial_metric_r[2][0]; qr[29] = spatial_metric_r[2][1]; qr[30] = spatial_metric_r[2][2];

        qr[31] = extrinsic_curvature_r[0][0]; qr[32] = extrinsic_curvature_r[0][1]; qr[33] = extrinsic_curvature_r[0][2];
        qr[34] = extrinsic_curvature_r[1][0]; qr[35] = extrinsic_curvature_r[1][1]; qr[36] = extrinsic_curvature_r[1][2];
        qr[37] = extrinsic_curvature_r[2][0]; qr[38] = extrinsic_curvature_r[2][1]; qr[39] = extrinsic_curvature_r[2][2];

        qr[40] = 1.0;

        qr[41] = lapse_der_r[0]; qr[42] = lapse_der_r[1]; qr[43] = lapse_der_r[2];
        qr[44] = shift_der_r[0][0]; qr[45] = shift_der_r[0][1]; qr[46] = shift_der_r[0][2];
        qr[47] = shift_der_r[1][0]; qr[48] = shift_der_r[1][1]; qr[49] = shift_der_r[1][2];
        qr[50] = shift_der_r[2][0]; qr[51] = shift_der_r[2][1]; qr[52] = shift_der_r[2][2];

        qr[53] = spatial_metric_der_r[0][0][0]; qr[54] = spatial_metric_der_r[0][0][1]; qr[55] = spatial_metric_der_r[0][0][2];
        qr[56] = spatial_metric_der_r[0][1][0]; qr[57] = spatial_metric_der_r[0][1][1]; qr[58] = spatial_metric_der_r[0][1][2];
        qr[59] = spatial_metric_der_r[0][2][0]; qr[60] = spatial_metric_der_r[0][2][1]; qr[61] = spatial_metric_der_r[0][2][2];

        qr[62] = spatial_metric_der_r[1][0][0]; qr[63] = spatial_metric_der_r[1][0][1]; qr[64] = spatial_metric_der_r[1][0][2];
        qr[65] = spatial_metric_der_r[1][1][0]; qr[66] = spatial_metric_der_r[1][1][1]; qr[67] = spatial_metric_der_r[1][1][2];
        qr[68] = spatial_metric_der_r[1][2][0]; qr[69] = spatial_metric_der_r[1][2][1]; qr[70] = spatial_metric_der_r[1][2][2];

        qr[71] = spatial_metric_der_r[2][0][0]; qr[72] = spatial_metric_der_r[2][0][1]; qr[73] = spatial_metric_der_r[2][0][2];
        qr[74] = spatial_metric_der_r[2][1][0]; qr[75] = spatial_metric_der_r[2][1][1]; qr[76] = spatial_metric_der_r[2][1][2];
        qr[77] = spatial_metric_der_r[2][2][0]; qr[78] = spatial_metric_der_r[2][2][1]; qr[79] = spatial_metric_der_r[2][2][2];

        qr[80] = 0.0;
        qr[81] = x + 0.5; qr[82] = y; qr[83] = 0.0;

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

        for (int d = 0; d < 3; d++) {
          double speeds[3], waves[3 * 84], waves_local[3 * 84];

          double ql_local[84], qr_local[84];
          gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[84];
          for (int i = 0; i < 84; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_twofluid_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[84], amdq_local[84];
          gkyl_wv_eqn_qfluct(gr_twofluid_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 3; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 84], &waves[i * 84]);
          }

          double apdq[84], amdq[84];
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local_sr[84], fr_local_sr[84];
          gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql_local, fl_local_sr);
          gkyl_gr_twofluid_tetrad_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr_local, fr_local_sr);

          double fl_local_gr[84], fr_local_gr[84];
          gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, ql_local, fl_local_sr, fl_local_gr);
          gkyl_gr_twofluid_tetrad_flux_correction(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, qr_local, fr_local_sr, fr_local_gr);

          double fl[84], fr[84];
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
          gkyl_wv_eqn_rotate_to_global(gr_twofluid_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

          for (int i = 0; i < 84; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_twofluid_tetrad_basic_minkowski_ho", test_gr_twofluid_tetrad_basic_minkowski_ho },
  { "gr_twofluid_tetrad_basic_schwarzschild_ho", test_gr_twofluid_tetrad_basic_schwarzschild_ho },
  { "gr_twofluid_tetrad_basic_kerr_ho", test_gr_twofluid_tetrad_basic_kerr_ho },
  { "gr_twofluid_tetrad_waves_minkowski_ho", test_gr_twofluid_tetrad_waves_minkowski_ho },
  { "gr_twofluid_tetrad_waves_schwarzschild_ho", test_gr_twofluid_tetrad_waves_schwarzschild_ho },
  { "gr_twofluid_tetrad_waves_kerr_ho", test_gr_twofluid_tetrad_waves_kerr_ho },
  { NULL, NULL },
};
