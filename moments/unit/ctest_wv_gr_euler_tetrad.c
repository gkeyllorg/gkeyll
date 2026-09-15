#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_euler_tetrad_basic_minkowski_ho()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_euler_tetrad = gkyl_wv_gr_euler_tetrad_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_euler_tetrad->num_equations == 71 );
  TEST_CHECK( gr_euler_tetrad->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;

      double spatial_det = 0.0, lapse = 0.0;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region = false;

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
      double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

      double q[71];
      q[0] = sqrt(spatial_det) * rho * W;
      q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
      q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
      q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
      q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

      q[5] = lapse;
      q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];

      q[9] = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
      q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
      q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];

      q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
      q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
      q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];

      q[27] = 1.0;

      q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
      q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
      q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
      q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];

      q[40] = spatial_metric_der[0][0][0]; q[41] = spatial_metric_der[0][0][1]; q[42] = spatial_metric_der[0][0][2];
      q[43] = spatial_metric_der[0][1][0]; q[44] = spatial_metric_der[0][1][1]; q[45] = spatial_metric_der[0][1][2];
      q[46] = spatial_metric_der[0][2][0]; q[47] = spatial_metric_der[0][2][1]; q[48] = spatial_metric_der[0][2][2];

      q[49] = spatial_metric_der[1][0][0]; q[50] = spatial_metric_der[1][0][1]; q[51] = spatial_metric_der[1][0][2];
      q[52] = spatial_metric_der[1][1][0]; q[53] = spatial_metric_der[1][1][1]; q[54] = spatial_metric_der[1][1][2];
      q[55] = spatial_metric_der[1][2][0]; q[56] = spatial_metric_der[1][2][1]; q[57] = spatial_metric_der[1][2][2];

      q[58] = spatial_metric_der[2][0][0]; q[59] = spatial_metric_der[2][0][1]; q[60] = spatial_metric_der[2][0][2];
      q[61] = spatial_metric_der[2][1][0]; q[62] = spatial_metric_der[2][1][1]; q[63] = spatial_metric_der[2][1][2];
      q[64] = spatial_metric_der[2][2][0]; q[65] = spatial_metric_der[2][2][1]; q[66] = spatial_metric_der[2][2][2];

      q[67] = 0.0;
      q[68] = x; q[69] = y; q[70] = 0.0;

      double prims[71] = {0};
      gkyl_gr_euler_tetrad_prim_vars(gas_gamma, q, prims);
      
      TEST_CHECK( gkyl_compare(prims[0], rho, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[1], u, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[2], v, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[3], w, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[4], p, 1e-13) );

      double fluxes[3][5] = {
        { (lapse * sqrt(spatial_det)) * (rho * W * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])) },
        { (lapse * sqrt(spatial_det)) * (rho * W * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])) },
        { (lapse * sqrt(spatial_det)) * (rho * W * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])) },
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

      double q_local[71] = {0}, flux_local_sr[71] = {0}, flux_local_gr[71] = {0}, flux[71] = {0};
      for (int d = 0; d < 3; d++) {
        gr_euler_tetrad->rotate_to_local_func(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_euler_tetrad_flux(gas_gamma, q_local, flux_local_sr);
        gkyl_gr_euler_tetrad_flux_correction(gas_gamma, q_local, flux_local_sr, flux_local_gr);
        gr_euler_tetrad->rotate_to_global_func(gr_euler_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

        for (int i = 0; i < 5; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-8) );
        }
      }

      double q_l[71] = {0}, q_g[71] = {0};
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 71; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
        }

        double w1[71] = {0}, q1[71] = {0};
        gr_euler_tetrad->cons_to_riem(gr_euler_tetrad, q_local, q_local, w1);
        gr_euler_tetrad->riem_to_cons(gr_euler_tetrad, q_local, w1, q1);

        for (int i = 0; i < 71; i++) {
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

  gkyl_wv_eqn_release(gr_euler_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_basic_schwarzschild_ho()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler_tetrad = gkyl_wv_gr_euler_tetrad_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_euler_tetrad->num_equations == 71 );
  TEST_CHECK( gr_euler_tetrad->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;

      double spatial_det = 0.0, lapse = 0.0;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region = false;

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
      double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region) {
        double q[71];
        q[0] = sqrt(spatial_det) * rho * W;
        q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
        q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
        q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
        q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

        q[5] = lapse;
        q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];

        q[9] = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
        q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
        q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];

        q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
        q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
        q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];

        q[27] = 1.0;

        q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
        q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
        q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
        q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];

        q[40] = spatial_metric_der[0][0][0]; q[41] = spatial_metric_der[0][0][1]; q[42] = spatial_metric_der[0][0][2];
        q[43] = spatial_metric_der[0][1][0]; q[44] = spatial_metric_der[0][1][1]; q[45] = spatial_metric_der[0][1][2];
        q[46] = spatial_metric_der[0][2][0]; q[47] = spatial_metric_der[0][2][1]; q[48] = spatial_metric_der[0][2][2];

        q[49] = spatial_metric_der[1][0][0]; q[50] = spatial_metric_der[1][0][1]; q[51] = spatial_metric_der[1][0][2];
        q[52] = spatial_metric_der[1][1][0]; q[53] = spatial_metric_der[1][1][1]; q[54] = spatial_metric_der[1][1][2];
        q[55] = spatial_metric_der[1][2][0]; q[56] = spatial_metric_der[1][2][1]; q[57] = spatial_metric_der[1][2][2];

        q[58] = spatial_metric_der[2][0][0]; q[59] = spatial_metric_der[2][0][1]; q[60] = spatial_metric_der[2][0][2];
        q[61] = spatial_metric_der[2][1][0]; q[62] = spatial_metric_der[2][1][1]; q[63] = spatial_metric_der[2][1][2];
        q[64] = spatial_metric_der[2][2][0]; q[65] = spatial_metric_der[2][2][1]; q[66] = spatial_metric_der[2][2][2];

        q[67] = 0.0;
        q[68] = x; q[69] = y; q[70] = 0.0;

        double prims[71] = {0};
        gkyl_gr_euler_tetrad_prim_vars(gas_gamma, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        double fluxes[3][5] = {
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])) },
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

        double q_local[71] = {0}, flux_local_sr[71] = {0}, flux_local_gr[71] = {0}, flux[71] = {0};
        for (int d = 0; d < 3; d++) {
          gr_euler_tetrad->rotate_to_local_func(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_euler_tetrad_flux(gas_gamma, q_local, flux_local_sr);
          gkyl_gr_euler_tetrad_flux_correction(gas_gamma, q_local, flux_local_sr, flux_local_gr);
          gr_euler_tetrad->rotate_to_global_func(gr_euler_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 5; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[71] = {0}, q_g[71] = {0};
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[71] = {0}, q1[71] = {0};
          gr_euler_tetrad->cons_to_riem(gr_euler_tetrad, q_local, q_local, w1);
          gr_euler_tetrad->riem_to_cons(gr_euler_tetrad, q_local, w1, q1);

          for (int i = 0; i < 71; i++) {
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

  gkyl_wv_eqn_release(gr_euler_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_basic_kerr_ho()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler_tetrad = gkyl_wv_gr_euler_tetrad_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_euler_tetrad->num_equations == 71 );
  TEST_CHECK( gr_euler_tetrad->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;

      double spatial_det = 0.0, lapse = 0.0;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region = false;

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
      double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region) {
        double q[71];
        q[0] = sqrt(spatial_det) * rho * W;
        q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
        q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
        q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
        q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

        q[5] = lapse;
        q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];

        q[9] = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
        q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
        q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];

        q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
        q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
        q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];

        q[27] = 1.0;

        q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
        q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
        q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
        q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];

        q[40] = spatial_metric_der[0][0][0]; q[41] = spatial_metric_der[0][0][1]; q[42] = spatial_metric_der[0][0][2];
        q[43] = spatial_metric_der[0][1][0]; q[44] = spatial_metric_der[0][1][1]; q[45] = spatial_metric_der[0][1][2];
        q[46] = spatial_metric_der[0][2][0]; q[47] = spatial_metric_der[0][2][1]; q[48] = spatial_metric_der[0][2][2];

        q[49] = spatial_metric_der[1][0][0]; q[50] = spatial_metric_der[1][0][1]; q[51] = spatial_metric_der[1][0][2];
        q[52] = spatial_metric_der[1][1][0]; q[53] = spatial_metric_der[1][1][1]; q[54] = spatial_metric_der[1][1][2];
        q[55] = spatial_metric_der[1][2][0]; q[56] = spatial_metric_der[1][2][1]; q[57] = spatial_metric_der[1][2][2];

        q[58] = spatial_metric_der[2][0][0]; q[59] = spatial_metric_der[2][0][1]; q[60] = spatial_metric_der[2][0][2];
        q[61] = spatial_metric_der[2][1][0]; q[62] = spatial_metric_der[2][1][1]; q[63] = spatial_metric_der[2][1][2];
        q[64] = spatial_metric_der[2][2][0]; q[65] = spatial_metric_der[2][2][1]; q[66] = spatial_metric_der[2][2][2];

        q[67] = 0.0;
        q[68] = x; q[69] = y; q[70] = 0.0;

        double prims[71] = {0};
        gkyl_gr_euler_tetrad_prim_vars(gas_gamma, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        double fluxes[3][5] = {
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])) },
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

        double q_local[71] = {0}, flux_local_sr[71] = {0}, flux_local_gr[71] = {0}, flux[71] = {0};
        for (int d = 0; d < 3; d++) {
          gr_euler_tetrad->rotate_to_local_func(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_euler_tetrad_flux(gas_gamma, q_local, flux_local_sr);
          gkyl_gr_euler_tetrad_flux_correction(gas_gamma, q_local, flux_local_sr, flux_local_gr);
          gr_euler_tetrad->rotate_to_global_func(gr_euler_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 5; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[71] = {0}, q_g[71] = {0};
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[71] = {0}, q1[71] = {0};
          gr_euler_tetrad->cons_to_riem(gr_euler_tetrad, q_local, q_local, w1);
          gr_euler_tetrad->riem_to_cons(gr_euler_tetrad, q_local, w1, q1);

          for (int i = 0; i < 71; i++) {
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

  gkyl_wv_eqn_release(gr_euler_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_waves_minkowski_ho()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_euler_tetrad = gkyl_wv_gr_euler_tetrad_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;

      double spatial_det_l = 0.0, spatial_det_r = 0.0;
      double lapse_l = 0.0, lapse_r = 0.0;
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
      double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
      double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

      double ql[71], qr[71];
      ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
      ql[1] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * u_l;
      ql[2] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * v_l;
      ql[3] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * w_l;
      ql[4] = sqrt(spatial_det_l) * ((rho_l * h_l * (W_l * W_l)) - p_l - (rho_l * W_l));

      ql[5] = lapse_l;
      ql[6] = shift_l[0]; ql[7] = shift_l[1]; ql[8] = shift_l[2];

      ql[9] = spatial_metric_l[0][0]; ql[10] = spatial_metric_l[0][1]; ql[11] = spatial_metric_l[0][2];
      ql[12] = spatial_metric_l[1][0]; ql[13] = spatial_metric_l[1][1]; ql[14] = spatial_metric_l[1][2];
      ql[15] = spatial_metric_l[2][0]; ql[16] = spatial_metric_l[2][1]; ql[17] = spatial_metric_l[2][2];

      ql[18] = extrinsic_curvature_l[0][0]; ql[19] = extrinsic_curvature_l[0][1]; ql[20] = extrinsic_curvature_l[0][2];
      ql[21] = extrinsic_curvature_l[1][0]; ql[22] = extrinsic_curvature_l[1][1]; ql[23] = extrinsic_curvature_l[1][2];
      ql[24] = extrinsic_curvature_l[2][0]; ql[25] = extrinsic_curvature_l[2][1]; ql[26] = extrinsic_curvature_l[2][2];

      ql[27] = 1.0;

      ql[28] = lapse_der_l[0]; ql[29] = lapse_der_l[1]; ql[30] = lapse_der_l[2];
      ql[31] = shift_der_l[0][0]; ql[32] = shift_der_l[0][1]; ql[33] = shift_der_l[0][2];
      ql[34] = shift_der_l[1][0]; ql[35] = shift_der_l[1][1]; ql[36] = shift_der_l[1][2];
      ql[37] = shift_der_l[2][0]; ql[38] = shift_der_l[2][1]; ql[39] = shift_der_l[2][2];

      ql[40] = spatial_metric_der_l[0][0][0]; ql[41] = spatial_metric_der_l[0][0][1]; ql[42] = spatial_metric_der_l[0][0][2];
      ql[43] = spatial_metric_der_l[0][1][0]; ql[44] = spatial_metric_der_l[0][1][1]; ql[45] = spatial_metric_der_l[0][1][2];
      ql[46] = spatial_metric_der_l[0][2][0]; ql[47] = spatial_metric_der_l[0][2][1]; ql[48] = spatial_metric_der_l[0][2][2];

      ql[49] = spatial_metric_der_l[1][0][0]; ql[50] = spatial_metric_der_l[1][0][1]; ql[51] = spatial_metric_der_l[1][0][2];
      ql[52] = spatial_metric_der_l[1][1][0]; ql[53] = spatial_metric_der_l[1][1][1]; ql[54] = spatial_metric_der_l[1][1][2];
      ql[55] = spatial_metric_der_l[1][2][0]; ql[56] = spatial_metric_der_l[1][2][1]; ql[57] = spatial_metric_der_l[1][2][2];

      ql[58] = spatial_metric_der_l[2][0][0]; ql[59] = spatial_metric_der_l[2][0][1]; ql[60] = spatial_metric_der_l[2][0][2];
      ql[61] = spatial_metric_der_l[2][1][0]; ql[62] = spatial_metric_der_l[2][1][1]; ql[63] = spatial_metric_der_l[2][1][2];
      ql[64] = spatial_metric_der_l[2][2][0]; ql[65] = spatial_metric_der_l[2][2][1]; ql[66] = spatial_metric_der_l[2][2][2];

      ql[67] = 0.0;
      ql[68] = x - 0.5; ql[69] = y; ql[70] = 0.0;

      qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
      qr[1] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * u_r;
      qr[2] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * v_r;
      qr[3] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * w_r;
      qr[4] = sqrt(spatial_det_r) * ((rho_r * h_r * (W_r * W_r)) - p_r * (rho_r * W_r));

      qr[5] = lapse_r;
      qr[6] = shift_r[0]; qr[7] = shift_r[1]; qr[8] = shift_r[2];

      qr[9] = spatial_metric_r[0][0]; qr[10] = spatial_metric_r[0][1]; qr[11] = spatial_metric_r[0][2];
      qr[12] = spatial_metric_r[1][0]; qr[13] = spatial_metric_r[1][1]; qr[14] = spatial_metric_r[1][2];
      qr[15] = spatial_metric_r[2][0]; qr[16] = spatial_metric_r[2][1]; qr[17] = spatial_metric_r[2][2];

      qr[18] = extrinsic_curvature_r[0][0]; qr[19] = extrinsic_curvature_r[0][1]; qr[20] = extrinsic_curvature_r[0][2];
      qr[21] = extrinsic_curvature_r[1][0]; qr[22] = extrinsic_curvature_r[1][1]; qr[23] = extrinsic_curvature_r[1][2];
      qr[24] = extrinsic_curvature_r[2][0]; qr[25] = extrinsic_curvature_r[2][1]; qr[26] = extrinsic_curvature_r[2][2];

      qr[27] = 1.0;

      qr[28] = lapse_der_r[0]; qr[29] = lapse_der_r[1]; qr[30] = lapse_der_r[2];
      qr[31] = shift_der_r[0][0]; qr[32] = shift_der_r[0][1]; qr[33] = shift_der_r[0][2];
      qr[34] = shift_der_r[1][0]; qr[35] = shift_der_r[1][1]; qr[36] = shift_der_r[1][2];
      qr[37] = shift_der_r[2][0]; qr[38] = shift_der_r[2][1]; qr[39] = shift_der_r[2][2];

      qr[40] = spatial_metric_der_r[0][0][0]; qr[41] = spatial_metric_der_r[0][0][1]; qr[42] = spatial_metric_der_r[0][0][2];
      qr[43] = spatial_metric_der_r[0][1][0]; qr[44] = spatial_metric_der_r[0][1][1]; qr[45] = spatial_metric_der_r[0][1][2];
      qr[46] = spatial_metric_der_r[0][2][0]; qr[47] = spatial_metric_der_r[0][2][1]; qr[48] = spatial_metric_der_r[0][2][2];

      qr[49] = spatial_metric_der_r[1][0][0]; qr[50] = spatial_metric_der_r[1][0][1]; qr[51] = spatial_metric_der_r[1][0][2];
      qr[52] = spatial_metric_der_r[1][1][0]; qr[53] = spatial_metric_der_r[1][1][1]; qr[54] = spatial_metric_der_r[1][1][2];
      qr[55] = spatial_metric_der_r[1][2][0]; qr[56] = spatial_metric_der_r[1][2][1]; qr[57] = spatial_metric_der_r[1][2][2];

      qr[58] = spatial_metric_der_r[2][0][0]; qr[59] = spatial_metric_der_r[2][0][1]; qr[60] = spatial_metric_der_r[2][0][2];
      qr[61] = spatial_metric_der_r[2][1][0]; qr[62] = spatial_metric_der_r[2][1][1]; qr[63] = spatial_metric_der_r[2][1][2];
      qr[64] = spatial_metric_der_r[2][2][0]; qr[65] = spatial_metric_der_r[2][2][1]; qr[66] = spatial_metric_der_r[2][2][2];

      qr[67] = 0.0;
      qr[68] = x + 0.5; qr[69] = y; qr[70] = 0.0;

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
        double speeds[3] = {0}, waves[3 * 71] = {0}, waves_local[3 * 71] = {0};

        double ql_local[71] = {0}, qr_local[71] = {0};
        gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
        gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

        double delta[71] = {0};
        for (int i = 0; i < 71; i++) {
          delta[i] = qr_local[i] - ql_local[i];
        }

        gkyl_wv_eqn_waves(gr_euler_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

        double apdq_local[71] = {0}, amdq_local[71] = {0};
        gkyl_wv_eqn_qfluct(gr_euler_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

        for (int i = 0; i < 3; i++) {
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 71], &waves[i * 71]);
        }

        double apdq[71] = {0}, amdq[71] = {0};
        gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
        gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

        double fl_local_sr[71] = {0}, fr_local_sr[71] = {0};
        gkyl_gr_euler_tetrad_flux(gas_gamma, ql_local, fl_local_sr);
        gkyl_gr_euler_tetrad_flux(gas_gamma, qr_local, fr_local_sr);

        double fl_local_gr[71], fr_local_gr[71];
        gkyl_gr_euler_tetrad_flux_correction(gas_gamma, ql_local, fl_local_sr, fl_local_gr);
        gkyl_gr_euler_tetrad_flux_correction(gas_gamma, qr_local, fr_local_sr, fr_local_gr);

        double fl[71], fr[71];
        gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
        gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

        for (int i = 0; i < 71; i++) {
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

  gkyl_wv_eqn_release(gr_euler_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_waves_schwarzschild_ho()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler_tetrad = gkyl_wv_gr_euler_tetrad_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;

      double spatial_det_l = 0.0, spatial_det_r = 0.0;
      double lapse_l = 0.0, lapse_r = 0.0;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l = false, in_excision_region_r = false;

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
      double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
      double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[71], qr[71];
        ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
        ql[1] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * u_l;
        ql[2] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * v_l;
        ql[3] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * w_l;
        ql[4] = sqrt(spatial_det_l) * ((rho_l * h_l * (W_l * W_l)) - p_l - (rho_l * W_l));

        ql[5] = lapse_l;
        ql[6] = shift_l[0]; ql[7] = shift_l[1]; ql[8] = shift_l[2];

        ql[9] = spatial_metric_l[0][0]; ql[10] = spatial_metric_l[0][1]; ql[11] = spatial_metric_l[0][2];
        ql[12] = spatial_metric_l[1][0]; ql[13] = spatial_metric_l[1][1]; ql[14] = spatial_metric_l[1][2];
        ql[15] = spatial_metric_l[2][0]; ql[16] = spatial_metric_l[2][1]; ql[17] = spatial_metric_l[2][2];

        ql[18] = extrinsic_curvature_l[0][0]; ql[19] = extrinsic_curvature_l[0][1]; ql[20] = extrinsic_curvature_l[0][2];
        ql[21] = extrinsic_curvature_l[1][0]; ql[22] = extrinsic_curvature_l[1][1]; ql[23] = extrinsic_curvature_l[1][2];
        ql[24] = extrinsic_curvature_l[2][0]; ql[25] = extrinsic_curvature_l[2][1]; ql[26] = extrinsic_curvature_l[2][2];

        ql[27] = 1.0;

        ql[28] = lapse_der_l[0]; ql[29] = lapse_der_l[1]; ql[30] = lapse_der_l[2];
        ql[31] = shift_der_l[0][0]; ql[32] = shift_der_l[0][1]; ql[33] = shift_der_l[0][2];
        ql[34] = shift_der_l[1][0]; ql[35] = shift_der_l[1][1]; ql[36] = shift_der_l[1][2];
        ql[37] = shift_der_l[2][0]; ql[38] = shift_der_l[2][1]; ql[39] = shift_der_l[2][2];

        ql[40] = spatial_metric_der_l[0][0][0]; ql[41] = spatial_metric_der_l[0][0][1]; ql[42] = spatial_metric_der_l[0][0][2];
        ql[43] = spatial_metric_der_l[0][1][0]; ql[44] = spatial_metric_der_l[0][1][1]; ql[45] = spatial_metric_der_l[0][1][2];
        ql[46] = spatial_metric_der_l[0][2][0]; ql[47] = spatial_metric_der_l[0][2][1]; ql[48] = spatial_metric_der_l[0][2][2];

        ql[49] = spatial_metric_der_l[1][0][0]; ql[50] = spatial_metric_der_l[1][0][1]; ql[51] = spatial_metric_der_l[1][0][2];
        ql[52] = spatial_metric_der_l[1][1][0]; ql[53] = spatial_metric_der_l[1][1][1]; ql[54] = spatial_metric_der_l[1][1][2];
        ql[55] = spatial_metric_der_l[1][2][0]; ql[56] = spatial_metric_der_l[1][2][1]; ql[57] = spatial_metric_der_l[1][2][2];

        ql[58] = spatial_metric_der_l[2][0][0]; ql[59] = spatial_metric_der_l[2][0][1]; ql[60] = spatial_metric_der_l[2][0][2];
        ql[61] = spatial_metric_der_l[2][1][0]; ql[62] = spatial_metric_der_l[2][1][1]; ql[63] = spatial_metric_der_l[2][1][2];
        ql[64] = spatial_metric_der_l[2][2][0]; ql[65] = spatial_metric_der_l[2][2][1]; ql[66] = spatial_metric_der_l[2][2][2];

        ql[67] = 0.0;
        ql[68] = x - 0.5; ql[69] = y; ql[70] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
        qr[1] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * u_r;
        qr[2] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * v_r;
        qr[3] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * w_r;
        qr[4] = sqrt(spatial_det_r) * ((rho_r * h_r * (W_r * W_r)) - p_r * (rho_r * W_r));

        qr[5] = lapse_r;
        qr[6] = shift_r[0]; qr[7] = shift_r[1]; qr[8] = shift_r[2];

        qr[9] = spatial_metric_r[0][0]; qr[10] = spatial_metric_r[0][1]; qr[11] = spatial_metric_r[0][2];
        qr[12] = spatial_metric_r[1][0]; qr[13] = spatial_metric_r[1][1]; qr[14] = spatial_metric_r[1][2];
        qr[15] = spatial_metric_r[2][0]; qr[16] = spatial_metric_r[2][1]; qr[17] = spatial_metric_r[2][2];

        qr[18] = extrinsic_curvature_r[0][0]; qr[19] = extrinsic_curvature_r[0][1]; qr[20] = extrinsic_curvature_r[0][2];
        qr[21] = extrinsic_curvature_r[1][0]; qr[22] = extrinsic_curvature_r[1][1]; qr[23] = extrinsic_curvature_r[1][2];
        qr[24] = extrinsic_curvature_r[2][0]; qr[25] = extrinsic_curvature_r[2][1]; qr[26] = extrinsic_curvature_r[2][2];

        qr[27] = 1.0;

        qr[28] = lapse_der_r[0]; qr[29] = lapse_der_r[1]; qr[30] = lapse_der_r[2];
        qr[31] = shift_der_r[0][0]; qr[32] = shift_der_r[0][1]; qr[33] = shift_der_r[0][2];
        qr[34] = shift_der_r[1][0]; qr[35] = shift_der_r[1][1]; qr[36] = shift_der_r[1][2];
        qr[37] = shift_der_r[2][0]; qr[38] = shift_der_r[2][1]; qr[39] = shift_der_r[2][2];

        qr[40] = spatial_metric_der_r[0][0][0]; qr[41] = spatial_metric_der_r[0][0][1]; qr[42] = spatial_metric_der_r[0][0][2];
        qr[43] = spatial_metric_der_r[0][1][0]; qr[44] = spatial_metric_der_r[0][1][1]; qr[45] = spatial_metric_der_r[0][1][2];
        qr[46] = spatial_metric_der_r[0][2][0]; qr[47] = spatial_metric_der_r[0][2][1]; qr[48] = spatial_metric_der_r[0][2][2];

        qr[49] = spatial_metric_der_r[1][0][0]; qr[50] = spatial_metric_der_r[1][0][1]; qr[51] = spatial_metric_der_r[1][0][2];
        qr[52] = spatial_metric_der_r[1][1][0]; qr[53] = spatial_metric_der_r[1][1][1]; qr[54] = spatial_metric_der_r[1][1][2];
        qr[55] = spatial_metric_der_r[1][2][0]; qr[56] = spatial_metric_der_r[1][2][1]; qr[57] = spatial_metric_der_r[1][2][2];

        qr[58] = spatial_metric_der_r[2][0][0]; qr[59] = spatial_metric_der_r[2][0][1]; qr[60] = spatial_metric_der_r[2][0][2];
        qr[61] = spatial_metric_der_r[2][1][0]; qr[62] = spatial_metric_der_r[2][1][1]; qr[63] = spatial_metric_der_r[2][1][2];
        qr[64] = spatial_metric_der_r[2][2][0]; qr[65] = spatial_metric_der_r[2][2][1]; qr[66] = spatial_metric_der_r[2][2][2];

        qr[67] = 0.0;
        qr[68] = x + 0.5; qr[69] = y; qr[70] = 0.0;

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
          double speeds[3] = {0}, waves[3 * 71] = {0}, waves_local[3 * 71] = {0};

          double ql_local[71] = {0}, qr_local[71] = {0};
          gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[71] = {0};
          for (int i = 0; i < 71; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_euler_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[71] = {0}, amdq_local[71] = {0};
          gkyl_wv_eqn_qfluct(gr_euler_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 3; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 71], &waves[i * 71]);
          }

          double apdq[71] = {0}, amdq[71] = {0};
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local_sr[71] = {0}, fr_local_sr[71] = {0};
          gkyl_gr_euler_tetrad_flux(gas_gamma, ql_local, fl_local_sr);
          gkyl_gr_euler_tetrad_flux(gas_gamma, qr_local, fr_local_sr);

          double fl_local_gr[71], fr_local_gr[71];
          gkyl_gr_euler_tetrad_flux_correction(gas_gamma, ql_local, fl_local_sr, fl_local_gr);
          gkyl_gr_euler_tetrad_flux_correction(gas_gamma, qr_local, fr_local_sr, fr_local_gr);

          double fl[71], fr[71];
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

          for (int i = 0; i < 71; i++) {
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

  gkyl_wv_eqn_release(gr_euler_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_waves_kerr_ho()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler_tetrad = gkyl_wv_gr_euler_tetrad_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;

      double spatial_det_l = 0.0, spatial_det_r = 0.0;
      double lapse_l = 0.0, lapse_r = 0.0;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l = false, in_excision_region_r = false;

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
      double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
      double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[71], qr[71];
        ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
        ql[1] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * u_l;
        ql[2] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * v_l;
        ql[3] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * w_l;
        ql[4] = sqrt(spatial_det_l) * ((rho_l * h_l * (W_l * W_l)) - p_l - (rho_l * W_l));

        ql[5] = lapse_l;
        ql[6] = shift_l[0]; ql[7] = shift_l[1]; ql[8] = shift_l[2];

        ql[9] = spatial_metric_l[0][0]; ql[10] = spatial_metric_l[0][1]; ql[11] = spatial_metric_l[0][2];
        ql[12] = spatial_metric_l[1][0]; ql[13] = spatial_metric_l[1][1]; ql[14] = spatial_metric_l[1][2];
        ql[15] = spatial_metric_l[2][0]; ql[16] = spatial_metric_l[2][1]; ql[17] = spatial_metric_l[2][2];

        ql[18] = extrinsic_curvature_l[0][0]; ql[19] = extrinsic_curvature_l[0][1]; ql[20] = extrinsic_curvature_l[0][2];
        ql[21] = extrinsic_curvature_l[1][0]; ql[22] = extrinsic_curvature_l[1][1]; ql[23] = extrinsic_curvature_l[1][2];
        ql[24] = extrinsic_curvature_l[2][0]; ql[25] = extrinsic_curvature_l[2][1]; ql[26] = extrinsic_curvature_l[2][2];

        ql[27] = 1.0;

        ql[28] = lapse_der_l[0]; ql[29] = lapse_der_l[1]; ql[30] = lapse_der_l[2];
        ql[31] = shift_der_l[0][0]; ql[32] = shift_der_l[0][1]; ql[33] = shift_der_l[0][2];
        ql[34] = shift_der_l[1][0]; ql[35] = shift_der_l[1][1]; ql[36] = shift_der_l[1][2];
        ql[37] = shift_der_l[2][0]; ql[38] = shift_der_l[2][1]; ql[39] = shift_der_l[2][2];

        ql[40] = spatial_metric_der_l[0][0][0]; ql[41] = spatial_metric_der_l[0][0][1]; ql[42] = spatial_metric_der_l[0][0][2];
        ql[43] = spatial_metric_der_l[0][1][0]; ql[44] = spatial_metric_der_l[0][1][1]; ql[45] = spatial_metric_der_l[0][1][2];
        ql[46] = spatial_metric_der_l[0][2][0]; ql[47] = spatial_metric_der_l[0][2][1]; ql[48] = spatial_metric_der_l[0][2][2];

        ql[49] = spatial_metric_der_l[1][0][0]; ql[50] = spatial_metric_der_l[1][0][1]; ql[51] = spatial_metric_der_l[1][0][2];
        ql[52] = spatial_metric_der_l[1][1][0]; ql[53] = spatial_metric_der_l[1][1][1]; ql[54] = spatial_metric_der_l[1][1][2];
        ql[55] = spatial_metric_der_l[1][2][0]; ql[56] = spatial_metric_der_l[1][2][1]; ql[57] = spatial_metric_der_l[1][2][2];

        ql[58] = spatial_metric_der_l[2][0][0]; ql[59] = spatial_metric_der_l[2][0][1]; ql[60] = spatial_metric_der_l[2][0][2];
        ql[61] = spatial_metric_der_l[2][1][0]; ql[62] = spatial_metric_der_l[2][1][1]; ql[63] = spatial_metric_der_l[2][1][2];
        ql[64] = spatial_metric_der_l[2][2][0]; ql[65] = spatial_metric_der_l[2][2][1]; ql[66] = spatial_metric_der_l[2][2][2];

        ql[67] = 0.0;
        ql[68] = x - 0.5; ql[69] = y; ql[70] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
        qr[1] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * u_r;
        qr[2] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * v_r;
        qr[3] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * w_r;
        qr[4] = sqrt(spatial_det_r) * ((rho_r * h_r * (W_r * W_r)) - p_r * (rho_r * W_r));

        qr[5] = lapse_r;
        qr[6] = shift_r[0]; qr[7] = shift_r[1]; qr[8] = shift_r[2];

        qr[9] = spatial_metric_r[0][0]; qr[10] = spatial_metric_r[0][1]; qr[11] = spatial_metric_r[0][2];
        qr[12] = spatial_metric_r[1][0]; qr[13] = spatial_metric_r[1][1]; qr[14] = spatial_metric_r[1][2];
        qr[15] = spatial_metric_r[2][0]; qr[16] = spatial_metric_r[2][1]; qr[17] = spatial_metric_r[2][2];

        qr[18] = extrinsic_curvature_r[0][0]; qr[19] = extrinsic_curvature_r[0][1]; qr[20] = extrinsic_curvature_r[0][2];
        qr[21] = extrinsic_curvature_r[1][0]; qr[22] = extrinsic_curvature_r[1][1]; qr[23] = extrinsic_curvature_r[1][2];
        qr[24] = extrinsic_curvature_r[2][0]; qr[25] = extrinsic_curvature_r[2][1]; qr[26] = extrinsic_curvature_r[2][2];

        qr[27] = 1.0;

        qr[28] = lapse_der_r[0]; qr[29] = lapse_der_r[1]; qr[30] = lapse_der_r[2];
        qr[31] = shift_der_r[0][0]; qr[32] = shift_der_r[0][1]; qr[33] = shift_der_r[0][2];
        qr[34] = shift_der_r[1][0]; qr[35] = shift_der_r[1][1]; qr[36] = shift_der_r[1][2];
        qr[37] = shift_der_r[2][0]; qr[38] = shift_der_r[2][1]; qr[39] = shift_der_r[2][2];

        qr[40] = spatial_metric_der_r[0][0][0]; qr[41] = spatial_metric_der_r[0][0][1]; qr[42] = spatial_metric_der_r[0][0][2];
        qr[43] = spatial_metric_der_r[0][1][0]; qr[44] = spatial_metric_der_r[0][1][1]; qr[45] = spatial_metric_der_r[0][1][2];
        qr[46] = spatial_metric_der_r[0][2][0]; qr[47] = spatial_metric_der_r[0][2][1]; qr[48] = spatial_metric_der_r[0][2][2];

        qr[49] = spatial_metric_der_r[1][0][0]; qr[50] = spatial_metric_der_r[1][0][1]; qr[51] = spatial_metric_der_r[1][0][2];
        qr[52] = spatial_metric_der_r[1][1][0]; qr[53] = spatial_metric_der_r[1][1][1]; qr[54] = spatial_metric_der_r[1][1][2];
        qr[55] = spatial_metric_der_r[1][2][0]; qr[56] = spatial_metric_der_r[1][2][1]; qr[57] = spatial_metric_der_r[1][2][2];

        qr[58] = spatial_metric_der_r[2][0][0]; qr[59] = spatial_metric_der_r[2][0][1]; qr[60] = spatial_metric_der_r[2][0][2];
        qr[61] = spatial_metric_der_r[2][1][0]; qr[62] = spatial_metric_der_r[2][1][1]; qr[63] = spatial_metric_der_r[2][1][2];
        qr[64] = spatial_metric_der_r[2][2][0]; qr[65] = spatial_metric_der_r[2][2][1]; qr[66] = spatial_metric_der_r[2][2][2];

        qr[67] = 0.0;
        qr[68] = x + 0.5; qr[69] = y; qr[70] = 0.0;

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
          double speeds[3] = {0}, waves[3 * 71] = {0}, waves_local[3 * 71] = {0};

          double ql_local[71] = {0}, qr_local[71] = {0};
          gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_euler_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[71] = {0};
          for (int i = 0; i < 71; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_euler_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[71] = {0}, amdq_local[71] = {0};
          gkyl_wv_eqn_qfluct(gr_euler_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 3; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 71], &waves[i * 71]);
          }

          double apdq[71] = {0}, amdq[71] = {0};
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local_sr[71] = {0}, fr_local_sr[71] = {0};
          gkyl_gr_euler_tetrad_flux(gas_gamma, ql_local, fl_local_sr);
          gkyl_gr_euler_tetrad_flux(gas_gamma, qr_local, fr_local_sr);

          double fl_local_gr[71], fr_local_gr[71];
          gkyl_gr_euler_tetrad_flux_correction(gas_gamma, ql_local, fl_local_sr, fl_local_gr);
          gkyl_gr_euler_tetrad_flux_correction(gas_gamma, qr_local, fr_local_sr, fr_local_gr);

          double fl[71], fr[71];
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
          gkyl_wv_eqn_rotate_to_global(gr_euler_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

          for (int i = 0; i < 71; i++) {
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

  gkyl_wv_eqn_release(gr_euler_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_euler_tetrad_basic_minkowski_ho", test_gr_euler_tetrad_basic_minkowski_ho },
  { "gr_euler_tetrad_basic_schwarzschild_ho", test_gr_euler_tetrad_basic_schwarzschild_ho },
  { "gr_euler_tetrad_basic_kerr_ho", test_gr_euler_tetrad_basic_kerr_ho },
  { "gr_euler_tetrad_waves_minkowski_ho", test_gr_euler_tetrad_waves_minkowski_ho },
  { "gr_euler_tetrad_waves_schwarzschild_ho", test_gr_euler_tetrad_waves_schwarzschild_ho },
  { "gr_euler_tetrad_waves_kerr_ho", test_gr_euler_tetrad_waves_kerr_ho },
  { NULL, NULL },
};
