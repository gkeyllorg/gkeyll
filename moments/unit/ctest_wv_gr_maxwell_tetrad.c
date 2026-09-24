#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_maxwell_tetrad.h>
#include <gkyl_wv_gr_maxwell_tetrad_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_maxwell_tetrad_basic_minkowski()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_maxwell_tetrad->num_equations == 26 );
  TEST_CHECK( gr_maxwell_tetrad->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);

      double q[26];
      q[0] = Dx; q[1] = Dy; q[2] = Dz;
      q[3] = Bx; q[4] = By; q[5] = Bz;
      q[6] = phi; q[7] = psi;

      q[8] = lapse;
      q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

      q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
      q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
      q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

      q[21] = 1.0;

      q[22] = 0.0;
      q[23] = x; q[24] = y; q[25] = 0.0;

      double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
      double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
      double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

      double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
      double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
      double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

      double fluxes[3][8] = {
        { e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
          -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
        { -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
          -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
        { (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
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

      double q_local[26], flux_local_sr[26], flux_local_gr[26], flux[26];
      for (int d = 0; d < 3; d++) {
        gr_maxwell_tetrad->rotate_to_local_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, q_local, flux_local_sr);
        gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
        gr_maxwell_tetrad->rotate_to_global_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-13) );
        }
      }

      double q_l[26], q_g[26];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
        }

        double w1[26], q1[26];
        gr_maxwell_tetrad->cons_to_riem(gr_maxwell_tetrad, q_local, q_local, w1);
        gr_maxwell_tetrad->riem_to_cons(gr_maxwell_tetrad, q_local, w1, q1);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(shift);
    }
  }

  gkyl_wv_eqn_release(gr_maxwell_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_basic_schwarzschild()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_maxwell_tetrad->num_equations == 26 );
  TEST_CHECK( gr_maxwell_tetrad->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);

      if (!in_excision_region) {
        double q[26];
        q[0] = Dx; q[1] = Dy; q[2] = Dz;
        q[3] = Bx; q[4] = By; q[5] = Bz;
        q[6] = phi; q[7] = psi;

        q[8] = lapse;
        q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

        q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
        q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
        q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

        q[21] = 1.0;

        q[22] = 0.0;
        q[23] = x; q[24] = y; q[25] = 0.0;

        double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
        double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
        double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

        double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
        double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
        double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

        double fluxes[3][8] = {
          { e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
            -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
          { -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
            -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
          { (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
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

        double q_local[26], flux_local_sr[26], flux_local_gr[26], flux[26];
        for (int d = 0; d < 3; d++) {
          gr_maxwell_tetrad->rotate_to_local_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, q_local, flux_local_sr);
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
          gr_maxwell_tetrad->rotate_to_global_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-13) );
          }
        }

        double q_l[26], q_g[26];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[26], q1[26];
          gr_maxwell_tetrad->cons_to_riem(gr_maxwell_tetrad, q_local, q_local, w1);
          gr_maxwell_tetrad->riem_to_cons(gr_maxwell_tetrad, q_local, w1, q1);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(shift);
    }
  }

  gkyl_wv_eqn_release(gr_maxwell_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_basic_kerr()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_maxwell_tetrad->num_equations == 26 );
  TEST_CHECK( gr_maxwell_tetrad->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);

      if (!in_excision_region) {
        double q[26];
        q[0] = Dx; q[1] = Dy; q[2] = Dz;
        q[3] = Bx; q[4] = By; q[5] = Bz;
        q[6] = phi; q[7] = psi;

        q[8] = lapse;
        q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

        q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
        q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
        q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

        q[21] = 1.0;

        q[22] = 0.0;
        q[23] = x; q[24] = y; q[25] = 0.0;

        double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
        double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
        double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

        double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
        double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
        double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

        double fluxes[3][8] = {
          { e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
            -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
          { -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
            -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
          { (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
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

        double q_local[26], flux_local_sr[26], flux_local_gr[26], flux[26];
        for (int d = 0; d < 3; d++) {
          gr_maxwell_tetrad->rotate_to_local_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, q_local, flux_local_sr);
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
          gr_maxwell_tetrad->rotate_to_global_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-13) );
          }
        }

        double q_l[26], q_g[26];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[26], q1[26];
          gr_maxwell_tetrad->cons_to_riem(gr_maxwell_tetrad, q_local, q_local, w1);
          gr_maxwell_tetrad->riem_to_cons(gr_maxwell_tetrad, q_local, w1, q1);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(shift);
    }
  }

  gkyl_wv_eqn_release(gr_maxwell_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_waves_minkowski()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx_l = 0.1, Dy_l = 0.2, Dz_l = 0.3;
      double Bx_l = 0.4, By_l = 0.5, Bz_l = 0.6;
      double phi_l = 0.0, psi_l = 0.0;

      double Dx_r = 1.0, Dy_r = 0.9, Dz_r = 0.8;
      double Bx_r = 0.7, By_r = 0.6, Bz_r = 0.5;
      double phi_r = 0.0, psi_r = 0.0;

      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      
      double ql[26], qr[26];
      ql[0] = Dx_l; ql[1] = Dy_l; ql[2] = Dz_l;
      ql[3] = Bx_l; ql[4] = By_l; ql[5] = Bz_l;
      ql[6] = phi_l; ql[7] = psi_l;

      ql[8] = lapse_l;
      ql[9] = shift_l[0]; ql[10] = shift_l[1]; ql[11] = shift_l[2];

      ql[12] = spatial_metric_l[0][0]; ql[13] = spatial_metric_l[0][1]; ql[14] = spatial_metric_l[0][2];
      ql[15] = spatial_metric_l[1][0]; ql[16] = spatial_metric_l[1][1]; ql[17] = spatial_metric_l[1][2];
      ql[18] = spatial_metric_l[2][0]; ql[19] = spatial_metric_l[2][1]; ql[20] = spatial_metric_l[2][2];

      ql[21] = 1.0;

      ql[22] = 0.0;
      ql[23] = x - 0.5; ql[24] = y; ql[25] = 0.0;

      qr[0] = Dx_r; qr[1] = Dy_r; qr[2] = Dz_r;
      qr[3] = Bx_r; qr[4] = By_r; qr[5] = Bz_r;
      qr[6] = phi_r; qr[7] = psi_r;

      qr[8] = lapse_r;
      qr[9] = shift_r[0]; qr[10] = shift_r[1]; qr[11] = shift_r[2];

      qr[12] = spatial_metric_r[0][0]; qr[13] = spatial_metric_r[0][1]; qr[14] = spatial_metric_r[0][2];
      qr[15] = spatial_metric_r[1][0]; qr[16] = spatial_metric_r[1][1]; qr[17] = spatial_metric_r[1][2];
      qr[18] = spatial_metric_r[2][0]; qr[19] = spatial_metric_r[2][1]; qr[20] = spatial_metric_r[2][2];

      qr[21] = 1.0;

      qr[22] = 0.0;
      qr[23] = x + 0.5; qr[24] = y; qr[25] = 0.0;

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
        double speeds[6], waves[6 * 26], waves_local[6 * 26];

        double ql_local[26], qr_local[26];
        gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
        gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

        double delta[26];
        for (int i = 0; i < 26; i++) {
          delta[i] = qr_local[i] - ql_local[i];
        }

        gkyl_wv_eqn_waves(gr_maxwell_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

        double apdq_local[26], amdq_local[26];
        gkyl_wv_eqn_qfluct(gr_maxwell_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

        for (int i = 0; i < 2; i++) {
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 26], &waves[i * 26]);
        }

        double apdq[26], amdq[26];
        gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
        gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

        double fl_local_sr[26], fr_local_sr[26];
        gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, ql_local, fl_local_sr);
        gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, qr_local, fr_local_sr);

        double fl_local_gr[26], fr_local_gr[26];
        gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, ql_local, fl_local_sr, fl_local_gr);
        gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, qr_local, fr_local_sr, fr_local_gr);

        double fl[26], fr[26];
        gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
        gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

        for (int i = 0; i < 26; i++) {
          TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
        }
      }
      
      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
    }
  }
  
  gkyl_wv_eqn_release(gr_maxwell_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_waves_schwarzschild()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx_l = 0.1, Dy_l = 0.2, Dz_l = 0.3;
      double Bx_l = 0.4, By_l = 0.5, Bz_l = 0.6;
      double phi_l = 0.0, psi_l = 0.0;

      double Dx_r = 1.0, Dy_r = 0.9, Dz_r = 0.8;
      double Bx_r = 0.7, By_r = 0.6, Bz_r = 0.5;
      double phi_r = 0.0, psi_r = 0.0;

      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[26], qr[26];
        ql[0] = Dx_l; ql[1] = Dy_l; ql[2] = Dz_l;
        ql[3] = Bx_l; ql[4] = By_l; ql[5] = Bz_l;
        ql[6] = phi_l; ql[7] = psi_l;

        ql[8] = lapse_l;
        ql[9] = shift_l[0]; ql[10] = shift_l[1]; ql[11] = shift_l[2];

        ql[12] = spatial_metric_l[0][0]; ql[13] = spatial_metric_l[0][1]; ql[14] = spatial_metric_l[0][2];
        ql[15] = spatial_metric_l[1][0]; ql[16] = spatial_metric_l[1][1]; ql[17] = spatial_metric_l[1][2];
        ql[18] = spatial_metric_l[2][0]; ql[19] = spatial_metric_l[2][1]; ql[20] = spatial_metric_l[2][2];

        ql[21] = 1.0;

        ql[22] = 0.0;
        ql[23] = x - 0.5; ql[24] = y; ql[25] = 0.0;

        qr[0] = Dx_r; qr[1] = Dy_r; qr[2] = Dz_r;
        qr[3] = Bx_r; qr[4] = By_r; qr[5] = Bz_r;
        qr[6] = phi_r; qr[7] = psi_r;

        qr[8] = lapse_r;
        qr[9] = shift_r[0]; qr[10] = shift_r[1]; qr[11] = shift_r[2];

        qr[12] = spatial_metric_r[0][0]; qr[13] = spatial_metric_r[0][1]; qr[14] = spatial_metric_r[0][2];
        qr[15] = spatial_metric_r[1][0]; qr[16] = spatial_metric_r[1][1]; qr[17] = spatial_metric_r[1][2];
        qr[18] = spatial_metric_r[2][0]; qr[19] = spatial_metric_r[2][1]; qr[20] = spatial_metric_r[2][2];

        qr[21] = 1.0;

        qr[22] = 0.0;
        qr[23] = x + 0.5; qr[24] = y; qr[25] = 0.0;

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
          double speeds[6], waves[6 * 26], waves_local[6 * 26];

          double ql_local[26], qr_local[26];
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[26];
          for (int i = 0; i < 26; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_maxwell_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[26], amdq_local[26];
          gkyl_wv_eqn_qfluct(gr_maxwell_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 2; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 26], &waves[i * 26]);
          }

          double apdq[26], amdq[26];
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local_sr[26], fr_local_sr[26];
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, ql_local, fl_local_sr);
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, qr_local, fr_local_sr);

          double fl_local_gr[26], fr_local_gr[26];
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, ql_local, fl_local_sr, fl_local_gr);
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, qr_local, fr_local_sr, fr_local_gr);

          double fl[26], fr[26];
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

          for (int i = 0; i < 26; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }
      
      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
    }
  }
  
  gkyl_wv_eqn_release(gr_maxwell_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_waves_kerr()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx_l = 0.1, Dy_l = 0.2, Dz_l = 0.3;
      double Bx_l = 0.4, By_l = 0.5, Bz_l = 0.6;
      double phi_l = 0.0, psi_l = 0.0;

      double Dx_r = 1.0, Dy_r = 0.9, Dz_r = 0.8;
      double Bx_r = 0.7, By_r = 0.6, Bz_r = 0.5;
      double phi_r = 0.0, psi_r = 0.0;

      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[26], qr[26];
        ql[0] = Dx_l; ql[1] = Dy_l; ql[2] = Dz_l;
        ql[3] = Bx_l; ql[4] = By_l; ql[5] = Bz_l;
        ql[6] = phi_l; ql[7] = psi_l;

        ql[8] = lapse_l;
        ql[9] = shift_l[0]; ql[10] = shift_l[1]; ql[11] = shift_l[2];

        ql[12] = spatial_metric_l[0][0]; ql[13] = spatial_metric_l[0][1]; ql[14] = spatial_metric_l[0][2];
        ql[15] = spatial_metric_l[1][0]; ql[16] = spatial_metric_l[1][1]; ql[17] = spatial_metric_l[1][2];
        ql[18] = spatial_metric_l[2][0]; ql[19] = spatial_metric_l[2][1]; ql[20] = spatial_metric_l[2][2];

        ql[21] = 1.0;

        ql[22] = 0.0;
        ql[23] = x - 0.5; ql[24] = y; ql[25] = 0.0;

        qr[0] = Dx_r; qr[1] = Dy_r; qr[2] = Dz_r;
        qr[3] = Bx_r; qr[4] = By_r; qr[5] = Bz_r;
        qr[6] = phi_r; qr[7] = psi_r;

        qr[8] = lapse_r;
        qr[9] = shift_r[0]; qr[10] = shift_r[1]; qr[11] = shift_r[2];

        qr[12] = spatial_metric_r[0][0]; qr[13] = spatial_metric_r[0][1]; qr[14] = spatial_metric_r[0][2];
        qr[15] = spatial_metric_r[1][0]; qr[16] = spatial_metric_r[1][1]; qr[17] = spatial_metric_r[1][2];
        qr[18] = spatial_metric_r[2][0]; qr[19] = spatial_metric_r[2][1]; qr[20] = spatial_metric_r[2][2];

        qr[21] = 1.0;

        qr[22] = 0.0;
        qr[23] = x + 0.5; qr[24] = y; qr[25] = 0.0;

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
          double speeds[6], waves[6 * 26], waves_local[6 * 26];

          double ql_local[26], qr_local[26];
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[26];
          for (int i = 0; i < 26; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_maxwell_tetrad, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[26], amdq_local[26];
          gkyl_wv_eqn_qfluct(gr_maxwell_tetrad, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 2; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], &waves_local[i * 26], &waves[i * 26]);
          }

          double apdq[26], amdq[26];
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local_sr[26], fr_local_sr[26];
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, ql_local, fl_local_sr);
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, qr_local, fr_local_sr);

          double fl_local_gr[26], fr_local_gr[26];
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, ql_local, fl_local_sr, fl_local_gr);
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, qr_local, fr_local_sr, fr_local_gr);

          double fl[26], fr[26];
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], fl_local_gr, fl);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], fr_local_gr, fr);

          for (int i = 0; i < 26; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
          }
        }
      }
      
      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
    }
  }
  
  gkyl_wv_eqn_release(gr_maxwell_tetrad);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_maxwell_tetrad_basic_minkowski", test_gr_maxwell_tetrad_basic_minkowski },
  { "gr_maxwell_tetrad_basic_schwarzschild", test_gr_maxwell_tetrad_basic_schwarzschild },
  { "gr_maxwell_tetrad_basic_kerr", test_gr_maxwell_tetrad_basic_kerr },
  { "gr_maxwell_tetrad_waves_minkowski", test_gr_maxwell_tetrad_waves_minkowski },
  { "gr_maxwell_tetrad_waves_schwarzschild", test_gr_maxwell_tetrad_waves_schwarzschild },
  { "gr_maxwell_tetrad_waves_kerr", test_gr_maxwell_tetrad_waves_kerr },
  { NULL, NULL },
};
