#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_medium.h>
#include <gkyl_wv_gr_medium_priv.h>

void
test_gr_medium_basic_ho()
{
  double gas_gamma = 5.0 / 3.0;
  double kappa = 8.0 * M_PI;
  struct gkyl_wv_eqn *gr_medium = gkyl_wv_gr_medium_new(gas_gamma, kappa, false);

  TEST_CHECK( gr_medium->num_equations == 15 );
  TEST_CHECK( gr_medium->num_waves == 2 );

  double exp_2a = 0.1;
  double a_dt = 0.2, a_dx = 0.3;
  double b_dt = 0.4, b_dx = 0.5;
  double c_dt = 0.6, c_dx = 0.7;

  double a_dt_dx = 0.1, a_dx_dx = 0.2;
  double b_dt_dx = 0.3, b_dx_dx = 0.4;
  double c_dt_dx = 0.5, c_dx_dx = 0.6;

  double rho = 1.0, vel = 0.2;
  double p = (gas_gamma - 1.0) * rho;
  double W = 1.0 / sqrt(1.0 - (vel * vel));

  double q[15] = { exp_2a, a_dt, a_dx, b_dt, b_dx, c_dt, c_dx, a_dt_dx, a_dx_dx, b_dt_dx, b_dx_dx, c_dt_dx, c_dx_dx,
    ((rho + p) * (W * W)) - p, (rho + p) * vel * (W * W) };

  double prims[15];
  gkyl_gr_medium_prim_vars(gas_gamma, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], exp_2a, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[1], a_dt, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[2], a_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[3], b_dt, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[4], b_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[5], c_dt, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[6], c_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[7], a_dt_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[8], a_dx_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[9], b_dt_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[10], b_dx_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[11], c_dt_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[12], c_dx_dx, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[13], rho, 1e-14) );
  TEST_CHECK( gkyl_compare(prims[14], vel, 1e-14) );

  double Etot = ((rho + p) * (W * W)) - p;
  double mom = (rho + p) * vel * (W * W);

  double fluxes[3][15] = {
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     -a_dx_dx + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))), -a_dt_dx,
     -b_dx_dx - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))), -b_dt_dx,
     -c_dx_dx, -c_dt_dx, mom, (mom * vel) + p },
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     -a_dx_dx + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))), -a_dt_dx,
     -b_dx_dx - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))), -b_dt_dx,
     -c_dx_dx, -c_dt_dx, mom, (mom * vel) + p },
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     -a_dx_dx + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))), -a_dt_dx,
     -b_dx_dx - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))), -b_dt_dx,
     -c_dx_dx, -c_dt_dx, mom, (mom * vel) + p },
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

  double q_local[15], flux_local[15], flux[15];
  for (int d = 0; d < 3; d++) {
    gr_medium->rotate_to_local_func(gr_medium, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_gr_medium_flux(gas_gamma, kappa, q_local, flux_local);
    gr_medium->rotate_to_global_func(gr_medium, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-14) );
    }
  }

  double q_l[15], q_g[15];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(gr_medium, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-14) );
    }

    double w1[15], q1[15];
    gr_medium->cons_to_riem(gr_medium, q_local, q_local, w1);
    gr_medium->riem_to_cons(gr_medium, q_local, w1, q1);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(gr_medium);
}

void
test_gr_medium_waves_ho()
{
  double gas_gamma = 5.0 / 3.0;
  double kappa = 8.0 * M_PI;
  struct gkyl_wv_eqn *gr_medium = gkyl_wv_gr_medium_new(gas_gamma, kappa, false);

  double exp_2a_l = 0.1;
  double a_dt_l = 0.2, a_dx_l = 0.3;
  double b_dt_l = 0.4, b_dx_l = 0.5;
  double c_dt_l = 0.6, c_dx_l = 0.7;

  double a_dt_dx_l = 0.1, a_dx_dx_l = 0.2;
  double b_dt_dx_l = 0.3, b_dx_dx_l = 0.4;
  double c_dt_dx_l = 0.5, c_dx_dx_l = 0.6;

  double rho_l = 1.0, vel_l = 0.2;
  double p_l = (gas_gamma - 1.0) * rho_l;
  double W_l = 1.0 / sqrt(1.0 - (vel_l * vel_l));

  double ql[15] = { exp_2a_l, a_dt_l, a_dx_l, b_dt_l, b_dx_l, c_dt_l, c_dx_l, a_dt_dx_l, a_dx_dx_l, b_dt_dx_l, b_dx_dx_l, c_dt_dx_l, c_dx_dx_l,
    ((rho_l + p_l) * (W_l * W_l)) - p_l, (rho_l + p_l) * vel_l * (W_l * W_l) };

  double exp_2a_r = 0.07;
  double a_dt_r = 0.06, a_dx_r = 0.05;
  double b_dt_r = 0.04, b_dx_r = 0.03;
  double c_dt_r = 0.02, c_dx_r = 0.01;

  double a_dt_dx_r = 0.06, a_dx_dx_r = 0.05;
  double b_dt_dx_r = 0.04, b_dx_dx_r = 0.03;
  double c_dt_dx_r = 0.02, c_dx_dx_r = 0.01;

  double rho_r = 0.1, vel_r = 0.4;
  double p_r = (gas_gamma - 1.0) * rho_r;
  double W_r = 1.0 / sqrt(1.0 - (vel_r * vel_r));

  double qr[15] = { exp_2a_r, a_dt_r, a_dx_r, b_dt_r, b_dx_r, c_dt_r, c_dx_r, a_dt_dx_r, a_dx_dx_r, b_dt_dx_r, b_dx_dx_r, c_dt_dx_r, c_dx_dx_r,
    ((rho_r + p_r) * (W_r * W_r)) - p_r, (rho_r + p_r) * vel_r * (W_r * W_r) };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };

  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 15], waves_local[2 * 15];

    double ql_local[15], qr_local[15];
    gkyl_wv_eqn_rotate_to_local(gr_medium, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(gr_medium, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[15];
    for (int i = 0; i < 15; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(gr_medium, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[15], amdq_local[15];
    gkyl_wv_eqn_qfluct(gr_medium, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], &waves_local[i * 15], &waves[i * 15]);
    }

    double apdq[15], amdq[15];
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[15], fr_local[15];
    gkyl_gr_medium_flux(gas_gamma, kappa, ql_local, fl_local);
    gkyl_gr_medium_flux(gas_gamma, kappa, qr_local, fr_local);

    double fl[15], fr[15];
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(gr_medium);
}

void
test_gr_medium_waves_2_ho()
{
  double gas_gamma = 5.0 / 3.0;
  double kappa = 8.0 * M_PI;
  struct gkyl_wv_eqn *gr_medium = gkyl_wv_gr_medium_new(gas_gamma, kappa, false);

  double exp_2a_l = 1.0;
  double a_dt_l = 2.0, a_dx_l = 3.0;
  double b_dt_l = 4.0, b_dx_l = 5.0;
  double c_dt_l = 6.0, c_dx_l = 7.0;

  double a_dt_dx_l = 1.0, a_dx_dx_l = 2.0;
  double b_dt_dx_l = 3.0, b_dx_dx_l = 4.0;
  double c_dt_dx_l = 5.0, c_dx_dx_l = 6.0;

  double rho_l = 10.0, vel_l = 0.99;
  double p_l = (gas_gamma - 1.0) * rho_l;
  double W_l = 1.0 / sqrt(1.0 - (vel_l * vel_l));

  double ql[15] = { exp_2a_l, a_dt_l, a_dx_l, b_dt_l, b_dx_l, c_dt_l, c_dx_l, a_dt_dx_l, a_dx_dx_l, b_dt_dx_l, b_dx_dx_l, c_dt_dx_l, c_dx_dx_l,
    ((rho_l + p_l) * (W_l * W_l)) - p_l, (rho_l + p_l) * vel_l * (W_l * W_l) };

  double exp_2a_r = 0.07;
  double a_dt_r = 0.06, a_dx_r = 0.05;
  double b_dt_r = 0.04, b_dx_r = 0.03;
  double c_dt_r = 0.02, c_dx_r = 0.01;

  double a_dt_dx_r = 0.06, a_dx_dx_r = 0.05;
  double b_dt_dx_r = 0.04, b_dx_dx_r = 0.03;
  double c_dt_dx_r = 0.02, c_dx_dx_r = 0.01;

  double rho_r = 0.01, vel_r = 0.3;
  double p_r = (gas_gamma - 1.0) * rho_r;
  double W_r = 1.0 / sqrt(1.0 - (vel_r * vel_r));

  double qr[15] = { exp_2a_r, a_dt_r, a_dx_r, b_dt_r, b_dx_r, c_dt_r, c_dx_r, a_dt_dx_r, a_dx_dx_r, b_dt_dx_r, b_dx_dx_r, c_dt_dx_r, c_dx_dx_r,
    ((rho_r + p_r) * (W_r * W_r)) - p_r, (rho_r + p_r) * vel_r * (W_r * W_r) };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };

  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 15], waves_local[2 * 15];

    double ql_local[15], qr_local[15];
    gkyl_wv_eqn_rotate_to_local(gr_medium, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(gr_medium, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[15];
    for (int i = 0; i < 15; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(gr_medium, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[15], amdq_local[15];
    gkyl_wv_eqn_qfluct(gr_medium, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], &waves_local[i * 15], &waves[i * 15]);
    }

    double apdq[15], amdq[15];
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[15], fr_local[15];
    gkyl_gr_medium_flux(gas_gamma, kappa, ql_local, fl_local);
    gkyl_gr_medium_flux(gas_gamma, kappa, qr_local, fr_local);

    double fl[15], fr[15];
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(gr_medium, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(gr_medium);
}

TEST_LIST = {
  { "gr_medium_basic_ho", test_gr_medium_basic_ho },
  { "gr_medium_waves_ho", test_gr_medium_waves_ho },
  { "gr_medium_waves_2_ho", test_gr_medium_waves_2_ho },
  { NULL, NULL },
};
