#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_euler_mixture.h>
#include <gkyl_wv_euler_mixture_priv.h>

void
test_euler_mixture_twocomponent_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(2, gas_gamma_s, false);

  TEST_CHECK( euler_mixture->num_equations == 8 );
  TEST_CHECK( euler_mixture->num_waves == 2 );

  double alpha1 = 0.75, rho1 = 1.0, rho2 = 2.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E_total = (alpha1 * E1) + ((1.0 - alpha1) * E2);

  double q[8] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total, rho_total * alpha1,
    alpha1 * rho1, (1.0 - alpha1) * rho2 };

  double prims[8];
  gkyl_euler_mixture_prim_vars(2, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed mixture pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical mixture rules into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho2, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed mixture pressure.

  double fluxes[3][8] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + p_total, rho_total * (vx_total * vy_total), rho_total * (vx_total * vz_total),
     (E_total * vx_total) + (vx_total * p_total), rho_total * vx_total * alpha1, alpha1 * (vx_total * rho1), (1.0 - alpha1) * (vx_total * rho2) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + p_total, rho_total * (vy_total * vz_total),
     (E_total * vy_total) + (vy_total * p_total), rho_total * vy_total * alpha1, alpha1 * (vy_total * rho1), (1.0 - alpha1) * (vy_total * rho2) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total), (rho_total * (vz_total * vz_total)) + p_total,
     (E_total * vz_total) + (vz_total * p_total), rho_total * vz_total * alpha1, alpha1 * (vz_total * rho1), (1.0 - alpha1) * (vz_total * rho2) },
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

  double q_local[8], flux_local[8], flux[8];
  for (int d = 0; d < 3; d++) {
    euler_mixture->rotate_to_local_func(euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_mixture_flux(2, gas_gamma_s, q_local, flux_local);
    euler_mixture->rotate_to_global_func(euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[8], q_g[8];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[8], q1[8];
    euler_mixture->cons_to_riem(euler_mixture, q_local, q_local, w1);
    euler_mixture->riem_to_cons(euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

void
test_euler_mixture_threecomponent_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(3, gas_gamma_s, false);

  TEST_CHECK( euler_mixture->num_equations == 10 );
  TEST_CHECK( euler_mixture->num_waves == 2 );

  double alpha1 = 0.5, alpha2 = 0.3, rho1 = 1.0, rho2 = 2.0, rho3 = 3.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double rho_total = (alpha1 * rho1) + (alpha2 * rho2) + ((1.0 - (alpha1 + alpha2)) * rho3);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E3 = (p_total / (gas_gamma3 - 1.0)) + (0.5 * rho3 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E_total = (alpha1 * E1) + (alpha2 * E2) + ((1.0 - (alpha1 + alpha2)) * rho3);

  double q[10] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total, rho_total * alpha1,
    rho_total * alpha2, alpha1 * rho1, alpha2 * rho2, (1.0 - (alpha1 + alpha2)) * rho3 };

  double prims[10];
  gkyl_euler_mixture_prim_vars(3, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed mixture pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical mixture rules into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], alpha2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[8], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[9], rho3, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed mixture pressure.

  double fluxes[3][10] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + p_total, rho_total * (vx_total * vy_total), rho_total * (vx_total * vz_total),
     (E_total * vx_total) + (vx_total * p_total), rho_total * vx_total * alpha1, rho_total * vx_total * alpha2, alpha1 * (vx_total * rho1),
     alpha2 * (vx_total * rho2), (1.0 - (alpha1 + alpha2)) * (vx_total * rho3) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + p_total, rho_total * (vy_total * vz_total),
     (E_total * vy_total) + (vy_total * p_total), rho_total * vy_total * alpha1, rho_total * vy_total * alpha2, alpha1 * (vy_total * rho1),
     alpha2 * (vy_total * rho2), (1.0 - (alpha1 + alpha2)) * (vy_total * rho3) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total), (rho_total * (vz_total * vz_total)) + p_total,
     (E_total * vz_total) + (vz_total * p_total), rho_total * vz_total * alpha1, rho_total * vz_total * alpha2, alpha1 * (vz_total * rho1),
     alpha2 * (vz_total * rho2), (1.0 - (alpha1 + alpha2)) * (vz_total * rho3) },
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

  double q_local[10], flux_local[10], flux[10];
  for (int d = 0; d < 3; d++) {
    euler_mixture->rotate_to_local_func(euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_mixture_flux(3, gas_gamma_s, q_local, flux_local);
    euler_mixture->rotate_to_global_func(euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
    }
  }

  double q_l[10], q_g[10];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[10], q1[10];
    euler_mixture->cons_to_riem(euler_mixture, q_local, q_local, w1);
    euler_mixture->riem_to_cons(euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

void
test_euler_mixture_twocomponent_waves()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(2, gas_gamma_s, false);
  
  double alpha1_l = 0.75, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (alpha1_l * rho1_l) + ((1.0 - alpha1_l) * rho2_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (alpha1_l * E1_l) + ((1.0 - alpha1_l) * E2_l);

  double alpha1_r = 0.25, rho1_r = 0.1, rho2_r = 0.2, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 0.15;
  double rho_total_r = (alpha1_r * rho1_r) + ((1.0 - alpha1_r) * rho2_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (alpha1_r * E1_r) + ((1.0 - alpha1_r) * E2_r);

  double ql[8] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * alpha1_l,
    alpha1_l * rho1_l, (1.0 - alpha1_l) * rho2_l };
  double qr[8] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * alpha1_r,
    alpha1_r * rho1_r, (1.0 - alpha1_r) * rho2_r };
  
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
    double speeds[2], waves[2 * 8], waves_local[2 * 8];

    double ql_local[8], qr_local[8];
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i = 0; i < 8; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[8], amdq_local[8];
    gkyl_wv_eqn_qfluct(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 8], &waves[i * 8]);
    }

    double apdq[8], amdq[8];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[8], fr_local[8];
    gkyl_euler_mixture_flux(2, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_mixture_flux(2, gas_gamma_s, qr_local, fr_local);

    double fl[8], fr[8];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

void
test_euler_mixture_twocomponent_waves_2()
{
  double gas_gamma1 = 1.2;
  double gas_gamma2 = 1.7;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(2, gas_gamma_s, false);
  
  double alpha1_l = 0.9, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (alpha1_l * rho1_l) + ((1.0 - alpha1_l) * rho2_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (alpha1_l * E1_l) + ((1.0 - alpha1_l) * E2_l);

  double alpha1_r = 0.1, rho1_r = 0.01, rho2_r = 0.02, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 15.0;
  double rho_total_r = (alpha1_r * rho1_r) + ((1.0 - alpha1_r) * rho2_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (alpha1_r * E1_r) + ((1.0 - alpha1_r) * E2_r);

  double ql[8] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * alpha1_l,
    alpha1_l * rho1_l, (1.0 - alpha1_l) * rho2_l };
  double qr[8] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * alpha1_r,
    alpha1_r * rho1_r, (1.0 - alpha1_r) * rho2_r };
  
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
    double speeds[2], waves[2 * 8], waves_local[2 * 8];

    double ql_local[8], qr_local[8];
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i = 0; i < 8; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[8], amdq_local[8];
    gkyl_wv_eqn_qfluct(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 8], &waves[i * 8]);
    }

    double apdq[8], amdq[8];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[8], fr_local[8];
    gkyl_euler_mixture_flux(2, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_mixture_flux(2, gas_gamma_s, qr_local, fr_local);

    double fl[8], fr[8];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 8; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

void
test_euler_mixture_threecomponent_waves()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(3, gas_gamma_s, false);
  
  double alpha1_l = 0.5, alpha2_l = 0.3, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (alpha1_l * rho1_l) + (alpha2_l * rho2_l) + ((1.0 - (alpha1_l + alpha2_l)) * rho3_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E3_l = (p_total_l / (gas_gamma3 - 1.0)) + (0.5 * rho3_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (alpha1_l * E1_l) + (alpha2_l * E2_l) + ((1.0 - (alpha1_l + alpha2_l)) * E3_l);

  double alpha1_r = 0.4, alpha2_r = 0.2, rho1_r = 0.1, rho2_r = 0.2, rho3_r = 0.3, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 0.15;
  double rho_total_r = (alpha1_r * rho1_r) + (alpha2_r * rho2_r) + ((1.0 - (alpha1_r + alpha2_r)) * rho3_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E3_r = (p_total_r / (gas_gamma3 - 1.0)) + (0.5 * rho3_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (alpha1_r * E1_r) + (alpha2_r * E2_r) + ((1.0 - (alpha1_r + alpha2_r)) * E3_r);

  double ql[10] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * alpha1_l,
    rho_total_l * alpha2_l, alpha1_l * rho1_l, alpha2_l * rho2_l, (1.0 - (alpha1_l + alpha2_l)) * rho3_l };
  double qr[10] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * alpha1_r,
    rho_total_r * alpha2_r, alpha1_r * rho1_r, alpha2_r * rho2_r, (1.0 - (alpha1_r + alpha2_r)) * rho3_r };
  
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
    double speeds[2], waves[2 * 10], waves_local[2 * 10];

    double ql_local[10], qr_local[10];
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[10];
    for (int i = 0; i < 10; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[10], amdq_local[10];
    gkyl_wv_eqn_qfluct(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 10], &waves[i * 10]);
    }

    double apdq[10], amdq[10];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[10], fr_local[10];
    gkyl_euler_mixture_flux(3, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_mixture_flux(3, gas_gamma_s, qr_local, fr_local);

    double fl[10], fr[10];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

void
test_euler_mixture_threecomponent_waves_2()
{
  double gas_gamma1 = 1.1;
  double gas_gamma2 = 1.5;
  double gas_gamma3 = 1.9;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(3, gas_gamma_s, false);
  
  double alpha1_l = 0.8, alpha2_l = 0.1, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (alpha1_l * rho1_l) + (alpha2_l * rho2_l) + ((1.0 - (alpha1_l + alpha2_l)) * rho3_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E3_l = (p_total_l / (gas_gamma3 - 1.0)) + (0.5 * rho3_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (alpha1_l * E1_l) + (alpha2_l * E2_l) + ((1.0 - (alpha1_l + alpha2_l)) * E3_l);

  double alpha1_r = 0.75, alpha2_r = 0.2, rho1_r = 0.01, rho2_r = 0.02, rho3_r = 0.03, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 15.0;
  double rho_total_r = (alpha1_r * rho1_r) + (alpha2_r * rho2_r) + ((1.0 - (alpha1_r + alpha2_r)) * rho3_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E3_r = (p_total_r / (gas_gamma3 - 1.0)) + (0.5 * rho3_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (alpha1_r * E1_r) + (alpha2_r * E2_r) + ((1.0 - (alpha1_r + alpha2_r)) * E3_r);

  double ql[10] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * alpha1_l,
    rho_total_l * alpha2_l, alpha1_l * rho1_l, alpha2_l * rho2_l, (1.0 - (alpha1_l + alpha2_l)) * rho3_l };
  double qr[10] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * alpha1_r,
    rho_total_r * alpha2_r, alpha1_r * rho1_r, alpha2_r * rho2_r, (1.0 - (alpha1_r + alpha2_r)) * rho3_r };
  
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
    double speeds[2], waves[2 * 10], waves_local[2 * 10];

    double ql_local[10], qr_local[10];
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[10];
    for (int i = 0; i < 10; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[10], amdq_local[10];
    gkyl_wv_eqn_qfluct(euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 10], &waves[i * 10]);
    }

    double apdq[10], amdq[10];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[10], fr_local[10];
    gkyl_euler_mixture_flux(3, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_mixture_flux(3, gas_gamma_s, qr_local, fr_local);

    double fl[10], fr[10];
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 10; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
    }
  }

  gkyl_wv_eqn_release(euler_mixture);
  gkyl_free(gas_gamma_s);
}

TEST_LIST = {
  { "euler_mixture_twocomponent_basic", test_euler_mixture_twocomponent_basic },
  { "euler_mixture_threecomponent_basic", test_euler_mixture_threecomponent_basic },
  { "euler_mixture_twocomponent_waves", test_euler_mixture_twocomponent_waves },
  { "euler_mixture_twocomponent_waves_2", test_euler_mixture_twocomponent_waves_2 },
  { "euler_mixture_threecomponent_waves", test_euler_mixture_threecomponent_waves },
  { "euler_mixture_threecomponent_waves_2", test_euler_mixture_threecomponent_waves_2 },
  { NULL, NULL },
};
