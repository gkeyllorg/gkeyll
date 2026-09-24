#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_euler_rgfm.h>
#include <gkyl_wv_euler_rgfm_priv.h>

void
test_euler_rgfm_twospecies_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_rgfm = gkyl_wv_euler_rgfm_new(2, gas_gamma_s, 0, false);

  TEST_CHECK( euler_rgfm->num_equations == 9 );
  TEST_CHECK( euler_rgfm->num_waves == 2 );

  double phi1 = 0.75, rho1 = 1.0, rho2 = 2.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double rho_total = (phi1 * rho1) + ((1.0 - phi1) * rho2);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E_total = (phi1 * E1) + ((1.0 - phi1) * E2);

  double q[9] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total, rho_total * phi1,
    phi1 * rho1, (1.0 - phi1) * rho2, 0.0 };

  double prims[9];
  gkyl_euler_rgfm_prim_vars(2, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed interface pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical boundary conditions into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], phi1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho2, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed interface pressure.

  double fluxes[3][9] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + p_total, rho_total * (vx_total * vy_total), rho_total * (vx_total * vz_total),
     (E_total * vx_total) + (vx_total * p_total), rho_total * vx_total * phi1, phi1 * (vx_total * rho1), (1.0 - phi1) * (vx_total * rho2), 0.0 },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + p_total, rho_total * (vy_total * vz_total),
     (E_total * vy_total) + (vy_total * p_total), rho_total * vy_total * phi1, phi1 * (vy_total * rho1), (1.0 - phi1) * (vy_total * rho2), 0.0 },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total), (rho_total * (vz_total * vz_total)) + p_total,
     (E_total * vz_total) + (vz_total * p_total), rho_total * vz_total * phi1, phi1 * (vz_total * rho1), (1.0 - phi1) * (vz_total * rho2), 0.0 },
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

  double q_local[9], flux_local[9], flux[9];
  for (int d = 0; d < 3; d++) {
    euler_rgfm->rotate_to_local_func(euler_rgfm, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_rgfm_flux(2, gas_gamma_s, q_local, flux_local);
    euler_rgfm->rotate_to_global_func(euler_rgfm, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[9], q_g[9];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[9], q1[9];
    euler_rgfm->cons_to_riem(euler_rgfm, q_local, q_local, w1);
    euler_rgfm->riem_to_cons(euler_rgfm, q_local, w1, q1);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(euler_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_euler_rgfm_threespecies_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_rgfm = gkyl_wv_euler_rgfm_new(3, gas_gamma_s, 0, false);

  TEST_CHECK( euler_rgfm->num_equations == 11 );
  TEST_CHECK( euler_rgfm->num_waves == 2 );

  double phi1 = 0.5, phi2 = 0.3, rho1 = 1.0, rho2 = 2.0, rho3 = 3.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double rho_total = (phi1 * rho1) + (phi2 * rho2) + ((1.0 - (phi1 + phi2)) * rho3);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E3 = (p_total / (gas_gamma3 - 1.0)) + (0.5 * rho3 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total)));
  double E_total = (phi1 * E1) + (phi2 * E2) + ((1.0 - (phi1 + phi2)) * rho3);

  double q[11] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total, rho_total * phi1,
    rho_total * phi2, phi1 * rho1, phi2 * rho2, (1.0 - (phi1 + phi2)) * rho3, 0.0 };

  double prims[11];
  gkyl_euler_rgfm_prim_vars(3, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed rgfm pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical rgfm rules into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], phi1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], phi2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[8], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[9], rho3, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed rgfm pressure.

  double fluxes[3][11] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + p_total, rho_total * (vx_total * vy_total), rho_total * (vx_total * vz_total),
     (E_total * vx_total) + (vx_total * p_total), rho_total * vx_total * phi1, rho_total * vx_total * phi2, phi1 * (vx_total * rho1),
     phi2 * (vx_total * rho2), (1.0 - (phi1 + phi2)) * (vx_total * rho3), 0.0 },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + p_total, rho_total * (vy_total * vz_total),
     (E_total * vy_total) + (vy_total * p_total), rho_total * vy_total * phi1, rho_total * vy_total * phi2, phi1 * (vy_total * rho1),
     phi2 * (vy_total * rho2), (1.0 - (phi1 + phi2)) * (vy_total * rho3), 0.0 },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total), (rho_total * (vz_total * vz_total)) + p_total,
     (E_total * vz_total) + (vz_total * p_total), rho_total * vz_total * phi1, rho_total * vz_total * phi2, phi1 * (vz_total * rho1),
     phi2 * (vz_total * rho2), (1.0 - (phi1 + phi2)) * (vz_total * rho3), 0.0 },
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

  double q_local[11], flux_local[11], flux[11];
  for (int d = 0; d < 3; d++) {
    euler_rgfm->rotate_to_local_func(euler_rgfm, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_rgfm_flux(3, gas_gamma_s, q_local, flux_local);
    euler_rgfm->rotate_to_global_func(euler_rgfm, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 11; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
    }
  }

  double q_l[11], q_g[11];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 11; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[11], q1[11];
    euler_rgfm->cons_to_riem(euler_rgfm, q_local, q_local, w1);
    euler_rgfm->riem_to_cons(euler_rgfm, q_local, w1, q1);

    for (int i = 0; i < 11; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(euler_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_euler_rgfm_twospecies_waves()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_rgfm = gkyl_wv_euler_rgfm_new(2, gas_gamma_s, 0, false);
  
  double phi1_l = 0.75, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (phi1_l * rho1_l) + ((1.0 - phi1_l) * rho2_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (phi1_l * E1_l) + ((1.0 - phi1_l) * E2_l);

  double phi1_r = 0.25, rho1_r = 0.1, rho2_r = 0.2, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 0.15;
  double rho_total_r = (phi1_r * rho1_r) + ((1.0 - phi1_r) * rho2_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (phi1_r * E1_r) + ((1.0 - phi1_r) * E2_r);

  double ql[9] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * phi1_l,
    phi1_l * rho1_l, (1.0 - phi1_l) * rho2_l, 0.0 };
  double qr[9] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * phi1_r,
    phi1_r * rho1_r, (1.0 - phi1_r) * rho2_r, 0.0 };
  
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
    double speeds[2], waves[2 * 9], waves_local[2 * 9];

    double ql_local[9], qr_local[9];
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i = 0; i < 9; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 8], &waves[i * 8]);
    }

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[9], fr_local[9];
    gkyl_euler_rgfm_flux(2, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_rgfm_flux(2, gas_gamma_s, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(euler_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_euler_rgfm_twospecies_waves_2()
{
  double gas_gamma1 = 1.2;
  double gas_gamma2 = 1.7;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *euler_rgfm = gkyl_wv_euler_rgfm_new(2, gas_gamma_s, 0, false);
  
  double phi1_l = 0.9, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (phi1_l * rho1_l) + ((1.0 - phi1_l) * rho2_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (phi1_l * E1_l) + ((1.0 - phi1_l) * E2_l);

  double phi1_r = 0.1, rho1_r = 0.01, rho2_r = 0.02, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 15.0;
  double rho_total_r = (phi1_r * rho1_r) + ((1.0 - phi1_r) * rho2_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (phi1_r * E1_r) + ((1.0 - phi1_r) * E2_r);

  double ql[9] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * phi1_l,
    phi1_l * rho1_l, (1.0 - phi1_l) * rho2_l, 0.0 };
  double qr[9] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * phi1_r,
    phi1_r * rho1_r, (1.0 - phi1_r) * rho2_r, 0.0 };
  
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
    double speeds[2], waves[2 * 9], waves_local[2 * 9];

    double ql_local[9], qr_local[9];
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i = 0; i < 9; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 8], &waves[i * 8]);
    }

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[9], fr_local[9];
    gkyl_euler_rgfm_flux(2, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_rgfm_flux(2, gas_gamma_s, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(euler_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_euler_rgfm_threespecies_waves()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_rgfm = gkyl_wv_euler_rgfm_new(3, gas_gamma_s, 0, false);
  
  double phi1_l = 0.5, phi2_l = 0.3, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (phi1_l * rho1_l) + (phi2_l * rho2_l) + ((1.0 - (phi1_l + phi2_l)) * rho3_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E3_l = (p_total_l / (gas_gamma3 - 1.0)) + (0.5 * rho3_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (phi1_l * E1_l) + (phi2_l * E2_l) + ((1.0 - (phi1_l + phi2_l)) * E3_l);

  double phi1_r = 0.4, phi2_r = 0.2, rho1_r = 0.1, rho2_r = 0.2, rho3_r = 0.3, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 0.15;
  double rho_total_r = (phi1_r * rho1_r) + (phi2_r * rho2_r) + ((1.0 - (phi1_r + phi2_r)) * rho3_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E3_r = (p_total_r / (gas_gamma3 - 1.0)) + (0.5 * rho3_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (phi1_r * E1_r) + (phi2_r * E2_r) + ((1.0 - (phi1_r + phi2_r)) * E3_r);

  double ql[11] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * phi1_l,
    rho_total_l * phi2_l, phi1_l * rho1_l, phi2_l * rho2_l, (1.0 - (phi1_l + phi2_l)) * rho3_l, 0.0 };
  double qr[11] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * phi1_r,
    rho_total_r * phi2_r, phi1_r * rho1_r, phi2_r * rho2_r, (1.0 - (phi1_r + phi2_r)) * rho3_r, 0.0 };
  
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
    double speeds[2], waves[2 * 11], waves_local[2 * 11];

    double ql_local[11], qr_local[11];
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[11];
    for (int i = 0; i < 11; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[11], amdq_local[11];
    gkyl_wv_eqn_qfluct(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 10], &waves[i * 10]);
    }

    double apdq[11], amdq[11];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[11], fr_local[11];
    gkyl_euler_rgfm_flux(3, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_rgfm_flux(3, gas_gamma_s, qr_local, fr_local);

    double fl[11], fr[11];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 11; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
    }
  }

  gkyl_wv_eqn_release(euler_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_euler_rgfm_threespecies_waves_2()
{
  double gas_gamma1 = 1.1;
  double gas_gamma2 = 1.5;
  double gas_gamma3 = 1.9;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *euler_rgfm = gkyl_wv_euler_rgfm_new(3, gas_gamma_s, 0, false);
  
  double phi1_l = 0.8, phi2_l = 0.1, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double rho_total_l = (phi1_l * rho1_l) + (phi2_l * rho2_l) + ((1.0 - (phi1_l + phi2_l)) * rho3_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E3_l = (p_total_l / (gas_gamma3 - 1.0)) + (0.5 * rho3_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l)));
  double E_total_l = (phi1_l * E1_l) + (phi2_l * E2_l) + ((1.0 - (phi1_l + phi2_l)) * E3_l);

  double phi1_r = 0.75, phi2_r = 0.2, rho1_r = 0.01, rho2_r = 0.02, rho3_r = 0.03, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 15.0;
  double rho_total_r = (phi1_r * rho1_r) + (phi2_r * rho2_r) + ((1.0 - (phi1_r + phi2_r)) * rho3_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E3_r = (p_total_r / (gas_gamma3 - 1.0)) + (0.5 * rho3_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r)));
  double E_total_r = (phi1_r * E1_r) + (phi2_r * E2_r) + ((1.0 - (phi1_r + phi2_r)) * E3_r);

  double ql[11] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l, rho_total_l * phi1_l,
    rho_total_l * phi2_l, phi1_l * rho1_l, phi2_l * rho2_l, (1.0 - (phi1_l + phi2_l)) * rho3_l, 0.0 };
  double qr[11] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r, rho_total_r * phi1_r,
    rho_total_r * phi2_r, phi1_r * rho1_r, phi2_r * rho2_r, (1.0 - (phi1_r + phi2_r)) * rho3_r, 0.0 };
  
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
    double speeds[2], waves[2 * 11], waves_local[2 * 11];

    double ql_local[11], qr_local[11];
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[11];
    for (int i = 0; i < 11; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[11], amdq_local[11];
    gkyl_wv_eqn_qfluct(euler_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 10], &waves[i * 10]);
    }

    double apdq[11], amdq[11];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[11], fr_local[11];
    gkyl_euler_rgfm_flux(3, gas_gamma_s, ql_local, fl_local);
    gkyl_euler_rgfm_flux(3, gas_gamma_s, qr_local, fr_local);

    double fl[11], fr[11];
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 11; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
    }
  }

  gkyl_wv_eqn_release(euler_rgfm);
  gkyl_free(gas_gamma_s);
}

TEST_LIST = {
  { "euler_rgfm_twospecies_basic", test_euler_rgfm_twospecies_basic },
  { "euler_rgfm_threespecies_basic", test_euler_rgfm_threespecies_basic },
  { "euler_rgfm_twospecies_waves", test_euler_rgfm_twospecies_waves },
  { "euler_rgfm_twospecies_waves_2", test_euler_rgfm_twospecies_waves_2 },
  { "euler_rgfm_threespecies_waves", test_euler_rgfm_threespecies_waves },
  { "euler_rgfm_threespecies_waves_2", test_euler_rgfm_threespecies_waves_2 },
  { NULL, NULL },
};
