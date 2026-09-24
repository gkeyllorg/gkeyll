#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_reactive_euler.h>
#include <gkyl_wv_reactive_euler_priv.h>

void
test_reactive_euler_basic()
{
  double gas_gamma = 1.4;
  double specific_heat_capacity = 2.5;
  double energy_of_formation = 1.0;
  double ignition_temperature = 0.25;
  double reaction_rate = 250.0;
  struct gkyl_wv_eqn *reactive_euler = gkyl_wv_reactive_euler_new(gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature,
    reaction_rate, false);

  TEST_CHECK( reactive_euler->num_equations == 6 );
  TEST_CHECK( reactive_euler->num_waves == 2 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5, reac = 0.5;
  double q[6] = { rho, rho * u, rho * v, rho * w, p / (gas_gamma - 1.0) + (0.5 * rho * ((u * u) + (v * v) + (w * w))) +
    (energy_of_formation * (reac - 1.0)), rho * reac };

  double prims[6];
  gkyl_reactive_euler_prim_vars(gas_gamma, energy_of_formation, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], u, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], v, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], w, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[4], p, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[5], reac, 1e-16) );

  double E = p / (gas_gamma - 1.0) + (0.5 * rho * ((u * u) + (v * v) + (w * w))) + (energy_of_formation * (reac - 1.0));

  double fluxes[3][6] = {
   { rho * u, (rho * (u * u)) + p, rho * (u * v), rho * (u * w), (E + p) * u, rho * (u * reac) },
   { rho * v, rho * (v * u), (rho * (v * v)) + p, rho * (v * w), (E + p) * v, rho * (v * reac) },
   { rho * w, rho * (w * u), rho * (w * v), (rho * (w * w)) + p, (E + p) * w, rho * (w * reac) },
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

  double q_local[6], flux_local[6], flux[6];
  for (int d = 0; d < 3; d++) {
    reactive_euler->rotate_to_local_func(reactive_euler, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, q_local, flux_local);
    reactive_euler->rotate_to_global_func(reactive_euler, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 6; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
    }
  }

  double q_l[6], q_g[6];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(reactive_euler, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 6; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[6], q1[6];
    reactive_euler->cons_to_riem(reactive_euler, q_local, q_local, w1);
    reactive_euler->riem_to_cons(reactive_euler, q_local, w1, q1);

    for (int i = 0; i < 6; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(reactive_euler);
}

void
test_reactive_euler_waves()
{
  double gas_gamma = 1.4;
  double specific_heat_capacity = 2.5;
  double energy_of_formation = 1.0;
  double ignition_temperature = 0.25;
  double reaction_rate = 250.0;
  struct gkyl_wv_eqn *reactive_euler = gkyl_wv_reactive_euler_new(gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature,
    reaction_rate, false);
  
  double rhol = 1.0, ul = 0.1, vl = 0.2, wl = 0.3, pl = 1.5, reacl = 0.75;
  double rhor = 0.1, ur = 1.0, vr = 2.0, wr = 3.0, pr = 0.15, reacr = 0.25;

  double ql[6] = { rhol, rhol * ul, rhol * vl, rhol * wl, pl / (gas_gamma - 1.0) + (0.5 * rhol * ((ul * ul) + (vl * vl) + (wl * wl))) +
    (energy_of_formation * (reacl - 1.0)), rhol * reacl };
  double qr[6] = { rhor, rhor * ur, rhor * vr, rhor * wr, pr / (gas_gamma - 1.0) + (0.5 * rhor * ((ur * ur) + (vr * vr) + (wr * wr))) +
    (energy_of_formation * (reacr - 1.0)), rhor * reacr };
  
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
    double speeds[2], waves[2 * 6], waves_local[2 * 6];

    double ql_local[6], qr_local[6];
    gkyl_wv_eqn_rotate_to_local(reactive_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(reactive_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[6];
    for (int i = 0; i < 6; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(reactive_euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[6], amdq_local[6];
    gkyl_wv_eqn_qfluct(reactive_euler, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], &waves_local[i * 6], &waves[i * 6]);
    }

    double apdq[6], amdq[6];
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[6], fr_local[6];
    gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, ql_local, fl_local);
    gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, qr_local, fr_local);

    double fl[6], fr[6];
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 6; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-15) );
    }
  }

  gkyl_wv_eqn_release(reactive_euler);
}

void
test_reactive_euler_waves_2()
{
  double gas_gamma = 1.4;
  double specific_heat_capacity = 25.0;
  double energy_of_formation = 10.0;
  double ignition_temperature = 0.025;
  double reaction_rate = 2500.0;
  struct gkyl_wv_eqn *reactive_euler = gkyl_wv_reactive_euler_new(gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature,
    reaction_rate, false);
  
  double rhol = 1.0, ul = 0.1, vl = 0.2, wl = 0.3, pl = 1.5, reacl = 0.9;
  double rhor = 0.01, ur = 1.0, vr = 2.0, wr = 3.0, pr = 15.0, reacr = 0.1;

  double ql[6] = { rhol, rhol * ul, rhol * vl, rhol * wl, pl / (gas_gamma - 1.0) + (0.5 * rhol * ((ul * ul) + (vl * vl) + (wl * wl))) +
    (energy_of_formation * (reacl - 1.0)), rhol * reacl };
  double qr[6] = { rhor, rhor * ur, rhor * vr, rhor * wr, pr / (gas_gamma - 1.0) + (0.5 * rhor * ((ur * ur) + (vr * vr) + (wr * wr))) +
    (energy_of_formation * (reacr - 1.0)), rhor * reacr };
  
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
    double speeds[2], waves[2 * 6], waves_local[2 * 6];

    double ql_local[6], qr_local[6];
    gkyl_wv_eqn_rotate_to_local(reactive_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(reactive_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[6];
    for (int i = 0; i < 6; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(reactive_euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[6], amdq_local[6];
    gkyl_wv_eqn_qfluct(reactive_euler, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], &waves_local[i * 6], &waves[i * 6]);
    }

    double apdq[6], amdq[6];
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[6], fr_local[6];
    gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, ql_local, fl_local);
    gkyl_reactive_euler_flux(gas_gamma, energy_of_formation, qr_local, fr_local);

    double fl[6], fr[6];
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(reactive_euler, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 6; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(reactive_euler);
}

TEST_LIST = {
  { "reactive_euler_basic", test_reactive_euler_basic },
  { "reactive_euler_waves", test_reactive_euler_waves },
  { "reactive_euler_waves_2", test_reactive_euler_waves_2 },
  { NULL, NULL },
};
