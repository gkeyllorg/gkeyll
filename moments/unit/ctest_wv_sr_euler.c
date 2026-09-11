#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_moment_prim_sr_euler.h>
#include <gkyl_wv_sr_euler.h>

void
calcq(double gas_gamma, const double pv[5], double q[5])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = gas_gamma * pr / (gas_gamma - 1)  + rho;
  q[0] = gamma*rho;
  q[1] = gamma*gamma*rhoh - pr;
  q[2] = gamma*gamma*rhoh*u;
  q[3] = gamma*gamma*rhoh*v;
  q[4] = gamma*gamma*rhoh*w;
  
}

void
test_sr_euler_prim1_ho()
{
  double gas_gamma = 1.333;
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(gas_gamma);

  TEST_CHECK( sr_euler->num_equations == 5 );
  TEST_CHECK( sr_euler->num_waves == 3 );

  double rho = 1.0, u = 0.9999, v = 0.0, w = 0.01, pr = 1.5;
  double q[5], pv2[5], pv[5] = { rho, u, v, w, pr };
  calcq(gas_gamma, pv, q);
  
  gkyl_sr_euler_prim_vars(gas_gamma, q, pv2);
  
  TEST_CHECK ( gkyl_compare(rho, pv2[0], 1e-10) );
  TEST_CHECK ( gkyl_compare(pr, pv2[1], 1e-10) );
  TEST_CHECK ( gkyl_compare(u, pv2[2], 1e-15) );
  TEST_CHECK ( gkyl_compare(v, pv2[3], 1e-15) );
  TEST_CHECK ( gkyl_compare(w, pv2[4], 1e-15) );

  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = gas_gamma * pr / (gas_gamma - 1)  + rho;

  double fluxes[3][5] = {
    { gamma*rho*u, gamma*gamma*rhoh*u, gamma*gamma*rhoh*u*u + pr, gamma*gamma*rhoh*u*v, gamma*gamma*rhoh*u*w },
    { gamma*rho*v, gamma*gamma*rhoh*v, gamma*gamma*rhoh*v*u, gamma*gamma*rhoh*v*v + pr, gamma*gamma*rhoh*v*w },
    { gamma*rho*w, gamma*gamma*rhoh*w, gamma*gamma*rhoh*w*u, gamma*gamma*rhoh*w*v, gamma*gamma*rhoh*w*w + pr },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 }
  };

  double q_local[5], flux_local[5], flux[5];
  for (int d=1; d<2; ++d) {
    sr_euler->rotate_to_local_func(sr_euler, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_sr_euler_flux(gas_gamma, q_local, flux_local);
    sr_euler->rotate_to_global_func(sr_euler, tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<5; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-10) );
  }

  double q_l[5], q_g[5];
  for (int d=1; d<3; ++d) {
    gkyl_wv_eqn_rotate_to_local(sr_euler, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<5; ++m)
      TEST_CHECK( gkyl_compare(q[m], q_g[m], 1e-12) );

    // check Riemann transform
    double w1[5], q1[5];
    sr_euler->cons_to_riem(sr_euler, q_local, q_local, w1);
    sr_euler->riem_to_cons(sr_euler, q_local, w1, q1);
    
    for (int m=0; m<5; ++m)
      TEST_CHECK( gkyl_compare_double(q_local[m], q1[m], 1e-14) );    
  }
  
  gkyl_wv_eqn_release(sr_euler);
  
}

void
test_sr_euler_waves_ho()
{
  double gas_gamma = 1.333;
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(gas_gamma);

  double vl[5] = { 1.0, 0.09, 0.02, 0.03, 1.5};
  double vr[5] = { 0.1, .9, 0.2, 0.3, 0.15};

  double ql[5], qr[5];
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*5], waves_local[3*5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(sr_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(sr_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[5];
    for (int i=0; i<5; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(sr_euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], &waves_local[mw*5], &waves[mw*5]);
    
    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(sr_euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);

    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[5], fr_local[5];
    gkyl_sr_euler_flux(gas_gamma, ql_local, fl_local);
    gkyl_sr_euler_flux(gas_gamma, qr_local, fr_local);

    double fl[5], fr[5];
    gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-12) );
  }
    
  gkyl_wv_eqn_release(sr_euler);
}

void
test_sr_euler_waves2_ho()
{
  double gas_gamma = 1.3333;
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(gas_gamma);

  double vl[5] = { 1.0, 0.0999, 0.02, 0.03, 1500.};
  double vr[5] = { 0.1, 0.999, 0.01, 0.0, 0.015};
  
  double ql[5], qr[5];
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*5], waves_local[3*5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(sr_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(sr_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[5];
    for (int i=0; i<5; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(sr_euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], &waves_local[mw*5], &waves[mw*5]);

    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(sr_euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[5], fr_local[5];
    gkyl_sr_euler_flux(gas_gamma, ql_local, fl_local);
    gkyl_sr_euler_flux(gas_gamma, qr_local, fr_local);

    double fl[5], fr[5];
    gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(sr_euler, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-10) );
  }

  gkyl_wv_eqn_release(sr_euler);
}


TEST_LIST = {
  { "sr_euler_prim1_ho", test_sr_euler_prim1_ho },
  { "test_sr_euler_waves_ho", test_sr_euler_waves_ho},
  { "test_sr_euler_waves2_ho", test_sr_euler_waves2_ho},
  { NULL, NULL },
};
