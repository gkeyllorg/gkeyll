#include <acutest.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_iso_euler_priv.h>

void
calcq(const double pv[4], double q[4])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
}

void
test_iso_euler_basic_ho()
{
  double vt = 1.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt, false);

  TEST_CHECK( iso_euler->num_equations == 4 );
  TEST_CHECK( iso_euler->num_waves == 3 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double q[4], q_local[4], pv[4] = { rho, u, v, w };
  calcq(pv, q);

  double fluxes[3][4] = {
    { rho*u, rho*(u*u+vt*vt), rho*u*v, rho*u*w },
    { rho*v, rho*u*v, rho*(v*v+vt*vt), rho*v*w },
    { rho*w, rho*u*w, rho*v*w, rho*(w*w+vt*vt) },
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

  double flux[4], flux_local[4];  
  for (int d=1; d<2; ++d) {
    iso_euler->rotate_to_local_func(iso_euler, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_flux(vt, q_local, flux_local);
    iso_euler->rotate_to_global_func(iso_euler, tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<4; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) );

    // check Riemann transform
    double w1[4], q1[4];
    iso_euler->cons_to_riem(iso_euler, q_local, q_local, w1);
    iso_euler->riem_to_cons(iso_euler, q_local, w1, q1);
    
    for (int m=0; m<4; ++m)
      TEST_CHECK( gkyl_compare_double(q_local[m], q1[m], 1e-14) );
  }  

  iso_euler->rotate_to_local_func(iso_euler, tau1[0], tau2[0], norm[0], q, q_local);
  gkyl_iso_euler_flux(vt, q_local, flux_local);
  iso_euler->rotate_to_global_func(iso_euler, tau1[0], tau2[0], norm[0], flux_local, flux);

  TEST_CHECK( flux[0] == rho*u );
  TEST_CHECK( flux[1] == rho*(u*u + vt*vt) );
  TEST_CHECK( flux[2] == rho*u*v );
  TEST_CHECK( flux[3] == rho*u*w );

  iso_euler->rotate_to_local_func(iso_euler, tau1[1], tau2[1], norm[1], q, q_local);
  gkyl_iso_euler_flux(vt, q_local, flux_local);
  iso_euler->rotate_to_global_func(iso_euler, tau1[1], tau2[1], norm[1], flux_local, flux);
  
  TEST_CHECK( flux[0] == rho*v );
  TEST_CHECK( flux[1] == rho*v*u );
  TEST_CHECK( flux[2] == rho*(v*v + vt*vt) );
  TEST_CHECK( flux[3] == rho*v*w );

  iso_euler->rotate_to_local_func(iso_euler, tau1[2], tau2[2], norm[2], q, q_local);
  gkyl_iso_euler_flux(vt, q_local, flux_local);
  iso_euler->rotate_to_global_func(iso_euler, tau1[2], tau2[2], norm[2], flux_local, flux);

  TEST_CHECK( flux[0] == rho*w );
  TEST_CHECK( flux[1] == rho*w*u );
  TEST_CHECK( flux[2] == rho*w*v );
  TEST_CHECK( flux[3] == rho*(w*w + vt*vt) );
  
  gkyl_wv_eqn_release(iso_euler);
}

void
test_iso_euler_waves_ho()
{
  double vt = 1.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt, false);

  double vl[4] = { 1.0, 0.1, 0.2, 0.3};
  double vr[4] = { 0.1, 1.0, 2.0, 3.0};

  double ql[4], qr[4];
  double ql_local[4], qr_local[4];
  calcq(vl, ql); calcq(vr, qr);

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

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*4], waves_local[3*4];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(iso_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(iso_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[4];
    for (int i=0; i<4; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(iso_euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(iso_euler, tau1[d], tau2[d], norm[d], &waves_local[mw*4], &waves[mw*4]);
    
    double apdq[4], amdq[4];
    gkyl_wv_eqn_qfluct(iso_euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[4], fr_local[4];
    gkyl_iso_euler_flux(vt, ql_local, fl_local);
    gkyl_iso_euler_flux(vt, qr_local, fr_local);
    
    double fl[4], fr[4];
    gkyl_wv_eqn_rotate_to_global(iso_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(iso_euler, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(iso_euler);
}

void
test_iso_euler_waves_2_ho()
{
  double vt = 10.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt, false);

  double vl[4] = { 1.0, 0.1, 0.2, 0.3};
  double vr[4] = { 0.01, 1.0, 2.0, 3.0};

  double ql[4], qr[4];
  double ql_local[4], qr_local[4];
  calcq(vl, ql); calcq(vr, qr);

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

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*4], waves_local[3*4];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(iso_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(iso_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[4];
    for (int i=0; i<4; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(iso_euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(iso_euler, tau1[d], tau2[d], norm[d], &waves_local[mw*4], &waves[mw*4]);

    double apdq[4], amdq[4];
    gkyl_wv_eqn_qfluct(iso_euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[4], fr_local[4];
    gkyl_iso_euler_flux(vt, ql_local, fl_local);
    gkyl_iso_euler_flux(vt, qr_local, fr_local);
    
    double fl[4], fr[4];
    gkyl_wv_eqn_rotate_to_global(iso_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(iso_euler, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
  
  gkyl_wv_eqn_release(iso_euler);
}

TEST_LIST = {
  { "iso_euler_basic_ho", test_iso_euler_basic_ho },
  { "iso_euler_waves_ho", test_iso_euler_waves_ho },
  { "iso_euler_waves_2_ho", test_iso_euler_waves_2_ho },
  { NULL, NULL },
};
