#include <acutest.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_ten_moment_priv.h>

static const int dir_u_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

static const int dir_p_shuffle[][6] = {
  {4, 5, 6, 7, 8, 9},
  {7, 8, 5, 9, 6, 4},
  {9, 6, 8, 4, 5, 7}
};

// Make indexing cleaner with the dir_shuffle
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

#define PXX dp[0]
#define PXY dp[1]
#define PXZ dp[2]
#define PYY dp[3]
#define PYZ dp[4]
#define PZZ dp[5]

void
calcq(const double pv[10], double q[10])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3];
  double pxx = pv[4], pxy = pv[5], pxz = pv[6], pyy = pv[7], pyz = pv[8], pzz = pv[9];

  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pxx + rho*u*u; q[5] = pxy + rho*u*v; q[6] = pxz + rho*u*w;
  q[7] = pyy + rho*v*v; q[8] = pyz + rho*v*w; q[9] = pzz + rho*w*w;
}

void
test_ten_moment_basic_ho()
{
  struct gkyl_wv_eqn *ten_moment = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  TEST_CHECK( ten_moment->num_equations == 10 );
  TEST_CHECK( ten_moment->num_waves == 5 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double pxx = 0.5, pxy = 0.1, pxz = 0.2, pyy = 1.0, pyz = 0.3, pzz = 1.5;
  double q[10], pv[10] = { rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz };
  calcq(pv, q);

  // new array for re-computed primitive variables
  double var[10];
  gkyl_ten_moment_primitive(q, var);
  // check implementation of primitive variable calculation
  TEST_CHECK ( var[0] == pv[0] );
  TEST_CHECK ( var[1] == pv[1] );
  TEST_CHECK ( var[2] == pv[2] );
  TEST_CHECK ( var[3] == pv[3] );
  TEST_CHECK ( var[4] == pv[4] );
  TEST_CHECK ( var[5] == pv[5] );
  TEST_CHECK ( var[6] == pv[6] );
  TEST_CHECK ( var[7] == pv[7] );
  TEST_CHECK ( var[8] == pv[8] );
  TEST_CHECK ( var[9] == pv[9] );

  double fluxes[3][10] = {
    { rho*u, rho*u*u + pxx, rho*u*v + pxy, rho*u*w + pxz, 
      rho*u*u*u + 3*u*pxx, rho*u*u*v + 2*u*pxy + v*pxx, rho*u*u*w + 2*u*pxz + w*pxx,
      rho*u*v*v + 2*v*pxy + u*pyy, rho*u*v*w + u*pyz + v*pxz + w*pxy, rho*u*w*w + 2*w*pxz + u*pzz },
    { rho*v, rho*u*v + pxy, rho*v*v + pyy, rho*v*w + pyz, 
      rho*v*u*u + 2*u*pxy + v*pxx, rho*u*v*v + 2*v*pxy + u*pyy, rho*u*v*w + u*pyz + v*pxz + w*pxy,
      rho*v*v*v + 3*v*pyy, rho*v*v*w + 2*v*pyz + w*pyy, rho*v*w*w + 2*w*pyz + v*pzz},
    { rho*w, rho*u*w + pxz, rho*v*w + pyz, rho*w*w + pzz, 
      rho*u*u*w + 2*u*pxz + w*pxx, rho*u*v*w + u*pyz + v*pxz + w*pxy, rho*u*w*w + 2*w*pxz + u*pzz,
      rho*v*v*w + 2*v*pyz + w*pyy, rho*v*w*w + 2*w*pyz + v*pzz, rho*w*w*w + 3*w*pzz },
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

  double q_local[10], flux_local[10], flux[10];
  for (int d=0; d<3; ++d) {
    ten_moment->rotate_to_local_func(ten_moment, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_ten_moment_flux(q_local, flux_local);
    ten_moment->rotate_to_global_func(ten_moment, tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<10; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) );
  }

  double q_l[10], q_g[10];
  for (int d=0; d<3; ++d) {
    gkyl_wv_eqn_rotate_to_local(ten_moment, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(ten_moment, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<10; ++m) TEST_CHECK( q[m] == q_g[m] );

    // check Riemann transform
    double w1[10], q1[10];
    ten_moment->cons_to_riem(ten_moment, q_local, q_local, w1);
    ten_moment->riem_to_cons(ten_moment, q_local, w1, q1);
    
    for (int m=0; m<10; ++m)
      TEST_CHECK( gkyl_compare_double(q_local[m], q1[m], 1e-14) );
  }
  
  gkyl_wv_eqn_release(ten_moment);
}

void
test_ten_moment_waves_ho()
{
  struct gkyl_wv_eqn *ten_moment = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double vl[10] = { 1.0, 0.1, 0.2, 0.3, 0.5, 0.0, 0.0, 1.0, 0.0, 1.5};
  double vr[10] = { 0.1, 1.0, 2.0, 3.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.3};

  double ql[10], qr[10];
  double ql_local[10], qr_local[10];
  calcq(vl, ql); calcq(vr, qr);

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
    double speeds[5], waves[5*10], waves_local[5*10];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(ten_moment, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(ten_moment, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[10];
    for (int i=0; i<10; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(ten_moment, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<5; ++mw)
      gkyl_wv_eqn_rotate_to_global(ten_moment, tau1[d], tau2[d], norm[d], &waves_local[mw*10], &waves[mw*10]);

    double apdq[10], amdq[10];
    gkyl_wv_eqn_qfluct(ten_moment, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[10], fr_local[10];
    gkyl_ten_moment_flux(ql_local, fl_local);
    gkyl_ten_moment_flux(qr_local, fr_local);

    double fl[10], fr[10];
    gkyl_wv_eqn_rotate_to_global(ten_moment, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(ten_moment, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<10; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(ten_moment);
}

#ifdef GKYL_HAVE_CUDA

int cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn);

void
test_wv_ten_moment_dev()
{
  double k0 = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(k0, false, false, 0, 0, true);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct wv_ten_moment *ten_moment = container_of(eqn, struct wv_ten_moment, eqn);

  TEST_CHECK( ten_moment->k0 == 1.0 ); 
  
  // call CUDA test
  int nfail = cu_wv_ten_moment_test(eqn->on_dev);

  TEST_CHECK( nfail == 0 );

  gkyl_wv_eqn_release(eqn);
}

#endif

TEST_LIST = {
  { "ten_moment_basic_ho", test_ten_moment_basic_ho },
  { "ten_moment_waves_ho", test_ten_moment_waves_ho },
#ifdef GKYL_HAVE_CUDA
  { "wv_ten_moment_dev", test_wv_ten_moment_dev },
#endif  
  { NULL, NULL },
};
