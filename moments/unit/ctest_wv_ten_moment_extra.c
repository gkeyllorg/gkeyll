#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_wv_ten_moment.h>

// Build a 10-moment conserved state from primitives:
//   density rho, velocity (u,v,w), pressure tensor (Pxx,Pxy,Pxz,Pyy,Pyz,Pzz).
//   q = {rho, rho u, rho v, rho w,
//        Pxx + rho u^2, Pxy + rho u v, Pxz + rho u w,
//        Pyy + rho v^2, Pyz + rho v w, Pzz + rho w^2}
static void
calcq(double rho, double u, double v, double w,
  double pxx, double pxy, double pxz, double pyy, double pyz, double pzz,
  double q[10])
{
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pxx + rho*u*u;
  q[5] = pxy + rho*u*v;
  q[6] = pxz + rho*u*w;
  q[7] = pyy + rho*v*v;
  q[8] = pyz + rho*v*w;
  q[9] = pzz + rho*w*w;
}

// Structural properties and accessor parameters.
void
test_ten_moment_struct()
{
  double k0 = 5.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(k0, false, false, 1, 0, false);

  TEST_CHECK( eqn->type == GKYL_EQN_TEN_MOMENT );
  TEST_CHECK( eqn->num_equations == 10 );
  TEST_CHECK( eqn->num_waves == 5 );
  TEST_CHECK( eqn->num_diag == 10 );

  TEST_CHECK( gkyl_compare(gkyl_wv_ten_moment_k0(eqn), k0, 1e-15) );
  TEST_CHECK( gkyl_wv_ten_moment_use_grad_closure(eqn) == false );
  TEST_CHECK( gkyl_wv_ten_moment_use_nn_closure(eqn) == false );
  TEST_CHECK( gkyl_wv_ten_moment_poly_order(eqn) == 1 );
  TEST_CHECK( gkyl_wv_ten_moment_ann(eqn) == 0 );

  gkyl_wv_eqn_release(eqn);
}

// Gradient-based closure flag and distinct k0/poly_order via inew.
void
test_ten_moment_grad_closure()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_inew(&(struct gkyl_wv_ten_moment_inp) {
      .k0 = 2.5, .use_grad_closure = true, .use_nn_closure = false,
      .poly_order = 3, .ann = 0, .embed_geo = 0, .use_gpu = false });

  TEST_CHECK( gkyl_compare(gkyl_wv_ten_moment_k0(eqn), 2.5, 1e-15) );
  TEST_CHECK( gkyl_wv_ten_moment_use_grad_closure(eqn) == true );
  TEST_CHECK( gkyl_wv_ten_moment_use_nn_closure(eqn) == false );
  TEST_CHECK( gkyl_wv_ten_moment_poly_order(eqn) == 3 );

  gkyl_wv_eqn_release(eqn);
}

// Max speed = |u| + sqrt(3 Pxx / rho).
void
test_ten_moment_max_speed()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double rho = 2.0, u = 1.5, pxx = 4.0;
  double q[10]; calcq(rho, u, 0.3, -0.2, pxx, 0.1, 0.2, 3.0, 0.05, 2.0, q);

  double expect = fabs(u) + sqrt(3.0*pxx/rho);
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q), expect, 1e-13) );

  gkyl_wv_eqn_release(eqn);
}

// check_inv: requires positive density and positive diagonal pressures.
void
test_ten_moment_check_inv()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double qgood[10]; calcq(1.0, 0.1, 0.2, 0.3, 2.0, 0.1, 0.0, 3.0, 0.0, 1.5, qgood);
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qgood) == true );

  // Negative density.
  double qbad_rho[10]; calcq(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, qbad_rho);
  qbad_rho[0] = -1.0;
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qbad_rho) == false );

  // Negative Pyy (set primitive pyy < 0).
  double qbad_p[10]; calcq(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -2.0, 0.0, 1.0, qbad_p);
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qbad_p) == false );

  gkyl_wv_eqn_release(eqn);
}

// Default cons_to_diag copies all 10 conserved components.
void
test_ten_moment_cons_to_diag()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double q[10]; calcq(1.3, 0.5, -0.4, 0.2, 2.0, 0.1, 0.05, 1.5, 0.02, 1.1, q);
  double diag[10];
  eqn->cons_to_diag(eqn, q, diag);
  for (int i=0; i<10; ++i) TEST_CHECK( diag[i] == q[i] );

  gkyl_wv_eqn_release(eqn);
}

// Rotation round-trip recovers the full 10-component state, and density is invariant.
void
test_ten_moment_rotate_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double norm[3] = { 0.0, 1.0, 0.0 };
  double tau1[3] = { -1.0, 0.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double q[10]; calcq(1.3, 0.5, -0.7, 1.1, 2.0, 0.3, 0.1, 1.7, 0.2, 1.4, q);
  double qlocal[10], qback[10];

  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, q, qlocal);
  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, qlocal, qback);

  for (int i=0; i<10; ++i) TEST_CHECK( gkyl_compare(qback[i], q[i], 1e-13) );
  TEST_CHECK( qlocal[0] == q[0] );

  gkyl_wv_eqn_release(eqn);
}

// Riemann round-trip recovers conserved state (identity transform).
void
test_ten_moment_riem_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double qstate[10]; calcq(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, qstate);
  double qin[10]; calcq(2.0, 1.0, -1.0, 0.5, 2.0, 0.2, 0.1, 1.5, 0.05, 1.2, qin);
  double w[10], qout[10];

  eqn->cons_to_riem(eqn, qstate, qin, w);
  eqn->riem_to_cons(eqn, qstate, w, qout);
  for (int i=0; i<10; ++i) TEST_CHECK( gkyl_compare(qout[i], qin[i], 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Roe solver: equal states => zero jump => zero fluctuations on both sides.
void
test_ten_moment_waves_zero_jump()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 1, 0, false);

  double q[10]; calcq(1.4, 0.5, 0.3, -0.2, 2.0, 0.1, 0.05, 1.6, 0.02, 1.3, q);
  double delta[10] = { 0 };
  double waves[5*10], speeds[5];

  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q, q, 0.0, 0.0, waves, speeds);

  double amdq[10], apdq[10];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, q, q, 0.0, 0.0, waves, speeds, amdq, apdq);
  for (int i=0; i<10; ++i) {
    TEST_CHECK( gkyl_compare(amdq[i], 0.0, 1e-12) );
    TEST_CHECK( gkyl_compare(apdq[i], 0.0, 1e-12) );
  }

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "ten_moment_struct", test_ten_moment_struct },
  { "ten_moment_grad_closure", test_ten_moment_grad_closure },
  { "ten_moment_max_speed", test_ten_moment_max_speed },
  { "ten_moment_check_inv", test_ten_moment_check_inv },
  { "ten_moment_cons_to_diag", test_ten_moment_cons_to_diag },
  { "ten_moment_rotate_roundtrip", test_ten_moment_rotate_roundtrip },
  { "ten_moment_riem_roundtrip", test_ten_moment_riem_roundtrip },
  { "ten_moment_waves_zero_jump", test_ten_moment_waves_zero_jump },
  { NULL, NULL },
};
