#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_coldfluid.h>

// Build conserved state {rho, rho*u, rho*v, rho*w} from primitives.
static void
calcq(double rho, double u, double v, double w, double q[4])
{
  q[0] = rho; q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
}

// Structural properties of the cold-fluid equation object.
void
test_coldfluid_basic()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  TEST_CHECK( eqn->num_equations == 4 );
  TEST_CHECK( eqn->num_waves == 2 );
  TEST_CHECK( eqn->num_diag == 5 ); // KE is final diagnostic component
  TEST_CHECK( eqn->type == GKYL_EQN_COLDFLUID );

  gkyl_wv_eqn_release(eqn);
}

// Max speed = |u| = |rho*u / rho|.
void
test_coldfluid_max_speed()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double q[4];
  calcq(2.0, -3.0, 1.0, 5.0, q);
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q), 3.0, 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// check_inv: valid iff density positive.
void
test_coldfluid_check_inv()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double qgood[4]; calcq(1.0, 0.1, 0.2, 0.3, qgood);
  double qbad[4]  = { -1.0, 0.0, 0.0, 0.0 };

  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qgood) == true );
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qbad) == false );

  gkyl_wv_eqn_release(eqn);
}

// Source term is zero (homogeneous).
void
test_coldfluid_source()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();
  double q[4]; calcq(1.0, 1.0, 1.0, 1.0, q);
  double s[4] = { 9, 9, 9, 9 };
  gkyl_wv_eqn_source(eqn, q, s);
  for (int i=0; i<4; ++i) TEST_CHECK( s[i] == 0.0 );
  gkyl_wv_eqn_release(eqn);
}

// Diagnostics: first 4 are conserved vars, 5th is kinetic energy density
// KE = 0.5*(|rho*v|^2)/rho.
void
test_coldfluid_cons_to_diag()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double rho = 2.0, u = 1.0, v = 2.0, w = 3.0;
  double q[4]; calcq(rho, u, v, w, q);
  double diag[5];
  eqn->cons_to_diag(eqn, q, diag);

  for (int i=0; i<4; ++i) TEST_CHECK( diag[i] == q[i] );
  double ke = 0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0];
  TEST_CHECK( gkyl_compare(diag[4], ke, 1e-14) );
  // KE = 0.5*rho*(u^2+v^2+w^2).
  TEST_CHECK( gkyl_compare(diag[4], 0.5*rho*(u*u+v*v+w*w), 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// Rotation round-trip recovers the global state for the 3-vector momentum.
void
test_coldfluid_rotate_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  // Orthonormal frame in the y direction.
  double norm[3] = { 0.0, 1.0, 0.0 };
  double tau1[3] = { -1.0, 0.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double q[4]; calcq(1.3, 0.5, -0.7, 1.1, q);
  double qlocal[4], qback[4];

  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, q, qlocal);
  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, qlocal, qback);

  for (int i=0; i<4; ++i) TEST_CHECK( gkyl_compare(qback[i], q[i], 1e-14) );

  // Density is rotation-invariant.
  TEST_CHECK( qlocal[0] == q[0] );
  // Normal momentum component in local frame is q . norm = rho*v.
  TEST_CHECK( gkyl_compare(qlocal[1], q[2], 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// Riemann round-trip recovers conserved state (identity transform here).
void
test_coldfluid_riem_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double qstate[4]; calcq(1.0, 0.0, 0.0, 0.0, qstate);
  double qin[4]; calcq(2.0, 1.0, -1.0, 0.5, qin);
  double w[4], qout[4];

  eqn->cons_to_riem(eqn, qstate, qin, w);
  eqn->riem_to_cons(eqn, qstate, w, qout);

  for (int i=0; i<4; ++i) TEST_CHECK( gkyl_compare(qout[i], qin[i], 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Flux jump in x: F = {rho*u, rho*u*u, rho*v*u, rho*w*u}.
void
test_coldfluid_flux_jump()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double ul = 2.0, vl = 1.0, wl = -1.0, rl = 1.5;
  double ur = 3.0, vr = -2.0, wr = 0.5, rr = 2.0;
  double ql[4]; calcq(rl, ul, vl, wl, ql);
  double qr[4]; calcq(rr, ur, vr, wr, qr);

  double fjump[4];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  double fl[4] = { rl*ul, rl*ul*ul, rl*vl*ul, rl*wl*ul };
  double fr[4] = { rr*ur, rr*ur*ur, rr*vr*ur, rr*wr*ur };
  for (int i=0; i<4; ++i)
    TEST_CHECK( gkyl_compare(fjump[i], fr[i]-fl[i], 1e-13) );

  gkyl_wv_eqn_release(eqn);
}

// Roe solver: for identical states the jump is zero, so all fluctuations vanish
// and the wave speed equals the common flow speed u.
void
test_coldfluid_waves_zero_jump()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double q[4]; calcq(1.7, 2.0, 0.3, -0.4, q);
  double delta[4] = { 0.0, 0.0, 0.0, 0.0 };
  double waves[8], speeds[2];

  double maxs = gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q, q, 0.0, 0.0, waves, speeds);

  // Both Roe speeds equal u = 2.0; max speed is |u|.
  TEST_CHECK( gkyl_compare(speeds[0], 2.0, 1e-14) );
  TEST_CHECK( gkyl_compare(speeds[1], 2.0, 1e-14) );
  TEST_CHECK( gkyl_compare(maxs, 2.0, 1e-14) );

  double amdq[4], apdq[4];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, q, q, 0.0, 0.0, waves, speeds, amdq, apdq);
  for (int i=0; i<4; ++i) {
    TEST_CHECK( gkyl_compare(amdq[i], 0.0, 1e-14) );
    TEST_CHECK( gkyl_compare(apdq[i], 0.0, 1e-14) );
  }

  gkyl_wv_eqn_release(eqn);
}

// Roe solver, both speeds positive: all fluctuation right-going and equal to the
// flux jump (conservation). Use states with the same positive velocity so the
// Roe average is exactly that velocity and waves carry the full delta.
void
test_coldfluid_waves_conservation()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  // Same positive velocity u=2 on both sides; differing density/transverse mom.
  double ql[4]; calcq(1.0, 2.0, 0.0, 0.0, ql);
  double qr[4]; calcq(3.0, 2.0, 1.0, -1.0, qr);
  double delta[4] = { qr[0]-ql[0], qr[1]-ql[1], qr[2]-ql[2], qr[3]-ql[3] };
  double waves[8], speeds[2];

  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  TEST_CHECK( speeds[0] > 0.0 );
  TEST_CHECK( speeds[1] > 0.0 );

  double amdq[4], apdq[4];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  double fjump[4];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  for (int i=0; i<4; ++i) {
    TEST_CHECK( gkyl_compare(amdq[i], 0.0, 1e-13) );      // all right-going
    TEST_CHECK( gkyl_compare(amdq[i]+apdq[i], fjump[i], 1e-12) );
  }

  gkyl_wv_eqn_release(eqn);
}

// f-fluctuations: amdq + apdq must equal the total f-wave content (sum of waves).
void
test_coldfluid_ffluct_conservation()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_coldfluid_new();

  double ql[4]; calcq(1.0, 2.0, 0.0, 0.0, ql);
  double qr[4]; calcq(3.0, 2.0, 1.0, -1.0, qr);
  double delta[4] = { qr[0]-ql[0], qr[1]-ql[1], qr[2]-ql[2], qr[3]-ql[3] };
  double waves[8], speeds[2];

  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double amdq[4], apdq[4];
  gkyl_wv_eqn_ffluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  for (int i=0; i<4; ++i) {
    double wsum = waves[i] + waves[4+i];
    TEST_CHECK( gkyl_compare(amdq[i]+apdq[i], wsum, 1e-13) );
  }

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "coldfluid_basic", test_coldfluid_basic },
  { "coldfluid_max_speed", test_coldfluid_max_speed },
  { "coldfluid_check_inv", test_coldfluid_check_inv },
  { "coldfluid_source", test_coldfluid_source },
  { "coldfluid_cons_to_diag", test_coldfluid_cons_to_diag },
  { "coldfluid_rotate_roundtrip", test_coldfluid_rotate_roundtrip },
  { "coldfluid_riem_roundtrip", test_coldfluid_riem_roundtrip },
  { "coldfluid_flux_jump", test_coldfluid_flux_jump },
  { "coldfluid_waves_zero_jump", test_coldfluid_waves_zero_jump },
  { "coldfluid_waves_conservation", test_coldfluid_waves_conservation },
  { "coldfluid_ffluct_conservation", test_coldfluid_ffluct_conservation },
  { NULL, NULL },
};
