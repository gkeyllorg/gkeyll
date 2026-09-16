#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_advect.h>
#include <gkyl_wv_advect_priv.h>

// Basic structural properties of the linear advection equation object.
void
test_advect_basic_roe()
{
  double a = 2.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(a, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->num_waves == 1 ); // Roe is default
  TEST_CHECK( eqn->num_diag == 1 );
  TEST_CHECK( eqn->type == GKYL_EQN_ADVECTION );

  gkyl_wv_eqn_release(eqn);
}

void
test_advect_basic_lax()
{
  double a = -1.3;
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_inew(&(struct gkyl_wv_advect_inp) {
      .a = a,
      .rp_type = WV_ADVECT_RP_LAX,
      .use_gpu = false,
    });

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->num_waves == 2 ); // Lax has 2 waves
  TEST_CHECK( eqn->num_diag == 1 );

  gkyl_wv_eqn_release(eqn);
}

// Flux F(q) = a*q, and flux derivative df/dq = a (the advection speed).
void
test_advect_flux()
{
  double a = 3.7;
  double q[1] = { 4.2 };
  double flux[1], flux_deriv[1];

  gkyl_advect_flux(a, q, flux);
  gkyl_advect_flux_deriv(a, q, flux_deriv);

  TEST_CHECK( gkyl_compare(flux[0], a*q[0], 1e-15) );
  TEST_CHECK( gkyl_compare(flux_deriv[0], a, 1e-15) );

  // Flux is linear: F(2q) = 2 F(q).
  double q2[1] = { 2.0*q[0] }, flux2[1];
  gkyl_advect_flux(a, q2, flux2);
  TEST_CHECK( gkyl_compare(flux2[0], 2.0*flux[0], 1e-14) );
}

// Max speed is |a| independent of state.
void
test_advect_max_speed()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(-4.0, false);

  double q1[1] = { 1.0 }, q2[1] = { -100.0 };
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q1), 4.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q2), 4.0, 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Rotation to local/global is the identity (scalar eqn). Round-trip recovers q.
void
test_advect_rotate()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(1.0, false);

  double norm[3] = { 0.0, 1.0, 0.0 };
  double tau1[3] = { 1.0, 0.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double q[1] = { 7.3 }, qlocal[1], qback[1];
  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, q, qlocal);
  TEST_CHECK( qlocal[0] == q[0] );
  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, qlocal, qback);
  TEST_CHECK( qback[0] == q[0] );

  gkyl_wv_eqn_release(eqn);
}

// Riemann variable round-trip (identity here) recovers state.
void
test_advect_riem_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(2.0, false);

  double qstate[1] = { 1.5 }, qin[1] = { 9.1 }, w[1], qout[1];
  eqn->cons_to_riem(eqn, qstate, qin, w);
  eqn->riem_to_cons(eqn, qstate, w, qout);
  TEST_CHECK( gkyl_compare(qout[0], qin[0], 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Flux jump must equal F(qr) - F(ql) = a*(qr - ql).
void
test_advect_flux_jump()
{
  double a = 1.7;
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(a, false);

  double ql[1] = { 2.0 }, qr[1] = { 5.0 }, fjump[1];
  double maxs = gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  TEST_CHECK( gkyl_compare(fjump[0], a*(qr[0]-ql[0]), 1e-14) );
  TEST_CHECK( gkyl_compare(maxs, fabs(a), 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// check_inv always true for linear advection.
void
test_advect_check_inv()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(1.0, false);
  double q[1] = { -3.0 };
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, q) == true );
  gkyl_wv_eqn_release(eqn);
}

// Source term is zero (homogeneous equation).
void
test_advect_source()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(1.0, false);
  double q[1] = { 3.0 }, s[1] = { 999.0 };
  gkyl_wv_eqn_source(eqn, q, s);
  TEST_CHECK( s[0] == 0.0 );
  gkyl_wv_eqn_release(eqn);
}

// cons_to_diag copies the state.
void
test_advect_cons_to_diag()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(1.0, false);
  double q[1] = { 6.6 }, diag[1];
  eqn->cons_to_diag(eqn, q, diag);
  TEST_CHECK( diag[0] == q[0] );
  gkyl_wv_eqn_release(eqn);
}

// Roe waves/qfluct: wave = delta, speed = a. Fluctuations split by sign of a.
// Consistency: amdq + apdq = flux jump (conservation of f-waves/q-waves through
// the Roe linearization, since wave*s = a*delta = F(qr)-F(ql)).
void
test_advect_waves_roe()
{
  double a = 2.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(a, false);

  double ql[1] = { 1.0 }, qr[1] = { 4.0 };
  double delta[1] = { qr[0]-ql[0] };
  double waves[1], speeds[1];

  double maxs = gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  TEST_CHECK( gkyl_compare(speeds[0], a, 1e-15) );
  TEST_CHECK( gkyl_compare(maxs, a, 1e-15) );
  TEST_CHECK( gkyl_compare(waves[0], delta[0], 1e-15) );

  double amdq[1], apdq[1];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  // a > 0 so all fluctuation is right-going.
  TEST_CHECK( amdq[0] == 0.0 );
  TEST_CHECK( gkyl_compare(apdq[0], a*delta[0], 1e-14) );

  // Conservation: amdq + apdq = F(qr) - F(ql).
  double fjump[1];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  TEST_CHECK( gkyl_compare(amdq[0]+apdq[0], fjump[0], 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// Roe with negative advection speed: all fluctuation should be left-going.
void
test_advect_waves_roe_negative()
{
  double a = -3.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_new(a, false);

  double ql[1] = { 5.0 }, qr[1] = { 2.0 };
  double delta[1] = { qr[0]-ql[0] };
  double waves[1], speeds[1];
  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double amdq[1], apdq[1];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  TEST_CHECK( apdq[0] == 0.0 );
  TEST_CHECK( gkyl_compare(amdq[0], a*delta[0], 1e-14) );

  double fjump[1];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  TEST_CHECK( gkyl_compare(amdq[0]+apdq[0], fjump[0], 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// Lax solver: two symmetric waves with speeds +-amax. Conservation must hold.
void
test_advect_waves_lax()
{
  double a = 2.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_advect_inew(&(struct gkyl_wv_advect_inp) {
      .a = a, .rp_type = WV_ADVECT_RP_LAX, .use_gpu = false });

  double ql[1] = { 1.0 }, qr[1] = { 4.0 };
  double delta[1] = { qr[0]-ql[0] };
  double waves[2], speeds[2];

  gkyl_wv_eqn_waves(eqn, GKYL_WV_LOW_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  TEST_CHECK( gkyl_compare(speeds[0], -fabs(a), 1e-15) );
  TEST_CHECK( gkyl_compare(speeds[1],  fabs(a), 1e-15) );

  double amdq[1], apdq[1];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_LOW_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  // Conservation: amdq + apdq = F(qr) - F(ql) = a*delta.
  TEST_CHECK( gkyl_compare(amdq[0]+apdq[0], a*delta[0], 1e-13) );

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "advect_basic_roe", test_advect_basic_roe },
  { "advect_basic_lax", test_advect_basic_lax },
  { "advect_flux", test_advect_flux },
  { "advect_max_speed", test_advect_max_speed },
  { "advect_rotate", test_advect_rotate },
  { "advect_riem_roundtrip", test_advect_riem_roundtrip },
  { "advect_flux_jump", test_advect_flux_jump },
  { "advect_check_inv", test_advect_check_inv },
  { "advect_source", test_advect_source },
  { "advect_cons_to_diag", test_advect_cons_to_diag },
  { "advect_waves_roe", test_advect_waves_roe },
  { "advect_waves_roe_negative", test_advect_waves_roe_negative },
  { "advect_waves_lax", test_advect_waves_lax },
  { NULL, NULL },
};
