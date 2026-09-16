#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_burgers.h>
#include <gkyl_wv_burgers_priv.h>

// Structural properties of the inviscid Burgers' equation object.
void
test_burgers_basic_roe()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->num_waves == 1 ); // Roe default
  TEST_CHECK( eqn->num_diag == 1 );
  TEST_CHECK( eqn->type == GKYL_EQN_BURGERS );

  gkyl_wv_eqn_release(eqn);
}

void
test_burgers_basic_lax()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_inew(&(struct gkyl_wv_burgers_inp) {
      .rp_type = WV_BURGERS_RP_LAX, .use_gpu = false });

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->num_waves == 2 );

  gkyl_wv_eqn_release(eqn);
}

// Flux F(u) = u^2/2 ; F'(u) = u.
void
test_burgers_flux()
{
  double q[1] = { 3.0 };
  double flux[1], flux_deriv[1];

  gkyl_burgers_flux(q, flux);
  gkyl_burgers_flux_deriv(q, flux_deriv);

  TEST_CHECK( gkyl_compare(flux[0], 0.5*q[0]*q[0], 1e-15) );
  TEST_CHECK( gkyl_compare(flux_deriv[0], q[0], 1e-15) );

  // Symmetry: F(u) == F(-u).
  double qm[1] = { -3.0 }, fm[1];
  gkyl_burgers_flux(qm, fm);
  TEST_CHECK( gkyl_compare(fm[0], flux[0], 1e-15) );
}

// Max speed = |u|.
void
test_burgers_max_speed()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  double q1[1] = { 2.5 }, q2[1] = { -7.0 };
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q1), 2.5, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q2), 7.0, 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Rotation is identity (scalar). Round-trip recovers state.
void
test_burgers_rotate()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  double norm[3] = { 0.0, 0.0, 1.0 };
  double tau1[3] = { 1.0, 0.0, 0.0 };
  double tau2[3] = { 0.0, 1.0, 0.0 };

  double q[1] = { 4.4 }, qlocal[1], qback[1];
  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, q, qlocal);
  TEST_CHECK( qlocal[0] == q[0] );
  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, qlocal, qback);
  TEST_CHECK( qback[0] == q[0] );

  gkyl_wv_eqn_release(eqn);
}

// Riemann round-trip recovers state.
void
test_burgers_riem_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  double qstate[1] = { 2.0 }, qin[1] = { 8.5 }, w[1], qout[1];
  eqn->cons_to_riem(eqn, qstate, qin, w);
  eqn->riem_to_cons(eqn, qstate, w, qout);
  TEST_CHECK( gkyl_compare(qout[0], qin[0], 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Flux jump = F(qr) - F(ql) = (qr^2 - ql^2)/2.
void
test_burgers_flux_jump()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  double ql[1] = { 1.0 }, qr[1] = { 3.0 }, fjump[1];
  double maxs = gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  TEST_CHECK( gkyl_compare(fjump[0], 0.5*(qr[0]*qr[0]-ql[0]*ql[0]), 1e-14) );
  TEST_CHECK( gkyl_compare(maxs, fmax(fabs(ql[0]), fabs(qr[0])), 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// check_inv always true.
void
test_burgers_check_inv()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);
  double q[1] = { -10.0 };
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, q) == true );
  gkyl_wv_eqn_release(eqn);
}

// Source term is zero.
void
test_burgers_source()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);
  double q[1] = { 5.0 }, s[1] = { 123.0 };
  gkyl_wv_eqn_source(eqn, q, s);
  TEST_CHECK( s[0] == 0.0 );
  gkyl_wv_eqn_release(eqn);
}

// Roe wave: wave = delta, speed = Roe average = (ul+ur)/2.
// Conservation: amdq + apdq = F(qr) - F(ql).
// For Burgers the Roe speed (ql+qr)/2 makes wave*s exactly the flux jump.
void
test_burgers_waves_roe()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  double ql[1] = { 1.0 }, qr[1] = { 3.0 };
  double delta[1] = { qr[0]-ql[0] };
  double waves[1], speeds[1];

  double maxs = gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double roe = 0.5*(ql[0]+qr[0]);
  TEST_CHECK( gkyl_compare(speeds[0], roe, 1e-15) );
  TEST_CHECK( gkyl_compare(maxs, roe, 1e-15) );
  TEST_CHECK( gkyl_compare(waves[0], delta[0], 1e-15) );

  // wave * speed should equal flux jump exactly (Roe consistency).
  double fjump[1];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  TEST_CHECK( gkyl_compare(waves[0]*speeds[0], fjump[0], 1e-14) );

  double amdq[1], apdq[1];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  // roe > 0 here so all right-going.
  TEST_CHECK( amdq[0] == 0.0 );
  TEST_CHECK( gkyl_compare(apdq[0]+amdq[0], fjump[0], 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// Roe with negative Roe speed -> left-going fluctuation.
void
test_burgers_waves_roe_negative()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new(false);

  double ql[1] = { -3.0 }, qr[1] = { -1.0 };
  double delta[1] = { qr[0]-ql[0] };
  double waves[1], speeds[1];
  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  TEST_CHECK( speeds[0] < 0.0 );

  double amdq[1], apdq[1];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  TEST_CHECK( apdq[0] == 0.0 );
  double fjump[1];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  TEST_CHECK( gkyl_compare(amdq[0]+apdq[0], fjump[0], 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// Lax solver: symmetric speeds +-amax, conservation holds.
void
test_burgers_waves_lax()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_inew(&(struct gkyl_wv_burgers_inp) {
      .rp_type = WV_BURGERS_RP_LAX, .use_gpu = false });

  double ql[1] = { 1.0 }, qr[1] = { 4.0 };
  double delta[1] = { qr[0]-ql[0] };
  double waves[2], speeds[2];

  gkyl_wv_eqn_waves(eqn, GKYL_WV_LOW_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double amax = fmax(fabs(ql[0]), fabs(qr[0]));
  TEST_CHECK( gkyl_compare(speeds[0], -amax, 1e-15) );
  TEST_CHECK( gkyl_compare(speeds[1],  amax, 1e-15) );

  double amdq[1], apdq[1];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_LOW_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  double fjump[1];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  TEST_CHECK( gkyl_compare(amdq[0]+apdq[0], fjump[0], 1e-13) );

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "burgers_basic_roe", test_burgers_basic_roe },
  { "burgers_basic_lax", test_burgers_basic_lax },
  { "burgers_flux", test_burgers_flux },
  { "burgers_max_speed", test_burgers_max_speed },
  { "burgers_rotate", test_burgers_rotate },
  { "burgers_riem_roundtrip", test_burgers_riem_roundtrip },
  { "burgers_flux_jump", test_burgers_flux_jump },
  { "burgers_check_inv", test_burgers_check_inv },
  { "burgers_source", test_burgers_source },
  { "burgers_waves_roe", test_burgers_waves_roe },
  { "burgers_waves_roe_negative", test_burgers_waves_roe_negative },
  { "burgers_waves_lax", test_burgers_waves_lax },
  { NULL, NULL },
};
