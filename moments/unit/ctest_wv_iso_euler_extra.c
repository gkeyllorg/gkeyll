#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_wv_iso_euler.h>

// Build conserved state {rho, rho*u, rho*v, rho*w}.
static void
calcq(double rho, double u, double v, double w, double q[4])
{
  q[0] = rho; q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
}

// Constructor wiring, equation/diag counts, type tag, and vt accessor.
void
test_iso_euler_struct()
{
  double vt = 3.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(vt, false);

  TEST_CHECK( eqn->type == GKYL_EQN_ISO_EULER );
  TEST_CHECK( eqn->num_equations == 4 );
  TEST_CHECK( eqn->num_diag == 4 );
  // Default RP type is Roe => 3 waves.
  TEST_CHECK( eqn->num_waves == 3 );
  TEST_CHECK( gkyl_compare(gkyl_wv_iso_euler_vt(eqn), vt, 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Lax RP selects 2 waves; Roe RP selects 3. inew honours rp_type.
void
test_iso_euler_rp_types()
{
  struct gkyl_wv_eqn *lax = gkyl_wv_iso_euler_inew(&(struct gkyl_wv_iso_euler_inp) {
      .vt = 1.0, .rp_type = WV_ISO_EULER_RP_LAX, .use_gpu = false });
  struct gkyl_wv_eqn *roe = gkyl_wv_iso_euler_inew(&(struct gkyl_wv_iso_euler_inp) {
      .vt = 1.0, .rp_type = WV_ISO_EULER_RP_ROE, .use_gpu = false });

  TEST_CHECK( lax->num_waves == 2 );
  TEST_CHECK( roe->num_waves == 3 );
  TEST_CHECK( gkyl_compare(gkyl_wv_iso_euler_vt(lax), 1.0, 1e-15) );

  gkyl_wv_eqn_release(lax);
  gkyl_wv_eqn_release(roe);
}

// Max speed = |u| + vt.
void
test_iso_euler_max_speed()
{
  double vt = 2.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(vt, false);

  double q[4]; calcq(2.0, -3.0, 1.0, 5.0, q);
  // fmax(|u-vt|, |u+vt|) = |u| + vt = 3 + 2 = 5.
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q), 5.0, 1e-14) );

  double q2[4]; calcq(1.0, 0.5, 0.0, 0.0, q2);
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q2), 2.5, 1e-14) );

  gkyl_wv_eqn_release(eqn);
}

// check_inv: valid iff density positive.
void
test_iso_euler_check_inv()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(1.0, false);

  double qgood[4]; calcq(1.0, 0.1, 0.2, 0.3, qgood);
  double qbad[4]  = { -0.5, 0.0, 0.0, 0.0 };
  double qzero[4] = { 0.0, 0.0, 0.0, 0.0 };

  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qgood) == true );
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qbad) == false );
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qzero) == false );

  gkyl_wv_eqn_release(eqn);
}

// Diagnostics are just the conserved variables (no extra KE component).
void
test_iso_euler_cons_to_diag()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(1.0, false);

  double q[4]; calcq(2.0, 1.0, -2.0, 3.0, q);
  double diag[4];
  eqn->cons_to_diag(eqn, q, diag);

  for (int i=0; i<4; ++i) TEST_CHECK( diag[i] == q[i] );

  gkyl_wv_eqn_release(eqn);
}

// Source term is homogeneous (zero).
void
test_iso_euler_source()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(1.0, false);

  double q[4]; calcq(1.5, 1.0, 1.0, 1.0, q);
  double s[4] = { 7, 7, 7, 7 };
  gkyl_wv_eqn_source(eqn, q, s);
  for (int i=0; i<4; ++i) TEST_CHECK( s[i] == 0.0 );

  gkyl_wv_eqn_release(eqn);
}

// Flux jump equals the analytic flux difference in the x direction.
// F = {rho*u, rho*u^2 + rho*vt^2, rho*u*v, rho*u*w}.
void
test_iso_euler_flux_jump()
{
  double vt = 1.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(vt, false);

  double rl = 1.5, ul = 2.0, vl = 1.0, wl = -1.0;
  double rr = 2.0, ur = 3.0, vr = -2.0, wr = 0.5;
  double ql[4]; calcq(rl, ul, vl, wl, ql);
  double qr[4]; calcq(rr, ur, vr, wr, qr);

  double fjump[4];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  double fl[4] = { rl*ul, rl*ul*ul + rl*vt*vt, rl*ul*vl, rl*ul*wl };
  double fr[4] = { rr*ur, rr*ur*ur + rr*vt*vt, rr*ur*vr, rr*ur*wr };
  for (int i=0; i<4; ++i)
    TEST_CHECK( gkyl_compare(fjump[i], fr[i]-fl[i], 1e-12) );

  gkyl_wv_eqn_release(eqn);
}

// Riemann round-trip recovers conserved state (identity transform here).
void
test_iso_euler_riem_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(1.0, false);

  double qstate[4]; calcq(1.0, 0.0, 0.0, 0.0, qstate);
  double qin[4]; calcq(2.0, 1.0, -1.0, 0.5, qin);
  double w[4], qout[4];

  eqn->cons_to_riem(eqn, qstate, qin, w);
  eqn->riem_to_cons(eqn, qstate, w, qout);

  for (int i=0; i<4; ++i) TEST_CHECK( gkyl_compare(qout[i], qin[i], 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Roe solver, equal states => zero jump => zero fluctuations.
void
test_iso_euler_waves_zero_jump()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_new(1.0, false);

  double q[4]; calcq(1.7, 0.6, 0.3, -0.4, q);
  double delta[4] = { 0.0, 0.0, 0.0, 0.0 };
  double waves[3*4], speeds[3];

  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q, q, 0.0, 0.0, waves, speeds);

  double amdq[4], apdq[4];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, q, q, 0.0, 0.0, waves, speeds, amdq, apdq);
  for (int i=0; i<4; ++i) {
    TEST_CHECK( gkyl_compare(amdq[i], 0.0, 1e-13) );
    TEST_CHECK( gkyl_compare(apdq[i], 0.0, 1e-13) );
  }

  gkyl_wv_eqn_release(eqn);
}

// Lax solver: fluctuations are conservative (amdq + apdq == flux jump) in x.
void
test_iso_euler_lax_conservation()
{
  double vt = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_inew(&(struct gkyl_wv_iso_euler_inp) {
      .vt = vt, .rp_type = WV_ISO_EULER_RP_LAX, .use_gpu = false });

  double ql[4]; calcq(1.0, 0.2, 0.1, -0.1, ql);
  double qr[4]; calcq(2.0, 0.3, -0.2, 0.4, qr);
  double delta[4];
  for (int i=0; i<4; ++i) delta[i] = qr[i]-ql[i];

  double waves[2*4], speeds[2];
  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double amdq[4], apdq[4];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  double fjump[4];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  for (int i=0; i<4; ++i)
    TEST_CHECK( gkyl_compare(amdq[i]+apdq[i], fjump[i], 1e-12) );

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "iso_euler_struct", test_iso_euler_struct },
  { "iso_euler_rp_types", test_iso_euler_rp_types },
  { "iso_euler_max_speed", test_iso_euler_max_speed },
  { "iso_euler_check_inv", test_iso_euler_check_inv },
  { "iso_euler_cons_to_diag", test_iso_euler_cons_to_diag },
  { "iso_euler_source", test_iso_euler_source },
  { "iso_euler_flux_jump", test_iso_euler_flux_jump },
  { "iso_euler_riem_roundtrip", test_iso_euler_riem_roundtrip },
  { "iso_euler_waves_zero_jump", test_iso_euler_waves_zero_jump },
  { "iso_euler_lax_conservation", test_iso_euler_lax_conservation },
  { NULL, NULL },
};
