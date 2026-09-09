#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_maxwell_priv.h>

// State layout: q[0..7] = {Ex, Ey, Ez, Bx, By, Bz, phi, psi}.

void
test_maxwell_extra_basic()
{
  double c = 1.0, e_fact = 1.0, b_fact = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, e_fact, b_fact, false);

  TEST_CHECK( eqn->num_equations == 8 );
  TEST_CHECK( eqn->num_waves == 6 ); // Roe default
  TEST_CHECK( eqn->num_diag == 6 );
  TEST_CHECK( eqn->type == GKYL_EQN_MAXWELL );

  gkyl_wv_eqn_release(eqn);
}

// Explicit flux for the perfectly hyperbolic Maxwell system.
void
test_maxwell_extra_flux()
{
  double c = 2.0, e_fact = 1.5, b_fact = 1.2;
  double q[8] = { 0.3, -0.4, 0.5, 0.6, -0.7, 0.8, 0.9, -1.0 };
  double flux[8];
  gkyl_maxwell_flux(c, e_fact, b_fact, q, flux);

  TEST_CHECK( gkyl_compare(flux[0], e_fact*c*c*q[6], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[1], c*c*q[5], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[2], -c*c*q[4], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[3], b_fact*q[7], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[4], -q[2], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[5], q[1], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[6], e_fact*q[0], 1e-14) );
  TEST_CHECK( gkyl_compare(flux[7], b_fact*c*c*q[3], 1e-14) );
}

// Max speed: with unit correction factors it is just the speed of light.
void
test_maxwell_extra_max_speed()
{
  double c = 3.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, 1.0, 1.0, false);
  double q[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q), c, 1e-14) );
  gkyl_wv_eqn_release(eqn);

  // With correction factor > 1, max speed scales by the larger factor.
  double e_fact = 2.0, b_fact = 0.5;
  struct gkyl_wv_eqn *eqn2 = gkyl_wv_maxwell_new(c, e_fact, b_fact, false);
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn2, q), c*e_fact, 1e-14) );
  gkyl_wv_eqn_release(eqn2);
}

// Diagnostics are the squared components of the first 6 (field) variables.
void
test_maxwell_extra_cons_to_diag()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(1.0, 1.0, 1.0, false);
  double q[8] = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, 8.0 };
  double diag[6];
  eqn->cons_to_diag(eqn, q, diag);
  for (int i=0; i<6; ++i) TEST_CHECK( gkyl_compare(diag[i], q[i]*q[i], 1e-14) );
  gkyl_wv_eqn_release(eqn);
}

// Rotation round-trip across a non-axis-aligned orthonormal frame recovers q.
void
test_maxwell_extra_rotate_roundtrip()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(1.0, 1.0, 1.0, false);

  // Build an orthonormal frame: norm = tau1 x tau2.
  double inv = 1.0/sqrt(2.0);
  double norm[3] = { inv, inv, 0.0 };
  double tau1[3] = { -inv, inv, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double q[8] = { 0.3, -0.4, 0.5, 0.6, -0.7, 0.8, 0.9, -1.0 };
  double qlocal[8], qback[8];

  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, q, qlocal);
  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, qlocal, qback);

  for (int i=0; i<8; ++i) TEST_CHECK( gkyl_compare(qback[i], q[i], 1e-13) );

  // Scalar potentials unchanged by rotation.
  TEST_CHECK( gkyl_compare(qlocal[6], q[6], 1e-14) );
  TEST_CHECK( gkyl_compare(qlocal[7], q[7], 1e-14) );

  // Rotation preserves the magnitude of the E and B field vectors.
  double E2g = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
  double E2l = qlocal[0]*qlocal[0]+qlocal[1]*qlocal[1]+qlocal[2]*qlocal[2];
  TEST_CHECK( gkyl_compare(E2g, E2l, 1e-13) );
  double B2g = q[3]*q[3]+q[4]*q[4]+q[5]*q[5];
  double B2l = qlocal[3]*qlocal[3]+qlocal[4]*qlocal[4]+qlocal[5]*qlocal[5];
  TEST_CHECK( gkyl_compare(B2g, B2l, 1e-13) );

  gkyl_wv_eqn_release(eqn);
}

// Flux jump equals F(qr) - F(ql) computed directly from the flux function.
void
test_maxwell_extra_flux_jump()
{
  double c = 1.5, e_fact = 1.0, b_fact = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, e_fact, b_fact, false);

  double ql[8] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  double qr[8] = { 1.1, -0.9, 0.2, -0.3, 0.8, -0.1, 0.4, 0.5 };

  double fjump[8];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  double fl[8], fr[8];
  gkyl_maxwell_flux(c, e_fact, b_fact, ql, fl);
  gkyl_maxwell_flux(c, e_fact, b_fact, qr, fr);

  for (int i=0; i<8; ++i)
    TEST_CHECK( gkyl_compare(fjump[i], fr[i]-fl[i], 1e-13) );

  gkyl_wv_eqn_release(eqn);
}

// Roe solver q-fluctuation conservation: amdq + apdq = F(qr) - F(ql).
void
test_maxwell_extra_waves_conservation()
{
  double c = 1.0, e_fact = 1.0, b_fact = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, e_fact, b_fact, false);

  double ql[8] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  double qr[8] = { 1.1, -0.9, 0.2, -0.3, 0.8, -0.1, 0.4, 0.5 };
  double delta[8];
  for (int i=0; i<8; ++i) delta[i] = qr[i]-ql[i];

  double waves[6*8], speeds[6];
  double maxs = gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);
  TEST_CHECK( maxs > 0.0 );

  double amdq[8], apdq[8];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  double fjump[8];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  for (int i=0; i<8; ++i)
    TEST_CHECK( gkyl_compare(amdq[i]+apdq[i], fjump[i], 1e-12) );

  gkyl_wv_eqn_release(eqn);
}

// The sum of all Roe waves must reconstruct the full jump delta.
void
test_maxwell_extra_wave_sum()
{
  double c = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, 1.0, 1.0, false);

  double ql[8] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  double qr[8] = { 1.1, -0.9, 0.2, -0.3, 0.8, -0.1, 0.4, 0.5 };
  double delta[8];
  for (int i=0; i<8; ++i) delta[i] = qr[i]-ql[i];

  double waves[6*8], speeds[6];
  gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  for (int i=0; i<8; ++i) {
    double sum = 0.0;
    for (int w=0; w<6; ++w) sum += waves[w*8 + i];
    TEST_CHECK( gkyl_compare(sum, delta[i], 1e-12) );
  }

  gkyl_wv_eqn_release(eqn);
}

// Lax solver: structural check and conservation of fluctuations.
void
test_maxwell_extra_waves_lax()
{
  double c = 1.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_inew(&(struct gkyl_wv_maxwell_inp) {
      .c = c, .e_fact = 1.0, .b_fact = 1.0,
      .rp_type = WV_MAXWELL_RP_LAX, .use_gpu = false });

  TEST_CHECK( eqn->num_waves == 2 );

  double ql[8] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  double qr[8] = { 1.1, -0.9, 0.2, -0.3, 0.8, -0.1, 0.4, 0.5 };
  double delta[8];
  for (int i=0; i<8; ++i) delta[i] = qr[i]-ql[i];

  double waves[2*8], speeds[2];
  gkyl_wv_eqn_waves(eqn, GKYL_WV_LOW_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double amdq[8], apdq[8];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_LOW_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  double fjump[8];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  for (int i=0; i<8; ++i)
    TEST_CHECK( gkyl_compare(amdq[i]+apdq[i], fjump[i], 1e-11) );

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "maxwell_extra_basic", test_maxwell_extra_basic },
  { "maxwell_extra_flux", test_maxwell_extra_flux },
  { "maxwell_extra_max_speed", test_maxwell_extra_max_speed },
  { "maxwell_extra_cons_to_diag", test_maxwell_extra_cons_to_diag },
  { "maxwell_extra_rotate_roundtrip", test_maxwell_extra_rotate_roundtrip },
  { "maxwell_extra_flux_jump", test_maxwell_extra_flux_jump },
  { "maxwell_extra_waves_conservation", test_maxwell_extra_waves_conservation },
  { "maxwell_extra_wave_sum", test_maxwell_extra_wave_sum },
  { "maxwell_extra_waves_lax", test_maxwell_extra_waves_lax },
  { NULL, NULL },
};
