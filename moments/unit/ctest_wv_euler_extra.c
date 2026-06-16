#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_priv.h>

// Conserved state q = {rho, rho*u, rho*v, rho*w, E} from primitives, where
// E = p/(gamma-1) + 0.5*rho*|v|^2.
static void
calcq(double gas_gamma, double rho, double u, double v, double w, double pr, double q[5])
{
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pr/(gas_gamma-1.0) + 0.5*rho*(u*u+v*v+w*w);
}

void
test_euler_extra_basic()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);

  TEST_CHECK( eqn->num_equations == 5 );
  TEST_CHECK( eqn->type == GKYL_EQN_EULER );
  TEST_CHECK( gkyl_compare(gkyl_wv_euler_gas_gamma(eqn), gas_gamma, 1e-15) );

  gkyl_wv_eqn_release(eqn);
}

// Pressure recovered from conserved variables matches the input pressure.
void
test_euler_extra_pressure()
{
  double gas_gamma = 1.4;
  double rho = 2.0, u = 1.0, v = -0.5, w = 0.3, pr = 3.0;
  double q[5];
  calcq(gas_gamma, rho, u, v, w, pr, q);

  TEST_CHECK( gkyl_compare(gkyl_euler_pressure(gas_gamma, q), pr, 1e-13) );
}

// Primitive-variable extraction is the exact inverse of calcq.
void
test_euler_extra_prim_vars()
{
  double gas_gamma = 5.0/3.0;
  double rho = 1.3, u = 0.7, v = 0.2, w = -0.9, pr = 2.5;
  double q[5], prim[5];
  calcq(gas_gamma, rho, u, v, w, pr, q);

  gkyl_euler_prim_vars(gas_gamma, q, prim);

  TEST_CHECK( gkyl_compare(prim[0], rho, 1e-13) );
  TEST_CHECK( gkyl_compare(prim[1], u, 1e-13) );
  TEST_CHECK( gkyl_compare(prim[2], v, 1e-13) );
  TEST_CHECK( gkyl_compare(prim[3], w, 1e-13) );
  TEST_CHECK( gkyl_compare(prim[4], pr, 1e-13) );
}

// Max abs speed = |velocity| + sound speed.
void
test_euler_extra_max_speed()
{
  double gas_gamma = 1.4;
  double rho = 1.0, u = 2.0, v = 0.0, w = 0.0, pr = 1.0;
  double q[5];
  calcq(gas_gamma, rho, u, v, w, pr, q);

  double cs = sqrt(gas_gamma*pr/rho);
  double expected = fabs(u) + cs;

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);
  TEST_CHECK( gkyl_compare(gkyl_wv_eqn_max_speed(eqn, q), expected, 1e-13) );
  TEST_CHECK( gkyl_compare(gkyl_euler_max_abs_speed(gas_gamma, q), expected, 1e-13) );
  gkyl_wv_eqn_release(eqn);
}

// Euler flux in the local (x) frame.
void
test_euler_extra_flux()
{
  double gas_gamma = 1.4;
  double rho = 1.5, u = 0.8, v = -0.3, w = 0.4, pr = 2.0;
  double q[5];
  calcq(gas_gamma, rho, u, v, w, pr, q);

  double flux[5];
  gkyl_euler_flux(gas_gamma, q, flux);

  double E = q[4];
  TEST_CHECK( gkyl_compare(flux[0], rho*u, 1e-13) );
  TEST_CHECK( gkyl_compare(flux[1], rho*u*u + pr, 1e-13) );
  TEST_CHECK( gkyl_compare(flux[2], rho*v*u, 1e-13) );
  TEST_CHECK( gkyl_compare(flux[3], rho*w*u, 1e-13) );
  TEST_CHECK( gkyl_compare(flux[4], (E+pr)*u, 1e-13) );
}

// Rotating the state to each axis frame, computing the local flux, and rotating
// back must reproduce the directional flux. The directional flux differs from
// the x-flux only by which momentum component carries the pressure.
void
test_euler_extra_flux_rotation()
{
  double gas_gamma = 1.4;
  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double q[5];
  calcq(gas_gamma, rho, u, v, w, pr, q);

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);

  double norm[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
  double tau1[3][3] = { {0,1,0}, {1,0,0}, {1,0,0} };
  double tau2[3][3] = { {0,0,1}, {0,0,-1}, {0,1,0} };

  // Expected directional fluxes (momentum eqns 1,2,3 map to x,y,z mom).
  double E = q[4];
  double fx[5] = { rho*u, rho*u*u+pr, rho*v*u, rho*w*u, (E+pr)*u };
  double fy[5] = { rho*v, rho*u*v, rho*v*v+pr, rho*w*v, (E+pr)*v };
  double fz[5] = { rho*w, rho*u*w, rho*v*w, rho*w*w+pr, (E+pr)*w };
  double *fexp[3] = { fx, fy, fz };

  for (int d=0; d<3; ++d) {
    double qloc[5], floc[5], fglob[5];
    eqn->rotate_to_local_func(eqn, tau1[d], tau2[d], norm[d], q, qloc);
    gkyl_euler_flux(gas_gamma, qloc, floc);
    eqn->rotate_to_global_func(eqn, tau1[d], tau2[d], norm[d], floc, fglob);
    for (int m=0; m<5; ++m)
      TEST_CHECK( gkyl_compare(fglob[m], fexp[d][m], 1e-13) );
  }

  gkyl_wv_eqn_release(eqn);
}

// Rotation round-trip recovers the state, and density/energy are invariant.
void
test_euler_extra_rotate_roundtrip()
{
  double gas_gamma = 1.4;
  double q[5];
  calcq(gas_gamma, 1.2, 0.5, -0.4, 0.9, 2.2, q);

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);

  double inv = 1.0/sqrt(2.0);
  double norm[3] = { inv, inv, 0.0 };
  double tau1[3] = { -inv, inv, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double qloc[5], qback[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, q, qloc);
  eqn->rotate_to_global_func(eqn, tau1, tau2, norm, qloc, qback);

  for (int m=0; m<5; ++m)
    TEST_CHECK( gkyl_compare(qback[m], q[m], 1e-12) );

  // Density and total energy are rotation-invariant scalars.
  TEST_CHECK( gkyl_compare(qloc[0], q[0], 1e-13) );
  TEST_CHECK( gkyl_compare(qloc[4], q[4], 1e-13) );
  // Pressure (a scalar) is preserved under rotation.
  TEST_CHECK( gkyl_compare(gkyl_euler_pressure(gas_gamma, qloc),
                           gkyl_euler_pressure(gas_gamma, q), 1e-12) );

  gkyl_wv_eqn_release(eqn);
}

// Flux jump equals F(qr) - F(ql).
void
test_euler_extra_flux_jump()
{
  double gas_gamma = 1.4;
  double ql[5], qr[5];
  calcq(gas_gamma, 1.0, 0.2, 0.1, 0.0, 1.0, ql);
  calcq(gas_gamma, 2.0, -0.3, 0.5, 0.2, 2.5, qr);

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);

  double fjump[5];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);

  double fl[5], fr[5];
  gkyl_euler_flux(gas_gamma, ql, fl);
  gkyl_euler_flux(gas_gamma, qr, fr);
  for (int m=0; m<5; ++m)
    TEST_CHECK( gkyl_compare(fjump[m], fr[m]-fl[m], 1e-12) );

  gkyl_wv_eqn_release(eqn);
}

// check_inv: positive density and pressure -> valid; negative pressure invalid.
void
test_euler_extra_check_inv()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);

  double qgood[5];
  calcq(gas_gamma, 1.0, 0.1, 0.0, 0.0, 1.0, qgood);
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qgood) == true );

  // Energy too small for the kinetic part -> negative pressure -> invalid.
  double qbad[5] = { 1.0, 5.0, 0.0, 0.0, 0.1 };
  TEST_CHECK( gkyl_euler_pressure(gas_gamma, qbad) < 0.0 );
  TEST_CHECK( gkyl_wv_eqn_check_inv(eqn, qbad) == false );

  gkyl_wv_eqn_release(eqn);
}

// High-order Roe waves must sum to reconstruct the full jump delta. This is a
// fundamental property of any wave-propagation flux: sum_p W^p = qr - ql.
void
test_euler_extra_wave_sum()
{
  double gas_gamma = 1.4;
  double ql[5], qr[5];
  calcq(gas_gamma, 1.0, 0.0, 0.0, 0.0, 1.0, ql);
  calcq(gas_gamma, 0.8, 0.1, 0.05, 0.0, 0.9, qr);

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, false);

  double delta[5];
  for (int i=0; i<5; ++i) delta[i] = qr[i]-ql[i];

  double waves[3*5], speeds[3];
  double maxs = gkyl_wv_eqn_waves(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);
  TEST_CHECK( maxs > 0.0 );

  for (int i=0; i<5; ++i) {
    double sum = 0.0;
    for (int w=0; w<3; ++w) sum += waves[w*5+i];
    TEST_CHECK( gkyl_compare(sum, delta[i], 1e-11) );
  }

  gkyl_wv_eqn_release(eqn);
}

// Low-order (Lax) flux fluctuation conservation: amdq + apdq = F(qr) - F(ql).
// The Lax-Friedrichs splitting is exactly conservative by construction.
void
test_euler_extra_waves_conservation_lax()
{
  double gas_gamma = 1.4;
  double ql[5], qr[5];
  calcq(gas_gamma, 1.0, 0.0, 0.0, 0.0, 1.0, ql);
  calcq(gas_gamma, 0.8, 0.1, 0.05, 0.0, 0.9, qr);

  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_inew(&(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma, .rp_type = WV_EULER_RP_LAX, .use_gpu = false });

  double delta[5];
  for (int i=0; i<5; ++i) delta[i] = qr[i]-ql[i];

  double waves[2*5], speeds[2];
  gkyl_wv_eqn_waves(eqn, GKYL_WV_LOW_ORDER_FLUX, delta, ql, qr, 0.0, 0.0, waves, speeds);

  double amdq[5], apdq[5];
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_LOW_ORDER_FLUX, ql, qr, 0.0, 0.0, waves, speeds, amdq, apdq);

  double fjump[5];
  gkyl_wv_eqn_flux_jump(eqn, ql, qr, fjump);
  for (int i=0; i<5; ++i)
    TEST_CHECK( gkyl_compare(amdq[i]+apdq[i], fjump[i], 1e-12) );

  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "euler_extra_basic", test_euler_extra_basic },
  { "euler_extra_pressure", test_euler_extra_pressure },
  { "euler_extra_prim_vars", test_euler_extra_prim_vars },
  { "euler_extra_max_speed", test_euler_extra_max_speed },
  { "euler_extra_flux", test_euler_extra_flux },
  { "euler_extra_flux_rotation", test_euler_extra_flux_rotation },
  { "euler_extra_rotate_roundtrip", test_euler_extra_rotate_roundtrip },
  { "euler_extra_flux_jump", test_euler_extra_flux_jump },
  { "euler_extra_check_inv", test_euler_extra_check_inv },
  { "euler_extra_wave_sum", test_euler_extra_wave_sum },
  { "euler_extra_waves_conservation_lax", test_euler_extra_waves_conservation_lax },
  { NULL, NULL },
};
