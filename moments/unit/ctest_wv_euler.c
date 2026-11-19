#include <gkyl_util.h>
#include <acutest.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_priv.h>

void
calcq(double gas_gamma, const double pv[5], double q[5])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho * u;
  q[2] = rho * v;
  q[3] = rho * w;
  q[4] = pr / (gas_gamma - 1) + 0.5 * rho * (u * u + v * v + w * w);
}

void
test_euler_basic_ho()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(gas_gamma, false);

  TEST_CHECK(euler->num_equations == 5);
  TEST_CHECK(euler->num_waves == 3);

  //double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 0.0;
  double q[5], pv[5] = {rho, u, v, w, pr};

  calcq(gas_gamma, pv, q);
  double E = q[4];

  double fluxes[3][5] = {
    {rho * u, rho * u * u + pr, rho * u * v, rho * u * w, (E + pr) * u},
    {rho * v, rho * u * v, rho * v * v + pr, rho * v * w, (E + pr) * v},
    {rho * w, rho * u * w, rho * v * w, rho * w * w, (E + pr) * w}
  };

  double norm[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  double tau1[3][3] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

  double tau2[3][3] = {{0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}};

  TEST_CHECK(pr == gkyl_euler_pressure(gas_gamma, q));

  double q_local[5], flux_local[5], flux[5];
  for (int d = 0; d < 3; ++d) {
    euler->rotate_to_local_func(euler, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_flux(gas_gamma, q_local, flux_local);
    euler->rotate_to_global_func(euler, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int m = 0; m < 5; ++m) {
      TEST_CHECK(gkyl_compare(flux[m], fluxes[d][m], 1e-15));
    }
  }

  double q_l[5], q_g[5];
  for (int d = 0; d < 3; ++d) {
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m = 0; m < 5; ++m) {
      TEST_CHECK(q[m] == q_g[m]);
    }

    // check Riemann transform
    double w1[5], q1[5];
    euler->cons_to_riem(euler, q_local, q_local, w1);
    euler->riem_to_cons(euler, q_local, w1, q1);

    for (int m = 0; m < 5; ++m) {
      TEST_CHECK(gkyl_compare_double(q_local[m], q1[m], 1e-14));
    }
  }

  gkyl_wv_eqn_release(euler);
}

void
test_euler_waves(enum gkyl_wv_flux_type ftype)
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(gas_gamma, false);

  double vl[5] = {1.0, 0.1, 0.2, 0.3, 1.5};
  double vr[5] = {0.1, 1.0, 2.0, 3.0, 0.15};

  double ql[5], qr[5];
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql);
  calcq(gas_gamma, vr, qr);

  double norm[3][3] = {{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}};

  double tau1[3][3] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

  double tau2[3][3] = {{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};

  for (int d = 0; d < 3; ++d) {
    double speeds[3], waves[3 * 5], waves_local[3 * 5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[5];
    for (int i = 0; i < 5; ++i) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler, ftype, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[5], amdq_local[5];
    gkyl_wv_eqn_qfluct(
      euler, ftype, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local
    );

    // rotate waves back to global frame
    for (int mw = 0; mw < 3; ++mw) {
      gkyl_wv_eqn_rotate_to_global(
        euler, tau1[d], tau2[d], norm[d], &waves_local[mw * 5], &waves[mw * 5]
      );
    }

    double apdq[5], amdq[5];
    // rotate fluctuations back to global frame
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[5], fr_local[5];
    gkyl_euler_flux(gas_gamma, ql_local, fl_local);
    gkyl_euler_flux(gas_gamma, qr_local, fr_local);

    double fl[5], fr[5];
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 5; ++i) {
      TEST_CHECK(gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14));
    }
  }

  gkyl_wv_eqn_release(euler);
}

void
test_euler_waves_hof_ho(void)
{
  test_euler_waves(GKYL_WV_HIGH_ORDER_FLUX);
}
void
test_euler_waves_lof_ho(void)
{
  test_euler_waves(GKYL_WV_LOW_ORDER_FLUX);
}

void
test_euler_waves_2(enum gkyl_wv_flux_type ftype, enum gkyl_wv_euler_rp rp_type)
{
  double gas_gamma = 1.4;
  struct gkyl_wv_euler_inp inp = {.gas_gamma = gas_gamma, .rp_type = rp_type, .use_gpu = false};
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_inew(&inp);

  double vl[5] = {1.0, 0.1, 0.2, 0.3, 1.5};
  double vr[5] = {0.01, 1.0, 2.0, 3.0, 15.0};

  double ql[5], qr[5];
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql);
  calcq(gas_gamma, vr, qr);

  double norm[3][3] = {{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}};

  double tau1[3][3] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

  double tau2[3][3] = {{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};

  for (int d = 0; d < 3; ++d) {
    double speeds[3], waves[3 * 5], waves_local[3 * 5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[5];
    for (int i = 0; i < 5; ++i) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(euler, ftype, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[5], amdq_local[5];
    gkyl_wv_eqn_qfluct(
      euler, ftype, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local
    );

    // rotate waves back to global frame
    for (int mw = 0; mw < 3; ++mw) {
      gkyl_wv_eqn_rotate_to_global(
        euler, tau1[d], tau2[d], norm[d], &waves_local[mw * 5], &waves[mw * 5]
      );
    }

    double apdq[5], amdq[5];
    // rotate fluctuations back to global frame
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[5], fr_local[5];
    gkyl_euler_flux(gas_gamma, ql_local, fl_local);
    gkyl_euler_flux(gas_gamma, qr_local, fr_local);

    double fl[5], fr[5];
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 5; ++i) {
      TEST_CHECK(gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12));
    }
  }

  gkyl_wv_eqn_release(euler);
}

#ifdef GKYL_HAVE_CUDA

int cu_wv_euler_test(const struct gkyl_wv_eqn *eqn);

void
test_wv_euler_dev()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *eqn = gkyl_wv_euler_new(gas_gamma, true);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);

  TEST_CHECK(euler->gas_gamma == 1.4);

  // call CUDA test
  int nfail = cu_wv_euler_test(eqn->on_dev);

  TEST_CHECK(nfail == 0);

  gkyl_wv_eqn_release(eqn);
}

#endif

// HLL and HLLC wave speeds: ordered, finite, and with a positive HLLC middle density, for
// colliding streams, two strong shocks and a near vacuum; the fluctuations sum to the flux jump,
// and an isolated contact is exact for HLLC.
static void
test_euler_hll_speeds(enum gkyl_wv_euler_rp rp_type)
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_inew(&(struct gkyl_wv_euler_inp){
    .gas_gamma = gas_gamma,
    .rp_type = rp_type,
    .use_gpu = false,
  });
  int num_waves = euler->num_waves;

  // primitive states (rho, u, v, w, p)
  double vl[][5] = {
    {1.0, 2.5, 0.0, 0.0, 1.0}, // colliding streams, Mach 2.1 each
    {1.0, 2.0897846388247809, 0.0, 0.0, 1.0}, // colliding streams, Mach 1.77 each
    {5.99924, 19.5975, 0.0, 0.0, 460.894}, // two strong shocks (Toro test 5)
    {5.99924, 19.5975, 0.0, 0.0, 460.894}, // two strong shocks, slower right stream
    {1.0, 0.0, 0.0, 0.0, 1.0}, // near vacuum on the right
    {1.4, 0.3, 0.2, -0.1, 1.0} // isolated moving contact with shear
  };
  double vr[][5] = {
    {1.0, -2.5, 0.0, 0.0, 1.0},
    {1.0, -2.0897846388247809, 0.0, 0.0, 1.0},
    {5.99242, -6.19633, 0.0, 0.0, 46.0950},
    {5.99242, -6.0170177200762591, 0.0, 0.0, 46.0950},
    {1.0e-8, 0.0, 0.0, 0.0, 1.0e-10},
    {1.0, 0.3, -0.4, 0.5, 1.0}
  };
  int num_cases = sizeof(vl) / sizeof(vl[0]);

  for (int n = 0; n < num_cases; ++n) {
    double ql[5], qr[5], fl[5], fr[5], delta[5];
    const double *prim[2] = {vl[n], vr[n]};
    double *cons[2] = {ql, qr}, *flux[2] = {fl, fr};
    for (int side = 0; side < 2; ++side) {
      const double *v = prim[side];
      double *q = cons[side];
      q[0] = v[0];
      q[1] = v[0] * v[1];
      q[2] = v[0] * v[2];
      q[3] = v[0] * v[3];
      q[4] = v[4] / (gas_gamma - 1) + 0.5 * v[0] * (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
      gkyl_euler_flux(gas_gamma, q, flux[side]);
    }
    for (int i = 0; i < 5; ++i) {
      delta[i] = qr[i] - ql[i];
    }

    double waves[15] = {0.0}, speeds[3] = {0.0}, amdq[5], apdq[5];
    gkyl_wv_eqn_waves(euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 1.0, 1.0, waves, speeds);
    gkyl_wv_eqn_qfluct(euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);

    // ordered speeds that bound the characteristic speeds of both sides
    for (int m = 0; m < num_waves; ++m) {
      TEST_CHECK(isfinite(speeds[m]));
      if (m > 0) {
        TEST_CHECK(speeds[m] > speeds[m - 1]);
      }
    }
    TEST_CHECK(speeds[0] <= vl[n][1] - sqrt(gas_gamma * vl[n][4] / vl[n][0]));
    TEST_CHECK(speeds[num_waves - 1] >= vr[n][1] + sqrt(gas_gamma * vr[n][4] / vr[n][0]));

    // positive density behind every wave, and waves that sum to the jump
    double rho = ql[0];
    for (int m = 0; m < num_waves - 1; ++m) {
      rho += waves[5 * m];
      TEST_CHECK(rho > 0.0);
    }
    for (int i = 0; i < 5; ++i) {
      double wsum = 0.0;
      for (int m = 0; m < num_waves; ++m) {
        wsum += waves[5 * m + i];
      }
      TEST_CHECK(gkyl_compare(wsum, delta[i], 1e-12));
      TEST_CHECK(isfinite(amdq[i]) && isfinite(apdq[i]));
      TEST_CHECK(gkyl_compare(amdq[i] + apdq[i], fr[i] - fl[i], 1e-12));
    }
  }

  // HLLC keeps the isolated contact (the last case): its flux is the physical flux of either side
  if (rp_type == WV_EULER_RP_HLLC) {
    int n = num_cases - 1;
    double ql[5] = {vl[n][0], vl[n][0] * vl[n][1], vl[n][0] * vl[n][2], vl[n][0] * vl[n][3], 0.0};
    double qr[5] = {vr[n][0], vr[n][0] * vr[n][1], vr[n][0] * vr[n][2], vr[n][0] * vr[n][3], 0.0};
    ql[4] =
      vl[n][4] / (gas_gamma - 1) + 0.5 * (ql[1] * ql[1] + ql[2] * ql[2] + ql[3] * ql[3]) / ql[0];
    qr[4] =
      vr[n][4] / (gas_gamma - 1) + 0.5 * (qr[1] * qr[1] + qr[2] * qr[2] + qr[3] * qr[3]) / qr[0];
    double delta[5], waves[15] = {0.0}, speeds[3] = {0.0}, amdq[5], apdq[5], fl[5];
    for (int i = 0; i < 5; ++i) {
      delta[i] = qr[i] - ql[i];
    }
    gkyl_euler_flux(gas_gamma, ql, fl);
    gkyl_wv_eqn_waves(euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql, qr, 1.0, 1.0, waves, speeds);
    gkyl_wv_eqn_qfluct(euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);
    // the contact moves to the right (u = 0.3): the interface flux is that of the left state
    for (int i = 0; i < 5; ++i) {
      TEST_CHECK(gkyl_compare(fl[i] + amdq[i], fl[i], 1e-13));
    }
  }

  gkyl_wv_eqn_release(euler);
}

void
test_euler_hll_speeds_hll(void)
{
  test_euler_hll_speeds(WV_EULER_RP_HLL);
}

void
test_euler_hll_speeds_hllc(void)
{
  test_euler_hll_speeds(WV_EULER_RP_HLLC);
}

void
test_euler_waves_2_hof_roe_ho(void)
{
  test_euler_waves_2(GKYL_WV_HIGH_ORDER_FLUX, WV_EULER_RP_ROE);
}
void
test_euler_waves_2_lof_roe_ho(void)
{
  test_euler_waves_2(GKYL_WV_LOW_ORDER_FLUX, WV_EULER_RP_ROE);
}

void
test_euler_waves_2_hof_hllc_ho(void)
{
  test_euler_waves_2(GKYL_WV_HIGH_ORDER_FLUX, WV_EULER_RP_HLLC);
}
void
test_euler_waves_2_lof_hllc_ho(void)
{
  test_euler_waves_2(GKYL_WV_LOW_ORDER_FLUX, WV_EULER_RP_HLLC);
}

void
test_euler_waves_2_hof_lax_ho(void)
{
  test_euler_waves_2(GKYL_WV_HIGH_ORDER_FLUX, WV_EULER_RP_LAX);
}
void
test_euler_waves_2_lof_lax_ho(void)
{
  test_euler_waves_2(GKYL_WV_LOW_ORDER_FLUX, WV_EULER_RP_LAX);
}

void
test_euler_waves_2_hof_hll_ho(void)
{
  test_euler_waves_2(GKYL_WV_HIGH_ORDER_FLUX, WV_EULER_RP_HLL);
}
void
test_euler_waves_2_lof_hll_ho(void)
{
  test_euler_waves_2(GKYL_WV_LOW_ORDER_FLUX, WV_EULER_RP_HLL);
}

TEST_LIST = {
  {"euler_basic_ho", test_euler_basic_ho},
  {"euler_waves_hof_ho", test_euler_waves_hof_ho},
  {"euler_waves_lof_ho", test_euler_waves_lof_ho},
  {"euler_waves_2_hof_roe_ho", test_euler_waves_2_hof_roe_ho},
  {"euler_waves_2_lof_roe_ho", test_euler_waves_2_lof_roe_ho},
  {"euler_waves_2_hof_hllc_ho", test_euler_waves_2_hof_hllc_ho},
  {"euler_waves_2_lof_hllc_ho", test_euler_waves_2_lof_hllc_ho},
  {"euler_waves_2_hof_lax_ho", test_euler_waves_2_hof_lax_ho},
  {"euler_waves_2_lof_lax_ho", test_euler_waves_2_lof_lax_ho},
  {"euler_hll_speeds_hll", test_euler_hll_speeds_hll},
  {"euler_hll_speeds_hllc", test_euler_hll_speeds_hllc},
  {"euler_waves_2_hof_hll_ho", test_euler_waves_2_hof_hll_ho},
  {"euler_waves_2_lof_hll_ho", test_euler_waves_2_lof_hll_ho},
#ifdef GKYL_HAVE_CUDA
  {"wv_euler_dev", test_wv_euler_dev},
#endif
  {NULL, NULL}
};
