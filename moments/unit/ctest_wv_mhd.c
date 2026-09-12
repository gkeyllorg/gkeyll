#include <acutest.h>
#include <float.h>
#include <gkyl_moment_prim_mhd.h>
#include <gkyl_wv_mhd.h>
#include <math.h>

void calcq(double gas_gamma, const double *pv, double *q) {
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho * u;
  q[2] = rho * v;
  q[3] = rho * w;
  q[5] = pv[5];
  q[6] = pv[6];
  q[7] = pv[7]; // B field
  double pb = 0.5 * (pv[5] * pv[5] + pv[6] * pv[6] +
                        pv[7] * pv[7]); // magnetic pressure
  q[4] = pr / (gas_gamma - 1) + 0.5 * rho * (u * u + v * v + w * w) + pb;
}

/**************************************/
/* CHECK FLUX FUNCTION IMPLEMENTATION */
/**************************************/
void test_mhd_basic_ho() {
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(&(struct gkyl_wv_mhd_inp){
      .gas_gamma = gas_gamma, .divergence_constraint = GKYL_MHD_DIVB_NONE});

  TEST_CHECK(mhd->num_equations == 8);
  TEST_CHECK(mhd->num_waves == 7);

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double bx = 0.11, by = 0.22, bz = 0.33;
  double q[8], pv[8] = {rho, u, v, w, pr, bx, by, bz};
  calcq(gas_gamma, pv, q);
  double pb = 0.5 * (bx * bx + by * by + bz * bz);
  double u_dot_b = u * bx + v * by + w * bz;
  double E = q[4];

  double fluxes[3][8] = {
      {rho * u, rho * u * u - bx * bx + pr + pb, rho * u * v - bx * by,
          rho * u * w - bx * bz, (E + pr + pb) * u - bx * u_dot_b, 0.0,
          u * by - v * bx, u * bz - w * bx},
      {rho * v, rho * v * u - bx * by, rho * v * v - by * by + pr + pb,
          rho * v * w - by * bz, (E + pr + pb) * v - by * u_dot_b,
          v * bx - u * by, 0.0, v * bz - w * by},
      {rho * w, rho * w * u - bx * bz, rho * w * v - by * bz,
          rho * w * w - bz * bz + pr + pb, (E + pr + pb) * w - bz * u_dot_b,
          w * bx - u * bz, w * by - v * bz, 0.0},
  };

  double norm[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  double tau1[3][3] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

  double tau2[3][3] = {{0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}};

  TEST_CHECK(gkyl_compare(pr, gkyl_mhd_pressure(gas_gamma, q), 1e-14));

  double q_local[8], flux_local[8], flux[8];
  for (int d = 1; d < 2; ++d) {
    mhd->rotate_to_local_func(mhd, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_mhd_flux(gas_gamma, q_local, flux_local);
    mhd->rotate_to_global_func(mhd, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int m = 0; m < 8; ++m)
      TEST_CHECK(gkyl_compare(flux[m], fluxes[d][m], 1e-14));

    // check Riemann transform
    double w1[8], q1[8];
    mhd->cons_to_riem(mhd, q_local, q_local, w1);
    mhd->riem_to_cons(mhd, q_local, w1, q1);

    for (int m = 0; m < 8; ++m)
      TEST_CHECK(gkyl_compare_double(q_local[m], q1[m], 1e-14));
  }

  double q_l[8], q_g[8];
  for (int d = 1; d < 3; ++d) {
    gkyl_wv_eqn_rotate_to_local(mhd, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m = 0; m < 8; ++m)
      TEST_CHECK(gkyl_compare(q[m], q_g[m], 1e-14));
  }

  gkyl_wv_eqn_release(mhd);
}

/*********************************************************************/
/* CHECK IF SUM OF LEFT/RIGHT GOING FLUCTUATIONS SUM TO JUMP IN FLUX */
/*********************************************************************/
void do_test_mhd_qfluct(enum gkyl_wv_mhd_rp rp_type,
    enum gkyl_wv_flux_type ftype, enum gkyl_wv_mhd_div_constraint divb,
    const double vl[], const double vr[], const int d, const double eps) {
  double gas_gamma = 5.0 / 3.0;
  double ch = 1.2345;
  struct gkyl_wv_eqn *eqn =
      gkyl_wv_mhd_new(&(struct gkyl_wv_mhd_inp){.rp_type = rp_type,
          .gas_gamma = gas_gamma,
          .divergence_constraint = divb,
          .glm_ch = ch});

  int meq = eqn->num_equations;
  int mwv = eqn->num_waves;

  double norm[3][3] = {{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}};

  double tau1[3][3] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

  double tau2[3][3] = {{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};

  double ql[meq], qr[meq];
  calcq(gas_gamma, vl, ql);
  calcq(gas_gamma, vr, qr);
  if (divb == GKYL_MHD_DIVB_GLM) {
    ql[8] = vl[8];
    qr[8] = vr[8];
  }

  double ql_local[meq], qr_local[meq];
  double speeds[mwv], waves[mwv * meq], waves_local[mwv * meq], delta[meq];
  double apdq[meq], amdq[meq], apdq_local[meq], amdq_local[meq];
  double fl[meq], fr[meq], fl_local[meq], fr_local[meq];

  // rotate left/right states to local tangent-normal frame
  gkyl_wv_eqn_rotate_to_local(eqn, tau1[d], tau2[d], norm[d], ql, ql_local);
  gkyl_wv_eqn_rotate_to_local(eqn, tau1[d], tau2[d], norm[d], qr, qr_local);

  // compute waves in local frame
  for (int i = 0; i < meq; ++i)
    delta[i] = qr_local[i] - ql_local[i];
  gkyl_wv_eqn_waves(eqn, ftype, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

  // compute left/right-going fluctuations in local frame
  gkyl_wv_eqn_qfluct(eqn, ftype, ql_local, qr_local, 1.0, 1.0, waves_local, speeds,
      amdq_local, apdq_local);

  // compute fluxes in local frame
  if (divb == GKYL_MHD_DIVB_GLM) {
    gkyl_glm_mhd_flux(gas_gamma, ch, ql_local, fl_local);
    gkyl_glm_mhd_flux(gas_gamma, ch, qr_local, fr_local);
  } else {
    gkyl_mhd_flux(gas_gamma, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, qr_local, fr_local);
  }

  // rotate local-frame waves back to global frame
  for (int mw = 0; mw < mwv; ++mw)
    gkyl_wv_eqn_rotate_to_global(eqn, tau1[d], tau2[d], norm[d],
        &waves_local[mw * meq], &waves[mw * meq]);

  // rotate local-frame fluctuations back to global frame
  gkyl_wv_eqn_rotate_to_global(
      eqn, tau1[d], tau2[d], norm[d], amdq_local, amdq);
  gkyl_wv_eqn_rotate_to_global(
      eqn, tau1[d], tau2[d], norm[d], apdq_local, apdq);

  // rotate local-frame fluxes back to global frame
  gkyl_wv_eqn_rotate_to_global(eqn, tau1[d], tau2[d], norm[d], fl_local, fl);
  gkyl_wv_eqn_rotate_to_global(eqn, tau1[d], tau2[d], norm[d], fr_local, fr);

  // check if sum of left/right-going fluctuations sum to jump in flux
  for (int i = 0; i < meq; ++i)
    TEST_CHECK(gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], eps));
  for (int i = 0; i < meq; ++i)
    TEST_CHECK(gkyl_compare(
        fr_local[i] - fl_local[i], amdq_local[i] + apdq_local[i], eps));

  gkyl_wv_eqn_release(eqn);
}

void test_mhd_qfluct_lax_ho() {
  // jumps in bx, by, and bz; checking all three directions
  double vl[8] = {1.0, 0.1, 0.2, 0.3, 1.5, 0.4, 0.43, 0.3};
  double vr[8] = {1.1, 0.13, 0.25, 0.34, 15.0, 0.42, 0.4, -0.3};
  int rp_type = WV_MHD_RP_LAX;
  int ftype = GKYL_WV_LOW_ORDER_FLUX;
  int divb = GKYL_MHD_DIVB_NONE;
  double eps = 1e-13;
  for (int d = 0; d < 3; ++d)
    do_test_mhd_qfluct(rp_type, ftype, divb, vl, vr, d, eps);
}

void test_mhd_qfluct_roe_ho() {
  // no jump in bx
  double vl[8] = {1.0, 0.1, 0.2, 0.3, 1.5, 0.4, 0.4, 0.3};
  double vr[8] = {1.1, 0.13, 0.25, 0.34, 1.54, 0.4, 0.44, 0.34};
  int rp_type = WV_MHD_RP_ROE;
  int ftype = GKYL_WV_HIGH_ORDER_FLUX;
  int divb = GKYL_MHD_DIVB_NONE;
  int d = 0;
  double eps = 1e-13;
  do_test_mhd_qfluct(WV_MHD_RP_ROE, ftype, divb, vl, vr, d, eps);
}

void test_mhd_qfluct_hlld_ho() {
  // no jump in bx
  double vl[8] = {1.0, 0.1, 0.2, 0.3, 1.5, 0.4, 0.4, 0.3};
  double vr[8] = {1.1, 0.13, 0.25, 0.34, 1.54, 0.4, 0.44, 0.34};
  int rp_type = WV_MHD_RP_HLLD;
  int ftype = GKYL_WV_HIGH_ORDER_FLUX;
  int divb = GKYL_MHD_DIVB_NONE;
  int d = 0;
  double eps = 1e-12; // FIXME 1e-13 fails
  do_test_mhd_qfluct(rp_type, ftype, divb, vl, vr, d, eps);
}

void test_glm_mhd_qfluct_lax_ho() {
  // jumps in bx, by, and bz; checking all three directions
  double vl[9] = {1.0, 0.1, 0.2, 0.3, 1.5, 0.4, 0.5, 0.2, 0.0};
  double vr[9] = {1.1, 0.13, 0.25, 0.34, 15.0, 0.42, 0.4, -0.3, 0.1};
  int rp_type = WV_MHD_RP_LAX;
  int ftype = GKYL_WV_LOW_ORDER_FLUX;
  int divb = GKYL_MHD_DIVB_GLM;
  double eps = 1e-13;
  for (int d = 0; d < 3; ++d)
    do_test_mhd_qfluct(rp_type, ftype, divb, vl, vr, d, eps);
}

void test_glm_mhd_qfluct_roe_ho() {
  // no jump in bx
  double vl[9] = {1.0, 0.1, 0.2, 0.3, 1.5, 0.4, 0.5, 0.2, 0.0};
  double vr[9] = {1.1, 0.13, 0.25, 0.34, 15.0, 0.4, 0.4, -0.3, 0.1};
  int rp_type = WV_MHD_RP_ROE;
  int ftype = GKYL_WV_HIGH_ORDER_FLUX;
  int divb = GKYL_MHD_DIVB_GLM;
  int d = 0;
  double eps = 1e-13;
  do_test_mhd_qfluct(rp_type, ftype, divb, vl, vr, d, eps);
}

void test_glm_mhd_qfluct_hlld_ho() {
  // no jump in bx
  double vl[9] = {1.0, 0.1, 0.2, 0.3, 1.5, 0.4, 0.5, 0.2, 0.0};
  double vr[9] = {1.1, 0.13, 0.25, 0.34, 15.0, 0.4, 0.4, -0.3, 0.1};
  int rp_type = WV_MHD_RP_HLLD;
  int ftype = GKYL_WV_HIGH_ORDER_FLUX;
  int divb = GKYL_MHD_DIVB_GLM;
  int d = 0;
  double eps = 1e-13;
  do_test_mhd_qfluct(rp_type, ftype, divb, vl, vr, d, eps);
}

TEST_LIST = {
    {"mhd_basic_ho", test_mhd_basic_ho},
    {"mhd_qfluct_lax_ho", test_mhd_qfluct_lax_ho},
    {"mhd_qfluct_roe_ho", test_mhd_qfluct_roe_ho},
    {"mhd_qfluct_hlld_ho", test_mhd_qfluct_hlld_ho},
    {"glm_mhd_qfluct_lax_ho", test_glm_mhd_qfluct_lax_ho},
    {"glm_mhd_qfluct_roe_ho", test_glm_mhd_qfluct_roe_ho},
    {"glm_mhd_qfluct_hlld_ho", test_glm_mhd_qfluct_hlld_ho},
    {NULL, NULL},
};
