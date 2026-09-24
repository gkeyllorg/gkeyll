#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_iso_euler_mixture.h>
#include <gkyl_wv_iso_euler_mixture_priv.h>

void
test_iso_euler_mixture_twocomponent_basic()
{
  double vt1 = 1.0;
  double vt2 = 10.0;

  double *vt_s = gkyl_malloc(sizeof(double[2]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(2, vt_s, false);

  TEST_CHECK( iso_euler_mixture->num_equations == 7 );
  TEST_CHECK( iso_euler_mixture->num_waves == 2 );

  double alpha1 = 0.75, rho1 = 1.0, rho2 = 2.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3;
  double rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2);
  double vt_total = (alpha1 * vt1) + ((1.0 - alpha1) * vt2);

  double q[7] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, rho_total * alpha1, alpha1 * rho1,
    (1.0 - alpha1) * rho2 };

  double prims[7];
  gkyl_iso_euler_mixture_prim_vars(2, vt_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) ); 
  TEST_CHECK( gkyl_compare(prims[4], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[5], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho2, 1e-16) );

  double fluxes[3][7] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + (rho_total * (vt_total * vt_total)), rho_total * (vx_total * vy_total),
     rho_total * (vx_total * vz_total), rho_total * vx_total * alpha1, alpha1 * (vx_total * rho1), (1.0 - alpha1) * (vx_total * rho2) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + (rho_total * (vt_total * vt_total)),
     rho_total * (vy_total * vz_total), rho_total * vy_total * alpha1, alpha1 * (vy_total * rho1), (1.0 - alpha1) * (vy_total * rho2) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total),
     (rho_total * (vz_total * vz_total)) + (rho_total * (vt_total * vt_total)), rho_total * vz_total * alpha1, alpha1 * (vz_total * rho1),
     (1.0 - alpha1) * (vz_total * rho2) },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[7], flux_local[7], flux[7];
  for (int d = 0; d < 3; d++) {
    iso_euler_mixture->rotate_to_local_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_mixture_flux(2, vt_s, q_local, flux_local);
    iso_euler_mixture->rotate_to_global_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[7], q_g[7];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[7], q1[7];
    iso_euler_mixture->cons_to_riem(iso_euler_mixture, q_local, q_local, w1);
    iso_euler_mixture->riem_to_cons(iso_euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

void
test_iso_euler_mixture_threecomponent_basic()
{
  double vt1 = 1.0;
  double vt2 = 10.0;
  double vt3 = 100.0;

  double *vt_s = gkyl_malloc(sizeof(double[3]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;
  vt_s[2] = vt3;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(3, vt_s, false);

  TEST_CHECK( iso_euler_mixture->num_equations == 9 );
  TEST_CHECK( iso_euler_mixture->num_waves == 2 );

  double alpha1 = 0.5, alpha2 = 0.3, rho1 = 1.0, rho2 = 2.0, rho3 = 3.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3;
  double rho_total = (alpha1 * rho1) + (alpha2 * rho2) + (1.0 - (alpha1 + alpha2)) * rho3;
  double vt_total = (alpha1 * vt1) + (alpha2 * vt2) + (1.0 - (alpha1 + alpha2)) * vt3;

  double q[9] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, rho_total * alpha1, rho_total * alpha2,
    alpha1 * rho1, alpha2 * rho2, (1.0 - (alpha1 + alpha2)) * rho3 };

  double prims[9];
  gkyl_iso_euler_mixture_prim_vars(3, vt_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) ); 
  TEST_CHECK( gkyl_compare(prims[4], alpha1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[5], alpha2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[8], rho3, 1e-16) );

  double fluxes[3][9] = {
   { rho_total * vx_total, (rho_total * (vx_total * vx_total)) + (rho_total * (vt_total * vt_total)), rho_total * (vx_total * vy_total),
     rho_total * (vx_total * vz_total), rho_total * vx_total * alpha1, rho_total * vx_total * alpha2, alpha1 * (vx_total * rho1),
     alpha2 * (vx_total * rho2), (1.0 - (alpha1 + alpha2)) * (vx_total * rho3) },
   { rho_total * vy_total, rho_total * (vy_total * vx_total), (rho_total * (vy_total * vy_total)) + (rho_total * (vt_total * vt_total)),
     rho_total * (vy_total * vz_total), rho_total * vy_total * alpha1, rho_total * vy_total * alpha2, alpha1 * (vy_total * rho1),
     alpha2 * (vy_total * rho2), (1.0 - (alpha1 + alpha2)) * (vy_total * rho3) },
   { rho_total * vz_total, rho_total * (vz_total * vx_total), rho_total * (vz_total * vy_total),
     (rho_total * (vz_total * vz_total)) + (rho_total * (vt_total * vt_total)), rho_total * vz_total * alpha1, rho_total * vz_total * alpha2,
     alpha1 * (vz_total * rho1), alpha2 * (vz_total * rho2), (1.0 - (alpha1 + alpha2)) * (vz_total * rho3) },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[9], flux_local[9], flux[9];
  for (int d = 0; d < 3; d++) {
    iso_euler_mixture->rotate_to_local_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_mixture_flux(3, vt_s, q_local, flux_local);
    iso_euler_mixture->rotate_to_global_func(iso_euler_mixture, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[9], q_g[9];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[9], q1[9];
    iso_euler_mixture->cons_to_riem(iso_euler_mixture, q_local, q_local, w1);
    iso_euler_mixture->riem_to_cons(iso_euler_mixture, q_local, w1, q1);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

void
test_iso_euler_mixture_twocomponent_waves()
{
  double vt1 = 1.0;
  double vt2 = 10.0;
  
  double *vt_s = gkyl_malloc(sizeof(double[2]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(2, vt_s, false);
  
  double alpha1_l = 0.75, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3;
  double rho_total_l = (alpha1_l * rho1_l) + ((1.0 - alpha1_l) * rho2_l);
  double vt_total_l = (alpha1_l * vt1) + ((1.0 - alpha1_l) * vt2);

  double alpha1_r = 0.25, rho1_r = 0.1, rho2_r = 0.2, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0;
  double rho_total_r = (alpha1_r * rho1_r) + ((1.0 - alpha1_r) * rho2_r);
  double vt_total_r = (alpha1_r * vt1) + ((1.0 - alpha1_r) * vt2);

  double ql[7] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, rho_total_l * alpha1_l,
    alpha1_l * rho1_l, (1.0 - alpha1_l) * rho2_l };
  double qr[7] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, rho_total_r * alpha1_r,
    alpha1_r * rho1_r, (1.0 - alpha1_r) * rho2_r };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 7], waves_local[2 * 7];

    double ql_local[7], qr_local[7];
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[7];
    for (int i = 0; i < 7; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[7], amdq_local[7];
    gkyl_wv_eqn_qfluct(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 7], &waves[i * 7]);
    }

    double apdq[7], amdq[7];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[7], fr_local[7];
    gkyl_iso_euler_mixture_flux(2, vt_s, ql_local, fl_local);
    gkyl_iso_euler_mixture_flux(2, vt_s, qr_local, fr_local);

    double fl[7], fr[7];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

void
test_iso_euler_mixture_twocomponent_waves_2()
{
  double vt1 = 0.5;
  double vt2 = 50.0;
  
  double *vt_s = gkyl_malloc(sizeof(double[2]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(2, vt_s, false);
  
  double alpha1_l = 0.9, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3;
  double rho_total_l = (alpha1_l * rho1_l) + ((1.0 - alpha1_l) * rho2_l);
  double vt_total_l = (alpha1_l * vt1) + ((1.0 - alpha1_l) * vt2);

  double alpha1_r = 0.1, rho1_r = 0.01, rho2_r = 0.02, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0;
  double rho_total_r = (alpha1_r * rho1_r) + ((1.0 - alpha1_r) * rho2_r);
  double vt_total_r = (alpha1_r * vt1) + ((1.0 - alpha1_r) * vt2);

  double ql[7] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, rho_total_l * alpha1_l,
    alpha1_l * rho1_l, (1.0 - alpha1_l) * rho2_l };
  double qr[7] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, rho_total_r * alpha1_r,
    alpha1_r * rho1_r, (1.0 - alpha1_r) * rho2_r };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 7], waves_local[2 * 7];

    double ql_local[7], qr_local[7];
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[7];
    for (int i = 0; i < 7; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[7], amdq_local[7];
    gkyl_wv_eqn_qfluct(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 7], &waves[i * 7]);
    }

    double apdq[7], amdq[7];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[7], fr_local[7];
    gkyl_iso_euler_mixture_flux(2, vt_s, ql_local, fl_local);
    gkyl_iso_euler_mixture_flux(2, vt_s, qr_local, fr_local);

    double fl[7], fr[7];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 7; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

void
test_iso_euler_mixture_threecomponent_waves()
{
  double vt1 = 1.0;
  double vt2 = 10.0;
  double vt3 = 100.0;
  
  double *vt_s = gkyl_malloc(sizeof(double[3]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;
  vt_s[2] = vt3;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(3, vt_s, false);
  
  double alpha1_l = 0.5, alpha2_l = 0.3, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3;
  double rho_total_l = (alpha1_l * rho1_l) + (alpha2_l * rho2_l) + ((1.0 - (alpha1_l + alpha2_l)) * rho3_l);
  double vt_total_l = (alpha1_l * vt1) + (alpha2_l * vt2) + ((1.0 - (alpha1_l + alpha2_l)) * vt3);

  double alpha1_r = 0.4, alpha2_r = 0.2, rho1_r = 0.1, rho2_r = 0.2, rho3_r = 0.3, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0;
  double rho_total_r = (alpha1_r * rho1_r) + (alpha2_r * rho2_r) + ((1.0 - (alpha1_r + alpha2_r)) * rho3_r);
  double vt_total_r = (alpha1_r * vt1) + (alpha2_r * vt2) + ((1.0 - (alpha1_r + alpha2_r)) * vt3);

  double ql[9] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, rho_total_l * alpha1_l,
    rho_total_l * alpha2_l, alpha1_l * rho1_l, alpha2_l * rho2_l, (1.0 - (alpha1_l + alpha2_l)) * rho3_l };
  double qr[9] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, rho_total_r * alpha1_r,
    rho_total_r * alpha2_r, alpha1_r * rho1_r, alpha2_r * rho2_r, (1.0 - (alpha1_r + alpha2_r)) * rho3_r };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 9], waves_local[2 * 9];

    double ql_local[9], qr_local[9];
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i = 0; i < 9; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 9], &waves[i * 9]);
    }

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[9], fr_local[9];
    gkyl_iso_euler_mixture_flux(3, vt_s, ql_local, fl_local);
    gkyl_iso_euler_mixture_flux(3, vt_s, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

void
test_iso_euler_mixture_threecomponent_waves_2()
{
  double vt1 = 0.5;
  double vt2 = 50.0;
  double vt3 = 800.0;
  
  double *vt_s = gkyl_malloc(sizeof(double[3]));
  vt_s[0] = vt1;
  vt_s[1] = vt2;
  vt_s[2] = vt3;

  struct gkyl_wv_eqn *iso_euler_mixture = gkyl_wv_iso_euler_mixture_new(3, vt_s, false);
  
  double alpha1_l = 0.8, alpha2_l = 0.1, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3;
  double rho_total_l = (alpha1_l * rho1_l) + (alpha2_l * rho2_l) + ((1.0 - (alpha1_l + alpha2_l)) * rho3_l);
  double vt_total_l = (alpha1_l * vt1) + (alpha2_l * vt2) + ((1.0 - (alpha1_l + alpha2_l)) * vt3);

  double alpha1_r = 0.75, alpha2_r = 0.2, rho1_r = 0.01, rho2_r = 0.02, rho3_r = 0.03, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0;
  double rho_total_r = (alpha1_r * rho1_r) + (alpha2_r * rho2_r) + ((1.0 - (alpha1_r + alpha2_r)) * rho3_r);
  double vt_total_r = (alpha1_r * vt1) + (alpha2_r * vt2) + ((1.0 - (alpha1_r + alpha2_r)) * vt3);

  double ql[9] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, rho_total_l * alpha1_l,
    rho_total_l * alpha2_l, alpha1_l * rho1_l, alpha2_l * rho2_l, (1.0 - (alpha1_l + alpha2_l)) * rho3_l };
  double qr[9] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, rho_total_r * alpha1_r,
    rho_total_r * alpha2_r, alpha1_r * rho1_r, alpha2_r * rho2_r, (1.0 - (alpha1_r + alpha2_r)) * rho3_r };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 9], waves_local[2 * 9];

    double ql_local[9], qr_local[9];
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(iso_euler_mixture, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i = 0; i < 9; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(iso_euler_mixture, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], &waves_local[i * 9], &waves[i * 9]);
    }

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[9], fr_local[9];
    gkyl_iso_euler_mixture_flux(3, vt_s, ql_local, fl_local);
    gkyl_iso_euler_mixture_flux(3, vt_s, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(iso_euler_mixture, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(iso_euler_mixture);
  gkyl_free(vt_s);
}

TEST_LIST = {
  { "iso_euler_mixture_twocomponent_basic", test_iso_euler_mixture_twocomponent_basic },
  { "iso_euler_mixture_threecomponent_basic", test_iso_euler_mixture_threecomponent_basic },
  { "iso_euler_mixture_twocomponent_waves", test_iso_euler_mixture_twocomponent_waves },
  { "iso_euler_mixture_twocomponent_waves_2", test_iso_euler_mixture_twocomponent_waves_2 },
  { "iso_euler_mixture_threecomponent_waves", test_iso_euler_mixture_threecomponent_waves },
  { "iso_euler_mixture_threecomponent_waves_2", test_iso_euler_mixture_threecomponent_waves_2 },
  { NULL, NULL },
};
