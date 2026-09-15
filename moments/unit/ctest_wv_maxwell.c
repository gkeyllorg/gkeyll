#include <acutest.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_maxwell_priv.h>

// Make indexing cleaner with the dir_shuffle
#define EX 0
#define EY 1
#define EZ 2
#define BX 3
#define BY 4
#define BZ 5

void
test_maxwell_basic_ho()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 2.0;
  double b_fact = 2.5;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, e_fact, b_fact, false);

  TEST_CHECK( maxwell->num_equations == 8 );
  TEST_CHECK( maxwell->num_waves == 6 );

  double Ex = 1.0, Ey = 0.1, Ez = 0.2;
  double Bx = 10.0, By = 10.1, Bz = 10.2;
  double phi = 0.01, psi = 0.02;
  double q[8] = { Ex, Ey, Ez, Bx, By, Bz, phi, psi };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 }
  };

  double fluxes[3][8] = {
    {
      e_fact*c2*q[6],
      c2*q[BZ],
      -c2*q[BY],
      
      b_fact*q[7],
      -q[EZ],
      q[EY],
      
      e_fact*q[EX],
      b_fact*c2*q[BX]
    },
    {
      -c2*q[BZ],
      e_fact*c2*q[6],
      c2*q[BX],
      
      q[EZ],
      b_fact*q[7],
      -q[EX],
      
      e_fact*q[EY],
      b_fact*c2*q[BY]
    },
    {
      c2*q[BY],
      -c2*q[BX],      
      e_fact*c2*q[6],
      
      -q[EY],
      q[EX],
      b_fact*q[7],
      
      e_fact*q[EZ],
      b_fact*c2*q[BZ]
    },

  };


  double q_local[8], flux_local[8], flux[8];

  for (int d=0; d<3; ++d) {
    maxwell->rotate_to_local_func(maxwell, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_maxwell_flux(c, e_fact, b_fact, q_local, flux_local);
    maxwell->rotate_to_global_func(maxwell, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int m=0; m<8; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) );

    // check Riemann transform
    double w1[8], q1[8];
    maxwell->cons_to_riem(maxwell, q_local, q_local, w1);
    maxwell->riem_to_cons(maxwell, q_local, w1, q1);
    
    for (int m=0; m<8; ++m)
      TEST_CHECK( gkyl_compare_double(q_local[m], q1[m], 1e-14) );    
  }

  gkyl_wv_eqn_release(maxwell);
}

void
test_maxwell_waves_ho()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 1.0;
  double b_fact = 1.0;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, e_fact, b_fact, false);

  double ql[8] = { 0.0, 1.0, 0.0, 1.0, -0.75, 0.0, 0.0, 0.0};
  double qr[8] = { 0.0, -1.0, 0.0, 1.0, 0.75, 0.0, 0.0, 0.0};

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<3; ++d) {
    double speeds[6], waves[6*8], waves_local[6*8];
    double ql_local[8], qr_local[8];
    
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(maxwell, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(maxwell, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i=0; i<8; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(maxwell, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[8], amdq_local[8];
    gkyl_wv_eqn_qfluct(maxwell, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds,
      amdq_local, apdq_local);

    // rotate waves back to global frame
    for (int mw=0; mw<6; ++mw)
      gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], &waves_local[mw*8], &waves[mw*8]);

    double apdq[8], amdq[8];
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], amdq_local, amdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[8], fr_local[8];
    gkyl_maxwell_flux(c, e_fact, b_fact, ql_local, fl_local);
    gkyl_maxwell_flux(c, e_fact, b_fact, qr_local, fr_local);
    
    double fl[8], fr[8];
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(maxwell, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<8; ++i) {
      //printf("%d: %g %g (%g)\n", i, fr[i]-fl[i], amdq[i]+apdq[i], (fr[i]-fl[i])-(amdq[i]+apdq[i]));
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
    }
  }
    
  gkyl_wv_eqn_release(maxwell);
}

#ifdef GKYL_HAVE_CUDA

int cu_wv_maxwell_test(const struct gkyl_wv_eqn *eqn);

void
test_wv_maxwell_dev()
{
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 2.0;
  double b_fact = 2.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, e_fact, b_fact, true);
  
  // call CUDA test
  int nfail = cu_wv_maxwell_test(eqn->on_dev);

  TEST_CHECK( nfail == 0 );

  gkyl_wv_eqn_release(eqn);
}

#endif

TEST_LIST = {
  { "maxwell_basic_ho", test_maxwell_basic_ho },
  { "maxwell_waves_ho", test_maxwell_waves_ho },
#ifdef GKYL_HAVE_CUDA
  { "wv_maxwell_dev", test_wv_maxwell_dev },
#endif  
  { NULL, NULL },
};
