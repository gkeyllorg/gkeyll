/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_ten_moment_priv.h>
int cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn);
}

__global__
void ker_cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 10, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 5, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_ten_moment *ten_moment = container_of(eqn, struct wv_ten_moment, eqn);

  GKYL_CU_CHECK( ten_moment->k0 == 1.0, nfail );

  double vl[10] = { 1.0, 0.1, 0.2, 0.3, 0.5, 0.0, 0.0, 1.0, 0.0, 1.5};
  double vr[10] = { 0.1, 1.0, 2.0, 3.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.3};

  double ql[10], qr[10];
  double ql_local[10], qr_local[10];

  ql[0] = vl[0];
  ql[1] = vl[0]*vl[1]; 
  ql[2] = vl[0]*vl[2]; 
  ql[3] = vl[0]*vl[3];
  ql[4] = vl[4] + vl[0]*vl[1]*vl[1]; 
  ql[5] = vl[5] + vl[0]*vl[1]*vl[2]; 
  ql[6] = vl[6] + vl[0]*vl[1]*vl[3];
  ql[7] = vl[7] + vl[0]*vl[2]*vl[2]; 
  ql[8] = vl[8] + vl[0]*vl[2]*vl[3]; 
  ql[9] = vl[9] + vl[0]*vl[3]*vl[3];

  qr[0] = vr[0];
  qr[1] = vr[0]*vr[1]; 
  qr[2] = vr[0]*vr[2]; 
  qr[3] = vr[0]*vr[3];
  qr[4] = vr[4] + vr[0]*vr[1]*vr[1]; 
  qr[5] = vr[5] + vr[0]*vr[1]*vr[2]; 
  qr[6] = vr[6] + vr[0]*vr[1]*vr[3];
  qr[7] = vr[7] + vr[0]*vr[2]*vr[2]; 
  qr[8] = vr[8] + vr[0]*vr[2]*vr[3]; 
  qr[9] = vr[9] + vr[0]*vr[3]*vr[3];

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  
  
  for (int d=0; d<3; ++d) {
    double speeds[5], waves[5*10], waves_local[5*10];
    // rotate to local tangent-normal frame
    eqn->rotate_to_local_func(eqn, tau1[d], tau2[d], norm[d], ql, ql_local);
    eqn->rotate_to_local_func(eqn, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[10];
    for (int i=0; i<10; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<5; ++mw)
      eqn->rotate_to_global_func(eqn, tau1[d], tau2[d], norm[d], &waves_local[mw*10], &waves[mw*10]);

    double apdq[10], amdq[10];
    eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, 1.0, 1.0, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[10], fr_local[10];
    gkyl_ten_moment_flux(ql_local, fl_local);
    gkyl_ten_moment_flux(qr_local, fr_local);

    double fl[10], fr[10];
    eqn->rotate_to_global_func(eqn, tau1[d], tau2[d], norm[d], fl_local, fl);
    eqn->rotate_to_global_func(eqn, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<10; ++i)
      GKYL_CU_CHECK( fabs((fr[i]-fl[i]) - (amdq[i]+apdq[i])) < 1e-14, nfail );
  }
}

int cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_wv_ten_moment_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
