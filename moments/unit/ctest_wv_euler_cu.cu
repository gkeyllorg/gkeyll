/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_priv.h>
int cu_wv_euler_test(const struct gkyl_wv_eqn *eqn);
}

__global__
void ker_cu_wv_euler_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 5, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 3, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);

  GKYL_CU_CHECK( euler->gas_gamma == 1.4, nfail );
  double gas_gamma = euler->gas_gamma;
  
  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 0.0;
  double q[5];

  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v+w*w);
  double E = q[4];

  double fluxes[3][5] = {
    { rho*u, rho*u*u+pr, rho*u*v, rho*u*w, (E+pr)*u },
    { rho*v, rho*u*v, rho*v*v+pr, rho*v*w, (E+pr)*v },
    { rho*w, rho*u*w, rho*v*w, rho*w*w, (E+pr)*w },
  };

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

  GKYL_CU_CHECK( pr == gkyl_euler_pressure(gas_gamma, q), nfail );

  double q_local[5], flux_local[5], flux[5];
  for (int d=0; d<3; ++d) {
    eqn->rotate_to_local_func(eqn, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_flux(gas_gamma, q_local, flux_local);
    eqn->rotate_to_global_func(eqn, tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<5; ++m)
      GKYL_CU_CHECK( flux[m] == fluxes[d][m], nfail );
  }

  double q_l[5], q_g[5];
  for (int d=0; d<3; ++d) {
    eqn->rotate_to_local_func(eqn, tau1[d], tau2[d], norm[d], q, q_l);
    eqn->rotate_to_global_func(eqn, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<5; ++m) GKYL_CU_CHECK( q[m] == q_g[m], nfail );
  }
}

int cu_wv_euler_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_wv_euler_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
