#include <gkyl_lbo_pkpm_nodal_kernels.h> 
GKYL_CU_DH void lbo_pkpm_nodal_vol_1v_p2(const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[1]:    Velocity cell centre. 
  // dxv[1]:  Velocity cell width. 
  // nu:      Collision frequency at the node. 
  // nu_u:    Collision frequency times the drag velocity at the node. 
  // nu_vtsq: Collision frequency times the thermal speed squared at the node. 
  // stride:  Offset between the blocks [F_0, G_1, M0_perp (3), N_perp (3)] of a cell. 
  // f:       Velocity expansions of the blocks at the node. 
  // out:     Incremented output. 
  const double rdv = 2.0/dxv[0]; 
  const double a0 = rdv*(nu*w[0] - nu_u), a1 = nu, d2 = nu_vtsq*rdv*rdv; 
  for (int b = 0; b < 8; ++b) { 
    const double *h = &f[b*stride]; 
    double *o = &out[b*stride]; 
    o[1] += -1.7320508075688772*a0*h[0] - 1.0*a1*h[1]; 
    o[2] += -3.872983346207417*a0*h[1] - 2.23606797749979*a1*h[0] - 2.0*a1*h[2] + 6.708203932499369*d2*h[0]; 
  } 
  const double *F_0 = &f[0], *G_1 = &f[stride]; 
  double *out_G_1 = &out[stride]; 
  for (int j = 0; j < 3; ++j) 
    out_G_1[j] += 2.0*(nu_vtsq*F_0[j] - nu*G_1[j]); 
  for (int k = 0; k < 3; ++k) { 
    const double *M0 = &f[(2 + k)*stride], *N = &f[(5 + k)*stride]; 
    double *out_M0 = &out[(2 + k)*stride], *out_N = &out[(5 + k)*stride]; 
    for (int j = 0; j < 3; ++j) { 
      out_M0[j] -= nu*M0[j]; 
      out_N[j] += 4.0*nu_vtsq*M0[j] - 3.0*nu*N[j]; 
    } 
  } 
} 
