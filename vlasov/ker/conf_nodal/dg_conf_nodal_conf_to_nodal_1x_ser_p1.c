#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_conf_to_nodal_1x_ser_p1(const double *fmodal, double* GKYL_RESTRICT fnodal) 
{ 
  // fmodal: p=1 configuration-space expansion (2 coefficients). 
  // fnodal: its values at the 2 interior Gauss-Legendre nodes (row-major). 
  const double t0[2] = {fmodal[0], fmodal[1]};
  const double t1[2] = {t0[0]-t0[1], t0[0]+t0[1]};
  fnodal[0] = 0.7071067811865475*t1[0]; 
  fnodal[1] = 0.7071067811865475*t1[1]; 
} 
