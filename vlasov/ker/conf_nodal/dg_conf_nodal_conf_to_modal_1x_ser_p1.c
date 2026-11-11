#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_conf_to_modal_1x_ser_p1(const double *fnodal, double* GKYL_RESTRICT fmodal) 
{ 
  // fnodal: values at the 2 interior Gauss-Legendre nodes. 
  // fmodal: the p=1 configuration-space expansion with those values. 
  const double t0[2] = {fnodal[0], fnodal[1]};
  const double t1[2] = {t0[0]+t0[1], t0[1]-t0[0]};
  fmodal[0] = 0.7071067811865475*t1[0]; 
  fmodal[1] = 0.7071067811865475*t1[1]; 
} 
