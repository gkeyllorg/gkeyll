#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_conf_to_modal_2x_ser_p1(const double *fnodal, double* GKYL_RESTRICT fmodal) 
{ 
  // fnodal: values at the 4 interior Gauss-Legendre nodes. 
  // fmodal: the p=1 configuration-space expansion with those values. 
  const double t0[4] = {fnodal[0], fnodal[1], fnodal[2], fnodal[3]};
  const double t1[4] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2]};
  const double t2[4] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1]};
  fmodal[0] = 0.5*t2[0]; 
  fmodal[1] = 0.5*t2[1]; 
  fmodal[2] = 0.5*t2[2]; 
  fmodal[3] = 0.5*t2[3]; 
} 
