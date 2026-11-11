#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_conf_to_modal_3x_ser_p1(const double *fnodal, double* GKYL_RESTRICT fmodal) 
{ 
  // fnodal: values at the 8 interior Gauss-Legendre nodes. 
  // fmodal: the p=1 configuration-space expansion with those values. 
  const double t0[8] = {fnodal[0], fnodal[1], fnodal[2], fnodal[3], fnodal[4], fnodal[5], fnodal[6], fnodal[7]};
  const double t1[8] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4], t0[6]+t0[7], t0[7]-t0[6]};
  const double t2[8] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1], t1[4]+t1[6], t1[5]+t1[7], t1[6]-t1[4], t1[7]-t1[5]};
  const double t3[8] = {t2[0]+t2[4], t2[1]+t2[5], t2[2]+t2[6], t2[3]+t2[7], t2[4]-t2[0], t2[5]-t2[1], t2[6]-t2[2], t2[7]-t2[3]};
  fmodal[0] = 0.3535533905932737*t3[0]; 
  fmodal[1] = 0.3535533905932737*t3[1]; 
  fmodal[2] = 0.3535533905932737*t3[2]; 
  fmodal[4] = 0.3535533905932737*t3[3]; 
  fmodal[3] = 0.3535533905932737*t3[4]; 
  fmodal[5] = 0.3535533905932737*t3[5]; 
  fmodal[6] = 0.3535533905932737*t3[6]; 
  fmodal[7] = 0.3535533905932737*t3[7]; 
} 
