#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_conf_to_nodal_2x_ser_p1(const double *fmodal, double* GKYL_RESTRICT fnodal) 
{ 
  // fmodal: p=1 configuration-space expansion (4 coefficients). 
  // fnodal: its values at the 4 interior Gauss-Legendre nodes (row-major). 
  const double t0[4] = {fmodal[0], fmodal[1], fmodal[2], fmodal[3]};
  const double t1[4] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3]};
  const double t2[4] = {t1[0]-t1[2], t1[1]-t1[3], t1[0]+t1[2], t1[1]+t1[3]};
  fnodal[0] = 0.5*t2[0]; 
  fnodal[1] = 0.5*t2[2]; 
  fnodal[2] = 0.5*t2[1]; 
  fnodal[3] = 0.5*t2[3]; 
} 
