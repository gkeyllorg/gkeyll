#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_phase_to_nodal_2x1v_hyb_p1(const double *fmodal, double* GKYL_RESTRICT fnodal) 
{ 
  // fmodal: phase-space expansion (12 coefficients). 
  // fnodal: its 3 velocity modes at each of the 4 configuration-space nodes, fnodal[q*3 + j]. 
  const double t0[12] = {fmodal[0], fmodal[1], fmodal[2], fmodal[4], fmodal[3], fmodal[5], fmodal[6], fmodal[8], fmodal[7], fmodal[9], fmodal[10], fmodal[11]};
  const double t1[12] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3], t0[4]-t0[5], t0[4]+t0[5], t0[6]-t0[7], t0[6]+t0[7], t0[8]-t0[9], t0[8]+t0[9], t0[10]-t0[11], t0[10]+t0[11]};
  const double t2[12] = {t1[0]-t1[2], t1[1]-t1[3], t1[0]+t1[2], t1[1]+t1[3], t1[4]-t1[6], t1[5]-t1[7], t1[4]+t1[6], t1[5]+t1[7], t1[8]-t1[10], t1[9]-t1[11], t1[8]+t1[10], t1[9]+t1[11]};
  fnodal[0] = 0.5*t2[0]; 
  fnodal[1] = 0.5*t2[4]; 
  fnodal[2] = 0.5*t2[8]; 
  fnodal[3] = 0.5*t2[1]; 
  fnodal[4] = 0.5*t2[5]; 
  fnodal[5] = 0.5*t2[9]; 
  fnodal[6] = 0.5*t2[2]; 
  fnodal[7] = 0.5*t2[6]; 
  fnodal[8] = 0.5*t2[10]; 
  fnodal[9] = 0.5*t2[3]; 
  fnodal[10] = 0.5*t2[7]; 
  fnodal[11] = 0.5*t2[11]; 
} 
