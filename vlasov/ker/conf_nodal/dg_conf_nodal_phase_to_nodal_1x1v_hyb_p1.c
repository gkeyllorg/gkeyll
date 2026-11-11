#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_phase_to_nodal_1x1v_hyb_p1(const double *fmodal, double* GKYL_RESTRICT fnodal) 
{ 
  // fmodal: phase-space expansion (6 coefficients). 
  // fnodal: its 3 velocity modes at each of the 2 configuration-space nodes, fnodal[q*3 + j]. 
  const double t0[6] = {fmodal[0], fmodal[1], fmodal[2], fmodal[3], fmodal[4], fmodal[5]};
  const double t1[6] = {t0[0]-t0[1], t0[0]+t0[1], t0[2]-t0[3], t0[2]+t0[3], t0[4]-t0[5], t0[4]+t0[5]};
  fnodal[0] = 0.7071067811865475*t1[0]; 
  fnodal[1] = 0.7071067811865475*t1[2]; 
  fnodal[2] = 0.7071067811865475*t1[4]; 
  fnodal[3] = 0.7071067811865475*t1[1]; 
  fnodal[4] = 0.7071067811865475*t1[3]; 
  fnodal[5] = 0.7071067811865475*t1[5]; 
} 
