#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_phase_to_modal_accumulate_2x1v_hyb_p1(const double *fnodal, double* GKYL_RESTRICT fmodal) 
{ 
  // fnodal: 3 velocity modes at each of the 4 configuration-space nodes, fnodal[q*3 + j]. 
  // fmodal: phase-space expansion (12 coefficients) to increment with that of fnodal. 
  const double t0[12] = {fnodal[0], fnodal[6], fnodal[3], fnodal[9], fnodal[1], fnodal[7], fnodal[4], fnodal[10], fnodal[2], fnodal[8], fnodal[5], fnodal[11]};
  const double t1[12] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4], t0[6]+t0[7], t0[7]-t0[6], t0[8]+t0[9], t0[9]-t0[8], t0[10]+t0[11], t0[11]-t0[10]};
  const double t2[12] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1], t1[4]+t1[6], t1[5]+t1[7], t1[6]-t1[4], t1[7]-t1[5], t1[8]+t1[10], t1[9]+t1[11], t1[10]-t1[8], t1[11]-t1[9]};
  fmodal[0] += 0.5*t2[0]; 
  fmodal[1] += 0.5*t2[1]; 
  fmodal[2] += 0.5*t2[2]; 
  fmodal[4] += 0.5*t2[3]; 
  fmodal[3] += 0.5*t2[4]; 
  fmodal[5] += 0.5*t2[5]; 
  fmodal[6] += 0.5*t2[6]; 
  fmodal[8] += 0.5*t2[7]; 
  fmodal[7] += 0.5*t2[8]; 
  fmodal[9] += 0.5*t2[9]; 
  fmodal[10] += 0.5*t2[10]; 
  fmodal[11] += 0.5*t2[11]; 
} 
