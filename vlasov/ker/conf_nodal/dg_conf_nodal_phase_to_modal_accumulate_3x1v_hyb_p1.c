#include <gkyl_dg_conf_nodal_kernels.h> 
GKYL_CU_DH void dg_conf_nodal_phase_to_modal_accumulate_3x1v_hyb_p1(const double *fnodal, double* GKYL_RESTRICT fmodal) 
{ 
  // fnodal: 3 velocity modes at each of the 8 configuration-space nodes, fnodal[q*3 + j]. 
  // fmodal: phase-space expansion (24 coefficients) to increment with that of fnodal. 
  const double t0[24] = {fnodal[0], fnodal[3], fnodal[6], fnodal[9], fnodal[12], fnodal[15], fnodal[18], fnodal[21], fnodal[1], fnodal[4], fnodal[7], fnodal[10], fnodal[13], fnodal[16], fnodal[19], fnodal[22], fnodal[2], fnodal[5], fnodal[8], fnodal[11], fnodal[14], fnodal[17], fnodal[20], fnodal[23]};
  const double t1[24] = {t0[0]+t0[1], t0[1]-t0[0], t0[2]+t0[3], t0[3]-t0[2], t0[4]+t0[5], t0[5]-t0[4], t0[6]+t0[7], t0[7]-t0[6], t0[8]+t0[9], t0[9]-t0[8], t0[10]+t0[11], t0[11]-t0[10], t0[12]+t0[13], t0[13]-t0[12], t0[14]+t0[15], t0[15]-t0[14], t0[16]+t0[17], t0[17]-t0[16], t0[18]+t0[19], t0[19]-t0[18], t0[20]+t0[21], t0[21]-t0[20], t0[22]+t0[23], t0[23]-t0[22]};
  const double t2[24] = {t1[0]+t1[2], t1[1]+t1[3], t1[2]-t1[0], t1[3]-t1[1], t1[4]+t1[6], t1[5]+t1[7], t1[6]-t1[4], t1[7]-t1[5], t1[8]+t1[10], t1[9]+t1[11], t1[10]-t1[8], t1[11]-t1[9], t1[12]+t1[14], t1[13]+t1[15], t1[14]-t1[12], t1[15]-t1[13], t1[16]+t1[18], t1[17]+t1[19], t1[18]-t1[16], t1[19]-t1[17], t1[20]+t1[22], t1[21]+t1[23], t1[22]-t1[20], t1[23]-t1[21]};
  const double t3[24] = {t2[0]+t2[4], t2[1]+t2[5], t2[2]+t2[6], t2[3]+t2[7], t2[4]-t2[0], t2[5]-t2[1], t2[6]-t2[2], t2[7]-t2[3], t2[8]+t2[12], t2[9]+t2[13], t2[10]+t2[14], t2[11]+t2[15], t2[12]-t2[8], t2[13]-t2[9], t2[14]-t2[10], t2[15]-t2[11], t2[16]+t2[20], t2[17]+t2[21], t2[18]+t2[22], t2[19]+t2[23], t2[20]-t2[16], t2[21]-t2[17], t2[22]-t2[18], t2[23]-t2[19]};
  fmodal[0] += 0.3535533905932737*t3[0]; 
  fmodal[1] += 0.3535533905932737*t3[1]; 
  fmodal[2] += 0.3535533905932737*t3[2]; 
  fmodal[5] += 0.3535533905932737*t3[3]; 
  fmodal[3] += 0.3535533905932737*t3[4]; 
  fmodal[6] += 0.3535533905932737*t3[5]; 
  fmodal[7] += 0.3535533905932737*t3[6]; 
  fmodal[12] += 0.3535533905932737*t3[7]; 
  fmodal[4] += 0.3535533905932737*t3[8]; 
  fmodal[8] += 0.3535533905932737*t3[9]; 
  fmodal[9] += 0.3535533905932737*t3[10]; 
  fmodal[13] += 0.3535533905932737*t3[11]; 
  fmodal[10] += 0.3535533905932737*t3[12]; 
  fmodal[14] += 0.3535533905932737*t3[13]; 
  fmodal[15] += 0.3535533905932737*t3[14]; 
  fmodal[19] += 0.3535533905932737*t3[15]; 
  fmodal[11] += 0.3535533905932737*t3[16]; 
  fmodal[16] += 0.3535533905932737*t3[17]; 
  fmodal[17] += 0.3535533905932737*t3[18]; 
  fmodal[20] += 0.3535533905932737*t3[19]; 
  fmodal[18] += 0.3535533905932737*t3[20]; 
  fmodal[21] += 0.3535533905932737*t3[21]; 
  fmodal[22] += 0.3535533905932737*t3[22]; 
  fmodal[23] += 0.3535533905932737*t3[23]; 
} 
