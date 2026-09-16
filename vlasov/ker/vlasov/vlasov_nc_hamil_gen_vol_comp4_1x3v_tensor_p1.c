#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_1x3v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double *jacob_vz = &jacob_vel[6]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  out[0] += 0.4330127018922193*alpha[0]*dv12; 
  out[1] += 0.4330127018922193*alpha[1]*dv12; 
  out[2] += 0.4330127018922193*alpha[2]*dv12; 
  out[3] += 0.4330127018922193*alpha[3]*dv12; 
  out[4] += 0.4330127018922193*alpha[4]*dv12; 
  out[5] += 0.4330127018922193*alpha[5]*dv12; 
  out[6] += 0.4330127018922193*alpha[6]*dv12; 
  out[7] += 0.4330127018922193*alpha[7]*dv12; 
  out[8] += 0.4330127018922193*alpha[8]*dv12; 
  out[9] += 0.4330127018922193*alpha[9]*dv12; 
  out[10] += 0.4330127018922193*alpha[10]*dv12; 
  out[11] += 0.4330127018922193*alpha[11]*dv12; 
  out[12] += 0.4330127018922193*alpha[12]*dv12; 
  out[13] += 0.4330127018922193*alpha[13]*dv12; 
  out[14] += 0.4330127018922193*alpha[14]*dv12; 
  out[15] += 0.4330127018922193*alpha[15]*dv12; 
  out[16] += 0.4330127018922193*alpha[16]*dv12; 
  out[17] += 0.4330127018922193*alpha[17]*dv12; 
  out[18] += 0.4330127018922193*alpha[18]*dv12; 
  out[19] += 0.4330127018922193*alpha[19]*dv12; 
  out[20] += 0.4330127018922193*alpha[20]*dv12; 
  out[21] += 0.4330127018922193*alpha[21]*dv12; 
  out[22] += 0.4330127018922193*alpha[22]*dv12; 
  out[23] += 0.4330127018922193*alpha[23]*dv12; 
  out[24] += 0.4330127018922193*alpha[24]*dv12; 
  out[25] += 0.4330127018922193*alpha[25]*dv12; 
  out[26] += 0.4330127018922193*alpha[26]*dv12; 
  out[27] += 0.4330127018922193*alpha[27]*dv12; 
  out[28] += 0.4330127018922193*alpha[28]*dv12; 
  out[29] += 0.4330127018922193*alpha[29]*dv12; 
  out[30] += 0.4330127018922193*alpha[30]*dv12; 
  out[31] += 0.4330127018922193*alpha[31]*dv12; 
  out[32] += 0.4330127018922193*alpha[32]*dv12; 
  out[33] += 0.4330127018922193*alpha[33]*dv12; 
  out[34] += 0.4330127018922193*alpha[34]*dv12; 
  out[35] += 0.4330127018922193*alpha[35]*dv12; 
  out[36] += 0.4330127018922193*alpha[36]*dv12; 
  out[37] += 0.4330127018922193*alpha[37]*dv12; 
  out[38] += 0.4330127018922193*alpha[38]*dv12; 
  out[39] += 0.4330127018922193*alpha[39]*dv12; 
  out[40] += 0.4330127018922193*alpha[40]*dv12; 
  out[41] += 0.4330127018922193*alpha[41]*dv12; 
  out[42] += 0.4330127018922193*alpha[42]*dv12; 
  out[43] += 0.4330127018922193*alpha[43]*dv12; 
  out[44] += 0.4330127018922193*alpha[44]*dv12; 
  out[45] += 0.4330127018922193*alpha[45]*dv12; 
  out[46] += 0.4330127018922193*alpha[46]*dv12; 
  out[47] += 0.4330127018922193*alpha[47]*dv12; 
  out[48] += 0.4330127018922193*alpha[48]*dv12; 
  out[49] += 0.4330127018922193*alpha[49]*dv12; 
  out[50] += 0.4330127018922193*alpha[50]*dv12; 
  out[51] += 0.4330127018922193*alpha[51]*dv12; 
  out[52] += 0.4330127018922193*alpha[52]*dv12; 
  out[53] += 0.4330127018922193*alpha[53]*dv12; 
  } 

} 
