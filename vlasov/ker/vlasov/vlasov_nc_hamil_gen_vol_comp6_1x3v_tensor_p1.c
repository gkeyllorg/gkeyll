#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_1x3v_tensor_p1(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.4330127018922193*alpha[3]*dx10*jacob_cx_inv; 
  out[1] += 0.4330127018922193*alpha[6]*dx10*jacob_cx_inv; 
  out[2] += 0.4330127018922193*alpha[7]*dx10*jacob_cx_inv; 
  out[3] += (0.38729833462074165*alpha[12]+0.4330127018922193*alpha[0])*dx10*jacob_cx_inv; 
  out[4] += 0.4330127018922193*alpha[10]*dx10*jacob_cx_inv; 
  out[5] += 0.4330127018922193*alpha[14]*dx10*jacob_cx_inv; 
  out[6] += (0.3872983346207417*alpha[20]+0.4330127018922193*alpha[1])*dx10*jacob_cx_inv; 
  out[7] += (0.3872983346207417*alpha[21]+0.4330127018922193*alpha[2])*dx10*jacob_cx_inv; 
  out[8] += 0.4330127018922193*alpha[16]*dx10*jacob_cx_inv; 
  out[9] += 0.4330127018922193*alpha[17]*dx10*jacob_cx_inv; 
  out[10] += (0.3872983346207417*alpha[23]+0.4330127018922193*alpha[4])*dx10*jacob_cx_inv; 
  out[11] += 0.43301270189221935*alpha[19]*dx10*jacob_cx_inv; 
  out[12] += 0.38729833462074165*alpha[3]*dx10*jacob_cx_inv; 
  out[13] += 0.43301270189221935*alpha[26]*dx10*jacob_cx_inv; 
  out[14] += (0.38729833462074165*alpha[29]+0.4330127018922193*alpha[5])*dx10*jacob_cx_inv; 
  out[15] += 0.4330127018922193*alpha[27]*dx10*jacob_cx_inv; 
  out[16] += (0.38729833462074165*alpha[32]+0.4330127018922193*alpha[8])*dx10*jacob_cx_inv; 
  out[17] += (0.38729833462074165*alpha[33]+0.4330127018922193*alpha[9])*dx10*jacob_cx_inv; 
  out[18] += 0.43301270189221935*alpha[28]*dx10*jacob_cx_inv; 
  out[19] += (0.3872983346207417*alpha[37]+0.43301270189221935*alpha[11])*dx10*jacob_cx_inv; 
  out[20] += 0.3872983346207417*alpha[6]*dx10*jacob_cx_inv; 
  out[21] += 0.3872983346207417*alpha[7]*dx10*jacob_cx_inv; 
  out[22] += 0.43301270189221935*alpha[31]*dx10*jacob_cx_inv; 
  out[23] += 0.3872983346207417*alpha[10]*dx10*jacob_cx_inv; 
  out[24] += 0.43301270189221935*alpha[35]*dx10*jacob_cx_inv; 
  out[25] += 0.43301270189221935*alpha[36]*dx10*jacob_cx_inv; 
  out[26] += (0.3872983346207417*alpha[39]+0.43301270189221935*alpha[13])*dx10*jacob_cx_inv; 
  out[27] += (0.3872983346207417*alpha[41]+0.4330127018922193*alpha[15])*dx10*jacob_cx_inv; 
  out[28] += (0.38729833462074165*alpha[43]+0.43301270189221935*alpha[18])*dx10*jacob_cx_inv; 
  out[29] += 0.38729833462074165*alpha[14]*dx10*jacob_cx_inv; 
  out[30] += 0.43301270189221935*alpha[40]*dx10*jacob_cx_inv; 
  out[31] += (0.38729833462074165*alpha[44]+0.43301270189221935*alpha[22])*dx10*jacob_cx_inv; 
  out[32] += 0.38729833462074165*alpha[16]*dx10*jacob_cx_inv; 
  out[33] += 0.38729833462074165*alpha[17]*dx10*jacob_cx_inv; 
  out[34] += 0.43301270189221935*alpha[42]*dx10*jacob_cx_inv; 
  out[35] += (0.38729833462074165*alpha[47]+0.43301270189221935*alpha[24])*dx10*jacob_cx_inv; 
  out[36] += (0.38729833462074165*alpha[48]+0.43301270189221935*alpha[25])*dx10*jacob_cx_inv; 
  out[37] += 0.3872983346207417*alpha[19]*dx10*jacob_cx_inv; 
  out[38] += 0.4330127018922193*alpha[46]*dx10*jacob_cx_inv; 
  out[39] += 0.3872983346207417*alpha[26]*dx10*jacob_cx_inv; 
  out[40] += (0.3872983346207417*alpha[49]+0.43301270189221935*alpha[30])*dx10*jacob_cx_inv; 
  out[41] += 0.3872983346207417*alpha[27]*dx10*jacob_cx_inv; 
  out[42] += (0.3872983346207417*alpha[51]+0.43301270189221935*alpha[34])*dx10*jacob_cx_inv; 
  out[43] += 0.38729833462074165*alpha[28]*dx10*jacob_cx_inv; 
  out[44] += 0.38729833462074165*alpha[31]*dx10*jacob_cx_inv; 
  out[45] += 0.4330127018922193*alpha[50]*dx10*jacob_cx_inv; 
  out[46] += (0.38729833462074165*alpha[52]+0.4330127018922193*alpha[38])*dx10*jacob_cx_inv; 
  out[47] += 0.38729833462074165*alpha[35]*dx10*jacob_cx_inv; 
  out[48] += 0.38729833462074165*alpha[36]*dx10*jacob_cx_inv; 
  out[49] += 0.3872983346207417*alpha[40]*dx10*jacob_cx_inv; 
  out[50] += (0.3872983346207417*alpha[53]+0.4330127018922193*alpha[45])*dx10*jacob_cx_inv; 
  out[51] += 0.3872983346207417*alpha[42]*dx10*jacob_cx_inv; 
  out[52] += 0.38729833462074165*alpha[46]*dx10*jacob_cx_inv; 
  out[53] += 0.3872983346207417*alpha[50]*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  out[0] += 0.4330127018922193*alpha[1]*dv11; 
  out[1] += 0.4330127018922193*alpha[0]*dv11; 
  out[2] += 0.4330127018922193*alpha[5]*dv11; 
  out[3] += 0.4330127018922193*alpha[6]*dv11; 
  out[4] += 0.4330127018922193*alpha[8]*dv11; 
  out[5] += 0.4330127018922193*alpha[2]*dv11; 
  out[6] += 0.4330127018922193*alpha[3]*dv11; 
  out[7] += 0.4330127018922193*alpha[14]*dv11; 
  out[8] += 0.4330127018922193*alpha[4]*dv11; 
  out[9] += 0.4330127018922193*alpha[15]*dv11; 
  out[10] += 0.4330127018922193*alpha[16]*dv11; 
  out[11] += 0.43301270189221935*alpha[18]*dv11; 
  out[12] += 0.43301270189221935*alpha[20]*dv11; 
  out[13] += 0.43301270189221935*alpha[24]*dv11; 
  out[14] += 0.4330127018922193*alpha[7]*dv11; 
  out[15] += 0.4330127018922193*alpha[9]*dv11; 
  out[16] += 0.4330127018922193*alpha[10]*dv11; 
  out[17] += 0.4330127018922193*alpha[27]*dv11; 
  out[18] += 0.43301270189221935*alpha[11]*dv11; 
  out[19] += 0.43301270189221935*alpha[28]*dv11; 
  out[20] += 0.43301270189221935*alpha[12]*dv11; 
  out[21] += 0.43301270189221935*alpha[29]*dv11; 
  out[22] += 0.43301270189221935*alpha[30]*dv11; 
  out[23] += 0.43301270189221935*alpha[32]*dv11; 
  out[24] += 0.43301270189221935*alpha[13]*dv11; 
  out[25] += 0.43301270189221935*alpha[34]*dv11; 
  out[26] += 0.43301270189221935*alpha[35]*dv11; 
  out[27] += 0.4330127018922193*alpha[17]*dv11; 
  out[28] += 0.43301270189221935*alpha[19]*dv11; 
  out[29] += 0.43301270189221935*alpha[21]*dv11; 
  out[30] += 0.43301270189221935*alpha[22]*dv11; 
  out[31] += 0.43301270189221935*alpha[40]*dv11; 
  out[32] += 0.43301270189221935*alpha[23]*dv11; 
  out[33] += 0.43301270189221935*alpha[41]*dv11; 
  out[34] += 0.43301270189221935*alpha[25]*dv11; 
  out[35] += 0.43301270189221935*alpha[26]*dv11; 
  out[36] += 0.43301270189221935*alpha[42]*dv11; 
  out[37] += 0.4330127018922193*alpha[43]*dv11; 
  out[38] += 0.4330127018922193*alpha[45]*dv11; 
  out[39] += 0.4330127018922193*alpha[47]*dv11; 
  out[40] += 0.43301270189221935*alpha[31]*dv11; 
  out[41] += 0.43301270189221935*alpha[33]*dv11; 
  out[42] += 0.43301270189221935*alpha[36]*dv11; 
  out[43] += 0.4330127018922193*alpha[37]*dv11; 
  out[44] += 0.4330127018922193*alpha[49]*dv11; 
  out[45] += 0.4330127018922193*alpha[38]*dv11; 
  out[46] += 0.4330127018922193*alpha[50]*dv11; 
  out[47] += 0.4330127018922193*alpha[39]*dv11; 
  out[48] += 0.4330127018922193*alpha[51]*dv11; 
  out[49] += 0.4330127018922193*alpha[44]*dv11; 
  out[50] += 0.4330127018922193*alpha[46]*dv11; 
  out[51] += 0.4330127018922193*alpha[48]*dv11; 
  out[52] += 0.43301270189221935*alpha[53]*dv11; 
  out[53] += 0.43301270189221935*alpha[52]*dv11; 
  } 

  if (dir == 3) { 
  } 

} 
