#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel[6]; 
  
  if (dir == 0) { 
  out[0] += 0.4330127018922193*alpha[3]*dx10*jacob_cx_inv; 
  out[1] += 0.4330127018922193*alpha[6]*dx10*jacob_cx_inv; 
  out[2] += 0.4330127018922193*alpha[7]*dx10*jacob_cx_inv; 
  out[3] += (0.38729833462074165*alpha[13]+0.4330127018922193*alpha[0])*dx10*jacob_cx_inv; 
  out[4] += 0.4330127018922193*alpha[10]*dx10*jacob_cx_inv; 
  out[5] += 0.4330127018922193*alpha[15]*dx10*jacob_cx_inv; 
  out[6] += (0.3872983346207417*alpha[23]+0.4330127018922193*alpha[1])*dx10*jacob_cx_inv; 
  out[7] += (0.3872983346207417*alpha[24]+0.4330127018922193*alpha[2])*dx10*jacob_cx_inv; 
  out[8] += 0.4330127018922193*alpha[17]*dx10*jacob_cx_inv; 
  out[9] += 0.4330127018922193*alpha[18]*dx10*jacob_cx_inv; 
  out[10] += (0.3872983346207417*alpha[27]+0.4330127018922193*alpha[4])*dx10*jacob_cx_inv; 
  out[11] += 0.43301270189221935*alpha[21]*dx10*jacob_cx_inv; 
  out[12] += 0.43301270189221935*alpha[22]*dx10*jacob_cx_inv; 
  out[13] += 0.38729833462074165*alpha[3]*dx10*jacob_cx_inv; 
  out[14] += 0.43301270189221935*alpha[30]*dx10*jacob_cx_inv; 
  out[15] += (0.38729833462074165*alpha[34]+0.4330127018922193*alpha[5])*dx10*jacob_cx_inv; 
  out[16] += 0.4330127018922193*alpha[31]*dx10*jacob_cx_inv; 
  out[17] += (0.38729833462074165*alpha[39]+0.4330127018922193*alpha[8])*dx10*jacob_cx_inv; 
  out[18] += (0.38729833462074165*alpha[40]+0.4330127018922193*alpha[9])*dx10*jacob_cx_inv; 
  out[19] += 0.43301270189221935*alpha[32]*dx10*jacob_cx_inv; 
  out[20] += 0.43301270189221935*alpha[33]*dx10*jacob_cx_inv; 
  out[21] += 0.43301270189221935*alpha[11]*dx10*jacob_cx_inv; 
  out[22] += 0.43301270189221935*alpha[12]*dx10*jacob_cx_inv; 
  out[23] += 0.3872983346207417*alpha[6]*dx10*jacob_cx_inv; 
  out[24] += 0.3872983346207417*alpha[7]*dx10*jacob_cx_inv; 
  out[25] += 0.43301270189221935*alpha[37]*dx10*jacob_cx_inv; 
  out[26] += 0.43301270189221935*alpha[38]*dx10*jacob_cx_inv; 
  out[27] += 0.3872983346207417*alpha[10]*dx10*jacob_cx_inv; 
  out[28] += 0.43301270189221935*alpha[42]*dx10*jacob_cx_inv; 
  out[29] += 0.43301270189221935*alpha[43]*dx10*jacob_cx_inv; 
  out[30] += 0.43301270189221935*alpha[14]*dx10*jacob_cx_inv; 
  out[31] += (0.3872983346207417*alpha[46]+0.4330127018922193*alpha[16])*dx10*jacob_cx_inv; 
  out[32] += 0.43301270189221935*alpha[19]*dx10*jacob_cx_inv; 
  out[33] += 0.43301270189221935*alpha[20]*dx10*jacob_cx_inv; 
  out[34] += 0.38729833462074165*alpha[15]*dx10*jacob_cx_inv; 
  out[35] += 0.43301270189221935*alpha[44]*dx10*jacob_cx_inv; 
  out[36] += 0.43301270189221935*alpha[45]*dx10*jacob_cx_inv; 
  out[37] += 0.43301270189221935*alpha[25]*dx10*jacob_cx_inv; 
  out[38] += 0.43301270189221935*alpha[26]*dx10*jacob_cx_inv; 
  out[39] += 0.38729833462074165*alpha[17]*dx10*jacob_cx_inv; 
  out[40] += 0.38729833462074165*alpha[18]*dx10*jacob_cx_inv; 
  out[41] += 0.43301270189221935*alpha[47]*dx10*jacob_cx_inv; 
  out[42] += 0.43301270189221935*alpha[28]*dx10*jacob_cx_inv; 
  out[43] += 0.43301270189221935*alpha[29]*dx10*jacob_cx_inv; 
  out[44] += 0.43301270189221935*alpha[35]*dx10*jacob_cx_inv; 
  out[45] += 0.43301270189221935*alpha[36]*dx10*jacob_cx_inv; 
  out[46] += 0.3872983346207417*alpha[31]*dx10*jacob_cx_inv; 
  out[47] += 0.43301270189221935*alpha[41]*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  out[0] += 0.4330127018922193*alpha[1]*dv11*jacob_vy_inv; 
  out[1] += (0.38729833462074165*alpha[11]+0.4330127018922193*alpha[0])*dv11*jacob_vy_inv; 
  out[2] += 0.4330127018922193*alpha[5]*dv11*jacob_vy_inv; 
  out[3] += 0.4330127018922193*alpha[6]*dv11*jacob_vy_inv; 
  out[4] += 0.4330127018922193*alpha[8]*dv11*jacob_vy_inv; 
  out[5] += (0.3872983346207417*alpha[19]+0.4330127018922193*alpha[2])*dv11*jacob_vy_inv; 
  out[6] += (0.3872983346207417*alpha[21]+0.4330127018922193*alpha[3])*dv11*jacob_vy_inv; 
  out[7] += 0.4330127018922193*alpha[15]*dv11*jacob_vy_inv; 
  out[8] += (0.3872983346207417*alpha[25]+0.4330127018922193*alpha[4])*dv11*jacob_vy_inv; 
  out[9] += 0.4330127018922193*alpha[16]*dv11*jacob_vy_inv; 
  out[10] += 0.4330127018922193*alpha[17]*dv11*jacob_vy_inv; 
  out[11] += 0.38729833462074165*alpha[1]*dv11*jacob_vy_inv; 
  out[12] += 0.43301270189221935*alpha[20]*dv11*jacob_vy_inv; 
  out[13] += 0.43301270189221935*alpha[23]*dv11*jacob_vy_inv; 
  out[14] += 0.43301270189221935*alpha[28]*dv11*jacob_vy_inv; 
  out[15] += (0.38729833462074165*alpha[32]+0.4330127018922193*alpha[7])*dv11*jacob_vy_inv; 
  out[16] += (0.38729833462074165*alpha[35]+0.4330127018922193*alpha[9])*dv11*jacob_vy_inv; 
  out[17] += (0.38729833462074165*alpha[37]+0.4330127018922193*alpha[10])*dv11*jacob_vy_inv; 
  out[18] += 0.4330127018922193*alpha[31]*dv11*jacob_vy_inv; 
  out[19] += 0.3872983346207417*alpha[5]*dv11*jacob_vy_inv; 
  out[20] += 0.43301270189221935*alpha[12]*dv11*jacob_vy_inv; 
  out[21] += 0.3872983346207417*alpha[6]*dv11*jacob_vy_inv; 
  out[22] += 0.43301270189221935*alpha[33]*dv11*jacob_vy_inv; 
  out[23] += 0.43301270189221935*alpha[13]*dv11*jacob_vy_inv; 
  out[24] += 0.43301270189221935*alpha[34]*dv11*jacob_vy_inv; 
  out[25] += 0.3872983346207417*alpha[8]*dv11*jacob_vy_inv; 
  out[26] += 0.43301270189221935*alpha[36]*dv11*jacob_vy_inv; 
  out[27] += 0.43301270189221935*alpha[39]*dv11*jacob_vy_inv; 
  out[28] += 0.43301270189221935*alpha[14]*dv11*jacob_vy_inv; 
  out[29] += 0.43301270189221935*alpha[41]*dv11*jacob_vy_inv; 
  out[30] += 0.43301270189221935*alpha[42]*dv11*jacob_vy_inv; 
  out[31] += (0.3872983346207417*alpha[44]+0.4330127018922193*alpha[18])*dv11*jacob_vy_inv; 
  out[32] += 0.38729833462074165*alpha[15]*dv11*jacob_vy_inv; 
  out[33] += 0.43301270189221935*alpha[22]*dv11*jacob_vy_inv; 
  out[34] += 0.43301270189221935*alpha[24]*dv11*jacob_vy_inv; 
  out[35] += 0.38729833462074165*alpha[16]*dv11*jacob_vy_inv; 
  out[36] += 0.43301270189221935*alpha[26]*dv11*jacob_vy_inv; 
  out[37] += 0.38729833462074165*alpha[17]*dv11*jacob_vy_inv; 
  out[38] += 0.43301270189221935*alpha[45]*dv11*jacob_vy_inv; 
  out[39] += 0.43301270189221935*alpha[27]*dv11*jacob_vy_inv; 
  out[40] += 0.43301270189221935*alpha[46]*dv11*jacob_vy_inv; 
  out[41] += 0.43301270189221935*alpha[29]*dv11*jacob_vy_inv; 
  out[42] += 0.43301270189221935*alpha[30]*dv11*jacob_vy_inv; 
  out[43] += 0.43301270189221935*alpha[47]*dv11*jacob_vy_inv; 
  out[44] += 0.3872983346207417*alpha[31]*dv11*jacob_vy_inv; 
  out[45] += 0.43301270189221935*alpha[38]*dv11*jacob_vy_inv; 
  out[46] += 0.43301270189221935*alpha[40]*dv11*jacob_vy_inv; 
  out[47] += 0.43301270189221935*alpha[43]*dv11*jacob_vy_inv; 
  } 

  if (dir == 3) { 
  } 

} 
