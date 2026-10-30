#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.4330127018922193*alpha[9]*dx10*jacob_cx_inv; 
  out[1] += 0.4330127018922193*alpha[16]*dx10*jacob_cx_inv; 
  out[2] += (0.3872983346207417*alpha[26]+0.4330127018922193*alpha[4])*dx10*jacob_cx_inv; 
  out[3] += 0.4330127018922193*alpha[18]*dx10*jacob_cx_inv; 
  out[4] += (0.3872983346207417*alpha[29]+0.4330127018922193*alpha[2])*dx10*jacob_cx_inv; 
  out[5] += (0.38729833462074165*alpha[36]+0.4330127018922193*alpha[8])*dx10*jacob_cx_inv; 
  out[6] += 0.4330127018922193*alpha[31]*dx10*jacob_cx_inv; 
  out[7] += (0.38729833462074165*alpha[38]+0.4330127018922193*alpha[10])*dx10*jacob_cx_inv; 
  out[8] += (0.38729833462074165*alpha[41]+0.4330127018922193*alpha[5])*dx10*jacob_cx_inv; 
  out[9] += (0.38729833462074165*alpha[14]+0.38729833462074165*alpha[12]+0.4330127018922193*alpha[0])*dx10*jacob_cx_inv; 
  out[10] += (0.38729833462074165*alpha[43]+0.4330127018922193*alpha[7])*dx10*jacob_cx_inv; 
  out[11] += 0.4330127018922193*alpha[35]*dx10*jacob_cx_inv; 
  out[12] += 0.38729833462074165*alpha[9]*dx10*jacob_cx_inv; 
  out[13] += 0.4330127018922193*alpha[40]*dx10*jacob_cx_inv; 
  out[14] += 0.38729833462074165*alpha[9]*dx10*jacob_cx_inv; 
  out[15] += (0.3872983346207417*alpha[45]+0.4330127018922193*alpha[17])*dx10*jacob_cx_inv; 
  out[16] += (0.3872983346207417*alpha[28]+0.3872983346207417*alpha[20]+0.4330127018922193*alpha[1])*dx10*jacob_cx_inv; 
  out[17] += (0.3872983346207417*alpha[47]+0.4330127018922193*alpha[15])*dx10*jacob_cx_inv; 
  out[18] += (0.3872983346207417*alpha[30]+0.3872983346207417*alpha[22]+0.4330127018922193*alpha[3])*dx10*jacob_cx_inv; 
  out[19] += 0.4330127018922193*alpha[25]*dx10*jacob_cx_inv; 
  out[20] += 0.3872983346207417*alpha[16]*dx10*jacob_cx_inv; 
  out[21] += 0.4330127018922193*alpha[44]*dx10*jacob_cx_inv; 
  out[22] += 0.3872983346207417*alpha[18]*dx10*jacob_cx_inv; 
  out[23] += 0.4330127018922193*alpha[46]*dx10*jacob_cx_inv; 
  out[24] += 0.4330127018922193*alpha[27]*dx10*jacob_cx_inv; 
  out[25] += 0.4330127018922193*alpha[19]*dx10*jacob_cx_inv; 
  out[26] += (0.34641016151377546*alpha[29]+0.3872983346207417*alpha[2])*dx10*jacob_cx_inv; 
  out[27] += 0.4330127018922193*alpha[24]*dx10*jacob_cx_inv; 
  out[28] += 0.3872983346207417*alpha[16]*dx10*jacob_cx_inv; 
  out[29] += (0.34641016151377546*alpha[26]+0.3872983346207417*alpha[4])*dx10*jacob_cx_inv; 
  out[30] += 0.3872983346207417*alpha[18]*dx10*jacob_cx_inv; 
  out[31] += (0.38729833462074165*alpha[42]+0.38729833462074165*alpha[33]+0.4330127018922193*alpha[6])*dx10*jacob_cx_inv; 
  out[32] += 0.4330127018922193*alpha[37]*dx10*jacob_cx_inv; 
  out[33] += 0.38729833462074165*alpha[31]*dx10*jacob_cx_inv; 
  out[34] += 0.4330127018922193*alpha[39]*dx10*jacob_cx_inv; 
  out[35] += 0.4330127018922193*alpha[11]*dx10*jacob_cx_inv; 
  out[36] += (0.34641016151377546*alpha[41]+0.38729833462074165*alpha[5])*dx10*jacob_cx_inv; 
  out[37] += 0.4330127018922193*alpha[32]*dx10*jacob_cx_inv; 
  out[38] += (0.34641016151377546*alpha[43]+0.38729833462074165*alpha[7])*dx10*jacob_cx_inv; 
  out[39] += 0.4330127018922193*alpha[34]*dx10*jacob_cx_inv; 
  out[40] += 0.4330127018922193*alpha[13]*dx10*jacob_cx_inv; 
  out[41] += (0.34641016151377546*alpha[36]+0.38729833462074165*alpha[8])*dx10*jacob_cx_inv; 
  out[42] += 0.38729833462074165*alpha[31]*dx10*jacob_cx_inv; 
  out[43] += (0.34641016151377546*alpha[38]+0.38729833462074165*alpha[10])*dx10*jacob_cx_inv; 
  out[44] += 0.4330127018922193*alpha[21]*dx10*jacob_cx_inv; 
  out[45] += (0.34641016151377546*alpha[47]+0.3872983346207417*alpha[15])*dx10*jacob_cx_inv; 
  out[46] += 0.4330127018922193*alpha[23]*dx10*jacob_cx_inv; 
  out[47] += (0.34641016151377546*alpha[45]+0.3872983346207417*alpha[17])*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  out[0] += 0.4330127018922193*alpha[8]*dv10*jacob_vx_inv; 
  out[1] += (0.3872983346207417*alpha[25]+0.4330127018922193*alpha[4])*dv10*jacob_vx_inv; 
  out[2] += 0.4330127018922193*alpha[16]*dv10*jacob_vx_inv; 
  out[3] += 0.4330127018922193*alpha[17]*dv10*jacob_vx_inv; 
  out[4] += (0.3872983346207417*alpha[28]+0.4330127018922193*alpha[1])*dv10*jacob_vx_inv; 
  out[5] += (0.38729833462074165*alpha[35]+0.4330127018922193*alpha[9])*dv10*jacob_vx_inv; 
  out[6] += (0.38729833462074165*alpha[37]+0.4330127018922193*alpha[10])*dv10*jacob_vx_inv; 
  out[7] += 0.4330127018922193*alpha[31]*dv10*jacob_vx_inv; 
  out[8] += (0.38729833462074165*alpha[14]+0.38729833462074165*alpha[11]+0.4330127018922193*alpha[0])*dv10*jacob_vx_inv; 
  out[9] += (0.38729833462074165*alpha[41]+0.4330127018922193*alpha[5])*dv10*jacob_vx_inv; 
  out[10] += (0.38729833462074165*alpha[42]+0.4330127018922193*alpha[6])*dv10*jacob_vx_inv; 
  out[11] += 0.38729833462074165*alpha[8]*dv10*jacob_vx_inv; 
  out[12] += 0.4330127018922193*alpha[36]*dv10*jacob_vx_inv; 
  out[13] += 0.4330127018922193*alpha[39]*dv10*jacob_vx_inv; 
  out[14] += 0.38729833462074165*alpha[8]*dv10*jacob_vx_inv; 
  out[15] += (0.3872983346207417*alpha[44]+0.4330127018922193*alpha[18])*dv10*jacob_vx_inv; 
  out[16] += (0.3872983346207417*alpha[29]+0.3872983346207417*alpha[19]+0.4330127018922193*alpha[2])*dv10*jacob_vx_inv; 
  out[17] += (0.3872983346207417*alpha[30]+0.3872983346207417*alpha[21]+0.4330127018922193*alpha[3])*dv10*jacob_vx_inv; 
  out[18] += (0.3872983346207417*alpha[47]+0.4330127018922193*alpha[15])*dv10*jacob_vx_inv; 
  out[19] += 0.3872983346207417*alpha[16]*dv10*jacob_vx_inv; 
  out[20] += 0.4330127018922193*alpha[26]*dv10*jacob_vx_inv; 
  out[21] += 0.3872983346207417*alpha[17]*dv10*jacob_vx_inv; 
  out[22] += 0.4330127018922193*alpha[45]*dv10*jacob_vx_inv; 
  out[23] += 0.4330127018922193*alpha[27]*dv10*jacob_vx_inv; 
  out[24] += 0.4330127018922193*alpha[46]*dv10*jacob_vx_inv; 
  out[25] += (0.34641016151377546*alpha[28]+0.3872983346207417*alpha[1])*dv10*jacob_vx_inv; 
  out[26] += 0.4330127018922193*alpha[20]*dv10*jacob_vx_inv; 
  out[27] += 0.4330127018922193*alpha[23]*dv10*jacob_vx_inv; 
  out[28] += (0.34641016151377546*alpha[25]+0.3872983346207417*alpha[4])*dv10*jacob_vx_inv; 
  out[29] += 0.3872983346207417*alpha[16]*dv10*jacob_vx_inv; 
  out[30] += 0.3872983346207417*alpha[17]*dv10*jacob_vx_inv; 
  out[31] += (0.38729833462074165*alpha[43]+0.38729833462074165*alpha[32]+0.4330127018922193*alpha[7])*dv10*jacob_vx_inv; 
  out[32] += 0.38729833462074165*alpha[31]*dv10*jacob_vx_inv; 
  out[33] += 0.4330127018922193*alpha[38]*dv10*jacob_vx_inv; 
  out[34] += 0.4330127018922193*alpha[40]*dv10*jacob_vx_inv; 
  out[35] += (0.34641016151377546*alpha[41]+0.38729833462074165*alpha[5])*dv10*jacob_vx_inv; 
  out[36] += 0.4330127018922193*alpha[12]*dv10*jacob_vx_inv; 
  out[37] += (0.34641016151377546*alpha[42]+0.38729833462074165*alpha[6])*dv10*jacob_vx_inv; 
  out[38] += 0.4330127018922193*alpha[33]*dv10*jacob_vx_inv; 
  out[39] += 0.4330127018922193*alpha[13]*dv10*jacob_vx_inv; 
  out[40] += 0.4330127018922193*alpha[34]*dv10*jacob_vx_inv; 
  out[41] += (0.34641016151377546*alpha[35]+0.38729833462074165*alpha[9])*dv10*jacob_vx_inv; 
  out[42] += (0.34641016151377546*alpha[37]+0.38729833462074165*alpha[10])*dv10*jacob_vx_inv; 
  out[43] += 0.38729833462074165*alpha[31]*dv10*jacob_vx_inv; 
  out[44] += (0.34641016151377546*alpha[47]+0.3872983346207417*alpha[15])*dv10*jacob_vx_inv; 
  out[45] += 0.4330127018922193*alpha[22]*dv10*jacob_vx_inv; 
  out[46] += 0.4330127018922193*alpha[24]*dv10*jacob_vx_inv; 
  out[47] += (0.34641016151377546*alpha[44]+0.3872983346207417*alpha[18])*dv10*jacob_vx_inv; 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  out[0] += 0.4330127018922193*alpha[5]*dv12*jacob_vz_inv; 
  out[1] += (0.3872983346207417*alpha[19]+0.4330127018922193*alpha[2])*dv12*jacob_vz_inv; 
  out[2] += (0.3872983346207417*alpha[20]+0.4330127018922193*alpha[1])*dv12*jacob_vz_inv; 
  out[3] += 0.4330127018922193*alpha[15]*dv12*jacob_vz_inv; 
  out[4] += 0.4330127018922193*alpha[16]*dv12*jacob_vz_inv; 
  out[5] += (0.38729833462074165*alpha[12]+0.38729833462074165*alpha[11]+0.4330127018922193*alpha[0])*dv12*jacob_vz_inv; 
  out[6] += (0.38729833462074165*alpha[32]+0.4330127018922193*alpha[7])*dv12*jacob_vz_inv; 
  out[7] += (0.38729833462074165*alpha[33]+0.4330127018922193*alpha[6])*dv12*jacob_vz_inv; 
  out[8] += (0.38729833462074165*alpha[35]+0.4330127018922193*alpha[9])*dv12*jacob_vz_inv; 
  out[9] += (0.38729833462074165*alpha[36]+0.4330127018922193*alpha[8])*dv12*jacob_vz_inv; 
  out[10] += 0.4330127018922193*alpha[31]*dv12*jacob_vz_inv; 
  out[11] += 0.38729833462074165*alpha[5]*dv12*jacob_vz_inv; 
  out[12] += 0.38729833462074165*alpha[5]*dv12*jacob_vz_inv; 
  out[13] += 0.4330127018922193*alpha[34]*dv12*jacob_vz_inv; 
  out[14] += 0.4330127018922193*alpha[41]*dv12*jacob_vz_inv; 
  out[15] += (0.3872983346207417*alpha[22]+0.3872983346207417*alpha[21]+0.4330127018922193*alpha[3])*dv12*jacob_vz_inv; 
  out[16] += (0.3872983346207417*alpha[26]+0.3872983346207417*alpha[25]+0.4330127018922193*alpha[4])*dv12*jacob_vz_inv; 
  out[17] += (0.3872983346207417*alpha[44]+0.4330127018922193*alpha[18])*dv12*jacob_vz_inv; 
  out[18] += (0.3872983346207417*alpha[45]+0.4330127018922193*alpha[17])*dv12*jacob_vz_inv; 
  out[19] += (0.34641016151377546*alpha[20]+0.3872983346207417*alpha[1])*dv12*jacob_vz_inv; 
  out[20] += (0.34641016151377546*alpha[19]+0.3872983346207417*alpha[2])*dv12*jacob_vz_inv; 
  out[21] += 0.3872983346207417*alpha[15]*dv12*jacob_vz_inv; 
  out[22] += 0.3872983346207417*alpha[15]*dv12*jacob_vz_inv; 
  out[23] += 0.4330127018922193*alpha[24]*dv12*jacob_vz_inv; 
  out[24] += 0.4330127018922193*alpha[23]*dv12*jacob_vz_inv; 
  out[25] += 0.3872983346207417*alpha[16]*dv12*jacob_vz_inv; 
  out[26] += 0.3872983346207417*alpha[16]*dv12*jacob_vz_inv; 
  out[27] += 0.4330127018922193*alpha[46]*dv12*jacob_vz_inv; 
  out[28] += 0.4330127018922193*alpha[29]*dv12*jacob_vz_inv; 
  out[29] += 0.4330127018922193*alpha[28]*dv12*jacob_vz_inv; 
  out[30] += 0.4330127018922193*alpha[47]*dv12*jacob_vz_inv; 
  out[31] += (0.38729833462074165*alpha[38]+0.38729833462074165*alpha[37]+0.4330127018922193*alpha[10])*dv12*jacob_vz_inv; 
  out[32] += (0.34641016151377546*alpha[33]+0.38729833462074165*alpha[6])*dv12*jacob_vz_inv; 
  out[33] += (0.34641016151377546*alpha[32]+0.38729833462074165*alpha[7])*dv12*jacob_vz_inv; 
  out[34] += 0.4330127018922193*alpha[13]*dv12*jacob_vz_inv; 
  out[35] += (0.34641016151377546*alpha[36]+0.38729833462074165*alpha[8])*dv12*jacob_vz_inv; 
  out[36] += (0.34641016151377546*alpha[35]+0.38729833462074165*alpha[9])*dv12*jacob_vz_inv; 
  out[37] += 0.38729833462074165*alpha[31]*dv12*jacob_vz_inv; 
  out[38] += 0.38729833462074165*alpha[31]*dv12*jacob_vz_inv; 
  out[39] += 0.4330127018922193*alpha[40]*dv12*jacob_vz_inv; 
  out[40] += 0.4330127018922193*alpha[39]*dv12*jacob_vz_inv; 
  out[41] += 0.4330127018922193*alpha[14]*dv12*jacob_vz_inv; 
  out[42] += 0.4330127018922193*alpha[43]*dv12*jacob_vz_inv; 
  out[43] += 0.4330127018922193*alpha[42]*dv12*jacob_vz_inv; 
  out[44] += (0.34641016151377546*alpha[45]+0.3872983346207417*alpha[17])*dv12*jacob_vz_inv; 
  out[45] += (0.34641016151377546*alpha[44]+0.3872983346207417*alpha[18])*dv12*jacob_vz_inv; 
  out[46] += 0.4330127018922193*alpha[27]*dv12*jacob_vz_inv; 
  out[47] += 0.4330127018922193*alpha[30]*dv12*jacob_vz_inv; 
  } 

} 
