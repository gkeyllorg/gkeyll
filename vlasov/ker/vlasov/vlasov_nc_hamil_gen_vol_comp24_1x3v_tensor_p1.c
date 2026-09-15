#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_1x3v_tensor_p1(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.43301270189221935*alpha[13]*dx10*jacob_cx_inv; 
  out[1] += 0.4330127018922193*alpha[24]*dx10*jacob_cx_inv; 
  out[2] += 0.4330127018922193*alpha[25]*dx10*jacob_cx_inv; 
  out[3] += 0.4330127018922193*alpha[26]*dx10*jacob_cx_inv; 
  out[4] += 0.3872983346207417*alpha[4]*dx10*jacob_cx_inv; 
  out[5] += 0.43301270189221935*alpha[34]*dx10*jacob_cx_inv; 
  out[6] += 0.43301270189221935*alpha[35]*dx10*jacob_cx_inv; 
  out[7] += 0.43301270189221935*alpha[36]*dx10*jacob_cx_inv; 
  out[8] += 0.3872983346207417*alpha[8]*dx10*jacob_cx_inv; 
  out[9] += 0.3872983346207417*alpha[9]*dx10*jacob_cx_inv; 
  out[10] += 0.3872983346207417*alpha[10]*dx10*jacob_cx_inv; 
  out[11] += 0.43301270189221935*alpha[38]*dx10*jacob_cx_inv; 
  out[12] += 0.43301270189221935*alpha[39]*dx10*jacob_cx_inv; 
  out[13] += (0.27664166758624403*alpha[13]+0.43301270189221935*alpha[0])*dx10*jacob_cx_inv; 
  out[14] += 0.4330127018922193*alpha[42]*dx10*jacob_cx_inv; 
  out[15] += 0.3872983346207417*alpha[15]*dx10*jacob_cx_inv; 
  out[16] += 0.3872983346207417*alpha[16]*dx10*jacob_cx_inv; 
  out[17] += 0.3872983346207417*alpha[17]*dx10*jacob_cx_inv; 
  out[18] += 0.4330127018922193*alpha[45]*dx10*jacob_cx_inv; 
  out[19] += 0.4330127018922193*alpha[46]*dx10*jacob_cx_inv; 
  out[20] += 0.4330127018922193*alpha[47]*dx10*jacob_cx_inv; 
  out[21] += 0.4330127018922193*alpha[48]*dx10*jacob_cx_inv; 
  out[22] += 0.3872983346207417*alpha[22]*dx10*jacob_cx_inv; 
  out[23] += 0.3872983346207417*alpha[23]*dx10*jacob_cx_inv; 
  out[24] += (0.27664166758624403*alpha[24]+0.4330127018922193*alpha[1])*dx10*jacob_cx_inv; 
  out[25] += (0.27664166758624403*alpha[25]+0.4330127018922193*alpha[2])*dx10*jacob_cx_inv; 
  out[26] += (0.27664166758624403*alpha[26]+0.4330127018922193*alpha[3])*dx10*jacob_cx_inv; 
  out[27] += 0.3872983346207417*alpha[27]*dx10*jacob_cx_inv; 
  out[28] += 0.43301270189221935*alpha[50]*dx10*jacob_cx_inv; 
  out[29] += 0.43301270189221935*alpha[51]*dx10*jacob_cx_inv; 
  out[30] += 0.3872983346207417*alpha[30]*dx10*jacob_cx_inv; 
  out[31] += 0.3872983346207417*alpha[31]*dx10*jacob_cx_inv; 
  out[32] += 0.3872983346207417*alpha[32]*dx10*jacob_cx_inv; 
  out[33] += 0.3872983346207417*alpha[33]*dx10*jacob_cx_inv; 
  out[34] += (0.27664166758624403*alpha[34]+0.43301270189221935*alpha[5])*dx10*jacob_cx_inv; 
  out[35] += (0.27664166758624403*alpha[35]+0.43301270189221935*alpha[6])*dx10*jacob_cx_inv; 
  out[36] += (0.27664166758624403*alpha[36]+0.43301270189221935*alpha[7])*dx10*jacob_cx_inv; 
  out[37] += 0.43301270189221935*alpha[52]*dx10*jacob_cx_inv; 
  out[38] += (0.27664166758624403*alpha[38]+0.43301270189221935*alpha[11])*dx10*jacob_cx_inv; 
  out[39] += (0.27664166758624403*alpha[39]+0.43301270189221935*alpha[12])*dx10*jacob_cx_inv; 
  out[40] += 0.3872983346207417*alpha[40]*dx10*jacob_cx_inv; 
  out[41] += 0.3872983346207417*alpha[41]*dx10*jacob_cx_inv; 
  out[42] += (0.27664166758624403*alpha[42]+0.4330127018922193*alpha[14])*dx10*jacob_cx_inv; 
  out[43] += 0.4330127018922193*alpha[53]*dx10*jacob_cx_inv; 
  out[44] += 0.3872983346207417*alpha[44]*dx10*jacob_cx_inv; 
  out[45] += (0.27664166758624403*alpha[45]+0.4330127018922193*alpha[18])*dx10*jacob_cx_inv; 
  out[46] += (0.27664166758624403*alpha[46]+0.4330127018922193*alpha[19])*dx10*jacob_cx_inv; 
  out[47] += (0.27664166758624403*alpha[47]+0.4330127018922193*alpha[20])*dx10*jacob_cx_inv; 
  out[48] += (0.27664166758624403*alpha[48]+0.4330127018922193*alpha[21])*dx10*jacob_cx_inv; 
  out[49] += 0.3872983346207417*alpha[49]*dx10*jacob_cx_inv; 
  out[50] += (0.27664166758624403*alpha[50]+0.43301270189221935*alpha[28])*dx10*jacob_cx_inv; 
  out[51] += (0.27664166758624403*alpha[51]+0.43301270189221935*alpha[29])*dx10*jacob_cx_inv; 
  out[52] += (0.27664166758624403*alpha[52]+0.43301270189221935*alpha[37])*dx10*jacob_cx_inv; 
  out[53] += (0.27664166758624403*alpha[53]+0.4330127018922193*alpha[43])*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  out[0] += 0.9682458365518543*alpha[8]*dv12; 
  out[1] += 0.9682458365518543*alpha[4]*dv12; 
  out[2] += 0.9682458365518543*alpha[15]*dv12; 
  out[3] += 0.9682458365518543*alpha[16]*dv12; 
  out[4] += (0.8660254037844386*alpha[24]+0.9682458365518543*alpha[1])*dv12; 
  out[5] += 0.9682458365518543*alpha[9]*dv12; 
  out[6] += 0.9682458365518543*alpha[10]*dv12; 
  out[7] += 0.9682458365518543*alpha[27]*dv12; 
  out[8] += (0.8660254037844387*alpha[13]+0.9682458365518543*alpha[0])*dv12; 
  out[9] += (0.8660254037844387*alpha[34]+0.9682458365518543*alpha[5])*dv12; 
  out[10] += (0.8660254037844387*alpha[35]+0.9682458365518543*alpha[6])*dv12; 
  out[11] += 0.9682458365518543*alpha[30]*dv12; 
  out[12] += 0.9682458365518543*alpha[32]*dv12; 
  out[13] += 0.8660254037844387*alpha[8]*dv12; 
  out[14] += 0.9682458365518543*alpha[17]*dv12; 
  out[15] += (0.8660254037844386*alpha[25]+0.9682458365518543*alpha[2])*dv12; 
  out[16] += (0.8660254037844386*alpha[26]+0.9682458365518543*alpha[3])*dv12; 
  out[17] += (0.8660254037844386*alpha[42]+0.9682458365518543*alpha[14])*dv12; 
  out[18] += 0.9682458365518543*alpha[22]*dv12; 
  out[19] += 0.9682458365518543*alpha[40]*dv12; 
  out[20] += 0.9682458365518543*alpha[23]*dv12; 
  out[21] += 0.9682458365518543*alpha[41]*dv12; 
  out[22] += (0.8660254037844386*alpha[45]+0.9682458365518543*alpha[18])*dv12; 
  out[23] += (0.8660254037844386*alpha[47]+0.9682458365518543*alpha[20])*dv12; 
  out[24] += 0.8660254037844386*alpha[4]*dv12; 
  out[25] += 0.8660254037844386*alpha[15]*dv12; 
  out[26] += 0.8660254037844386*alpha[16]*dv12; 
  out[27] += (0.8660254037844387*alpha[36]+0.9682458365518543*alpha[7])*dv12; 
  out[28] += 0.9682458365518543*alpha[31]*dv12; 
  out[29] += 0.9682458365518543*alpha[33]*dv12; 
  out[30] += (0.8660254037844387*alpha[38]+0.9682458365518543*alpha[11])*dv12; 
  out[31] += (0.8660254037844387*alpha[50]+0.9682458365518543*alpha[28])*dv12; 
  out[32] += (0.8660254037844387*alpha[39]+0.9682458365518543*alpha[12])*dv12; 
  out[33] += (0.8660254037844387*alpha[51]+0.9682458365518543*alpha[29])*dv12; 
  out[34] += 0.8660254037844387*alpha[9]*dv12; 
  out[35] += 0.8660254037844387*alpha[10]*dv12; 
  out[36] += 0.8660254037844387*alpha[27]*dv12; 
  out[37] += 0.9682458365518543*alpha[49]*dv12; 
  out[38] += 0.8660254037844387*alpha[30]*dv12; 
  out[39] += 0.8660254037844387*alpha[32]*dv12; 
  out[40] += (0.8660254037844386*alpha[46]+0.9682458365518543*alpha[19])*dv12; 
  out[41] += (0.8660254037844386*alpha[48]+0.9682458365518543*alpha[21])*dv12; 
  out[42] += 0.8660254037844386*alpha[17]*dv12; 
  out[43] += 0.9682458365518543*alpha[44]*dv12; 
  out[44] += (0.8660254037844386*alpha[53]+0.9682458365518543*alpha[43])*dv12; 
  out[45] += 0.8660254037844386*alpha[22]*dv12; 
  out[46] += 0.8660254037844386*alpha[40]*dv12; 
  out[47] += 0.8660254037844386*alpha[23]*dv12; 
  out[48] += 0.8660254037844386*alpha[41]*dv12; 
  out[49] += (0.8660254037844387*alpha[52]+0.9682458365518543*alpha[37])*dv12; 
  out[50] += 0.8660254037844387*alpha[31]*dv12; 
  out[51] += 0.8660254037844387*alpha[33]*dv12; 
  out[52] += 0.8660254037844387*alpha[49]*dv12; 
  out[53] += 0.8660254037844386*alpha[44]*dv12; 
  } 

} 
