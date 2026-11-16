#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[3]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  out[0] += 0.43301270189221935*alpha[13]*dx11*jacob_cy_inv; 
  out[1] += 0.4330127018922193*alpha[23]*dx11*jacob_cy_inv; 
  out[2] += 0.4330127018922193*alpha[24]*dx11*jacob_cy_inv; 
  out[3] += 0.3872983346207417*alpha[3]*dx11*jacob_cy_inv; 
  out[4] += 0.4330127018922193*alpha[27]*dx11*jacob_cy_inv; 
  out[5] += 0.43301270189221935*alpha[34]*dx11*jacob_cy_inv; 
  out[6] += 0.3872983346207417*alpha[6]*dx11*jacob_cy_inv; 
  out[7] += 0.3872983346207417*alpha[7]*dx11*jacob_cy_inv; 
  out[8] += 0.43301270189221935*alpha[39]*dx11*jacob_cy_inv; 
  out[9] += 0.43301270189221935*alpha[40]*dx11*jacob_cy_inv; 
  out[10] += 0.3872983346207417*alpha[10]*dx11*jacob_cy_inv; 
  out[13] += (0.27664166758624403*alpha[13]+0.43301270189221935*alpha[0])*dx11*jacob_cy_inv; 
  out[15] += 0.3872983346207417*alpha[15]*dx11*jacob_cy_inv; 
  out[16] += 0.4330127018922193*alpha[46]*dx11*jacob_cy_inv; 
  out[17] += 0.3872983346207417*alpha[17]*dx11*jacob_cy_inv; 
  out[18] += 0.3872983346207417*alpha[18]*dx11*jacob_cy_inv; 
  out[21] += 0.3872983346207417*alpha[21]*dx11*jacob_cy_inv; 
  out[22] += 0.3872983346207417*alpha[22]*dx11*jacob_cy_inv; 
  out[23] += (0.27664166758624403*alpha[23]+0.4330127018922193*alpha[1])*dx11*jacob_cy_inv; 
  out[24] += (0.27664166758624403*alpha[24]+0.4330127018922193*alpha[2])*dx11*jacob_cy_inv; 
  out[27] += (0.27664166758624403*alpha[27]+0.4330127018922193*alpha[4])*dx11*jacob_cy_inv; 
  out[30] += 0.3872983346207417*alpha[30]*dx11*jacob_cy_inv; 
  out[31] += 0.3872983346207417*alpha[31]*dx11*jacob_cy_inv; 
  out[32] += 0.3872983346207417*alpha[32]*dx11*jacob_cy_inv; 
  out[33] += 0.3872983346207417*alpha[33]*dx11*jacob_cy_inv; 
  out[34] += (0.27664166758624403*alpha[34]+0.43301270189221935*alpha[5])*dx11*jacob_cy_inv; 
  out[37] += 0.3872983346207417*alpha[37]*dx11*jacob_cy_inv; 
  out[38] += 0.3872983346207417*alpha[38]*dx11*jacob_cy_inv; 
  out[39] += (0.27664166758624403*alpha[39]+0.43301270189221935*alpha[8])*dx11*jacob_cy_inv; 
  out[40] += (0.27664166758624403*alpha[40]+0.43301270189221935*alpha[9])*dx11*jacob_cy_inv; 
  out[42] += 0.3872983346207417*alpha[42]*dx11*jacob_cy_inv; 
  out[43] += 0.3872983346207417*alpha[43]*dx11*jacob_cy_inv; 
  out[44] += 0.3872983346207417*alpha[44]*dx11*jacob_cy_inv; 
  out[45] += 0.3872983346207417*alpha[45]*dx11*jacob_cy_inv; 
  out[46] += (0.27664166758624403*alpha[46]+0.4330127018922193*alpha[16])*dx11*jacob_cy_inv; 
  out[47] += 0.3872983346207417*alpha[47]*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  out[0] += 0.9682458365518543*alpha[7]*dv10*jacob_vx_inv; 
  out[1] += 0.9682458365518543*alpha[15]*dv10*jacob_vx_inv; 
  out[2] += (0.8660254037844386*alpha[22]+0.9682458365518543*alpha[3])*dv10*jacob_vx_inv; 
  out[3] += (0.8660254037844386*alpha[24]+0.9682458365518543*alpha[2])*dv10*jacob_vx_inv; 
  out[4] += 0.9682458365518543*alpha[18]*dv10*jacob_vx_inv; 
  out[5] += (0.8660254037844387*alpha[33]+0.9682458365518543*alpha[6])*dv10*jacob_vx_inv; 
  out[6] += (0.8660254037844387*alpha[34]+0.9682458365518543*alpha[5])*dv10*jacob_vx_inv; 
  out[7] += (0.8660254037844387*alpha[13]+0.8660254037844387*alpha[12]+0.9682458365518543*alpha[0])*dv10*jacob_vx_inv; 
  out[8] += 0.9682458365518543*alpha[31]*dv10*jacob_vx_inv; 
  out[9] += (0.8660254037844387*alpha[38]+0.9682458365518543*alpha[10])*dv10*jacob_vx_inv; 
  out[10] += (0.8660254037844387*alpha[40]+0.9682458365518543*alpha[9])*dv10*jacob_vx_inv; 
  out[11] += 0.9682458365518543*alpha[32]*dv10*jacob_vx_inv; 
  out[12] += 0.8660254037844387*alpha[7]*dv10*jacob_vx_inv; 
  out[13] += 0.8660254037844387*alpha[7]*dv10*jacob_vx_inv; 
  out[14] += 0.9682458365518543*alpha[43]*dv10*jacob_vx_inv; 
  out[15] += (0.8660254037844386*alpha[23]+0.8660254037844386*alpha[20]+0.9682458365518543*alpha[1])*dv10*jacob_vx_inv; 
  out[16] += (0.8660254037844386*alpha[45]+0.9682458365518543*alpha[17])*dv10*jacob_vx_inv; 
  out[17] += (0.8660254037844386*alpha[46]+0.9682458365518543*alpha[16])*dv10*jacob_vx_inv; 
  out[18] += (0.8660254037844386*alpha[27]+0.8660254037844386*alpha[26]+0.9682458365518543*alpha[4])*dv10*jacob_vx_inv; 
  out[19] += 0.9682458365518543*alpha[21]*dv10*jacob_vx_inv; 
  out[20] += 0.8660254037844386*alpha[15]*dv10*jacob_vx_inv; 
  out[21] += 0.9682458365518543*alpha[19]*dv10*jacob_vx_inv; 
  out[22] += (0.7745966692414834*alpha[24]+0.8660254037844386*alpha[2])*dv10*jacob_vx_inv; 
  out[23] += 0.8660254037844386*alpha[15]*dv10*jacob_vx_inv; 
  out[24] += (0.7745966692414834*alpha[22]+0.8660254037844386*alpha[3])*dv10*jacob_vx_inv; 
  out[25] += 0.9682458365518543*alpha[44]*dv10*jacob_vx_inv; 
  out[26] += 0.8660254037844386*alpha[18]*dv10*jacob_vx_inv; 
  out[27] += 0.8660254037844386*alpha[18]*dv10*jacob_vx_inv; 
  out[28] += 0.9682458365518543*alpha[47]*dv10*jacob_vx_inv; 
  out[29] += 0.9682458365518543*alpha[30]*dv10*jacob_vx_inv; 
  out[30] += 0.9682458365518543*alpha[29]*dv10*jacob_vx_inv; 
  out[31] += (0.8660254037844387*alpha[39]+0.8660254037844387*alpha[36]+0.9682458365518543*alpha[8])*dv10*jacob_vx_inv; 
  out[32] += 0.9682458365518543*alpha[11]*dv10*jacob_vx_inv; 
  out[33] += (0.7745966692414834*alpha[34]+0.8660254037844387*alpha[5])*dv10*jacob_vx_inv; 
  out[34] += (0.7745966692414834*alpha[33]+0.8660254037844387*alpha[6])*dv10*jacob_vx_inv; 
  out[35] += 0.9682458365518543*alpha[37]*dv10*jacob_vx_inv; 
  out[36] += 0.8660254037844387*alpha[31]*dv10*jacob_vx_inv; 
  out[37] += 0.9682458365518543*alpha[35]*dv10*jacob_vx_inv; 
  out[38] += (0.7745966692414834*alpha[40]+0.8660254037844387*alpha[9])*dv10*jacob_vx_inv; 
  out[39] += 0.8660254037844387*alpha[31]*dv10*jacob_vx_inv; 
  out[40] += (0.7745966692414834*alpha[38]+0.8660254037844387*alpha[10])*dv10*jacob_vx_inv; 
  out[41] += 0.9682458365518543*alpha[42]*dv10*jacob_vx_inv; 
  out[42] += 0.9682458365518543*alpha[41]*dv10*jacob_vx_inv; 
  out[43] += 0.9682458365518543*alpha[14]*dv10*jacob_vx_inv; 
  out[44] += 0.9682458365518543*alpha[25]*dv10*jacob_vx_inv; 
  out[45] += (0.7745966692414834*alpha[46]+0.8660254037844386*alpha[16])*dv10*jacob_vx_inv; 
  out[46] += (0.7745966692414834*alpha[45]+0.8660254037844386*alpha[17])*dv10*jacob_vx_inv; 
  out[47] += 0.9682458365518543*alpha[28]*dv10*jacob_vx_inv; 
  } 

  if (dir == 3) { 
  } 

} 
