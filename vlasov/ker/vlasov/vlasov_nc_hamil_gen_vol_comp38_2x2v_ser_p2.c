#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp38_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.9682458365518543*alpha[18]*dx11*jacob_cy_inv; 
  out[1] += 0.9682458365518543*alpha[31]*dx11*jacob_cy_inv; 
  out[2] += (0.8660254037844386*alpha[38]+0.9682458365518543*alpha[10])*dx11*jacob_cy_inv; 
  out[3] += (0.8660254037844386*alpha[40]+0.9682458365518543*alpha[9])*dx11*jacob_cy_inv; 
  out[4] += (0.8660254037844386*alpha[43]+0.9682458365518543*alpha[7])*dx11*jacob_cy_inv; 
  out[5] += (0.8660254037844387*alpha[45]+0.9682458365518543*alpha[17])*dx11*jacob_cy_inv; 
  out[6] += (0.8660254037844387*alpha[46]+0.9682458365518543*alpha[16])*dx11*jacob_cy_inv; 
  out[7] += (0.8660254037844387*alpha[27]+0.8660254037844387*alpha[26]+0.9682458365518543*alpha[4])*dx11*jacob_cy_inv; 
  out[8] += (0.8660254037844387*alpha[47]+0.9682458365518543*alpha[15])*dx11*jacob_cy_inv; 
  out[9] += (0.8660254037844387*alpha[30]+0.8660254037844387*alpha[22]+0.9682458365518543*alpha[3])*dx11*jacob_cy_inv; 
  out[10] += (0.8660254037844387*alpha[29]+0.8660254037844387*alpha[24]+0.9682458365518543*alpha[2])*dx11*jacob_cy_inv; 
  out[11] += 0.9682458365518543*alpha[44]*dx11*jacob_cy_inv; 
  out[12] += 0.8660254037844386*alpha[18]*dx11*jacob_cy_inv; 
  out[13] += 0.8660254037844386*alpha[18]*dx11*jacob_cy_inv; 
  out[14] += 0.8660254037844386*alpha[18]*dx11*jacob_cy_inv; 
  out[15] += (0.8660254037844386*alpha[39]+0.8660254037844386*alpha[36]+0.9682458365518543*alpha[8])*dx11*jacob_cy_inv; 
  out[16] += (0.8660254037844386*alpha[42]+0.8660254037844386*alpha[33]+0.9682458365518543*alpha[6])*dx11*jacob_cy_inv; 
  out[17] += (0.8660254037844386*alpha[41]+0.8660254037844386*alpha[34]+0.9682458365518543*alpha[5])*dx11*jacob_cy_inv; 
  out[18] += (0.8660254037844386*alpha[14]+0.8660254037844386*alpha[13]+0.8660254037844386*alpha[12]+0.9682458365518543*alpha[0])*dx11*jacob_cy_inv; 
  out[19] += 0.9682458365518543*alpha[37]*dx11*jacob_cy_inv; 
  out[20] += 0.8660254037844387*alpha[31]*dx11*jacob_cy_inv; 
  out[21] += 0.9682458365518543*alpha[35]*dx11*jacob_cy_inv; 
  out[22] += (0.7745966692414834*alpha[40]+0.8660254037844387*alpha[9])*dx11*jacob_cy_inv; 
  out[23] += 0.8660254037844387*alpha[31]*dx11*jacob_cy_inv; 
  out[24] += (0.7745966692414834*alpha[38]+0.8660254037844387*alpha[10])*dx11*jacob_cy_inv; 
  out[25] += 0.9682458365518543*alpha[32]*dx11*jacob_cy_inv; 
  out[26] += (0.7745966692414834*alpha[43]+0.8660254037844387*alpha[7])*dx11*jacob_cy_inv; 
  out[27] += (0.7745966692414834*alpha[43]+0.8660254037844387*alpha[7])*dx11*jacob_cy_inv; 
  out[28] += 0.8660254037844387*alpha[31]*dx11*jacob_cy_inv; 
  out[29] += (0.7745966692414834*alpha[38]+0.8660254037844387*alpha[10])*dx11*jacob_cy_inv; 
  out[30] += (0.7745966692414834*alpha[40]+0.8660254037844387*alpha[9])*dx11*jacob_cy_inv; 
  out[31] += (0.8660254037844387*alpha[28]+0.8660254037844387*alpha[23]+0.8660254037844387*alpha[20]+0.9682458365518543*alpha[1])*dx11*jacob_cy_inv; 
  out[32] += 0.9682458365518543*alpha[25]*dx11*jacob_cy_inv; 
  out[33] += (0.7745966692414834*alpha[46]+0.8660254037844386*alpha[16])*dx11*jacob_cy_inv; 
  out[34] += (0.7745966692414834*alpha[45]+0.8660254037844386*alpha[17])*dx11*jacob_cy_inv; 
  out[35] += 0.9682458365518543*alpha[21]*dx11*jacob_cy_inv; 
  out[36] += (0.7745966692414834*alpha[47]+0.8660254037844386*alpha[15])*dx11*jacob_cy_inv; 
  out[37] += 0.9682458365518543*alpha[19]*dx11*jacob_cy_inv; 
  out[38] += (0.7745966692414834*alpha[29]+0.7745966692414834*alpha[24]+0.8660254037844386*alpha[2])*dx11*jacob_cy_inv; 
  out[39] += (0.7745966692414834*alpha[47]+0.8660254037844386*alpha[15])*dx11*jacob_cy_inv; 
  out[40] += (0.7745966692414834*alpha[30]+0.7745966692414834*alpha[22]+0.8660254037844386*alpha[3])*dx11*jacob_cy_inv; 
  out[41] += (0.7745966692414834*alpha[45]+0.8660254037844386*alpha[17])*dx11*jacob_cy_inv; 
  out[42] += (0.7745966692414834*alpha[46]+0.8660254037844386*alpha[16])*dx11*jacob_cy_inv; 
  out[43] += (0.7745966692414834*alpha[27]+0.7745966692414834*alpha[26]+0.8660254037844386*alpha[4])*dx11*jacob_cy_inv; 
  out[44] += 0.9682458365518543*alpha[11]*dx11*jacob_cy_inv; 
  out[45] += (0.7745966692414834*alpha[41]+0.7745966692414834*alpha[34]+0.8660254037844387*alpha[5])*dx11*jacob_cy_inv; 
  out[46] += (0.7745966692414834*alpha[42]+0.7745966692414834*alpha[33]+0.8660254037844387*alpha[6])*dx11*jacob_cy_inv; 
  out[47] += (0.7745966692414834*alpha[39]+0.7745966692414834*alpha[36]+0.8660254037844387*alpha[8])*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  out[0] += 0.43301270189221935*alpha[26]*dv10*jacob_vx_inv; 
  out[1] += 0.4330127018922193*alpha[36]*dv10*jacob_vx_inv; 
  out[2] += 0.38729833462074165*alpha[9]*dv10*jacob_vx_inv; 
  out[3] += 0.4330127018922193*alpha[38]*dv10*jacob_vx_inv; 
  out[4] += 0.4330127018922193*alpha[12]*dv10*jacob_vx_inv; 
  out[5] += 0.38729833462074165*alpha[16]*dv10*jacob_vx_inv; 
  out[6] += 0.43301270189221935*alpha[45]*dv10*jacob_vx_inv; 
  out[7] += 0.38729833462074165*alpha[18]*dv10*jacob_vx_inv; 
  out[8] += 0.43301270189221935*alpha[20]*dv10*jacob_vx_inv; 
  out[9] += (0.3464101615137754*alpha[29]+0.38729833462074165*alpha[2])*dv10*jacob_vx_inv; 
  out[10] += 0.43301270189221935*alpha[22]*dv10*jacob_vx_inv; 
  out[12] += (0.27664166758624403*alpha[26]+0.4330127018922193*alpha[4])*dv10*jacob_vx_inv; 
  out[14] += 0.3872983346207417*alpha[26]*dv10*jacob_vx_inv; 
  out[15] += 0.38729833462074165*alpha[31]*dv10*jacob_vx_inv; 
  out[16] += (0.34641016151377546*alpha[41]+0.38729833462074165*alpha[5])*dv10*jacob_vx_inv; 
  out[17] += 0.4330127018922193*alpha[33]*dv10*jacob_vx_inv; 
  out[18] += (0.34641016151377546*alpha[43]+0.38729833462074165*alpha[7])*dv10*jacob_vx_inv; 
  out[19] += 0.3872983346207417*alpha[35]*dv10*jacob_vx_inv; 
  out[20] += (0.27664166758624403*alpha[36]+0.43301270189221935*alpha[8])*dv10*jacob_vx_inv; 
  out[22] += (0.27664166758624403*alpha[38]+0.43301270189221935*alpha[10])*dv10*jacob_vx_inv; 
  out[24] += 0.3872983346207417*alpha[40]*dv10*jacob_vx_inv; 
  out[26] += (0.3872983346207417*alpha[14]+0.27664166758624403*alpha[12]+0.43301270189221935*alpha[0])*dv10*jacob_vx_inv; 
  out[28] += 0.3872983346207417*alpha[36]*dv10*jacob_vx_inv; 
  out[29] += 0.3464101615137754*alpha[9]*dv10*jacob_vx_inv; 
  out[30] += 0.3872983346207417*alpha[38]*dv10*jacob_vx_inv; 
  out[31] += (0.3464101615137754*alpha[47]+0.38729833462074165*alpha[15])*dv10*jacob_vx_inv; 
  out[32] += 0.3872983346207417*alpha[44]*dv10*jacob_vx_inv; 
  out[33] += (0.27664166758624403*alpha[45]+0.4330127018922193*alpha[17])*dv10*jacob_vx_inv; 
  out[34] += 0.3872983346207417*alpha[46]*dv10*jacob_vx_inv; 
  out[35] += 0.3872983346207417*alpha[19]*dv10*jacob_vx_inv; 
  out[36] += (0.3872983346207417*alpha[28]+0.27664166758624403*alpha[20]+0.4330127018922193*alpha[1])*dv10*jacob_vx_inv; 
  out[38] += (0.3872983346207417*alpha[30]+0.27664166758624403*alpha[22]+0.4330127018922193*alpha[3])*dv10*jacob_vx_inv; 
  out[40] += 0.3872983346207417*alpha[24]*dv10*jacob_vx_inv; 
  out[41] += 0.34641016151377546*alpha[16]*dv10*jacob_vx_inv; 
  out[42] += 0.3872983346207417*alpha[45]*dv10*jacob_vx_inv; 
  out[43] += 0.34641016151377546*alpha[18]*dv10*jacob_vx_inv; 
  out[44] += 0.3872983346207417*alpha[32]*dv10*jacob_vx_inv; 
  out[45] += (0.3872983346207417*alpha[42]+0.27664166758624403*alpha[33]+0.43301270189221935*alpha[6])*dv10*jacob_vx_inv; 
  out[46] += 0.3872983346207417*alpha[34]*dv10*jacob_vx_inv; 
  out[47] += 0.3464101615137754*alpha[31]*dv10*jacob_vx_inv; 
  } 

  if (dir == 3) { 
  out[0] += 0.43301270189221935*alpha[22]*dv11*jacob_vy_inv; 
  out[1] += 0.4330127018922193*alpha[33]*dv11*jacob_vy_inv; 
  out[2] += 0.38729833462074165*alpha[7]*dv11*jacob_vy_inv; 
  out[3] += 0.4330127018922193*alpha[12]*dv11*jacob_vy_inv; 
  out[4] += 0.4330127018922193*alpha[38]*dv11*jacob_vy_inv; 
  out[5] += 0.38729833462074165*alpha[15]*dv11*jacob_vy_inv; 
  out[6] += 0.43301270189221935*alpha[20]*dv11*jacob_vy_inv; 
  out[7] += (0.3464101615137754*alpha[24]+0.38729833462074165*alpha[2])*dv11*jacob_vy_inv; 
  out[8] += 0.43301270189221935*alpha[45]*dv11*jacob_vy_inv; 
  out[9] += 0.38729833462074165*alpha[18]*dv11*jacob_vy_inv; 
  out[10] += 0.43301270189221935*alpha[26]*dv11*jacob_vy_inv; 
  out[12] += (0.27664166758624403*alpha[22]+0.4330127018922193*alpha[3])*dv11*jacob_vy_inv; 
  out[13] += 0.3872983346207417*alpha[22]*dv11*jacob_vy_inv; 
  out[15] += (0.34641016151377546*alpha[34]+0.38729833462074165*alpha[5])*dv11*jacob_vy_inv; 
  out[16] += 0.38729833462074165*alpha[31]*dv11*jacob_vy_inv; 
  out[17] += 0.4330127018922193*alpha[36]*dv11*jacob_vy_inv; 
  out[18] += (0.34641016151377546*alpha[40]+0.38729833462074165*alpha[9])*dv11*jacob_vy_inv; 
  out[19] += 0.3872983346207417*alpha[32]*dv11*jacob_vy_inv; 
  out[20] += (0.27664166758624403*alpha[33]+0.43301270189221935*alpha[6])*dv11*jacob_vy_inv; 
  out[22] += (0.3872983346207417*alpha[13]+0.27664166758624403*alpha[12]+0.43301270189221935*alpha[0])*dv11*jacob_vy_inv; 
  out[23] += 0.3872983346207417*alpha[33]*dv11*jacob_vy_inv; 
  out[24] += 0.3464101615137754*alpha[7]*dv11*jacob_vy_inv; 
  out[26] += (0.27664166758624403*alpha[38]+0.43301270189221935*alpha[10])*dv11*jacob_vy_inv; 
  out[27] += 0.3872983346207417*alpha[38]*dv11*jacob_vy_inv; 
  out[29] += 0.3872983346207417*alpha[43]*dv11*jacob_vy_inv; 
  out[31] += (0.3464101615137754*alpha[46]+0.38729833462074165*alpha[16])*dv11*jacob_vy_inv; 
  out[32] += 0.3872983346207417*alpha[19]*dv11*jacob_vy_inv; 
  out[33] += (0.3872983346207417*alpha[23]+0.27664166758624403*alpha[20]+0.4330127018922193*alpha[1])*dv11*jacob_vy_inv; 
  out[34] += 0.34641016151377546*alpha[15]*dv11*jacob_vy_inv; 
  out[35] += 0.3872983346207417*alpha[44]*dv11*jacob_vy_inv; 
  out[36] += (0.27664166758624403*alpha[45]+0.4330127018922193*alpha[17])*dv11*jacob_vy_inv; 
  out[38] += (0.3872983346207417*alpha[27]+0.27664166758624403*alpha[26]+0.4330127018922193*alpha[4])*dv11*jacob_vy_inv; 
  out[39] += 0.3872983346207417*alpha[45]*dv11*jacob_vy_inv; 
  out[40] += 0.34641016151377546*alpha[18]*dv11*jacob_vy_inv; 
  out[41] += 0.3872983346207417*alpha[47]*dv11*jacob_vy_inv; 
  out[43] += 0.3872983346207417*alpha[29]*dv11*jacob_vy_inv; 
  out[44] += 0.3872983346207417*alpha[35]*dv11*jacob_vy_inv; 
  out[45] += (0.3872983346207417*alpha[39]+0.27664166758624403*alpha[36]+0.43301270189221935*alpha[8])*dv11*jacob_vy_inv; 
  out[46] += 0.3464101615137754*alpha[31]*dv11*jacob_vy_inv; 
  out[47] += 0.3872983346207417*alpha[41]*dv11*jacob_vy_inv; 
  } 

} 
