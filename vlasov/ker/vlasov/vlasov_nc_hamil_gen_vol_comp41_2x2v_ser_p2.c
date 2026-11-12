#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp41_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
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
  out[0] += 0.43301270189221935*alpha[29]*dx10*jacob_cx_inv; 
  out[1] += 0.4330127018922193*alpha[41]*dx10*jacob_cx_inv; 
  out[2] += 0.4330127018922193*alpha[14]*dx10*jacob_cx_inv; 
  out[3] += 0.4330127018922193*alpha[43]*dx10*jacob_cx_inv; 
  out[4] += 0.38729833462074165*alpha[9]*dx10*jacob_cx_inv; 
  out[5] += 0.43301270189221935*alpha[28]*dx10*jacob_cx_inv; 
  out[6] += 0.43301270189221935*alpha[47]*dx10*jacob_cx_inv; 
  out[7] += 0.43301270189221935*alpha[30]*dx10*jacob_cx_inv; 
  out[8] += 0.38729833462074165*alpha[16]*dx10*jacob_cx_inv; 
  out[9] += (0.3464101615137754*alpha[26]+0.38729833462074165*alpha[4])*dx10*jacob_cx_inv; 
  out[10] += 0.38729833462074165*alpha[18]*dx10*jacob_cx_inv; 
  out[12] += 0.3872983346207417*alpha[29]*dx10*jacob_cx_inv; 
  out[14] += (0.27664166758624403*alpha[29]+0.4330127018922193*alpha[2])*dx10*jacob_cx_inv; 
  out[15] += 0.4330127018922193*alpha[42]*dx10*jacob_cx_inv; 
  out[16] += (0.34641016151377546*alpha[36]+0.38729833462074165*alpha[8])*dx10*jacob_cx_inv; 
  out[17] += 0.38729833462074165*alpha[31]*dx10*jacob_cx_inv; 
  out[18] += (0.34641016151377546*alpha[38]+0.38729833462074165*alpha[10])*dx10*jacob_cx_inv; 
  out[20] += 0.3872983346207417*alpha[41]*dx10*jacob_cx_inv; 
  out[22] += 0.3872983346207417*alpha[43]*dx10*jacob_cx_inv; 
  out[25] += 0.3872983346207417*alpha[35]*dx10*jacob_cx_inv; 
  out[26] += 0.3464101615137754*alpha[9]*dx10*jacob_cx_inv; 
  out[27] += 0.3872983346207417*alpha[40]*dx10*jacob_cx_inv; 
  out[28] += (0.27664166758624403*alpha[41]+0.43301270189221935*alpha[5])*dx10*jacob_cx_inv; 
  out[29] += (0.27664166758624403*alpha[14]+0.3872983346207417*alpha[12]+0.43301270189221935*alpha[0])*dx10*jacob_cx_inv; 
  out[30] += (0.27664166758624403*alpha[43]+0.43301270189221935*alpha[7])*dx10*jacob_cx_inv; 
  out[31] += (0.3464101615137754*alpha[45]+0.38729833462074165*alpha[17])*dx10*jacob_cx_inv; 
  out[33] += 0.3872983346207417*alpha[47]*dx10*jacob_cx_inv; 
  out[35] += 0.3872983346207417*alpha[25]*dx10*jacob_cx_inv; 
  out[36] += 0.34641016151377546*alpha[16]*dx10*jacob_cx_inv; 
  out[37] += 0.3872983346207417*alpha[44]*dx10*jacob_cx_inv; 
  out[38] += 0.34641016151377546*alpha[18]*dx10*jacob_cx_inv; 
  out[39] += 0.3872983346207417*alpha[46]*dx10*jacob_cx_inv; 
  out[40] += 0.3872983346207417*alpha[27]*dx10*jacob_cx_inv; 
  out[41] += (0.27664166758624403*alpha[28]+0.3872983346207417*alpha[20]+0.4330127018922193*alpha[1])*dx10*jacob_cx_inv; 
  out[42] += (0.27664166758624403*alpha[47]+0.4330127018922193*alpha[15])*dx10*jacob_cx_inv; 
  out[43] += (0.27664166758624403*alpha[30]+0.3872983346207417*alpha[22]+0.4330127018922193*alpha[3])*dx10*jacob_cx_inv; 
  out[44] += 0.3872983346207417*alpha[37]*dx10*jacob_cx_inv; 
  out[45] += 0.3464101615137754*alpha[31]*dx10*jacob_cx_inv; 
  out[46] += 0.3872983346207417*alpha[39]*dx10*jacob_cx_inv; 
  out[47] += (0.27664166758624403*alpha[42]+0.3872983346207417*alpha[33]+0.43301270189221935*alpha[6])*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  out[0] += 0.43301270189221935*alpha[28]*dx11*jacob_cy_inv; 
  out[1] += 0.4330127018922193*alpha[14]*dx11*jacob_cy_inv; 
  out[2] += 0.4330127018922193*alpha[41]*dx11*jacob_cy_inv; 
  out[3] += 0.4330127018922193*alpha[42]*dx11*jacob_cy_inv; 
  out[4] += 0.38729833462074165*alpha[8]*dx11*jacob_cy_inv; 
  out[5] += 0.43301270189221935*alpha[29]*dx11*jacob_cy_inv; 
  out[6] += 0.43301270189221935*alpha[30]*dx11*jacob_cy_inv; 
  out[7] += 0.43301270189221935*alpha[47]*dx11*jacob_cy_inv; 
  out[8] += (0.3464101615137754*alpha[25]+0.38729833462074165*alpha[4])*dx11*jacob_cy_inv; 
  out[9] += 0.38729833462074165*alpha[16]*dx11*jacob_cy_inv; 
  out[10] += 0.38729833462074165*alpha[17]*dx11*jacob_cy_inv; 
  out[11] += 0.3872983346207417*alpha[28]*dx11*jacob_cy_inv; 
  out[14] += (0.27664166758624403*alpha[28]+0.4330127018922193*alpha[1])*dx11*jacob_cy_inv; 
  out[15] += 0.4330127018922193*alpha[43]*dx11*jacob_cy_inv; 
  out[16] += (0.34641016151377546*alpha[35]+0.38729833462074165*alpha[9])*dx11*jacob_cy_inv; 
  out[17] += (0.34641016151377546*alpha[37]+0.38729833462074165*alpha[10])*dx11*jacob_cy_inv; 
  out[18] += 0.38729833462074165*alpha[31]*dx11*jacob_cy_inv; 
  out[19] += 0.3872983346207417*alpha[41]*dx11*jacob_cy_inv; 
  out[21] += 0.3872983346207417*alpha[42]*dx11*jacob_cy_inv; 
  out[25] += 0.3464101615137754*alpha[8]*dx11*jacob_cy_inv; 
  out[26] += 0.3872983346207417*alpha[36]*dx11*jacob_cy_inv; 
  out[27] += 0.3872983346207417*alpha[39]*dx11*jacob_cy_inv; 
  out[28] += (0.27664166758624403*alpha[14]+0.3872983346207417*alpha[11]+0.43301270189221935*alpha[0])*dx11*jacob_cy_inv; 
  out[29] += (0.27664166758624403*alpha[41]+0.43301270189221935*alpha[5])*dx11*jacob_cy_inv; 
  out[30] += (0.27664166758624403*alpha[42]+0.43301270189221935*alpha[6])*dx11*jacob_cy_inv; 
  out[31] += (0.3464101615137754*alpha[44]+0.38729833462074165*alpha[18])*dx11*jacob_cy_inv; 
  out[32] += 0.3872983346207417*alpha[47]*dx11*jacob_cy_inv; 
  out[35] += 0.34641016151377546*alpha[16]*dx11*jacob_cy_inv; 
  out[36] += 0.3872983346207417*alpha[26]*dx11*jacob_cy_inv; 
  out[37] += 0.34641016151377546*alpha[17]*dx11*jacob_cy_inv; 
  out[38] += 0.3872983346207417*alpha[45]*dx11*jacob_cy_inv; 
  out[39] += 0.3872983346207417*alpha[27]*dx11*jacob_cy_inv; 
  out[40] += 0.3872983346207417*alpha[46]*dx11*jacob_cy_inv; 
  out[41] += (0.27664166758624403*alpha[29]+0.3872983346207417*alpha[19]+0.4330127018922193*alpha[2])*dx11*jacob_cy_inv; 
  out[42] += (0.27664166758624403*alpha[30]+0.3872983346207417*alpha[21]+0.4330127018922193*alpha[3])*dx11*jacob_cy_inv; 
  out[43] += (0.27664166758624403*alpha[47]+0.4330127018922193*alpha[15])*dx11*jacob_cy_inv; 
  out[44] += 0.3464101615137754*alpha[31]*dx11*jacob_cy_inv; 
  out[45] += 0.3872983346207417*alpha[38]*dx11*jacob_cy_inv; 
  out[46] += 0.3872983346207417*alpha[40]*dx11*jacob_cy_inv; 
  out[47] += (0.27664166758624403*alpha[43]+0.3872983346207417*alpha[32]+0.43301270189221935*alpha[7])*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  out[0] += 0.9682458365518543*alpha[16]*dv11*jacob_vy_inv; 
  out[1] += (0.8660254037844386*alpha[35]+0.9682458365518543*alpha[9])*dv11*jacob_vy_inv; 
  out[2] += (0.8660254037844386*alpha[36]+0.9682458365518543*alpha[8])*dv11*jacob_vy_inv; 
  out[3] += 0.9682458365518543*alpha[31]*dv11*jacob_vy_inv; 
  out[4] += (0.8660254037844386*alpha[41]+0.9682458365518543*alpha[5])*dv11*jacob_vy_inv; 
  out[5] += (0.8660254037844387*alpha[26]+0.8660254037844387*alpha[25]+0.9682458365518543*alpha[4])*dv11*jacob_vy_inv; 
  out[6] += (0.8660254037844387*alpha[44]+0.9682458365518543*alpha[18])*dv11*jacob_vy_inv; 
  out[7] += (0.8660254037844387*alpha[45]+0.9682458365518543*alpha[17])*dv11*jacob_vy_inv; 
  out[8] += (0.8660254037844387*alpha[29]+0.8660254037844387*alpha[19]+0.9682458365518543*alpha[2])*dv11*jacob_vy_inv; 
  out[9] += (0.8660254037844387*alpha[28]+0.8660254037844387*alpha[20]+0.9682458365518543*alpha[1])*dv11*jacob_vy_inv; 
  out[10] += (0.8660254037844387*alpha[47]+0.9682458365518543*alpha[15])*dv11*jacob_vy_inv; 
  out[11] += 0.8660254037844386*alpha[16]*dv11*jacob_vy_inv; 
  out[12] += 0.8660254037844386*alpha[16]*dv11*jacob_vy_inv; 
  out[13] += 0.9682458365518543*alpha[46]*dv11*jacob_vy_inv; 
  out[14] += 0.8660254037844386*alpha[16]*dv11*jacob_vy_inv; 
  out[15] += (0.8660254037844386*alpha[38]+0.8660254037844386*alpha[37]+0.9682458365518543*alpha[10])*dv11*jacob_vy_inv; 
  out[16] += (0.8660254037844386*alpha[14]+0.8660254037844386*alpha[12]+0.8660254037844386*alpha[11]+0.9682458365518543*alpha[0])*dv11*jacob_vy_inv; 
  out[17] += (0.8660254037844386*alpha[43]+0.8660254037844386*alpha[32]+0.9682458365518543*alpha[7])*dv11*jacob_vy_inv; 
  out[18] += (0.8660254037844386*alpha[42]+0.8660254037844386*alpha[33]+0.9682458365518543*alpha[6])*dv11*jacob_vy_inv; 
  out[19] += (0.7745966692414834*alpha[36]+0.8660254037844387*alpha[8])*dv11*jacob_vy_inv; 
  out[20] += (0.7745966692414834*alpha[35]+0.8660254037844387*alpha[9])*dv11*jacob_vy_inv; 
  out[21] += 0.8660254037844387*alpha[31]*dv11*jacob_vy_inv; 
  out[22] += 0.8660254037844387*alpha[31]*dv11*jacob_vy_inv; 
  out[23] += 0.9682458365518543*alpha[40]*dv11*jacob_vy_inv; 
  out[24] += 0.9682458365518543*alpha[39]*dv11*jacob_vy_inv; 
  out[25] += (0.7745966692414834*alpha[41]+0.8660254037844387*alpha[5])*dv11*jacob_vy_inv; 
  out[26] += (0.7745966692414834*alpha[41]+0.8660254037844387*alpha[5])*dv11*jacob_vy_inv; 
  out[27] += 0.9682458365518543*alpha[34]*dv11*jacob_vy_inv; 
  out[28] += (0.7745966692414834*alpha[35]+0.8660254037844387*alpha[9])*dv11*jacob_vy_inv; 
  out[29] += (0.7745966692414834*alpha[36]+0.8660254037844387*alpha[8])*dv11*jacob_vy_inv; 
  out[30] += 0.8660254037844387*alpha[31]*dv11*jacob_vy_inv; 
  out[31] += (0.8660254037844387*alpha[30]+0.8660254037844387*alpha[22]+0.8660254037844387*alpha[21]+0.9682458365518543*alpha[3])*dv11*jacob_vy_inv; 
  out[32] += (0.7745966692414834*alpha[45]+0.8660254037844386*alpha[17])*dv11*jacob_vy_inv; 
  out[33] += (0.7745966692414834*alpha[44]+0.8660254037844386*alpha[18])*dv11*jacob_vy_inv; 
  out[34] += 0.9682458365518543*alpha[27]*dv11*jacob_vy_inv; 
  out[35] += (0.7745966692414834*alpha[28]+0.7745966692414834*alpha[20]+0.8660254037844386*alpha[1])*dv11*jacob_vy_inv; 
  out[36] += (0.7745966692414834*alpha[29]+0.7745966692414834*alpha[19]+0.8660254037844386*alpha[2])*dv11*jacob_vy_inv; 
  out[37] += (0.7745966692414834*alpha[47]+0.8660254037844386*alpha[15])*dv11*jacob_vy_inv; 
  out[38] += (0.7745966692414834*alpha[47]+0.8660254037844386*alpha[15])*dv11*jacob_vy_inv; 
  out[39] += 0.9682458365518543*alpha[24]*dv11*jacob_vy_inv; 
  out[40] += 0.9682458365518543*alpha[23]*dv11*jacob_vy_inv; 
  out[41] += (0.7745966692414834*alpha[26]+0.7745966692414834*alpha[25]+0.8660254037844386*alpha[4])*dv11*jacob_vy_inv; 
  out[42] += (0.7745966692414834*alpha[44]+0.8660254037844386*alpha[18])*dv11*jacob_vy_inv; 
  out[43] += (0.7745966692414834*alpha[45]+0.8660254037844386*alpha[17])*dv11*jacob_vy_inv; 
  out[44] += (0.7745966692414834*alpha[42]+0.7745966692414834*alpha[33]+0.8660254037844387*alpha[6])*dv11*jacob_vy_inv; 
  out[45] += (0.7745966692414834*alpha[43]+0.7745966692414834*alpha[32]+0.8660254037844387*alpha[7])*dv11*jacob_vy_inv; 
  out[46] += 0.9682458365518543*alpha[13]*dv11*jacob_vy_inv; 
  out[47] += (0.7745966692414834*alpha[38]+0.7745966692414834*alpha[37]+0.8660254037844387*alpha[10])*dv11*jacob_vy_inv; 
  } 

} 
