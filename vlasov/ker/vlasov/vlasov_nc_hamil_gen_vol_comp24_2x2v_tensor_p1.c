#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_2x2v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double *jacob_vy = &jacob_vel[3]; 
  
  if (dir == 0) { 
  out[0] += 0.43301270189221935*alpha[18]*dx10*jacob_cx_inv; 
  out[1] += 0.4330127018922193*alpha[24]*dx10*jacob_cx_inv; 
  out[2] += 0.4330127018922193*alpha[11]*dx10*jacob_cx_inv; 
  out[3] += 0.38729833462074165*alpha[7]*dx10*jacob_cx_inv; 
  out[4] += 0.4330127018922193*alpha[26]*dx10*jacob_cx_inv; 
  out[5] += 0.43301270189221935*alpha[17]*dx10*jacob_cx_inv; 
  out[6] += 0.38729833462074165*alpha[13]*dx10*jacob_cx_inv; 
  out[7] += 0.38729833462074165*alpha[3]*dx10*jacob_cx_inv; 
  out[8] += 0.43301270189221935*alpha[31]*dx10*jacob_cx_inv; 
  out[9] += 0.43301270189221935*alpha[19]*dx10*jacob_cx_inv; 
  out[10] += 0.38729833462074165*alpha[16]*dx10*jacob_cx_inv; 
  out[11] += (0.27664166758624403*alpha[18]+0.4330127018922193*alpha[2])*dx10*jacob_cx_inv; 
  out[12] += 0.4330127018922193*alpha[34]*dx10*jacob_cx_inv; 
  out[13] += 0.38729833462074165*alpha[6]*dx10*jacob_cx_inv; 
  out[14] += 0.4330127018922193*alpha[25]*dx10*jacob_cx_inv; 
  out[15] += 0.38729833462074165*alpha[23]*dx10*jacob_cx_inv; 
  out[16] += 0.38729833462074165*alpha[10]*dx10*jacob_cx_inv; 
  out[17] += (0.27664166758624403*alpha[24]+0.43301270189221935*alpha[5])*dx10*jacob_cx_inv; 
  out[18] += (0.27664166758624403*alpha[11]+0.43301270189221935*alpha[0])*dx10*jacob_cx_inv; 
  out[19] += (0.27664166758624403*alpha[26]+0.43301270189221935*alpha[9])*dx10*jacob_cx_inv; 
  out[20] += 0.43301270189221935*alpha[35]*dx10*jacob_cx_inv; 
  out[21] += 0.43301270189221935*alpha[30]*dx10*jacob_cx_inv; 
  out[22] += 0.3872983346207417*alpha[29]*dx10*jacob_cx_inv; 
  out[23] += 0.38729833462074165*alpha[15]*dx10*jacob_cx_inv; 
  out[24] += (0.27664166758624403*alpha[17]+0.4330127018922193*alpha[1])*dx10*jacob_cx_inv; 
  out[25] += (0.27664166758624403*alpha[31]+0.4330127018922193*alpha[14])*dx10*jacob_cx_inv; 
  out[26] += (0.27664166758624403*alpha[19]+0.4330127018922193*alpha[4])*dx10*jacob_cx_inv; 
  out[27] += 0.4330127018922193*alpha[33]*dx10*jacob_cx_inv; 
  out[28] += 0.3872983346207417*alpha[32]*dx10*jacob_cx_inv; 
  out[29] += 0.3872983346207417*alpha[22]*dx10*jacob_cx_inv; 
  out[30] += (0.27664166758624403*alpha[34]+0.43301270189221935*alpha[21])*dx10*jacob_cx_inv; 
  out[31] += (0.27664166758624403*alpha[25]+0.43301270189221935*alpha[8])*dx10*jacob_cx_inv; 
  out[32] += 0.3872983346207417*alpha[28]*dx10*jacob_cx_inv; 
  out[33] += (0.27664166758624403*alpha[35]+0.4330127018922193*alpha[27])*dx10*jacob_cx_inv; 
  out[34] += (0.27664166758624403*alpha[30]+0.4330127018922193*alpha[12])*dx10*jacob_cx_inv; 
  out[35] += (0.27664166758624403*alpha[33]+0.43301270189221935*alpha[20])*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  out[0] += 0.43301270189221935*alpha[17]*dx11*jacob_cy_inv; 
  out[1] += 0.4330127018922193*alpha[11]*dx11*jacob_cy_inv; 
  out[2] += 0.4330127018922193*alpha[24]*dx11*jacob_cy_inv; 
  out[3] += 0.38729833462074165*alpha[6]*dx11*jacob_cy_inv; 
  out[4] += 0.4330127018922193*alpha[25]*dx11*jacob_cy_inv; 
  out[5] += 0.43301270189221935*alpha[18]*dx11*jacob_cy_inv; 
  out[6] += 0.38729833462074165*alpha[3]*dx11*jacob_cy_inv; 
  out[7] += 0.38729833462074165*alpha[13]*dx11*jacob_cy_inv; 
  out[8] += 0.43301270189221935*alpha[19]*dx11*jacob_cy_inv; 
  out[9] += 0.43301270189221935*alpha[31]*dx11*jacob_cy_inv; 
  out[10] += 0.38729833462074165*alpha[15]*dx11*jacob_cy_inv; 
  out[11] += (0.27664166758624403*alpha[17]+0.4330127018922193*alpha[1])*dx11*jacob_cy_inv; 
  out[12] += 0.4330127018922193*alpha[33]*dx11*jacob_cy_inv; 
  out[13] += 0.38729833462074165*alpha[7]*dx11*jacob_cy_inv; 
  out[14] += 0.4330127018922193*alpha[26]*dx11*jacob_cy_inv; 
  out[15] += 0.38729833462074165*alpha[10]*dx11*jacob_cy_inv; 
  out[16] += 0.38729833462074165*alpha[23]*dx11*jacob_cy_inv; 
  out[17] += (0.27664166758624403*alpha[11]+0.43301270189221935*alpha[0])*dx11*jacob_cy_inv; 
  out[18] += (0.27664166758624403*alpha[24]+0.43301270189221935*alpha[5])*dx11*jacob_cy_inv; 
  out[19] += (0.27664166758624403*alpha[25]+0.43301270189221935*alpha[8])*dx11*jacob_cy_inv; 
  out[20] += 0.43301270189221935*alpha[30]*dx11*jacob_cy_inv; 
  out[21] += 0.43301270189221935*alpha[35]*dx11*jacob_cy_inv; 
  out[22] += 0.3872983346207417*alpha[28]*dx11*jacob_cy_inv; 
  out[23] += 0.38729833462074165*alpha[16]*dx11*jacob_cy_inv; 
  out[24] += (0.27664166758624403*alpha[18]+0.4330127018922193*alpha[2])*dx11*jacob_cy_inv; 
  out[25] += (0.27664166758624403*alpha[19]+0.4330127018922193*alpha[4])*dx11*jacob_cy_inv; 
  out[26] += (0.27664166758624403*alpha[31]+0.4330127018922193*alpha[14])*dx11*jacob_cy_inv; 
  out[27] += 0.4330127018922193*alpha[34]*dx11*jacob_cy_inv; 
  out[28] += 0.3872983346207417*alpha[22]*dx11*jacob_cy_inv; 
  out[29] += 0.3872983346207417*alpha[32]*dx11*jacob_cy_inv; 
  out[30] += (0.27664166758624403*alpha[33]+0.43301270189221935*alpha[20])*dx11*jacob_cy_inv; 
  out[31] += (0.27664166758624403*alpha[26]+0.43301270189221935*alpha[9])*dx11*jacob_cy_inv; 
  out[32] += 0.3872983346207417*alpha[29]*dx11*jacob_cy_inv; 
  out[33] += (0.27664166758624403*alpha[30]+0.4330127018922193*alpha[12])*dx11*jacob_cy_inv; 
  out[34] += (0.27664166758624403*alpha[35]+0.4330127018922193*alpha[27])*dx11*jacob_cy_inv; 
  out[35] += (0.27664166758624403*alpha[34]+0.43301270189221935*alpha[21])*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  out[0] += 0.9682458365518543*alpha[13]*dv10; 
  out[1] += 0.9682458365518543*alpha[7]*dv10; 
  out[2] += 0.9682458365518543*alpha[6]*dv10; 
  out[3] += (0.8660254037844386*alpha[24]+0.9682458365518543*alpha[5])*dv10; 
  out[4] += 0.9682458365518543*alpha[23]*dv10; 
  out[5] += 0.9682458365518543*alpha[3]*dv10; 
  out[6] += (0.8660254037844387*alpha[18]+0.9682458365518543*alpha[2])*dv10; 
  out[7] += (0.8660254037844387*alpha[17]+0.9682458365518543*alpha[1])*dv10; 
  out[8] += 0.9682458365518543*alpha[16]*dv10; 
  out[9] += 0.9682458365518543*alpha[15]*dv10; 
  out[10] += (0.8660254037844387*alpha[31]+0.9682458365518543*alpha[14])*dv10; 
  out[11] += 0.8660254037844386*alpha[13]*dv10; 
  out[12] += 0.9682458365518543*alpha[32]*dv10; 
  out[13] += (0.8660254037844386*alpha[11]+0.9682458365518543*alpha[0])*dv10; 
  out[14] += 0.9682458365518543*alpha[10]*dv10; 
  out[15] += (0.8660254037844386*alpha[26]+0.9682458365518543*alpha[9])*dv10; 
  out[16] += (0.8660254037844386*alpha[25]+0.9682458365518543*alpha[8])*dv10; 
  out[17] += 0.8660254037844387*alpha[7]*dv10; 
  out[18] += 0.8660254037844387*alpha[6]*dv10; 
  out[19] += 0.8660254037844387*alpha[23]*dv10; 
  out[20] += 0.9682458365518543*alpha[29]*dv10; 
  out[21] += 0.9682458365518543*alpha[28]*dv10; 
  out[22] += (0.8660254037844387*alpha[35]+0.9682458365518543*alpha[27])*dv10; 
  out[23] += (0.8660254037844387*alpha[19]+0.9682458365518543*alpha[4])*dv10; 
  out[24] += 0.8660254037844386*alpha[3]*dv10; 
  out[25] += 0.8660254037844386*alpha[16]*dv10; 
  out[26] += 0.8660254037844386*alpha[15]*dv10; 
  out[27] += 0.9682458365518543*alpha[22]*dv10; 
  out[28] += (0.8660254037844386*alpha[34]+0.9682458365518543*alpha[21])*dv10; 
  out[29] += (0.8660254037844386*alpha[33]+0.9682458365518543*alpha[20])*dv10; 
  out[30] += 0.8660254037844387*alpha[32]*dv10; 
  out[31] += 0.8660254037844387*alpha[10]*dv10; 
  out[32] += (0.8660254037844387*alpha[30]+0.9682458365518543*alpha[12])*dv10; 
  out[33] += 0.8660254037844386*alpha[29]*dv10; 
  out[34] += 0.8660254037844386*alpha[28]*dv10; 
  out[35] += 0.8660254037844387*alpha[22]*dv10; 
  } 

  if (dir == 3) { 
  } 

} 
