#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp20_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[2]; 
  const double jacob_vz_inv = 1.0/jacob_vel[4]; 
  
  if (dir == 0) { 
  out[0] += 0.3061862178478971*alpha[13]*dx10*jacob_cx_inv; 
  out[1] += 0.3061862178478971*alpha[20]*dx10*jacob_cx_inv; 
  out[2] += 0.3061862178478971*alpha[5]*dx10*jacob_cx_inv; 
  out[3] += 0.3061862178478971*alpha[22]*dx10*jacob_cx_inv; 
  out[4] += 0.3061862178478971*alpha[24]*dx10*jacob_cx_inv; 
  out[5] += 0.3061862178478971*alpha[2]*dx10*jacob_cx_inv; 
  out[6] += 0.3061862178478971*alpha[12]*dx10*jacob_cx_inv; 
  out[7] += 0.3061862178478971*alpha[27]*dx10*jacob_cx_inv; 
  out[8] += 0.3061862178478971*alpha[14]*dx10*jacob_cx_inv; 
  out[9] += 0.3061862178478971*alpha[28]*dx10*jacob_cx_inv; 
  out[10] += 0.3061862178478971*alpha[15]*dx10*jacob_cx_inv; 
  out[11] += 0.3061862178478971*alpha[30]*dx10*jacob_cx_inv; 
  out[12] += 0.3061862178478971*alpha[6]*dx10*jacob_cx_inv; 
  out[13] += 0.3061862178478971*alpha[0]*dx10*jacob_cx_inv; 
  out[14] += 0.3061862178478971*alpha[8]*dx10*jacob_cx_inv; 
  out[15] += 0.3061862178478971*alpha[10]*dx10*jacob_cx_inv; 
  out[16] += 0.3061862178478971*alpha[21]*dx10*jacob_cx_inv; 
  out[17] += 0.3061862178478971*alpha[23]*dx10*jacob_cx_inv; 
  out[18] += 0.3061862178478971*alpha[31]*dx10*jacob_cx_inv; 
  out[19] += 0.3061862178478971*alpha[25]*dx10*jacob_cx_inv; 
  out[20] += 0.3061862178478971*alpha[1]*dx10*jacob_cx_inv; 
  out[21] += 0.3061862178478971*alpha[16]*dx10*jacob_cx_inv; 
  out[22] += 0.3061862178478971*alpha[3]*dx10*jacob_cx_inv; 
  out[23] += 0.3061862178478971*alpha[17]*dx10*jacob_cx_inv; 
  out[24] += 0.3061862178478971*alpha[4]*dx10*jacob_cx_inv; 
  out[25] += 0.3061862178478971*alpha[19]*dx10*jacob_cx_inv; 
  out[26] += 0.3061862178478971*alpha[29]*dx10*jacob_cx_inv; 
  out[27] += 0.3061862178478971*alpha[7]*dx10*jacob_cx_inv; 
  out[28] += 0.3061862178478971*alpha[9]*dx10*jacob_cx_inv; 
  out[29] += 0.3061862178478971*alpha[26]*dx10*jacob_cx_inv; 
  out[30] += 0.3061862178478971*alpha[11]*dx10*jacob_cx_inv; 
  out[31] += 0.3061862178478971*alpha[18]*dx10*jacob_cx_inv; 
  } 

  if (dir == 1) { 
  out[0] += 0.3061862178478971*alpha[12]*dx11*jacob_cy_inv; 
  out[1] += 0.3061862178478971*alpha[5]*dx11*jacob_cy_inv; 
  out[2] += 0.3061862178478971*alpha[20]*dx11*jacob_cy_inv; 
  out[3] += 0.3061862178478971*alpha[21]*dx11*jacob_cy_inv; 
  out[4] += 0.3061862178478971*alpha[23]*dx11*jacob_cy_inv; 
  out[5] += 0.3061862178478971*alpha[1]*dx11*jacob_cy_inv; 
  out[6] += 0.3061862178478971*alpha[13]*dx11*jacob_cy_inv; 
  out[7] += 0.3061862178478971*alpha[14]*dx11*jacob_cy_inv; 
  out[8] += 0.3061862178478971*alpha[27]*dx11*jacob_cy_inv; 
  out[9] += 0.3061862178478971*alpha[15]*dx11*jacob_cy_inv; 
  out[10] += 0.3061862178478971*alpha[28]*dx11*jacob_cy_inv; 
  out[11] += 0.3061862178478971*alpha[29]*dx11*jacob_cy_inv; 
  out[12] += 0.3061862178478971*alpha[0]*dx11*jacob_cy_inv; 
  out[13] += 0.3061862178478971*alpha[6]*dx11*jacob_cy_inv; 
  out[14] += 0.3061862178478971*alpha[7]*dx11*jacob_cy_inv; 
  out[15] += 0.3061862178478971*alpha[9]*dx11*jacob_cy_inv; 
  out[16] += 0.3061862178478971*alpha[22]*dx11*jacob_cy_inv; 
  out[17] += 0.3061862178478971*alpha[24]*dx11*jacob_cy_inv; 
  out[18] += 0.3061862178478971*alpha[25]*dx11*jacob_cy_inv; 
  out[19] += 0.3061862178478971*alpha[31]*dx11*jacob_cy_inv; 
  out[20] += 0.3061862178478971*alpha[2]*dx11*jacob_cy_inv; 
  out[21] += 0.3061862178478971*alpha[3]*dx11*jacob_cy_inv; 
  out[22] += 0.3061862178478971*alpha[16]*dx11*jacob_cy_inv; 
  out[23] += 0.3061862178478971*alpha[4]*dx11*jacob_cy_inv; 
  out[24] += 0.3061862178478971*alpha[17]*dx11*jacob_cy_inv; 
  out[25] += 0.3061862178478971*alpha[18]*dx11*jacob_cy_inv; 
  out[26] += 0.3061862178478971*alpha[30]*dx11*jacob_cy_inv; 
  out[27] += 0.3061862178478971*alpha[8]*dx11*jacob_cy_inv; 
  out[28] += 0.3061862178478971*alpha[10]*dx11*jacob_cy_inv; 
  out[29] += 0.3061862178478971*alpha[11]*dx11*jacob_cy_inv; 
  out[30] += 0.3061862178478971*alpha[26]*dx11*jacob_cy_inv; 
  out[31] += 0.3061862178478971*alpha[19]*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  out[0] += 0.3061862178478971*alpha[6]*dv12*jacob_vz_inv; 
  out[1] += 0.3061862178478971*alpha[2]*dv12*jacob_vz_inv; 
  out[2] += 0.3061862178478971*alpha[1]*dv12*jacob_vz_inv; 
  out[3] += 0.3061862178478971*alpha[16]*dv12*jacob_vz_inv; 
  out[4] += 0.3061862178478971*alpha[17]*dv12*jacob_vz_inv; 
  out[5] += 0.3061862178478971*alpha[20]*dv12*jacob_vz_inv; 
  out[6] += 0.3061862178478971*alpha[0]*dv12*jacob_vz_inv; 
  out[7] += 0.3061862178478971*alpha[8]*dv12*jacob_vz_inv; 
  out[8] += 0.3061862178478971*alpha[7]*dv12*jacob_vz_inv; 
  out[9] += 0.3061862178478971*alpha[10]*dv12*jacob_vz_inv; 
  out[10] += 0.3061862178478971*alpha[9]*dv12*jacob_vz_inv; 
  out[11] += 0.3061862178478971*alpha[26]*dv12*jacob_vz_inv; 
  out[12] += 0.3061862178478971*alpha[13]*dv12*jacob_vz_inv; 
  out[13] += 0.3061862178478971*alpha[12]*dv12*jacob_vz_inv; 
  out[14] += 0.3061862178478971*alpha[27]*dv12*jacob_vz_inv; 
  out[15] += 0.3061862178478971*alpha[28]*dv12*jacob_vz_inv; 
  out[16] += 0.3061862178478971*alpha[3]*dv12*jacob_vz_inv; 
  out[17] += 0.3061862178478971*alpha[4]*dv12*jacob_vz_inv; 
  out[18] += 0.3061862178478971*alpha[19]*dv12*jacob_vz_inv; 
  out[19] += 0.3061862178478971*alpha[18]*dv12*jacob_vz_inv; 
  out[20] += 0.3061862178478971*alpha[5]*dv12*jacob_vz_inv; 
  out[21] += 0.3061862178478971*alpha[22]*dv12*jacob_vz_inv; 
  out[22] += 0.3061862178478971*alpha[21]*dv12*jacob_vz_inv; 
  out[23] += 0.3061862178478971*alpha[24]*dv12*jacob_vz_inv; 
  out[24] += 0.3061862178478971*alpha[23]*dv12*jacob_vz_inv; 
  out[25] += 0.3061862178478971*alpha[31]*dv12*jacob_vz_inv; 
  out[26] += 0.3061862178478971*alpha[11]*dv12*jacob_vz_inv; 
  out[27] += 0.3061862178478971*alpha[14]*dv12*jacob_vz_inv; 
  out[28] += 0.3061862178478971*alpha[15]*dv12*jacob_vz_inv; 
  out[29] += 0.3061862178478971*alpha[30]*dv12*jacob_vz_inv; 
  out[30] += 0.3061862178478971*alpha[29]*dv12*jacob_vz_inv; 
  out[31] += 0.3061862178478971*alpha[25]*dv12*jacob_vz_inv; 
  } 

} 
