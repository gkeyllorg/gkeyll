#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
   const double *jacob_pos, const double *jacob_vel,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dx12 = 2.0/dxv[2]; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double jacob_cz_inv = 1.0/jacob_pos[4]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[2]; 
  const double jacob_vz_inv = 1.0/jacob_vel[4]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  out[0] += 0.21650635094610965*alpha[0]*dx11*jacob_cy_inv; 
  out[1] += 0.21650635094610965*alpha[1]*dx11*jacob_cy_inv; 
  out[2] += 0.21650635094610965*alpha[2]*dx11*jacob_cy_inv; 
  out[3] += 0.21650635094610965*alpha[3]*dx11*jacob_cy_inv; 
  out[4] += 0.21650635094610965*alpha[4]*dx11*jacob_cy_inv; 
  out[5] += 0.21650635094610965*alpha[5]*dx11*jacob_cy_inv; 
  out[6] += 0.21650635094610965*alpha[6]*dx11*jacob_cy_inv; 
  out[7] += 0.21650635094610965*alpha[7]*dx11*jacob_cy_inv; 
  out[8] += 0.21650635094610965*alpha[8]*dx11*jacob_cy_inv; 
  out[9] += 0.21650635094610965*alpha[9]*dx11*jacob_cy_inv; 
  out[10] += 0.21650635094610965*alpha[10]*dx11*jacob_cy_inv; 
  out[11] += 0.21650635094610965*alpha[11]*dx11*jacob_cy_inv; 
  out[12] += 0.21650635094610965*alpha[12]*dx11*jacob_cy_inv; 
  out[13] += 0.21650635094610965*alpha[13]*dx11*jacob_cy_inv; 
  out[14] += 0.21650635094610965*alpha[14]*dx11*jacob_cy_inv; 
  out[15] += 0.21650635094610965*alpha[15]*dx11*jacob_cy_inv; 
  out[16] += 0.21650635094610965*alpha[16]*dx11*jacob_cy_inv; 
  out[17] += 0.21650635094610965*alpha[17]*dx11*jacob_cy_inv; 
  out[18] += 0.21650635094610965*alpha[18]*dx11*jacob_cy_inv; 
  out[19] += 0.21650635094610965*alpha[19]*dx11*jacob_cy_inv; 
  out[20] += 0.21650635094610965*alpha[20]*dx11*jacob_cy_inv; 
  out[21] += 0.21650635094610965*alpha[21]*dx11*jacob_cy_inv; 
  out[22] += 0.21650635094610965*alpha[22]*dx11*jacob_cy_inv; 
  out[23] += 0.21650635094610965*alpha[23]*dx11*jacob_cy_inv; 
  out[24] += 0.21650635094610965*alpha[24]*dx11*jacob_cy_inv; 
  out[25] += 0.21650635094610965*alpha[25]*dx11*jacob_cy_inv; 
  out[26] += 0.21650635094610965*alpha[26]*dx11*jacob_cy_inv; 
  out[27] += 0.21650635094610965*alpha[27]*dx11*jacob_cy_inv; 
  out[28] += 0.21650635094610965*alpha[28]*dx11*jacob_cy_inv; 
  out[29] += 0.21650635094610965*alpha[29]*dx11*jacob_cy_inv; 
  out[30] += 0.21650635094610965*alpha[30]*dx11*jacob_cy_inv; 
  out[31] += 0.21650635094610965*alpha[31]*dx11*jacob_cy_inv; 
  out[32] += 0.21650635094610965*alpha[32]*dx11*jacob_cy_inv; 
  out[33] += 0.21650635094610965*alpha[33]*dx11*jacob_cy_inv; 
  out[34] += 0.21650635094610965*alpha[34]*dx11*jacob_cy_inv; 
  out[35] += 0.21650635094610965*alpha[35]*dx11*jacob_cy_inv; 
  out[36] += 0.21650635094610965*alpha[36]*dx11*jacob_cy_inv; 
  out[37] += 0.21650635094610965*alpha[37]*dx11*jacob_cy_inv; 
  out[38] += 0.21650635094610965*alpha[38]*dx11*jacob_cy_inv; 
  out[39] += 0.21650635094610965*alpha[39]*dx11*jacob_cy_inv; 
  out[40] += 0.21650635094610965*alpha[40]*dx11*jacob_cy_inv; 
  out[41] += 0.21650635094610965*alpha[41]*dx11*jacob_cy_inv; 
  out[42] += 0.21650635094610965*alpha[42]*dx11*jacob_cy_inv; 
  out[43] += 0.21650635094610965*alpha[43]*dx11*jacob_cy_inv; 
  out[44] += 0.21650635094610965*alpha[44]*dx11*jacob_cy_inv; 
  out[45] += 0.21650635094610965*alpha[45]*dx11*jacob_cy_inv; 
  out[46] += 0.21650635094610965*alpha[46]*dx11*jacob_cy_inv; 
  out[47] += 0.21650635094610965*alpha[47]*dx11*jacob_cy_inv; 
  out[48] += 0.21650635094610965*alpha[48]*dx11*jacob_cy_inv; 
  out[49] += 0.21650635094610965*alpha[49]*dx11*jacob_cy_inv; 
  out[50] += 0.21650635094610965*alpha[50]*dx11*jacob_cy_inv; 
  out[51] += 0.21650635094610965*alpha[51]*dx11*jacob_cy_inv; 
  out[52] += 0.21650635094610965*alpha[52]*dx11*jacob_cy_inv; 
  out[53] += 0.21650635094610965*alpha[53]*dx11*jacob_cy_inv; 
  out[54] += 0.21650635094610965*alpha[54]*dx11*jacob_cy_inv; 
  out[55] += 0.21650635094610965*alpha[55]*dx11*jacob_cy_inv; 
  out[56] += 0.21650635094610965*alpha[56]*dx11*jacob_cy_inv; 
  out[57] += 0.21650635094610965*alpha[57]*dx11*jacob_cy_inv; 
  out[58] += 0.21650635094610965*alpha[58]*dx11*jacob_cy_inv; 
  out[59] += 0.21650635094610965*alpha[59]*dx11*jacob_cy_inv; 
  out[60] += 0.21650635094610965*alpha[60]*dx11*jacob_cy_inv; 
  out[61] += 0.21650635094610965*alpha[61]*dx11*jacob_cy_inv; 
  out[62] += 0.21650635094610965*alpha[62]*dx11*jacob_cy_inv; 
  out[63] += 0.21650635094610965*alpha[63]*dx11*jacob_cy_inv; 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  } 

  if (dir == 5) { 
  } 

} 
