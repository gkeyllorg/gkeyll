#include <gkyl_nc_hamil_gen_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp13_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
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
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  out[0] += 0.9682458365518543*alpha[3]*dv11*jacob_vy_inv; 
  out[1] += 0.9682458365518543*alpha[6]*dv11*jacob_vy_inv; 
  out[2] += 0.9682458365518543*alpha[7]*dv11*jacob_vy_inv; 
  out[3] += (0.8660254037844386*alpha[13]+0.9682458365518543*alpha[0])*dv11*jacob_vy_inv; 
  out[4] += 0.9682458365518543*alpha[10]*dv11*jacob_vy_inv; 
  out[5] += 0.9682458365518543*alpha[15]*dv11*jacob_vy_inv; 
  out[6] += (0.8660254037844387*alpha[23]+0.9682458365518543*alpha[1])*dv11*jacob_vy_inv; 
  out[7] += (0.8660254037844387*alpha[24]+0.9682458365518543*alpha[2])*dv11*jacob_vy_inv; 
  out[8] += 0.9682458365518543*alpha[17]*dv11*jacob_vy_inv; 
  out[9] += 0.9682458365518543*alpha[18]*dv11*jacob_vy_inv; 
  out[10] += (0.8660254037844387*alpha[27]+0.9682458365518543*alpha[4])*dv11*jacob_vy_inv; 
  out[11] += 0.9682458365518543*alpha[21]*dv11*jacob_vy_inv; 
  out[12] += 0.9682458365518543*alpha[22]*dv11*jacob_vy_inv; 
  out[13] += 0.8660254037844386*alpha[3]*dv11*jacob_vy_inv; 
  out[14] += 0.9682458365518543*alpha[30]*dv11*jacob_vy_inv; 
  out[15] += (0.8660254037844386*alpha[34]+0.9682458365518543*alpha[5])*dv11*jacob_vy_inv; 
  out[16] += 0.9682458365518543*alpha[31]*dv11*jacob_vy_inv; 
  out[17] += (0.8660254037844386*alpha[39]+0.9682458365518543*alpha[8])*dv11*jacob_vy_inv; 
  out[18] += (0.8660254037844386*alpha[40]+0.9682458365518543*alpha[9])*dv11*jacob_vy_inv; 
  out[19] += 0.9682458365518543*alpha[32]*dv11*jacob_vy_inv; 
  out[20] += 0.9682458365518543*alpha[33]*dv11*jacob_vy_inv; 
  out[21] += 0.9682458365518543*alpha[11]*dv11*jacob_vy_inv; 
  out[22] += 0.9682458365518543*alpha[12]*dv11*jacob_vy_inv; 
  out[23] += 0.8660254037844387*alpha[6]*dv11*jacob_vy_inv; 
  out[24] += 0.8660254037844387*alpha[7]*dv11*jacob_vy_inv; 
  out[25] += 0.9682458365518543*alpha[37]*dv11*jacob_vy_inv; 
  out[26] += 0.9682458365518543*alpha[38]*dv11*jacob_vy_inv; 
  out[27] += 0.8660254037844387*alpha[10]*dv11*jacob_vy_inv; 
  out[28] += 0.9682458365518543*alpha[42]*dv11*jacob_vy_inv; 
  out[29] += 0.9682458365518543*alpha[43]*dv11*jacob_vy_inv; 
  out[30] += 0.9682458365518543*alpha[14]*dv11*jacob_vy_inv; 
  out[31] += (0.8660254037844387*alpha[46]+0.9682458365518543*alpha[16])*dv11*jacob_vy_inv; 
  out[32] += 0.9682458365518543*alpha[19]*dv11*jacob_vy_inv; 
  out[33] += 0.9682458365518543*alpha[20]*dv11*jacob_vy_inv; 
  out[34] += 0.8660254037844386*alpha[15]*dv11*jacob_vy_inv; 
  out[35] += 0.9682458365518543*alpha[44]*dv11*jacob_vy_inv; 
  out[36] += 0.9682458365518543*alpha[45]*dv11*jacob_vy_inv; 
  out[37] += 0.9682458365518543*alpha[25]*dv11*jacob_vy_inv; 
  out[38] += 0.9682458365518543*alpha[26]*dv11*jacob_vy_inv; 
  out[39] += 0.8660254037844386*alpha[17]*dv11*jacob_vy_inv; 
  out[40] += 0.8660254037844386*alpha[18]*dv11*jacob_vy_inv; 
  out[41] += 0.9682458365518543*alpha[47]*dv11*jacob_vy_inv; 
  out[42] += 0.9682458365518543*alpha[28]*dv11*jacob_vy_inv; 
  out[43] += 0.9682458365518543*alpha[29]*dv11*jacob_vy_inv; 
  out[44] += 0.9682458365518543*alpha[35]*dv11*jacob_vy_inv; 
  out[45] += 0.9682458365518543*alpha[36]*dv11*jacob_vy_inv; 
  out[46] += 0.8660254037844387*alpha[31]*dv11*jacob_vy_inv; 
  out[47] += 0.9682458365518543*alpha[41]*dv11*jacob_vy_inv; 
  } 

  if (dir == 3) { 
  } 

} 
