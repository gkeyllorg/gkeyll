#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[62]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[78]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[2]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[8]*subvol; 
  out[6] += 0.25*win[0]*fin[10]*subvol; 
  out[7] += 0.25*win[0]*fin[11]*subvol; 
  out[8] += 0.25*win[0]*fin[13]*subvol; 
  out[9] += 0.25*win[0]*fin[14]*subvol; 
  out[10] += 0.25*win[0]*fin[15]*subvol; 
  out[11] += 0.25*win[0]*fin[19]*subvol; 
  out[12] += 0.25*win[0]*fin[22]*subvol; 
  out[13] += 0.25*win[0]*fin[24]*subvol; 
  out[14] += 0.25*win[0]*fin[25]*subvol; 
  out[15] += 0.25*win[0]*fin[30]*subvol; 
  out[16] += 0.25*win[0]*fin[32]*subvol; 
  out[17] += 0.25*win[0]*fin[34]*subvol; 
  out[18] += 0.25*win[0]*fin[35]*subvol; 
  out[19] += 0.25*win[0]*fin[36]*subvol; 
  out[20] += 0.25*win[0]*fin[39]*subvol; 
  out[21] += 0.25*win[0]*fin[41]*subvol; 
  out[22] += 0.25*win[0]*fin[42]*subvol; 
  out[23] += 0.25*win[0]*fin[46]*subvol; 
  out[24] += 0.25*win[0]*fin[48]*subvol; 
  out[25] += 0.25*win[0]*fin[50]*subvol; 
  out[26] += 0.25*win[0]*fin[51]*subvol; 
  out[27] += 0.25*win[0]*fin[52]*subvol; 
  out[28] += 0.25*win[0]*fin[55]*subvol; 
  out[29] += 0.25*win[0]*fin[57]*subvol; 
  out[30] += 0.25*win[0]*fin[58]*subvol; 
  out[31] += 0.25*win[0]*fin[62]*subvol; 
  out[32] += 0.25*win[0]*fin[64]*subvol; 
  out[33] += 0.25*win[0]*fin[66]*subvol; 
  out[34] += 0.25*win[0]*fin[67]*subvol; 
  out[35] += 0.25*win[0]*fin[68]*subvol; 
  out[36] += 0.25*win[0]*fin[71]*subvol; 
  out[37] += 0.25*win[0]*fin[73]*subvol; 
  out[38] += 0.25*win[0]*fin[74]*subvol; 
  out[39] += 0.25*win[0]*fin[78]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[22], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[30], 0.3535533905932737*win[0]*fin[30]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[35], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[39], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[41], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[46], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[55], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[57], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[62], 0.3535533905932737*win[0]*fin[62]*subvol); 
  atomicAdd(&out[64], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[66], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[67], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[68], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[71], 0.3535533905932737*win[0]*fin[71]*subvol); 
  atomicAdd(&out[73], 0.3535533905932737*win[0]*fin[73]*subvol); 
  atomicAdd(&out[74], 0.3535533905932737*win[0]*fin[74]*subvol); 
  atomicAdd(&out[78], 0.3535533905932737*win[0]*fin[78]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[22] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[24] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[25] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[30] += 0.3535533905932737*win[0]*fin[30]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[35] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[39] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[41] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[42] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[46] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[50] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[51] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[52] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[55] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[57] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[58] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[62] += 0.3535533905932737*win[0]*fin[62]*subvol; 
  out[64] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[66] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[67] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[68] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[71] += 0.3535533905932737*win[0]*fin[71]*subvol; 
  out[73] += 0.3535533905932737*win[0]*fin[73]*subvol; 
  out[74] += 0.3535533905932737*win[0]*fin[74]*subvol; 
  out[78] += 0.3535533905932737*win[0]*fin[78]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[61]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[77]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[3]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[7]*subvol; 
  out[6] += 0.25*win[0]*fin[9]*subvol; 
  out[7] += 0.25*win[0]*fin[11]*subvol; 
  out[8] += 0.25*win[0]*fin[12]*subvol; 
  out[9] += 0.25*win[0]*fin[14]*subvol; 
  out[10] += 0.25*win[0]*fin[15]*subvol; 
  out[11] += 0.25*win[0]*fin[18]*subvol; 
  out[12] += 0.25*win[0]*fin[21]*subvol; 
  out[13] += 0.25*win[0]*fin[23]*subvol; 
  out[14] += 0.25*win[0]*fin[25]*subvol; 
  out[15] += 0.25*win[0]*fin[29]*subvol; 
  out[16] += 0.25*win[0]*fin[32]*subvol; 
  out[17] += 0.25*win[0]*fin[33]*subvol; 
  out[18] += 0.25*win[0]*fin[35]*subvol; 
  out[19] += 0.25*win[0]*fin[36]*subvol; 
  out[20] += 0.25*win[0]*fin[38]*subvol; 
  out[21] += 0.25*win[0]*fin[40]*subvol; 
  out[22] += 0.25*win[0]*fin[42]*subvol; 
  out[23] += 0.25*win[0]*fin[45]*subvol; 
  out[24] += 0.25*win[0]*fin[48]*subvol; 
  out[25] += 0.25*win[0]*fin[49]*subvol; 
  out[26] += 0.25*win[0]*fin[51]*subvol; 
  out[27] += 0.25*win[0]*fin[52]*subvol; 
  out[28] += 0.25*win[0]*fin[54]*subvol; 
  out[29] += 0.25*win[0]*fin[56]*subvol; 
  out[30] += 0.25*win[0]*fin[58]*subvol; 
  out[31] += 0.25*win[0]*fin[61]*subvol; 
  out[32] += 0.25*win[0]*fin[64]*subvol; 
  out[33] += 0.25*win[0]*fin[65]*subvol; 
  out[34] += 0.25*win[0]*fin[67]*subvol; 
  out[35] += 0.25*win[0]*fin[68]*subvol; 
  out[36] += 0.25*win[0]*fin[70]*subvol; 
  out[37] += 0.25*win[0]*fin[72]*subvol; 
  out[38] += 0.25*win[0]*fin[74]*subvol; 
  out[39] += 0.25*win[0]*fin[77]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[21], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[23], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[25], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[29], 0.3535533905932737*win[0]*fin[29]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[38], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[40], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[42], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[45], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[51], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[54], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[56], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[58], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[61], 0.3535533905932737*win[0]*fin[61]*subvol); 
  atomicAdd(&out[64], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[67], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[68], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[70], 0.3535533905932737*win[0]*fin[70]*subvol); 
  atomicAdd(&out[72], 0.3535533905932737*win[0]*fin[72]*subvol); 
  atomicAdd(&out[74], 0.3535533905932737*win[0]*fin[74]*subvol); 
  atomicAdd(&out[77], 0.3535533905932737*win[0]*fin[77]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[21] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[23] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[25] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[29] += 0.3535533905932737*win[0]*fin[29]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[35] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[38] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[40] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[42] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[45] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[49] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[51] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[52] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[54] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[56] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[58] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[61] += 0.3535533905932737*win[0]*fin[61]*subvol; 
  out[64] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[65] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[67] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[68] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[70] += 0.3535533905932737*win[0]*fin[70]*subvol; 
  out[72] += 0.3535533905932737*win[0]*fin[72]*subvol; 
  out[74] += 0.3535533905932737*win[0]*fin[74]*subvol; 
  out[77] += 0.3535533905932737*win[0]*fin[77]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[74]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[74]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
  atomicAdd(&out[25], 0.7071067811865475*win[0]*fin[25]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[35], 0.7071067811865475*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[42], 0.7071067811865475*win[0]*fin[42]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[51], 0.7071067811865475*win[0]*fin[51]*subvol); 
  atomicAdd(&out[52], 0.7071067811865475*win[0]*fin[52]*subvol); 
  atomicAdd(&out[58], 0.7071067811865475*win[0]*fin[58]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[67], 0.7071067811865475*win[0]*fin[67]*subvol); 
  atomicAdd(&out[68], 0.7071067811865475*win[0]*fin[68]*subvol); 
  atomicAdd(&out[74], 0.7071067811865475*win[0]*fin[74]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[25] += 0.7071067811865475*win[0]*fin[25]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[35] += 0.7071067811865475*win[0]*fin[35]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[42] += 0.7071067811865475*win[0]*fin[42]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[51] += 0.7071067811865475*win[0]*fin[51]*subvol; 
  out[52] += 0.7071067811865475*win[0]*fin[52]*subvol; 
  out[58] += 0.7071067811865475*win[0]*fin[58]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[67] += 0.7071067811865475*win[0]*fin[67]*subvol; 
  out[68] += 0.7071067811865475*win[0]*fin[68]*subvol; 
  out[74] += 0.7071067811865475*win[0]*fin[74]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[15]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[60]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[76]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[4]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[6]*subvol; 
  out[6] += 0.25*win[0]*fin[9]*subvol; 
  out[7] += 0.25*win[0]*fin[10]*subvol; 
  out[8] += 0.25*win[0]*fin[12]*subvol; 
  out[9] += 0.25*win[0]*fin[13]*subvol; 
  out[10] += 0.25*win[0]*fin[15]*subvol; 
  out[11] += 0.25*win[0]*fin[17]*subvol; 
  out[12] += 0.25*win[0]*fin[20]*subvol; 
  out[13] += 0.25*win[0]*fin[23]*subvol; 
  out[14] += 0.25*win[0]*fin[24]*subvol; 
  out[15] += 0.25*win[0]*fin[28]*subvol; 
  out[16] += 0.25*win[0]*fin[48]*subvol; 
  out[17] += 0.25*win[0]*fin[49]*subvol; 
  out[18] += 0.25*win[0]*fin[50]*subvol; 
  out[19] += 0.25*win[0]*fin[52]*subvol; 
  out[20] += 0.25*win[0]*fin[53]*subvol; 
  out[21] += 0.25*win[0]*fin[56]*subvol; 
  out[22] += 0.25*win[0]*fin[57]*subvol; 
  out[23] += 0.25*win[0]*fin[60]*subvol; 
  out[24] += 0.25*win[0]*fin[64]*subvol; 
  out[25] += 0.25*win[0]*fin[65]*subvol; 
  out[26] += 0.25*win[0]*fin[66]*subvol; 
  out[27] += 0.25*win[0]*fin[68]*subvol; 
  out[28] += 0.25*win[0]*fin[69]*subvol; 
  out[29] += 0.25*win[0]*fin[72]*subvol; 
  out[30] += 0.25*win[0]*fin[73]*subvol; 
  out[31] += 0.25*win[0]*fin[76]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[23], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[28], 0.3535533905932737*win[0]*fin[28]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[52], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[53], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[56], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[60], 0.3535533905932737*win[0]*fin[60]*subvol); 
  atomicAdd(&out[64], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[66], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[68], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[69], 0.3535533905932737*win[0]*fin[69]*subvol); 
  atomicAdd(&out[72], 0.3535533905932737*win[0]*fin[72]*subvol); 
  atomicAdd(&out[73], 0.3535533905932737*win[0]*fin[73]*subvol); 
  atomicAdd(&out[76], 0.3535533905932737*win[0]*fin[76]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[23] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[24] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[28] += 0.3535533905932737*win[0]*fin[28]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[49] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[50] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[52] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[53] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[56] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[57] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[60] += 0.3535533905932737*win[0]*fin[60]*subvol; 
  out[64] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[65] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[66] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[68] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[69] += 0.3535533905932737*win[0]*fin[69]*subvol; 
  out[72] += 0.3535533905932737*win[0]*fin[72]*subvol; 
  out[73] += 0.3535533905932737*win[0]*fin[73]*subvol; 
  out[76] += 0.3535533905932737*win[0]*fin[76]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[73]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[73]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
  atomicAdd(&out[24], 0.7071067811865475*win[0]*fin[24]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.7071067811865475*win[0]*fin[50]*subvol); 
  atomicAdd(&out[52], 0.7071067811865475*win[0]*fin[52]*subvol); 
  atomicAdd(&out[57], 0.7071067811865475*win[0]*fin[57]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[66], 0.7071067811865475*win[0]*fin[66]*subvol); 
  atomicAdd(&out[68], 0.7071067811865475*win[0]*fin[68]*subvol); 
  atomicAdd(&out[73], 0.7071067811865475*win[0]*fin[73]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[24] += 0.7071067811865475*win[0]*fin[24]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[50] += 0.7071067811865475*win[0]*fin[50]*subvol; 
  out[52] += 0.7071067811865475*win[0]*fin[52]*subvol; 
  out[57] += 0.7071067811865475*win[0]*fin[57]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[66] += 0.7071067811865475*win[0]*fin[66]*subvol; 
  out[68] += 0.7071067811865475*win[0]*fin[68]*subvol; 
  out[73] += 0.7071067811865475*win[0]*fin[73]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[72]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[72]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
  atomicAdd(&out[23], 0.7071067811865475*win[0]*fin[23]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.7071067811865475*win[0]*fin[49]*subvol); 
  atomicAdd(&out[52], 0.7071067811865475*win[0]*fin[52]*subvol); 
  atomicAdd(&out[56], 0.7071067811865475*win[0]*fin[56]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.7071067811865475*win[0]*fin[65]*subvol); 
  atomicAdd(&out[68], 0.7071067811865475*win[0]*fin[68]*subvol); 
  atomicAdd(&out[72], 0.7071067811865475*win[0]*fin[72]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[23] += 0.7071067811865475*win[0]*fin[23]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[49] += 0.7071067811865475*win[0]*fin[49]*subvol; 
  out[52] += 0.7071067811865475*win[0]*fin[52]*subvol; 
  out[56] += 0.7071067811865475*win[0]*fin[56]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[65] += 0.7071067811865475*win[0]*fin[65]*subvol; 
  out[68] += 0.7071067811865475*win[0]*fin[68]*subvol; 
  out[72] += 0.7071067811865475*win[0]*fin[72]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[64]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[52]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[68]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[4]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[15]*subvol; 
  out[4] += 0.5*win[0]*fin[48]*subvol; 
  out[5] += 0.5*win[0]*fin[64]*subvol; 
  out[6] += 0.5*win[0]*fin[52]*subvol; 
  out[7] += 0.5*win[0]*fin[68]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[15], 1.4142135623730951*win[0]*fin[15]*subvol); 
  atomicAdd(&out[48], 1.4142135623730951*win[0]*fin[48]*subvol); 
  atomicAdd(&out[52], 1.4142135623730951*win[0]*fin[52]*subvol); 
  atomicAdd(&out[64], 1.4142135623730951*win[0]*fin[64]*subvol); 
  atomicAdd(&out[68], 1.4142135623730951*win[0]*fin[68]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[15] += 1.4142135623730951*win[0]*fin[15]*subvol; 
  out[48] += 1.4142135623730951*win[0]*fin[48]*subvol; 
  out[52] += 1.4142135623730951*win[0]*fin[52]*subvol; 
  out[64] += 1.4142135623730951*win[0]*fin[64]*subvol; 
  out[68] += 1.4142135623730951*win[0]*fin[68]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[12]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[13]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[14]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[75]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[5]*subvol; 
  out[5] += 0.25*win[0]*fin[6]*subvol; 
  out[6] += 0.25*win[0]*fin[7]*subvol; 
  out[7] += 0.25*win[0]*fin[8]*subvol; 
  out[8] += 0.25*win[0]*fin[12]*subvol; 
  out[9] += 0.25*win[0]*fin[13]*subvol; 
  out[10] += 0.25*win[0]*fin[14]*subvol; 
  out[11] += 0.25*win[0]*fin[16]*subvol; 
  out[12] += 0.25*win[0]*fin[20]*subvol; 
  out[13] += 0.25*win[0]*fin[21]*subvol; 
  out[14] += 0.25*win[0]*fin[22]*subvol; 
  out[15] += 0.25*win[0]*fin[27]*subvol; 
  out[16] += 0.25*win[0]*fin[32]*subvol; 
  out[17] += 0.25*win[0]*fin[33]*subvol; 
  out[18] += 0.25*win[0]*fin[34]*subvol; 
  out[19] += 0.25*win[0]*fin[36]*subvol; 
  out[20] += 0.25*win[0]*fin[37]*subvol; 
  out[21] += 0.25*win[0]*fin[40]*subvol; 
  out[22] += 0.25*win[0]*fin[41]*subvol; 
  out[23] += 0.25*win[0]*fin[44]*subvol; 
  out[24] += 0.25*win[0]*fin[64]*subvol; 
  out[25] += 0.25*win[0]*fin[65]*subvol; 
  out[26] += 0.25*win[0]*fin[66]*subvol; 
  out[27] += 0.25*win[0]*fin[67]*subvol; 
  out[28] += 0.25*win[0]*fin[69]*subvol; 
  out[29] += 0.25*win[0]*fin[70]*subvol; 
  out[30] += 0.25*win[0]*fin[71]*subvol; 
  out[31] += 0.25*win[0]*fin[75]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[27], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[37], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[40], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[41], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[44], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[64], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[66], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[67], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[69], 0.3535533905932737*win[0]*fin[69]*subvol); 
  atomicAdd(&out[70], 0.3535533905932737*win[0]*fin[70]*subvol); 
  atomicAdd(&out[71], 0.3535533905932737*win[0]*fin[71]*subvol); 
  atomicAdd(&out[75], 0.3535533905932737*win[0]*fin[75]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[21] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[22] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[27] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[37] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[40] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[41] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[44] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[64] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[65] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[66] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[67] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[69] += 0.3535533905932737*win[0]*fin[69]*subvol; 
  out[70] += 0.3535533905932737*win[0]*fin[70]*subvol; 
  out[71] += 0.3535533905932737*win[0]*fin[71]*subvol; 
  out[75] += 0.3535533905932737*win[0]*fin[75]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[71]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[71]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[22], 0.7071067811865475*win[0]*fin[22]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[41], 0.7071067811865475*win[0]*fin[41]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[66], 0.7071067811865475*win[0]*fin[66]*subvol); 
  atomicAdd(&out[67], 0.7071067811865475*win[0]*fin[67]*subvol); 
  atomicAdd(&out[71], 0.7071067811865475*win[0]*fin[71]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[22] += 0.7071067811865475*win[0]*fin[22]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[41] += 0.7071067811865475*win[0]*fin[41]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[66] += 0.7071067811865475*win[0]*fin[66]*subvol; 
  out[67] += 0.7071067811865475*win[0]*fin[67]*subvol; 
  out[71] += 0.7071067811865475*win[0]*fin[71]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[70]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[70]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[21], 0.7071067811865475*win[0]*fin[21]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[40], 0.7071067811865475*win[0]*fin[40]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.7071067811865475*win[0]*fin[65]*subvol); 
  atomicAdd(&out[67], 0.7071067811865475*win[0]*fin[67]*subvol); 
  atomicAdd(&out[70], 0.7071067811865475*win[0]*fin[70]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[21] += 0.7071067811865475*win[0]*fin[21]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[40] += 0.7071067811865475*win[0]*fin[40]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[65] += 0.7071067811865475*win[0]*fin[65]*subvol; 
  out[67] += 0.7071067811865475*win[0]*fin[67]*subvol; 
  out[70] += 0.7071067811865475*win[0]*fin[70]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[64]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[36]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[67]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[14]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[64]*subvol; 
  out[6] += 0.5*win[0]*fin[36]*subvol; 
  out[7] += 0.5*win[0]*fin[67]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[14], 1.4142135623730951*win[0]*fin[14]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[36], 1.4142135623730951*win[0]*fin[36]*subvol); 
  atomicAdd(&out[64], 1.4142135623730951*win[0]*fin[64]*subvol); 
  atomicAdd(&out[67], 1.4142135623730951*win[0]*fin[67]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[14] += 1.4142135623730951*win[0]*fin[14]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[36] += 1.4142135623730951*win[0]*fin[36]*subvol; 
  out[64] += 1.4142135623730951*win[0]*fin[64]*subvol; 
  out[67] += 1.4142135623730951*win[0]*fin[67]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[69]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[69]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.7071067811865475*win[0]*fin[65]*subvol); 
  atomicAdd(&out[66], 0.7071067811865475*win[0]*fin[66]*subvol); 
  atomicAdd(&out[69], 0.7071067811865475*win[0]*fin[69]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[65] += 0.7071067811865475*win[0]*fin[65]*subvol; 
  out[66] += 0.7071067811865475*win[0]*fin[66]*subvol; 
  out[69] += 0.7071067811865475*win[0]*fin[69]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[64]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[66]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[13]*subvol; 
  out[4] += 0.5*win[0]*fin[64]*subvol; 
  out[5] += 0.5*win[0]*fin[66]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[13], 1.4142135623730951*win[0]*fin[13]*subvol); 
  atomicAdd(&out[64], 1.4142135623730951*win[0]*fin[64]*subvol); 
  atomicAdd(&out[66], 1.4142135623730951*win[0]*fin[66]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[13] += 1.4142135623730951*win[0]*fin[13]*subvol; 
  out[64] += 1.4142135623730951*win[0]*fin[64]*subvol; 
  out[66] += 1.4142135623730951*win[0]*fin[66]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[64]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[65]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[12]*subvol; 
  out[4] += 0.5*win[0]*fin[64]*subvol; 
  out[5] += 0.5*win[0]*fin[65]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[12], 1.4142135623730951*win[0]*fin[12]*subvol); 
  atomicAdd(&out[64], 1.4142135623730951*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 1.4142135623730951*win[0]*fin[65]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[12] += 1.4142135623730951*win[0]*fin[12]*subvol; 
  out[64] += 1.4142135623730951*win[0]*fin[64]*subvol; 
  out[65] += 1.4142135623730951*win[0]*fin[65]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[64]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[64]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[5], 2.8284271247461907*win[0]*fin[5]*subvol); 
  atomicAdd(&out[64], 2.8284271247461907*win[0]*fin[64]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[5] += 2.8284271247461907*win[0]*fin[5]*subvol; 
  out[64] += 2.8284271247461907*win[0]*fin[64]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.25*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.25*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.25*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.25*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.25*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.25*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.25*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.25*win[0]*fin[8]*subvol); 
  atomicAdd(&out[8], 0.25*win[0]*fin[9]*subvol); 
  atomicAdd(&out[9], 0.25*win[0]*fin[10]*subvol); 
  atomicAdd(&out[10], 0.25*win[0]*fin[11]*subvol); 
  atomicAdd(&out[11], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[59]*subvol); 
#else
  out[0] += 0.25*fin[0]*win[0]*subvol; 
  out[1] += 0.25*win[0]*fin[1]*subvol; 
  out[2] += 0.25*win[0]*fin[2]*subvol; 
  out[3] += 0.25*win[0]*fin[3]*subvol; 
  out[4] += 0.25*win[0]*fin[4]*subvol; 
  out[5] += 0.25*win[0]*fin[6]*subvol; 
  out[6] += 0.25*win[0]*fin[7]*subvol; 
  out[7] += 0.25*win[0]*fin[8]*subvol; 
  out[8] += 0.25*win[0]*fin[9]*subvol; 
  out[9] += 0.25*win[0]*fin[10]*subvol; 
  out[10] += 0.25*win[0]*fin[11]*subvol; 
  out[11] += 0.25*win[0]*fin[16]*subvol; 
  out[12] += 0.25*win[0]*fin[17]*subvol; 
  out[13] += 0.25*win[0]*fin[18]*subvol; 
  out[14] += 0.25*win[0]*fin[19]*subvol; 
  out[15] += 0.25*win[0]*fin[26]*subvol; 
  out[16] += 0.25*win[0]*fin[32]*subvol; 
  out[17] += 0.25*win[0]*fin[33]*subvol; 
  out[18] += 0.25*win[0]*fin[34]*subvol; 
  out[19] += 0.25*win[0]*fin[35]*subvol; 
  out[20] += 0.25*win[0]*fin[37]*subvol; 
  out[21] += 0.25*win[0]*fin[38]*subvol; 
  out[22] += 0.25*win[0]*fin[39]*subvol; 
  out[23] += 0.25*win[0]*fin[43]*subvol; 
  out[24] += 0.25*win[0]*fin[48]*subvol; 
  out[25] += 0.25*win[0]*fin[49]*subvol; 
  out[26] += 0.25*win[0]*fin[50]*subvol; 
  out[27] += 0.25*win[0]*fin[51]*subvol; 
  out[28] += 0.25*win[0]*fin[53]*subvol; 
  out[29] += 0.25*win[0]*fin[54]*subvol; 
  out[30] += 0.25*win[0]*fin[55]*subvol; 
  out[31] += 0.25*win[0]*fin[59]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[26], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[35], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[37], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[38], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[43], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[53], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[54], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[55], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[59], 0.3535533905932737*win[0]*fin[59]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[26] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[35] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[37] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[38] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[39] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[43] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[49] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[50] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[51] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[53] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[54] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[55] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[59] += 0.3535533905932737*win[0]*fin[59]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[55]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[55]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[35], 0.7071067811865475*win[0]*fin[35]*subvol); 
  atomicAdd(&out[39], 0.7071067811865475*win[0]*fin[39]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.7071067811865475*win[0]*fin[50]*subvol); 
  atomicAdd(&out[51], 0.7071067811865475*win[0]*fin[51]*subvol); 
  atomicAdd(&out[55], 0.7071067811865475*win[0]*fin[55]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[35] += 0.7071067811865475*win[0]*fin[35]*subvol; 
  out[39] += 0.7071067811865475*win[0]*fin[39]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[50] += 0.7071067811865475*win[0]*fin[50]*subvol; 
  out[51] += 0.7071067811865475*win[0]*fin[51]*subvol; 
  out[55] += 0.7071067811865475*win[0]*fin[55]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[54]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[54]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 0.7071067811865475*win[0]*fin[35]*subvol); 
  atomicAdd(&out[38], 0.7071067811865475*win[0]*fin[38]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.7071067811865475*win[0]*fin[49]*subvol); 
  atomicAdd(&out[51], 0.7071067811865475*win[0]*fin[51]*subvol); 
  atomicAdd(&out[54], 0.7071067811865475*win[0]*fin[54]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[35] += 0.7071067811865475*win[0]*fin[35]*subvol; 
  out[38] += 0.7071067811865475*win[0]*fin[38]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[49] += 0.7071067811865475*win[0]*fin[49]*subvol; 
  out[51] += 0.7071067811865475*win[0]*fin[51]*subvol; 
  out[54] += 0.7071067811865475*win[0]*fin[54]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[35]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[51]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[11]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[48]*subvol; 
  out[6] += 0.5*win[0]*fin[35]*subvol; 
  out[7] += 0.5*win[0]*fin[51]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[11], 1.4142135623730951*win[0]*fin[11]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[35], 1.4142135623730951*win[0]*fin[35]*subvol); 
  atomicAdd(&out[48], 1.4142135623730951*win[0]*fin[48]*subvol); 
  atomicAdd(&out[51], 1.4142135623730951*win[0]*fin[51]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[11] += 1.4142135623730951*win[0]*fin[11]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[35] += 1.4142135623730951*win[0]*fin[35]*subvol; 
  out[48] += 1.4142135623730951*win[0]*fin[48]*subvol; 
  out[51] += 1.4142135623730951*win[0]*fin[51]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[53]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[53]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.7071067811865475*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.7071067811865475*win[0]*fin[50]*subvol); 
  atomicAdd(&out[53], 0.7071067811865475*win[0]*fin[53]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[49] += 0.7071067811865475*win[0]*fin[49]*subvol; 
  out[50] += 0.7071067811865475*win[0]*fin[50]*subvol; 
  out[53] += 0.7071067811865475*win[0]*fin[53]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[50]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[10]*subvol; 
  out[4] += 0.5*win[0]*fin[48]*subvol; 
  out[5] += 0.5*win[0]*fin[50]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[10], 1.4142135623730951*win[0]*fin[10]*subvol); 
  atomicAdd(&out[48], 1.4142135623730951*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 1.4142135623730951*win[0]*fin[50]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[10] += 1.4142135623730951*win[0]*fin[10]*subvol; 
  out[48] += 1.4142135623730951*win[0]*fin[48]*subvol; 
  out[50] += 1.4142135623730951*win[0]*fin[50]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[48]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[49]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[9]*subvol; 
  out[4] += 0.5*win[0]*fin[48]*subvol; 
  out[5] += 0.5*win[0]*fin[49]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[9], 1.4142135623730951*win[0]*fin[9]*subvol); 
  atomicAdd(&out[48], 1.4142135623730951*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 1.4142135623730951*win[0]*fin[49]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[9] += 1.4142135623730951*win[0]*fin[9]*subvol; 
  out[48] += 1.4142135623730951*win[0]*fin[48]*subvol; 
  out[49] += 1.4142135623730951*win[0]*fin[49]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[48]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[48]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.8284271247461907*win[0]*fin[4]*subvol); 
  atomicAdd(&out[48], 2.8284271247461907*win[0]*fin[48]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[4] += 2.8284271247461907*win[0]*fin[4]*subvol; 
  out[48] += 2.8284271247461907*win[0]*fin[48]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[37]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[37]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[37], 0.7071067811865475*win[0]*fin[37]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[37] += 0.7071067811865475*win[0]*fin[37]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[34]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[8]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[34]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[8], 1.4142135623730951*win[0]*fin[8]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 1.4142135623730951*win[0]*fin[34]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[8] += 1.4142135623730951*win[0]*fin[8]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[34] += 1.4142135623730951*win[0]*fin[34]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[33]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[7]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[33]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[7], 1.4142135623730951*win[0]*fin[7]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 1.4142135623730951*win[0]*fin[33]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[7] += 1.4142135623730951*win[0]*fin[7]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[33] += 1.4142135623730951*win[0]*fin[33]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[32]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[32]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.8284271247461907*win[0]*fin[3]*subvol); 
  atomicAdd(&out[32], 2.8284271247461907*win[0]*fin[32]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[3] += 2.8284271247461907*win[0]*fin[3]*subvol; 
  out[32] += 2.8284271247461907*win[0]*fin[32]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[6]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[6]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[6], 1.4142135623730951*win[0]*fin[6]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[6] += 1.4142135623730951*win[0]*fin[6]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.8284271247461907*win[0]*fin[2]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[2] += 2.8284271247461907*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgyvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.8284271247461907*win[0]*fin[1]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[1] += 2.8284271247461907*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_hyb_p1_avgxyvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 5.656854249492382*fin[0]*win[0]*subvol); 
#else
  out[0] += 5.656854249492382*fin[0]*win[0]*subvol; 
#endif
} 

