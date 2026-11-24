#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[11], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[55]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[62]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[97]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[99]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[102]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[106]*subvol); 
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
  out[11] += 0.25*win[0]*fin[17]*subvol; 
  out[12] += 0.25*win[0]*fin[18]*subvol; 
  out[13] += 0.25*win[0]*fin[19]*subvol; 
  out[14] += 0.25*win[0]*fin[20]*subvol; 
  out[15] += 0.25*win[0]*fin[24]*subvol; 
  out[16] += 0.25*win[0]*fin[27]*subvol; 
  out[17] += 0.25*win[0]*fin[29]*subvol; 
  out[18] += 0.25*win[0]*fin[30]*subvol; 
  out[19] += 0.25*win[0]*fin[34]*subvol; 
  out[20] += 0.25*win[0]*fin[36]*subvol; 
  out[21] += 0.25*win[0]*fin[38]*subvol; 
  out[22] += 0.25*win[0]*fin[39]*subvol; 
  out[23] += 0.25*win[0]*fin[41]*subvol; 
  out[24] += 0.25*win[0]*fin[42]*subvol; 
  out[25] += 0.25*win[0]*fin[44]*subvol; 
  out[26] += 0.25*win[0]*fin[45]*subvol; 
  out[27] += 0.25*win[0]*fin[46]*subvol; 
  out[28] += 0.25*win[0]*fin[48]*subvol; 
  out[29] += 0.25*win[0]*fin[49]*subvol; 
  out[30] += 0.25*win[0]*fin[50]*subvol; 
  out[31] += 0.25*win[0]*fin[55]*subvol; 
  out[32] += 0.25*win[0]*fin[62]*subvol; 
  out[33] += 0.25*win[0]*fin[64]*subvol; 
  out[34] += 0.25*win[0]*fin[67]*subvol; 
  out[35] += 0.25*win[0]*fin[71]*subvol; 
  out[36] += 0.25*win[0]*fin[73]*subvol; 
  out[37] += 0.25*win[0]*fin[75]*subvol; 
  out[38] += 0.25*win[0]*fin[76]*subvol; 
  out[39] += 0.25*win[0]*fin[78]*subvol; 
  out[40] += 0.25*win[0]*fin[79]*subvol; 
  out[41] += 0.25*win[0]*fin[82]*subvol; 
  out[42] += 0.25*win[0]*fin[84]*subvol; 
  out[43] += 0.25*win[0]*fin[85]*subvol; 
  out[44] += 0.25*win[0]*fin[97]*subvol; 
  out[45] += 0.25*win[0]*fin[99]*subvol; 
  out[46] += 0.25*win[0]*fin[102]*subvol; 
  out[47] += 0.25*win[0]*fin[106]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[24], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[27], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[29], 0.3535533905932737*win[0]*fin[29]*subvol); 
  atomicAdd(&out[30], 0.3535533905932737*win[0]*fin[30]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[38], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[41], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[44], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[45], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[55], 0.3535533905932737*win[0]*fin[55]*subvol); 
  atomicAdd(&out[62], 0.3535533905932737*win[0]*fin[62]*subvol); 
  atomicAdd(&out[64], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[67], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[71], 0.3535533905932737*win[0]*fin[71]*subvol); 
  atomicAdd(&out[73], 0.3535533905932737*win[0]*fin[73]*subvol); 
  atomicAdd(&out[75], 0.3535533905932737*win[0]*fin[75]*subvol); 
  atomicAdd(&out[76], 0.3535533905932737*win[0]*fin[76]*subvol); 
  atomicAdd(&out[78], 0.3535533905932737*win[0]*fin[78]*subvol); 
  atomicAdd(&out[79], 0.3535533905932737*win[0]*fin[79]*subvol); 
  atomicAdd(&out[82], 0.3535533905932737*win[0]*fin[82]*subvol); 
  atomicAdd(&out[84], 0.3535533905932737*win[0]*fin[84]*subvol); 
  atomicAdd(&out[85], 0.3535533905932737*win[0]*fin[85]*subvol); 
  atomicAdd(&out[97], 0.3535533905932737*win[0]*fin[97]*subvol); 
  atomicAdd(&out[99], 0.3535533905932737*win[0]*fin[99]*subvol); 
  atomicAdd(&out[102], 0.3535533905932737*win[0]*fin[102]*subvol); 
  atomicAdd(&out[106], 0.3535533905932737*win[0]*fin[106]*subvol); 
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
  out[17] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[24] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[27] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[29] += 0.3535533905932737*win[0]*fin[29]*subvol; 
  out[30] += 0.3535533905932737*win[0]*fin[30]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[38] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[39] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[41] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[42] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[44] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[45] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[46] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[49] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[50] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[55] += 0.3535533905932737*win[0]*fin[55]*subvol; 
  out[62] += 0.3535533905932737*win[0]*fin[62]*subvol; 
  out[64] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[67] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[71] += 0.3535533905932737*win[0]*fin[71]*subvol; 
  out[73] += 0.3535533905932737*win[0]*fin[73]*subvol; 
  out[75] += 0.3535533905932737*win[0]*fin[75]*subvol; 
  out[76] += 0.3535533905932737*win[0]*fin[76]*subvol; 
  out[78] += 0.3535533905932737*win[0]*fin[78]*subvol; 
  out[79] += 0.3535533905932737*win[0]*fin[79]*subvol; 
  out[82] += 0.3535533905932737*win[0]*fin[82]*subvol; 
  out[84] += 0.3535533905932737*win[0]*fin[84]*subvol; 
  out[85] += 0.3535533905932737*win[0]*fin[85]*subvol; 
  out[97] += 0.3535533905932737*win[0]*fin[97]*subvol; 
  out[99] += 0.3535533905932737*win[0]*fin[99]*subvol; 
  out[102] += 0.3535533905932737*win[0]*fin[102]*subvol; 
  out[106] += 0.3535533905932737*win[0]*fin[106]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[11], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[30]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[54]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[61]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[63]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[76]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[79]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[85]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[96]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[98]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[101]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[105]*subvol); 
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
  out[11] += 0.25*win[0]*fin[16]*subvol; 
  out[12] += 0.25*win[0]*fin[18]*subvol; 
  out[13] += 0.25*win[0]*fin[19]*subvol; 
  out[14] += 0.25*win[0]*fin[20]*subvol; 
  out[15] += 0.25*win[0]*fin[23]*subvol; 
  out[16] += 0.25*win[0]*fin[26]*subvol; 
  out[17] += 0.25*win[0]*fin[28]*subvol; 
  out[18] += 0.25*win[0]*fin[30]*subvol; 
  out[19] += 0.25*win[0]*fin[33]*subvol; 
  out[20] += 0.25*win[0]*fin[35]*subvol; 
  out[21] += 0.25*win[0]*fin[37]*subvol; 
  out[22] += 0.25*win[0]*fin[39]*subvol; 
  out[23] += 0.25*win[0]*fin[40]*subvol; 
  out[24] += 0.25*win[0]*fin[42]*subvol; 
  out[25] += 0.25*win[0]*fin[43]*subvol; 
  out[26] += 0.25*win[0]*fin[45]*subvol; 
  out[27] += 0.25*win[0]*fin[46]*subvol; 
  out[28] += 0.25*win[0]*fin[47]*subvol; 
  out[29] += 0.25*win[0]*fin[49]*subvol; 
  out[30] += 0.25*win[0]*fin[50]*subvol; 
  out[31] += 0.25*win[0]*fin[54]*subvol; 
  out[32] += 0.25*win[0]*fin[61]*subvol; 
  out[33] += 0.25*win[0]*fin[63]*subvol; 
  out[34] += 0.25*win[0]*fin[66]*subvol; 
  out[35] += 0.25*win[0]*fin[70]*subvol; 
  out[36] += 0.25*win[0]*fin[72]*subvol; 
  out[37] += 0.25*win[0]*fin[74]*subvol; 
  out[38] += 0.25*win[0]*fin[76]*subvol; 
  out[39] += 0.25*win[0]*fin[77]*subvol; 
  out[40] += 0.25*win[0]*fin[79]*subvol; 
  out[41] += 0.25*win[0]*fin[81]*subvol; 
  out[42] += 0.25*win[0]*fin[83]*subvol; 
  out[43] += 0.25*win[0]*fin[85]*subvol; 
  out[44] += 0.25*win[0]*fin[96]*subvol; 
  out[45] += 0.25*win[0]*fin[98]*subvol; 
  out[46] += 0.25*win[0]*fin[101]*subvol; 
  out[47] += 0.25*win[0]*fin[105]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[23], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[26], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[28], 0.3535533905932737*win[0]*fin[28]*subvol); 
  atomicAdd(&out[30], 0.3535533905932737*win[0]*fin[30]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[37], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[39], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[40], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[42], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[43], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[45], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[49], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[54], 0.3535533905932737*win[0]*fin[54]*subvol); 
  atomicAdd(&out[61], 0.3535533905932737*win[0]*fin[61]*subvol); 
  atomicAdd(&out[63], 0.3535533905932737*win[0]*fin[63]*subvol); 
  atomicAdd(&out[66], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[70], 0.3535533905932737*win[0]*fin[70]*subvol); 
  atomicAdd(&out[72], 0.3535533905932737*win[0]*fin[72]*subvol); 
  atomicAdd(&out[74], 0.3535533905932737*win[0]*fin[74]*subvol); 
  atomicAdd(&out[76], 0.3535533905932737*win[0]*fin[76]*subvol); 
  atomicAdd(&out[77], 0.3535533905932737*win[0]*fin[77]*subvol); 
  atomicAdd(&out[79], 0.3535533905932737*win[0]*fin[79]*subvol); 
  atomicAdd(&out[81], 0.3535533905932737*win[0]*fin[81]*subvol); 
  atomicAdd(&out[83], 0.3535533905932737*win[0]*fin[83]*subvol); 
  atomicAdd(&out[85], 0.3535533905932737*win[0]*fin[85]*subvol); 
  atomicAdd(&out[96], 0.3535533905932737*win[0]*fin[96]*subvol); 
  atomicAdd(&out[98], 0.3535533905932737*win[0]*fin[98]*subvol); 
  atomicAdd(&out[101], 0.3535533905932737*win[0]*fin[101]*subvol); 
  atomicAdd(&out[105], 0.3535533905932737*win[0]*fin[105]*subvol); 
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
  out[16] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[23] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[26] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[28] += 0.3535533905932737*win[0]*fin[28]*subvol; 
  out[30] += 0.3535533905932737*win[0]*fin[30]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[35] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[37] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[39] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[40] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[42] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[43] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[45] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[46] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[47] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[49] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[50] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[54] += 0.3535533905932737*win[0]*fin[54]*subvol; 
  out[61] += 0.3535533905932737*win[0]*fin[61]*subvol; 
  out[63] += 0.3535533905932737*win[0]*fin[63]*subvol; 
  out[66] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[70] += 0.3535533905932737*win[0]*fin[70]*subvol; 
  out[72] += 0.3535533905932737*win[0]*fin[72]*subvol; 
  out[74] += 0.3535533905932737*win[0]*fin[74]*subvol; 
  out[76] += 0.3535533905932737*win[0]*fin[76]*subvol; 
  out[77] += 0.3535533905932737*win[0]*fin[77]*subvol; 
  out[79] += 0.3535533905932737*win[0]*fin[79]*subvol; 
  out[81] += 0.3535533905932737*win[0]*fin[81]*subvol; 
  out[83] += 0.3535533905932737*win[0]*fin[83]*subvol; 
  out[85] += 0.3535533905932737*win[0]*fin[85]*subvol; 
  out[96] += 0.3535533905932737*win[0]*fin[96]*subvol; 
  out[98] += 0.3535533905932737*win[0]*fin[98]*subvol; 
  out[101] += 0.3535533905932737*win[0]*fin[101]*subvol; 
  out[105] += 0.3535533905932737*win[0]*fin[105]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[30]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[76]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[79]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[85]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[30]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[76]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[79]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[85]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[30], 0.7071067811865475*win[0]*fin[30]*subvol); 
  atomicAdd(&out[39], 0.7071067811865475*win[0]*fin[39]*subvol); 
  atomicAdd(&out[42], 0.7071067811865475*win[0]*fin[42]*subvol); 
  atomicAdd(&out[45], 0.7071067811865475*win[0]*fin[45]*subvol); 
  atomicAdd(&out[46], 0.7071067811865475*win[0]*fin[46]*subvol); 
  atomicAdd(&out[49], 0.7071067811865475*win[0]*fin[49]*subvol); 
  atomicAdd(&out[50], 0.7071067811865475*win[0]*fin[50]*subvol); 
  atomicAdd(&out[76], 0.7071067811865475*win[0]*fin[76]*subvol); 
  atomicAdd(&out[79], 0.7071067811865475*win[0]*fin[79]*subvol); 
  atomicAdd(&out[85], 0.7071067811865475*win[0]*fin[85]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[30] += 0.7071067811865475*win[0]*fin[30]*subvol; 
  out[39] += 0.7071067811865475*win[0]*fin[39]*subvol; 
  out[42] += 0.7071067811865475*win[0]*fin[42]*subvol; 
  out[45] += 0.7071067811865475*win[0]*fin[45]*subvol; 
  out[46] += 0.7071067811865475*win[0]*fin[46]*subvol; 
  out[49] += 0.7071067811865475*win[0]*fin[49]*subvol; 
  out[50] += 0.7071067811865475*win[0]*fin[50]*subvol; 
  out[76] += 0.7071067811865475*win[0]*fin[76]*subvol; 
  out[79] += 0.7071067811865475*win[0]*fin[79]*subvol; 
  out[85] += 0.7071067811865475*win[0]*fin[85]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[11], 0.25*win[0]*fin[16]*subvol); 
  atomicAdd(&out[12], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[19]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[28]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[29]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[46]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[50]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[53]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[60]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[74]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[75]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[77]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[78]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[83]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[84]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[94]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[95]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[100]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[104]*subvol); 
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
  out[11] += 0.25*win[0]*fin[16]*subvol; 
  out[12] += 0.25*win[0]*fin[17]*subvol; 
  out[13] += 0.25*win[0]*fin[19]*subvol; 
  out[14] += 0.25*win[0]*fin[20]*subvol; 
  out[15] += 0.25*win[0]*fin[22]*subvol; 
  out[16] += 0.25*win[0]*fin[25]*subvol; 
  out[17] += 0.25*win[0]*fin[28]*subvol; 
  out[18] += 0.25*win[0]*fin[29]*subvol; 
  out[19] += 0.25*win[0]*fin[31]*subvol; 
  out[20] += 0.25*win[0]*fin[32]*subvol; 
  out[21] += 0.25*win[0]*fin[37]*subvol; 
  out[22] += 0.25*win[0]*fin[38]*subvol; 
  out[23] += 0.25*win[0]*fin[40]*subvol; 
  out[24] += 0.25*win[0]*fin[41]*subvol; 
  out[25] += 0.25*win[0]*fin[43]*subvol; 
  out[26] += 0.25*win[0]*fin[44]*subvol; 
  out[27] += 0.25*win[0]*fin[46]*subvol; 
  out[28] += 0.25*win[0]*fin[47]*subvol; 
  out[29] += 0.25*win[0]*fin[48]*subvol; 
  out[30] += 0.25*win[0]*fin[50]*subvol; 
  out[31] += 0.25*win[0]*fin[53]*subvol; 
  out[32] += 0.25*win[0]*fin[59]*subvol; 
  out[33] += 0.25*win[0]*fin[60]*subvol; 
  out[34] += 0.25*win[0]*fin[65]*subvol; 
  out[35] += 0.25*win[0]*fin[68]*subvol; 
  out[36] += 0.25*win[0]*fin[69]*subvol; 
  out[37] += 0.25*win[0]*fin[74]*subvol; 
  out[38] += 0.25*win[0]*fin[75]*subvol; 
  out[39] += 0.25*win[0]*fin[77]*subvol; 
  out[40] += 0.25*win[0]*fin[78]*subvol; 
  out[41] += 0.25*win[0]*fin[80]*subvol; 
  out[42] += 0.25*win[0]*fin[83]*subvol; 
  out[43] += 0.25*win[0]*fin[84]*subvol; 
  out[44] += 0.25*win[0]*fin[94]*subvol; 
  out[45] += 0.25*win[0]*fin[95]*subvol; 
  out[46] += 0.25*win[0]*fin[100]*subvol; 
  out[47] += 0.25*win[0]*fin[104]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[22], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[25], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[28], 0.3535533905932737*win[0]*fin[28]*subvol); 
  atomicAdd(&out[29], 0.3535533905932737*win[0]*fin[29]*subvol); 
  atomicAdd(&out[31], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[37], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[38], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[40], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[41], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[43], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[44], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[46], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[53], 0.3535533905932737*win[0]*fin[53]*subvol); 
  atomicAdd(&out[59], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[60], 0.3535533905932737*win[0]*fin[60]*subvol); 
  atomicAdd(&out[65], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[68], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[69], 0.3535533905932737*win[0]*fin[69]*subvol); 
  atomicAdd(&out[74], 0.3535533905932737*win[0]*fin[74]*subvol); 
  atomicAdd(&out[75], 0.3535533905932737*win[0]*fin[75]*subvol); 
  atomicAdd(&out[77], 0.3535533905932737*win[0]*fin[77]*subvol); 
  atomicAdd(&out[78], 0.3535533905932737*win[0]*fin[78]*subvol); 
  atomicAdd(&out[80], 0.3535533905932737*win[0]*fin[80]*subvol); 
  atomicAdd(&out[83], 0.3535533905932737*win[0]*fin[83]*subvol); 
  atomicAdd(&out[84], 0.3535533905932737*win[0]*fin[84]*subvol); 
  atomicAdd(&out[94], 0.3535533905932737*win[0]*fin[94]*subvol); 
  atomicAdd(&out[95], 0.3535533905932737*win[0]*fin[95]*subvol); 
  atomicAdd(&out[100], 0.3535533905932737*win[0]*fin[100]*subvol); 
  atomicAdd(&out[104], 0.3535533905932737*win[0]*fin[104]*subvol); 
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
  out[16] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[22] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[25] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[28] += 0.3535533905932737*win[0]*fin[28]*subvol; 
  out[29] += 0.3535533905932737*win[0]*fin[29]*subvol; 
  out[31] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[37] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[38] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[40] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[41] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[43] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[44] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[46] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[47] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[50] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[53] += 0.3535533905932737*win[0]*fin[53]*subvol; 
  out[59] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[60] += 0.3535533905932737*win[0]*fin[60]*subvol; 
  out[65] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[68] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[69] += 0.3535533905932737*win[0]*fin[69]*subvol; 
  out[74] += 0.3535533905932737*win[0]*fin[74]*subvol; 
  out[75] += 0.3535533905932737*win[0]*fin[75]*subvol; 
  out[77] += 0.3535533905932737*win[0]*fin[77]*subvol; 
  out[78] += 0.3535533905932737*win[0]*fin[78]*subvol; 
  out[80] += 0.3535533905932737*win[0]*fin[80]*subvol; 
  out[83] += 0.3535533905932737*win[0]*fin[83]*subvol; 
  out[84] += 0.3535533905932737*win[0]*fin[84]*subvol; 
  out[94] += 0.3535533905932737*win[0]*fin[94]*subvol; 
  out[95] += 0.3535533905932737*win[0]*fin[95]*subvol; 
  out[100] += 0.3535533905932737*win[0]*fin[100]*subvol; 
  out[104] += 0.3535533905932737*win[0]*fin[104]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[29]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[75]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[78]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[84]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[29]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[75]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[78]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[84]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[29], 0.7071067811865475*win[0]*fin[29]*subvol); 
  atomicAdd(&out[38], 0.7071067811865475*win[0]*fin[38]*subvol); 
  atomicAdd(&out[41], 0.7071067811865475*win[0]*fin[41]*subvol); 
  atomicAdd(&out[44], 0.7071067811865475*win[0]*fin[44]*subvol); 
  atomicAdd(&out[46], 0.7071067811865475*win[0]*fin[46]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[50], 0.7071067811865475*win[0]*fin[50]*subvol); 
  atomicAdd(&out[75], 0.7071067811865475*win[0]*fin[75]*subvol); 
  atomicAdd(&out[78], 0.7071067811865475*win[0]*fin[78]*subvol); 
  atomicAdd(&out[84], 0.7071067811865475*win[0]*fin[84]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[29] += 0.7071067811865475*win[0]*fin[29]*subvol; 
  out[38] += 0.7071067811865475*win[0]*fin[38]*subvol; 
  out[41] += 0.7071067811865475*win[0]*fin[41]*subvol; 
  out[44] += 0.7071067811865475*win[0]*fin[44]*subvol; 
  out[46] += 0.7071067811865475*win[0]*fin[46]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[50] += 0.7071067811865475*win[0]*fin[50]*subvol; 
  out[75] += 0.7071067811865475*win[0]*fin[75]*subvol; 
  out[78] += 0.7071067811865475*win[0]*fin[78]*subvol; 
  out[84] += 0.7071067811865475*win[0]*fin[84]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[28]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[46]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[50]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[74]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[77]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[83]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[28]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[46]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[50]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[74]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[77]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[83]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[28], 0.7071067811865475*win[0]*fin[28]*subvol); 
  atomicAdd(&out[37], 0.7071067811865475*win[0]*fin[37]*subvol); 
  atomicAdd(&out[40], 0.7071067811865475*win[0]*fin[40]*subvol); 
  atomicAdd(&out[43], 0.7071067811865475*win[0]*fin[43]*subvol); 
  atomicAdd(&out[46], 0.7071067811865475*win[0]*fin[46]*subvol); 
  atomicAdd(&out[47], 0.7071067811865475*win[0]*fin[47]*subvol); 
  atomicAdd(&out[50], 0.7071067811865475*win[0]*fin[50]*subvol); 
  atomicAdd(&out[74], 0.7071067811865475*win[0]*fin[74]*subvol); 
  atomicAdd(&out[77], 0.7071067811865475*win[0]*fin[77]*subvol); 
  atomicAdd(&out[83], 0.7071067811865475*win[0]*fin[83]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[28] += 0.7071067811865475*win[0]*fin[28]*subvol; 
  out[37] += 0.7071067811865475*win[0]*fin[37]*subvol; 
  out[40] += 0.7071067811865475*win[0]*fin[40]*subvol; 
  out[43] += 0.7071067811865475*win[0]*fin[43]*subvol; 
  out[46] += 0.7071067811865475*win[0]*fin[46]*subvol; 
  out[47] += 0.7071067811865475*win[0]*fin[47]*subvol; 
  out[50] += 0.7071067811865475*win[0]*fin[50]*subvol; 
  out[74] += 0.7071067811865475*win[0]*fin[74]*subvol; 
  out[77] += 0.7071067811865475*win[0]*fin[77]*subvol; 
  out[83] += 0.7071067811865475*win[0]*fin[83]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[46]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[50]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[4]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[15]*subvol; 
  out[4] += 0.5*win[0]*fin[19]*subvol; 
  out[5] += 0.5*win[0]*fin[20]*subvol; 
  out[6] += 0.5*win[0]*fin[46]*subvol; 
  out[7] += 0.5*win[0]*fin[50]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[15], 1.4142135623730951*win[0]*fin[15]*subvol); 
  atomicAdd(&out[19], 1.4142135623730951*win[0]*fin[19]*subvol); 
  atomicAdd(&out[20], 1.4142135623730951*win[0]*fin[20]*subvol); 
  atomicAdd(&out[46], 1.4142135623730951*win[0]*fin[46]*subvol); 
  atomicAdd(&out[50], 1.4142135623730951*win[0]*fin[50]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[15] += 1.4142135623730951*win[0]*fin[15]*subvol; 
  out[19] += 1.4142135623730951*win[0]*fin[19]*subvol; 
  out[20] += 1.4142135623730951*win[0]*fin[20]*subvol; 
  out[46] += 1.4142135623730951*win[0]*fin[46]*subvol; 
  out[50] += 1.4142135623730951*win[0]*fin[50]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[12], 0.25*win[0]*fin[17]*subvol); 
  atomicAdd(&out[13], 0.25*win[0]*fin[18]*subvol); 
  atomicAdd(&out[14], 0.25*win[0]*fin[20]*subvol); 
  atomicAdd(&out[15], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[25]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[26]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[27]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[43]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[44]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[45]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[47]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[48]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[49]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[52]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[68]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[69]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[70]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[71]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[72]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[73]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[80]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[81]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[82]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[91]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[92]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[93]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[103]*subvol); 
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
  out[12] += 0.25*win[0]*fin[17]*subvol; 
  out[13] += 0.25*win[0]*fin[18]*subvol; 
  out[14] += 0.25*win[0]*fin[20]*subvol; 
  out[15] += 0.25*win[0]*fin[21]*subvol; 
  out[16] += 0.25*win[0]*fin[25]*subvol; 
  out[17] += 0.25*win[0]*fin[26]*subvol; 
  out[18] += 0.25*win[0]*fin[27]*subvol; 
  out[19] += 0.25*win[0]*fin[31]*subvol; 
  out[20] += 0.25*win[0]*fin[32]*subvol; 
  out[21] += 0.25*win[0]*fin[33]*subvol; 
  out[22] += 0.25*win[0]*fin[34]*subvol; 
  out[23] += 0.25*win[0]*fin[35]*subvol; 
  out[24] += 0.25*win[0]*fin[36]*subvol; 
  out[25] += 0.25*win[0]*fin[43]*subvol; 
  out[26] += 0.25*win[0]*fin[44]*subvol; 
  out[27] += 0.25*win[0]*fin[45]*subvol; 
  out[28] += 0.25*win[0]*fin[47]*subvol; 
  out[29] += 0.25*win[0]*fin[48]*subvol; 
  out[30] += 0.25*win[0]*fin[49]*subvol; 
  out[31] += 0.25*win[0]*fin[52]*subvol; 
  out[32] += 0.25*win[0]*fin[56]*subvol; 
  out[33] += 0.25*win[0]*fin[57]*subvol; 
  out[34] += 0.25*win[0]*fin[58]*subvol; 
  out[35] += 0.25*win[0]*fin[68]*subvol; 
  out[36] += 0.25*win[0]*fin[69]*subvol; 
  out[37] += 0.25*win[0]*fin[70]*subvol; 
  out[38] += 0.25*win[0]*fin[71]*subvol; 
  out[39] += 0.25*win[0]*fin[72]*subvol; 
  out[40] += 0.25*win[0]*fin[73]*subvol; 
  out[41] += 0.25*win[0]*fin[80]*subvol; 
  out[42] += 0.25*win[0]*fin[81]*subvol; 
  out[43] += 0.25*win[0]*fin[82]*subvol; 
  out[44] += 0.25*win[0]*fin[91]*subvol; 
  out[45] += 0.25*win[0]*fin[92]*subvol; 
  out[46] += 0.25*win[0]*fin[93]*subvol; 
  out[47] += 0.25*win[0]*fin[103]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[21], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[25], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[26], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[27], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[31], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[35], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[43], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[44], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[45], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[47], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[52], 0.3535533905932737*win[0]*fin[52]*subvol); 
  atomicAdd(&out[56], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[68], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[69], 0.3535533905932737*win[0]*fin[69]*subvol); 
  atomicAdd(&out[70], 0.3535533905932737*win[0]*fin[70]*subvol); 
  atomicAdd(&out[71], 0.3535533905932737*win[0]*fin[71]*subvol); 
  atomicAdd(&out[72], 0.3535533905932737*win[0]*fin[72]*subvol); 
  atomicAdd(&out[73], 0.3535533905932737*win[0]*fin[73]*subvol); 
  atomicAdd(&out[80], 0.3535533905932737*win[0]*fin[80]*subvol); 
  atomicAdd(&out[81], 0.3535533905932737*win[0]*fin[81]*subvol); 
  atomicAdd(&out[82], 0.3535533905932737*win[0]*fin[82]*subvol); 
  atomicAdd(&out[91], 0.3535533905932737*win[0]*fin[91]*subvol); 
  atomicAdd(&out[92], 0.3535533905932737*win[0]*fin[92]*subvol); 
  atomicAdd(&out[93], 0.3535533905932737*win[0]*fin[93]*subvol); 
  atomicAdd(&out[103], 0.3535533905932737*win[0]*fin[103]*subvol); 
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
  out[17] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[20] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[21] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[25] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[26] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[27] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[31] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[35] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[43] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[44] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[45] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[47] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[48] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[49] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[52] += 0.3535533905932737*win[0]*fin[52]*subvol; 
  out[56] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[57] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[58] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[68] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[69] += 0.3535533905932737*win[0]*fin[69]*subvol; 
  out[70] += 0.3535533905932737*win[0]*fin[70]*subvol; 
  out[71] += 0.3535533905932737*win[0]*fin[71]*subvol; 
  out[72] += 0.3535533905932737*win[0]*fin[72]*subvol; 
  out[73] += 0.3535533905932737*win[0]*fin[73]*subvol; 
  out[80] += 0.3535533905932737*win[0]*fin[80]*subvol; 
  out[81] += 0.3535533905932737*win[0]*fin[81]*subvol; 
  out[82] += 0.3535533905932737*win[0]*fin[82]*subvol; 
  out[91] += 0.3535533905932737*win[0]*fin[91]*subvol; 
  out[92] += 0.3535533905932737*win[0]*fin[92]*subvol; 
  out[93] += 0.3535533905932737*win[0]*fin[93]*subvol; 
  out[103] += 0.3535533905932737*win[0]*fin[103]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[71]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[73]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[82]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[71]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[73]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[82]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[27], 0.7071067811865475*win[0]*fin[27]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[44], 0.7071067811865475*win[0]*fin[44]*subvol); 
  atomicAdd(&out[45], 0.7071067811865475*win[0]*fin[45]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[49], 0.7071067811865475*win[0]*fin[49]*subvol); 
  atomicAdd(&out[71], 0.7071067811865475*win[0]*fin[71]*subvol); 
  atomicAdd(&out[73], 0.7071067811865475*win[0]*fin[73]*subvol); 
  atomicAdd(&out[82], 0.7071067811865475*win[0]*fin[82]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[27] += 0.7071067811865475*win[0]*fin[27]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[44] += 0.7071067811865475*win[0]*fin[44]*subvol; 
  out[45] += 0.7071067811865475*win[0]*fin[45]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[49] += 0.7071067811865475*win[0]*fin[49]*subvol; 
  out[71] += 0.7071067811865475*win[0]*fin[71]*subvol; 
  out[73] += 0.7071067811865475*win[0]*fin[73]*subvol; 
  out[82] += 0.7071067811865475*win[0]*fin[82]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[45]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[49]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[70]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[72]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[81]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[45]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[49]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[70]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[72]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[81]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[26], 0.7071067811865475*win[0]*fin[26]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 0.7071067811865475*win[0]*fin[35]*subvol); 
  atomicAdd(&out[43], 0.7071067811865475*win[0]*fin[43]*subvol); 
  atomicAdd(&out[45], 0.7071067811865475*win[0]*fin[45]*subvol); 
  atomicAdd(&out[47], 0.7071067811865475*win[0]*fin[47]*subvol); 
  atomicAdd(&out[49], 0.7071067811865475*win[0]*fin[49]*subvol); 
  atomicAdd(&out[70], 0.7071067811865475*win[0]*fin[70]*subvol); 
  atomicAdd(&out[72], 0.7071067811865475*win[0]*fin[72]*subvol); 
  atomicAdd(&out[81], 0.7071067811865475*win[0]*fin[81]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[26] += 0.7071067811865475*win[0]*fin[26]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[35] += 0.7071067811865475*win[0]*fin[35]*subvol; 
  out[43] += 0.7071067811865475*win[0]*fin[43]*subvol; 
  out[45] += 0.7071067811865475*win[0]*fin[45]*subvol; 
  out[47] += 0.7071067811865475*win[0]*fin[47]*subvol; 
  out[49] += 0.7071067811865475*win[0]*fin[49]*subvol; 
  out[70] += 0.7071067811865475*win[0]*fin[70]*subvol; 
  out[72] += 0.7071067811865475*win[0]*fin[72]*subvol; 
  out[81] += 0.7071067811865475*win[0]*fin[81]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[45]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[49]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[14]*subvol; 
  out[4] += 0.5*win[0]*fin[18]*subvol; 
  out[5] += 0.5*win[0]*fin[20]*subvol; 
  out[6] += 0.5*win[0]*fin[45]*subvol; 
  out[7] += 0.5*win[0]*fin[49]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[14], 1.4142135623730951*win[0]*fin[14]*subvol); 
  atomicAdd(&out[18], 1.4142135623730951*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 1.4142135623730951*win[0]*fin[20]*subvol); 
  atomicAdd(&out[45], 1.4142135623730951*win[0]*fin[45]*subvol); 
  atomicAdd(&out[49], 1.4142135623730951*win[0]*fin[49]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[14] += 1.4142135623730951*win[0]*fin[14]*subvol; 
  out[18] += 1.4142135623730951*win[0]*fin[18]*subvol; 
  out[20] += 1.4142135623730951*win[0]*fin[20]*subvol; 
  out[45] += 1.4142135623730951*win[0]*fin[45]*subvol; 
  out[49] += 1.4142135623730951*win[0]*fin[49]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[20]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[43]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[44]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[47]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[48]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[68]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[69]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[80]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[20]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[43]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[44]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[47]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[48]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[68]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[69]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[80]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[20], 0.7071067811865475*win[0]*fin[20]*subvol); 
  atomicAdd(&out[25], 0.7071067811865475*win[0]*fin[25]*subvol); 
  atomicAdd(&out[31], 0.7071067811865475*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[43], 0.7071067811865475*win[0]*fin[43]*subvol); 
  atomicAdd(&out[44], 0.7071067811865475*win[0]*fin[44]*subvol); 
  atomicAdd(&out[47], 0.7071067811865475*win[0]*fin[47]*subvol); 
  atomicAdd(&out[48], 0.7071067811865475*win[0]*fin[48]*subvol); 
  atomicAdd(&out[68], 0.7071067811865475*win[0]*fin[68]*subvol); 
  atomicAdd(&out[69], 0.7071067811865475*win[0]*fin[69]*subvol); 
  atomicAdd(&out[80], 0.7071067811865475*win[0]*fin[80]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
  out[25] += 0.7071067811865475*win[0]*fin[25]*subvol; 
  out[31] += 0.7071067811865475*win[0]*fin[31]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[43] += 0.7071067811865475*win[0]*fin[43]*subvol; 
  out[44] += 0.7071067811865475*win[0]*fin[44]*subvol; 
  out[47] += 0.7071067811865475*win[0]*fin[47]*subvol; 
  out[48] += 0.7071067811865475*win[0]*fin[48]*subvol; 
  out[68] += 0.7071067811865475*win[0]*fin[68]*subvol; 
  out[69] += 0.7071067811865475*win[0]*fin[69]*subvol; 
  out[80] += 0.7071067811865475*win[0]*fin[80]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[44]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[48]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[13]*subvol; 
  out[4] += 0.5*win[0]*fin[17]*subvol; 
  out[5] += 0.5*win[0]*fin[20]*subvol; 
  out[6] += 0.5*win[0]*fin[44]*subvol; 
  out[7] += 0.5*win[0]*fin[48]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[13], 1.4142135623730951*win[0]*fin[13]*subvol); 
  atomicAdd(&out[17], 1.4142135623730951*win[0]*fin[17]*subvol); 
  atomicAdd(&out[20], 1.4142135623730951*win[0]*fin[20]*subvol); 
  atomicAdd(&out[44], 1.4142135623730951*win[0]*fin[44]*subvol); 
  atomicAdd(&out[48], 1.4142135623730951*win[0]*fin[48]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[13] += 1.4142135623730951*win[0]*fin[13]*subvol; 
  out[17] += 1.4142135623730951*win[0]*fin[17]*subvol; 
  out[20] += 1.4142135623730951*win[0]*fin[20]*subvol; 
  out[44] += 1.4142135623730951*win[0]*fin[44]*subvol; 
  out[48] += 1.4142135623730951*win[0]*fin[48]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[20]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[43]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[47]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[12]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[20]*subvol; 
  out[6] += 0.5*win[0]*fin[43]*subvol; 
  out[7] += 0.5*win[0]*fin[47]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[12], 1.4142135623730951*win[0]*fin[12]*subvol); 
  atomicAdd(&out[16], 1.4142135623730951*win[0]*fin[16]*subvol); 
  atomicAdd(&out[20], 1.4142135623730951*win[0]*fin[20]*subvol); 
  atomicAdd(&out[43], 1.4142135623730951*win[0]*fin[43]*subvol); 
  atomicAdd(&out[47], 1.4142135623730951*win[0]*fin[47]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[12] += 1.4142135623730951*win[0]*fin[12]*subvol; 
  out[16] += 1.4142135623730951*win[0]*fin[16]*subvol; 
  out[20] += 1.4142135623730951*win[0]*fin[20]*subvol; 
  out[43] += 1.4142135623730951*win[0]*fin[43]*subvol; 
  out[47] += 1.4142135623730951*win[0]*fin[47]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[20]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[20]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[5], 2.8284271247461907*win[0]*fin[5]*subvol); 
  atomicAdd(&out[20], 2.8284271247461907*win[0]*fin[20]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[5] += 2.8284271247461907*win[0]*fin[5]*subvol; 
  out[20] += 2.8284271247461907*win[0]*fin[20]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[15], 0.25*win[0]*fin[21]*subvol); 
  atomicAdd(&out[16], 0.25*win[0]*fin[22]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[23]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[24]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[31]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[35]*subvol); 
  atomicAdd(&out[24], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[25], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[26], 0.25*win[0]*fin[38]*subvol); 
  atomicAdd(&out[27], 0.25*win[0]*fin[39]*subvol); 
  atomicAdd(&out[28], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[29], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[30], 0.25*win[0]*fin[42]*subvol); 
  atomicAdd(&out[31], 0.25*win[0]*fin[51]*subvol); 
  atomicAdd(&out[32], 0.25*win[0]*fin[56]*subvol); 
  atomicAdd(&out[33], 0.25*win[0]*fin[57]*subvol); 
  atomicAdd(&out[34], 0.25*win[0]*fin[58]*subvol); 
  atomicAdd(&out[35], 0.25*win[0]*fin[59]*subvol); 
  atomicAdd(&out[36], 0.25*win[0]*fin[60]*subvol); 
  atomicAdd(&out[37], 0.25*win[0]*fin[61]*subvol); 
  atomicAdd(&out[38], 0.25*win[0]*fin[62]*subvol); 
  atomicAdd(&out[39], 0.25*win[0]*fin[63]*subvol); 
  atomicAdd(&out[40], 0.25*win[0]*fin[64]*subvol); 
  atomicAdd(&out[41], 0.25*win[0]*fin[65]*subvol); 
  atomicAdd(&out[42], 0.25*win[0]*fin[66]*subvol); 
  atomicAdd(&out[43], 0.25*win[0]*fin[67]*subvol); 
  atomicAdd(&out[44], 0.25*win[0]*fin[87]*subvol); 
  atomicAdd(&out[45], 0.25*win[0]*fin[88]*subvol); 
  atomicAdd(&out[46], 0.25*win[0]*fin[89]*subvol); 
  atomicAdd(&out[47], 0.25*win[0]*fin[90]*subvol); 
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
  out[15] += 0.25*win[0]*fin[21]*subvol; 
  out[16] += 0.25*win[0]*fin[22]*subvol; 
  out[17] += 0.25*win[0]*fin[23]*subvol; 
  out[18] += 0.25*win[0]*fin[24]*subvol; 
  out[19] += 0.25*win[0]*fin[31]*subvol; 
  out[20] += 0.25*win[0]*fin[32]*subvol; 
  out[21] += 0.25*win[0]*fin[33]*subvol; 
  out[22] += 0.25*win[0]*fin[34]*subvol; 
  out[23] += 0.25*win[0]*fin[35]*subvol; 
  out[24] += 0.25*win[0]*fin[36]*subvol; 
  out[25] += 0.25*win[0]*fin[37]*subvol; 
  out[26] += 0.25*win[0]*fin[38]*subvol; 
  out[27] += 0.25*win[0]*fin[39]*subvol; 
  out[28] += 0.25*win[0]*fin[40]*subvol; 
  out[29] += 0.25*win[0]*fin[41]*subvol; 
  out[30] += 0.25*win[0]*fin[42]*subvol; 
  out[31] += 0.25*win[0]*fin[51]*subvol; 
  out[32] += 0.25*win[0]*fin[56]*subvol; 
  out[33] += 0.25*win[0]*fin[57]*subvol; 
  out[34] += 0.25*win[0]*fin[58]*subvol; 
  out[35] += 0.25*win[0]*fin[59]*subvol; 
  out[36] += 0.25*win[0]*fin[60]*subvol; 
  out[37] += 0.25*win[0]*fin[61]*subvol; 
  out[38] += 0.25*win[0]*fin[62]*subvol; 
  out[39] += 0.25*win[0]*fin[63]*subvol; 
  out[40] += 0.25*win[0]*fin[64]*subvol; 
  out[41] += 0.25*win[0]*fin[65]*subvol; 
  out[42] += 0.25*win[0]*fin[66]*subvol; 
  out[43] += 0.25*win[0]*fin[67]*subvol; 
  out[44] += 0.25*win[0]*fin[87]*subvol; 
  out[45] += 0.25*win[0]*fin[88]*subvol; 
  out[46] += 0.25*win[0]*fin[89]*subvol; 
  out[47] += 0.25*win[0]*fin[90]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[21], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[22], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[23], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[24], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[31], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[35], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[37], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[38], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[40], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[41], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[51], 0.3535533905932737*win[0]*fin[51]*subvol); 
  atomicAdd(&out[56], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.3535533905932737*win[0]*fin[58]*subvol); 
  atomicAdd(&out[59], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[60], 0.3535533905932737*win[0]*fin[60]*subvol); 
  atomicAdd(&out[61], 0.3535533905932737*win[0]*fin[61]*subvol); 
  atomicAdd(&out[62], 0.3535533905932737*win[0]*fin[62]*subvol); 
  atomicAdd(&out[63], 0.3535533905932737*win[0]*fin[63]*subvol); 
  atomicAdd(&out[64], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[65], 0.3535533905932737*win[0]*fin[65]*subvol); 
  atomicAdd(&out[66], 0.3535533905932737*win[0]*fin[66]*subvol); 
  atomicAdd(&out[67], 0.3535533905932737*win[0]*fin[67]*subvol); 
  atomicAdd(&out[87], 0.3535533905932737*win[0]*fin[87]*subvol); 
  atomicAdd(&out[88], 0.3535533905932737*win[0]*fin[88]*subvol); 
  atomicAdd(&out[89], 0.3535533905932737*win[0]*fin[89]*subvol); 
  atomicAdd(&out[90], 0.3535533905932737*win[0]*fin[90]*subvol); 
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
  out[21] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[22] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[23] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[24] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[31] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[35] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[37] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[38] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[39] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[40] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[41] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[42] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[51] += 0.3535533905932737*win[0]*fin[51]*subvol; 
  out[56] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[57] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[58] += 0.3535533905932737*win[0]*fin[58]*subvol; 
  out[59] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[60] += 0.3535533905932737*win[0]*fin[60]*subvol; 
  out[61] += 0.3535533905932737*win[0]*fin[61]*subvol; 
  out[62] += 0.3535533905932737*win[0]*fin[62]*subvol; 
  out[63] += 0.3535533905932737*win[0]*fin[63]*subvol; 
  out[64] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[65] += 0.3535533905932737*win[0]*fin[65]*subvol; 
  out[66] += 0.3535533905932737*win[0]*fin[66]*subvol; 
  out[67] += 0.3535533905932737*win[0]*fin[67]*subvol; 
  out[87] += 0.3535533905932737*win[0]*fin[87]*subvol; 
  out[88] += 0.3535533905932737*win[0]*fin[88]*subvol; 
  out[89] += 0.3535533905932737*win[0]*fin[89]*subvol; 
  out[90] += 0.3535533905932737*win[0]*fin[90]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[62]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[64]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[67]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[62]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[64]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[67]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[24], 0.7071067811865475*win[0]*fin[24]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[38], 0.7071067811865475*win[0]*fin[38]*subvol); 
  atomicAdd(&out[39], 0.7071067811865475*win[0]*fin[39]*subvol); 
  atomicAdd(&out[41], 0.7071067811865475*win[0]*fin[41]*subvol); 
  atomicAdd(&out[42], 0.7071067811865475*win[0]*fin[42]*subvol); 
  atomicAdd(&out[62], 0.7071067811865475*win[0]*fin[62]*subvol); 
  atomicAdd(&out[64], 0.7071067811865475*win[0]*fin[64]*subvol); 
  atomicAdd(&out[67], 0.7071067811865475*win[0]*fin[67]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[24] += 0.7071067811865475*win[0]*fin[24]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[38] += 0.7071067811865475*win[0]*fin[38]*subvol; 
  out[39] += 0.7071067811865475*win[0]*fin[39]*subvol; 
  out[41] += 0.7071067811865475*win[0]*fin[41]*subvol; 
  out[42] += 0.7071067811865475*win[0]*fin[42]*subvol; 
  out[62] += 0.7071067811865475*win[0]*fin[62]*subvol; 
  out[64] += 0.7071067811865475*win[0]*fin[64]*subvol; 
  out[67] += 0.7071067811865475*win[0]*fin[67]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[23]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[39]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[42]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[61]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[63]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[66]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[39]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[42]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[61]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[63]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[66]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[23], 0.7071067811865475*win[0]*fin[23]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 0.7071067811865475*win[0]*fin[35]*subvol); 
  atomicAdd(&out[37], 0.7071067811865475*win[0]*fin[37]*subvol); 
  atomicAdd(&out[39], 0.7071067811865475*win[0]*fin[39]*subvol); 
  atomicAdd(&out[40], 0.7071067811865475*win[0]*fin[40]*subvol); 
  atomicAdd(&out[42], 0.7071067811865475*win[0]*fin[42]*subvol); 
  atomicAdd(&out[61], 0.7071067811865475*win[0]*fin[61]*subvol); 
  atomicAdd(&out[63], 0.7071067811865475*win[0]*fin[63]*subvol); 
  atomicAdd(&out[66], 0.7071067811865475*win[0]*fin[66]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[23] += 0.7071067811865475*win[0]*fin[23]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[35] += 0.7071067811865475*win[0]*fin[35]*subvol; 
  out[37] += 0.7071067811865475*win[0]*fin[37]*subvol; 
  out[39] += 0.7071067811865475*win[0]*fin[39]*subvol; 
  out[40] += 0.7071067811865475*win[0]*fin[40]*subvol; 
  out[42] += 0.7071067811865475*win[0]*fin[42]*subvol; 
  out[61] += 0.7071067811865475*win[0]*fin[61]*subvol; 
  out[63] += 0.7071067811865475*win[0]*fin[63]*subvol; 
  out[66] += 0.7071067811865475*win[0]*fin[66]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[39]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[42]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[11]*subvol; 
  out[4] += 0.5*win[0]*fin[18]*subvol; 
  out[5] += 0.5*win[0]*fin[19]*subvol; 
  out[6] += 0.5*win[0]*fin[39]*subvol; 
  out[7] += 0.5*win[0]*fin[42]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[11], 1.4142135623730951*win[0]*fin[11]*subvol); 
  atomicAdd(&out[18], 1.4142135623730951*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 1.4142135623730951*win[0]*fin[19]*subvol); 
  atomicAdd(&out[39], 1.4142135623730951*win[0]*fin[39]*subvol); 
  atomicAdd(&out[42], 1.4142135623730951*win[0]*fin[42]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[11] += 1.4142135623730951*win[0]*fin[11]*subvol; 
  out[18] += 1.4142135623730951*win[0]*fin[18]*subvol; 
  out[19] += 1.4142135623730951*win[0]*fin[19]*subvol; 
  out[39] += 1.4142135623730951*win[0]*fin[39]*subvol; 
  out[42] += 1.4142135623730951*win[0]*fin[42]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[38]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[59]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[60]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[65]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[38]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[59]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[60]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[65]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.7071067811865475*win[0]*fin[19]*subvol); 
  atomicAdd(&out[22], 0.7071067811865475*win[0]*fin[22]*subvol); 
  atomicAdd(&out[31], 0.7071067811865475*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[37], 0.7071067811865475*win[0]*fin[37]*subvol); 
  atomicAdd(&out[38], 0.7071067811865475*win[0]*fin[38]*subvol); 
  atomicAdd(&out[40], 0.7071067811865475*win[0]*fin[40]*subvol); 
  atomicAdd(&out[41], 0.7071067811865475*win[0]*fin[41]*subvol); 
  atomicAdd(&out[59], 0.7071067811865475*win[0]*fin[59]*subvol); 
  atomicAdd(&out[60], 0.7071067811865475*win[0]*fin[60]*subvol); 
  atomicAdd(&out[65], 0.7071067811865475*win[0]*fin[65]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[19] += 0.7071067811865475*win[0]*fin[19]*subvol; 
  out[22] += 0.7071067811865475*win[0]*fin[22]*subvol; 
  out[31] += 0.7071067811865475*win[0]*fin[31]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[37] += 0.7071067811865475*win[0]*fin[37]*subvol; 
  out[38] += 0.7071067811865475*win[0]*fin[38]*subvol; 
  out[40] += 0.7071067811865475*win[0]*fin[40]*subvol; 
  out[41] += 0.7071067811865475*win[0]*fin[41]*subvol; 
  out[59] += 0.7071067811865475*win[0]*fin[59]*subvol; 
  out[60] += 0.7071067811865475*win[0]*fin[60]*subvol; 
  out[65] += 0.7071067811865475*win[0]*fin[65]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[38]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[41]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[10]*subvol; 
  out[4] += 0.5*win[0]*fin[17]*subvol; 
  out[5] += 0.5*win[0]*fin[19]*subvol; 
  out[6] += 0.5*win[0]*fin[38]*subvol; 
  out[7] += 0.5*win[0]*fin[41]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[10], 1.4142135623730951*win[0]*fin[10]*subvol); 
  atomicAdd(&out[17], 1.4142135623730951*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 1.4142135623730951*win[0]*fin[19]*subvol); 
  atomicAdd(&out[38], 1.4142135623730951*win[0]*fin[38]*subvol); 
  atomicAdd(&out[41], 1.4142135623730951*win[0]*fin[41]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[10] += 1.4142135623730951*win[0]*fin[10]*subvol; 
  out[17] += 1.4142135623730951*win[0]*fin[17]*subvol; 
  out[19] += 1.4142135623730951*win[0]*fin[19]*subvol; 
  out[38] += 1.4142135623730951*win[0]*fin[38]*subvol; 
  out[41] += 1.4142135623730951*win[0]*fin[41]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[37]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[40]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[9]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[19]*subvol; 
  out[6] += 0.5*win[0]*fin[37]*subvol; 
  out[7] += 0.5*win[0]*fin[40]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[9], 1.4142135623730951*win[0]*fin[9]*subvol); 
  atomicAdd(&out[16], 1.4142135623730951*win[0]*fin[16]*subvol); 
  atomicAdd(&out[19], 1.4142135623730951*win[0]*fin[19]*subvol); 
  atomicAdd(&out[37], 1.4142135623730951*win[0]*fin[37]*subvol); 
  atomicAdd(&out[40], 1.4142135623730951*win[0]*fin[40]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[9] += 1.4142135623730951*win[0]*fin[9]*subvol; 
  out[16] += 1.4142135623730951*win[0]*fin[16]*subvol; 
  out[19] += 1.4142135623730951*win[0]*fin[19]*subvol; 
  out[37] += 1.4142135623730951*win[0]*fin[37]*subvol; 
  out[40] += 1.4142135623730951*win[0]*fin[40]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvxvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[19]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[19]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvxvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.8284271247461907*win[0]*fin[4]*subvol); 
  atomicAdd(&out[19], 2.8284271247461907*win[0]*fin[19]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[4] += 2.8284271247461907*win[0]*fin[4]*subvol; 
  out[19] += 2.8284271247461907*win[0]*fin[19]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[31]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[16], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[17], 0.3535533905932737*win[0]*fin[56]*subvol); 
  atomicAdd(&out[18], 0.3535533905932737*win[0]*fin[57]*subvol); 
  atomicAdd(&out[19], 0.3535533905932737*win[0]*fin[58]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[31]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[16] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[17] += 0.3535533905932737*win[0]*fin[56]*subvol; 
  out[18] += 0.3535533905932737*win[0]*fin[57]*subvol; 
  out[19] += 0.3535533905932737*win[0]*fin[58]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[17], 0.7071067811865475*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.7071067811865475*win[0]*fin[18]*subvol); 
  atomicAdd(&out[21], 0.7071067811865475*win[0]*fin[21]*subvol); 
  atomicAdd(&out[31], 0.7071067811865475*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[35], 0.7071067811865475*win[0]*fin[35]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[56], 0.7071067811865475*win[0]*fin[56]*subvol); 
  atomicAdd(&out[57], 0.7071067811865475*win[0]*fin[57]*subvol); 
  atomicAdd(&out[58], 0.7071067811865475*win[0]*fin[58]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[18] += 0.7071067811865475*win[0]*fin[18]*subvol; 
  out[21] += 0.7071067811865475*win[0]*fin[21]*subvol; 
  out[31] += 0.7071067811865475*win[0]*fin[31]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[35] += 0.7071067811865475*win[0]*fin[35]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[56] += 0.7071067811865475*win[0]*fin[56]*subvol; 
  out[57] += 0.7071067811865475*win[0]*fin[57]*subvol; 
  out[58] += 0.7071067811865475*win[0]*fin[58]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[34]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[36]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[8]*subvol; 
  out[4] += 0.5*win[0]*fin[17]*subvol; 
  out[5] += 0.5*win[0]*fin[18]*subvol; 
  out[6] += 0.5*win[0]*fin[34]*subvol; 
  out[7] += 0.5*win[0]*fin[36]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[8], 1.4142135623730951*win[0]*fin[8]*subvol); 
  atomicAdd(&out[17], 1.4142135623730951*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 1.4142135623730951*win[0]*fin[18]*subvol); 
  atomicAdd(&out[34], 1.4142135623730951*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 1.4142135623730951*win[0]*fin[36]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[8] += 1.4142135623730951*win[0]*fin[8]*subvol; 
  out[17] += 1.4142135623730951*win[0]*fin[17]*subvol; 
  out[18] += 1.4142135623730951*win[0]*fin[18]*subvol; 
  out[34] += 1.4142135623730951*win[0]*fin[34]*subvol; 
  out[36] += 1.4142135623730951*win[0]*fin[36]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[33]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[35]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[7]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[18]*subvol; 
  out[6] += 0.5*win[0]*fin[33]*subvol; 
  out[7] += 0.5*win[0]*fin[35]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[7], 1.4142135623730951*win[0]*fin[7]*subvol); 
  atomicAdd(&out[16], 1.4142135623730951*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 1.4142135623730951*win[0]*fin[18]*subvol); 
  atomicAdd(&out[33], 1.4142135623730951*win[0]*fin[33]*subvol); 
  atomicAdd(&out[35], 1.4142135623730951*win[0]*fin[35]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[7] += 1.4142135623730951*win[0]*fin[7]*subvol; 
  out[16] += 1.4142135623730951*win[0]*fin[16]*subvol; 
  out[18] += 1.4142135623730951*win[0]*fin[18]*subvol; 
  out[33] += 1.4142135623730951*win[0]*fin[33]*subvol; 
  out[35] += 1.4142135623730951*win[0]*fin[35]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[18]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[18]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.8284271247461907*win[0]*fin[3]*subvol); 
  atomicAdd(&out[18], 2.8284271247461907*win[0]*fin[18]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[3] += 2.8284271247461907*win[0]*fin[3]*subvol; 
  out[18] += 2.8284271247461907*win[0]*fin[18]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[31]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[32]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[6]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[17]*subvol; 
  out[6] += 0.5*win[0]*fin[31]*subvol; 
  out[7] += 0.5*win[0]*fin[32]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[6], 1.4142135623730951*win[0]*fin[6]*subvol); 
  atomicAdd(&out[16], 1.4142135623730951*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 1.4142135623730951*win[0]*fin[17]*subvol); 
  atomicAdd(&out[31], 1.4142135623730951*win[0]*fin[31]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[6] += 1.4142135623730951*win[0]*fin[6]*subvol; 
  out[16] += 1.4142135623730951*win[0]*fin[16]*subvol; 
  out[17] += 1.4142135623730951*win[0]*fin[17]*subvol; 
  out[31] += 1.4142135623730951*win[0]*fin[31]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[17]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[17]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.8284271247461907*win[0]*fin[2]*subvol); 
  atomicAdd(&out[17], 2.8284271247461907*win[0]*fin[17]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[2] += 2.8284271247461907*win[0]*fin[2]*subvol; 
  out[17] += 2.8284271247461907*win[0]*fin[17]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[16]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgyvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.8284271247461907*win[0]*fin[1]*subvol); 
  atomicAdd(&out[16], 2.8284271247461907*win[0]*fin[16]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[1] += 2.8284271247461907*win[0]*fin[1]*subvol; 
  out[16] += 2.8284271247461907*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvxvyvz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x3v_ser_p2_avgxyvxvyvz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 5.656854249492382*fin[0]*win[0]*subvol); 
#else
  out[0] += 5.656854249492382*fin[0]*win[0]*subvol; 
#endif
} 

