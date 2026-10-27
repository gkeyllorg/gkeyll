#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[15]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[35]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[42]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[35]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[42]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[16], 0.25*win[0]*fin[32]*subvol); 
  atomicAdd(&out[17], 0.25*win[0]*fin[33]*subvol); 
  atomicAdd(&out[18], 0.25*win[0]*fin[34]*subvol); 
  atomicAdd(&out[19], 0.25*win[0]*fin[36]*subvol); 
  atomicAdd(&out[20], 0.25*win[0]*fin[37]*subvol); 
  atomicAdd(&out[21], 0.25*win[0]*fin[40]*subvol); 
  atomicAdd(&out[22], 0.25*win[0]*fin[41]*subvol); 
  atomicAdd(&out[23], 0.25*win[0]*fin[44]*subvol); 
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
  out[16] += 0.25*win[0]*fin[32]*subvol; 
  out[17] += 0.25*win[0]*fin[33]*subvol; 
  out[18] += 0.25*win[0]*fin[34]*subvol; 
  out[19] += 0.25*win[0]*fin[36]*subvol; 
  out[20] += 0.25*win[0]*fin[37]*subvol; 
  out[21] += 0.25*win[0]*fin[40]*subvol; 
  out[22] += 0.25*win[0]*fin[41]*subvol; 
  out[23] += 0.25*win[0]*fin[44]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[32], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[37], 0.3535533905932737*win[0]*fin[37]*subvol); 
  atomicAdd(&out[40], 0.3535533905932737*win[0]*fin[40]*subvol); 
  atomicAdd(&out[41], 0.3535533905932737*win[0]*fin[41]*subvol); 
  atomicAdd(&out[44], 0.3535533905932737*win[0]*fin[44]*subvol); 
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
  out[32] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[33] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[34] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[36] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[37] += 0.3535533905932737*win[0]*fin[37]*subvol; 
  out[40] += 0.3535533905932737*win[0]*fin[40]*subvol; 
  out[41] += 0.3535533905932737*win[0]*fin[41]*subvol; 
  out[44] += 0.3535533905932737*win[0]*fin[44]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[41]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[41]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[41], 0.7071067811865475*win[0]*fin[41]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[24] += 0.7071067811865475*win[0]*fin[24]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[41] += 0.7071067811865475*win[0]*fin[41]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[36]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[40]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[15]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[23]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[36]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[40]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[36], 0.7071067811865475*win[0]*fin[36]*subvol); 
  atomicAdd(&out[40], 0.7071067811865475*win[0]*fin[40]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
  out[23] += 0.7071067811865475*win[0]*fin[23]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[36] += 0.7071067811865475*win[0]*fin[36]*subvol; 
  out[40] += 0.7071067811865475*win[0]*fin[40]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[15]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[36]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[4]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[15]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[36]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyz_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[15], 1.4142135623730951*win[0]*fin[15]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[36], 1.4142135623730951*win[0]*fin[36]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[15] += 1.4142135623730951*win[0]*fin[15]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[36] += 1.4142135623730951*win[0]*fin[36]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[22]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[22] += 0.7071067811865475*win[0]*fin[22]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[21]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[21] += 0.7071067811865475*win[0]*fin[21]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[14]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[14]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[14], 1.4142135623730951*win[0]*fin[14]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[14] += 1.4142135623730951*win[0]*fin[14]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgzvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[20]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgzvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[20] += 0.7071067811865475*win[0]*fin[20]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxzvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[13]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[13]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxzvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[13], 1.4142135623730951*win[0]*fin[13]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[13] += 1.4142135623730951*win[0]*fin[13]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyzvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[12]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[5]*subvol; 
  out[3] += 0.5*win[0]*fin[12]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyzvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[5], 1.4142135623730951*win[0]*fin[5]*subvol); 
  atomicAdd(&out[12], 1.4142135623730951*win[0]*fin[12]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[5] += 1.4142135623730951*win[0]*fin[5]*subvol; 
  out[12] += 1.4142135623730951*win[0]*fin[12]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyzvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[5]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[5]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyzvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[5], 2.8284271247461907*win[0]*fin[5]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[5] += 2.8284271247461907*win[0]*fin[5]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgymu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgymu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxymu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[35]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[11]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[35]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxymu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[11], 1.4142135623730951*win[0]*fin[11]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[35], 1.4142135623730951*win[0]*fin[35]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[11] += 1.4142135623730951*win[0]*fin[11]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[35] += 1.4142135623730951*win[0]*fin[35]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgzmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[32]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[33]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[34]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[37]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[32]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[33]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[34]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[37]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgzmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
  atomicAdd(&out[32], 0.7071067811865475*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 0.7071067811865475*win[0]*fin[33]*subvol); 
  atomicAdd(&out[34], 0.7071067811865475*win[0]*fin[34]*subvol); 
  atomicAdd(&out[37], 0.7071067811865475*win[0]*fin[37]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
  out[17] += 0.7071067811865475*win[0]*fin[17]*subvol; 
  out[32] += 0.7071067811865475*win[0]*fin[32]*subvol; 
  out[33] += 0.7071067811865475*win[0]*fin[33]*subvol; 
  out[34] += 0.7071067811865475*win[0]*fin[34]*subvol; 
  out[37] += 0.7071067811865475*win[0]*fin[37]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxzmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[34]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[10]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[34]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxzmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[10], 1.4142135623730951*win[0]*fin[10]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[34], 1.4142135623730951*win[0]*fin[34]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[10] += 1.4142135623730951*win[0]*fin[10]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[34] += 1.4142135623730951*win[0]*fin[34]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyzmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[32]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[33]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[9]*subvol; 
  out[4] += 0.5*win[0]*fin[32]*subvol; 
  out[5] += 0.5*win[0]*fin[33]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyzmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], 1.4142135623730951*win[0]*fin[4]*subvol); 
  atomicAdd(&out[9], 1.4142135623730951*win[0]*fin[9]*subvol); 
  atomicAdd(&out[32], 1.4142135623730951*win[0]*fin[32]*subvol); 
  atomicAdd(&out[33], 1.4142135623730951*win[0]*fin[33]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[4] += 1.4142135623730951*win[0]*fin[4]*subvol; 
  out[9] += 1.4142135623730951*win[0]*fin[9]*subvol; 
  out[32] += 1.4142135623730951*win[0]*fin[32]*subvol; 
  out[33] += 1.4142135623730951*win[0]*fin[33]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyzmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[32]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[32]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyzmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.8284271247461907*win[0]*fin[4]*subvol); 
  atomicAdd(&out[32], 2.8284271247461907*win[0]*fin[32]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[4] += 2.8284271247461907*win[0]*fin[4]*subvol; 
  out[32] += 2.8284271247461907*win[0]*fin[32]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[8]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[8], 1.4142135623730951*win[0]*fin[8]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[8] += 1.4142135623730951*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[7]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[7]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[7], 1.4142135623730951*win[0]*fin[7]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[7] += 1.4142135623730951*win[0]*fin[7]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.8284271247461907*win[0]*fin[3]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[3] += 2.8284271247461907*win[0]*fin[3]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgzvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgzvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxzvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxzvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.8284271247461907*win[0]*fin[2]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[2] += 2.8284271247461907*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyzvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgyzvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.8284271247461907*win[0]*fin[1]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
  out[1] += 2.8284271247461907*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyzvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_3x2v_gkhyb_p1_avgxyzvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 5.656854249492382*fin[0]*win[0]*subvol); 
#else
  out[0] += 5.656854249492382*fin[0]*win[0]*subvol; 
#endif
} 

