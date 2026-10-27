#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[14]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[22]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[30]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[14]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[22]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[30]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[14], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[22], 0.5*win[0]*fin[22]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[26], 0.5*win[0]*fin[26]*subvol); 
  atomicAdd(&out[27], 0.5*win[0]*fin[27]*subvol); 
  atomicAdd(&out[30], 0.5*win[0]*fin[30]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[14] += 0.5*win[0]*fin[14]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[22] += 0.5*win[0]*fin[22]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[26] += 0.5*win[0]*fin[26]*subvol; 
  out[27] += 0.5*win[0]*fin[27]*subvol; 
  out[30] += 0.5*win[0]*fin[30]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[10]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[13]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[19]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[21]*subvol); 
  atomicAdd(&out[12], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[13], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[14], 0.3535533905932737*win[0]*fin[27]*subvol); 
  atomicAdd(&out[15], 0.3535533905932737*win[0]*fin[29]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[10]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[13]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[19]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[21]*subvol; 
  out[12] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[13] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[14] += 0.3535533905932737*win[0]*fin[27]*subvol; 
  out[15] += 0.3535533905932737*win[0]*fin[29]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[13], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[19], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[21], 0.5*win[0]*fin[21]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], 0.5*win[0]*fin[25]*subvol); 
  atomicAdd(&out[27], 0.5*win[0]*fin[27]*subvol); 
  atomicAdd(&out[29], 0.5*win[0]*fin[29]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[10] += 0.5*win[0]*fin[10]*subvol; 
  out[13] += 0.5*win[0]*fin[13]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[19] += 0.5*win[0]*fin[19]*subvol; 
  out[21] += 0.5*win[0]*fin[21]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[25] += 0.5*win[0]*fin[25]*subvol; 
  out[27] += 0.5*win[0]*fin[27]*subvol; 
  out[29] += 0.5*win[0]*fin[29]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[19]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[27]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[10]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[24]*subvol; 
  out[6] += 0.5*win[0]*fin[19]*subvol; 
  out[7] += 0.5*win[0]*fin[27]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[10], win[0]*fin[10]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[19], win[0]*fin[19]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[27], win[0]*fin[27]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[10] += win[0]*fin[10]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[19] += win[0]*fin[19]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[27] += win[0]*fin[27]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[9]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[12]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[24]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[25]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[26]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[28]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[12]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[24]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[25]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[26]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[28]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[8], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[12], 0.5*win[0]*fin[12]*subvol); 
  atomicAdd(&out[24], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], 0.5*win[0]*fin[25]*subvol); 
  atomicAdd(&out[26], 0.5*win[0]*fin[26]*subvol); 
  atomicAdd(&out[28], 0.5*win[0]*fin[28]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
  out[24] += 0.5*win[0]*fin[24]*subvol; 
  out[25] += 0.5*win[0]*fin[25]*subvol; 
  out[26] += 0.5*win[0]*fin[26]*subvol; 
  out[28] += 0.5*win[0]*fin[28]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[26]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[9]*subvol; 
  out[4] += 0.5*win[0]*fin[24]*subvol; 
  out[5] += 0.5*win[0]*fin[26]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[9], win[0]*fin[9]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[26], win[0]*fin[26]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[9] += win[0]*fin[9]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[26] += win[0]*fin[26]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[24]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[25]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[8]*subvol; 
  out[4] += 0.5*win[0]*fin[24]*subvol; 
  out[5] += 0.5*win[0]*fin[25]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[8], win[0]*fin[8]*subvol); 
  atomicAdd(&out[24], win[0]*fin[24]*subvol); 
  atomicAdd(&out[25], win[0]*fin[25]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[8] += win[0]*fin[8]*subvol; 
  out[24] += win[0]*fin[24]*subvol; 
  out[25] += win[0]*fin[25]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxyvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[24]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[24]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxyvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
  atomicAdd(&out[24], 2.0*win[0]*fin[24]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
  out[24] += 2.0*win[0]*fin[24]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.3535533905932737*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.3535533905932737*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.3535533905932737*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.3535533905932737*win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], 0.3535533905932737*win[0]*fin[5]*subvol); 
  atomicAdd(&out[5], 0.3535533905932737*win[0]*fin[6]*subvol); 
  atomicAdd(&out[6], 0.3535533905932737*win[0]*fin[7]*subvol); 
  atomicAdd(&out[7], 0.3535533905932737*win[0]*fin[11]*subvol); 
  atomicAdd(&out[8], 0.3535533905932737*win[0]*fin[16]*subvol); 
  atomicAdd(&out[9], 0.3535533905932737*win[0]*fin[17]*subvol); 
  atomicAdd(&out[10], 0.3535533905932737*win[0]*fin[18]*subvol); 
  atomicAdd(&out[11], 0.3535533905932737*win[0]*fin[20]*subvol); 
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[3]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[6]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[7]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[11]*subvol; 
  out[8] += 0.3535533905932737*win[0]*fin[16]*subvol; 
  out[9] += 0.3535533905932737*win[0]*fin[17]*subvol; 
  out[10] += 0.3535533905932737*win[0]*fin[18]*subvol; 
  out[11] += 0.3535533905932737*win[0]*fin[20]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[11], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[16], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], 0.5*win[0]*fin[17]*subvol); 
  atomicAdd(&out[18], 0.5*win[0]*fin[18]*subvol); 
  atomicAdd(&out[20], 0.5*win[0]*fin[20]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[3]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[6] += 0.5*win[0]*fin[6]*subvol; 
  out[7] += 0.5*win[0]*fin[7]*subvol; 
  out[11] += 0.5*win[0]*fin[11]*subvol; 
  out[16] += 0.5*win[0]*fin[16]*subvol; 
  out[17] += 0.5*win[0]*fin[17]*subvol; 
  out[18] += 0.5*win[0]*fin[18]*subvol; 
  out[20] += 0.5*win[0]*fin[20]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[18]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[7]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[18]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[7], win[0]*fin[7]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[18], win[0]*fin[18]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[7] += win[0]*fin[7]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[18] += win[0]*fin[18]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[17]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[6]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[17]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], win[0]*fin[6]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[17], win[0]*fin[17]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[6] += win[0]*fin[6]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[17] += win[0]*fin[17]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxyvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[16]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxyvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 2.0*win[0]*fin[3]*subvol); 
  atomicAdd(&out[16], 2.0*win[0]*fin[16]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[3] += 2.0*win[0]*fin[3]*subvol; 
  out[16] += 2.0*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[5]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[5]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[5], win[0]*fin[5]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[5] += win[0]*fin[5]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxyvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_hyb_p1_avgxyvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
#endif
} 

