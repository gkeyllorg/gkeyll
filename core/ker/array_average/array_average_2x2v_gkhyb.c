#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[10]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[16]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[19]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[3]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[10]*subvol; 
  out[4] += 0.5*win[0]*fin[16]*subvol; 
  out[5] += 0.5*win[0]*fin[19]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], win[0]*fin[3]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[10], win[0]*fin[10]*subvol); 
  atomicAdd(&out[16], win[0]*fin[16]*subvol); 
  atomicAdd(&out[19], win[0]*fin[19]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[3] += win[0]*fin[3]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[10] += win[0]*fin[10]*subvol; 
  out[16] += win[0]*fin[16]*subvol; 
  out[19] += win[0]*fin[19]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.3535533905932737*fin[0]*win[0]*subvol; 
  out[1] += 0.3535533905932737*win[0]*fin[1]*subvol; 
  out[2] += 0.3535533905932737*win[0]*fin[2]*subvol; 
  out[3] += 0.3535533905932737*win[0]*fin[4]*subvol; 
  out[4] += 0.3535533905932737*win[0]*fin[5]*subvol; 
  out[5] += 0.3535533905932737*win[0]*fin[8]*subvol; 
  out[6] += 0.3535533905932737*win[0]*fin[9]*subvol; 
  out[7] += 0.3535533905932737*win[0]*fin[12]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[4] += 0.5*win[0]*fin[4]*subvol; 
  out[5] += 0.5*win[0]*fin[5]*subvol; 
  out[8] += 0.5*win[0]*fin[8]*subvol; 
  out[9] += 0.5*win[0]*fin[9]*subvol; 
  out[12] += 0.5*win[0]*fin[12]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[9]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[9], win[0]*fin[9]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[9] += win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgyvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[8]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[4]*subvol; 
  out[3] += 0.5*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgyvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
  atomicAdd(&out[8], win[0]*fin[8]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
  out[8] += win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxyvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[4]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[4]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxyvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[4], 2.0*win[0]*fin[4]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[4] += 2.0*win[0]*fin[4]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgymu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgymu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxymu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxymu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 2.0*win[0]*fin[2]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[2] += 2.0*win[0]*fin[2]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgyvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgyvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 2.0*win[0]*fin[1]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
  out[1] += 2.0*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxyvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x2v_gkhyb_p1_avgxyvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 4.0*fin[0]*win[0]*subvol); 
#else
  out[0] += 4.0*fin[0]*win[0]*subvol; 
#endif
} 

