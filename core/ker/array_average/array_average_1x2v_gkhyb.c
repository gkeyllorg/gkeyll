#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[10]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[6]*subvol; 
  out[4] += 0.5*win[0]*fin[8]*subvol; 
  out[5] += 0.5*win[0]*fin[10]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[10], 0.7071067811865475*win[0]*fin[10]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[10] += 0.7071067811865475*win[0]*fin[10]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[5]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[5]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgxvpar(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgxvpar_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[9]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[4]*subvol; 
  out[4] += 0.5*win[0]*fin[8]*subvol; 
  out[5] += 0.5*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgxmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[8]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgxmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 1.4142135623730951*win[0]*fin[2]*subvol); 
  atomicAdd(&out[8], 1.4142135623730951*win[0]*fin[8]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[2] += 1.4142135623730951*win[0]*fin[2]*subvol; 
  out[8] += 1.4142135623730951*win[0]*fin[8]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgxvparmu(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_gkhyb_p1_avgxvparmu_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
#endif
} 

