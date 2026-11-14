#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[6]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[14]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[16]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[2]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[6]*subvol; 
  out[4] += 0.5*win[0]*fin[8]*subvol; 
  out[5] += 0.5*win[0]*fin[9]*subvol; 
  out[6] += 0.5*win[0]*fin[14]*subvol; 
  out[7] += 0.5*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[6], 0.7071067811865475*win[0]*fin[6]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[14], 0.7071067811865475*win[0]*fin[14]*subvol); 
  atomicAdd(&out[16], 0.7071067811865475*win[0]*fin[16]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[6] += 0.7071067811865475*win[0]*fin[6]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[14] += 0.7071067811865475*win[0]*fin[14]*subvol; 
  out[16] += 0.7071067811865475*win[0]*fin[16]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[3]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[5]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[9]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[13]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[15]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[3]*subvol; 
  out[3] += 0.5*win[0]*fin[5]*subvol; 
  out[4] += 0.5*win[0]*fin[7]*subvol; 
  out[5] += 0.5*win[0]*fin[9]*subvol; 
  out[6] += 0.5*win[0]*fin[13]*subvol; 
  out[7] += 0.5*win[0]*fin[15]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[3], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[5], 0.7071067811865475*win[0]*fin[5]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[9], 0.7071067811865475*win[0]*fin[9]*subvol); 
  atomicAdd(&out[13], 0.7071067811865475*win[0]*fin[13]*subvol); 
  atomicAdd(&out[15], 0.7071067811865475*win[0]*fin[15]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[3] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[5] += 0.7071067811865475*win[0]*fin[5]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[9] += 0.7071067811865475*win[0]*fin[9]*subvol; 
  out[13] += 0.7071067811865475*win[0]*fin[13]*subvol; 
  out[15] += 0.7071067811865475*win[0]*fin[15]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgxvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[3]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[9]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[3]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgxvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[3], 1.4142135623730951*win[0]*fin[3]*subvol); 
  atomicAdd(&out[9], 1.4142135623730951*win[0]*fin[9]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[3] += 1.4142135623730951*win[0]*fin[3]*subvol; 
  out[9] += 1.4142135623730951*win[0]*fin[9]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.5*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.5*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.5*win[0]*fin[2]*subvol); 
  atomicAdd(&out[3], 0.5*win[0]*fin[4]*subvol); 
  atomicAdd(&out[4], 0.5*win[0]*fin[7]*subvol); 
  atomicAdd(&out[5], 0.5*win[0]*fin[8]*subvol); 
  atomicAdd(&out[6], 0.5*win[0]*fin[11]*subvol); 
  atomicAdd(&out[7], 0.5*win[0]*fin[12]*subvol); 
#else
  out[0] += 0.5*fin[0]*win[0]*subvol; 
  out[1] += 0.5*win[0]*fin[1]*subvol; 
  out[2] += 0.5*win[0]*fin[2]*subvol; 
  out[3] += 0.5*win[0]*fin[4]*subvol; 
  out[4] += 0.5*win[0]*fin[7]*subvol; 
  out[5] += 0.5*win[0]*fin[8]*subvol; 
  out[6] += 0.5*win[0]*fin[11]*subvol; 
  out[7] += 0.5*win[0]*fin[12]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], 0.7071067811865475*win[0]*fin[4]*subvol); 
  atomicAdd(&out[7], 0.7071067811865475*win[0]*fin[7]*subvol); 
  atomicAdd(&out[8], 0.7071067811865475*win[0]*fin[8]*subvol); 
  atomicAdd(&out[11], 0.7071067811865475*win[0]*fin[11]*subvol); 
  atomicAdd(&out[12], 0.7071067811865475*win[0]*fin[12]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[4] += 0.7071067811865475*win[0]*fin[4]*subvol; 
  out[7] += 0.7071067811865475*win[0]*fin[7]*subvol; 
  out[8] += 0.7071067811865475*win[0]*fin[8]*subvol; 
  out[11] += 0.7071067811865475*win[0]*fin[11]*subvol; 
  out[12] += 0.7071067811865475*win[0]*fin[12]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[7]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[7]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 1.4142135623730951*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 1.4142135623730951*win[0]*fin[1]*subvol); 
  atomicAdd(&out[7], 1.4142135623730951*win[0]*fin[7]*subvol); 
#else
  out[0] += 1.4142135623730951*fin[0]*win[0]*subvol; 
  out[1] += 1.4142135623730951*win[0]*fin[1]*subvol; 
  out[7] += 1.4142135623730951*win[0]*fin[7]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgxvxvy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x2v_ser_p2_avgxvxvy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.8284271247461907*fin[0]*win[0]*subvol); 
#else
  out[0] += 2.8284271247461907*fin[0]*win[0]*subvol; 
#endif
} 

