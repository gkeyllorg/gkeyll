#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[1]*win[1]+fin[0]*win[0])*subvol); 
#else
  out[0] += (fin[1]*win[1]+fin[0]*win[0])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (1.4142135623730951*fin[1]*win[1]+1.4142135623730951*fin[0]*win[0])*subvol); 
#else
  out[0] += (1.4142135623730951*fin[1]*win[1]+1.4142135623730951*fin[0]*win[0])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
#else
  out[0] += (fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x_ser_p2_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (1.4142135623730951*fin[2]*win[2]+1.4142135623730951*fin[1]*win[1]+1.4142135623730951*fin[0]*win[0])*subvol); 
#else
  out[0] += (1.4142135623730951*fin[2]*win[2]+1.4142135623730951*fin[1]*win[1]+1.4142135623730951*fin[0]*win[0])*subvol; 
#endif
} 

