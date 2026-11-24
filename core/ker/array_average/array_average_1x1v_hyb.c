#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x1v_hyb_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[2]*subvol); 
  atomicAdd(&out[2], 0.7071067811865475*win[0]*fin[4]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[2]*subvol; 
  out[2] += 0.7071067811865475*win[0]*fin[4]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x1v_hyb_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[2], win[0]*fin[2]*subvol); 
  atomicAdd(&out[4], win[0]*fin[4]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[2] += win[0]*fin[2]*subvol; 
  out[4] += win[0]*fin[4]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x1v_hyb_p1_avgvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 0.7071067811865475*fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], 0.7071067811865475*win[0]*fin[1]*subvol); 
#else
  out[0] += 0.7071067811865475*fin[0]*win[0]*subvol; 
  out[1] += 0.7071067811865475*win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x1v_hyb_p1_avgvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
  atomicAdd(&out[1], win[0]*fin[1]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
  out[1] += win[0]*fin[1]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x1v_hyb_p1_avgxvx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], fin[0]*win[0]*subvol); 
#else
  out[0] += fin[0]*win[0]*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_1x1v_hyb_p1_avgxvx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], 2.0*fin[0]*win[0]*subvol); 
#else
  out[0] += 2.0*fin[0]*win[0]*subvol; 
#endif
} 

