#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol); 
  atomicAdd(&out[1], (0.7071067811865475*fin[1]*win[3]+0.7071067811865475*win[1]*fin[3]+0.7071067811865475*fin[0]*win[2]+0.7071067811865475*win[0]*fin[2])*subvol); 
#else
  out[0] += (0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[1]*win[3]+0.7071067811865475*win[1]*fin[3]+0.7071067811865475*fin[0]*win[2]+0.7071067811865475*win[0]*fin[2])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
  atomicAdd(&out[2], (fin[1]*win[3]+win[1]*fin[3]+fin[0]*win[2]+win[0]*fin[2])*subvol); 
#else
  out[0] += (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
  out[2] += (fin[1]*win[3]+win[1]*fin[3]+fin[0]*win[2]+win[0]*fin[2])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol); 
  atomicAdd(&out[1], (0.7071067811865475*fin[2]*win[3]+0.7071067811865475*win[2]*fin[3]+0.7071067811865475*fin[0]*win[1]+0.7071067811865475*win[0]*fin[1])*subvol); 
#else
  out[0] += (0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[2]*win[3]+0.7071067811865475*win[2]*fin[3]+0.7071067811865475*fin[0]*win[1]+0.7071067811865475*win[0]*fin[1])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
  atomicAdd(&out[1], (fin[2]*win[3]+win[2]*fin[3]+fin[0]*win[1]+win[0]*fin[1])*subvol); 
#else
  out[0] += (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
  out[1] += (fin[2]*win[3]+win[2]*fin[3]+fin[0]*win[1]+win[0]*fin[1])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
#else
  out[0] += (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (2.0*fin[3]*win[3]+2.0*fin[2]*win[2]+2.0*fin[1]*win[1]+2.0*fin[0]*win[0])*subvol); 
#else
  out[0] += (2.0*fin[3]*win[3]+2.0*fin[2]*win[2]+2.0*fin[1]*win[1]+2.0*fin[0]*win[0])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol); 
  atomicAdd(&out[1], (0.632455532033676*fin[3]*win[7]+0.632455532033676*win[3]*fin[7]+0.7071067811865475*fin[4]*win[6]+0.7071067811865475*win[4]*fin[6]+0.6324555320336759*fin[2]*win[5]+0.6324555320336759*win[2]*fin[5]+0.7071067811865475*fin[1]*win[3]+0.7071067811865475*win[1]*fin[3]+0.7071067811865475*fin[0]*win[2]+0.7071067811865475*win[0]*fin[2])*subvol); 
  atomicAdd(&out[2], (0.45175395145262565*fin[7]*win[7]+0.7071067811865475*fin[1]*win[7]+0.7071067811865475*win[1]*fin[7]+0.6324555320336759*fin[6]*win[6]+0.45175395145262565*fin[5]*win[5]+0.7071067811865475*fin[0]*win[5]+0.7071067811865475*win[0]*fin[5]+0.6324555320336759*fin[3]*win[3]+0.6324555320336759*fin[2]*win[2])*subvol); 
#else
  out[0] += (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.632455532033676*fin[3]*win[7]+0.632455532033676*win[3]*fin[7]+0.7071067811865475*fin[4]*win[6]+0.7071067811865475*win[4]*fin[6]+0.6324555320336759*fin[2]*win[5]+0.6324555320336759*win[2]*fin[5]+0.7071067811865475*fin[1]*win[3]+0.7071067811865475*win[1]*fin[3]+0.7071067811865475*fin[0]*win[2]+0.7071067811865475*win[0]*fin[2])*subvol; 
  out[2] += (0.45175395145262565*fin[7]*win[7]+0.7071067811865475*fin[1]*win[7]+0.7071067811865475*win[1]*fin[7]+0.6324555320336759*fin[6]*win[6]+0.45175395145262565*fin[5]*win[5]+0.7071067811865475*fin[0]*win[5]+0.7071067811865475*win[0]*fin[5]+0.6324555320336759*fin[3]*win[3]+0.6324555320336759*fin[2]*win[2])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgx_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
  atomicAdd(&out[2], (0.8944271909999161*fin[3]*win[7]+0.8944271909999161*win[3]*fin[7]+1.0000000000000002*fin[4]*win[6]+1.0000000000000002*win[4]*fin[6]+0.8944271909999159*fin[2]*win[5]+0.8944271909999159*win[2]*fin[5]+fin[1]*win[3]+win[1]*fin[3]+fin[0]*win[2]+win[0]*fin[2])*subvol); 
  atomicAdd(&out[5], (0.6388765649999399*fin[7]*win[7]+1.0000000000000002*fin[1]*win[7]+1.0000000000000002*win[1]*fin[7]+0.8944271909999159*fin[6]*win[6]+0.6388765649999399*fin[5]*win[5]+fin[0]*win[5]+win[0]*fin[5]+0.8944271909999159*fin[3]*win[3]+0.8944271909999159*fin[2]*win[2])*subvol); 
#else
  out[0] += (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
  out[2] += (0.8944271909999161*fin[3]*win[7]+0.8944271909999161*win[3]*fin[7]+1.0000000000000002*fin[4]*win[6]+1.0000000000000002*win[4]*fin[6]+0.8944271909999159*fin[2]*win[5]+0.8944271909999159*win[2]*fin[5]+fin[1]*win[3]+win[1]*fin[3]+fin[0]*win[2]+win[0]*fin[2])*subvol; 
  out[5] += (0.6388765649999399*fin[7]*win[7]+1.0000000000000002*fin[1]*win[7]+1.0000000000000002*win[1]*fin[7]+0.8944271909999159*fin[6]*win[6]+0.6388765649999399*fin[5]*win[5]+fin[0]*win[5]+win[0]*fin[5]+0.8944271909999159*fin[3]*win[3]+0.8944271909999159*fin[2]*win[2])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol); 
  atomicAdd(&out[1], (0.7071067811865475*fin[5]*win[7]+0.7071067811865475*win[5]*fin[7]+0.632455532033676*fin[3]*win[6]+0.632455532033676*win[3]*fin[6]+0.6324555320336759*fin[1]*win[4]+0.6324555320336759*win[1]*fin[4]+0.7071067811865475*fin[2]*win[3]+0.7071067811865475*win[2]*fin[3]+0.7071067811865475*fin[0]*win[1]+0.7071067811865475*win[0]*fin[1])*subvol); 
  atomicAdd(&out[2], (0.6324555320336759*fin[7]*win[7]+0.45175395145262565*fin[6]*win[6]+0.7071067811865475*fin[2]*win[6]+0.7071067811865475*win[2]*fin[6]+0.45175395145262565*fin[4]*win[4]+0.7071067811865475*fin[0]*win[4]+0.7071067811865475*win[0]*fin[4]+0.6324555320336759*fin[3]*win[3]+0.6324555320336759*fin[1]*win[1])*subvol); 
#else
  out[0] += (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[5]*win[7]+0.7071067811865475*win[5]*fin[7]+0.632455532033676*fin[3]*win[6]+0.632455532033676*win[3]*fin[6]+0.6324555320336759*fin[1]*win[4]+0.6324555320336759*win[1]*fin[4]+0.7071067811865475*fin[2]*win[3]+0.7071067811865475*win[2]*fin[3]+0.7071067811865475*fin[0]*win[1]+0.7071067811865475*win[0]*fin[1])*subvol; 
  out[2] += (0.6324555320336759*fin[7]*win[7]+0.45175395145262565*fin[6]*win[6]+0.7071067811865475*fin[2]*win[6]+0.7071067811865475*win[2]*fin[6]+0.45175395145262565*fin[4]*win[4]+0.7071067811865475*fin[0]*win[4]+0.7071067811865475*win[0]*fin[4]+0.6324555320336759*fin[3]*win[3]+0.6324555320336759*fin[1]*win[1])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
  atomicAdd(&out[1], (1.0000000000000002*fin[5]*win[7]+1.0000000000000002*win[5]*fin[7]+0.8944271909999161*fin[3]*win[6]+0.8944271909999161*win[3]*fin[6]+0.8944271909999159*fin[1]*win[4]+0.8944271909999159*win[1]*fin[4]+fin[2]*win[3]+win[2]*fin[3]+fin[0]*win[1]+win[0]*fin[1])*subvol); 
  atomicAdd(&out[4], (0.8944271909999159*fin[7]*win[7]+0.6388765649999399*fin[6]*win[6]+1.0000000000000002*fin[2]*win[6]+1.0000000000000002*win[2]*fin[6]+0.6388765649999399*fin[4]*win[4]+fin[0]*win[4]+win[0]*fin[4]+0.8944271909999159*fin[3]*win[3]+0.8944271909999159*fin[1]*win[1])*subvol); 
#else
  out[0] += (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
  out[1] += (1.0000000000000002*fin[5]*win[7]+1.0000000000000002*win[5]*fin[7]+0.8944271909999161*fin[3]*win[6]+0.8944271909999161*win[3]*fin[6]+0.8944271909999159*fin[1]*win[4]+0.8944271909999159*win[1]*fin[4]+fin[2]*win[3]+win[2]*fin[3]+fin[0]*win[1]+win[0]*fin[1])*subvol; 
  out[4] += (0.8944271909999159*fin[7]*win[7]+0.6388765649999399*fin[6]*win[6]+1.0000000000000002*fin[2]*win[6]+1.0000000000000002*win[2]*fin[6]+0.6388765649999399*fin[4]*win[4]+fin[0]*win[4]+win[0]*fin[4]+0.8944271909999159*fin[3]*win[3]+0.8944271909999159*fin[1]*win[1])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol); 
#else
  out[0] += (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
#endif
} 

GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgxy_fulldim(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
#ifdef __CUDA_ARCH__
  atomicAdd(&out[0], (2.0*fin[7]*win[7]+2.0*fin[6]*win[6]+2.0*fin[5]*win[5]+2.0*fin[4]*win[4]+2.0*fin[3]*win[3]+2.0*fin[2]*win[2]+2.0*fin[1]*win[1]+2.0*fin[0]*win[0])*subvol); 
#else
  out[0] += (2.0*fin[7]*win[7]+2.0*fin[6]*win[6]+2.0*fin[5]*win[5]+2.0*fin[4]*win[4]+2.0*fin[3]*win[3]+2.0*fin[2]*win[2]+2.0*fin[1]*win[1]+2.0*fin[0]*win[0])*subvol; 
#endif
} 

