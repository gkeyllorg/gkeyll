#pragma once

#include <gkyl_array.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
#include <gkyl_ref_count.h>


#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#ifdef GKYL_HAVE_CUDA
#if !defined(CUBLAS_V2_H_)
#include <cublas_v2.h>
#endif
#endif


struct gkyl_mat_mm_array_mem {

  // info for alpha*matrix_multiplication(A,B) + Beta*C = C 
  // using mat_mm_array
  bool on_gpu; // flag to indicate if we are on GPU
  double alpha;
  double beta;
  enum gkyl_mat_trans transa;
  enum gkyl_mat_trans transb;
  struct gkyl_mat *A;

#ifdef GKYL_HAVE_CUDA
  cublasHandle_t cuh; // cublas handle
#endif  
};