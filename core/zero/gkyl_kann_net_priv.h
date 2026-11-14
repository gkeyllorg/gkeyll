#pragma once

#include <gkyl_kann_net.h>
#include <gkyl_alloc_flags_priv.h>

#ifdef GKYL_HAVE_CUDA
#include <cublas_v2.h>
#endif

struct gkyl_kann_net {
  kann_t *ann; // underlying KANN network

  int n_var; // number of trainable variable floats
  int n_const; // number of constant floats
  int n_in; // input dimension (per sample)
  int n_out; // output dimension (per sample)

  uint32_t flags;
  struct gkyl_ref_count ref_count;

#ifdef GKYL_HAVE_CUDA
  cublasHandle_t cublas_h; // persistent cuBLAS handle (GPU only)
  struct kann_cu_graph *cg; // flattened graph on device (GPU only)
#endif
};
