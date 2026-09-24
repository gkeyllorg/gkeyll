// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_binop_div_ser.h>
#include <gkyl_binop_mul_ser.h>
#include <gkyl_binop_cross_mul_ser.h>
#include <gkyl_binop_cross_mul_tensor.h>
#include <gkyl_binop_cross_mul_hyb.h>
#include <gkyl_binop_cross_mul_gkhyb.h>
#include <gkyl_basis_ser_1x_p1_inv.h>
#include <gkyl_basis_ser_2x_p1_inv.h>
#include <gkyl_basis_ser_3x_p1_inv.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

enum gkyl_dg_op { GKYL_DG_OP_MEAN, GKYL_DG_OP_MEAN_L2 };

// Memory for use in the bin ops
struct gkyl_dg_bin_op_mem {
  bool on_gpu; // flag to indicate if we are on GPU  
  size_t batch_sz; // number of elements in batch
  size_t nrows, ncols; // number of rows and colsx
  struct gkyl_nmat *As, *xs; // data for matrices needed in division
  gkyl_nmat_mem *lu_mem; // data for use in LU solve
};

// Function pointer type for multiplication
typedef void (*mul_op_t)(const double *f, const double *g, double *fg);
typedef void (*mul_accumulate_op_t)(double a, const double *f, const double *g, double *fg);
typedef void (*mul_accumulate_comp_par_op_t)(double a, const double *f, const double *g, double *fg, int linc2);

typedef struct gkyl_kern_op_count (*mul_op_count_t)(void);

// Function pointer type for setting matrices for division
typedef void (*div_set_op_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);

// Function pointer type for inverting a DG field.
typedef void (*inv_op_t)(const double *A, double *A_inv);

// for use in kernel tables
typedef struct { mul_op_t kernels[4]; } mul_op_kern_list;
typedef struct { mul_accumulate_op_t kernels[4]; } mul_accumulate_op_kern_list;
typedef struct { mul_accumulate_comp_par_op_t kernels[4]; } mul_accumulate_comp_par_op_kern_list;
typedef struct { mul_op_kern_list list[3]; } cross_mul_op_kern_list;
typedef struct { mul_accumulate_op_kern_list list[3]; } cross_mul_accumulate_op_kern_list;
typedef struct { mul_accumulate_comp_par_op_kern_list list[3]; } cross_mul_accumulate_comp_par_op_kern_list;
typedef struct { mul_op_count_t kernels[4]; } mul_op_count_kern_list;
typedef struct { div_set_op_t kernels[4]; } div_set_op_kern_list;
typedef struct { inv_op_t kernels[4]; } inv_op_kern_list;

// Serendipity multiplication kernels
GKYL_CU_D
static const mul_op_kern_list ser_mul_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_mul_1d_ser_p0, binop_mul_1d_ser_p1, binop_mul_1d_ser_p2, binop_mul_1d_ser_p3 },
  { binop_mul_2d_ser_p0, binop_mul_2d_ser_p1, binop_mul_2d_ser_p2, binop_mul_2d_ser_p3 },
  { binop_mul_3d_ser_p0, binop_mul_3d_ser_p1, binop_mul_3d_ser_p2, binop_mul_3d_ser_p3 }
};

// Tensor multiplication kernels
GKYL_CU_D
static const mul_op_kern_list ten_mul_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_mul_1d_ser_p0, binop_mul_1d_ser_p1, binop_mul_1d_ser_p2, binop_mul_1d_ser_p3 },
  { binop_mul_2d_ser_p0, binop_mul_2d_ser_p1, binop_mul_2d_tensor_p2, NULL },
  { binop_mul_3d_ser_p0, binop_mul_3d_ser_p1, binop_mul_3d_tensor_p2, NULL }
};

// Serendipity conf*phase multiplication kernels
GKYL_CU_D
static const cross_mul_op_kern_list ser_cross_mul_list[] = {
  // pdim=2
  { .list = {{ NULL, binop_cross_mul_1d_2d_ser_p1, binop_cross_mul_1d_2d_ser_p2, binop_cross_mul_1d_2d_ser_p3 },
             { NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=3
  { .list = {{ NULL, binop_cross_mul_1d_3d_ser_p1, binop_cross_mul_1d_3d_ser_p2, binop_cross_mul_1d_3d_ser_p3 },
             { NULL, binop_cross_mul_2d_3d_ser_p1, binop_cross_mul_2d_3d_ser_p2, binop_cross_mul_2d_3d_ser_p3 },
             { NULL, NULL, NULL, NULL },} },
  // pdim=4
  { .list = {{ NULL, binop_cross_mul_1d_4d_ser_p1, binop_cross_mul_1d_4d_ser_p2, binop_cross_mul_1d_4d_ser_p3 },
             { NULL, binop_cross_mul_2d_4d_ser_p1, binop_cross_mul_2d_4d_ser_p2, binop_cross_mul_2d_4d_ser_p3 },
             { NULL, binop_cross_mul_3d_4d_ser_p1, binop_cross_mul_3d_4d_ser_p2, binop_cross_mul_3d_4d_ser_p3 },} },
  // pdim=5
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_2d_5d_ser_p1, binop_cross_mul_2d_5d_ser_p2, NULL },
             { NULL, binop_cross_mul_3d_5d_ser_p1, binop_cross_mul_3d_5d_ser_p2, NULL },} },
  // pdim=6
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_3d_6d_ser_p1, NULL, NULL },} },
};

// Serendipity conf*phase multiplication with accumulation to output kernels
GKYL_CU_D
static const cross_mul_accumulate_op_kern_list ser_cross_mul_accumulate_list[] = {
  // pdim=2
  { .list = {{ NULL, binop_cross_mul_accumulate_1d_2d_ser_p1, binop_cross_mul_accumulate_1d_2d_ser_p2, binop_cross_mul_accumulate_1d_2d_ser_p3 },
             { NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=3
  { .list = {{ NULL, binop_cross_mul_accumulate_1d_3d_ser_p1, binop_cross_mul_accumulate_1d_3d_ser_p2, binop_cross_mul_accumulate_1d_3d_ser_p3 },
             { NULL, binop_cross_mul_accumulate_2d_3d_ser_p1, binop_cross_mul_accumulate_2d_3d_ser_p2, binop_cross_mul_accumulate_2d_3d_ser_p3 },
             { NULL, NULL, NULL, NULL },} },
  // pdim=4
  { .list = {{ NULL, binop_cross_mul_accumulate_1d_4d_ser_p1, binop_cross_mul_accumulate_1d_4d_ser_p2, binop_cross_mul_accumulate_1d_4d_ser_p3 },
             { NULL, binop_cross_mul_accumulate_2d_4d_ser_p1, binop_cross_mul_accumulate_2d_4d_ser_p2, binop_cross_mul_accumulate_2d_4d_ser_p3 },
             { NULL, binop_cross_mul_accumulate_3d_4d_ser_p1, binop_cross_mul_accumulate_3d_4d_ser_p2, binop_cross_mul_accumulate_3d_4d_ser_p3 },} },
  // pdim=5
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_2d_5d_ser_p1, binop_cross_mul_accumulate_2d_5d_ser_p2, NULL },
             { NULL, binop_cross_mul_accumulate_3d_5d_ser_p1, binop_cross_mul_accumulate_3d_5d_ser_p2, NULL },} },
  // pdim=6
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_3d_6d_ser_p1, NULL, NULL },} },
};

// Serendipity conf*phase multiplication with accumulation to output kernels and parallelization over components
GKYL_CU_D
static const cross_mul_accumulate_comp_par_op_kern_list ser_cross_mul_accumulate_comp_par_list[] = {
  // pdim=2
  { .list = {{ NULL, binop_cross_mul_accumulate_comp_par_1d_2d_ser_p1, binop_cross_mul_accumulate_comp_par_1d_2d_ser_p2, binop_cross_mul_accumulate_comp_par_1d_2d_ser_p3 },
             { NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=3
  { .list = {{ NULL, binop_cross_mul_accumulate_comp_par_1d_3d_ser_p1, binop_cross_mul_accumulate_comp_par_1d_3d_ser_p2, binop_cross_mul_accumulate_comp_par_1d_3d_ser_p3 },
             { NULL, binop_cross_mul_accumulate_comp_par_2d_3d_ser_p1, binop_cross_mul_accumulate_comp_par_2d_3d_ser_p2, binop_cross_mul_accumulate_comp_par_2d_3d_ser_p3 },
             { NULL, NULL, NULL, NULL },} },
  // pdim=4
  { .list = {{ NULL, binop_cross_mul_accumulate_comp_par_1d_4d_ser_p1, binop_cross_mul_accumulate_comp_par_1d_4d_ser_p2, binop_cross_mul_accumulate_comp_par_1d_4d_ser_p3 },
             { NULL, binop_cross_mul_accumulate_comp_par_2d_4d_ser_p1, binop_cross_mul_accumulate_comp_par_2d_4d_ser_p2, binop_cross_mul_accumulate_comp_par_2d_4d_ser_p3 },
             { NULL, binop_cross_mul_accumulate_comp_par_3d_4d_ser_p1, binop_cross_mul_accumulate_comp_par_3d_4d_ser_p2, binop_cross_mul_accumulate_comp_par_3d_4d_ser_p3 },} },
  // pdim=5
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_2d_5d_ser_p1, binop_cross_mul_accumulate_comp_par_2d_5d_ser_p2, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_3d_5d_ser_p1, binop_cross_mul_accumulate_comp_par_3d_5d_ser_p2, NULL },} },
  // pdim=6
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_3d_6d_ser_p1, NULL, NULL },} },
};

// Tensor conf*phase multiplication kernels
GKYL_CU_D
static const cross_mul_op_kern_list ten_cross_mul_list[] = {
  // pdim=2
  { .list = {{ NULL, binop_cross_mul_1d_2d_ser_p1, binop_cross_mul_1d_2d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=3
  { .list = {{ NULL, binop_cross_mul_1d_3d_ser_p1, binop_cross_mul_1d_3d_tensor_p2, NULL },
             { NULL, binop_cross_mul_2d_3d_ser_p1, binop_cross_mul_2d_3d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=4
  { .list = {{ NULL, binop_cross_mul_1d_4d_ser_p1, binop_cross_mul_1d_4d_tensor_p2, NULL },
             { NULL, binop_cross_mul_2d_4d_ser_p1, binop_cross_mul_2d_4d_tensor_p2, NULL },
             { NULL, binop_cross_mul_3d_4d_ser_p1, NULL, NULL },} },
  // pdim=5
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_2d_5d_ser_p1, binop_cross_mul_2d_5d_tensor_p2, NULL },
             { NULL, binop_cross_mul_3d_5d_ser_p1, NULL, NULL },} },
  // pdim=6
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_3d_6d_ser_p1, NULL, NULL },} },
};

// Tensor conf*phase multiplication with accumulation to output kernels 
GKYL_CU_D
static const cross_mul_accumulate_op_kern_list ten_cross_mul_accumulate_list[] = {
  // pdim=2
  { .list = {{ NULL, binop_cross_mul_accumulate_1d_2d_ser_p1, binop_cross_mul_accumulate_1d_2d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=3
  { .list = {{ NULL, binop_cross_mul_accumulate_1d_3d_ser_p1, binop_cross_mul_accumulate_1d_3d_tensor_p2, NULL },
             { NULL, binop_cross_mul_accumulate_2d_3d_ser_p1, binop_cross_mul_accumulate_2d_3d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=4
  { .list = {{ NULL, binop_cross_mul_accumulate_1d_4d_ser_p1, binop_cross_mul_accumulate_1d_4d_tensor_p2, NULL },
             { NULL, binop_cross_mul_accumulate_2d_4d_ser_p1, binop_cross_mul_accumulate_2d_4d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=5
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_2d_5d_ser_p1, binop_cross_mul_accumulate_2d_5d_tensor_p2, NULL },
             { NULL, binop_cross_mul_accumulate_3d_5d_ser_p1, NULL, NULL },} },
  // pdim=6
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_3d_6d_ser_p1, NULL, NULL },} },
};

// Tensor conf*phase multiplication with accumulation to output kernels and parallelization over components
GKYL_CU_D
static const cross_mul_accumulate_comp_par_op_kern_list ten_cross_mul_accumulate_comp_par_list[] = {
  // pdim=2
  { .list = {{ NULL, binop_cross_mul_accumulate_comp_par_1d_2d_ser_p1, binop_cross_mul_accumulate_comp_par_1d_2d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=3
  { .list = {{ NULL, binop_cross_mul_accumulate_comp_par_1d_3d_ser_p1, binop_cross_mul_accumulate_comp_par_1d_3d_tensor_p2, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_2d_3d_ser_p1, binop_cross_mul_accumulate_comp_par_2d_3d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=4
  { .list = {{ NULL, binop_cross_mul_accumulate_comp_par_1d_4d_ser_p1, binop_cross_mul_accumulate_comp_par_1d_4d_tensor_p2, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_2d_4d_ser_p1, binop_cross_mul_accumulate_comp_par_2d_4d_tensor_p2, NULL },
             { NULL, NULL, NULL, NULL },} },
  // pdim=5
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_2d_5d_ser_p1, binop_cross_mul_accumulate_comp_par_2d_5d_tensor_p2, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_3d_5d_ser_p1, NULL, NULL },} },
  // pdim=6
  { .list = {{ NULL, NULL, NULL, NULL },
             { NULL, NULL, NULL, NULL },
             { NULL, binop_cross_mul_accumulate_comp_par_3d_6d_ser_p1, NULL, NULL },} },
};

// Hybrid basis conf*phase multiplication kernels 
GKYL_CU_D
static const mul_op_kern_list hyb_cross_mul_list[] = {
  { binop_cross_mul_1x1v_hyb_p1, binop_cross_mul_1x2v_hyb_p1, binop_cross_mul_1x3v_hyb_p1 },
  { binop_cross_mul_2x1v_hyb_p1, binop_cross_mul_2x2v_hyb_p1, binop_cross_mul_2x3v_hyb_p1 },
  { binop_cross_mul_3x1v_hyb_p1, binop_cross_mul_3x2v_hyb_p1, binop_cross_mul_3x3v_hyb_p1 },
};

// Hybrid basis conf*phase multiplication with accumulation to output kernels 
GKYL_CU_D    
static const mul_accumulate_op_kern_list hyb_cross_mul_accumulate_list[] = {
  { binop_cross_mul_accumulate_1x1v_hyb_p1, binop_cross_mul_accumulate_1x2v_hyb_p1, binop_cross_mul_accumulate_1x3v_hyb_p1 },
  { binop_cross_mul_accumulate_2x1v_hyb_p1, binop_cross_mul_accumulate_2x2v_hyb_p1, binop_cross_mul_accumulate_2x3v_hyb_p1 },
  { binop_cross_mul_accumulate_3x1v_hyb_p1, binop_cross_mul_accumulate_3x2v_hyb_p1, binop_cross_mul_accumulate_3x3v_hyb_p1 },
};

// Hybrid basis conf*phase multiplication with accumulation to output kernels and parallelization over components
GKYL_CU_D    
static const mul_accumulate_comp_par_op_kern_list hyb_cross_mul_accumulate_comp_par_list[] = {
  { binop_cross_mul_accumulate_comp_par_1x1v_hyb_p1, binop_cross_mul_accumulate_comp_par_1x2v_hyb_p1, binop_cross_mul_accumulate_comp_par_1x3v_hyb_p1 },
  { binop_cross_mul_accumulate_comp_par_2x1v_hyb_p1, binop_cross_mul_accumulate_comp_par_2x2v_hyb_p1, binop_cross_mul_accumulate_comp_par_2x3v_hyb_p1 },
  { binop_cross_mul_accumulate_comp_par_3x1v_hyb_p1, binop_cross_mul_accumulate_comp_par_3x2v_hyb_p1, binop_cross_mul_accumulate_comp_par_3x3v_hyb_p1 },
};

// GK hybrid basis conf*phase multiplication kernels 
GKYL_CU_D
static const mul_op_kern_list gkhyb_cross_mul_list[] = {
  { binop_cross_mul_1x1v_gkhyb_p1, binop_cross_mul_1x2v_gkhyb_p1 },
  { NULL, binop_cross_mul_2x2v_gkhyb_p1 },
  { NULL, binop_cross_mul_3x2v_gkhyb_p1 },
};

// GK hybrid basis conf*phase multiplication with accumulation to output kernels 
GKYL_CU_D
static const mul_accumulate_op_kern_list gkhyb_cross_mul_accumulate_list[] = {
  { binop_cross_mul_accumulate_1x1v_gkhyb_p1, binop_cross_mul_accumulate_1x2v_gkhyb_p1 },
  { NULL, binop_cross_mul_accumulate_2x2v_gkhyb_p1 },
  { NULL, binop_cross_mul_accumulate_3x2v_gkhyb_p1 },
};

// GK hybrid basis conf*phase multiplication with accumulation to output kernels and parallelization over components
GKYL_CU_D
static const mul_accumulate_comp_par_op_kern_list gkhyb_cross_mul_accumulate_comp_par_list[] = {
  { binop_cross_mul_accumulate_comp_par_1x1v_gkhyb_p1, binop_cross_mul_accumulate_comp_par_1x2v_gkhyb_p1 },
  { NULL, binop_cross_mul_accumulate_comp_par_2x2v_gkhyb_p1 },
  { NULL, binop_cross_mul_accumulate_comp_par_3x2v_gkhyb_p1 },
};

static const mul_op_count_kern_list ser_mul_op_count_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { op_count_binop_mul_1d_ser_p0, op_count_binop_mul_1d_ser_p1, op_count_binop_mul_1d_ser_p2, op_count_binop_mul_1d_ser_p3 },
  { op_count_binop_mul_2d_ser_p0, op_count_binop_mul_2d_ser_p1, op_count_binop_mul_2d_ser_p2, op_count_binop_mul_2d_ser_p3 },
  { op_count_binop_mul_3d_ser_p0, op_count_binop_mul_3d_ser_p1, op_count_binop_mul_3d_ser_p2, op_count_binop_mul_3d_ser_p3 }
};

// Serendipity division kernels
GKYL_CU_D
static const div_set_op_kern_list ser_div_set_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_div_set_1d_ser_p0, binop_div_set_1d_ser_p1, binop_div_set_1d_ser_p2, binop_div_set_1d_ser_p3 },
  { binop_div_set_2d_ser_p0, binop_div_set_2d_ser_p1, binop_div_set_2d_ser_p2, binop_div_set_2d_ser_p3 },
  { binop_div_set_3d_ser_p0, binop_div_set_3d_ser_p1, binop_div_set_3d_ser_p2, binop_div_set_3d_ser_p3 } 
};

// Tensor division kernels
GKYL_CU_D
static const div_set_op_kern_list ten_div_set_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_div_set_1d_ser_p0, binop_div_set_1d_ser_p1, binop_div_set_1d_ser_p2, binop_div_set_1d_ser_p3 },
  { binop_div_set_2d_ser_p0, binop_div_set_2d_ser_p1, binop_div_set_2d_tensor_p2, NULL },
  { binop_div_set_3d_ser_p0, binop_div_set_3d_ser_p1, binop_div_set_3d_tensor_p2, NULL } 
};

// Serendipity inv kernels
GKYL_CU_D
static const inv_op_kern_list ser_inv_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, ser_1x_p1_inv, NULL, NULL },
  { NULL, ser_2x_p1_inv, NULL, NULL },
  { NULL, ser_3x_p1_inv, NULL, NULL } 
};

GKYL_CU_D
static mul_op_t
choose_ser_mul_kern(int dim, int poly_order)
{
  assert(dim < 4);
  return ser_mul_list[dim].kernels[poly_order];
}

GKYL_CU_D
static mul_op_t
choose_ten_mul_kern(int dim, int poly_order)
{
  assert(dim < 4);
  return ten_mul_list[dim].kernels[poly_order];
}

GKYL_CU_D
static mul_op_t
choose_mul_conf_phase_kern(enum gkyl_basis_type btype, int cdim, int vdim, int poly_order)
{
  int pdim = cdim+vdim;
  
  switch (btype) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_cross_mul_list[pdim-2].list[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      return hyb_cross_mul_list[cdim-1].kernels[vdim-1];
      break;
    case GKYL_BASIS_MODAL_GKHYBRID:
      return gkhyb_cross_mul_list[cdim-1].kernels[vdim-1];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_cross_mul_list[pdim-2].list[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_D
static mul_accumulate_op_t
choose_mul_conf_phase_accumulate_kern(enum gkyl_basis_type btype, int cdim, int vdim, int poly_order)
{
  int pdim = cdim+vdim;
  
  switch (btype) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_cross_mul_accumulate_list[pdim-2].list[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      return hyb_cross_mul_accumulate_list[cdim-1].kernels[vdim-1];
      break;
    case GKYL_BASIS_MODAL_GKHYBRID:
      return gkhyb_cross_mul_accumulate_list[cdim-1].kernels[vdim-1];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_cross_mul_accumulate_list[pdim-2].list[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_D
static mul_accumulate_comp_par_op_t
choose_mul_conf_phase_accumulate_comp_par_kern(enum gkyl_basis_type btype, int cdim, int vdim, int poly_order)
{
  int pdim = cdim+vdim;
  
  switch (btype) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_cross_mul_accumulate_comp_par_list[pdim-2].list[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      return hyb_cross_mul_accumulate_comp_par_list[cdim-1].kernels[vdim-1];
      break;
    case GKYL_BASIS_MODAL_GKHYBRID:
      return gkhyb_cross_mul_accumulate_comp_par_list[cdim-1].kernels[vdim-1];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_cross_mul_accumulate_comp_par_list[pdim-2].list[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

static mul_op_count_t
choose_ser_mul_op_count_kern(int dim, int poly_order)
{
  return ser_mul_op_count_list[dim].kernels[poly_order];
}

GKYL_CU_D
static div_set_op_t
choose_ser_div_set_kern(int dim, int poly_order)
{
  return ser_div_set_list[dim].kernels[poly_order];
}

GKYL_CU_D
static div_set_op_t
choose_ten_div_set_kern(int dim, int poly_order)
{
  return ten_div_set_list[dim].kernels[poly_order];
}

GKYL_CU_D
static inv_op_t
choose_ser_inv_kern(int dim, int poly_order)
{
  return ser_inv_list[dim].kernels[poly_order];
}

GKYL_CU_D
static inline double
dg_cell_mean(int nc, const double *f)
{
  return f[0];
}

GKYL_CU_D
static inline double
dg_cell_mean_l2(int nb, const double *f)
{
  double sum = 0.0;
  for (int i=0; i<nb; ++i)
    sum += f[i]*f[i];
  return sum;
}

typedef double (*dp_op_t)(int nb, const double *f);

GKYL_CU_D
static dp_op_t
dg_get_op_func(enum gkyl_dg_op op)
{
  if (op == GKYL_DG_OP_MEAN)
    return dg_cell_mean;
  return dg_cell_mean_l2;
}

void gkyl_dg_calc_op_range(const struct gkyl_basis *basis, int c_oop,
  struct gkyl_array *out, int c_iop,
  const struct gkyl_array *iop,
  struct gkyl_range range, enum gkyl_dg_op op);

void gkyl_dg_calc_op_range_cu(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop,
  struct gkyl_range range, enum gkyl_dg_op op);

