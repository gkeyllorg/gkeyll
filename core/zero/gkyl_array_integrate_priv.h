#pragma once

// Private header for array_integrate updater, not for direct use by user.

#include <gkyl_util.h>
#include <math.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_integrate_kernels.h>
#include <assert.h>

// Function pointer type for array_integrate kernels.
typedef void (*array_integrate_t)(double *dxSq, double vol, int num_comp,
  int num_basis, const double *weight, const double *fIn, double *out);

// For use in kernel tables.
typedef struct { array_integrate_t kernels[2]; } array_integrate_none_kern_list;
typedef struct { array_integrate_t kernels[2]; } array_integrate_abs_kern_list;
typedef struct { array_integrate_t kernels[2]; } array_integrate_sq_kern_list;
typedef struct { array_integrate_t kernels[2]; } array_integrate_sq_weighted_kern_list;
typedef struct { array_integrate_t kernels[2]; } array_integrate_gradsq_kern_list;
typedef struct { array_integrate_t kernels[2]; } array_integrate_gradperpsq_kern_list;
typedef struct { array_integrate_t kernels[2]; } array_integrate_epsgradperpsq_kern_list;

GKYL_CU_D
static const array_integrate_none_kern_list gkyl_array_integrate_none_ker_list_ser[] = {
  {gkyl_array_integrate_op_none_1x_ser_p1, gkyl_array_integrate_op_none_1x_ser_p2},
  {gkyl_array_integrate_op_none_2x_ser_p1, gkyl_array_integrate_op_none_2x_ser_p2},
  {gkyl_array_integrate_op_none_3x_ser_p1, gkyl_array_integrate_op_none_3x_ser_p2},
};

GKYL_CU_D
static const array_integrate_abs_kern_list gkyl_array_integrate_abs_ker_list_ser[] = {
  {gkyl_array_integrate_op_abs_1x_ser_p1, gkyl_array_integrate_op_abs_1x_ser_p2},
  {gkyl_array_integrate_op_abs_2x_ser_p1, gkyl_array_integrate_op_abs_2x_ser_p2},
  {gkyl_array_integrate_op_abs_3x_ser_p1, gkyl_array_integrate_op_abs_3x_ser_p2},
};

GKYL_CU_D
static const array_integrate_sq_kern_list gkyl_array_integrate_sq_ker_list_ser[] = {
  {gkyl_array_integrate_op_sq_1x_ser_p1, gkyl_array_integrate_op_sq_1x_ser_p2},
  {gkyl_array_integrate_op_sq_2x_ser_p1, gkyl_array_integrate_op_sq_2x_ser_p2},
  {gkyl_array_integrate_op_sq_3x_ser_p1, gkyl_array_integrate_op_sq_3x_ser_p2},
};

GKYL_CU_D
static const array_integrate_sq_weighted_kern_list gkyl_array_integrate_sq_weighted_ker_list_ser[] = {
  {gkyl_array_integrate_op_sq_weighted_1x_ser_p1, gkyl_array_integrate_op_sq_weighted_1x_ser_p2},
  {gkyl_array_integrate_op_sq_weighted_2x_ser_p1, gkyl_array_integrate_op_sq_weighted_2x_ser_p2},
  {gkyl_array_integrate_op_sq_weighted_3x_ser_p1, gkyl_array_integrate_op_sq_weighted_3x_ser_p2},
};

GKYL_CU_D
static const array_integrate_sq_weighted_kern_list gkyl_array_integrate_sq_weighted_ker_list_gkhyb[] = {
  {NULL, NULL},
  {gkyl_array_integrate_op_sq_weighted_1x1v_gkhyb_p1, NULL},
  {gkyl_array_integrate_op_sq_weighted_1x2v_gkhyb_p1, NULL},
  {gkyl_array_integrate_op_sq_weighted_2x2v_gkhyb_p1, NULL},
  {gkyl_array_integrate_op_sq_weighted_3x2v_gkhyb_p1, NULL},
};

GKYL_CU_D
static const array_integrate_gradsq_kern_list gkyl_array_integrate_gradsq_ker_list[] = {
  {gkyl_array_integrate_op_grad_sq_1x_ser_p1, gkyl_array_integrate_op_grad_sq_1x_ser_p2},
  {gkyl_array_integrate_op_grad_sq_2x_ser_p1, gkyl_array_integrate_op_grad_sq_2x_ser_p2},
  {gkyl_array_integrate_op_grad_sq_3x_ser_p1, NULL},
};

GKYL_CU_D
static const array_integrate_gradperpsq_kern_list gkyl_array_integrate_gradperpsq_ker_list[] = {
  {NULL, NULL},
  {gkyl_array_integrate_op_grad_sq_2x_ser_p1, gkyl_array_integrate_op_grad_sq_2x_ser_p2},
  {gkyl_array_integrate_op_gradperp_sq_3x_ser_p1, gkyl_array_integrate_op_gradperp_sq_3x_ser_p2},
};

GKYL_CU_D
static const array_integrate_epsgradperpsq_kern_list gkyl_array_integrate_epsgradperpsq_ker_list[] = {
  {NULL, NULL},
  {gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p1, gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p2},
  {gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p1, gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p2},
};

// Primary struct in this updater.
struct gkyl_array_integrate {
  enum gkyl_array_integrate_op op; // Operator type.
  int num_basis, num_comp;
  bool use_gpu;
  double dxSq[GKYL_MAX_DIM];
  double vol;  // Single-cell volume factor.
  array_integrate_t kernel;  // Single cell integration kernel.
  struct gkyl_array_integrate *on_dev;  // Pointer to itself on device.
};

GKYL_CU_D static
void gkyl_array_integrate_choose_kernel(enum gkyl_array_integrate_op op,
  const struct gkyl_basis *basis, struct gkyl_array_integrate *up)
{
  int ndim = basis->ndim, poly_order = basis->poly_order;

  if (op == GKYL_ARRAY_INTEGRATE_OP_NONE) {
    up->kernel = gkyl_array_integrate_none_ker_list_ser[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_ABS) {
    up->kernel = gkyl_array_integrate_abs_ker_list_ser[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_SQ) {
    up->kernel = gkyl_array_integrate_sq_ker_list_ser[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED) {
    if (basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY)
      up->kernel = gkyl_array_integrate_sq_weighted_ker_list_ser[ndim-1].kernels[poly_order-1];
    else if (basis->b_type == GKYL_BASIS_MODAL_GKHYBRID)
      up->kernel = gkyl_array_integrate_sq_weighted_ker_list_gkhyb[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ) {
    up->kernel = gkyl_array_integrate_gradsq_ker_list[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_GRADPERP_SQ) {
    up->kernel = gkyl_array_integrate_gradperpsq_ker_list[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ) {
    up->kernel = gkyl_array_integrate_epsgradperpsq_ker_list[ndim-1].kernels[poly_order-1];
  }
  else {
    assert(false);
  }
  assert(up->kernel);
}

struct gkyl_array_integrate*
gkyl_array_integrate_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_comp, enum gkyl_array_integrate_op op);

void gkyl_array_integrate_advance_cu(gkyl_array_integrate *up, const struct gkyl_array *arr,
  double factor, const struct gkyl_array *weight, const struct gkyl_range *range, const struct gkyl_range *weight_range, double *out);
