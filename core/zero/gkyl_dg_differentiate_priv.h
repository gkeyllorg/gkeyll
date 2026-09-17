// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_differentiate_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Function pointer type for differentiation kernels.
// Computes the DG derivative of fin in a single direction and at a given
// derivative order, writing the result into out.  Both arrays have
// num_basis entries.  dx is the cell size in the differentiation direction.
typedef void (*differentiate_op_t)(double dx, const double *fin, double *out);

// for use in kernel tables
typedef struct { differentiate_op_t kernels[4]; } differentiate_op_kern_list;

// Serendipity differentiation kernels.
// [ndim-1][dir][diff_order-1].kernels[poly_order-1]
//   ndim  : 1..3
//   dir   : 0..ndim-1
//   order : 1..2
GKYL_CU_D
static const differentiate_op_kern_list ser_differentiate_list[3][3][2] = {
  // ndim = 1
  {
    // dir0
    {
      { dg_differentiate_local_ser_1x_p1_dir0_order1, dg_differentiate_local_ser_1x_p2_dir0_order1, NULL }, // order1
      { dg_differentiate_local_ser_1x_p1_dir0_order2, dg_differentiate_local_ser_1x_p2_dir0_order2, NULL }, // order2
    },
    // dir1: not applicable for 1D
    {
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
    },
    // dir2: not applicable for 1D
    {
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
    },
  },
  // ndim = 2
  {
    // dir0
    {
      { dg_differentiate_local_ser_2x_p1_dir0_order1, dg_differentiate_local_ser_2x_p2_dir0_order1, NULL }, // order1
      { dg_differentiate_local_ser_2x_p1_dir0_order2, dg_differentiate_local_ser_2x_p2_dir0_order2, NULL }, // order2
    },
    // dir1
    {
      { dg_differentiate_local_ser_2x_p1_dir1_order1, dg_differentiate_local_ser_2x_p2_dir1_order1, NULL }, // order1
      { dg_differentiate_local_ser_2x_p1_dir1_order2, dg_differentiate_local_ser_2x_p2_dir1_order2, NULL }, // order2
    },
    // dir2: not applicable for 2D
    {
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
    },
  },
  // ndim = 3
  {
    // dir0
    {
      { dg_differentiate_local_ser_3x_p1_dir0_order1, NULL, NULL }, // order1
      { dg_differentiate_local_ser_3x_p1_dir0_order2, NULL, NULL }, // order2
    },
    // dir1
    {
      { dg_differentiate_local_ser_3x_p1_dir1_order1, NULL, NULL }, // order1
      { dg_differentiate_local_ser_3x_p1_dir1_order2, NULL, NULL }, // order2
    },
    // dir2
    {
      { dg_differentiate_local_ser_3x_p1_dir2_order1, NULL, NULL }, // order1
      { dg_differentiate_local_ser_3x_p1_dir2_order2, NULL, NULL }, // order2
    },
  },
};

// Tensor differentiation kernels.
// p1 tensor == p1 serendipity; tensor-specific kernels exist only for p2.
// No tensor kernels for 3D.
// [ndim-1][dir][diff_order-1].kernels[poly_order]
GKYL_CU_D
static const differentiate_op_kern_list ten_differentiate_list[3][3][2] = {
  // ndim = 1
  {
    // dir0
    {
      { dg_differentiate_local_ser_1x_p1_dir0_order1, dg_differentiate_local_tensor_1x_p2_dir0_order1, NULL }, // order1
      { dg_differentiate_local_ser_1x_p1_dir0_order2, dg_differentiate_local_tensor_1x_p2_dir0_order2, NULL }, // order2
    },
    // dir1: not applicable for 1D
    {
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
    },
    // dir2: not applicable for 1D
    {
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
    },
  },
  // ndim = 2
  {
    // dir0
    {
      { dg_differentiate_local_ser_2x_p1_dir0_order1, dg_differentiate_local_tensor_2x_p2_dir0_order1, NULL }, // order1
      { dg_differentiate_local_ser_2x_p1_dir0_order2, dg_differentiate_local_tensor_2x_p2_dir0_order2, NULL }, // order2
    },
    // dir1
    {
      { dg_differentiate_local_ser_2x_p1_dir1_order1, dg_differentiate_local_tensor_2x_p2_dir1_order1, NULL }, // order1
      { dg_differentiate_local_ser_2x_p1_dir1_order2, dg_differentiate_local_tensor_2x_p2_dir1_order2, NULL }, // order2
    },
    // dir2: not applicable for 2D
    {
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
    },
  },
  // ndim = 3: no tensor kernels available
  {
    { { NULL, NULL, NULL }, { NULL, NULL, NULL } },
    { { NULL, NULL, NULL }, { NULL, NULL, NULL } },
    { { NULL, NULL, NULL }, { NULL, NULL, NULL } },
  },
};

GKYL_CU_D
static differentiate_op_t
choose_ser_differentiate_kern(int ndim, int dir, int poly_order, int diff_order)
{
  assert(ndim >= 1 && ndim <= 3);
  assert(dir >= 0 && dir < ndim);
  assert(diff_order >= 1 && diff_order <= 2);
  return ser_differentiate_list[ndim-1][dir][diff_order-1].kernels[poly_order-1];
}

GKYL_CU_D
static differentiate_op_t
choose_ten_differentiate_kern(int ndim, int dir, int poly_order, int diff_order)
{
  assert(ndim >= 1 && ndim <= 3);
  assert(dir >= 0 && dir < ndim);
  assert(diff_order >= 1 && diff_order <= 2);
  return ten_differentiate_list[ndim-1][dir][diff_order-1].kernels[poly_order-1];
}

/**
 * Host-side wrappers for GPU implementations.
 */
void gkyl_dg_differentiate_op_local_cu(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp);

void gkyl_dg_differentiate_op_local_range_cu(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp,
  const struct gkyl_range *range);
