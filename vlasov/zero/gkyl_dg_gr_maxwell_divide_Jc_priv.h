// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*divide_Jc_t)(
  const double *J, const double *Jfields, double *GKYL_RESTRICT fields_no_J
);
typedef void (*rescale_Jc_t)(
  const double *J, const double *fields_no_J, double *GKYL_RESTRICT Jfields
);

// for use in kernel tables
typedef struct {
  divide_Jc_t kernels[4];
} gkyl_dg_divide_Jc_kern_list;
typedef struct {
  rescale_Jc_t kernels[4];
} gkyl_dg_rescale_Jc_kern_list;

// Configuration-space Jacobian division kernels (Serendipity basis).
GKYL_CU_D static const gkyl_dg_divide_Jc_kern_list ser_divide_Jc_kernels[] = {
  // 1x kernels
  {NULL, dg_gr_maxwell_divide_Jc_1x_ser_p1, dg_gr_maxwell_divide_Jc_1x_ser_p2, NULL}, // 0
  // 2x kernels
  {NULL, dg_gr_maxwell_divide_Jc_2x_ser_p1, dg_gr_maxwell_divide_Jc_2x_ser_p2, NULL}, // 1
  // 3x kernels
  {NULL, dg_gr_maxwell_divide_Jc_3x_ser_p1, dg_gr_maxwell_divide_Jc_3x_ser_p2, NULL} // 2
};

// Configuration-space Jacobian division kernels (Tensor basis).
GKYL_CU_D static const gkyl_dg_divide_Jc_kern_list tensor_divide_Jc_kernels[] = {
  // 1x kernels
  {NULL, dg_gr_maxwell_divide_Jc_1x_tensor_p1, dg_gr_maxwell_divide_Jc_1x_tensor_p2, NULL}, // 0
  // 2x kernels
  {NULL, dg_gr_maxwell_divide_Jc_2x_tensor_p1, dg_gr_maxwell_divide_Jc_2x_tensor_p2, NULL}, // 1
  // 3x kernels
  {NULL, dg_gr_maxwell_divide_Jc_3x_tensor_p1, dg_gr_maxwell_divide_Jc_3x_tensor_p2, NULL} // 2
};

// Configuration-space Jacobian rescale kernels (Serendipity basis).
GKYL_CU_D static const gkyl_dg_rescale_Jc_kern_list ser_rescale_Jc_kernels[] = {
  // 1x kernels
  {NULL, dg_gr_maxwell_rescale_Jc_1x_ser_p1, dg_gr_maxwell_rescale_Jc_1x_ser_p2, NULL}, // 0
  // 2x kernels
  {NULL, dg_gr_maxwell_rescale_Jc_2x_ser_p1, dg_gr_maxwell_rescale_Jc_2x_ser_p2, NULL}, // 1
  // 3x kernels
  {NULL, dg_gr_maxwell_rescale_Jc_3x_ser_p1, dg_gr_maxwell_rescale_Jc_3x_ser_p2, NULL} // 2
};

// Configuration-space Jacobian rescale kernels (Tensor basis).
GKYL_CU_D static const gkyl_dg_rescale_Jc_kern_list tensor_rescale_Jc_kernels[] = {
  // 1x kernels
  {NULL, dg_gr_maxwell_rescale_Jc_1x_tensor_p1, dg_gr_maxwell_rescale_Jc_1x_tensor_p2, NULL}, // 0
  // 2x kernels
  {NULL, dg_gr_maxwell_rescale_Jc_2x_tensor_p1, dg_gr_maxwell_rescale_Jc_2x_tensor_p2, NULL}, // 1
  // 3x kernels
  {NULL, dg_gr_maxwell_rescale_Jc_3x_tensor_p1, dg_gr_maxwell_rescale_Jc_3x_tensor_p2, NULL} // 2
};

GKYL_CU_D static divide_Jc_t
choose_ser_divide_Jc_kern(int cdim, int poly_order)
{
  return ser_divide_Jc_kernels[cdim - 1].kernels[poly_order];
}

GKYL_CU_D static divide_Jc_t
choose_tensor_divide_Jc_kern(int cdim, int poly_order)
{
  return tensor_divide_Jc_kernels[cdim - 1].kernels[poly_order];
}

GKYL_CU_D static rescale_Jc_t
choose_ser_rescale_Jc_kern(int cdim, int poly_order)
{
  return ser_rescale_Jc_kernels[cdim - 1].kernels[poly_order];
}

GKYL_CU_D static rescale_Jc_t
choose_tensor_rescale_Jc_kern(int cdim, int poly_order)
{
  return tensor_rescale_Jc_kernels[cdim - 1].kernels[poly_order];
}