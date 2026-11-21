// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_maxwell_kernels.h>
#include <gkyl_pkpm_vars_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef int (*pkpm_bvar_t)(const double *em, double *GKYL_RESTRICT b);
typedef int (*pkpm_dtB_over_magB_t)(
  const double *em, const double *dt_em, double *GKYL_RESTRICT out
);

typedef void (*pkpm_bvar_grad_t)(
  const double *dxv, const double *jacob_pos, const double *b_l, const double *b_c,
  const double *b_r, double *GKYL_RESTRICT grad_b
);

typedef struct {
  pkpm_bvar_grad_t kernels[3];
} gkyl_dg_pkpm_bvar_grad_kern_list;

struct gkyl_dg_calc_pkpm_bvar {
  int cdim; // Configuration space dimensionality
  pkpm_bvar_t pkpm_bvar; // kernel computing b = B/|B| at the Gauss-Legendre nodes
  pkpm_dtB_over_magB_t pkpm_dtB_over_magB; // kernel computing d_t B/|B| at the same nodes
  pkpm_bvar_grad_t pkpm_grad[3]; // weak (surface-averaged) gradient of b in each direction
  double dx[GKYL_MAX_CDIM]; // configuration-space cell spacing (for the gradient)
  const struct gkyl_vlasov_position_map
    *pos_map; // Position map (acquired; host pointer, 0 on the device copy)
  const struct gkyl_array
    *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device)
  int num_null_b; // number of null-field nodes in the last advance (host copy)
  int *num_null_b_dev; // device counter (0 on host)

  uint32_t flags;
  struct gkyl_dg_calc_pkpm_bvar *on_dev; // pointer to itself or device data
};

// Weak gradient of b (indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_bvar_grad_kern_list ser_pkpm_bvar_grad_kernels[] = {
  {pkpm_bvar_grad_x_1x_ser_p1, NULL, NULL}, // 0
  {pkpm_bvar_grad_x_2x_ser_p1, pkpm_bvar_grad_y_2x_ser_p1, NULL}, // 1
  {pkpm_bvar_grad_x_3x_ser_p1, pkpm_bvar_grad_y_3x_ser_p1, pkpm_bvar_grad_z_3x_ser_p1} // 2
};

GKYL_CU_D static pkpm_bvar_grad_t
choose_pkpm_bvar_grad_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  if (poly_order == 1 &&
      (b_type == GKYL_BASIS_MODAL_SERENDIPITY || b_type == GKYL_BASIS_MODAL_TENSOR)) {
    return ser_pkpm_bvar_grad_kernels[cdim - 1].kernels[dir];
  }
  return 0;
}

// Direct nodal normalization kernels (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_bvar_t ser_pkpm_bvar_kernels[] = {
  em_bvar_nodal_1x_ser_p1, em_bvar_nodal_2x_ser_p1, em_bvar_nodal_3x_ser_p1
};

// d_t B/|B| at the nodes of b (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_dtB_over_magB_t ser_pkpm_dtB_over_magB_kernels[] = {
  pkpm_dtB_over_magB_1x_ser_p1, pkpm_dtB_over_magB_2x_ser_p1, pkpm_dtB_over_magB_3x_ser_p1
};

GKYL_CU_D static pkpm_dtB_over_magB_t
choose_pkpm_dtB_over_magB_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 &&
      (b_type == GKYL_BASIS_MODAL_SERENDIPITY || b_type == GKYL_BASIS_MODAL_TENSOR)) {
    return ser_pkpm_dtB_over_magB_kernels[cdim - 1];
  }
  return 0;
}

GKYL_CU_D static pkpm_bvar_t
choose_pkpm_bvar_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
    case GKYL_BASIS_MODAL_TENSOR: // the same basis at polynomial order 1
      return ser_pkpm_bvar_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}
