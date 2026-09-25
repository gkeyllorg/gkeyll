// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_canonical_pb_vars_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>
#include <gkyl_vlasov_hyb_build.h>

typedef void (*canonical_pb_m1i_contra_to_cov_t)(
  const double *h_ij, const double *v_i, const double *nv_i, double *GKYL_RESTRICT v_i_cov,
  double *GKYL_RESTRICT nv_i_cov
);
typedef void (*canonical_pb_pressure_t)(
  const double *h_ij_inv, const double *MEnergy, const double *v_i, const double *nv_i,
  double *GKYL_RESTRICT d_Jv_P
);

// for use in kernel tables
typedef struct {
  canonical_pb_m1i_contra_to_cov_t kernels[3];
} gkyl_dg_canonical_pb_m1i_contra_to_cov_kern_list;
typedef struct {
  canonical_pb_pressure_t kernels[3];
} gkyl_dg_canonical_pb_pressure_kern_list;

struct gkyl_dg_calc_canonical_pb_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  canonical_pb_m1i_contra_to_cov_t
    canonical_pb_covariant_u_i; // Canonical pb covariant u_i components
  canonical_pb_pressure_t canonical_pb_pressure; // Canonical pb pressure
  uint32_t flags;
  struct gkyl_dg_calc_canonical_pb_vars *on_dev; // pointer to itself or device data
};

// The cv_index[cd].vdim[cd] is used to index the various list of
// kernels below
static struct {
  int vdim[4];
} cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1, 0, 1, 2}, // 1x kernel indices
  {-1, -1, 3, 4}, // 2x kernel indices
  {-1, -1, -1, 5} // 3x kernel indices
};

//
// Serendipity surface kernels general geometry
//
//
// Serendipity surface kernels general geometry
//

// canonical_pb contravaraint to covariant conversion (Serendipity kernels)
// (Jnu_i = h_{ij}Jnu^i and u_i = h_{ij}u^j)
GKYL_CU_D static const gkyl_dg_canonical_pb_m1i_contra_to_cov_kern_list
  ser_canonical_pb_m1i_contra_to_cov_kernels[] = {
    // 1x kernels
    {NULL, canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p1,
     canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p2}, // 0
    {NULL, canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p1,
     canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p2}, // 1
    {NULL, canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p1,
     canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p2}, // 2
    // 2x kernels
    {NULL, canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p1,
     canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p2}, // 3
    {NULL, canonical_pb_vars_m1i_contra_to_cov_2x3v_ser_p1, NULL}, // 4
    // 3x kernels
    {NULL, NULL, NULL} // 5
};

// canonical_pb Pressure (d*P*Jv = 2*E - n*h^{ij}*u_i*u_j) (Serendipity kernels)
GKYL_CU_D static const gkyl_dg_canonical_pb_pressure_kern_list ser_canonical_pb_pressure_kernels[] =
  {
    // 1x kernels
    {NULL, canonical_pb_vars_pressure_1x1v_ser_p1, canonical_pb_vars_pressure_1x1v_ser_p2}, // 0
    {NULL, canonical_pb_vars_pressure_1x2v_ser_p1, canonical_pb_vars_pressure_1x2v_ser_p2}, // 1
    {NULL, canonical_pb_vars_pressure_1x3v_ser_p1, canonical_pb_vars_pressure_1x3v_ser_p2}, // 2
    // 2x kernels
    {NULL, canonical_pb_vars_pressure_2x2v_ser_p1, canonical_pb_vars_pressure_2x2v_ser_p2}, // 3
    {NULL, canonical_pb_vars_pressure_2x3v_ser_p1, NULL}, // 4
    // 3x kernels
    {NULL, NULL, NULL} // 5
};

//
// Tensor surface kernels general geometry
//
//
// Tensor surface kernels general geometry
//

// canonical_pb contravaraint to covariant conversion (Tensor kernels)
// (Jnu_i = h_{ij}Jnu^i and u_i = h_{ij}u^j)
GKYL_CU_D static const gkyl_dg_canonical_pb_m1i_contra_to_cov_kern_list
  tensor_canonical_pb_m1i_contra_to_cov_kernels[] = {
    // 1x kernels
    {NULL, canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p1,
     canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p2}, // 0
    {NULL, canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p1,
     canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p2}, // 1
    {NULL, canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p1,
     canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p2}, // 2
    // 2x kernels
    {NULL, canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p1,
     canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p2}, // 3
    {NULL, GKYL_HYB_2X3V(canonical_pb_vars_m1i_contra_to_cov_2x3v_tensor_p1), NULL}, //4
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(canonical_pb_vars_m1i_contra_to_cov_3x3v_tensor_p1), NULL} // 5
};

// canonical_pb Pressure (d*P*Jv = h^{ij}*M2_{ij} - n*h^{ij}*u_i*u_j) (Tensor kernels)
GKYL_CU_D static const gkyl_dg_canonical_pb_pressure_kern_list
  tensor_canonical_pb_pressure_kernels[] = {
    // 1x kernels
    {NULL, canonical_pb_vars_pressure_1x1v_tensor_p1, canonical_pb_vars_pressure_1x1v_tensor_p2
    }, // 0
    {NULL, canonical_pb_vars_pressure_1x2v_tensor_p1, canonical_pb_vars_pressure_1x2v_tensor_p2
    }, // 1
    {NULL, canonical_pb_vars_pressure_1x3v_tensor_p1, canonical_pb_vars_pressure_1x3v_tensor_p2
    }, // 2
    // 2x kernels
    {NULL, canonical_pb_vars_pressure_2x2v_tensor_p1, canonical_pb_vars_pressure_2x2v_tensor_p2
    }, // 3
    {NULL, GKYL_HYB_2X3V(canonical_pb_vars_pressure_2x3v_tensor_p1), NULL}, //4
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(canonical_pb_vars_pressure_3x3v_tensor_p1), NULL} // 5
};

GKYL_CU_D static canonical_pb_m1i_contra_to_cov_t
choose_canonical_pb_m1i_contra_to_cov_kern(
  enum gkyl_basis_type b_type, int cv_index, int cdim, int poly_order
)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      return ser_canonical_pb_m1i_contra_to_cov_kernels[cv_index].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      // Tensor p=1 hybrid: conf-space content is tensor p=1, which shares the
      // tensor tables (the ser tables lack the 3x3v row).
      assert(poly_order == 1);
      return tensor_canonical_pb_m1i_contra_to_cov_kernels[cv_index].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_m1i_contra_to_cov_kernels[cv_index].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D static canonical_pb_pressure_t
choose_canonical_pb_pressure_kern(
  enum gkyl_basis_type b_type, int cv_index, int cdim, int poly_order
)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      return ser_canonical_pb_pressure_kernels[cv_index].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      // Tensor p=1 hybrid: conf-space content is tensor p=1, which shares the
      // tensor tables (the ser tables lack the 3x3v row).
      assert(poly_order == 1);
      return tensor_canonical_pb_pressure_kernels[cv_index].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_pressure_kernels[cv_index].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}
