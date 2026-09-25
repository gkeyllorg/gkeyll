// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_sr_Gamma_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>
#include <gkyl_vlasov_hyb_build.h>

typedef void (*sr_n_set_t)(
  int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i
);

typedef void (*sr_n_copy_t)(
  int count, struct gkyl_nmat *x, const double *M0, double *GKYL_RESTRICT n
);

typedef void (*sr_GammaV_t)(
  const double *u_i, double *GKYL_RESTRICT u_i_sq, double *GKYL_RESTRICT GammaV,
  double *GKYL_RESTRICT GammaV_sq
);

typedef void (*sr_pressure_t)(
  const double *w, const double *dxv, const double *vmap, const double *gamma,
  const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV,
  const double *GammaV_sq, const double *f, double *GKYL_RESTRICT sr_pressure
);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D static struct {
  int vdim[4];
} cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1, 0, 1, 2}, // 1x kernel indices
  {-1, 3, 4, 5}, // 2x kernel indices
  {-1, -1, -1, 6} // 3x kernel indices
};

// for use in kernel tables
typedef struct {
  sr_n_set_t kernels[4];
} gkyl_dg_sr_vars_n_set_kern_list;
typedef struct {
  sr_n_copy_t kernels[4];
} gkyl_dg_sr_vars_n_copy_kern_list;
typedef struct {
  sr_GammaV_t kernels[4];
} gkyl_dg_sr_vars_GammaV_kern_list;
typedef struct {
  sr_pressure_t kernels[4];
} gkyl_dg_sr_vars_pressure_kern_list;

struct gkyl_dg_calc_sr_vars {
  struct gkyl_rect_grid phase_grid; // Phase-space grid for cell spacing and cell center
    // in pressure velocity moment computation.
  struct gkyl_rect_grid
    vel_grid; // Momentum (four-velocity)-space grid for cell spacing and cell center
  // in gamma and 1/gamma computation.
  struct gkyl_range vel_range; // Momentum (four-velocity)-space range.

  bool use_vmap; // bool to determine if we are using a mapped Momentum (four-velocity)-space grid
  const struct gkyl_array *vmap; // mapping for Momentum (four-velocity)-space

  int poly_order; // polynomial order (determines whether we solve linear system or use basis_inv method).
  struct gkyl_range mem_range; // Configuration-space range for linear solve.

  struct gkyl_nmat *As,
    *xs; // matrices for LHS and RHS for V_drift solve to find rest-frame density.
  gkyl_nmat_mem
    *mem; // memory for use in batched linear solve for V_drift solve to find rest-frame density.
  int Ncomp; // number of components in the linear solve (vdim components of V_drift from weak division).

  sr_n_set_t sr_n_set; // kernel for setting matrices for linear solve for rest-frame density.
  sr_n_copy_t
    sr_n_copy; // kernel for copying solution for V_drift from weak division and computing rest-frame density.
  sr_GammaV_t sr_GammaV; // kernel for computing bulk four-velocity derived quantities such as GammaV.
  sr_pressure_t
    sr_pressure; // kernel for computing rest-frame pressure as a velocity moment of distribution function.

  uint32_t flags;
  struct gkyl_dg_calc_sr_vars *on_dev; // pointer to itself or device data.
};

// Set matrices for computing rest-frame density kernel list (Serendipity kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_n_set_kern_list ser_sr_vars_n_set_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_n_set_1x1v_ser_p1, sr_vars_n_set_1x1v_ser_p2, sr_vars_n_set_1x1v_ser_p3}, // 0
  {NULL, sr_vars_n_set_1x2v_ser_p1, sr_vars_n_set_1x2v_ser_p2, NULL}, // 1
  {NULL, sr_vars_n_set_1x3v_ser_p1, sr_vars_n_set_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_n_set_2x1v_ser_p1, sr_vars_n_set_2x1v_ser_p2, sr_vars_n_set_2x1v_ser_p3}, // 3
  {NULL, sr_vars_n_set_2x2v_ser_p1, sr_vars_n_set_2x2v_ser_p2, NULL}, // 4
  {NULL, sr_vars_n_set_2x3v_ser_p1, sr_vars_n_set_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_n_set_3x3v_ser_p1, NULL, NULL} // 6
};

// Set matrices for computing rest-frame density kernel list (Tensor kernels).
// These are configuration-space kernels: at p=1 the tensor conf basis is
// identical to Serendipity p=1, so the p=1 slots alias the ser kernels.
GKYL_CU_D static const gkyl_dg_sr_vars_n_set_kern_list tensor_sr_vars_n_set_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_n_set_1x1v_ser_p1, sr_vars_n_set_1x1v_tensor_p2, sr_vars_n_set_1x1v_tensor_p3
  }, // 0
  {NULL, sr_vars_n_set_1x2v_ser_p1, sr_vars_n_set_1x2v_tensor_p2, NULL}, // 1
  {NULL, sr_vars_n_set_1x3v_ser_p1, sr_vars_n_set_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_n_set_2x1v_ser_p1, sr_vars_n_set_2x1v_tensor_p2, sr_vars_n_set_2x1v_tensor_p3
  }, // 3
  {NULL, sr_vars_n_set_2x2v_ser_p1, sr_vars_n_set_2x2v_tensor_p2, NULL}, // 4
  {NULL, sr_vars_n_set_2x3v_ser_p1, sr_vars_n_set_2x3v_tensor_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_n_set_3x3v_ser_p1, NULL, NULL} // 6
};

// Copy solution for computing rest-frame density kernel list (Serendipity kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_n_copy_kern_list ser_sr_vars_n_copy_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_n_copy_1x1v_ser_p1, sr_vars_n_copy_1x1v_ser_p2, sr_vars_n_copy_1x1v_ser_p3}, // 0
  {NULL, sr_vars_n_copy_1x2v_ser_p1, sr_vars_n_copy_1x2v_ser_p2, NULL}, // 1
  {NULL, sr_vars_n_copy_1x3v_ser_p1, sr_vars_n_copy_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_n_copy_2x1v_ser_p1, sr_vars_n_copy_2x1v_ser_p2, sr_vars_n_copy_2x1v_ser_p3}, // 3
  {NULL, sr_vars_n_copy_2x2v_ser_p1, sr_vars_n_copy_2x2v_ser_p2, NULL}, // 4
  {NULL, sr_vars_n_copy_2x3v_ser_p1, sr_vars_n_copy_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_n_copy_3x3v_ser_p1, NULL, NULL} // 6
};

// Copy solution for computing rest-frame density kernel list (Tensor kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_n_copy_kern_list tensor_sr_vars_n_copy_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_n_copy_1x1v_ser_p1, sr_vars_n_copy_1x1v_tensor_p2, sr_vars_n_copy_1x1v_tensor_p3
  }, // 0
  {NULL, sr_vars_n_copy_1x2v_ser_p1, sr_vars_n_copy_1x2v_tensor_p2, NULL}, // 1
  {NULL, sr_vars_n_copy_1x3v_ser_p1, sr_vars_n_copy_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_n_copy_2x1v_ser_p1, sr_vars_n_copy_2x1v_tensor_p2, sr_vars_n_copy_2x1v_tensor_p3
  }, // 3
  {NULL, sr_vars_n_copy_2x2v_ser_p1, sr_vars_n_copy_2x2v_tensor_p2, NULL}, // 4
  {NULL, sr_vars_n_copy_2x3v_ser_p1, sr_vars_n_copy_2x3v_tensor_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_n_copy_3x3v_ser_p1, NULL, NULL} // 6
};

// Compute bulk four-velocity derived quantities kernel list (Serendipity kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_GammaV_kern_list ser_sr_vars_GammaV_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_GammaV_1x1v_ser_p1, sr_vars_GammaV_1x1v_ser_p2, sr_vars_GammaV_1x1v_ser_p3}, // 0
  {NULL, sr_vars_GammaV_1x2v_ser_p1, sr_vars_GammaV_1x2v_ser_p2, NULL}, // 1
  {NULL, sr_vars_GammaV_1x3v_ser_p1, sr_vars_GammaV_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_GammaV_2x1v_ser_p1, sr_vars_GammaV_2x1v_ser_p2, sr_vars_GammaV_2x1v_ser_p3}, // 3
  {NULL, sr_vars_GammaV_2x2v_ser_p1, sr_vars_GammaV_2x2v_ser_p2, NULL}, // 4
  {NULL, sr_vars_GammaV_2x3v_ser_p1, sr_vars_GammaV_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_GammaV_3x3v_ser_p1, NULL, NULL} // 6
};

// Compute bulk four-velocity derived quantities kernel list (Tensor kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_GammaV_kern_list tensor_sr_vars_GammaV_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_GammaV_1x1v_ser_p1, sr_vars_GammaV_1x1v_tensor_p2, sr_vars_GammaV_1x1v_tensor_p3
  }, // 0
  {NULL, sr_vars_GammaV_1x2v_ser_p1, sr_vars_GammaV_1x2v_tensor_p2, NULL}, // 1
  {NULL, sr_vars_GammaV_1x3v_ser_p1, sr_vars_GammaV_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_GammaV_2x1v_ser_p1, sr_vars_GammaV_2x1v_tensor_p2, sr_vars_GammaV_2x1v_tensor_p3
  }, // 3
  {NULL, sr_vars_GammaV_2x2v_ser_p1, sr_vars_GammaV_2x2v_tensor_p2, NULL}, // 4
  {NULL, sr_vars_GammaV_2x3v_ser_p1, sr_vars_GammaV_2x3v_tensor_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_GammaV_3x3v_ser_p1, NULL, NULL} // 6
};

// Compute rest-frame pressure kernel list (Serendipity kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_pressure_kern_list ser_sr_vars_pressure_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_pressure_1x1v_ser_p1, sr_vars_pressure_1x1v_ser_p2, sr_vars_pressure_1x1v_ser_p3
  }, // 0
  {NULL, sr_vars_pressure_1x2v_ser_p1, sr_vars_pressure_1x2v_ser_p2, NULL}, // 1
  {NULL, sr_vars_pressure_1x3v_ser_p1, sr_vars_pressure_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_pressure_2x1v_ser_p1, sr_vars_pressure_2x1v_ser_p2, sr_vars_pressure_2x1v_ser_p3
  }, // 3
  {NULL, sr_vars_pressure_2x2v_ser_p1, sr_vars_pressure_2x2v_ser_p2, NULL}, // 4
  {NULL, sr_vars_pressure_2x3v_ser_p1, sr_vars_pressure_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, sr_vars_pressure_3x3v_ser_p1, NULL, NULL} // 6
};

// Compute rest-frame pressure kernel list (Tensor kernels).
GKYL_CU_D static const gkyl_dg_sr_vars_pressure_kern_list tensor_sr_vars_pressure_kernels[] = {
  // 1x kernels
  {NULL, sr_vars_pressure_1x1v_tensor_p1, sr_vars_pressure_1x1v_tensor_p2,
   sr_vars_pressure_1x1v_tensor_p3}, // 0
  {NULL, sr_vars_pressure_1x2v_tensor_p1, sr_vars_pressure_1x2v_tensor_p2, NULL}, // 1
  {NULL, sr_vars_pressure_1x3v_tensor_p1, sr_vars_pressure_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, sr_vars_pressure_2x1v_tensor_p1, sr_vars_pressure_2x1v_tensor_p2,
   sr_vars_pressure_2x1v_tensor_p3}, // 3
  {NULL, sr_vars_pressure_2x2v_tensor_p1, sr_vars_pressure_2x2v_tensor_p2, NULL}, // 4
  {NULL, GKYL_HYB_2X3V(sr_vars_pressure_2x3v_tensor_p1), sr_vars_pressure_2x3v_tensor_p2, NULL
  }, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(sr_vars_pressure_3x3v_tensor_p1), NULL, NULL} // 6
};

GKYL_CU_D static sr_n_set_t
choose_sr_vars_n_set_kern(enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sr_vars_n_set_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sr_vars_n_set_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D static sr_n_copy_t
choose_sr_vars_n_copy_kern(enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sr_vars_n_copy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sr_vars_n_copy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D static sr_GammaV_t
choose_sr_vars_GammaV_kern(enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sr_vars_GammaV_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sr_vars_GammaV_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D static sr_pressure_t
choose_sr_vars_pressure_kern(enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sr_vars_pressure_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sr_vars_pressure_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}
