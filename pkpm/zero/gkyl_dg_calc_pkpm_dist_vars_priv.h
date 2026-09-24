// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_pkpm_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*pkpm_dist_mirror_force_t)(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 

typedef void (*pkpm_dist_div_ppar_t)(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 

// for use in kernel tables
typedef struct { pkpm_dist_mirror_force_t kernels[3]; } gkyl_dg_pkpm_dist_mirror_force_kern_list;
typedef struct { pkpm_dist_div_ppar_t kernels[3]; } gkyl_dg_pkpm_dist_div_ppar_kern_list;

struct gkyl_dg_calc_pkpm_dist_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  pkpm_dist_mirror_force_t pkpm_dist_mirror_force; // kernel for computing distribution function sources
  pkpm_dist_div_ppar_t pkpm_dist_div_ppar[3];  // kernel for computing consistent div(p_par b)

  uint32_t flags;
  struct gkyl_dg_calc_pkpm_dist_vars *on_dev; // pointer to itself or device data
};

// PKPM distribution function source in mirror force and vperp characteristics (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_mirror_force_kern_list ser_pkpm_dist_mirror_force_kernels[] = {
  { NULL, pkpm_dist_mirror_force_1x1v_ser_p1, pkpm_dist_mirror_force_1x1v_ser_p2 }, // 0
  { NULL, pkpm_dist_mirror_force_2x1v_ser_p1, NULL }, // 1
  { NULL, pkpm_dist_mirror_force_3x1v_ser_p1, NULL }, // 2
};

// PKPM distribution function source in mirror force and vperp characteristics (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_mirror_force_kern_list ten_pkpm_dist_mirror_force_kernels[] = {
  { NULL, pkpm_dist_mirror_force_1x1v_ser_p1, pkpm_dist_mirror_force_1x1v_tensor_p2 }, // 0
  { NULL, pkpm_dist_mirror_force_2x1v_ser_p1, pkpm_dist_mirror_force_2x1v_tensor_p2 }, // 1
  { NULL, pkpm_dist_mirror_force_3x1v_ser_p1, NULL }, // 2
};

// PKPM consistent div(p_par b) (in x) kernels (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_div_ppar_kern_list ser_pkpm_dist_div_ppar_x_kernels[] = {
  { NULL, pkpm_dist_div_ppar_x_1x1v_ser_p1, pkpm_dist_div_ppar_x_1x1v_ser_p2 }, // 0
  { NULL, pkpm_dist_div_ppar_x_2x1v_ser_p1, NULL }, // 1
  { NULL, pkpm_dist_div_ppar_x_3x1v_ser_p1, NULL }, // 2
};

// PKPM consistent div(p_par b) (in y) kernels (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_div_ppar_kern_list ser_pkpm_dist_div_ppar_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, pkpm_dist_div_ppar_y_2x1v_ser_p1, NULL }, // 1
  { NULL, pkpm_dist_div_ppar_y_3x1v_ser_p1, NULL }, // 2
};

// PKPM consistent div(p_par b) (in z) kernels (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_div_ppar_kern_list ser_pkpm_dist_div_ppar_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, pkpm_dist_div_ppar_z_3x1v_ser_p1, NULL }, // 2
};

// PKPM consistent div(p_par b) (in x) kernels (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_div_ppar_kern_list ten_pkpm_dist_div_ppar_x_kernels[] = {
  { NULL, pkpm_dist_div_ppar_x_1x1v_ser_p1, pkpm_dist_div_ppar_x_1x1v_tensor_p2 }, // 0
  { NULL, pkpm_dist_div_ppar_x_2x1v_ser_p1, pkpm_dist_div_ppar_x_2x1v_tensor_p2 }, // 1
  { NULL, pkpm_dist_div_ppar_x_3x1v_ser_p1, NULL }, // 2
};

// PKPM consistent div(p_par b) (in y) kernels (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_div_ppar_kern_list ten_pkpm_dist_div_ppar_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, pkpm_dist_div_ppar_y_2x1v_ser_p1, pkpm_dist_div_ppar_y_2x1v_tensor_p2 }, // 1
  { NULL, pkpm_dist_div_ppar_y_3x1v_ser_p1, NULL }, // 2
};

// PKPM consistent div(p_par b) (in z) kernels (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_div_ppar_kern_list ten_pkpm_dist_div_ppar_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, pkpm_dist_div_ppar_z_3x1v_ser_p1, NULL }, // 2
};

GKYL_CU_D
static pkpm_dist_mirror_force_t
choose_pkpm_dist_mirror_force_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_dist_mirror_force_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_dist_mirror_force_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_dist_div_ppar_t
choose_pkpm_dist_div_ppar_kern(int dir, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_pkpm_dist_div_ppar_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_pkpm_dist_div_ppar_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ser_pkpm_dist_div_ppar_z_kernels[cdim-1].kernels[poly_order];
      else 
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return ten_pkpm_dist_div_ppar_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ten_pkpm_dist_div_ppar_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ten_pkpm_dist_div_ppar_z_kernels[cdim-1].kernels[poly_order];
      else 
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}
