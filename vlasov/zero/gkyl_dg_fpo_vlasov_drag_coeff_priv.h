// Private header: not for direct use
#pragma once


#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_fpo_vlasov_kernels.h>

// Kernel function pointers
typedef void (*fpo_drag_coeff_t)(const double *dxv, const double *gamma, 
    const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, 
    double *drag_coeff_surf);

typedef void (*fpo_sgn_drag_coeff_t)(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf);

// For use in kernel tables
typedef struct { fpo_drag_coeff_t kernels[3]; } gkyl_dg_fpo_drag_coeff_kern_list;
typedef struct { gkyl_dg_fpo_drag_coeff_kern_list list[3]; } gkyl_dg_fpo_drag_coeff_stencil_list;

typedef struct { fpo_sgn_drag_coeff_t kernels[3]; } gkyl_dg_fpo_sgn_drag_coeff_kern_list;
typedef struct { gkyl_dg_fpo_sgn_drag_coeff_kern_list list[3]; } gkyl_dg_fpo_sgn_drag_coeff_stencil_list;

// drag coefficient kernel lists
GKYL_CU_D
static const gkyl_dg_fpo_drag_coeff_stencil_list ser_fpo_drag_coeff_vx_kernels[] = {
  {{ // 1x3v
    {NULL, NULL, NULL},
    {fpo_drag_coeff_1x3v_vx_ser_p1_invx, fpo_drag_coeff_1x3v_vx_ser_p1_lovx, fpo_drag_coeff_1x3v_vx_ser_p1_upvx},
    {fpo_drag_coeff_1x3v_vx_ser_p2_invx, fpo_drag_coeff_1x3v_vx_ser_p2_lovx, fpo_drag_coeff_1x3v_vx_ser_p2_upvx},
  }},
  {{ // 2x3v
    {NULL, NULL, NULL},
    {fpo_drag_coeff_2x3v_vx_ser_p1_invx, fpo_drag_coeff_2x3v_vx_ser_p1_lovx, fpo_drag_coeff_2x3v_vx_ser_p1_upvx},
    {NULL, NULL, NULL},
  }}
};

GKYL_CU_D
static const gkyl_dg_fpo_drag_coeff_stencil_list ser_fpo_drag_coeff_vy_kernels[] = {
  {{ // 1x3v
    {NULL, NULL, NULL},
    {fpo_drag_coeff_1x3v_vy_ser_p1_invy, fpo_drag_coeff_1x3v_vy_ser_p1_lovy, fpo_drag_coeff_1x3v_vy_ser_p1_upvy},
    {fpo_drag_coeff_1x3v_vy_ser_p2_invy, fpo_drag_coeff_1x3v_vy_ser_p2_lovy, fpo_drag_coeff_1x3v_vy_ser_p2_upvy},
  }},
  {{ // 2x3v
    {NULL, NULL, NULL},
    {fpo_drag_coeff_2x3v_vy_ser_p1_invy, fpo_drag_coeff_2x3v_vy_ser_p1_lovy, fpo_drag_coeff_2x3v_vy_ser_p1_upvy},
    {NULL, NULL, NULL},
  }}
};

GKYL_CU_D
static const gkyl_dg_fpo_drag_coeff_stencil_list ser_fpo_drag_coeff_vz_kernels[] = {
  {{ // 1x3v
    {NULL, NULL, NULL},
    {fpo_drag_coeff_1x3v_vz_ser_p1_invz, fpo_drag_coeff_1x3v_vz_ser_p1_lovz, fpo_drag_coeff_1x3v_vz_ser_p1_upvz},
    {fpo_drag_coeff_1x3v_vz_ser_p2_invz, fpo_drag_coeff_1x3v_vz_ser_p2_lovz, fpo_drag_coeff_1x3v_vz_ser_p2_upvz},
  }},
  {{ // 2x3v
    {NULL, NULL, NULL},
    {fpo_drag_coeff_2x3v_vz_ser_p1_invz, fpo_drag_coeff_2x3v_vz_ser_p1_lovz, fpo_drag_coeff_2x3v_vz_ser_p1_upvz},
    {NULL, NULL, NULL},
  }}
};

// sgn drag coefficient kernel lists
GKYL_CU_D
static const gkyl_dg_fpo_sgn_drag_coeff_stencil_list ser_fpo_sgn_drag_coeff_vx_kernels[] = {
  {{ // 1x3v
    {NULL, NULL, NULL},
    {fpo_sgn_drag_coeff_1x3v_vx_ser_p1_invx, fpo_sgn_drag_coeff_1x3v_vx_ser_p1_lovx, fpo_sgn_drag_coeff_1x3v_vx_ser_p1_upvx},
    {fpo_sgn_drag_coeff_1x3v_vx_ser_p2_invx, fpo_sgn_drag_coeff_1x3v_vx_ser_p2_lovx, fpo_sgn_drag_coeff_1x3v_vx_ser_p2_upvx},
  }},
  {{ // 2x3v
    {NULL, NULL, NULL},
    {fpo_sgn_drag_coeff_2x3v_vx_ser_p1_invx, fpo_sgn_drag_coeff_2x3v_vx_ser_p1_lovx, fpo_sgn_drag_coeff_2x3v_vx_ser_p1_upvx},
    {NULL, NULL, NULL}
  }}
};

GKYL_CU_D
static const gkyl_dg_fpo_sgn_drag_coeff_stencil_list ser_fpo_sgn_drag_coeff_vy_kernels[] = {
  {{ // 1x3v
    {NULL, NULL, NULL},
    {fpo_sgn_drag_coeff_1x3v_vy_ser_p1_invy, fpo_sgn_drag_coeff_1x3v_vy_ser_p1_lovy, fpo_sgn_drag_coeff_1x3v_vy_ser_p1_upvy},
    {fpo_sgn_drag_coeff_1x3v_vy_ser_p2_invy, fpo_sgn_drag_coeff_1x3v_vy_ser_p2_lovy, fpo_sgn_drag_coeff_1x3v_vy_ser_p2_upvy},
  }},
  {{ // 2x3v
    {NULL, NULL, NULL},
    {fpo_sgn_drag_coeff_2x3v_vy_ser_p1_invy, fpo_sgn_drag_coeff_2x3v_vy_ser_p1_lovy, fpo_sgn_drag_coeff_2x3v_vy_ser_p1_upvy},
    {NULL, NULL, NULL}
  }}
};

GKYL_CU_D
static const gkyl_dg_fpo_sgn_drag_coeff_stencil_list ser_fpo_sgn_drag_coeff_vz_kernels[] = {
  {{ // 1x3v
    {NULL, NULL, NULL},
    {fpo_sgn_drag_coeff_1x3v_vz_ser_p1_invz, fpo_sgn_drag_coeff_1x3v_vz_ser_p1_lovz, fpo_sgn_drag_coeff_1x3v_vz_ser_p1_upvz},
    {fpo_sgn_drag_coeff_1x3v_vz_ser_p2_invz, fpo_sgn_drag_coeff_1x3v_vz_ser_p2_lovz, fpo_sgn_drag_coeff_1x3v_vz_ser_p2_upvz},
  }},
  {{ // 2x3v
    {NULL, NULL, NULL},
    {fpo_sgn_drag_coeff_2x3v_vz_ser_p1_invz, fpo_sgn_drag_coeff_2x3v_vz_ser_p1_lovz, fpo_sgn_drag_coeff_2x3v_vz_ser_p1_upvz},
    {NULL, NULL, NULL}
  }}
};

GKYL_CU_D
static fpo_drag_coeff_t
choose_ser_fpo_drag_coeff_recovery_kern(int dir, int cdim, int poly_order, int stencil_idx)
{
  if (dir == 0)
    return ser_fpo_drag_coeff_vx_kernels[cdim-1].list[poly_order].kernels[stencil_idx];
  else if (dir == 1)
    return ser_fpo_drag_coeff_vy_kernels[cdim-1].list[poly_order].kernels[stencil_idx];
  else if (dir == 2)
    return ser_fpo_drag_coeff_vz_kernels[cdim-1].list[poly_order].kernels[stencil_idx];
  else
    return NULL;
};

GKYL_CU_D
static fpo_sgn_drag_coeff_t
choose_ser_fpo_sgn_drag_coeff_recovery_kern(int dir, int cdim, int poly_order, int stencil_idx)
{
  if (dir == 0)
    return ser_fpo_sgn_drag_coeff_vx_kernels[cdim-1].list[poly_order].kernels[stencil_idx];
  else if (dir == 1)
    return ser_fpo_sgn_drag_coeff_vy_kernels[cdim-1].list[poly_order].kernels[stencil_idx];
  else if (dir == 2)
    return ser_fpo_sgn_drag_coeff_vz_kernels[cdim-1].list[poly_order].kernels[stencil_idx];
  else
    return NULL;
};

