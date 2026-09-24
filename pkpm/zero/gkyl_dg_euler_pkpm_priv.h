#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_euler_pkpm_kernels.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>

// private header for use in Euler DG equation object (for parallel-kinetic-perpendicular-moment (pkpm) model) creation
// functions

// Types for various kernels
typedef double (*euler_pkpm_surf_t)(const double *w, const double *dxv, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
  const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r,
  const double *p_ij_l, const double *p_ij_c, const double *p_ij_r,
  const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r,
  const double *pkpm_lax_l, const double *pkpm_lax_r, 
  const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
  double* GKYL_RESTRICT out); 

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_euler_pkpm_vol_kern_list;
typedef struct { euler_pkpm_surf_t kernels[3]; } gkyl_dg_euler_pkpm_surf_kern_list;

struct dg_euler_pkpm {
  struct gkyl_dg_eqn eqn; // Base object  
  euler_pkpm_surf_t surf[3]; // pointers to surface kernels
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_euler_pkpm_auxfields auxfields; // Auxiliary fields.
};

//
// Volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_1x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_1x_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_2x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_2x_tensor_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_2x_tensor_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_3x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx),
    qIn, qRhsOut);
}

// PKPM Fluid Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_euler_pkpm_vol_1x_ser_p1, kernel_euler_pkpm_vol_1x_ser_p2 }, // 0
  { NULL, kernel_euler_pkpm_vol_2x_ser_p1, NULL }, // 1
  { NULL, kernel_euler_pkpm_vol_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_vol_kern_list ten_vol_kernels[] = {
  { NULL, kernel_euler_pkpm_vol_1x_ser_p1, kernel_euler_pkpm_vol_1x_ser_p2 }, // 0
  { NULL, kernel_euler_pkpm_vol_2x_ser_p1, kernel_euler_pkpm_vol_2x_tensor_p2 }, // 1
  { NULL, kernel_euler_pkpm_vol_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, euler_pkpm_surfx_1x_ser_p1, euler_pkpm_surfx_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_surfx_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_surfx_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, euler_pkpm_surfx_1x_ser_p1, euler_pkpm_surfx_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_surfx_2x_ser_p1, euler_pkpm_surfx_2x_tensor_p2 }, // 1
  { NULL, euler_pkpm_surfx_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, euler_pkpm_surfy_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_surfy_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, euler_pkpm_surfy_2x_ser_p1, euler_pkpm_surfy_2x_tensor_p2 }, // 1
  { NULL, euler_pkpm_surfy_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, euler_pkpm_surfz_3x_ser_p1, NULL }, // 2
};

// PKPM Fluid Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, euler_pkpm_surfz_3x_ser_p1, NULL }, // 2
};

/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_euler_pkpm_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);

  long cidx_l = gkyl_range_idx(&euler_pkpm->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&euler_pkpm->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&euler_pkpm->conf_range, idxR);

  return euler_pkpm->surf[dir](xcC, dxC, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.vlasov_pkpm_moms, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.vlasov_pkpm_moms, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.vlasov_pkpm_moms, cidx_r),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim_surf, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim_surf, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_prim_surf, cidx_r),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_p_ij, cidx_r),
    qInL, qInC, qInR, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_lax, cidx_c), 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_lax, cidx_r), 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_penalization, cidx_c), 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pkpm_penalization, cidx_r), 
    qRhsOut);
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  return 0.;
}
