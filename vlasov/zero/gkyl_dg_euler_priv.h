#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>

// private header for use in euler DG equation object creation
// functions

// Types for various kernels
typedef double (*euler_surf_t)(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
  const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
  const double *u_surf_l, const double *u_surf_c, const double *u_surf_r,
  const double *p_surf_l, const double *p_surf_c, const double *p_surf_r, 
  const double *fluid_l, const double *fluid_c, const double *fluid_r, 
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[4]; } gkyl_dg_euler_vol_kern_list;
typedef struct { euler_surf_t kernels[4]; } gkyl_dg_euler_surf_kern_list;

struct dg_euler {
  struct gkyl_dg_eqn eqn; // Base object  
  euler_surf_t surf[3]; // pointers to surface kernels
  enum gkyl_eqn_type eqn_type; // Equation type
  const struct gkyl_wv_eqn *wv_eqn; // wave equation object for Roe solve
  const struct gkyl_wave_geom *geom; // wave geometry
  double gas_gamma; // adiabatic index
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_euler_auxfields auxfields; // Auxiliary fields.
};

GKYL_CU_DH
static double
kernel_euler_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_1x_ser_p1(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_1x_ser_p2(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_1x_ser_p3(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_1x_ser_p3(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_2x_ser_p1(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_2x_tensor_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_2x_tensor_p2(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_3x_ser_p1(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p, cidx),
    qIn, qRhsOut);
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_euler_vol_1x_ser_p1, kernel_euler_vol_1x_ser_p2, kernel_euler_vol_1x_ser_p3 }, // 0
  { NULL, kernel_euler_vol_2x_ser_p1, NULL, NULL }, // 1
  { NULL, kernel_euler_vol_3x_ser_p1, NULL, NULL }, // 2
};

// Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, euler_surfx_1x_ser_p1, euler_surfx_1x_ser_p2, euler_surfx_1x_ser_p3 }, // 0
  { NULL, euler_surfx_2x_ser_p1, NULL, NULL }, // 1
  { NULL, euler_surfx_3x_ser_p1, NULL, NULL }, // 2
};

// Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, euler_surfy_2x_ser_p1, NULL, NULL }, // 1
  { NULL, euler_surfy_3x_ser_p1, NULL, NULL }, // 2
};

// Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, euler_surfz_3x_ser_p1, NULL, NULL }, // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_vol_kern_list ten_vol_kernels[] = {
  { NULL, kernel_euler_vol_1x_ser_p1, kernel_euler_vol_1x_ser_p2, kernel_euler_vol_1x_ser_p3 }, // 0
  { NULL, kernel_euler_vol_2x_ser_p1, kernel_euler_vol_2x_tensor_p2, NULL }, // 1
  { NULL, kernel_euler_vol_3x_ser_p1, NULL, NULL }, // 2
};

// Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, euler_surfx_1x_ser_p1, euler_surfx_1x_ser_p2, euler_surfx_1x_ser_p3 }, // 0
  { NULL, euler_surfx_2x_ser_p1, euler_surfx_2x_tensor_p2, NULL }, // 1
  { NULL, euler_surfx_3x_ser_p1, NULL, NULL }, // 2
};

// Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, euler_surfy_2x_ser_p1, euler_surfy_2x_tensor_p2, NULL }, // 1
  { NULL, euler_surfy_3x_ser_p1, NULL, NULL }, // 2
};

// Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, euler_surfz_3x_ser_p1, NULL, NULL }, // 2
};

/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_dg_euler_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);

  long cidx_l = gkyl_range_idx(&euler->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&euler->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&euler->conf_range, idxR);

  const struct gkyl_wave_cell_geom *geom_l = gkyl_wave_geom_get(euler->geom, idxC);
  const struct gkyl_wave_cell_geom *geom_r = gkyl_wave_geom_get(euler->geom, idxR);

  return euler->surf[dir](xcC, dxC, 
    euler->wv_eqn, geom_l, geom_r, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_surf, cidx_l),
    (const double*) gkyl_array_cfetch(euler->auxfields.u_surf, cidx_c),
    (const double*) gkyl_array_cfetch(euler->auxfields.u_surf, cidx_r), 
    (const double*) gkyl_array_cfetch(euler->auxfields.p_surf, cidx_l),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_surf, cidx_c),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_surf, cidx_r), 
    qInL, qInC, qInR, qRhsOut);
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
