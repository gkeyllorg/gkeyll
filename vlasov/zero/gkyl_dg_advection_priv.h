#pragma once

#include <gkyl_advection_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in advection DG equation object creation
// functions

// Types for various kernels
typedef double (*advection_surf_t)(const double *w, const double *dx,
  const double *ul, const double *uc, const double *ur,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_advection_vol_kern_list;
typedef struct { advection_surf_t kernels[3]; } gkyl_dg_advection_surf_kern_list;

struct dg_advection {
  struct gkyl_dg_eqn eqn; // Base object    
  advection_surf_t surf[3]; // pointers to surface kernels
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_advection_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_advection_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection_vol_1x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_advection_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection_vol_1x_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_advection_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection_vol_2x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_advection_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection_vol_2x_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_advection_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection_vol_3x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_advection_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection_vol_3x_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx), qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_advection_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_advection_vol_1x_ser_p1, kernel_advection_vol_1x_ser_p2 }, // 0
  { NULL, kernel_advection_vol_2x_ser_p1, kernel_advection_vol_2x_ser_p2 }, // 1
  { NULL, kernel_advection_vol_3x_ser_p1, kernel_advection_vol_3x_ser_p2 }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, advection_surfx_1x_ser_p1, advection_surfx_1x_ser_p2 }, // 0
  { NULL, advection_surfx_2x_ser_p1, advection_surfx_2x_ser_p2 }, // 1
  { NULL, advection_surfx_3x_ser_p1, advection_surfx_3x_ser_p2 }, // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, advection_surfy_2x_ser_p1, advection_surfy_2x_ser_p2 }, // 1
  { NULL, advection_surfy_3x_ser_p1, advection_surfy_3x_ser_p2 }, // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, advection_surfz_3x_ser_p1, advection_surfz_3x_ser_p2 }, // 2
};

/**
 * Free advection equation object
 *
 * @param ref Reference counter for advection equation
 */
void gkyl_advection_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);

  long cidx_l = gkyl_range_idx(&advection->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&advection->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&advection->conf_range, idxR);

  return advection->surf[dir](xcC, dxC, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx_l),
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx_c),
    (const double*) gkyl_array_cfetch(advection->auxfields.u_i, cidx_r), 
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
