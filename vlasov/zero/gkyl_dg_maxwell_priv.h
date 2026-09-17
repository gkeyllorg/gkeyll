#pragma once

#include <gkyl_maxwell_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in Maxwell DG equation object creation
// functions

// Types for various kernels
typedef double (*maxwell_surf_t)(const gkyl_maxwell_inp *meq, const double *w, const double *dx,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[4]; } gkyl_dg_maxwell_vol_kern_list;
typedef struct { maxwell_surf_t kernels[4]; } gkyl_dg_maxwell_surf_kern_list;

struct dg_maxwell {
  struct gkyl_dg_eqn eqn; // Base object    
  gkyl_maxwell_inp maxwell_data; // Parameters needed by kernels
  maxwell_surf_t surf[3]; // pointers to surface kernels
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_maxwell_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_1x_ser_p1(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_1x_ser_p2(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_1x_ser_p3(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_1x_ser_p3(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_2x_ser_p1(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_2x_ser_p2(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_ser_p3(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_2x_ser_p3(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_tensor_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_2x_tensor_p2(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_3x_ser_p1(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell_vol_3x_ser_p2(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_maxwell_vol_1x_ser_p1, kernel_maxwell_vol_1x_ser_p2, kernel_maxwell_vol_1x_ser_p3 }, // 0
  { NULL, kernel_maxwell_vol_2x_ser_p1, kernel_maxwell_vol_2x_ser_p2, kernel_maxwell_vol_2x_ser_p3 }, // 1
  { NULL, kernel_maxwell_vol_3x_ser_p1, kernel_maxwell_vol_3x_ser_p2, NULL },              // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ten_vol_kernels[] = {
  { NULL, kernel_maxwell_vol_1x_ser_p1, kernel_maxwell_vol_1x_ser_p2, kernel_maxwell_vol_1x_ser_p3 }, // 0
  { NULL, kernel_maxwell_vol_2x_ser_p1, kernel_maxwell_vol_2x_tensor_p2, NULL }, // 1
  { NULL, kernel_maxwell_vol_3x_ser_p1, NULL, NULL },              // 2
};

// Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2, maxwell_surfx_1x_ser_p3 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_ser_p2, maxwell_surfx_2x_ser_p3 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, maxwell_surfx_3x_ser_p2, NULL },                 // 2
};

// Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2, maxwell_surfx_1x_ser_p3 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_tensor_p2, NULL }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL, NULL },                 // 2
};

// Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_ser_p2, maxwell_surfy_2x_ser_p3 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, maxwell_surfy_3x_ser_p2, NULL },                 // 2
};

// Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_tensor_p2, NULL }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL, NULL },                 // 2
};

// Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, maxwell_surfz_3x_ser_p2, NULL }, // 2
};

// Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, NULL, NULL }, // 2
};

/**
 * Free Maxwell equation object
 *
 * @param ref Reference counter for Maxwell equation
 */
void gkyl_maxwell_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell->surf[dir](&maxwell->maxwell_data, xcC, dxC,
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
