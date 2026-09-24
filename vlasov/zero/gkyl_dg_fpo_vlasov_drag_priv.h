#pragma once

// Private header, not for direct use in user code

#include <gkyl_fpo_vlasov_kernels.h>

// Types for various kernels
typedef double (*fpo_vlasov_drag_surf_t)(const double *w, const double *dxv,
  const double *h, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*fpo_vlasov_drag_boundary_surf_t)(const double *w, const double *dxv,
  const double *h, 
  const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_vol_kern_list;
typedef struct { fpo_vlasov_drag_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_surf_kern_list;
typedef struct { fpo_vlasov_drag_boundary_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list;

struct dg_fpo_vlasov_drag {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  fpo_vlasov_drag_surf_t surf[3]; // Surface terms for acceleration
  fpo_vlasov_drag_boundary_surf_t boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range phase_range; // Configuration space range.
  struct gkyl_dg_fpo_vlasov_drag_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  return fpo_vlasov_drag_vol_1x3v_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  return fpo_vlasov_drag_vol_1x3v_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
    qIn, qRhsOut);
}

// GKYL_CU_DH
// static double
// kernel_fpo_vlasov_drag_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
//   const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
// {
//   struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
//   long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
//   return fpo_vlasov_drag_vol_2x3v_ser_p1(xc, dx, 
//     (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
//     qIn, qRhsOut);
// }

// GKYL_CU_DH
// static double
// kernel_fpo_vlasov_drag_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
//   const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
// {
//   struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
//   long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
//   return fpo_vlasov_drag_vol_2x3v_ser_p2(xc, dx, 
//     (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
//     qIn, qRhsOut);
// }

// GKYL_CU_DH
// static double
// kernel_fpo_vlasov_drag_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
//   const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
// {
//   struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
//   long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
//   return fpo_vlasov_drag_vol_3x3v_ser_p1(xc, dx, 
//     (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
//     qIn, qRhsOut);
// }

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_vol_kern_list ser_vol_kernels[] = {
  // { NULL, kernel_fpo_vlasov_drag_vol_1x3v_ser_p1, kernel_fpo_vlasov_drag_vol_1x3v_ser_p2 }, // 0
  // { NULL, kernel_fpo_vlasov_drag_vol_2x3v_ser_p1, kernel_fpo_vlasov_drag_vol_2x3v_ser_p2 }, // 1
  // { NULL, kernel_fpo_vlasov_drag_vol_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_vx_kernels[] = {
  // { NULL, fpo_vlasov_drag_surfvx_1x3v_ser_p1, fpo_vlasov_drag_surfvx_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_drag_surfvx_2x3v_ser_p1, fpo_vlasov_drag_surfvx_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_drag_surfvx_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_vy_kernels[] = {
  // { NULL, fpo_vlasov_drag_surfvy_1x3v_ser_p1, fpo_vlasov_drag_surfvy_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_drag_surfvy_2x3v_ser_p1, fpo_vlasov_drag_surfvy_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_drag_surfvy_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_vz_kernels[] = {
  // { NULL, fpo_vlasov_drag_surfvz_1x3v_ser_p1, fpo_vlasov_drag_surfvz_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_drag_surfvz_2x3v_ser_p1, fpo_vlasov_drag_surfvz_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_drag_surfvz_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vx_kernels[] = {
  // { NULL, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_drag_boundary_surfvx_2x3v_ser_p1, fpo_vlasov_drag_boundary_surfvx_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_drag_boundary_surfvx_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vy_kernels[] = {
  // { NULL, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_drag_boundary_surfvy_2x3v_ser_p1, fpo_vlasov_drag_boundary_surfvy_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_drag_boundary_surfvy_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vz_kernels[] = {
  // { NULL, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_drag_boundary_surfvz_2x3v_ser_p1, fpo_vlasov_drag_boundary_surfvz_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_drag_boundary_surfvz_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

void gkyl_fpo_vlasov_drag_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxC);
  if (dir >= fpo_vlasov_drag->cdim) {
    return fpo_vlasov_drag->surf[dir-fpo_vlasov_drag->cdim](xcC, dxC, 
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
      qInL, qInC, qInR, qRhsOut);
  }
  return 0.;
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
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxSkin);
  if (dir >= fpo_vlasov_drag->cdim) {
    return fpo_vlasov_drag->boundary_surf[dir-fpo_vlasov_drag->cdim](xcSkin, dxSkin, 
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.h, pidx), 
      edge, qInSkin, qInEdge, qRhsOut);
  }
  return 0.;
}

