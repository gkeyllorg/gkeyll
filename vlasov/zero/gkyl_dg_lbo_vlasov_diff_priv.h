#pragma once

// Private header, not for direct use in user code

#include <gkyl_lbo_vlasov_kernels.h>

// Types for various kernels
typedef double (*lbo_vlasov_diff_surf_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuPrimMomsSum,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*lbo_vlasov_diff_boundary_surf_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuPrimMomsSum,
  const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_lbo_vlasov_diff_vol_kern_list;
typedef struct { lbo_vlasov_diff_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_diff_surf_kern_list;
typedef struct { lbo_vlasov_diff_boundary_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_diff_boundary_surf_kern_list;

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst, cdim, vd, poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_lbo_vlasov_diff {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int vdim; // Velocity-space dimensions.
  int pdim; // Phase-space dimensions
  lbo_vlasov_diff_surf_t surf[3]; // Surface terms for acceleration
  lbo_vlasov_diff_boundary_surf_t boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_lbo_vlasov_diff_auxfields auxfields; // Auxiliary fields.
  double viMax[3], vMaxSq;
  int num_cbasis;
};

GKYL_CU_DH
static inline bool
checkPrimMomCross(struct dg_lbo_vlasov_diff *lbo_vlasov_diff,
  const double* nuSum_p, const double* nuUSum_p, const double* nuVtSqSum_p) {
  bool noPrimMomCross = true;
  for (int d=0; d<lbo_vlasov_diff->vdim; d++) {
    if (fabs(nuUSum_p[d*lbo_vlasov_diff->num_cbasis]/nuSum_p[0]) > lbo_vlasov_diff->viMax[d]) {
       noPrimMomCross = false;
       break;
    }
  }
  noPrimMomCross = noPrimMomCross && ((nuVtSqSum_p[0]>0.)
    && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_vlasov_diff->vMaxSq));
  return noPrimMomCross;
}

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_1x1v_ser_p1(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_1x1v_ser_p2(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_1x2v_ser_p1(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_1x2v_ser_p2(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_1x3v_ser_p1(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_1x3v_ser_p2(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_2x2v_ser_p1(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_2x2v_ser_p2(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_2x3v_ser_p1(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_2x3v_ser_p2(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_diff_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_vlasov_diff_vol_3x3v_ser_p1(xc, dx,
        nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_vlasov_diff_vol_1x1v_ser_p1, kernel_lbo_vlasov_diff_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_lbo_vlasov_diff_vol_1x2v_ser_p1, kernel_lbo_vlasov_diff_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_lbo_vlasov_diff_vol_1x3v_ser_p1, kernel_lbo_vlasov_diff_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_lbo_vlasov_diff_vol_2x2v_ser_p1, kernel_lbo_vlasov_diff_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_lbo_vlasov_diff_vol_2x3v_ser_p1, kernel_lbo_vlasov_diff_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_lbo_vlasov_diff_vol_3x3v_ser_p1, NULL               }, // 5
};

// Constant nu surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_surf_kern_list ser_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_diff_surfvx_1x1v_ser_p1, lbo_vlasov_diff_surfvx_1x1v_ser_p2 }, // 0
  { NULL, lbo_vlasov_diff_surfvx_1x2v_ser_p1, lbo_vlasov_diff_surfvx_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_diff_surfvx_1x3v_ser_p1, lbo_vlasov_diff_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_diff_surfvx_2x2v_ser_p1, lbo_vlasov_diff_surfvx_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_diff_surfvx_2x3v_ser_p1, lbo_vlasov_diff_surfvx_2x3v_ser_p2 }, // 
  // 3x kernels
  { NULL, lbo_vlasov_diff_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_surf_kern_list ser_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_vlasov_diff_surfvy_1x2v_ser_p1, lbo_vlasov_diff_surfvy_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_diff_surfvy_1x3v_ser_p1, lbo_vlasov_diff_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_diff_surfvy_2x2v_ser_p1, lbo_vlasov_diff_surfvy_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_diff_surfvy_2x3v_ser_p1, lbo_vlasov_diff_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_diff_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_surf_kern_list ser_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, lbo_vlasov_diff_surfvz_1x3v_ser_p1, lbo_vlasov_diff_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, lbo_vlasov_diff_surfvz_2x3v_ser_p1, lbo_vlasov_diff_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_diff_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_diff_boundary_surfvx_1x1v_ser_p1, lbo_vlasov_diff_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, lbo_vlasov_diff_boundary_surfvx_1x2v_ser_p1, lbo_vlasov_diff_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_diff_boundary_surfvx_1x3v_ser_p1, lbo_vlasov_diff_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_diff_boundary_surfvx_2x2v_ser_p1, lbo_vlasov_diff_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_diff_boundary_surfvx_2x3v_ser_p1, lbo_vlasov_diff_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_diff_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_vlasov_diff_boundary_surfvy_1x2v_ser_p1, lbo_vlasov_diff_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_diff_boundary_surfvy_1x3v_ser_p1, lbo_vlasov_diff_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_diff_boundary_surfvy_2x2v_ser_p1, lbo_vlasov_diff_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_diff_boundary_surfvy_2x3v_ser_p1, lbo_vlasov_diff_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_diff_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, lbo_vlasov_diff_boundary_surfvz_1x3v_ser_p1, lbo_vlasov_diff_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, lbo_vlasov_diff_boundary_surfvz_2x3v_ser_p1, lbo_vlasov_diff_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_diff_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

void gkyl_lbo_vlasov_diff_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idxC);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if ((dir >= lbo_vlasov_diff->cdim) && (noPrimMomCross)) {
    return lbo_vlasov_diff->surf[dir-lbo_vlasov_diff->cdim](xcC, dxC, 
      nuSum_p, nuPrimMomsSum_p, qInL, qInC, qInR, qRhsOut);
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
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_diff->conf_range, idxSkin);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_vlasov_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_vlasov_diff->vdim*lbo_vlasov_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_diff, nuSum_p, nuUSum_p, nuVtSqSum_p);
  if ((dir >= lbo_vlasov_diff->cdim) && (noPrimMomCross)) {
    return lbo_vlasov_diff->boundary_surf[dir-lbo_vlasov_diff->cdim](xcSkin, dxSkin, 
      nuSum_p, nuPrimMomsSum_p, edge, qInSkin, qInEdge, qRhsOut);
  }
  return 0.;
}

