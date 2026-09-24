#pragma once

// Private header, not for direct use in user code

#include <gkyl_lbo_pkpm_kernels.h>

// Types for various kernels
typedef double (*lbo_pkpm_diff_surf_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuPrimMomsSum,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*lbo_pkpm_diff_boundary_surf_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuPrimMomsSum,
  const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);


// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_lbo_pkpm_diff_vol_kern_list;
typedef struct { lbo_pkpm_diff_surf_t kernels[3]; } gkyl_dg_lbo_pkpm_diff_surf_kern_list;
typedef struct { lbo_pkpm_diff_boundary_surf_t kernels[3]; } gkyl_dg_lbo_pkpm_diff_boundary_surf_kern_list;

struct dg_lbo_pkpm_diff {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  lbo_pkpm_diff_surf_t surf; // Surface terms for acceleration
  lbo_pkpm_diff_boundary_surf_t boundary_surf; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_lbo_pkpm_diff_auxfields auxfields; // Auxiliary fields.
  double vMaxSq;
  int num_cbasis;
};

GKYL_CU_DH
static inline bool
checkPrimMomCross(struct dg_lbo_pkpm_diff *lbo_pkpm_diff,
  const double* nuSum_p, const double* nuVtSqSum_p) {
  bool noPrimMomCross = true;
  noPrimMomCross = noPrimMomCross && ((nuVtSqSum_p[0]>0.)
    && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_pkpm_diff->vMaxSq));
  return noPrimMomCross;
}

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_lbo_pkpm_diff_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_pkpm_diff_vol_1x1v_ser_p1(xc, dx, nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_pkpm_diff_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_pkpm_diff_vol_1x1v_ser_p2(xc, dx, nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_pkpm_diff_vol_1x1v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_pkpm_diff_vol_1x1v_tensor_p2(xc, dx, nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_pkpm_diff_vol_2x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_pkpm_diff_vol_2x1v_ser_p1(xc, dx, nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_pkpm_diff_vol_2x1v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_pkpm_diff_vol_2x1v_tensor_p2(xc, dx, nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_pkpm_diff_vol_3x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if (noPrimMomCross) {
    return lbo_pkpm_diff_vol_3x1v_ser_p1(xc, dx, nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_lbo_pkpm_diff_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_pkpm_diff_vol_1x1v_ser_p1, kernel_lbo_pkpm_diff_vol_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_lbo_pkpm_diff_vol_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, kernel_lbo_pkpm_diff_vol_3x1v_ser_p1, NULL }, // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_lbo_pkpm_diff_vol_kern_list ten_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_pkpm_diff_vol_1x1v_ser_p1, kernel_lbo_pkpm_diff_vol_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, kernel_lbo_pkpm_diff_vol_2x1v_ser_p1, kernel_lbo_pkpm_diff_vol_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, kernel_lbo_pkpm_diff_vol_3x1v_ser_p1, NULL }, // 2
};

// Constant nu surface kernel list: vpar-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_lbo_pkpm_diff_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_pkpm_diff_surfvpar_1x1v_ser_p1, lbo_pkpm_diff_surfvpar_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, lbo_pkpm_diff_surfvpar_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, lbo_pkpm_diff_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Constant nu surface kernel list: vpar-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_lbo_pkpm_diff_surf_kern_list ten_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_pkpm_diff_surfvpar_1x1v_ser_p1, lbo_pkpm_diff_surfvpar_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, lbo_pkpm_diff_surfvpar_2x1v_ser_p1, lbo_pkpm_diff_surfvpar_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, lbo_pkpm_diff_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vpar-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_lbo_pkpm_diff_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_pkpm_diff_boundary_surfvpar_1x1v_ser_p1, lbo_pkpm_diff_boundary_surfvpar_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, lbo_pkpm_diff_boundary_surfvpar_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, lbo_pkpm_diff_boundary_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vpar-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_lbo_pkpm_diff_boundary_surf_kern_list ten_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_pkpm_diff_boundary_surfvpar_1x1v_ser_p1, lbo_pkpm_diff_boundary_surfvpar_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, lbo_pkpm_diff_boundary_surfvpar_2x1v_ser_p1, lbo_pkpm_diff_boundary_surfvpar_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, lbo_pkpm_diff_boundary_surfvpar_3x1v_ser_p1, NULL }, // 2
};

void gkyl_lbo_pkpm_diff_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idxC);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if ((dir >= lbo_pkpm_diff->cdim) && (noPrimMomCross)) {
    return lbo_pkpm_diff->surf(xcC, dxC, 
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
  struct dg_lbo_pkpm_diff *lbo_pkpm_diff = container_of(eqn, struct dg_lbo_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_pkpm_diff->conf_range, idxSkin);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_pkpm_diff->auxfields.nuPrimMomsSum, cidx);
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_pkpm_diff->num_cbasis];
  bool noPrimMomCross = checkPrimMomCross(lbo_pkpm_diff, nuSum_p, nuVtSqSum_p);
  if ((dir >= lbo_pkpm_diff->cdim) && (noPrimMomCross)) {
    return lbo_pkpm_diff->boundary_surf(xcSkin, dxSkin, 
      nuSum_p, nuPrimMomsSum_p, edge, qInSkin, qInEdge, qRhsOut);
  }
  return 0.;
}

