#pragma once

// Private header for the gyrokinetic LBO drag equation object.
// Not for direct use in user code.

#include <gkyl_gk_geometry.h>
#include <gkyl_lbo_gyrokinetic_kernels.h>

// Types for various kernels
typedef double (*lbo_gyrokinetic_drag_surf_t)(const double *dxv, const double *vmap,
  const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, const double m_,
  const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*lbo_gyrokinetic_drag_boundary_surf_t)(const double *dxv, const double *vmap,
  const double *vmap_prime_edge, const double *vmap_prime_skin, const double m_,
  const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum,
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list;
typedef struct { lbo_gyrokinetic_drag_surf_t kernels[3]; } gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list;
typedef struct { lbo_gyrokinetic_drag_boundary_surf_t kernels[3]; } gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list;

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst, cdim, vd, poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_lbo_gyrokinetic_drag {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  lbo_gyrokinetic_drag_surf_t surf[2]; // Surface terms for acceleration.
  lbo_gyrokinetic_drag_boundary_surf_t boundary_surf[2]; // Surface terms for acceleration.
  struct gkyl_range conf_range; // Configuration space range.
  double mass; // Species mass.
  const struct gk_geometry *gk_geom; // Pointer to geometry struct
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxfields; // Auxiliary fields.
  double vparMax, vparMaxSq;
  int num_cbasis;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  int vel_idx[2];
  for (int d=lbo->cdim; d<lbo->pdim; d++) vel_idx[d-lbo->cdim] = idx[d];

  long cidx = gkyl_range_idx(&lbo->conf_range, idx);
  long vidx = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo->vparMaxSq) &&
      (m2self_p[0]>0.)) {
    return lbo_gyrokinetic_drag_vol_1x1v_ser_p1(dx,
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap, vidx),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidx), lbo->mass, 
      (const double*) gkyl_array_cfetch(lbo->gk_geom->geo_corn.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  int vel_idx[2];
  for (int d=lbo->cdim; d<lbo->pdim; d++) vel_idx[d-lbo->cdim] = idx[d];

  long cidx = gkyl_range_idx(&lbo->conf_range, idx);
  long vidx = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_1x2v_ser_p1(dx,
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap, vidx),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidx), lbo->mass, 
      (const double*) gkyl_array_cfetch(lbo->gk_geom->geo_corn.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  int vel_idx[2];
  for (int d=lbo->cdim; d<lbo->pdim; d++) vel_idx[d-lbo->cdim] = idx[d];

  long cidx = gkyl_range_idx(&lbo->conf_range, idx);
  long vidx = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_2x2v_ser_p1(dx, 
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap, vidx),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidx), lbo->mass, 
      (const double*) gkyl_array_cfetch(lbo->gk_geom->geo_corn.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  int vel_idx[2];
  for (int d=lbo->cdim; d<lbo->pdim; d++) vel_idx[d-lbo->cdim] = idx[d];

  long cidx = gkyl_range_idx(&lbo->conf_range, idx);
  long vidx = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idx);

  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_3x2v_ser_p1(dx,
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap, vidx),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidx), lbo->mass, 
      (const double*) gkyl_array_cfetch(lbo->gk_geom->geo_corn.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_gyrokinetic_drag_vol_1x1v_ser_p1, NULL }, // 0
  { NULL, kernel_lbo_gyrokinetic_drag_vol_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, kernel_lbo_gyrokinetic_drag_vol_2x2v_ser_p1, NULL }, // 3
  // 3x kernels
  { NULL, kernel_lbo_gyrokinetic_drag_vol_3x2v_ser_p1, NULL }, // 4
};

// Surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_gyrokinetic_drag_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, lbo_gyrokinetic_drag_surfvpar_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: mu-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list ser_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_gyrokinetic_drag_surfmu_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_surfmu_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_surfmu_3x2v_ser_p1, NULL }, // 3
};

// Boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Constant nu boundary surface kernel (zero-flux BCs) list: mu-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list ser_boundary_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p1, NULL }, // 3
};

void gkyl_lbo_gyrokinetic_drag_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  long cidx = gkyl_range_idx(&lbo->conf_range, idxC);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo->num_cbasis];
  if ((dir >= lbo->cdim) &&
      (fabs(nuUSum_p[0]/nuSum_p[0]) < lbo->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo->vparMaxSq) &&
      (m2self_p[0]>0.))
  {
    int vel_idxL[2], vel_idxC[2], vel_idxR[2];
    for (int d=lbo->cdim; d<lbo->pdim; d++) {
      vel_idxL[d-lbo->cdim] = idxL[d];
      vel_idxC[d-lbo->cdim] = idxC[d];
      vel_idxR[d-lbo->cdim] = idxR[d];
    }
    long vidxL = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idxL);
    long vidxC = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idxC);
    long vidxR = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idxR);

    return lbo->surf[dir-lbo->cdim](dxC,
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap, vidxC),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidxL),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidxC),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidxR), lbo->mass,
      (const double*) gkyl_array_cfetch(lbo->gk_geom->geo_corn.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qInL, qInC, qInR, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  long cidx = gkyl_range_idx(&lbo->conf_range, idxSkin);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo->auxfields.nuPrimMomsSum, cidx); 
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo->num_cbasis];
  if ((dir >= lbo->cdim) &&
      (fabs(nuUSum_p[0]/nuSum_p[0]) < lbo->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo->vparMaxSq) &&
      (m2self_p[0]>0.))
  {
    int vel_idxEdge[2], vel_idxSkin[2];
    for (int d=lbo->cdim; d<lbo->pdim; d++) {
      vel_idxEdge[d-lbo->cdim] = idxEdge[d];
      vel_idxSkin[d-lbo->cdim] = idxSkin[d];
    }
    long vidxEdge = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idxEdge);
    long vidxSkin = gkyl_range_idx(&lbo->vel_map->local_vel, vel_idxSkin);

    return lbo->boundary_surf[dir-lbo->cdim](dxSkin, 
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap, vidxSkin),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidxEdge),
      (const double*) gkyl_array_cfetch(lbo->vel_map->vmap_prime, vidxSkin), lbo->mass,
      (const double*) gkyl_array_cfetch(lbo->gk_geom->geo_corn.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, edge, qInEdge, qInSkin, qRhsOut);
  }
  return 0.;
}

#ifdef GKYL_HAVE_CUDA

/**
 * Create a new LBO equation object that lives on NV-GPU
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map);

/**
 * CUDA device function to set auxiliary fields needed in updating the drag flux term.
 */
void gkyl_lbo_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxin);

#endif

