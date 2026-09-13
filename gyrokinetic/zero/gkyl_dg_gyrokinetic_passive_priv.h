// Private header, not for direct use in user code.
#pragma once

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_dg_gyrokinetic_passive_kernels.h>
#include <gkyl_dg_gyrokinetic_passive.h>

typedef double (*dg_gyrokinetic_passive_surf_t)(const double *w, const double *dxv,
  const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
  const double *flux_surf_l, const double *flux_surf_r, double *GKYL_RESTRICT out);

typedef double (*dg_gyrokinetic_passive_boundary_surf_t)(const double *w, const double *dxv,
  const double *vmap_prime_edge, const double *vmap_prime_skin,
  const double *flux_surf_edge, const double *flux_surf_skin,
  const int edge, double *GKYL_RESTRICT out);

// The cv_index_gkp[cd].vdim[vd] maps (cdim, vdim) to a kernel-list index.
static struct { int vdim[3]; } cv_index_gkp[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// For use in kernel tables.
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_gyrokinetic_passive_vol_kern_list;
typedef struct { dg_gyrokinetic_passive_surf_t kernels[3]; } gkyl_dg_gyrokinetic_passive_surf_kern_list;
typedef struct { dg_gyrokinetic_passive_boundary_surf_t kernels[3]; } gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list;

// Private equation struct.
struct dg_gyrokinetic_passive {
  struct gkyl_dg_eqn eqn;           // Base object (must be first).
  int cdim;                         // Config-space dimensions.
  int pdim;                         // Phase-space dimensions.
  dg_gyrokinetic_passive_surf_t surf[GKYL_MAX_CDIM]; // Conf-space surface kernels.
  dg_gyrokinetic_passive_boundary_surf_t boundary_surf[GKYL_MAX_CDIM];
  struct gkyl_range conf_range;     // Configuration space range.
  struct gkyl_range phase_range;    // Phase space range.
  double charge, mass;
  const struct gk_geometry *gk_geom;
  const struct gkyl_velocity_map *vel_map;
  struct gkyl_dg_gyrokinetic_passive_auxfields auxfields;
};

//
// Volume kernel wrappers — these are called by hyper_dg via the vol_term
// function pointer in gkyl_dg_eqn.
//

GKYL_CU_DH
static double
kernel_dg_gyrokinetic_passive_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn,
  const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);

  int vel_idx[2];
  for (int d = gkp->cdim; d < gkp->pdim; d++) vel_idx[d-gkp->cdim] = idx[d];

  long cidx = gkyl_range_idx(&gkp->conf_range, idx);
  long vidx = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idx);

  return dg_gyrokinetic_passive_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap, vidx),
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_sq, vidx),
    gkp->charge, gkp->mass,
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_corn.bmag, cidx),
    (const double*) gkyl_array_cfetch(gkp->auxfields.speeds, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.dualcurlbhatoverB, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.rtg33inv, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.bioverJB, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_gyrokinetic_passive_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn,
  const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);

  int vel_idx[2];
  for (int d = gkp->cdim; d < gkp->pdim; d++) vel_idx[d-gkp->cdim] = idx[d];

  long cidx = gkyl_range_idx(&gkp->conf_range, idx);
  long vidx = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idx);

  return dg_gyrokinetic_passive_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap, vidx),
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_sq, vidx),
    gkp->charge, gkp->mass,
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_corn.bmag, cidx),
    (const double*) gkyl_array_cfetch(gkp->auxfields.speeds, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.dualcurlbhatoverB, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.rtg33inv, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.bioverJB, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_gyrokinetic_passive_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn,
  const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);

  int vel_idx[2];
  for (int d = gkp->cdim; d < gkp->pdim; d++) vel_idx[d-gkp->cdim] = idx[d];

  long cidx = gkyl_range_idx(&gkp->conf_range, idx);
  long vidx = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idx);

  return dg_gyrokinetic_passive_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap, vidx),
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_sq, vidx),
    gkp->charge, gkp->mass,
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_corn.bmag, cidx),
    (const double*) gkyl_array_cfetch(gkp->auxfields.speeds, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.dualcurlbhatoverB, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.rtg33inv, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.bioverJB, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_gyrokinetic_passive_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn,
  const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);

  int vel_idx[2];
  for (int d = gkp->cdim; d < gkp->pdim; d++) vel_idx[d-gkp->cdim] = idx[d];

  long cidx = gkyl_range_idx(&gkp->conf_range, idx);
  long vidx = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idx);

  return dg_gyrokinetic_passive_vol_3x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap, vidx),
    (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_sq, vidx),
    gkp->charge, gkp->mass,
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_corn.bmag, cidx),
    (const double*) gkyl_array_cfetch(gkp->auxfields.speeds, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.dualcurlbhatoverB, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.rtg33inv, cidx),
    (const double*) gkyl_array_cfetch(gkp->gk_geom->geo_int.bioverJB, cidx),
    qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_vol_kern_list ser_passive_vol_kernels[] = {
  { kernel_dg_gyrokinetic_passive_vol_1x1v_ser_p1, NULL },
  { kernel_dg_gyrokinetic_passive_vol_1x2v_ser_p1, NULL },
  { kernel_dg_gyrokinetic_passive_vol_2x2v_ser_p1, NULL },
  { kernel_dg_gyrokinetic_passive_vol_3x2v_ser_p1, NULL },
};

//
// Surface kernel lists — conf-space directions only.
//

// x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_surf_kern_list ser_passive_surf_x_kernels[] = {
  { dg_gyrokinetic_passive_surfx_1x1v_ser_p1, NULL },
  { dg_gyrokinetic_passive_surfx_1x2v_ser_p1, NULL },
  { dg_gyrokinetic_passive_surfx_2x2v_ser_p1, NULL },
  { dg_gyrokinetic_passive_surfx_3x2v_ser_p1, NULL },
};

// y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_surf_kern_list ser_passive_surf_y_kernels[] = {
  { NULL, NULL },
  { NULL, NULL },
  { dg_gyrokinetic_passive_surfy_2x2v_ser_p1, NULL },
  { dg_gyrokinetic_passive_surfy_3x2v_ser_p1, NULL },
};

// z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_surf_kern_list ser_passive_surf_z_kernels[] = {
  { NULL, NULL },
  { NULL, NULL },
  { NULL, NULL },
  { dg_gyrokinetic_passive_surfz_3x2v_ser_p1, NULL },
};

//
// Boundary-surface kernel lists.
//

// x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list ser_passive_boundary_surf_x_kernels[] = {
  { dg_gyrokinetic_passive_boundary_surfx_1x1v_ser_p1, NULL },
  { dg_gyrokinetic_passive_boundary_surfx_1x2v_ser_p1, NULL },
  { dg_gyrokinetic_passive_boundary_surfx_2x2v_ser_p1, NULL },
  { dg_gyrokinetic_passive_boundary_surfx_3x2v_ser_p1, NULL },
};

// y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list ser_passive_boundary_surf_y_kernels[] = {
  { NULL, NULL },
  { NULL, NULL },
  { dg_gyrokinetic_passive_boundary_surfy_2x2v_ser_p1, NULL },
  { dg_gyrokinetic_passive_boundary_surfy_3x2v_ser_p1, NULL },
};

// z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_passive_boundary_surf_kern_list ser_passive_boundary_surf_z_kernels[] = {
  { NULL, NULL },
  { NULL, NULL },
  { NULL, NULL },
  { dg_gyrokinetic_passive_boundary_surfz_3x2v_ser_p1, NULL },
};

// Macro to select a kernel from a list.
#define CK_PASSIVE(lst,cdim,vd,poly_order) lst[cv_index_gkp[cdim].vdim[vd]].kernels[poly_order-1]

/**
 * Free passive gyrokinetic equation object.
 */
void gkyl_gyrokinetic_passive_free(const struct gkyl_ref_count *ref);

// Surface dispatch: called by hyper_dg for each surface.
GKYL_CU_D
static double
surf_passive(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR,
  const int *idxL, const int *idxC, const int *idxR,
  const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);

  // Only update conf-space directions (dir < cdim).
  if (dir < gkp->cdim) {
    int vel_idxL[2], vel_idxC[2], vel_idxR[2];
    for (int d = gkp->cdim; d < gkp->pdim; d++) {
      vel_idxL[d-gkp->cdim] = idxL[d];
      vel_idxC[d-gkp->cdim] = idxC[d];
      vel_idxR[d-gkp->cdim] = idxR[d];
    }
    long vidxL = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idxL);
    long vidxC = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idxC);
    long vidxR = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idxR);

    long pidxC = gkyl_range_idx(&gkp->phase_range, idxC);
    long pidxR = gkyl_range_idx(&gkp->phase_range, idxR);

    return gkp->surf[dir](xcC, dxC,
      (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_prime, vidxL),
      (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_prime, vidxC),
      (const double*) gkyl_array_cfetch(gkp->vel_map->vmap_prime, vidxR),
      (const double*) gkyl_array_cfetch(gkp->auxfields.flux_surf, pidxC),
      (const double*) gkyl_array_cfetch(gkp->auxfields.flux_surf, pidxR),
      qRhsOut);
  }
  return 0.0;
}

// Boundary-surface dispatch.
GKYL_CU_D
static double
boundary_surf_passive(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin,
  const int *idxEdge, const int *idxSkin, const int edge,
  const double *qInEdge, const double *qInSkin,
  double *GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);

  if (dir < gkp->cdim) {
    int vel_idxEdge[2], vel_idxSkin[2];
    for (int d = gkp->cdim; d < gkp->pdim; d++) {
      vel_idxEdge[d-gkp->cdim] = idxEdge[d];
      vel_idxSkin[d-gkp->cdim] = idxSkin[d];
    }
    long vidxEdge = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idxEdge);
    long vidxSkin = gkyl_range_idx(&gkp->vel_map->local_vel, vel_idxSkin);

    long pidxEdge = gkyl_range_idx(&gkp->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&gkp->phase_range, idxSkin);

    return gkp->boundary_surf[dir](xcSkin, dxSkin,
      (const double *) gkyl_array_cfetch(gkp->vel_map->vmap_prime, vidxEdge),
      (const double *) gkyl_array_cfetch(gkp->vel_map->vmap_prime, vidxSkin),
      (const double *) gkyl_array_cfetch(gkp->auxfields.flux_surf, pidxEdge),
      (const double *) gkyl_array_cfetch(gkp->auxfields.flux_surf, pidxSkin),
      edge, qRhsOut);
  }
  return 0.0;
}

// boundary_diag mirrors boundary_surf.
GKYL_CU_D
static double
boundary_diag_passive(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin,
  const int *idxEdge, const int *idxSkin, const int edge,
  const double *qInEdge, const double *qInSkin,
  double *GKYL_RESTRICT qRhsOut)
{
  return boundary_surf_passive(eqn, dir, xcEdge, xcSkin, dxEdge, dxSkin,
    idxEdge, idxSkin, edge, qInEdge, qInSkin, qRhsOut);
}

#ifdef GKYL_HAVE_CUDA
// GPU-only functions (defined in dg_gyrokinetic_passive_cu.cu).
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_passive_cu_dev_new(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map);

void gkyl_gyrokinetic_passive_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_dg_gyrokinetic_passive_auxfields auxin);
#endif
