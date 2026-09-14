#pragma once

// Private header, not for direct use in user code.

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_dg_gyrokinetic_kernels.h>
#include <gkyl_dg_gyrokinetic.h>

// Types for various kernels.

typedef double (*dg_gyrokinetic_vol_es_t)(const double *w, const double *dxv,
  const double *vmap, const double *vmapSq, const double q_, const double m_,
  const double *bmag, const double *phi,
  const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
  const double *fin, double* GKYL_RESTRICT out); 

typedef double (*dg_gyrokinetic_vol_add_apar_t)(const double *w, const double *dxv, 
  const double *vmap, const double *vmapSq, const double q_, const double m_, 
  const double *bmag, const double *jacobtot_inv, const double *dualcurlbhatoverB, const double *bioverJB, const double *b_i, 
  const double *phi, const double *apar, const double *fin, double* GKYL_RESTRICT out); 

typedef double (*dg_gyrokinetic_vol_add_apardot_t)(const double *vmap, 
  const double q_, const double m_,
  const double *apardot, const double *fin, double* GKYL_RESTRICT out); 

typedef double (*dg_gyrokinetic_surf_t)( const double *w, const double *dxv,
  const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
  const double *flux_surf_l, const double *flux_surf_r, double* GKYL_RESTRICT out
); 

typedef double (*dg_gyrokinetic_boundary_surf_t)(const double *w, const double *dxv,
  const double *vmap_prime_edge, const double *vmap_prime_skin,
  const double *flux_surf_edge, const double *flux_surf_skin, 
  const int edge, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_gyrokinetic_vol_kern_list;
typedef struct { dg_gyrokinetic_vol_es_t kernels[3]; } gkyl_dg_gyrokinetic_vol_es_kern_list;
typedef struct { dg_gyrokinetic_vol_add_apar_t kernels[3]; } gkyl_dg_gyrokinetic_vol_add_apar_kern_list;
typedef struct { dg_gyrokinetic_vol_add_apardot_t kernels[3]; } gkyl_dg_gyrokinetic_vol_add_apardot_kern_list;
typedef struct { dg_gyrokinetic_surf_t kernels[3]; } gkyl_dg_gyrokinetic_surf_kern_list;
typedef struct { dg_gyrokinetic_boundary_surf_t kernels[3]; } gkyl_dg_gyrokinetic_boundary_surf_kern_list;

struct dg_gyrokinetic {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  dg_gyrokinetic_vol_es_t vol_es_kernel; // Electrostatic contribution to volume term.
  dg_gyrokinetic_vol_add_apar_t vol_add_apar_kernel; // Additional volume term for Apar contribution.
  dg_gyrokinetic_vol_add_apardot_t vol_add_apardot_kernel; // Additional volume term for Apar contribution.
  dg_gyrokinetic_surf_t surf[4]; // Surface terms.
  dg_gyrokinetic_boundary_surf_t boundary_surf[4]; // Surface terms for velocity boundary.
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_range phase_range; // Phase space range.
  double charge, mass;
  const struct gk_geometry *gk_geom; // Pointer to geometry struct
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  struct gkyl_dg_gyrokinetic_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels general geometry
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_dg_gyrokinetic_vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  int vel_idx[2];
  for (int d=gyrokinetic->cdim; d<gyrokinetic->pdim; d++) vel_idx[d-gyrokinetic->cdim] = idx[d];

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long vidx = gkyl_range_idx(&gyrokinetic->vel_map->local_vel, vel_idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  double cfl = 0.0;

  const double *phi = (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx);
  const double *apar = (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx);
  const double *apardot = (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx);
  const double *vm = (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap, vidx);
  const double *vm_sq = (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap_sq, vidx);
  const double *bmag = (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->geo_corn.bmag, cidx);
  // The jacobtot_inv is not used anymore.
  const double *jacobtot_inv = (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->geo_int.jacobtot_inv, cidx);
  const double *b_i = (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->geo_int.b_i, cidx);
  const double *dualcurlbhatoverB = (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->geo_int.dualcurlbhatoverB, cidx);
  const double *rtg33inv = (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->geo_int.rtg33inv, cidx);
  const double *bioverJB = (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->geo_int.bioverJB, cidx);
  
  // Electrostatic contribution to volume term.
  cfl += gyrokinetic->vol_es_kernel(xc, dx, vm, vm_sq, gyrokinetic->charge, gyrokinetic->mass, 
    bmag, phi, dualcurlbhatoverB, rtg33inv, bioverJB, qIn, qRhsOut);
  // Add Apar contribution to volume term if needed.
  cfl += gyrokinetic->vol_add_apar_kernel(xc, dx, vm, vm_sq, gyrokinetic->charge, gyrokinetic->mass, 
    bmag, jacobtot_inv, dualcurlbhatoverB, bioverJB, b_i, phi, apar, qIn, qRhsOut);
  // Add Apardot contribution to volume term if needed.
  cfl += gyrokinetic->vol_add_apardot_kernel( vm, gyrokinetic->charge, gyrokinetic->mass, 
    apardot, qIn, qRhsOut);
  return cfl;
}

// To turn off terms (used in EM).

GKYL_CU_DH 
static double 
dg_gyrokinetic_vol_none(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *phi,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out)
{ 
  return 0.; 
}

GKYL_CU_DH 
static double 
dg_gyrokinetic_add_apar_vol_none(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv, 
    const double *dualcurlbhatoverB, const double *bioverJB, const double *b_i, const double *phi, 
    const double *apar, const double *fin, double* GKYL_RESTRICT out) 
{ 
  return 0.; 
}

GKYL_CU_DH 
static double 
dg_gyrokinetic_add_apardot_vol_none(const double *vmap, const double q_, const double m_,
            const double *apardot, const double *fin, double* GKYL_RESTRICT out)
{ 
  return 0.; 
} 


// Volume electrostatic kernel list.
GKYL_CU_D
static const gkyl_dg_gyrokinetic_vol_es_kern_list ser_vol_es_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_vol_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_vol_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_vol_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_vol_3x2v_ser_p1, NULL }, // 3
};

//
// Serendipity volume kernels general geometry, no toroidal field (by=0)
// Need to be separated like this for GPU build
//
GKYL_CU_D
static const gkyl_dg_gyrokinetic_vol_es_kern_list ser_no_by_vol_es_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_vol_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_vol_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_no_by_vol_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_no_by_vol_3x2v_ser_p1, NULL }, // 3
};

// Electromagnetic serendipity volume kernels general geometry.
GKYL_CU_D
static const gkyl_dg_gyrokinetic_vol_add_apar_kern_list ser_add_apar_vol_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_add_apar_vol_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_add_apar_vol_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_add_apar_vol_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_add_apar_vol_3x2v_ser_p1, NULL }, // 3
};

// Electromagnetic volume kernel list (apardot).
GKYL_CU_D
static const gkyl_dg_gyrokinetic_vol_add_apardot_kern_list ser_add_apardot_vol_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_add_apardot_vol_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_add_apardot_vol_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_add_apardot_vol_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_add_apardot_vol_3x2v_ser_p1, NULL }, // 3
};

//
// Serendipity surface kernels general geometry
//
// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_x_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_surfx_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_surfx_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_surfx_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_surfy_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_surfy_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_surfz_3x2v_ser_p1, NULL }, // 3
};

// Acceleration surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_surfvpar_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_surfvpar_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_boundary_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_boundary_surfx_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_boundary_surfx_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_boundary_surfx_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_boundary_surfy_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_boundary_surfy_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_boundary_surfz_3x2v_ser_p1, NULL }, // 3
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, dg_gyrokinetic_boundary_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, dg_gyrokinetic_boundary_surfvpar_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, dg_gyrokinetic_boundary_surfvpar_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, dg_gyrokinetic_boundary_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

/**
 * Free gyrokinetic eqn object.
 *
 * @param ref Reference counter for gyrokinetic eqn
 */
void gkyl_gyrokinetic_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcL, const double* xcC, const double* xcR,
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  // Only in x,y,z,vpar directions.
  if (dir <= gyrokinetic->cdim) {
    int vel_idxL[2], vel_idxC[2], vel_idxR[2];
    for (int d=gyrokinetic->cdim; d<gyrokinetic->pdim; d++) {
      vel_idxL[d-gyrokinetic->cdim] = idxL[d];
      vel_idxC[d-gyrokinetic->cdim] = idxC[d];
      vel_idxR[d-gyrokinetic->cdim] = idxR[d];
    }
    long vidxL = gkyl_range_idx(&gyrokinetic->vel_map->local_vel, vel_idxL);
    long vidxC = gkyl_range_idx(&gyrokinetic->vel_map->local_vel, vel_idxC);
    long vidxR = gkyl_range_idx(&gyrokinetic->vel_map->local_vel, vel_idxR);

    // Each cell owns the *lower* edge surface alpha
    // Since alpha is continuous, fetch alpha_surf in center cell for lower edge
    // and fetch alpha_surf in right cell for upper edge
    long pidxC = gkyl_range_idx(&gyrokinetic->phase_range, idxC);
    long pidxR = gkyl_range_idx(&gyrokinetic->phase_range, idxR);
    return gyrokinetic->surf[dir](xcC, dxC, 
      (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap_prime, vidxL),
      (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap_prime, vidxC),
      (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap_prime, vidxR),
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.flux_surf, pidxC), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.flux_surf, pidxR), 
      qRhsOut);
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
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  // Only in x,y,z,vpar directions.
  if (dir <= gyrokinetic->cdim) {
    int vel_idxEdge[2], vel_idxSkin[2];
    for (int d=gyrokinetic->cdim; d<gyrokinetic->pdim; d++) {
      vel_idxEdge[d-gyrokinetic->cdim] = idxEdge[d];
      vel_idxSkin[d-gyrokinetic->cdim] = idxSkin[d];
    }
    long vidxEdge = gkyl_range_idx(&gyrokinetic->vel_map->local_vel, vel_idxEdge);
    long vidxSkin = gkyl_range_idx(&gyrokinetic->vel_map->local_vel, vel_idxSkin);

    // Each cell owns the *lower* edge surface alpha
    long pidxEdge = gkyl_range_idx(&gyrokinetic->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&gyrokinetic->phase_range, idxSkin);
    return gyrokinetic->boundary_surf[dir](xcSkin, dxSkin, 
      (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap_prime, vidxEdge),
      (const double*) gkyl_array_cfetch(gyrokinetic->vel_map->vmap_prime, vidxSkin),
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.flux_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.flux_surf, pidxSkin), 
      edge, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_diag(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return boundary_surf(eqn, dir, xcEdge, xcSkin, dxEdge, dxSkin, idxEdge, idxSkin, edge, qInEdge, qInSkin, qRhsOut);
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create new Gyrokinetic equation object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_cu_dev_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const double charge, const double mass, enum gkyl_gk_collisionless_type collless_type,
  const bool no_by, const bool complete_em, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map);

/**
 * CUDA device function to set the auxiliary fields (e.g. geometry & EM fields)
 * needed in computing gyrokinetic updates.
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gyrokinetic_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin);
#endif
