#pragma once

// Private header for the gyrokinetic RAD drag equation object.
// Not for direct use in user code.

#include <gkyl_rad_gyrokinetic_kernels.h>

// Types for various kernels
typedef double (*rad_gyrokinetic_surf_t)(const double *w, const double *dxv, const double *vmap,
  const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
  const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
  const double *fl, const double *fc, const double *fr, 
  double* GKYL_RESTRICT out);

typedef double (*rad_gyrokinetic_boundary_surf_t)(const double *w, const double *dxv, const double *vmap,
  const double *vmap_prime_edge, const double *vmap_prime_skin,
  const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_rad_gyrokinetic_vol_kern_list;
typedef struct { rad_gyrokinetic_surf_t kernels[3]; } gkyl_dg_rad_gyrokinetic_surf_kern_list;
typedef struct { rad_gyrokinetic_boundary_surf_t kernels[3]; } gkyl_dg_rad_gyrokinetic_boundary_surf_kern_list;

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst, cdim, vd, poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_rad_gyrokinetic_drag {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  double cellav_norm_conf; // Norm factor to multiply 0th DG coeff to get cell ave.
  rad_gyrokinetic_surf_t surf[2]; // Surface terms for acceleration.
  rad_gyrokinetic_boundary_surf_t boundary_surf[2]; // Surface terms for acceleration.
  struct gkyl_range phase_range; // Phase-space range.
  struct gkyl_range conf_range; // Configuration-space range.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  struct gkyl_dg_rad_gyrokinetic_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *grad_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);

  long cidx = gkyl_range_idx(&grad_drag->conf_range, idx);
  int vel_idx[2];
  for (int d=grad_drag->cdim; d<grad_drag->pdim; d++) vel_idx[d-grad_drag->cdim] = idx[d];
  
  long vidx = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idx);
  long pidx = gkyl_range_idx(&grad_drag->phase_range, idx);
  
  return rad_gyrokinetic_vol_1x2v_ser_p1(xc, dx,
     (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidx),
     (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu, pidx), 
     (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu, pidx), 
      qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *grad_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);

  long cidx = gkyl_range_idx(&grad_drag->conf_range, idx);
  int vel_idx[2];
  for (int d=grad_drag->cdim; d<grad_drag->pdim; d++) vel_idx[d-grad_drag->cdim] = idx[d];

  long vidx = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idx);
  long pidx = gkyl_range_idx(&grad_drag->phase_range, idx);

  return rad_gyrokinetic_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidx),
    (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu, pidx), 
    (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu, pidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *grad_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);

  long cidx = gkyl_range_idx(&grad_drag->conf_range, idx);
  int vel_idx[2];
  for (int d=grad_drag->cdim; d<grad_drag->pdim; d++) vel_idx[d-grad_drag->cdim] = idx[d];

  long vidx = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idx);
  long pidx = gkyl_range_idx(&grad_drag->phase_range, idx);

  return rad_gyrokinetic_vol_3x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidx),
    (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu, pidx), 
    (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu, pidx), 
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_rad_gyrokinetic_vol_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, kernel_rad_gyrokinetic_vol_2x2v_ser_p1, NULL }, // 3
  // 3x kernels
  { NULL, kernel_rad_gyrokinetic_vol_3x2v_ser_p1, NULL }, // 4
};

// Surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_surfvpar_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_surfvpar_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: mu-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_surf_kern_list ser_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_surfmu_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_surfmu_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_surfmu_3x2v_ser_p1, NULL }, // 3
};

// Boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_boundary_surfvpar_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_boundary_surfvpar_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_boundary_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Constant nu boundary surface kernel (zero-flux BCs) list: mu-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_boundary_surfmu_1x2v_ser_p1, NULL }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_boundary_surfmu_2x2v_ser_p1, NULL }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_boundary_surfmu_3x2v_ser_p1, NULL }, // 3
};

void gkyl_rad_gyrokinetic_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *grad_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  // Only in vpar,mu directions.
  if (dir >= grad_drag->cdim) {
    int vel_idxL[2], vel_idxC[2], vel_idxR[2];
    for (int d=grad_drag->cdim; d<grad_drag->pdim; d++) {
      vel_idxL[d-grad_drag->cdim] = idxL[d];
      vel_idxC[d-grad_drag->cdim] = idxC[d];
      vel_idxR[d-grad_drag->cdim] = idxR[d];
    }
    long vidxL = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idxL);
    long vidxC = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idxC);
    long vidxR = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idxR);

    // Each cell owns the *lower* edge surface alpha
    // Since alpha is continuous, fetch alpha_surf in center cell for lower edge
    // and fetch alpha_surf in right cell for upper edge
    long pidxC = gkyl_range_idx(&grad_drag->phase_range, idxC);
    long pidxR = gkyl_range_idx(&grad_drag->phase_range, idxR);
    long cidx = gkyl_range_idx(&grad_drag->conf_range, idxR);
    return grad_drag->surf[dir-grad_drag->cdim](xcC, dxC,
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap, vidxC),
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidxL),
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidxC),
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidxR),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu_surf, pidxC),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu_surf, pidxR),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu_surf, pidxC),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu_surf, pidxR),
      qInL, qInC, qInR, qRhsOut);

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
  struct dg_rad_gyrokinetic_drag *grad_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);

  // Only in vpar,mu directions.
  if (dir >= grad_drag->cdim) {
    long cidx = gkyl_range_idx(&grad_drag->conf_range, idxSkin);
    int vel_idxEdge[2], vel_idxSkin[2];
    for (int d=grad_drag->cdim; d<grad_drag->pdim; d++) {
      vel_idxEdge[d-grad_drag->cdim] = idxEdge[d];
      vel_idxSkin[d-grad_drag->cdim] = idxSkin[d];
    }
    long vidxEdge = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idxEdge);
    long vidxSkin = gkyl_range_idx(&grad_drag->vel_map->local_vel, vel_idxSkin);
    
    // Each cell owns the *lower* edge surface alpha
    long pidxEdge = gkyl_range_idx(&grad_drag->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&grad_drag->phase_range, idxSkin);
    
    return grad_drag->boundary_surf[dir-grad_drag->cdim](xcSkin, dxSkin,
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap, vidxSkin),
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidxEdge),
      (const double*) gkyl_array_cfetch(grad_drag->vel_map->vmap_prime, vidxSkin),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu_surf, pidxEdge),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvnu_surf, pidxSkin),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu_surf, pidxEdge),
      (const double*) gkyl_array_cfetch(grad_drag->auxfields.nvsqnu_surf, pidxSkin),
      edge, qInEdge, qInSkin, qRhsOut);
    
  }
  return 0.;
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create a new device-side gyrokinetic RAD drag term equation object.
 *
 * @param conf_basis Configuration-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param phase_range Phase-space range for use in indexing drag coefficients
 * @param conf_range Configuration-space range for use in indexing temperature
 * @param vel_map Velocity space mapping object.
 * @return Pointer to RAD equation object
 */
struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* conf_basis,
  const struct gkyl_basis* phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_velocity_map *vel_map);

/**
 * CUDA device function to set auxiliary fields needed in updating the drag flux term.
 */
void gkyl_rad_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_auxfields auxin);
#endif
