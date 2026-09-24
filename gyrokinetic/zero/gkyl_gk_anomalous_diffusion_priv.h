#pragma once

#include <gkyl_gk_anomalous_diffusion.h>
#include <gkyl_gk_anomalous_diffusion_kernels.h>
#include <gkyl_ref_count.h>

// Types for various kernels
typedef double (*gk_anom_diff_surf_t)(const double *wc, const double *dxc,
  const double *nul, const double *nuc, const double *nur,
  const double *jacobgeo_invl, const double *jacobgeo_invc, const double *jacobgeo_invr,
  const double *Jfl, const double *Jfc, const double *Jfr, double* GKYL_RESTRICT out);

typedef double (*gk_anom_diff_boundary_surf_t)(const double *wSkin, const double *dxSkin,
  const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin,
  int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out);

struct gk_anomalous_diffusion {
  struct gkyl_dg_eqn eqn;
  gk_anom_diff_surf_t surf;
  gk_anom_diff_boundary_surf_t boundary_surf[2]; // 2=lower,upper.
  gk_anom_diff_boundary_surf_t boundary_diag[2]; // 2=lower,upper.
  struct gkyl_range conf_range;
  struct gkyl_gk_anomalous_diffusion_auxfields auxfields;
  int num_basis;
};

#define _cfnu(idx) (const double *) gkyl_array_cfetch(gkad->auxfields.nu, gkyl_range_idx(&gkad->conf_range, idx))

#define _cfJacInv(idx) (const double *) gkyl_array_cfetch(gkad->auxfields.jacobgeo_inv, gkyl_range_idx(&gkad->conf_range, idx))

// for use in kernel tables
typedef struct { vol_termf_t kernels[2]; } gkyl_gk_anomalous_diffusion_vol_kern_list;
typedef struct { gk_anom_diff_surf_t kernels[2]; } gkyl_gk_anomalous_diffusion_surf_kern_list;
typedef struct { gk_anom_diff_boundary_surf_t kernels[2]; } gkyl_gk_anomalous_diffusion_boundary_surf_kern_list;

// ............... Inhomogeneous (spatially varying) diffusion coefficient ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 2x2v
GKYL_CU_DH static double ker_gk_anomalous_diffusion_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_anomalous_diffusion* gkad = container_of(eqn, struct gk_anomalous_diffusion, eqn);
  return gk_anomalous_diffusion_vol_2x2v_ser_p1(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_gk_anomalous_diffusion_vol_kern_list ser_vol_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  {ker_gk_anomalous_diffusion_vol_2x2v_ser_p1,NULL}, // 2x2v
  { NULL, NULL }, // 3x2v
};

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_surf_kern_list ser_gyrokinetic_surfx_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_surfx_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};

// Zero-flux boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_lower_zeroflux_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_surfx_lower_zero_flux_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_upper_zeroflux_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_surfx_upper_zero_flux_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};

// Boundary-local boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_lower_boundlocal_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_surfx_lower_bound_local_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_upper_boundlocal_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_surfx_upper_bound_local_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};

// Bound-local boundary diagnostic kernel list: x-direction
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diagx_lower_boundlocal_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_diagx_lower_bound_local_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diagx_upper_boundlocal_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_diagx_upper_bound_local_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};

// Bound-recovery boundary diagnostic kernel list: x-direction
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diagx_lower_boundrecovery_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_diagx_lower_bound_recovery_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};
GKYL_CU_D static const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diagx_upper_boundrecovery_kernels[] = {
  { NULL, NULL }, // 1x1v
  { NULL, NULL }, // 1x2v
  { gk_anomalous_diffusion_boundary_diagx_upper_bound_recovery_2x2v_ser_p1, NULL }, // 2x2v
  { NULL, NULL }, // 3x2v
};

// Macro for choosing volume and surface kernels.
#define CKVOL(lst,pdim,poly_order) lst[pdim-2].kernels[poly_order-1]
#define CKSURF(lst,pdim,poly_order) lst[pdim-2].kernels[poly_order-1]

GKYL_CU_D static double surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct gk_anomalous_diffusion* gkad = container_of(eqn, struct gk_anomalous_diffusion, eqn);

  if (dir == 0) {
    gkad->surf(xcC, dxC, _cfnu(idxL), _cfnu(idxC), _cfnu(idxR), _cfJacInv(idxL), _cfJacInv(idxC), _cfJacInv(idxR), qInL, qInC, qInR, qRhsOut);
  }
  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin, const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  struct gk_anomalous_diffusion* gkad = container_of(eqn, struct gk_anomalous_diffusion, eqn);
  
  if (dir == 0) {
    if (edge == -1) {
      gkad->boundary_surf[0](xcSkin, dxSkin, _cfnu(idxEdge), _cfnu(idxSkin), _cfJacInv(idxEdge), _cfJacInv(idxSkin), edge, qInEdge, qInSkin, qRhsOut);
    }
    else {
      gkad->boundary_surf[1](xcSkin, dxSkin, _cfnu(idxEdge), _cfnu(idxSkin), _cfJacInv(idxEdge), _cfJacInv(idxSkin), edge, qInEdge, qInSkin, qRhsOut);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_diag(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcSkin, const double* xcGhost, const double* dxSkin, const double* dxGhost,
  const int* idxSkin, const int* idxGhost, const int edge,
  const double* qInSkin, const double* qInGhost, double* GKYL_RESTRICT qRhsGhost)
{ 
  // This function is based on boundary_surf above, but notice we use Skin
  // where the boundary_surf used Ghost, because we assume this kernel is called
  // in the ghost range (e.g. by the boundary_flux updater).
  struct gk_anomalous_diffusion* gkad = container_of(eqn, struct gk_anomalous_diffusion, eqn);

  if (dir == 0) {
    if (edge == -1) {
      gkad->boundary_diag[0](xcSkin, dxSkin, _cfnu(idxSkin), _cfnu(idxGhost), _cfJacInv(idxSkin), _cfJacInv(idxGhost), edge, qInSkin, qInGhost, qRhsGhost);
    }
    else {
      gkad->boundary_diag[1](xcSkin, dxSkin, _cfnu(idxSkin), _cfnu(idxGhost), _cfJacInv(idxSkin), _cfJacInv(idxGhost), edge, qInSkin, qInGhost, qRhsGhost);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

#undef _cfnu
#undef _cfJacInv

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_gk_anomalous_diffusion_free(const struct gkyl_ref_count* ref);

#ifdef GKYL_HAVE_CUDA
/**
 * Create a new gyrokinetic diffusion equation object on the device.
 *
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param conf_range Conf-space range object.
 * @param bc_x_lower Boundary condition at lower x boundary.
 * @param bc_x_upper Boundary condition at upper x boundary.
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn*
gkyl_gk_anomalous_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *conf_range, enum gkyl_gyrokinetic_bc_type bc_x_lower, enum gkyl_gyrokinetic_bc_type bc_x_upper);

/**
 * CUDA device function to set auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gk_anomalous_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_gk_anomalous_diffusion_auxfields auxin);

#endif
