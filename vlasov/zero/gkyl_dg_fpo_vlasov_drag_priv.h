#pragma once

// Private header, not for direct use in user code
#include <gkyl_fpo_vlasov_kernels.h>
#include <string.h>

// Types for various kernels
typedef double (*fpo_vlasov_drag_surf_t)(const double* dxv, 
  const double *alpha_surf_L, const double *alpha_surf_R,
  const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R,
  const int *const_sgn_alpha_L, const int *const_sgn_alpha_R,
  const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

typedef double (*fpo_vlasov_drag_boundary_surf_t)(const double* dxv,
  const double *alpha_surf_Edge, const double *alpha_surf_Skin,
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin,
  const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin,
  const int Edge, const double *fEdge, const double *fSkin,
  double* GKYL_RESTRICT out); 

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_vol_kern_list;
typedef struct { fpo_vlasov_drag_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_surf_kern_list;
typedef struct { fpo_vlasov_drag_boundary_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list;

struct dg_fpo_vlasov_drag {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  fpo_vlasov_drag_surf_t surf[3]; // Surface terms
  fpo_vlasov_drag_boundary_surf_t boundary_surf[3];
  struct gkyl_range phase_range; // Configuration space range.
  struct gkyl_dg_fpo_vlasov_drag_auxfields auxfields; // Auxiliary fields.
};

// Serendipity volume kernels
GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  return fpo_vlasov_drag_vol_1x3v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  return fpo_vlasov_drag_vol_1x3v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  return fpo_vlasov_drag_vol_2x3v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff, pidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_fpo_vlasov_drag_vol_1x3v_ser_p1, kernel_fpo_vlasov_drag_vol_1x3v_ser_p2 },
  // 2x kernels
  { NULL, kernel_fpo_vlasov_drag_vol_2x3v_ser_p1, NULL },
  // 3x kernels
  { NULL, NULL, NULL },
};

// Surface kernel lists
// vx-direction surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, fpo_vlasov_drag_surfvx_1x3v_ser_p1, fpo_vlasov_drag_surfvx_1x3v_ser_p2 },
  // 2x kernels
  { NULL, fpo_vlasov_drag_surfvx_2x3v_ser_p1, NULL},
  // 3x kernels
  { NULL, NULL, NULL },
};

// vy-direction surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, fpo_vlasov_drag_surfvy_1x3v_ser_p1, fpo_vlasov_drag_surfvy_1x3v_ser_p2 },
  // 2x kernels
  { NULL, fpo_vlasov_drag_surfvy_2x3v_ser_p1, NULL},
  // 3x kernels
  { NULL, NULL, NULL },
};

// vz-direction surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, fpo_vlasov_drag_surfvz_1x3v_ser_p1, fpo_vlasov_drag_surfvz_1x3v_ser_p2 },
  // 2x kernels
  { NULL, fpo_vlasov_drag_surfvz_2x3v_ser_p1, NULL},
  // 3x kernels
  { NULL, NULL, NULL },
};

// vx-direction boundary surf kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p2 },
  // 2x kernels
  { NULL, fpo_vlasov_drag_boundary_surfvx_2x3v_ser_p1, NULL},
  // 3x kernels
  { NULL, NULL, NULL },
};

// vy-direction boundary surf kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p2 },
  // 2x kernels
  { NULL, fpo_vlasov_drag_boundary_surfvy_2x3v_ser_p1, NULL},
  // 3x kernels
  { NULL, NULL, NULL },
};

// vz-direction boundary surf kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p2 },
  // 2x kernels
  { NULL, fpo_vlasov_drag_boundary_surfvz_2x3v_ser_p1, NULL},
  // 3x kernels
  { NULL, NULL, NULL },
};

/* Free fpo_vlasov_diff equation object
 *
 * @param ref Reference counter for constant fpo_vlasov_diff equation
*/
void gkyl_fpo_vlasov_drag_free(const struct gkyl_ref_count* ref);

// Surface term called by gkyl_hyper_dg_advance
GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);

  long linl = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxL);
  long linc = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxC);
  long linr = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxR);

  if (dir >= fpo_vlasov_drag->cdim) {
    return fpo_vlasov_drag->surf[dir-fpo_vlasov_drag->cdim](dxC,
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff_surf, linc),
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff_surf, linr),
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf, linc),
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf, linr),
      (const int*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf, linc),
      (const int*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf, linr),
      qInL, qInC, qInR, qRhsOut);
  }
  return 0.0;
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
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);

  long lin_edge = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxEdge);
  long lin_skin = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idxSkin);

  if (dir >= fpo_vlasov_drag->cdim) {
    return fpo_vlasov_drag->boundary_surf[dir-fpo_vlasov_drag->cdim](dxEdge,
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff_surf, lin_edge),
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff_surf, lin_skin),
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf, lin_edge),
      (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf, lin_skin),
      (const int*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf, lin_edge),
      (const int*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf, lin_skin),
      edge, qInEdge, qInSkin, qRhsOut);
  }

  return 0.0;
}

