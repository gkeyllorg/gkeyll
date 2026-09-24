#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Types for various kernels
typedef double (*canonical_pb_fluid_surf_t)(const double *w, const double *dxv,
  const double *phi, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_canonical_pb_fluid_vol_kern_list;
typedef struct { canonical_pb_fluid_surf_t kernels[3]; } gkyl_dg_canonical_pb_fluid_surf_kern_list;

struct dg_canonical_pb_fluid {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  canonical_pb_fluid_surf_t surf[3]; // Surface terms for streaming
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_canonical_pb_fluid_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_canonical_pb_fluid_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  long cidx = gkyl_range_idx(&can_pb_fluid->conf_range, idx);

  return canonical_pb_vol_2x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.phi, cidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_fluid_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  long cidx = gkyl_range_idx(&can_pb_fluid->conf_range, idx);

  return canonical_pb_vol_2x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.phi, cidx), 
    qIn, qRhsOut);  
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, kernel_canonical_pb_fluid_vol_2x_ser_p1, kernel_canonical_pb_fluid_vol_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_vol_kern_list tensor_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, kernel_canonical_pb_fluid_vol_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

//
// Serendipity volume kernels for two-fluid canonical PB fluid system
// such as Hasegawa-Wakatani
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_canonical_pb_two_fluid_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  long cidx = gkyl_range_idx(&can_pb_fluid->conf_range, idx);

  return canonical_pb_two_fluid_vol_2x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.phi, cidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_two_fluid_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  long cidx = gkyl_range_idx(&can_pb_fluid->conf_range, idx);

  return canonical_pb_two_fluid_vol_2x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.phi, cidx), 
    qIn, qRhsOut);  
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_vol_kern_list ser_two_fluid_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, kernel_canonical_pb_two_fluid_vol_2x_ser_p1, kernel_canonical_pb_two_fluid_vol_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_vol_kern_list tensor_two_fluid_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, kernel_canonical_pb_two_fluid_vol_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list ser_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfx_2x_ser_p1, canonical_pb_surfx_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list tensor_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfx_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list ser_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfy_2x_ser_p1, canonical_pb_surfy_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list tensor_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfy_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface two fluid kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list ser_two_fluid_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_two_fluid_surfx_2x_ser_p1, canonical_pb_two_fluid_surfx_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface two fluid kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list tensor_two_fluid_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_two_fluid_surfx_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface two fluid kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list ser_two_fluid_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_two_fluid_surfy_2x_ser_p1, canonical_pb_two_fluid_surfy_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Surface two fluid kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_surf_kern_list tensor_two_fluid_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_two_fluid_surfy_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

/**
 * Free canonical PB for fluids eqn object.
 *
 * @param ref Reference counter for canonical PB for fluids eqn
 */
void gkyl_canonical_pb_fluid_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{    
  // Each cell owns the *lower* edge surface alpha
  // Since alpha is continuous, fetch alpha_surf in center cell for lower edge
  // and fetch alpha_surf in right cell for upper edge
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  long cidxC = gkyl_range_idx(&can_pb_fluid->conf_range, idxC);
  long cidxR = gkyl_range_idx(&can_pb_fluid->conf_range, idxR);
  if (dir < can_pb_fluid->cdim) {
    return can_pb_fluid->surf[dir](xcC, dxC,
        (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.phi, cidxC),
        (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.alpha_surf, cidxC), 
        (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.alpha_surf, cidxR), 
        (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.sgn_alpha_surf, cidxC), 
        (const double*) gkyl_array_cfetch(can_pb_fluid->auxfields.sgn_alpha_surf, cidxR), 
        (const int*) gkyl_array_cfetch(can_pb_fluid->auxfields.const_sgn_alpha, cidxC), 
        (const int*) gkyl_array_cfetch(can_pb_fluid->auxfields.const_sgn_alpha, cidxR), 
        qInL, qInC, qInR, qRhsOut);    
  }
  return 0.;
}
