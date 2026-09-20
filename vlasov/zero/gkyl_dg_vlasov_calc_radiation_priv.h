// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*calc_radiation_t)(const double *w, const double *dv, 
  const double *vmap, double t_cool, double p0, double* GKYL_RESTRICT rad);

// for use in kernel tables
typedef struct { calc_radiation_t kernels[4]; } gkyl_dg_calc_radiation_kern_list;

// Compton radiation, rad = -1/t_cool*(p - p0) (Serendipity kernels).
GKYL_CU_D
static const gkyl_dg_calc_radiation_kern_list ser_rad_compton_kernels[] = {
  // 1x kernels
  { NULL, rad_compton_1v_ser_p1, rad_compton_1v_ser_p2, rad_compton_1v_ser_p3 }, // 0
  { NULL, rad_compton_2v_ser_p1, rad_compton_2v_ser_p2, NULL }, // 1
  { NULL, rad_compton_3v_ser_p1, rad_compton_3v_ser_p2, NULL }, // 2
};

// Curvature radiation, rad = -1/t_cool*sign(p)*p^4 = -1/t_cool*fabs(p)*p^3 (Serendipity kernels).
GKYL_CU_D
static const gkyl_dg_calc_radiation_kern_list ser_rad_curvature_kernels[] = {
  // 1x kernels
  { NULL, rad_curvature_1v_ser_p1, rad_curvature_1v_ser_p2, rad_curvature_1v_ser_p3 }, // 0
  { NULL, rad_curvature_2v_ser_p1, rad_curvature_2v_ser_p2, NULL }, // 1
  { NULL, rad_curvature_3v_ser_p1, rad_curvature_3v_ser_p2, NULL }, // 2
};

// Compton radiation, rad = -1/t_cool*(p - p0) (Tensor kernels).
GKYL_CU_D
static const gkyl_dg_calc_radiation_kern_list tensor_rad_compton_kernels[] = {
  // 1x kernels
  { NULL, rad_compton_1v_tensor_p1, rad_compton_1v_tensor_p2, rad_compton_1v_tensor_p3 }, // 0
  { NULL, rad_compton_2v_tensor_p1, rad_compton_2v_tensor_p2, NULL }, // 1
  { NULL, rad_compton_3v_tensor_p1, rad_compton_3v_tensor_p2, NULL }, // 2
};

// Curvature radiation, rad = -1/t_cool*sign(p)*p^4 = -1/t_cool*fabs(p)*p^3 (Tensor kernels).
GKYL_CU_D
static const gkyl_dg_calc_radiation_kern_list tensor_rad_curvature_kernels[] = {
  // 1x kernels
  { NULL, rad_curvature_1v_tensor_p1, rad_curvature_1v_tensor_p2, rad_curvature_1v_tensor_p3 }, // 0
  { NULL, rad_curvature_2v_tensor_p1, rad_curvature_2v_tensor_p2, NULL }, // 1
  { NULL, rad_curvature_3v_tensor_p1, rad_curvature_3v_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static calc_radiation_t
choose_ser_radiation_kern(enum gkyl_vlasov_radiation_id radiation_id, int vdim, int poly_order)
{
  if (radiation_id == GKYL_VM_COMPTON_RADIATION) {
    return ser_rad_compton_kernels[vdim-1].kernels[poly_order];
  }
  else if (radiation_id == GKYL_VM_CURVATURE_RADIATION) {
    return ser_rad_curvature_kernels[vdim-1].kernels[poly_order];
  }
  else {
    assert(false); 
    return 0;   
  }
  return 0;
}

GKYL_CU_D
static calc_radiation_t
choose_tensor_radiation_kern(enum gkyl_vlasov_radiation_id radiation_id, int vdim, int poly_order)
{
  if (radiation_id == GKYL_VM_COMPTON_RADIATION) {
    return tensor_rad_compton_kernels[vdim-1].kernels[poly_order];
  }
  else if (radiation_id == GKYL_VM_CURVATURE_RADIATION) {
    return tensor_rad_curvature_kernels[vdim-1].kernels[poly_order];
  }
  else {
    assert(false); 
    return 0;   
  }
  return 0;
}