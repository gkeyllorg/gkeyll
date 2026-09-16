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

typedef void (*calc_hamil_t)(const double *w, const double *dv, 
  const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv);

// for use in kernel tables
typedef struct { calc_hamil_t kernels[4]; } gkyl_dg_calc_hamil_kern_list;

// Default Hamiltonian = v^2/2 (Serendipity kernels).
GKYL_CU_D
static const gkyl_dg_calc_hamil_kern_list ser_hamil_default_kernels[] = {
  // 1x kernels
  { NULL, hamil_default_1v_ser_p1, hamil_default_1v_ser_p2, hamil_default_1v_ser_p3 }, // 0
  { NULL, hamil_default_2v_ser_p1, hamil_default_2v_ser_p2, NULL }, // 1
  { NULL, hamil_default_3v_ser_p1, hamil_default_3v_ser_p2, NULL }, // 2
};

// Relativistic Hamiltonian = sqrt(1 + u^2) (Serendipity kernels).
GKYL_CU_D
static const gkyl_dg_calc_hamil_kern_list ser_hamil_sr_kernels[] = {
  // 1x kernels
  { NULL, hamil_sr_1v_ser_p1, hamil_sr_1v_ser_p2, hamil_sr_1v_ser_p3 }, // 0
  { NULL, hamil_sr_2v_ser_p1, hamil_sr_2v_ser_p2, NULL }, // 1
  { NULL, hamil_sr_3v_ser_p1, hamil_sr_3v_ser_p2, NULL }, // 2
};

// Default Hamiltonian = v^2/2 (Tensor kernels).
GKYL_CU_D
static const gkyl_dg_calc_hamil_kern_list tensor_hamil_default_kernels[] = {
  // 1x kernels
  { NULL, hamil_default_1v_tensor_p1, hamil_default_1v_tensor_p2, hamil_default_1v_tensor_p3 }, // 0
  { NULL, hamil_default_2v_tensor_p1, hamil_default_2v_tensor_p2, NULL }, // 1
  { NULL, hamil_default_3v_tensor_p1, hamil_default_3v_tensor_p2, NULL }, // 2
};

// Relativistic Hamiltonian = sqrt(1 + u^2) (Tensor kernels).
GKYL_CU_D
static const gkyl_dg_calc_hamil_kern_list tensor_hamil_sr_kernels[] = {
  // 1x kernels
  { NULL, hamil_sr_1v_tensor_p1, hamil_sr_1v_tensor_p2, hamil_sr_1v_tensor_p3 }, // 0
  { NULL, hamil_sr_2v_tensor_p1, hamil_sr_2v_tensor_p2, NULL }, // 1
  { NULL, hamil_sr_3v_tensor_p1, hamil_sr_3v_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static calc_hamil_t
choose_ser_hamil_kern(enum gkyl_model_id model_id, int vdim, int poly_order)
{
  if (model_id == GKYL_MODEL_DEFAULT) {
    return ser_hamil_default_kernels[vdim-1].kernels[poly_order];
  }
  else if (model_id == GKYL_MODEL_SR) {
    return ser_hamil_sr_kernels[vdim-1].kernels[poly_order];
  }
  else {
    assert(false); 
    return 0;   
  }
  return 0;
}

GKYL_CU_D
static calc_hamil_t
choose_tensor_hamil_kern(enum gkyl_model_id model_id, int vdim, int poly_order)
{
  if (model_id == GKYL_MODEL_DEFAULT) {
    return tensor_hamil_default_kernels[vdim-1].kernels[poly_order];
  }
  else if (model_id == GKYL_MODEL_SR) {
    return tensor_hamil_sr_kernels[vdim-1].kernels[poly_order];
  }
  else {
    assert(false); 
    return 0;   
  }
  return 0;
}