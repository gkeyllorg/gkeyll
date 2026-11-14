// Private header: not for direct use
#pragma once

#include <gkyl_binop_cross_mul_hyb.h>
#include <gkyl_dg_conf_nodal_kernels.h>
#include <gkyl_pkpm_layout.h>

#define PKPM_BGK_MAX_PHASE_BASIS (PKPM_NUM_VEL_MODES * PKPM_MAX_CONF_BASIS)

// A conversion between node values and an expansion (configuration field, or one block of a
// phase-space cell).
typedef void (*pkpm_bgk_conv_t)(const double *in, double *GKYL_RESTRICT out);
// Weak product of a configuration-space expansion and one block of a phase-space cell.
typedef void (*pkpm_bgk_mul_t)(const double *conf, const double *phase, double *out);
// The same product times a factor, accumulated.
typedef void (*pkpm_bgk_mul_accumulate_t)(
  double fac, const double *conf, const double *phase, double *out
);

GKYL_CU_D static const pkpm_bgk_conv_t pkpm_bgk_conf_to_modal_list[] = {
  dg_conf_nodal_conf_to_modal_1x_ser_p1, dg_conf_nodal_conf_to_modal_2x_ser_p1,
  dg_conf_nodal_conf_to_modal_3x_ser_p1
};

GKYL_CU_D static const pkpm_bgk_conv_t pkpm_bgk_phase_to_modal_list[] = {
  dg_conf_nodal_phase_to_modal_accumulate_1x1v_hyb_p1,
  dg_conf_nodal_phase_to_modal_accumulate_2x1v_hyb_p1,
  dg_conf_nodal_phase_to_modal_accumulate_3x1v_hyb_p1
};

GKYL_CU_D static const pkpm_bgk_mul_t pkpm_bgk_conf_phase_mul_list[] = {
  binop_cross_mul_1x1v_hyb_p1, binop_cross_mul_2x1v_hyb_p1, binop_cross_mul_3x1v_hyb_p1
};

GKYL_CU_D static const pkpm_bgk_mul_accumulate_t pkpm_bgk_conf_phase_mul_accumulate_list[] = {
  binop_cross_mul_accumulate_1x1v_hyb_p1, binop_cross_mul_accumulate_2x1v_hyb_p1,
  binop_cross_mul_accumulate_3x1v_hyb_p1
};

struct gkyl_pkpm_bgk_collisions {
  int cdim; // configuration-space dimensions
  int num_nodes; // nodes per configuration-space cell
  int num_phase_basis; // values per block of a phase-space cell
  pkpm_bgk_conv_t conf_to_modal; // node values of a configuration field to its expansion
  pkpm_bgk_conv_t phase_to_modal; // nodal form of a block to its expansion (accumulated)
  pkpm_bgk_mul_t conf_phase_mul; // weak product of a configuration field and a block
  pkpm_bgk_mul_accumulate_t conf_phase_mul_accumulate; // the same, scaled and accumulated
};
