// Private header: not for direct use
#pragma once

#include <gkyl_binop_cross_mul_hyb.h>
#include <gkyl_dg_conf_nodal_kernels.h>

#define PKPM_STAGE_WEIGHTS_MAX_NODES 8

// Node values of a configuration field to its expansion.
typedef void (*pkpm_stage_weights_conv_t)(const double *in, double *GKYL_RESTRICT out);
// Weak product of a configuration-space expansion and one block of a phase-space cell, times a
// factor, accumulated.
typedef void (*pkpm_stage_weights_mul_t)(
  double fac, const double *conf, const double *phase, double *out
);

GKYL_CU_D static const pkpm_stage_weights_conv_t pkpm_stage_weights_conf_to_modal_list[] = {
  dg_conf_nodal_conf_to_modal_1x_ser_p1, dg_conf_nodal_conf_to_modal_2x_ser_p1,
  dg_conf_nodal_conf_to_modal_3x_ser_p1
};

GKYL_CU_D static const pkpm_stage_weights_mul_t pkpm_stage_weights_mul_accumulate_list[] = {
  binop_cross_mul_accumulate_1x1v_hyb_p1, binop_cross_mul_accumulate_2x1v_hyb_p1,
  binop_cross_mul_accumulate_3x1v_hyb_p1
};

struct gkyl_pkpm_stage_weights {
  int num_nodes; // nodes per configuration-space cell
  int num_phase_basis; // values per block of a phase-space cell
  pkpm_stage_weights_conv_t conf_to_modal; // node values of a configuration field to its expansion
  pkpm_stage_weights_mul_t conf_phase_mul_accumulate; // weak product, scaled and accumulated
};
