// Private header: not for direct use
#pragma once

#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_pkpm_layout.h>

struct gkyl_pkpm_stage_weights {
  int num_nodes; // nodes per configuration-space cell
  int num_phase_basis; // values per block of a phase-space cell
  dg_conf_nodal_kern_t conf_to_modal; // node values of a configuration field to its expansion
  mul_accumulate_op_t conf_phase_mul_accumulate; // weak product, scaled and accumulated
};
