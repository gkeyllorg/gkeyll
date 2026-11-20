// Private header: not for direct use
#pragma once

#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_pkpm_layout.h>

struct gkyl_pkpm_bgk_collisions {
  int cdim; // configuration-space dimensions
  int num_nodes; // nodes per configuration-space cell
  int num_phase_basis; // values per block of a phase-space cell
  dg_conf_nodal_kern_t conf_to_modal; // node values of a configuration field to its expansion
  dg_conf_nodal_kern_t phase_to_modal; // nodal form of a block to its expansion (accumulated)
  mul_op_t conf_phase_mul; // weak product of a configuration field and a block
  mul_accumulate_op_t conf_phase_mul_accumulate; // the same, scaled and accumulated
};
