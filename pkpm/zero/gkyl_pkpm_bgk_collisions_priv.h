// Private header: not for direct use
#pragma once

#define PKPM_BGK_NUM_BLOCKS 8
#define PKPM_BGK_NUM_MODES 3
#define PKPM_BGK_MAX_NODES 8

struct gkyl_pkpm_bgk_collisions {
  int cdim; // configuration-space dimensions
  int num_nodes; // nodes per configuration-space cell
  int num_phase_basis; // values per block of a phase-space cell
};
