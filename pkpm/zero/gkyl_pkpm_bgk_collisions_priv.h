// Private header: not for direct use
#pragma once

#include <gkyl_pkpm_coll_nodal.h>

#define PKPM_BGK_NUM_BLOCKS 8
#define PKPM_BGK_NUM_MODES 3
#define PKPM_BGK_MAX_NODES 8

typedef void (*pkpm_bgk_phase_map_t)(const double *in, double *GKYL_RESTRICT out);

GKYL_CU_D static const pkpm_bgk_phase_map_t pkpm_bgk_phase_to_nodal_list[] = {
  pkpm_coll_phase_to_nodal_1x1v, pkpm_coll_phase_to_nodal_2x1v, pkpm_coll_phase_to_nodal_3x1v
};

GKYL_CU_D static const pkpm_bgk_phase_map_t pkpm_bgk_nodal_to_phase_list[] = {
  pkpm_coll_nodal_to_phase_accumulate_1x1v, pkpm_coll_nodal_to_phase_accumulate_2x1v,
  pkpm_coll_nodal_to_phase_accumulate_3x1v
};

struct gkyl_pkpm_bgk_collisions {
  int cdim; // configuration-space dimensions
  int num_nodes; // nodes per configuration-space cell
  pkpm_bgk_phase_map_t phase_to_nodal; // modal to nodal, for this cdim
  pkpm_bgk_phase_map_t nodal_to_phase; // accumulate a nodal increment on modal, for this cdim
};
