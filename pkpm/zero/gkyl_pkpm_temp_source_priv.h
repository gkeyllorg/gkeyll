// Private header: not for direct use
#pragma once

#include <gkyl_pkpm_coll_nodal.h>

#define PKPM_TEMP_SOURCE_MAX_NODES 8
#define PKPM_TEMP_SOURCE_NUM_FLUID 5 // blocks of the fluid state [rho, M (3), E]
#define PKPM_TEMP_SOURCE_ENERGY 4 // block of the energy

typedef void (*pkpm_temp_source_conf_map_t)(const double *in, double *GKYL_RESTRICT out);

GKYL_CU_D static const pkpm_temp_source_conf_map_t pkpm_temp_source_conf_to_nodal_list[] = {
  pkpm_coll_conf_to_nodal_1x, pkpm_coll_conf_to_nodal_2x, pkpm_coll_conf_to_nodal_3x
};

GKYL_CU_D static const pkpm_temp_source_conf_map_t pkpm_temp_source_nodal_to_conf_list[] = {
  pkpm_coll_nodal_to_conf_1x, pkpm_coll_nodal_to_conf_2x, pkpm_coll_nodal_to_conf_3x
};

struct gkyl_pkpm_temp_source {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  pkpm_temp_source_conf_map_t conf_to_nodal; // modal to nodal, for this cdim
  pkpm_temp_source_conf_map_t nodal_to_conf; // nodal to modal, for this cdim
};
