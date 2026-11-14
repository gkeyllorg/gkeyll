// Private header: not for direct use
#pragma once

#include <gkyl_dg_conf_nodal_kernels.h>
#include <gkyl_pkpm_layout.h>

typedef void (*pkpm_temp_source_conf_map_t)(const double *in, double *GKYL_RESTRICT out);

GKYL_CU_D static const pkpm_temp_source_conf_map_t pkpm_temp_source_conf_to_nodal_list[] = {
  dg_conf_nodal_conf_to_nodal_1x_ser_p1, dg_conf_nodal_conf_to_nodal_2x_ser_p1,
  dg_conf_nodal_conf_to_nodal_3x_ser_p1
};

GKYL_CU_D static const pkpm_temp_source_conf_map_t pkpm_temp_source_nodal_to_conf_list[] = {
  dg_conf_nodal_conf_to_modal_1x_ser_p1, dg_conf_nodal_conf_to_modal_2x_ser_p1,
  dg_conf_nodal_conf_to_modal_3x_ser_p1
};

struct gkyl_pkpm_temp_source {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  pkpm_temp_source_conf_map_t conf_to_nodal; // modal to nodal, for this cdim
  pkpm_temp_source_conf_map_t nodal_to_conf; // nodal to modal, for this cdim
};
