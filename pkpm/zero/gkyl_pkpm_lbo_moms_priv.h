// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_dg_conf_nodal_kernels.h>
#include <gkyl_range.h>

// Largest number of nodes per configuration cell (3x).
#define PKPM_LBO_MOMS_MAX_NODES 8

typedef void (*pkpm_lbo_moms_conf_t)(const double *in, double *GKYL_RESTRICT out);

GKYL_CU_D static const pkpm_lbo_moms_conf_t pkpm_lbo_moms_conf_to_nodal_list[] = {
  dg_conf_nodal_conf_to_nodal_1x_ser_p1, dg_conf_nodal_conf_to_nodal_2x_ser_p1,
  dg_conf_nodal_conf_to_nodal_3x_ser_p1
};

struct gkyl_pkpm_lbo_moms {
  int cdim; // configuration-space dimension
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double vel_bounds[2]; // lower and upper bound of the velocity grid
  pkpm_lbo_moms_conf_t conf_to_nodal; // modal-to-nodal conversion of a configuration field
};
