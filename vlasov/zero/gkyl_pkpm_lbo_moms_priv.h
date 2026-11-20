// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_range.h>

struct gkyl_pkpm_lbo_moms {
  int cdim; // configuration-space dimension
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double vel_bounds[2]; // lower and upper bound of the velocity grid
  double edge_modes[2][PKPM_NUM_VEL_MODES]; // velocity modes at the lower and upper edge of a cell
  dg_conf_nodal_kern_t conf_to_nodal; // expansion to node values
};
