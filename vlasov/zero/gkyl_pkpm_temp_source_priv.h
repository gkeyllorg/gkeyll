// Private header: not for direct use
#pragma once

#include <stdbool.h>

#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_pkpm_layout.h>

struct gkyl_pkpm_temp_source {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double energy_fac; // thermal energy per unit of rho vtsq, 1/(gas_gamma - 1)
  dg_conf_nodal_kern_t conf_to_nodal; // expansion to node values
  dg_conf_nodal_kern_t conf_to_modal; // node values to expansion
};
