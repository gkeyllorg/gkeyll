// Private header: not for direct use
#pragma once

#include <gkyl_dg_vlasov_pressure.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_range.h>

struct gkyl_dg_vlasov_pressure
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_range conf_range; // Configuration-space range
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim; // Number of velocity dimensions
  double mass; 
  bool use_gpu;

  struct gkyl_array *V_drift_dot_M1i;
  struct gkyl_array *V_drift_M1i;
  struct gkyl_array *V_drift_M1i_transpose;
};
