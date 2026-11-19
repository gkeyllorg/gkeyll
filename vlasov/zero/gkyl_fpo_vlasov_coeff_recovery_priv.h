#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_array.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>

struct gkyl_fpo_vlasov_coeff_recovery {
  int cdim;
  int vdim;
  int pdim;
  int poly_order;

  // Relative offsets in each velocity direction for 3- and 9-cell stencils ( = 36 total offsets)
  long offsets[GKYL_HYPER_DG_GEN_NUM_OFFSETS];

  // Drag coefficient kernel pointers.
  fpo_drag_coeff_t drag_coeff_recovery_stencil[3][3];
  fpo_sgn_drag_coeff_t sgn_drag_coeff_stencil[3][3];

  // Diffusion coefficient kernel pointers.
  fpo_diff_coeff_diag_t diff_coeff_diag_recovery_stencil[3][3];
  fpo_diff_coeff_cross_t diff_coeff_cross_recovery_stencil[3][3][9];
  fpo_diff_coeff_surf_t diff_coeff_surf_recovery[3];

  uint32_t flags;
  struct gkyl_fpo_vlasov_coeff_recovery *on_dev;
};
