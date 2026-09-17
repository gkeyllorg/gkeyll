#pragma once

#include <gkyl_array.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>

struct gkyl_fpo_vlasov_coeff_recovery {
  int cdim;
  int vdim;
  int pdim;
  int poly_order;

  // Relative offsets in each velocity direction for 3- and 9-cell stencils ( = 36 total offsets)
  long offsets[36];

  // Drag coefficient kernel pointers.
  fpo_drag_coeff_t drag_coeff_recovery_stencil[3][3];
  fpo_sgn_drag_coeff_t sgn_drag_coeff_stencil[3][3];

  // Diffusion coefficient kernel pointers.
  fpo_diff_coeff_diag_t diff_coeff_diag_recovery_stencil[3][3];
  fpo_diff_coeff_cross_t diff_coeff_cross_recovery_stencil[3][3][9];
  fpo_diff_coeff_surf_t diff_coeff_surf_recovery[3];

  uint32_t flags;
  struct gkyl_fpo_vlasov_coeff_recovery* on_dev;
};

static void
create_offsets(const int num_up_dirs, const int update_dirs[2], 
  const struct gkyl_range *range, const int idxc[GKYL_MAX_DIM], long offsets[9])
{
  
  // Check if we're at an upper or lower edge in each direction
  bool is_edge_upper[2], is_edge_lower[2];
  for (int i=0; i<num_up_dirs; ++i) {
    is_edge_lower[i] = idxc[update_dirs[i]] == range->lower[update_dirs[i]];
    is_edge_upper[i] = idxc[update_dirs[i]] == range->upper[update_dirs[i]];
  }

  // Construct the offsets *only* in the directions being updated.
  // No need to load the neighbors that are not needed for the update.
  int lower_offset[GKYL_MAX_DIM] = {0};
  int upper_offset[GKYL_MAX_DIM] = {0};
  for (int d=0; d<num_up_dirs; ++d) {
    int dir = update_dirs[d];
    lower_offset[dir] = -1 + is_edge_lower[d];
    upper_offset[dir] = 1 - is_edge_upper[d];
  }  

  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, lower_offset, upper_offset);
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}
