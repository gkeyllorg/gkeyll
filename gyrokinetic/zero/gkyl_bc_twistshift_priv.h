#pragma once

// Private header for the bc_twistshift orchestrator, not for direct use in
// user code.

#include <gkyl_bc_twistshift.h>
#include <gkyl_twistshift_dg.h>
#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>

// Primary struct in this updater.
struct gkyl_bc_twistshift {
  bool use_gpu; // Whether to apply the BC on the GPU.

  struct gkyl_twistshift_dg *ts; // Pure twist-shift updater.

  void (*advance_func)(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar);

  int filter_half_width; // Filter stencil half-width M in cells of the field's
                         // own grid (0 = no filter).
  int half_width_fine; // The same stencil measured in supersampled cells.
  double filter_cutoff_wavelength; // Filter cutoff wavelength.
  int upsample_factor; // Supersampling factor along shear_dir.

  struct gkyl_dg_lowpass_filter *filter; // Post-shift filter along shear_dir.
  struct gkyl_rect_grid ts_grid; // Grid refined along shear_dir.
  struct gkyl_range ts_ext_r; // ffine and filt_buff range.
  struct gkyl_range ts_update_r; // Update range on ts_grid.
  struct gkyl_range ghost_r; // Ghost plane the twist-shift fills, on ts_grid.
  struct gkyl_range coarse_ghost_r; // Same plane on the field's own grid.
  struct gkyl_array *ffine; // Ghost plane on the refined grid.
  struct gkyl_array *filt_buff; // Buffer for the filter (shaped like ffine).
  struct gkyl_array *shift_dg_fine; // Input shift refined onto the fine shear grid.

  struct gkyl_dg_interpolate *refine; // Coarse ghost plane -> fine ghost plane.
  struct gkyl_dg_interpolate *coarsen; // Fine ghost plane -> coarse ghost plane.
  void (*refine_func)(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo);
  void (*coarsen_func)(struct gkyl_bc_twistshift *up, struct gkyl_array *ftar);
};
