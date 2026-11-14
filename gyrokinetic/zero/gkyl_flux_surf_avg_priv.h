#pragma once

// Private header for the flux_surf_avg updater, not for direct use in user code.

#include <gkyl_array_average.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_flux_surf_avg.h>
#include <gkyl_translate_dim.h>

struct gkyl_flux_surf_avg {
  struct gkyl_basis basis; // Configuration-space basis.
  struct gkyl_basis basis_x; // 1D basis in x.
  struct gkyl_range local; // Local configuration-space range.
  struct gkyl_range local_x; // 1D range in x.
  struct gkyl_comm *comm; // Communicator for the reduction over the averaged directions (or NULL).
  bool use_gpu; // Whether to run on the GPU.

  struct gkyl_array *jacobgeo; // Jacobian J weighting the average.
  struct gkyl_array_average *avg; // (y,z) average.
  struct gkyl_array *jf; // Buffer for J*f.
  struct gkyl_array *avg_jf, *avg_jf_red; // (y,z) average of J*f, local and reduced over ranks.
  struct gkyl_array *avg_j; // Global (y,z) average of J.
  struct gkyl_array *avg_buff; // 1D buffer for the average before its extension (if !is_deflated).
  gkyl_dg_bin_op_mem *div_mem; // Memory for the 1D weak division.

  struct gkyl_range local_xy, local_xy_ext; // 2D (x,y) ranges (3x only).
  struct gkyl_array *inflate_buff; // 2D buffer for the 1D->3D extension (3x only).
  struct gkyl_translate_dim *inflate_lo, *inflate_up; // 1D->2D and 2D->3D extension updaters.
  // Sums the local (y,z) averages over the ranks of comm (a copy if there is no comm).
  void (*reduce_func)(
    struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg_jf, struct gkyl_array *avg_jf_red
  );
  // Computes the average into a 1D field (is_deflated) or into a cdim field (!is_deflated).
  void (*advance_func)(
    struct gkyl_flux_surf_avg *up, const struct gkyl_array *fin, struct gkyl_array *out
  );
  // Extends a function of x to a cdim field.
  void (*inflate_func)(
    struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg, struct gkyl_array *out
  );
};
