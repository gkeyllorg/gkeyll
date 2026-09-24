#pragma once

#include <gkyl_mom_type.h>
#include <gkyl_rect_grid.h>

struct gkyl_mom_calc_bcorr {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_basis; // number of basis functions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  int space;
  const struct gkyl_mom_type *momt; // moment type object

  uint32_t flags;
  struct gkyl_mom_calc_bcorr *on_dev;
};

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to update boundary corrections on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_calc_bcorr* 
gkyl_mom_calc_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Compute boundary correction moments.
 *
 * @param bcorr Boundary correction updater object
 * @param phase_rng Phase space range on which to compute.
 * @param conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param out Output
 */
void gkyl_mom_calc_bcorr_advance_cu(const struct gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out);
#endif
