#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_mom_calc_bcorr gkyl_mom_calc_bcorr;

/**
 * Create new updater to update boundary corrections.
 *
 * @param grid Grid object
 * @param momt Moment type object for boundary correction
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_calc_bcorr* 
gkyl_mom_calc_bcorr_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt, bool use_gpu);

/**
 * Compute boundary correction moments.
 *
 * @param bcorr Boundary correction updater object
 * @param phase_rng Phase space range on which to compute.
 * @param conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param out Output
 */
void gkyl_mom_calc_bcorr_advance(const struct gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out);
  
/**
 * Delete updater.
 *
 * @param bcorr Updater to delete.
 */
void gkyl_mom_calc_bcorr_release(struct gkyl_mom_calc_bcorr* up);

// "derived" class constructors
struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, bool use_gpu);

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_pkpm_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, double mass, bool use_gpu);

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  double mass, const struct gkyl_velocity_map *vel_map,
  bool use_gpu);
