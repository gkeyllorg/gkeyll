#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_triad_geom.h>

/* Structure needed for dg-gr-maxwell geometry */
struct gkyl_dg_gr_maxwell_geom_ctx {
  double mass_bh;
  double spin_bh;
};

/**
 * Preset function for the lapse
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_dg_gr_maxwell_preset_lapse(enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for the shift
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_dg_gr_maxwell_preset_shift(enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for contravariant geometric source factors.
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_dg_gr_maxwell_preset_geom_factor_con(enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for the h_ij
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_dg_gr_maxwell_preset_h_ij(enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for the h_ij_inv
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_dg_gr_maxwell_preset_h_ij_inv(enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for the det_h
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_dg_gr_maxwell_preset_det_h(enum gkyl_triad_preset_geom_type preset_geom_type);
