#pragma once

#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_array.h>
#include <gkyl_util.h>

// Object type.
typedef struct gkyl_positivity_shift_gyrokinetic gkyl_positivity_shift_gyrokinetic;

// Configuration-space regions, specified in z-coordinates, in which to restrict
// the positivity shift. The shift is only applied where the coordinate along
// the field line (z, the last configuration-space direction) lies within one of
// the regions. Region r covers the interval [lower[r], upper[r]]. If num_regions
// is 0 (the default) the shift is applied everywhere.
struct gkyl_positivity_shift_gyrokinetic_regions {
  int num_regions; // Number of regions; 0 applies the shift everywhere.
  double lower[GKYL_MAX_POSITIVITY_SHIFT_REGIONS]; // Lower z-coordinate of each region.
  double upper[GKYL_MAX_POSITIVITY_SHIFT_REGIONS]; // Upper z-coordinate of each region.
};

/**
 * Create a new upater which shifts f up to a floor wherever f<0, and computes
 * contribution to the integrated moments of the shift in f in each cell (the
 * final reduction needs to be performed outside of this updater, as it may
 * involve a communicator).
 *
 * @param cbasis Conf-space basis.
 * @param pbasis Phase-space basis.
 * @param pgrid Phase-space grid.
 * @param mass Mass of species.
 * @param gk_geom Geometry object.
 * @param vel_map Velocity mapping object.
 * @param conf_rng_ext Extended configuration space range.
 * @param shift_regions Conf-space regions (in z-coordinates) to restrict the shift to.
 * @param use_gpu bool to determine if on GPU.
 * @return New positivity shift updater pointer.
 */
struct gkyl_positivity_shift_gyrokinetic *gkyl_positivity_shift_gyrokinetic_new(
  struct gkyl_basis cbasis, struct gkyl_basis pbasis, struct gkyl_rect_grid grid, double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map,
  const struct gkyl_range *conf_rng_ext,
  struct gkyl_positivity_shift_gyrokinetic_regions shift_regions, bool use_gpu
);

/**
 * Run the positivity shift updater in the indicated range.
 *
 * @param up Positivity shift updater.
 * @param conf_rng Config-space range.
 * @param phase_rng Phase-space range.
 * @param distf Distribution function array.
 * @param m0 Output M0 moment array.
 * @param delta_m0 M0 moment of the shift in f.
 */
void gkyl_positivity_shift_gyrokinetic_advance(
  gkyl_positivity_shift_gyrokinetic *up, const struct gkyl_range *conf_rng,
  const struct gkyl_range *phase_rng, struct gkyl_array *GKYL_RESTRICT distf,
  struct gkyl_array *GKYL_RESTRICT m0, struct gkyl_array *GKYL_RESTRICT delta_m0
);

/**
 * Scale the distribution function in the indicated range after applying the
 * positivity shift.
 *
 * @param up Positivity shift updater.
 * @param conf_rng Config-space range.
 * @param phase_rng Phase-space range.
 * @param delta_m0s M0 moment of the shift in fs.
 * @param delta_m0s_tot Sum of M0 moment of the shift in species with same charge sign.
 * @param delta_m0r_tot Sum of M0 moment of the shift in species with opposite charge sign.
 * @param m0s M0 moment of fs.
 * @param fs Distribution function of the species we wish to scale.
 */
void gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(
  gkyl_positivity_shift_gyrokinetic *up, const struct gkyl_range *conf_rng,
  const struct gkyl_range *phase_rng, const struct gkyl_array *GKYL_RESTRICT delta_m0s,
  const struct gkyl_array *GKYL_RESTRICT delta_m0s_tot,
  const struct gkyl_array *GKYL_RESTRICT delta_m0r_tot, const struct gkyl_array *GKYL_RESTRICT m0s,
  struct gkyl_array *GKYL_RESTRICT fs
);

/**
 * Release the memory associated with this positivity shift updater.
 *
 * @param up Positivity shift updater.
 */
void gkyl_positivity_shift_gyrokinetic_release(gkyl_positivity_shift_gyrokinetic *up);
