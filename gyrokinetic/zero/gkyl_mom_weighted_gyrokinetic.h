#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

enum gkyl_mom_weight_type {
  GKYL_F_MOMENT_WEIGHT_NONE = 0, // No weight.
  GKYL_F_MOMENT_WEIGHT_CONF, // Weight is a function of conf-space.
  GKYL_F_MOMENT_WEIGHT_VEL, // Weight is a function of v-space.
  GKYL_F_MOMENT_WEIGHT_PHASE // Weight is a function of phase-space.
};

// Object type.
typedef struct gkyl_mom_weighted_gyrokinetic gkyl_mom_weighted_gyrokinetic;

/**
 * Create a new updater to compute weighted gyrokinetic moments.
 *
 * @param mass Positive species mass.
 * @param charge Species charge.
 * @param cbasis Configuration-space serendipity basis.
 * @param pbasis Phase-space gkhybrid (p1) or serendipity (p2) basis.
 * @param phase_grid Phase-space grid; metadata is copied.
 * @param vel_map Velocity map; acquired by the updater.
 * @param gk_geom Geometry; acquired by the updater.
 * @param mom_type Name of moment to compute.
 * @param wgt_type CONF or VEL; NONE and PHASE are rejected.
 * @param is_integrated Whether to compute the volume integrated moment.
 * @param use_gpu bool to determine if on GPU.
 * @return A new weighted moment updater.
 */
struct gkyl_mom_weighted_gyrokinetic *gkyl_mom_weighted_gyrokinetic_new(
  double mass, double charge, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, enum gkyl_distribution_moments mom_type,
  enum gkyl_mom_weight_type wgt_type, bool is_integrated, bool use_gpu
);

/**
 * Fetch the number of moments
 *
 * @param up Moment updater.
 * 
 * @return number of moments
 */
int gkyl_mom_weighted_gyrokinetic_num_mom(const struct gkyl_mom_weighted_gyrokinetic *up);

/**
 * Compute a weighted moment of the gyrokinetic distribution function.
 * The ranges must retain the strides of their backing arrays. conf_rng must lie
 * within the configuration extent of phase_rng and use the geometry/phi layout.
 * Only conf_rng is overwritten; other cells are untouched. Integrated output has
 * one scalar per moment per configuration cell (no global reduction). All arrays
 * must use the updater's backend. The updater owns no input or output arrays.
 *
 * @param up Weighted moment updater.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param wgt_rng Weight array range, covering the requested configuration or velocity cells.
 * @param wgt Weight in cbasis (CONF) or gkhybrid_vel p1/serendipity p2 (VEL).
 * @param phi Optional potential, unused by the supported raw weighted moments; NULL allowed.
 * @param fin Input distribution function array.
 * @param mout Output array: num_mom*cbasis.num_basis components, or num_mom if integrated.
 */
void gkyl_mom_weighted_gyrokinetic_advance(
  struct gkyl_mom_weighted_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_range *wgt_rng,
  const struct gkyl_array *GKYL_RESTRICT wgt, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
);

/**
 * Free memory associated with the weighted moment calculation.
 *
 * @param up Weighted moment updater.
 */
void gkyl_mom_weighted_gyrokinetic_release(struct gkyl_mom_weighted_gyrokinetic *up);
