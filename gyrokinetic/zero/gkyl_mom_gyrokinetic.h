#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

// Object type.
typedef struct gkyl_mom_gyrokinetic gkyl_mom_gyrokinetic;

/**
 * Create new updater to compute gyrokinetic moments.
 *
 * @param mass Positive species mass.
 * @param charge Species charge.
 * @param cbasis Configuration-space serendipity basis.
 * @param pbasis Phase-space gkhybrid (p1) or serendipity (p2) basis.
 * @param phase_grid Phase-space grid; metadata is copied.
 * @param vel_map Velocity map; acquired by the updater.
 * @param gk_geom Geometry; acquired by the updater.
 * @param mom_type Name of moment to compute.
 * @param is_integrated Whether to compute the volume integrated moment.
 * @param use_gpu bool to determine if on GPU.
 * @return A new moment updater.
 */
struct gkyl_mom_gyrokinetic *gkyl_mom_gyrokinetic_new(
  double mass, double charge, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, enum gkyl_distribution_moments mom_type, bool is_integrated,
  bool use_gpu
);

/**
 * Fetch the number of moments
 *
 * @param up Moment updater.
 *
 * @return number of moments
 */
int gkyl_mom_gyrokinetic_num_mom(const struct gkyl_mom_gyrokinetic *up);

/**
 * Compute a moment of the gyrokinetic distribution function.
 * The ranges must retain the strides of their backing arrays. conf_rng must lie
 * within the configuration extent of phase_rng and use the geometry/phi layout.
 * Only conf_rng is overwritten; other cells are untouched. Integrated output has
 * one scalar per moment per configuration cell (no global reduction). All arrays
 * must use the updater's backend. The updater owns no input or output arrays.
 *
 * @param up Moment updater.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param phi Current potential; required for Hamiltonian moments, otherwise NULL allowed.
 * @param fin Input distribution function array.
 * @param mout Output array: num_mom*cbasis.num_basis components, or num_mom if integrated.
 */
void gkyl_mom_gyrokinetic_advance(
  struct gkyl_mom_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
);

/**
 * Free memory associated with the moment calculation.
 *
 * @param up Moment updater.
 */
void gkyl_mom_gyrokinetic_release(struct gkyl_mom_gyrokinetic *up);
