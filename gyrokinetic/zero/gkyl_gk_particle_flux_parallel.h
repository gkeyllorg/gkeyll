#pragma once

#include <gkyl_basis.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>
#include <gkyl_velocity_map.h>

struct gkyl_gk_particle_flux_parallel;

struct gkyl_gk_particle_flux_parallel_inp {
  enum gkyl_edge_loc edge; // Lower or upper parallel boundary.
  const struct gkyl_basis *cbasis; // Configuration-space basis.
  const struct gkyl_basis *pbasis; // Phase-space basis.
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_range *conf_skin; // Configuration-space skin range.
  const struct gkyl_range *phase_skin; // Phase-space skin range.
  double charge; // Species charge.
  double mass; // Species mass.
  const struct gk_geometry *gk_geom; // Gyrokinetic geometry object.
  const struct gkyl_dg_geom *dg_geom; // Nodal metric geometry.
  const struct gkyl_gk_dg_geom *gk_dg_geom; // Nodal magnetic geometry.
  const struct gkyl_velocity_map *vel_map; // Computational-to-physical velocity map.
  bool use_gpu; // Run on the GPU.
};

/**
 * Construct an updater that just computes the parallel boundary flux, assuming
 * zero drifts and absorbing boundary (so nothing comes in from the ghost cell).
 *
 * @param inp Constructor inputs; metadata is copied.
 * @return New updater.
 */
struct gkyl_gk_particle_flux_parallel *gkyl_gk_particle_flux_parallel_inew(
  const struct gkyl_gk_particle_flux_parallel_inp *inp
);

/**
 * Compute outward particle flux per transverse coordinate area from the skin
 * distribution, with zero incoming/ghost distribution.
 * The output has conf_skin.volume cells, each with 2^(cdim-1) surface
 * Gauss-Legender nodal values in collisionless-kernel order.
 *
 * @param up Updater.
 * @param fin Distribution with the phase_skin backing-array layout.
 * @param flux Surface nodal output on the updater's backend.
 */
void gkyl_gk_particle_flux_parallel_advance(
  struct gkyl_gk_particle_flux_parallel *up, const struct gkyl_array *fin, struct gkyl_array *flux
);

/**
 * Get (coordinate face area * |b dot n|) at surface Gauss-Legendre nodes.
 *
 * @param up Updater.
 * @return Read-only array with the same layout as the flux output.
 */
const struct gkyl_array *gkyl_gk_particle_flux_parallel_area(
  const struct gkyl_gk_particle_flux_parallel *up
);

/**
 * Release the updater and acquired references.
 *
 * @param up Updater.
 */
void gkyl_gk_particle_flux_parallel_release(struct gkyl_gk_particle_flux_parallel *up);
