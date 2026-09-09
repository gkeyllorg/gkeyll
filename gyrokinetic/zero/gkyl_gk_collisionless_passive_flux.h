#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_bc_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>

// Object type.
typedef struct gkyl_gk_collisionless_passive_flux gkyl_gk_collisionless_passive_flux;

/**
 * Create new updater to compute the surface expansion of the phase-space
 * flux for passive (prescribed) conf-space advection.
 *
 * @param phase_grid Phase space grid.
 * @param conf_basis Configuration space basis functions.
 * @param phase_basis Phase space basis functions.
 * @param passive_speeds Conf-space advection speeds.
 * @param charge Species charge.
 * @param mass Species mass.
 * @param gk_geom Gyrokinetic geometry object.
 * @param dg_geom DG geometry object.
 * @param gk_dg_geom GK DG geometry object.
 * @param vel_map Velocity space mapping object.
 * @param bctype_conf Type of conf-space BCs (for choosing edge kernels).
 * @param use_gpu Whether to allocate on GPU.
 * @return New updater pointer.
 */
gkyl_gk_collisionless_passive_flux*
gkyl_gk_collisionless_passive_flux_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_array *passive_speeds,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom,
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf, bool use_gpu);

/**
 * Compute the surface expansion of the passive phase-space flux.
 *
 * Each cell stores the *lower* surface expansion.
 *
 * @param up Updater.
 * @param conf_range Configuration space range (local only; geometry defined here).
 * @param phase_range Phase space range.
 * @param conf_ext_range Extended configuration space range.
 * @param phase_ext_range Extended phase space range.
 * @param fin Distribution function.
 * @param flux_surf Output surface expansion (lower edge in each direction).
 * @param cflrate Output CFL rate.
 */
void gkyl_gk_collisionless_passive_flux_surf(gkyl_gk_collisionless_passive_flux *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range,
  const struct gkyl_array *fin,
  struct gkyl_array *flux_surf, struct gkyl_array *cflrate);

/**
 * Delete pointer to passive flux updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_collisionless_passive_flux_release(gkyl_gk_collisionless_passive_flux *up);
