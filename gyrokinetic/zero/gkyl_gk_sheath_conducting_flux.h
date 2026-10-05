#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_collisionless_flux.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Object type
typedef struct gkyl_gk_sheath_conducting_flux gkyl_gk_sheath_conducting_flux;

struct gkyl_gk_sheath_conducting_flux_inp {
  int dir; // Configuration-space direction normal to the sheath boundary.
  int cdim; // Number of configuration-space dimensions.
  enum gkyl_edge_loc edge; // Lower or upper boundary edge.
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *phase_basis; // Host phase-space basis.
  const struct gkyl_basis *phase_basis_cu; // Device phase-space basis, when use_gpu is true.
  const struct gkyl_range *conf_range; // Local configuration-space range.
  const struct gkyl_range *conf_ext_range; // Extended local configuration-space range.
  const struct gkyl_range *phase_ext_range; // Extended local phase-space range.
  const struct gkyl_range *skin_range; // Phase-space range of boundary skin cells.
  const struct gkyl_gk_collisionless_flux *flux_op; // Collisionless phase-space flux updater.
  double charge; // Species charge.
  double mass; // Species mass.
  bool use_gpu; // Whether to use GPU acceleration.
};

/**
 * Create a new updater to add flux-balanced conducting sheath return flux.
 *
 * @param inp gk_sheath_conducting_flux_inp struct containing updater inputs.
 * @return New updater pointer.
 */
gkyl_gk_sheath_conducting_flux *gkyl_gk_sheath_conducting_flux_inew(
  const struct gkyl_gk_sheath_conducting_flux_inp *inp
);

/**
 * Add flux-balanced conducting sheath return flux to the collisionless surface flux.
 *
 * @param up Conducting sheath flux updater.
 * @param phi Electrostatic potential used to set the sheath energy cutoff.
 * @param gyro_phi Gyroaveraged electrostatic potential.
 * @param phi_wall Wall potential at the sheath boundary.
 * @param fin Input distribution function.
 * @param flux_surf Collisionless phase-space surface flux modified in place.
 */
void gkyl_gk_sheath_conducting_flux_advance(
  gkyl_gk_sheath_conducting_flux *up, const struct gkyl_array *phi,
  const struct gkyl_array *gyro_phi, const struct gkyl_array *phi_wall,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
);

/**
 * Delete conducting sheath flux updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_sheath_conducting_flux_release(gkyl_gk_sheath_conducting_flux *up);
