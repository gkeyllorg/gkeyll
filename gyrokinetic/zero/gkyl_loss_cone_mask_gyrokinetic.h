#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type.
typedef struct gkyl_loss_cone_mask_gyrokinetic gkyl_loss_cone_mask_gyrokinetic;

/** Boundary behavior used by the orbit classifier. */
enum gkyl_loss_cone_boundary_type {
  GKYL_LOSS_CONE_BC_OPEN = 0, // A particle that reaches the domain edge is lost.
  GKYL_LOSS_CONE_BC_SHEATH, // The material-wall potential may reflect the particle.
  GKYL_LOSS_CONE_BC_CLOSED // No particle can be lost through this edge.
};

// Inputs packaged as a struct.
struct gkyl_loss_cone_mask_gyrokinetic_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  bool use_gpu; // Flag to indicate if GPU should be used.
  double mass; // Species mass.
  double charge; // Species charge.
  enum gkyl_loss_cone_boundary_type lower_boundary; // Lower parallel-boundary behavior.
  enum gkyl_loss_cone_boundary_type upper_boundary; // Upper parallel-boundary behavior.
};

/**
 * Create new updater that populates an array with the masking function
 * based on the escape-barrier criterion. For sheath boundaries the material
 * wall is a virtual endpoint outside the plasma domain, so
 *   EB_L = max(max_{s in [z_L,z]} U(s,mu), mu*B(z_L)+q*phi_wall_L)
 *   EB_R = max(max_{s in [z,z_R]} U(s,mu), mu*B(z_R)+q*phi_wall_R)
 * with U = mu*B + q*phi. Open boundaries omit the virtual endpoint and
 * closed boundaries cannot be escaped.
 * A node is trapped if H < EB.
 *
 * @param inp Input parameters defined in gkyl_loss_cone_mask_gyrokinetic_inp struct.
 * @return New updater pointer.
 */
struct gkyl_loss_cone_mask_gyrokinetic *
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp);

/**
 * Compute the loss-cone mask on phase-space cell nodes.
 *
 * The caller supplies the magnetic field magnitude and electrostatic
 * potential arrays. This keeps the updater free of any communication logic;
 * the app is responsible for assembling global data when needed.
 *
 * @param up Project on basis updater to run.
 * @param phase_rng Phase-space range.
 * @param conf_rng Configuration-space range.
 * @param bmag Magnetic field magnitude (cdim DG expansion).
 * @param phi Electrostatic potential at the plasma/sheath entrance.
 * @param phi_wall_lo Lower material-wall potential. Required for a lower sheath boundary.
 * @param phi_wall_up Upper material-wall potential. Required for an upper sheath boundary.
 * @param mask_out Cellwise mask: 1 if every cell corner is trapped, 0 if any
 *                 corner can escape through either boundary.
 */
void gkyl_loss_cone_mask_gyrokinetic_advance(
  gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic *up);
