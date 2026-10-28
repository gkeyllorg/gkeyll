#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

typedef struct gkyl_orbit_average_gyrokinetic gkyl_orbit_average_gyrokinetic;

/** Result of constructing a frozen orbit projector. */
enum gkyl_orbit_average_gyrokinetic_status {
  GKYL_ORBIT_AVERAGE_SUCCESS = 0,
  GKYL_ORBIT_AVERAGE_UNSUPPORTED,
  GKYL_ORBIT_AVERAGE_INVALID_INPUT,
  GKYL_ORBIT_AVERAGE_SINGULAR
};

/**
 * Inputs for a frozen-field, 1x2v P1 gyrokinetic orbit projector.
 * All numeric inputs and array coefficients must be finite. Grid spacings,
 * mass, B, and geometry/velocity Jacobians must be positive; mapped mu and
 * separatrix_margin must be nonnegative. Evaluated fields and intermediate
 * energies and metric weights must remain representable as finite doubles.
 */
struct gkyl_orbit_average_gyrokinetic_inp {
  const struct gkyl_rect_grid *grid;
  const struct gkyl_basis *conf_basis, *phase_basis;
  const struct gkyl_range *phase_range, *conf_range;
  const struct gkyl_velocity_map *vel_map;
  const struct gkyl_array *bmag, *phi, *jacobgeo;
  const struct gkyl_array *mask; // Cellwise candidate mask; zero cells remain unchanged.
  const struct gkyl_array *phi_wall_lo, *phi_wall_up;
  enum gkyl_gk_trapped_passing_orbit_type lower_orbit, upper_orbit;
  double mass, charge; // Finite mass > 0 and finite signed charge.
  double separatrix_margin; // Energy buffer at escape barriers and internal saddles, >=0.
  int num_energy_nodes; // At least two. Uniform base knots are split at critical energies.
  bool frequency_is_infinite; // Use the projection limit; ignore the advance frequency.
  bool positivity; // Limit relaxation at volume quadrature and tensor cell corners.
  bool use_gpu;
};

/** Diagnostics for one implicit update. A failed update does not change f. */
struct gkyl_orbit_average_gyrokinetic_result {
  bool success;
  double min_blend; // Smallest actual correction fraction, in [0,1].
  long num_limited_blocks;
};

/** Frozen setup statistics. */
struct gkyl_orbit_average_gyrokinetic_info {
  long num_selected_cells;
  long num_rejected_cells;
  long num_modes; // Total rank after removing dependent projected invariant columns.
  int num_blocks; // Independent magnetic-moment cells.
};

/**
 * Construct a DG-consistent invariant-space projector. Fields and geometry
 * are frozen at construction; release and reconstruct after any change.
 * The supplied range must cover complete parallel orbits. MPI communication
 * belongs to the caller. Input arrays can be on host or device: only frozen
 * setup data are copied to the host, never the evolving distribution.
 *
 * The stored distribution is g=J_geo B J_vel f. Physical f is reconstructed
 * at Gauss (4,4,2) nodes and projected in the corresponding invariant-volume
 * metric. No generated DG kernels are introduced.
 *
 * @param inp Frozen geometry, ranges, candidate mask, and numerical options.
 * @param status Optional construction status; NULL is allowed.
 * @return New updater, or NULL for unsupported or invalid input.
 */
gkyl_orbit_average_gyrokinetic *gkyl_orbit_average_gyrokinetic_inew(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp,
  enum gkyl_orbit_average_gyrokinetic_status *status
);

/**
 * Apply f <- f + beta(Pf-f), beta=dt*frequency/(1+dt*frequency), in place.
 * A single positivity factor is used per coupled magnetic-moment block.
 * dt must be finite and nonnegative. With frequency_is_infinite=false,
 * frequency must also be finite and nonnegative; frequency=0 is an exact no-op.
 * With frequency_is_infinite=true, frequency is ignored and dt>0 requests
 * the projection limit. dt=0 is always an exact no-op.
 * Distribution coefficients and their evaluated values must be finite, with
 * representable intermediate reductions. Positivity limiting additionally
 * requires nonnegative values at monitored nodes. Finiteness is not checked.
 * No allocation or geometry reconstruction occurs during this update.
 *
 * @param up Frozen projector.
 * @param dt Accepted time step.
 * @param frequency Orbit relaxation rate.
 * @param f Jacobian-weighted DG distribution on the constructor's phase range.
 * @return Success and limiter diagnostics; negative monitored input fails
 *         transactionally if positivity limiting is enabled.
 */
struct gkyl_orbit_average_gyrokinetic_result gkyl_orbit_average_gyrokinetic_advance(
  gkyl_orbit_average_gyrokinetic *up, double dt, double frequency, struct gkyl_array *f
);

/** Return statistics for the frozen projector. */
struct gkyl_orbit_average_gyrokinetic_info gkyl_orbit_average_gyrokinetic_get_info(
  const gkyl_orbit_average_gyrokinetic *up
);

/** Free all projector and work arrays. */
void gkyl_orbit_average_gyrokinetic_release(gkyl_orbit_average_gyrokinetic *up);
