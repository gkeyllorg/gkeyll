#pragma once

#include <gkyl_gen_miller_eqdsk.h>

/**
 * Resolve flux or normalized-radius endpoints for either EQDSK flux orientation.
 * Uses signed flux differences; the axis and zero-width domains are excluded.
 *
 * @param type PSI_POL or RHO_POL representation; other values fail.
 * @param bounds Non-NULL finite, numerically increasing input endpoints.
 * @param psi_axis Magnetic-axis flux [Wb/rad].
 * @param psi_lcfs LCFS flux [Wb/rad], distinct from psi_axis.
 * @param psi_bounds Non-NULL output endpoints, sorted numerically; valid on success.
 * @return True on success, false for invalid or unrepresentable bounds.
 */
bool gen_miller_flux_bounds(
  enum gkyl_gk_geometry_radial_coord type, const double bounds[2], double psi_axis, double psi_lcfs,
  double psi_bounds[2]
);
