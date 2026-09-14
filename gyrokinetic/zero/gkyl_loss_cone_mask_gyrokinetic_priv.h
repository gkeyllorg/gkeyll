// Private header: not for direct use
#pragma once

#include <float.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct gkyl_loss_cone_mask_gyrokinetic {
  int cdim; // Configuration-space dimension.
  int num_basis_conf; // Number of configuration-space basis functions.
  struct gkyl_range conf_corner_range; // Tensor-product P1 corner indices (0 or 1 per direction).
  bool use_gpu;

  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  double mass; // Species mass.
  double charge; // Species charge.
  enum gkyl_gk_trapped_passing_orbit_type lower_orbit; // Lower-z orbit type.
  enum gkyl_gk_trapped_passing_orbit_type upper_orbit; // Upper-z orbit type.

  struct gkyl_array *basis_at_corners_conf; // Basis functions at configuration-space P1 corners.
};

/** Return 0 or 1 for the lower/upper z side of a range-linear P1 corner. */
GKYL_CU_DH
static inline int
gk_lcm_conf_corner_z_side(const struct gkyl_range *conf_corner_range,
  int conf_corner, int zdim)
{
  int corner_idx[GKYL_MAX_CDIM] = { 0 };
  gkyl_range_inv_idx(conf_corner_range, conf_corner, corner_idx);
  return corner_idx[zdim];
}

/**
 * Return the P1 corner with the anchor's transverse coordinates and the
 * requested lower (0) or upper (1) side in zdim.
 */
GKYL_CU_DH
static inline int
gk_lcm_conf_corner_with_z_side(const struct gkyl_range *conf_corner_range,
  int anchor_corner, int zdim, int z_side)
{
  int corner_idx[GKYL_MAX_CDIM] = { 0 };
  gkyl_range_inv_idx(conf_corner_range, anchor_corner, corner_idx);
  corner_idx[zdim] = z_side;
  return gkyl_range_idx(conf_corner_range, corner_idx);
}

/** Evaluate a configuration-space DG field at a range-linear P1 cell corner. */
GKYL_CU_DH
static inline double
gk_lcm_field_corner_value(const struct gkyl_array *arr,
  const struct gkyl_array *basis_at_corners, int num_basis, long linidx,
  int corner)
{
  const double *arr_d = (const double *)gkyl_array_cfetch(arr, linidx);
  const double *basis_d = (const double *)gkyl_array_cfetch(basis_at_corners, corner);

  double val = 0.0;
  for (int k = 0; k < num_basis; ++k) {
    val += arr_d[k] * basis_d[k];
  }

  return val;
}

/**
 * Return true when the particle Hamiltonian is strictly below one directional
 * escape barrier. Equality means that the particle reaches the loss endpoint
 * with zero parallel kinetic energy.
 */
GKYL_CU_DH
static inline bool
gk_lcm_hamiltonian_below_barrier(double kinetic_energy, double magnetic_energy,
  double electric_energy, double barrier)
{
  if (barrier == DBL_MAX) {
    return true;
  }

  double hamiltonian = kinetic_energy + magnetic_energy + electric_energy;
  double energy_scale = fabs(kinetic_energy) + fabs(magnetic_energy)
    + fabs(electric_energy) + fabs(barrier);
  double roundoff_tol = 32.0 * DBL_EPSILON * energy_scale;

  // Treat equality, including equality to within modal-evaluation roundoff,
  // as reaching the boundary rather than reflecting from it.
  return hamiltonian < barrier - roundoff_tol;
}

/**
 * Compute the maximum effective potential U=mu*B+q*phi along the lower- and
 * upper-z orbits from one configuration-space P1 corner. base_idx fixes the
 * transverse cell indices, z_cell is the anchor cell, and a sheath-trapped
 * orbit includes the material-wall potential.
 */
GKYL_CU_DH
static inline void
gk_lcm_escape_barriers(int cdim, int num_basis_conf,
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_corner_range,
  const struct gkyl_array *basis_at_corners_conf, const struct gkyl_array *phi,
  const struct gkyl_array *bmag, const struct gkyl_array *phi_wall_lo,
  const struct gkyl_array *phi_wall_up, const int *base_idx, int z_cell,
  int anchor_conf_corner, double mu, double charge,
  enum gkyl_gk_trapped_passing_orbit_type lower_orbit,
  enum gkyl_gk_trapped_passing_orbit_type upper_orbit,
  double *barrier_left, double *barrier_right)
{
  int zdim = cdim - 1;

  int z_side = gk_lcm_conf_corner_z_side(conf_corner_range, anchor_conf_corner, zdim);
  int z_upper_corner = gk_lcm_conf_corner_with_z_side(conf_corner_range,
    anchor_conf_corner, zdim, 1);
  int z_lower_corner = gk_lcm_conf_corner_with_z_side(conf_corner_range,
    anchor_conf_corner, zdim, 0);

  int scan_idx[GKYL_MAX_CDIM];
  for (int d = 0; d < cdim; ++d) {
    scan_idx[d] = base_idx[d];
  }

  *barrier_left = -DBL_MAX;
  *barrier_right = -DBL_MAX;

  for (int iz = conf_range->lower[zdim]; iz <= conf_range->upper[zdim]; ++iz) {
    scan_idx[zdim] = iz;
    long linidx = gkyl_range_idx(conf_range, scan_idx);

    double phi_lower = gk_lcm_field_corner_value(phi, basis_at_corners_conf,
      num_basis_conf, linidx, z_lower_corner);
    double phi_upper = gk_lcm_field_corner_value(phi, basis_at_corners_conf,
      num_basis_conf, linidx, z_upper_corner);
    double bmag_lower = gk_lcm_field_corner_value(bmag, basis_at_corners_conf,
      num_basis_conf, linidx, z_lower_corner);
    double bmag_upper = gk_lcm_field_corner_value(bmag, basis_at_corners_conf,
      num_basis_conf, linidx, z_upper_corner);
    double u_lower = mu * bmag_lower + charge * phi_lower;
    double u_upper = mu * bmag_upper + charge * phi_upper;

    // A complete cell contributes both endpoint traces. In the anchor cell,
    // only the portion between the anchor corner and the relevant wall belongs
    // to that directional path. Keeping both traces also handles discontinuous
    // input fields conservatively.
    if (iz < z_cell) {
      *barrier_left = GKYL_MAX2(*barrier_left, GKYL_MAX2(u_lower, u_upper));
    }
    else if (iz == z_cell) {
      *barrier_left = GKYL_MAX2(*barrier_left, u_lower);
      *barrier_right = GKYL_MAX2(*barrier_right, u_upper);
      if (z_side == 1) {
        *barrier_left = GKYL_MAX2(*barrier_left, u_upper);
      }
      else {
        *barrier_right = GKYL_MAX2(*barrier_right, u_lower);
      }
    }
    else {
      *barrier_right = GKYL_MAX2(*barrier_right, GKYL_MAX2(u_lower, u_upper));
    }
  }

  if (lower_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL) {
    *barrier_left = DBL_MAX;
  }
  else if (lower_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH) {
    scan_idx[zdim] = conf_range->lower[zdim];
    long linidx = gkyl_range_idx(conf_range, scan_idx);
    double bmag_wall = gk_lcm_field_corner_value(bmag, basis_at_corners_conf,
      num_basis_conf, linidx, z_lower_corner);
    double phi_wall = gk_lcm_field_corner_value(phi_wall_lo, basis_at_corners_conf,
      num_basis_conf, linidx, z_lower_corner);
    *barrier_left = GKYL_MAX2(*barrier_left, mu * bmag_wall + charge * phi_wall);
  }

  if (upper_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL) {
    *barrier_right = DBL_MAX;
  }
  else if (upper_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH) {
    scan_idx[zdim] = conf_range->upper[zdim];
    long linidx = gkyl_range_idx(conf_range, scan_idx);
    double bmag_wall = gk_lcm_field_corner_value(bmag, basis_at_corners_conf,
      num_basis_conf, linidx, z_upper_corner);
    double phi_wall = gk_lcm_field_corner_value(phi_wall_up, basis_at_corners_conf,
      num_basis_conf, linidx, z_upper_corner);
    *barrier_right = GKYL_MAX2(*barrier_right, mu * bmag_wall + charge * phi_wall);
  }
}

#ifdef GKYL_HAVE_CUDA
void gkyl_loss_cone_mask_gyrokinetic_advance_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out);
#endif
