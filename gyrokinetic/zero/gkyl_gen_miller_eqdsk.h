#pragma once

#include <stdbool.h>

typedef struct gkyl_gen_miller_eqdsk gkyl_gen_miller_eqdsk;

struct gkyl_gen_miller_eqdsk_inp {
  double Rmajor_axis, Zvert_axis; // Major radius and vertical position of magnetic axis.
  double bmag_axis; // Signed toroidal field at the magnetic axis [T]; nonzero, sign sets B_phi.
  double Rmajor_LCFS; // Major radius of LCFS at OMP.
  double kappa; // Elongation (=0 sets kappa=1).
  double delta; // Triangularity.
  double a_shift; // Strength of Shafranov shift.
  double (*qsafety)(double R_omp, void *ctx); // Safety factor (>0) as a function of major radius.
  void *qsafety_ctx;
  bool has_Rmajor_omp_bounds;
  double Rmajor_omp_bounds[2]; // Optional requested simulation surfaces, as major radii [m].
  int num_nodes_Rmajor; // Number of major-radius mesh nodes; odd, 5..999 (0 selects 257).
  int num_nodes_Zvert; // Number of vertical mesh nodes; odd, 5..999 (0 selects 257).
  double rtol, guard_fraction; // Defaults 1e-10 and 0.05.
  bool has_RmajorZvert_bounds;
  double RmajorZvert_bounds[4]; // Rmin, Rmax, Zmin, Zmax [m].
};

struct gkyl_gen_miller_eqdsk_info {
  double psi_axis; // Poloidal flux divided by 2*pi at the magnetic axis [Wb/rad], fixed to zero.
  double psi_boundary; // Poloidal flux divided by 2*pi at the LCFS [Wb/rad].
  double psi_guard; // Poloidal flux divided by 2*pi at the outer physical guard [Wb/rad].
  double rminor_LCFS; // Miller minor radius of the LCFS [m].
  double rminor_required; // Larger of LCFS and outer requested Miller minor radii [m].
  double rminor_guard; // Outer physical minor radius: rminor_required*(1+guard_fraction) [m].
  double RmajorZvert_bounds[4]; // Resolved rectangular mesh bounds: Rmin, Rmax, Zmin, Zmax [m].
  // Maximum midpoint interpolation error in g=dpsi/d(r^2), normalized by max(g,1e-8*max(g_nodes)).
  double radial_error;
  // Minimum sampled LCFS/outer-required contour distance to guard or box, in R-Z mesh spacings.
  // Populated by advance; at least 3 after a successful advance.
  double clearance_cells;
  int num_nodes_radial; // Number of nodes in the adaptive s=r^2 table used to integrate poloidal flux.
  int num_nodes_Rmajor; // Resolved number of major-radius nodes in the rectangular GEQDSK mesh.
  int num_nodes_Zvert; // Resolved number of vertical nodes in the rectangular GEQDSK mesh.
};

/**
 * Construct the host generator and radial flux table; report failures to stderr.
 *
 * @param inp Physical and numerical inputs; callback context must outlive the generator.
 * @return New generator, or NULL on failure.
 */
gkyl_gen_miller_eqdsk *gkyl_gen_miller_eqdsk_inew(const struct gkyl_gen_miller_eqdsk_inp *inp);

/**
 * Sample and validate the rectangular mesh, using storage allocated at construction.
 *
 * @param up Generator.
 * @return True on success, false on failure with a diagnostic on stderr.
 */
bool gkyl_gen_miller_eqdsk_advance(gkyl_gen_miller_eqdsk *up);

/**
 * Write a previously sampled GEQDSK atomically in its destination directory.
 *
 * @param up Generator after a successful advance.
 * @param filename Destination path.
 * @return True on success, false on failure with a diagnostic on stderr.
 */
bool gkyl_gen_miller_eqdsk_write(const gkyl_gen_miller_eqdsk *up, const char *filename);

/**
 * Get generator metadata.
 *
 * @param up Generator.
 * @return Pointer to metadata owned by up; valid until up is released.
 */
const struct gkyl_gen_miller_eqdsk_info *gkyl_gen_miller_eqdsk_get_info(
  const gkyl_gen_miller_eqdsk *up
);

/**
 * Evaluate flux at a Miller minor radius.
 *
 * @param up Generator.
 * @param radius Miller minor radius [m].
 * @return Poloidal flux divided by 2*pi [Wb/rad], or NaN outside physical coverage.
 */
double gkyl_gen_miller_eqdsk_psi_at_r(const gkyl_gen_miller_eqdsk *up, double radius);

/**
 * Invert poloidal flux to Miller minor radius.
 *
 * @param up Generator.
 * @param psi Poloidal flux divided by 2*pi [Wb/rad].
 * @return Miller minor radius [m], or NaN outside physical coverage.
 */
double gkyl_gen_miller_eqdsk_r_at_psi(const gkyl_gen_miller_eqdsk *up, double psi);

/**
 * Convert normalized poloidal arc length to Miller poloidal angle.
 * This uses the default IWL full-contour convention: zero at the outboard midplane,
 * +/-pi at the inboard midplane, including SOL contours cut at the inboard midplane.
 * Custom plates, position maps, and other tokamak block conventions are excluded.
 * The conversion uses the analytic Miller contour, so numerical GEQDSK interpolation
 * introduces a small discrepancy when comparing with the tokamak geometry.
 *
 * @param up Generator with constant shaping parameters.
 * @param theta_tok Tokamak normalized poloidal arc coordinate in [-pi, pi].
 * @return Miller angle [rad], or NaN outside the stated interval.
 */
double gkyl_gen_miller_eqdsk_theta_at_theta_tok(const gkyl_gen_miller_eqdsk *up, double theta_tok);

/**
 * Evaluate flux at a physical outboard midplane major radius.
 *
 * @param up Generator.
 * @param R_omp Outboard midplane major radius in meters.
 * @return Poloidal flux divided by 2*pi [Wb/rad], or NaN outside physical coverage.
 */
double gkyl_gen_miller_eqdsk_psi_at_R_omp(const gkyl_gen_miller_eqdsk *up, double R_omp);

/**
 * Invert physical flux to outboard midplane major radius.
 *
 * @param up Generator.
 * @param psi Poloidal flux divided by 2*pi [Wb/rad].
 * @return Outboard midplane major radius in meters, or NaN outside physical coverage.
 */
double gkyl_gen_miller_eqdsk_R_omp_at_psi(const gkyl_gen_miller_eqdsk *up, double psi);

/**
 * Evaluate the Miller surface at minor radius r and poloidal angle theta.
 *
 * @param up Generator.
 * @param r Miller minor radius in meters.
 * @param theta Poloidal angle in radians; zero is the outboard midplane.
 * @param rz Output array containing cylindrical R,Z coordinates in meters.
 */
void gkyl_gen_miller_eqdsk_RZ_at_rtheta(
  const gkyl_gen_miller_eqdsk *up, double r, double theta, double rz[2]
);

/**
 * Release a generator and all its host storage.
 *
 * @param up Generator to release, or NULL.
 */
void gkyl_gen_miller_eqdsk_release(gkyl_gen_miller_eqdsk *up);
