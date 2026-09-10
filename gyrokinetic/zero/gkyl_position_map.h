#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>
#include <gkyl_util.h>

#include <math.h>

enum gkyl_position_map_id {
  GKYL_PMAP_USER_INPUT = 0, // Function projection. User specified. Default
  GKYL_PMAP_USER_INPUT_W_DERIVATIVE, // Function projection and derivative user specified
  GKYL_PMAP_CONSTANT_DB_POLYNOMIAL, // Makes a uniform dB in each cell. Polynomial approximation, assuming 2 local maxima in Bmag
  GKYL_PMAP_CONSTANT_DB_NUMERIC, // Makes a uniform dB in each cell, but calculates the dB numerically
  GKYL_PMAP_XPT_COMPRESSION, // Compresses cells near X-point (For use in MB Tokamaks)
};

typedef void (*mc2nu_t)(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx);

struct gkyl_position_map_inp {
  enum gkyl_position_map_id id;
  mc2nu_t maps[3]; // Position mapping in each position direction. This is defined in full 3x, 
  // not in deflated coordinates.
  mc2nu_t map_derivs[3]; // Derivative of mapping in each position direction. This is defined in full 3x, 
  // not in deflated coordinates.
  void *ctxs[3]; // Context for each position mapping function.
  double map_strength; // Zero is uniform mapping, one is fully nonuniform mapping. How strong the nonuniformity is
  // Call map_strength = s, xc computational coordinate, and xnu the nonuniform coordinate
  // xnu' = xnu * s + xc * (1-s)
  double maximum_slope_at_min_B; // The maximum slope of the mapping at a magnetic field minimum. A number > 1. Hard limits on cell sizes
  double maximum_slope_at_max_B; // The maximum slope of the mapping at a magnetic field maximum. A number > 1. Hard limits on cell sizes
  double gaussian_std; // The width of the moving average for the map to smooth it. Units of normalized field line length
  double gaussian_max_integration_width; // The maximum width to integrate the Gaussian filter. Units of normalized field line length
  double compression_factor; // For PMAP_XPT_Compression. Specifies how much smaller the cells are
                             // near the X-point
  double radial_compression_factor; // Factor by which cells are compressed radially near xpt
  bool compress_divertor; // Whether to apply compression at divertor plates for PMAP_XPT compression
};
struct gkyl_position_map_inew_inp {
  struct gkyl_position_map_inp pmap_info;
  struct gkyl_rect_grid grid; // Position space grid.
  struct gkyl_range local, local_ext; // Local & extended local position-space range.
  struct gkyl_range global, global_ext; // Global & extended global position-space range.
  struct gkyl_basis basis;  // Basis for position mapping.
};

struct gkyl_position_map {
  enum gkyl_position_map_id id;
  mc2nu_t maps[3]; // Position mapping in each position direction.
  mc2nu_t map_derivs[3]; // Derivative of position mapping in each position direction.
  void *ctxs[3]; // Context for each position mapping function.

  double cdim; // Number of computational dimensions.
  struct gkyl_rect_grid grid; // Position space grid.
  struct gkyl_range local, local_ext, global, global_ext; // Local & extended local position-space range.
  struct gkyl_basis basis;  // Basis for position mapping.
  struct gkyl_array *mc2nu; // Position mapping in each position direction.
  struct gkyl_ref_count ref_count;
  bool to_optimize; // Whether to optimize the position map for constant B mapping.
  bool use_map_derivs; // Whether to use analytical derivatives of the mapping

  // Stuff for constant B mapping
  struct gkyl_bmag_ctx *bmag_ctx; // Context for magnetic field calculation
  struct gkyl_position_map_const_B_ctx *constB_ctx; // Context for constant B mapping
  struct gkyl_position_map_xpt_ctx *xpt_ctx; // Context for X-point compression mapping

  struct gkyl_array *mc2nu_dev; // Device mirror of mc2nu (only with a device copy, see
                                // gkyl_position_map_make_cu_dev; kept in sync by set_mc2nu).
  struct gkyl_position_map *on_dev; // Device copy of itself (points to itself on CPU).
};

struct gkyl_position_map_const_B_ctx {
  mc2nu_t maps_backup[3]; // Backup of the position mapping functions.
  void *ctxs_backup[3]; // Backup of the context for each position mapping function.
  
  double psi, alpha; // The psi and alpha values for the middle flux surface to identify the 1D line we are optimizing 
  double psi_min, psi_max; // The max and min psi values for the simulation
  double alpha_min, alpha_max; // The max and min alpha values for the simulation
  double theta_min, theta_max; // The max and min theta values for the simulation
  double map_strength; // Zero is uniform mapping, one is fully nonuniform mapping. How strong the nonuniformity is
  bool enable_maximum_slope_limits_at_min_B; // Whether to enable the maximum slope limits at a magnetic field minimum
  bool enable_maximum_slope_limits_at_max_B; // Whether to enable the maximum slope limits at a magnetic field maximum
  double maximum_slope_at_min_B; // The maximum slope of the mapping at a magnetic field minimum
  double maximum_slope_at_max_B; // The maximum slope of the mapping at a magnetic field maximum
  double gaussian_std; // The standard deviation of the Gaussian filter used for smoothing the mapping
  double gaussian_max_integration_width; // The maximum width to integrate the Gaussian filter

  // Polynomial-based mapping
  double theta_throat, Bmag_throat; // The theta and Bmag values at the throat of the magnetic field
  int map_order_center, map_order_expander; // The polynomial order of the center and expander maps

  // Constant B mapping
  int N_theta_boundaries; // Number of times dB/dz changes sign
  int num_extrema; // Number of extrema in the magnetic field
  double *theta_extrema; // The theta values of the extrema
  double *bmag_extrema; // The Bmag values of the extrema
  bool *min_or_max; // Whether the extrema is a minima or maxima. 1 is maxima, 0 is minima
  double dB_cell; // The change in Bmag per cell
};

struct gkyl_position_map_xpt_ctx {
  mc2nu_t maps_backup[3]; // Backup of the position mapping functions.
  void *ctxs_backup[3]; // Backup of the context for each position mapping function.
  double compression_factor; // Factor by which cells near X-point are compressed
  double radial_compression_factor; // Factor by which cells are compressed in radial direction at separatrix
  bool compress_divertor; // Whether to apply compression at divertor plates
  double zcut; // Half-wavelength of sinusoidal mapping
  double zcenter; // Location of largest cells
  double w; // Radial width of domain in psi
  double psisep; // Separatrix psi value
};

/**
 * Create a new position map object. A position map is a function that maps 
 * uniform computational coordinates to non-uniform coordinates in the 
 * same coordinate space as computational coordinates. (e.g. uniform field
 * aligned -> non-uniform field aligned).
 *
 * @param pmap_info Comp. to phys. mapping input object (see definition above).
 * @param grid Position space grid.
 * @param local Local position range.
 * @param local_ext Local extended position range.
 * @param global Global position range.
 * @param use_gpu Whether to create a device copy of this new object.
 * @return New position map object.
 */
struct gkyl_position_map* gkyl_position_map_new(struct gkyl_position_map_inp pmap_info,
  struct gkyl_rect_grid grid, struct gkyl_range local, struct gkyl_range local_ext, 
  struct gkyl_range global, struct gkyl_range global_ext, struct gkyl_basis basis);

/**
 * Create a new position map object using the input structure.
 * @param inp Input structure (see definition of gkyl_position_map_inew_inp).
 * @return New position map object.
 */
struct gkyl_position_map*
gkyl_position_map_inew(struct gkyl_position_map_inew_inp inp);


/** Create a new null position map object. This is a position map that does nothing.
 *  All maps are identity maps.
 * @return New null position map object.
 */
struct gkyl_position_map*
gkyl_position_map_null_new();

/**
 * Set the position map object. Copy the non-uniform map array to the position map object.
 * 
 * @param gpm Position map object.
 * @param mc2nu Position map array.
 * 
 * @note This function is used to set the position map array in the position map object.
 */
void gkyl_position_map_set_mc2nu(struct gkyl_position_map* gpm, struct gkyl_array* mc2nu);

/**
 * Set the magnetic field array in the position map object. This is used to set the magnetic field
 * array in the position map object.
 *
 * @param gpm Position map object.
 * @param comm Communicator object.
 * @param bmag Magnetic field array.
 */
void
gkyl_position_map_set_bmag(struct gkyl_position_map* gpm, struct gkyl_comm* comm,
  struct gkyl_array* bmag);

/**
 * Set the function paramters for the map object.
 * 
 * @param gpm Position map object.
 * @param zcut half wavelength of sinusoidal mapping.
 * @param zcenter location of largest cells.
 * @param w radial width of domain in psi
 * @param psisep separatrix psi value.
 */
void
gkyl_position_map_set_compression(struct gkyl_position_map* gpm, double zcut, 
    double zcenter, double w, double psisep);

/**
 * Create a device (GPU) copy of the position map object, stored in
 * gpm->on_dev, for use inside device kernels with
 * gkyl_position_map_eval_mc2nu. Call it after the geometry has set the
 * mc2nu array (gkyl_position_map_set_mc2nu keeps the device mirror in sync
 * afterwards). Only available in GPU builds.
 *
 * @param gpm Position map object.
 */
void
gkyl_position_map_make_cu_dev(struct gkyl_position_map* gpm);

/**
 * Evaluate the position mapping at a specific computational (position) coordinate.
 * Callable on the host with the host object, or inside device kernels with the
 * device object (gpm->on_dev). On GPU builds do not call it on the host with
 * the device object or vice-versa.
 *
 * @param gpm Gkyl position map object.
 * @param x_comp Computational position coordinates.
 * @param x_fa Resulting non-uniform position coordinates.
 */
GKYL_CU_DH static inline void
gkyl_position_map_eval_mc2nu(const struct gkyl_position_map* gpm, const double *x_comp, double *x_fa)
{
  int cidx[GKYL_MAX_CDIM];
  for (int i=0; i<gpm->grid.ndim; i++) {
    int idxtemp = gpm->global.lower[i] + (int) floor((x_comp[i] - (gpm->grid.lower[i]) )/gpm->grid.dx[i]);
    idxtemp = GKYL_MAX2(GKYL_MIN2(idxtemp, gpm->local.upper[i]), gpm->local.lower[i]);
    cidx[i] = idxtemp;
  }
  long lidx = gkyl_range_idx(&gpm->local, cidx);
  const double *pmap_coeffs = (const double *) gkyl_array_cfetch(gpm->mc2nu, lidx);
  double cxc[GKYL_MAX_CDIM];
  double x_log[GKYL_MAX_CDIM];
  gkyl_rect_grid_cell_center(&gpm->grid, cidx, cxc);
  for (int i=0; i<gpm->grid.ndim; i++) {
    x_log[i] = (x_comp[i]-cxc[i])/(gpm->grid.dx[i]*0.5);
  }
  double xyz_fa[3];
  for (int i=0; i<3; i++) {
    xyz_fa[i] = gpm->basis.eval_expand(x_log, &pmap_coeffs[i*gpm->basis.num_basis]);
  }
  for (int i=0; i<gpm->grid.ndim; i++) {
    x_fa[i] = xyz_fa[i];
  }
  x_fa[gpm->grid.ndim-1] = xyz_fa[2];
}

/**
 * Evaluate the slope of the position mapping at a specific computational (position) coordinate.
 * 
 * @param gpm Gkyl position map object.
 * @param ix_map Index of the map to evaluate. Calls gpm->maps[index].
 * @param x Computational position coordinates.
 * @param dx Computational position increment to use for finite difference.
 * @param ix_comp Index in the geometry loop of which cell we are discussing
 * @param nrange Range of the computational coordinates.
 * @return Slope of the position mapping.
 */
double
gkyl_position_map_slope(const struct gkyl_position_map* gpm, int ix_map,
  double x, double dx, int ix_comp, const struct gkyl_range *nrange);

/**
 * Create a new pointer to the position map object.
 * Release it with gkyl_position_map_release.
 *
 * @param gpm Position map object.
 */
struct gkyl_position_map* gkyl_position_map_acquire(const struct gkyl_position_map* gpm);

/**
 * Optimize the position map object for constant B mapping.
 * 
 * @param gpm Position map object.
 * @param grid 3D Position space grid.
 * @param global 3D Global position range.
 */
void gkyl_position_map_optimize(struct gkyl_position_map* gpm, struct gkyl_rect_grid grid,
  struct gkyl_range global);


/**
 * Release pointer to (and eventually memory associated with)
 * the position map object. 
 *
 * @param Position map object.
 */
void gkyl_position_map_release(const struct gkyl_position_map *gpm);
