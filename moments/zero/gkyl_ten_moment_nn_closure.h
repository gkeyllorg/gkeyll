#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

struct gkyl_ten_moment_nn_closure_inp
{
  const struct gkyl_rect_grid *grid; // Grid on which to solve equations.
  int poly_order; // Polynomial order of learned DG coefficients.
  double k0; // Damping coefficient.
  struct gkyl_kann_net* ann; // Neural network architecture.
};

// Object type.
typedef struct gkyl_ten_moment_nn_closure gkyl_ten_moment_nn_closure;

// Cached magnetic-field-aligned geometry produced by the closure's geometry
// pass and consumed when constructing the heat-flux divergence after
// inference. Computed once per cell so the gather and scatter phases share it.
struct gkyl_ten_moment_nn_closure_geom {
  double local_mag[3];    // Unit magnetic field direction b = B / |B|.
  double local_mag_dx[3]; // d(b) / dx.
  double local_mag_dy[3]; // d(b) / dy (2D only; zero otherwise).
  double B_avg[3];        // Cell-averaged magnetic field.
  double rho_avg;         // Cell-averaged mass density.
  double p_avg[6];        // Cell-averaged pressure tensor (xx, xy, xz, yy, yz, zz).
};

/**
 * Number of neural-network input features for the given closure configuration.
 *
 * @param nnclosure Neural network closure updater object.
 * @return Input feature count (6 for 1D-p1, 9 for 1D-p2, 12 for 2D-p1).
 */
int
gkyl_ten_moment_nn_closure_n_in(const gkyl_ten_moment_nn_closure *nnclosure);

/**
 * Number of neural-network output predictions for the given closure configuration.
 *
 * @param nnclosure Neural network closure updater object.
 * @return Output count (4 for 1D-p1, 6 for 1D-p2, 8 for 2D-p1).
 */
int
gkyl_ten_moment_nn_closure_n_out(const gkyl_ten_moment_nn_closure *nnclosure);

/**
 * Compute the neural-network input feature vector (of length n_in) and cache
 * the magnetic-field-aligned geometry for a single cell from its fluid/EM
 * stencil. The geometry is computed exactly once here and reused when
 * constructing the heat-flux divergence after inference. 
 *
 * @param nnclosure Neural network closure updater object.
 * @param fluid_d Stencil of fluid variable pointers (length set by dim/poly_order).
 * @param em_tot_d Stencil of total EM field pointers (length set by dim/poly_order).
 * @param input_data Output: input feature vector (caller-allocated, length n_in).
 * @param geom Output: cached aligned geometry for this cell.
 */
void
gkyl_ten_moment_nn_closure_geom_calc(const gkyl_ten_moment_nn_closure *nnclosure, const double *fluid_d[], const double *em_tot_d[], float *input_data,
  struct gkyl_ten_moment_nn_closure_geom *geom);

/**
 * Construct the magnetized heat-flux tensor from a (network) prediction and
 * the cached aligned geometry, and write the ten-moment closure source into
 * rhs: d(P_ij)/dt = ... - d(q_ijk)/dx_k. 
 *
 * @param nnclosure Neural network closure updater object.
 * @param geom Cached aligned geometry for the cell (b = B/|B|, its gradients, B_avg).
 * @param output_data_predicted Predicted outputs (q and its gradients), length n_out.
 * @param rhs Output: closure source for the conserved variables (length 10).
 */
void
gkyl_ten_moment_nn_closure_construct(const gkyl_ten_moment_nn_closure *nnclosure, const struct gkyl_ten_moment_nn_closure_geom *geom,
  const float *output_data_predicted, double *rhs);

/**
 * Create a new updater to update the pressure tensor in the ten moment equation system using a machine-learned magnetized closure trained on PKPM simulations,
 * i.e. where the heat flux tensor q_ijk and its derivatives are calculated from predictions for q_par and q_perp on the basis of rho, p_par and p_perp.
 * Returns RHS for accumulation in a forward-Euler method.
 *
 * @param inp Input parameters for updater.
 * @return Pointer to updater.
 */
gkyl_ten_moment_nn_closure*
gkyl_ten_moment_nn_closure_new(struct gkyl_ten_moment_nn_closure_inp inp);

/**
 * Compute the right-hand-side contribution to the ten moment equation system from a machine-learned magnetized closure trained on PKPM simulations.
 * The update_rng MUST be a sub-range of the range on which the array is defined. That is, it must be either the same range as the arary range, or one
 * created using the gkyl_sub_range_init method.
 *
 * @param nnclosure Neural network closure updater object.
 * @param heat_flux_rng Range on which to compute the heat flux tensor (at cell nodes).
 * @param update_rng Range on which to compute update.
 * @param fluid Input array of fluid variables.
 * @param em_tot Total electromagnetic field variables (plasma + external).
 * @param heat_flux Array for storing intermediate computation of heat flux tensor (at cell nodes).
 * @param rhs Right-hand-side output.
 */
void
gkyl_ten_moment_nn_closure_advance(const gkyl_ten_moment_nn_closure *nnclosure, const struct gkyl_range *heat_flux_rng, const struct gkyl_range *update_rng,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot, struct gkyl_array *heat_flux, struct gkyl_array *rhs);

/**
 * Delete updater.
 *
 * @param nnclosure Updater to delete.
 */
 void
 gkyl_ten_moment_nn_closure_release(gkyl_ten_moment_nn_closure *nnclosure);