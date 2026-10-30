#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_flux_surf_avg.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_fem_poisson_perp_adiabatic gkyl_fem_poisson_perp_adiabatic;

// Input parameters for the perpendicular Poisson solver with an adiabatic response.
struct gkyl_fem_poisson_perp_adiabatic_inp {
  const struct gkyl_rect_grid *grid; // Configuration-space grid (2x or 3x).
  struct gkyl_basis basis; // Configuration-space basis (p=1).
  const struct gkyl_range *local; // Local range in which to solve.
  const struct gkyl_range *local_ext; // Local extended range.
  struct gkyl_basis basis_x; // 1D basis in x of the flux-surface averages.
  const struct gkyl_range *local_x; // 1D range in x of the flux-surface averages.
  const struct gkyl_range *local_x_ext; // 1D extended range in x of the flux-surface averages.
  struct gkyl_flux_surf_avg *fsa; // Flux-surface average (R) and extension (E) operators.
  struct gkyl_poisson_bc *bcs; // Perpendicular boundary conditions.
  struct gkyl_poisson_bias_line_list *bias_line_list; // Biased lines (NULL if none).
  struct gkyl_array *epsilon; // Polarization weight (permittivity) tensor.
  const struct gkyl_array *kJ; // Adiabatic response coefficient K.
  double lap_fac; // Factor alpha of the Laplacian term in the response (0 for none).
  // Optional in-place parallel projection (e.g. a smoothing along the field line) of every
  // right-hand side before the perpendicular solve (NULL for none), and its context.
  void (*fem_parproj_op)(struct gkyl_array *rhs, void *ctx);
  void *fem_parproj_ctx;
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater solving the perpendicular Poisson problem with an adiabatic response
 *   - nabla_perp . (epsilon nabla_perp phi) + M (phi - <phi>) = rho,
 *   M = K - alpha nabla_perp . (epsilon nabla_perp),
 * where <phi> is the flux-surface average and alpha = lap_fac (the Laplacian term in the
 * response arises with Pade FLR effects). With H = -(1+alpha) nabla_perp . (epsilon nabla_perp)
 * + K, R the average and E the extension of a function of x to a field constant along y and z,
 * the Woodbury identity gives
 *   phi = phi1 + H^{-1} K E psi + alpha E psi,  H phi1 = rho,  (I - G) psi = R phi1,
 * with G = R H^{-1} K E. Nonzero Dirichlet values add -beta R phi_0 to R phi1 and -beta phi_0 to
 * phi, with beta = alpha/(1+alpha) and phi_0 the solution with rho = 0.
 * The m x m matrix I - G (m = number of zonal degrees of freedom) is assembled with m Helmholtz
 * solves and inverted once here, so each solve costs two Helmholtz solves, one average and an
 * m x m matrix-vector product. Free using gkyl_fem_poisson_perp_adiabatic_release.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_fem_poisson_perp_adiabatic *gkyl_fem_poisson_perp_adiabatic_inew(
  const struct gkyl_fem_poisson_perp_adiabatic_inp *inp
);

/**
 * Solve for phi given rho. rho is overwritten by its projection if fem_parproj_op is set.
 *
 * @param up Updater object.
 * @param rho Right-hand side (DG, on local).
 * @param phi Output potential (DG, on local).
 */
void gkyl_fem_poisson_perp_adiabatic_advance(
  struct gkyl_fem_poisson_perp_adiabatic *up, struct gkyl_array *rho, struct gkyl_array *phi
);

/**
 * 1-norm condition number of the zonal matrix I - G.
 *
 * @param up Updater object.
 * @return Condition number.
 */
double gkyl_fem_poisson_perp_adiabatic_cond(const struct gkyl_fem_poisson_perp_adiabatic *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_perp_adiabatic_release(struct gkyl_fem_poisson_perp_adiabatic *up);
