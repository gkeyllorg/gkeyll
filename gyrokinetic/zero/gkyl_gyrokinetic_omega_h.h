#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

struct gkyl_gyrokinetic_omega_h;

/**
 * Cache omegaH_geo_fac squared at Gauss-Legendre nodes and the cell center.
 * omegaH_geo_fac is the geometry/polarization factor, excluding |q|/sqrt(m).
 *
 * In 3x, let (x,y,z) denote the computational coordinates. The formula applies
 * to both MAPC2P and TOKAMAK geometry;
 *
 *   omegaH_geo_fac = (C/(J B)) k_z_max
 *     / sqrt(k_x_min^2 epsilon_xx + 2 k_x_min k_y_min epsilon_xy
 *            + k_y_min^2 epsilon_yy).
 *
 * Here J is the configuration-space Jacobian, B = |B|, C = cmag, and epsilon_ij
 * are the stored polarization coefficients, including J. For field-aligned
 * coordinates, C = J B/sqrt(g_zz), where g_zz is the covariant metric component,
 * so (C/(J B)) k_z_max estimates the parallel wavenumber in inverse metres.
 * The computational wavenumbers are
 *
 *   k_z_max = pi (p+1)/Delta_z, k_x_min = pi/L_x, k_y_min = pi/L_y,
 *
 * with configuration-space polynomial order p, cell width Delta_z, and domain
 * widths L_x and L_y. In 2x omit the y terms; in 1x epsilon already includes
 * the prescribed k_perp^2. The implementation uses B, not B_parallel^*.
 *
 * The supplied modal factor is constructed using DG products and a projected
 * inverse square root. Its sampled values are squared directly; subsequent
 * products with density are not projected back into the DG basis.
 *
 * @param basis Configuration-space basis on the host.
 * @param range Interior configuration range, indexing omegaH_geo_fac and later density arrays.
 * @param omegaH_geo_fac Scalar modal geometry/polarization factor; device array
 * when use_gpu is true. Not retained.
 * @param use_gpu Whether subsequent reductions run on the GPU.
 * @return New updater; release with gkyl_gyrokinetic_omega_h_release.
 */
struct gkyl_gyrokinetic_omega_h *gkyl_gyrokinetic_omega_h_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  const struct gkyl_array *omegaH_geo_fac, bool use_gpu
);

/**
 * Return max(max(Jn,0)*omegaH_geo_fac^2) over the cached range and evaluation points.
 * Multiply its square root by |q|/sqrt(m) to obtain the species frequency
 * estimate. This is a sampled local estimate, not a discrete spectral bound.
 *
 * @param up Updater containing cached geometry and basis evaluations.
 * @param density Scalar modal Jacobian-weighted density Jn, with the constructor's
 * range indexing and memory location. Nonpositive sampled density contributes zero.
 * @return Local maximum on the host; zero for an empty range or nonpositive density.
 */
double gkyl_gyrokinetic_omega_h_advance(
  struct gkyl_gyrokinetic_omega_h *up, const struct gkyl_array *density
);

/**
 * Release the cached geometry, basis evaluations and reduction storage.
 *
 * @param up Updater to release; NULL is accepted.
 */
void gkyl_gyrokinetic_omega_h_release(struct gkyl_gyrokinetic_omega_h *up);
