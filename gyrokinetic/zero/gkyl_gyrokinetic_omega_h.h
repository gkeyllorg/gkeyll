#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

struct gkyl_gyrokinetic_omega_h;

/**
 * Cache G squared at Gauss-Legendre nodes and the cell center for the local
 * omega_H estimate. G is the geometry/polarization factor, excluding q/sqrt(m).
 * Products are evaluated pointwise without projecting them back to a DG basis.
 *
 * @param basis Configuration-space basis on the host.
 * @param range Interior configuration range, indexing geo_factor and later density arrays.
 * @param geo_factor Scalar modal G; device array when use_gpu is true. Not retained.
 * @param use_gpu Whether subsequent reductions run on the GPU.
 * @return New updater; release with gkyl_gyrokinetic_omega_h_release.
 */
struct gkyl_gyrokinetic_omega_h *gkyl_gyrokinetic_omega_h_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  const struct gkyl_array *geo_factor, bool use_gpu
);

/**
 * Return max(max(Jn,0)*G^2) over the cached range and evaluation points.
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
