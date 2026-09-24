#pragma once

#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_gyrokinetic_pol_density gkyl_gyrokinetic_pol_density;

/**
 * Create a new updater to compute the total polarization charge density
 *   sum_s d/dx_i jacobgeo*(n_s*m_s/B^2) g^{ij} d/dx_j phi,
 * where i,j \in {1,2}.
 *
 * Currently assumes phi is represented with a polynomial tensor basis of one
 * higher order, and the calculation is done in a single volume term without
 * recovery. It is only meant for ICs.
 *
 * @param cbasis Configuration-space basis.
 * @param cgrid Phase-space grid.
 * @param use_gpu bool to determine if on GPU.
 * @return New polarization density updater pointer.
 */
struct gkyl_gyrokinetic_pol_density*
gkyl_gyrokinetic_pol_density_new(struct gkyl_basis cbasis, struct gkyl_rect_grid cgrid,
  bool use_gpu);

/**
 * Run the polarization density updater in the indicated range.
 *
 * @param up Polarization density updater.
 * @param conf_rng Config-space range.
 * @param pol_weight Polarization density weight (factor inside differential op).
 * @param phi Electrostatic potential (represented with a p+1 tensor basis).
 * @param npol Polarization density.
 */
void
gkyl_gyrokinetic_pol_density_advance(gkyl_gyrokinetic_pol_density* up,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT pol_weight,
  const struct gkyl_array *GKYL_RESTRICT phi, struct gkyl_array *GKYL_RESTRICT npol);

/**
 * Release the memory associated with this polarization density updater.
 *
 * @param up Positivity shift updater.
 */
void
gkyl_gyrokinetic_pol_density_release(gkyl_gyrokinetic_pol_density* up);
