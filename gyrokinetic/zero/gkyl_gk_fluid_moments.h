#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_gk_fluid_moments;

struct gkyl_gk_fluid_moments_inp {
  const struct gkyl_basis *basis; // Configuration-space basis.
  const struct gkyl_rect_grid *grid; // Configuration-space grid.
  const struct gkyl_range *range, *range_ext; // Interior and extended ranges.
  const struct gkyl_array *bmag; // Magnetic magnitude (acquired by updater).
  double mass, charge; // Species mass and charge.
  int num_moments; // Number of conserved fields in the fluid state.
  enum gkyl_distribution_moments moment; // Moment to calculate.
  bool integrated; // Return cell integrals instead of DG coefficients.
  bool use_gpu; // Allocate scratch on the GPU.
};

/**
 * Create weak-DG conversions and Maxwellian closures for J*(rho,Mpar,E,pperp/B).
 * Always includes parallel and perpendicular velocity moments (2v).
 *
 * @param inp Configuration and geometry; arrays must reside on the selected device.
 */
struct gkyl_gk_fluid_moments *gkyl_gk_fluid_moments_inew(const struct gkyl_gk_fluid_moments_inp *inp
);

/**
 * Return the number of output moments.
 * @param up Fluid moment updater.
 */
int gkyl_gk_fluid_moments_num_mom(const struct gkyl_gk_fluid_moments *up);

/**
 * Calculate moments, retaining J on densities but not on primitive u and T/m.
 *
 * @param up Fluid moment updater.
 * @param range Configuration range used at initialization.
 * @param state Conservative state J*(rho,Mpar,E,pperp/B).
 * @param phi Potential, required only for Hamiltonian moments.
 * @param out Output DG moments or cell integrals.
 */
void gkyl_gk_fluid_moments_advance(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi, struct gkyl_array *out
);

/**
 * Convert primitive fields using weighted-first weak products.
 *
 * @param up Fluid moment updater.
 * @param range Configuration range used at initialization.
 * @param prim Fields (n,upar,Tpar/m,Tperp/m), without J.
 * @param jacobgeo Configuration-space Jacobian.
 * @param state Output J*(rho,Mpar,E,pperp/B); must not alias prim.
 */
void gkyl_gk_fluid_moments_from_prim(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *prim, const struct gkyl_array *jacobgeo, struct gkyl_array *state
);

/**
 * Scale all conserved fields by target/current density, then set rho exactly.
 *
 * @param up Fluid moment updater.
 * @param range Configuration range used at initialization.
 * @param density Target J*n, not modified.
 * @param state Conservative state to rescale in place.
 */
void gkyl_gk_fluid_moments_set_density(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *density, struct gkyl_array *state
);

/**
 * Release updater and scratch.
 *
 * @param up Fluid moment updater.
 */
void gkyl_gk_fluid_moments_release(struct gkyl_gk_fluid_moments *up);
