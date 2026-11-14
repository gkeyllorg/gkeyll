#pragma once

#include <gkyl_gk_fluid_moments.h>
#include <gkyl_dg_bin_ops.h>

struct gkyl_gk_fluid_moments {
  struct gkyl_basis basis; // Configuration-space basis.
  int num_mom_out; // Number of output moments.
  int num_mom_in; // Number of conserved fields in the fluid state.
  int component; // Conservative-state component used for M0, M1, or M2.
  long num_cells; // Batch size of the weak division solves.
  double mass, charge; // Species mass and charge.
  double volnorm_fac; // Converts the constant DG coefficient to a cell integral.
  double rmass_fac; // State-to-moment factor: 1/m for M0 and M1, 2/m for M2.
  struct gkyl_array *bmag; // Acquired magnetic-field magnitude.
  struct gkyl_array *upar; // Parallel-velocity buffer, also reused for conversion and scaling.
  struct gkyl_array *ppar; // Jacobian-weighted parallel-pressure buffer and intermediate results.
  struct gkyl_array *pperp; // Jacobian-weighted perpendicular-pressure buffer.
  struct gkyl_array *result; // DG coefficients of calculated moments, before output or integration.
  gkyl_dg_bin_op_mem *mem; // Preallocated workspace for weak-DG division.
  // Methods selected at initialization:
  // Calculate the requested moments into result.
  void (*calc_func)(const struct gkyl_gk_fluid_moments *, const struct gkyl_range *, const struct gkyl_array *, const struct gkyl_array *);
  // Copy DG moments to output or compute their cell integrals.
  void (*final_func)(const struct gkyl_gk_fluid_moments *, const struct gkyl_range *, struct gkyl_array *);
};

#ifdef GKYL_HAVE_CUDA
/** Extract constant coefficients as cell integrals on the GPU.
 * @param up Fluid updater with DG results in its buffer array.
 * @param range Configuration subrange.
 * @param out Output cell integrals.
 */
void gkyl_gk_fluid_moments_final_integrate_cu(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range, struct gkyl_array *out
);
#endif
