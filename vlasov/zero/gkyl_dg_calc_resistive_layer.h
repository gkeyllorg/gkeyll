#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_resistive_layer gkyl_dg_calc_resistive_layer;

// Map from computational to physical coordinates.
typedef void (*dg_calc_resistive_layer_c2p_t)(const double *xcomp, double *xphys, void *ctx);

// Input packaged as a struct
struct gkyl_dg_calc_resistive_layer_inp {
  const struct gkyl_rect_grid *grid; // Configuration-space grid.
  const struct gkyl_basis *basis; // Configuration-space basis.
  const struct gkyl_range *range; // Range on which sigma is sampled.
  evalf_t eval; // Resistive layer sigma(t, x) >= 0, the damping rate of E (sampled at t = 0).
  void *ctx; // Context for eval.
  dg_calc_resistive_layer_c2p_t c2p_func; // Computational to physical coordinates (0: identity).
  void *c2p_func_ctx; // Context for c2p_func.
  struct gkyl_array *sigma; // Output: sigma at the (p+1)^cdim quadrature points of each cell.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create an updater for a resistive layer, dE/dt = -sigma E. The constructor samples sigma at the
 * (p+1)^cdim Gauss-Legendre points of each cell of the range (row-major order, the last
 * direction fastest, as in the projection updaters) into inp->sigma. The increment of an expansion
 * E is the quadrature projection
 *
 *   (sigma E)_k = sum_q w_q b_k(x_q) sigma(x_q) E(x_q),
 *
 * so sum_k E_k (sigma E)_k = sum_q w_q sigma(x_q) |E(x_q)|^2: for any p the layer only removes
 * energy when sigma >= 0 at the points. At p = 1 this is the exact weak product of E and the p=1
 * interpolant of sigma, and the values are those the nodal fluid-EM source solve uses.
 *
 * @param inp Input parameters defined in gkyl_dg_calc_resistive_layer_inp struct.
 * @return New updater pointer.
 */
struct gkyl_dg_calc_resistive_layer *gkyl_dg_calc_resistive_layer_inew(
  const struct gkyl_dg_calc_resistive_layer_inp *inp
);

/**
 * Accumulate -(sigma E) onto the electric field rows of a right-hand side.
 *
 * @param up Updater.
 * @param range Range to update.
 * @param em Input fields: E in the first three expansions.
 * @param rhs Right-hand side: its first three expansions are incremented.
 */
void gkyl_dg_calc_resistive_layer_advance(
  const struct gkyl_dg_calc_resistive_layer *up, const struct gkyl_range *range,
  const struct gkyl_array *em, struct gkyl_array *rhs
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_resistive_layer_release(struct gkyl_dg_calc_resistive_layer *up);
