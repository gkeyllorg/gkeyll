#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

typedef struct gkyl_fem_proj gkyl_fem_proj;

// No slope or flux conditions are needed for a mass projection.
enum gkyl_fem_proj_bc_type {
  GKYL_FEM_PROJ_NONE = 0,
  GKYL_FEM_PROJ_PERIODIC,
  GKYL_FEM_PROJ_DIRICHLET,
  GKYL_FEM_PROJ_DIRICHLET_VARYING
};

struct gkyl_fem_proj_bc {
  enum gkyl_fem_proj_bc_type lo_type[GKYL_MAX_CDIM], up_type[GKYL_MAX_CDIM];
  double lo_value[GKYL_MAX_CDIM], up_value[GKYL_MAX_CDIM];
};

// Constant Dirichlet data on a grid-aligned plane spanning the solve range.
struct gkyl_fem_proj_bias_plane {
  int dir; // Direction normal to the plane.
  double loc; // Physical coordinate in the supplied grid.
  double val; // Prescribed solution on the plane.
};

struct gkyl_fem_proj_bias_plane_list {
  int num_bias_plane; // Number of planes supplied.
  const struct gkyl_fem_proj_bias_plane *bp; // Plane data, copied during construction.
};

struct gkyl_fem_proj_inp {
  const struct gkyl_rect_grid *grid;
  const struct gkyl_range *solve_range; // May be a subrange of an extended range.
  const struct gkyl_basis *basis;
  struct gkyl_fem_proj_bc bcs;
  // Scalar DG weights in the input basis/layout; NULL independently means unity.
  // Weights must remain unchanged until release. weight_left must be positive.
  const struct gkyl_array *weight_left, *weight_right;
  // Optional fixed biases. Requires solve indices within grid cells [1,cells[d]].
  // Planes outside the solve range are ignored. Intersecting prescriptions,
  // including varying boundary data and periodic copies, must agree.
  const struct gkyl_fem_proj_bias_plane_list *bias_plane_list;
  bool use_gpu;
};

/**
 * Project a DG field onto the FEM basis solving:
 *   integral(v*wL*u) = integral(v*wR*f).
 * Equal weights give weighted L2 projection. The output is continuous in every
 * coordinate and represented as DG coefficients.
 *
 * @param inp Constructor inputs, including boundary conditions and fixed weights.
 */
gkyl_fem_proj *gkyl_fem_proj_inew(const struct gkyl_fem_proj_inp *inp);

/**
 * Assemble a new RHS vector without changing the matrix.
 * 
 * @param up Projection updater.
 * @param rhsin Scalar DG input, defined over the parent of the solve range.
 * @param boundary_values DG field in the same layout, evaluated on skin-cell
 * boundary nodes for DIRICHLET_VARYING. NULL if no varying BC is used. All
 * prescribed values at shared nodes must agree, including periodic copies.
 */
void gkyl_fem_proj_set_rhs(
  gkyl_fem_proj *up, const struct gkyl_array *rhsin, const struct gkyl_array *boundary_values
);

/**
 * Solve and reconstruct the DG output, assigning only the solve range.
 *
 * @param up Projection updater with an assembled RHS.
 * @param out Scalar DG output in the input layout. May alias the RHS input.
 */
void gkyl_fem_proj_solve(gkyl_fem_proj *up, struct gkyl_array *out);

/**
 * Free memory allocated by this updater.
 * 
 * @param up Updater.
 */
void gkyl_fem_proj_release(gkyl_fem_proj *up);
