#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_comm.h>

// Object type
typedef struct gkyl_bc_updown_tok_core gkyl_bc_updown_tok_core;

// Inputs to the up-down tokamak core BC updater.
struct gkyl_bc_updown_tok_core_inp {
  int dir; // Direction in which to apply the BC (must be 0, i.e. x).
  enum gkyl_edge_loc edge; // Edge at which to apply the BC (must be GKYL_LOWER_EDGE).
  int cdim; // Number of configuration space dimensions (>=2, z is cdim-1).
  const struct gkyl_rect_grid *grid; // Phase-space grid (to check z symmetry).
  const struct gkyl_basis *basis; // Phase-space basis (a device pointer if use_gpu=true).
  const struct gkyl_range *skin_r; // Local phase-space skin range.
  const struct gkyl_range *ghost_r; // Local phase-space ghost range.
  const struct gkyl_rect_decomp *decomp; // Conf-space decomposition (only z may be cut).
  struct gkyl_comm *comm; // Communicator. We should not pass it to updaters but for now we do that for simplicity.
  bool avg_y; // Average the donor over y (only meaningful for cdim=3).
  bool use_gpu; // Whether to run on GPU.
};

/**
 * Create a new updater that fills the ghost cells of the lower radial boundary with
 * the mirror image of the skin cells of an up-down symmetric equilibrium:
 *   f_ghost(x_g, y, z, vpar, mu) = < f_skin(x_s, y', -z, vpar, mu) >_y'
 * This stands in for the absent core: particles leaving through the boundary
 * cross the core and come back at -z with the same vpar and mu.
 *
 * Each rank packs its skin, mirrored in z (within the local z extent) and
 * averaged over y, into a buffer, sends it to the rank owning the
 * mirrored z extent, and unpacks the buffer it receives into its ghost cells
 * (broadcasting in y). The ghost coefficients are also reflected in x, so the
 * ghost trace at the face equals the donor trace and the upwinded fluxes at
 * (z,vpar) and (-z,vpar) cancel exactly (the radial magnetic drifts are odd
 * in z and even in vpar). Since the z decomposition must be uniform
 * (Nz % cuts_z == 0), reversing z within the local skin is the global mirror.
 *
 * Valid for up-down symmetric equilibria with a z grid symmetric about z=0
 * (a nonuniform z map must be up-down symmetric too).
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_bc_updown_tok_core*
gkyl_bc_updown_tok_core_new(const struct gkyl_bc_updown_tok_core_inp *inp);

/**
 * Apply the up-down tokamak core BC. Must be called by all ranks.
 *
 * @param up BC updater.
 * @param distf Distribution function to apply the BC to.
 */
void gkyl_bc_updown_tok_core_advance(struct gkyl_bc_updown_tok_core *up, struct gkyl_array *distf);

/**
 * Free memory associated with bc_updown_tok_core updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_updown_tok_core_release(struct gkyl_bc_updown_tok_core *up);
