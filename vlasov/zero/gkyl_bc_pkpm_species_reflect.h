#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_bc_pkpm_species_reflect gkyl_bc_pkpm_species_reflect;

// Input packaged as a struct
struct gkyl_bc_pkpm_species_reflect_inp {
  int dir; // Configuration-space direction of the wall normal.
  enum gkyl_edge_loc edge; // Lower or upper wall.
  int cdim; // Configuration-space dimensions (the normalized parallel velocity is phase dimension cdim).
  const struct gkyl_rect_grid
    *phase_grid; // Phase-space grid (its velocity extent must be symmetric about w = 0).
  const struct gkyl_basis *phase_basis; // Phase-space basis (host object, for its size).
  const struct gkyl_basis *
    phase_basis_on_dev; // Phase-space basis whose flip functions are called (device object when use_gpu).
  int conf_num_basis; // Number of configuration-space basis functions (layout of bvar).
  int num_blocks; // Number of blocks of num_basis components in the distribution array.
  const struct gkyl_range *skin_r; // Phase-space skin range of the wall (interior velocity cells).
  const struct gkyl_range *ghost_r; // Phase-space ghost range of the wall (interior velocity cells).
  const struct gkyl_range *phase_range_ext; // Extended phase-space range indexing the distributions.
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range indexing bvar.
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater applying a specular (reflecting) wall to the normalized PKPM
 * distributions [F, G, M0_perp (x, y, z), N_perp (x, y, z)]. The ghost state is the mirror
 * image f'(v) = f(R v), R = I - 2 n n, expressed in the ghost cell's own field-aligned
 * coordinates. If the ghost field direction is b' = s R b (s = +1 for a PEC wall, whose normal B
 * is odd; s = -1 for a symmetry wall, whose normal B is even), the image is exact in the
 * first-harmonic model for any orientation of b:
 *   F'(w) = F(s w), G'(w) = G(s w), M0_perp'(w) = R M0_perp(s w), N_perp'(w) = R N_perp(s w),
 * each also mirrored in the wall direction. The sign follows the parity of the field across the
 * wall, per ghost cell from the cell averages of the ghost and skin b at every application: w is
 * flipped when the normal field is even (b'_n b_n > 0) or the image is the reversed mirror
 * (b' . R b < 0, which decides when b_n = 0), and kept otherwise. For PEC and symmetry walls this
 * is the exact image. For a copied field (b' = b) at oblique b the exact image is not gyrotropic
 * about b'; flipping w then keeps the mass flux through the wall exactly zero; the energy and
 * tangential momentum fluxes are not zero in general.
 *
 * The image at s w is taken cell for cell from the mirrored velocity cell, so the velocity grid
 * must be symmetric about w = 0.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_bc_pkpm_species_reflect *gkyl_bc_pkpm_species_reflect_inew(
  const struct gkyl_bc_pkpm_species_reflect_inp *inp
);

/**
 * Create the updater on NV-GPU. See the inew method for documentation.
 */
struct gkyl_bc_pkpm_species_reflect *gkyl_bc_pkpm_species_reflect_cu_dev_inew(
  const struct gkyl_bc_pkpm_species_reflect_inp *inp
);

/**
 * Fill the ghost cells of the wall with the mirror image of the skin cells.
 *
 * @param up Updater.
 * @param bvar Magnetic field unit vector [bx, by, bz] in the skin and ghost cells (after the
 *             field boundary conditions).
 * @param f Distributions; the ghost cells of the wall are overwritten.
 */
void gkyl_bc_pkpm_species_reflect_advance(
  const struct gkyl_bc_pkpm_species_reflect *up, const struct gkyl_array *bvar, struct gkyl_array *f
);

/**
 * Host-side wrapper of the advance method on NV-GPU.
 */
void gkyl_bc_pkpm_species_reflect_advance_cu(
  const struct gkyl_bc_pkpm_species_reflect *up, const struct gkyl_array *bvar, struct gkyl_array *f
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_bc_pkpm_species_reflect_release(struct gkyl_bc_pkpm_species_reflect *up);
