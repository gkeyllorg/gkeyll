#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

#include <stdbool.h>

// Object type.
typedef struct gkyl_gk_dens_mu_cells gkyl_gk_dens_mu_cells;

/**
 * Create a new updater that computes the density of a gyrokinetic species
 * resolved in the magnetic moment mu, i.e. the density of the particles in
 * each mu cell j,
 *
 *   N_j(x) = (2*pi/m) int_{mu cell j} int f(x,vpar,mu) dvpar dmu,
 *
 * so that sum_j N_j is the M0 moment of f. These are the weights that a
 * function of (x,mu) that is constant in each mu cell is integrated against
 * when multiplied by f and integrated over velocity space.
 *
 * @param phase_grid Phase space grid.
 * @param conf_basis Configuration space basis.
 * @param phase_basis Phase space basis.
 * @param mass Species mass.
 * @param use_gpu Whether to run on the GPU.
 * @return New updater pointer.
 */
struct gkyl_gk_dens_mu_cells *gkyl_gk_dens_mu_cells_new(
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, double mass, bool use_gpu
);

/**
 * Compute the mu-resolved density.
 *
 * @param up Updater.
 * @param phase_range Phase space range.
 * @param conf_range Configuration space range.
 * @param fin Distribution function.
 * @param out Density in each mu cell (num_mu*num_basis components, the
 *            coefficients for the j-th mu cell starting at j*num_basis, with
 *            num_mu the number of mu cells and num_basis the number of
 *            configuration space basis functions).
 */
void gkyl_gk_dens_mu_cells_advance(
  const struct gkyl_gk_dens_mu_cells *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *fin, struct gkyl_array *out
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_dens_mu_cells_release(struct gkyl_gk_dens_mu_cells *up);
