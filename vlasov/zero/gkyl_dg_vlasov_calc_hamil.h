#pragma once

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_velocity_map.h>

/**
 * Construct the velocity-space Hamiltonian on the velocity-space grid
 * given the input model_id, either GKYL_MODEL_DEFAULT (H = v^2/2) or
 * GKYL_MODEL_SR (H = sqrt(1 + u^2)).
 * Note that these Hamiltonians do not have factors of mass in them.
 *
 * @param vel_grid Velocity-space grid.
 * @param vel_basis Velocity-space basis.
 * @param vel_range Velocity-space range.
 * @param model_id Model indicator for which velocity-space Hamiltonian to project.
 * @param vel_map Velocity-space mapping object for computing Hamiltonian on mapped grids. NULL => uniform velocity grid.
 * @param hamil Output array storing the velocity-space Hamiltonian.
 * @param hamil_inv Output array storing the inverse velocity-space Hamiltonian (1/gamma for relativity).
 * @param use_gpu bool to determine if on GPU
 */
void gkyl_dg_vlasov_calc_hamil(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_model_id model_id,
  const struct gkyl_vlasov_velocity_map *vel_map, struct gkyl_array *hamil,
  struct gkyl_array *hamil_inv, bool use_gpu
);

/**
 * Construct the velocity-space Hamiltonian on the velocity-space grid on
 * NV-GPU. See new() method for documentation.
 */
void gkyl_dg_vlasov_calc_hamil_cu(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_model_id model_id,
  const struct gkyl_vlasov_velocity_map *vel_map, struct gkyl_array *hamil,
  struct gkyl_array *hamil_inv
);