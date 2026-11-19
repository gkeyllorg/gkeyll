#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_pkpm_mom_diffusion gkyl_dg_pkpm_mom_diffusion;

// Input packaged as a struct
struct gkyl_dg_pkpm_mom_diffusion_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid (computational).
  const struct gkyl_basis *conf_basis; // Configuration-space basis (p=1 serendipity).
  const struct gkyl_vlasov_position_map *pos_map; // Position map (REQUIRED; the identity map gives the uniform-grid operator).
  int order; // Order of the diffusion: 2 (grad^2) or 4 (grad^4).
  double D[GKYL_MAX_CDIM]; // Diffusion coefficient per configuration direction (0 switches a direction off).
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater for the momentum diffusion of the PKPM fluid on the position map,
 *
 *   d_t(J M) = J sum_d (-1)^(order/2+1) D_d d^order M / dx_d^order,
 *
 * with J the per-cell constant Jacobian of the position map and M = (J M)/J the physical momentum
 * (the momentum blocks of the J-weighted fluid state [rho, M_x, M_y, M_z, E]). Recovery DG in
 * physical coordinates: the face values of the cubic recovered from the
 * two cells sharing a face (their physical widths J_d dxi_d), shared by both cells, so the momentum
 * is conserved; on the identity map the operator is that of gkyl_dg_updater_diffusion_fluid.
 * p=1 serendipity only.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_dg_pkpm_mom_diffusion *gkyl_dg_pkpm_mom_diffusion_inew(
  const struct gkyl_dg_pkpm_mom_diffusion_inp *inp
);

/**
 * Create the updater on NV-GPU. See the inew method for documentation.
 */
struct gkyl_dg_pkpm_mom_diffusion *gkyl_dg_pkpm_mom_diffusion_cu_dev_inew(
  const struct gkyl_dg_pkpm_mom_diffusion_inp *inp
);

/**
 * Accumulate the momentum diffusion into the fluid right-hand side and its time-step rate into
 * cflrate. Reads the ghost cells of fluid in every direction with a nonzero coefficient.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range of the cells to update.
 * @param fluid J-weighted fluid state [rho, M_x, M_y, M_z, E] (ghost cells filled).
 * @param cflrate Configuration-space time-step rates, incremented.
 * @param rhs Fluid right-hand side, incremented in the momentum blocks.
 */
void gkyl_dg_pkpm_mom_diffusion_advance(
  const struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Host-side wrapper of the advance method on NV-GPU.
 */
void gkyl_dg_pkpm_mom_diffusion_advance_cu(
  const struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_mom_diffusion_release(struct gkyl_dg_pkpm_mom_diffusion *up);
