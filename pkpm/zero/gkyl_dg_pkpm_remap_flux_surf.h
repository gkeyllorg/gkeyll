#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_pkpm_remap_flux_surf gkyl_dg_pkpm_remap_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_remap_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis
    *phase_basis; // Phase-space basis functions (p=1 hybrid, one velocity dimension).
  const struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space position map (per-direction Jacobians of the cells).
  int num_face_nodes; // Gauss-Legendre nodes per transverse direction of a face (>= 3).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the configuration-space face fluxes of PKPM: the kinetic flux of
 * the eight blocks [F, G, M0_perp, N_perp], the fluid flux [F_rho, F_x, F_y, F_z, F_E] and the
 * time-step rates of both, in one pass (physical_jump_flux.tex, section "The exact remap";
 * PKPM_REMAP_FLUX_PLAN.md).
 *
 * At every transverse Gauss node of a face, rho, u, xi and b are averaged across the face and the
 * two sides' blocks are mapped exactly to that face frame, on the supermesh of both sides' mapped
 * velocity grids and of the kinks of the HLL weights. The HLL flux is formed there; it is mapped
 * back and projected exactly onto each side's velocity cells, so the two sides of a face carry
 * different kinetic fluxes. The fluid flux is the face-frame moment of the same flux over the
 * union of the two grids.
 *
 * The result is the blend (1 - blend) remap + blend label-paired flux, where the label-paired
 * flux pairs the two sides at equal label w. blend is the largest share of the flux outside the
 * receiving side's velocity grid, over the momentum and energy shares and over both sides.
 *
 * Each cell owns its lower face in each direction, and the upper ghost cell owns the upper
 * boundary face.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_remap_flux_surf_inp struct.
 * @return Pointer to the updater.
 */
struct gkyl_dg_pkpm_remap_flux_surf *gkyl_dg_pkpm_remap_flux_surf_inew(
  const struct gkyl_dg_pkpm_remap_flux_surf_inp *inp
);

/**
 * Compute the face fluxes and time-step rates. A face with a nonpositive trace of rho or xi at
 * any node, where the frame map does not exist, has blend = 1.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range of the interior cells.
 * @param conf_range_ext Extended configuration-space range (indexes the configuration arrays).
 * @param phase_range Phase-space range of the interior cells.
 * @param phase_range_ext Extended phase-space range (indexes the phase-space arrays).
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param fluid Input physical fluid state [rho, rho ux, rho uy, rho uz, E] (first component used).
 * @param u Input flow velocity [ux, uy, uz].
 * @param xi Input thermal speed [xi, xi^2, p] (first component used).
 * @param theta Input T_perp/T.
 * @param fin Input Jacobian-weighted normalized distributions [F, G, M0_perp, N_perp].
 * @param conf_flux_surf Output kinetic fluxes, [dir][side][block][surface mode] on the owning
 *   cell of the face: side 0 is the flux the lower cell receives through its upper face, side 1
 *   the flux the upper cell receives through its lower face.
 * @param fluid_flux_surf Output fluid fluxes, [dir][component][surface mode].
 * @param cflrate_f Accumulated kinetic rate of the owning cell.
 * @param cflrate_fluid Accumulated fluid rate of the owning cell.
 * @param blend Output weight of the label-paired flux of the lower face per direction, [dir].
 */
void gkyl_dg_pkpm_remap_flux_surf_advance(
  struct gkyl_dg_pkpm_remap_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar,
  const struct gkyl_array *fluid, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *fluid_flux_surf, struct gkyl_array *cflrate_f,
  struct gkyl_array *cflrate_fluid, struct gkyl_array *blend
);

/**
 * Number of surface modes per block, direction and side of the kinetic flux (the array needs
 * 2 * 8 * cdim times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_remap_flux_surf *up);

/**
 * Number of surface modes per component and direction of the fluid flux (the array needs
 * 5 * cdim times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(const struct gkyl_dg_pkpm_remap_flux_surf *up
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_remap_flux_surf_release(struct gkyl_dg_pkpm_remap_flux_surf *up);
