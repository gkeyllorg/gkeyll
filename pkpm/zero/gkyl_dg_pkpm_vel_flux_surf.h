#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_pkpm_vel_flux_surf gkyl_dg_pkpm_vel_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_vel_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 hybrid, one velocity dimension).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the v_par surface fluxes of the PKPM kinetic
 * equations (velocity-space pass) of the eight blocks U = [F, G, M0_perp, N_perp].
 * Per configuration cell the w flux is formed at the 5^cdim face nodes once from the raw values of
 * the w-flux coefficients there (written with the stored operator of the stage by
 * gkyl_dg_pkpm_harm_vol_coef), theta, b and the closures: the speed
 * alpha = alpha_0 + (c - K) w_face - (b.grad xi) w_face^2, K = bb:grad u, and the coupling
 * T: F <- t.M0 + beta G, G <- t.N + beta H, M0 <- G t + R_B1 N, N <- H t + R_B1 O, with
 * H = 4 theta G - 2 theta^2 F and O = 6 theta N - 6 theta^2 M0; the coefficients t, beta and
 * R_B1 are listed in gkyl_pkpm_harm_slots.h. On the lower v_par face of every velocity cell the
 * HLL flux of alpha U + T(U) with the speeds alpha + the bounds of the spectrum of T is formed at
 * the nodes (pkpm_harm_node_wflux; exact upwinding when the speeds share a sign) and projected
 * onto the surface modal basis. Each cell owns its *lower* face, stored as [block][surface mode];
 * the faces at the two velocity domain edges are the closure's.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_vel_flux_surf_inp struct.
 * @return Pointer to v_par surface flux updater.
 */
struct gkyl_dg_pkpm_vel_flux_surf *gkyl_dg_pkpm_vel_flux_surf_inew(
  const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
);

/**
 * Create a new updater computing the v_par surface fluxes on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_pkpm_vel_flux_surf *gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(
  const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
);

/**
 * Compute the v_par surface fluxes on the INTERIOR faces (the lower face of every cell but the
 * lowest; the domain-edge faces are set by the closure updater gkyl_dg_calc_pkpm_vars_closure and
 * left untouched here).
 *
 * @param up v_par surface flux updater.
 * @param conf_range Configuration-space range for indexing the configuration fields.
 * @param phase_range Phase-space range to compute fluxes over.
 * @param wnode Input raw w-flux values at the face nodes (gkyl_dg_pkpm_harm_vol_coef).
 * @param theta Input T_perp/T.
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param alpha_0 Input pressure force a_P/xi closing the first moment.
 * @param c_h Input heating rate closing the scalar second moment.
 * @param fin Input distribution functions [F, G, M0_perp, N_perp].
 * @param vel_flux_surf Output modal v_par fluxes.
 *
 * No time-step rate: that of the w fluxes is the column rate of gkyl_dg_calc_pkpm_vars_closure.
 */
void gkyl_dg_pkpm_vel_flux_surf_advance(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *bvar, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *vel_flux_surf
);

/**
 * Host-side wrapper for computing v_par fluxes on device.
 */
void gkyl_dg_pkpm_vel_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *bvar, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *vel_flux_surf
);

/**
 * Number of surface modal basis functions per component in the flux array this
 * updater fills (the array needs 8 times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_vel_flux_surf *up);




/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_vel_flux_surf_release(struct gkyl_dg_pkpm_vel_flux_surf *up);
