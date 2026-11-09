#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_pkpm_harm_vol gkyl_dg_pkpm_harm_vol;

// Input packaged as a struct
struct gkyl_dg_pkpm_harm_vol_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis (p=1 hybrid, one velocity dimension).
  double qbym; // Charge to mass ratio of the species (gyration of the first harmonic).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space position map (REQUIRED; the flux divergence carries 1/J_d).
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater for the volume integrals of the known velocity-space terms of the
 * configuration-space pass (stage A) of the normalized PKPM kinetic update, on the eight blocks
 * [F, G, M0_perp (x, y, z), N_perp (x, y, z)]: the w-flux volume terms int d_w phi Phi, the
 * sources int phi S and the configuration-space flux couplings int d_x phi . T of the first
 * harmonic (the streaming of all eight blocks stays in the configuration-pass equation object).
 *
 * The stored operator (gkyl_pkpm_harm_ops.h) is built once per stage into a configuration-space
 * array (the coef method, generated kernels pkpm_harm_build_*): every coefficient is a product of
 * the precomputed p=1 fields, stored as its exact projection onto degree <= 2 per direction (by
 * Gauss-Legendre quadrature exact for its degree): that is all that any Galerkin integral against
 * (test function) x (distribution) sees. The advance method forms the matrices of each
 * configuration cell (pkpm_harm_mat_*) and applies them per phase cell (pkpm_harm_vol_*). p=1
 * serendipity only.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_dg_pkpm_harm_vol *gkyl_dg_pkpm_harm_vol_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp);

/**
 * Create the updater on NV-GPU. See the inew method for documentation.
 */
struct gkyl_dg_pkpm_harm_vol *gkyl_dg_pkpm_harm_vol_cu_dev_inew(
  const struct gkyl_dg_pkpm_harm_vol_inp *inp
);

/**
 * Number of components per configuration cell of the stored operator (entries x 3^cdim modes).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_vol_num_coef(const struct gkyl_dg_pkpm_harm_vol *up);

/**
 * Number of components per configuration cell of the raw w-flux values at the w-face nodes
 * (PKPM_HS_W_END values at each of the 5^cdim nodes; gkyl_pkpm_harm_slots.h).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_vol_num_wnode(const struct gkyl_dg_pkpm_harm_vol *up);

/**
 * Build the stored operator of the configuration cells of conf_range, and the raw values of the
 * w-flux coefficients (the same pointwise products, before projection) at the w-face nodes: the
 * w fluxes (interior faces and the closure's velocity-domain edges) are formed there from the
 * actual local operator.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range of the cells to build.
 * @param bvar Magnetic field unit vector [bx, by, bz].
 * @param grad_b Weak gradient of b, [d_i bx, d_i by, d_i bz] per direction i.
 * @param u Flow velocity [ux, uy, uz].
 * @param grad_u Gradient of u of the sources (in-cell plus the blended face lift).
 * @param xi Thermal speed [xi, xi^2, p] (the first expansion is used).
 * @param grad_xi Weak gradient of [xi, 0, 0].
 * @param grad_ln_n grad(ln n), cdim expansions.
 * @param theta T_perp/T.
 * @param dtB_over_magB d_t B/|B| of the stage.
 * @param em Total physical electromagnetic field of the stage (B at component 3).
 * @param coef Output coefficients (gkyl_dg_pkpm_harm_vol_num_coef components).
 * @param wnode Output w-flux values at the w-face nodes (gkyl_dg_pkpm_harm_vol_num_wnode components).
 */
void gkyl_dg_pkpm_harm_vol_coef(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef,
  struct gkyl_array *wnode
);

/**
 * Host-side wrapper of the coef method on NV-GPU.
 */
void gkyl_dg_pkpm_harm_vol_coef_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef,
  struct gkyl_array *wnode
);

/**
 * Accumulate the first-harmonic volume terms into rhs, and the bound of their source rate over each
 * column (absolute coefficients over the configuration cell, at the larger |w| edge; plan 17.27),
 * plus the bound of the gyration rate |(q/m) B| over the cell (plan 17.44), into cflrate in every
 * cell of the column.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range indexing the coefficients.
 * @param phase_range Phase-space range of the cells to update.
 * @param coef Stored coefficients of the stage (gkyl_dg_pkpm_harm_vol_coef).
 * @param fin Distributions [F, G, M0_perp, N_perp].
 * @param cflrate Phase-space time-step rates, incremented.
 * @param rhs Right-hand side, incremented.
 */
void gkyl_dg_pkpm_harm_vol_advance(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Host-side wrapper of the advance method on NV-GPU.
 */
void gkyl_dg_pkpm_harm_vol_advance_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_harm_vol_release(struct gkyl_dg_pkpm_harm_vol *up);
