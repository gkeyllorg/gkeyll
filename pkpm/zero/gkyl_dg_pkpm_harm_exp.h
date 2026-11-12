#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_pkpm_harm_exp gkyl_dg_pkpm_harm_exp;

// Input packaged as a struct
struct gkyl_dg_pkpm_harm_exp_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid (the velocity cell width).
  const struct gkyl_basis *conf_basis; // Configuration-space basis (p=1 serendipity).
  double qbym; // Charge to mass ratio of the species.
  double mix; // Gyrophase-mixing rate in units of the local |Omega| (0: no mixing).
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater for the gyrophase-mixing integrating factor of the first-harmonic blocks
 * M0_perp and N_perp of the normalized PKPM distributions [F, G, M0_perp, N_perp] (plan 17.44). The
 * gyration (q/m) U x B is explicit in the stage right-hand side, so every block is advanced by plain
 * SSP-RK3 (|Omega| is part of the time-step rate, gkyl_dg_pkpm_harm_vol_advance); only the mixing
 * damping, which can be stiff, is integrated exactly.
 *
 * Gyrophase mixing (plan 17.24): agyrotropy beyond what one harmonic can represent relaxes to that
 * bound, d(M0_perp, N_perp)/dt = -nu ((M0_perp, N_perp) - P(M0_perp, N_perp)) with nu = mix |Omega|
 * and P the projection onto r <= 1 at each phase node, r = max(sqrt2 |M0_perp|/sqrt(F G),
 * |N_perp| sqrt(F)/G^(3/2)). As a rate on the first harmonic it is nu (1 - 1/r)_+, frozen at the
 * start of the step (capped at ln(r)/h, which takes the unforced state to the bound in one step), on
 * the plane perpendicular to b_n, the field direction at the start of the step. With L = -nu (1 - 1/r)_+
 * on that plane, SSP-RK3 advances V = e^(-(t - t_n) L) U with the full stage right-hand sides N_i (the
 * mixing is not part of them). In the Shu-Osher form of the driver (stages at t_n, t_n + h, t_n + h/2):
 *
 *   U_2 = e^(hL) (M + h N_1),
 *   U_3 = 3/4 e^(hL/2) M + 1/4 e^(-hL/2) (U_2 + h N_2),
 *   M^(n+1) = 1/3 e^(hL) M + 2/3 e^(hL/2) (U_3 + h N_3).
 *
 * The free relaxation is exact at any rate (to the bound in one step at the cap) and the gyration
 * keeps its SSP-RK3 error; the forced response has an error that grows with the rate times dt (stage 2
 * amplifies its forcing by e^(nu h/2) <= sqrt(r)). Where no node of a phase cell is damped the update
 * is SSP-RK3 itself and the cell is not touched. The mixing changes int M0_perp dw, which
 * gkyl_dg_pkpm_harm_exp_restore restores. F and G are not touched. p=1 serendipity only.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_dg_pkpm_harm_exp *gkyl_dg_pkpm_harm_exp_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp);

/**
 * Create the updater on NV-GPU. See the new method for documentation.
 */
struct gkyl_dg_pkpm_harm_exp *gkyl_dg_pkpm_harm_exp_cu_dev_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp);

/**
 * Number of components of the configuration-space coefficient array ((q/m) B_n at the
 * configuration nodes).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_exp_num_coef(const struct gkyl_dg_pkpm_harm_exp *up);

/**
 * Number of components of the phase-space array of the gyrophase-mixing rate (one per phase node).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_exp_num_rate(const struct gkyl_dg_pkpm_harm_exp *up);

/**
 * Freeze the field of the mixing at the start of the step: (q/m) B_n at the configuration nodes.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param em Total physical electromagnetic field at the start of the step.
 * @param coef Coefficients at the configuration nodes (output).
 */
void gkyl_dg_pkpm_harm_exp_coef(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *coef
);

/**
 * Replace the first-harmonic blocks of the forward Euler output of stage s (1, 2 or 3), out = U_s +
 * h rhs_s, at the damped phase nodes by the value X_s whose SSP-RK3 combination in the driver with the
 * state M at the start of the step gives the integrating-factor stage: X_1 = U_2, X_2 = 4 U_3 - 3 M
 * (3/4 M + 1/4 X_2 = U_3) and X_3 = (3 M^(n+1) - M)/2 (1/3 M + 2/3 X_3 = M^(n+1)), i.e. with Y = out,
 *
 *   X_1 = e^(hL) Y,  X_2 = 3 (e^(hL/2) - 1) M + e^(-hL/2) Y,  X_3 = (e^(hL) - 1) M / 2 + e^(hL/2) Y.
 *
 * At stage 1 (U_1 = M) the mixing rate of every phase node is set first, and act is set to 1 in the
 * configuration cells with a damped node (act is not cleared here). Without mixing (mix <= 0) nothing
 * is done. F and G of out are not touched.
 *
 * The factors 3 and 1/2 of M are the ratios c_M / c_Y of the weights of the combinations. If the
 * driver combines with the weights of gkyl_pkpm_stage_weights instead, w_M / w_Y = c_M / (c_Y lam),
 * pass lam: the factors are divided by it, and the weighted combination of M and X_s is the
 * integrating-factor stage with those weights. Not available on the device.
 *
 * @param up Updater.
 * @param stage Stage (1, 2 or 3).
 * @param h Time step.
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range.
 * @param coef Frozen coefficients.
 * @param f0 State M at the start of the step (eight blocks).
 * @param lam Temperature factor of the stage weights at the configuration nodes (ordered as in
 *   gkyl_dg_conf_nodal), or 0 for the SSP-RK3 weights.
 * @param nu Mixing rate at the phase nodes (written at stage 1).
 * @param act Configuration-space flag of damped columns (set at stage 1).
 * @param out Forward Euler output of the stage (eight blocks; blocks 2-7 replaced at damped nodes).
 */
void gkyl_dg_pkpm_harm_exp_advance(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  const struct gkyl_array *lam, struct gkyl_array *nu, struct gkyl_array *act, struct gkyl_array *out
);

/**
 * Number of components of the configuration-space array of the w integrals of the restoration.
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_exp_num_col(const struct gkyl_dg_pkpm_harm_exp *up);

/**
 * Restore int M0_perp dw = 0 in the configuration cells flagged by act (those with a damped node):
 * with m = int M0_perp dw / int F dw at the configuration nodes, M0_perp -= F m and
 * N_perp -= 2 G m at every w. (F, 2G) dv is the first-harmonic change of a perpendicular frame shift
 * dv, the direction of the perpendicular closure force, so the momentum the mixing removes is
 * returned as a frame shift; F and G (and with them the other constraints and the energy) are not
 * changed. The products are exact at the phase nodes, so int M0_perp dw = 0 holds to round-off. Applied
 * to a stage output it also holds for the SSP-RK3 combinations of the driver (they are linear).
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range (its velocity extent is the column).
 * @param act Configuration-space flag of damped columns.
 * @param col Configuration-space work array (gkyl_dg_pkpm_harm_exp_num_col components).
 * @param out Stage output (eight blocks), restored in place.
 */
void gkyl_dg_pkpm_harm_exp_restore(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *act, struct gkyl_array *col,
  struct gkyl_array *out
);

/**
 * Host-side wrappers of the device kernels (see the host methods for documentation).
 */
void gkyl_dg_pkpm_harm_exp_coef_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *coef
);
void gkyl_dg_pkpm_harm_exp_advance_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  struct gkyl_array *nu, struct gkyl_array *act, struct gkyl_array *out
);
void gkyl_dg_pkpm_harm_exp_restore_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *act, struct gkyl_array *col,
  struct gkyl_array *out
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_harm_exp_release(struct gkyl_dg_pkpm_harm_exp *up);
