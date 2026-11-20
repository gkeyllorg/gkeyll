#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_vars gkyl_dg_calc_pkpm_vars;

// Input packaged as a struct
struct gkyl_dg_calc_pkpm_vars_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space position map (REQUIRED; the gradients carry 1/J_d).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create new updater for the configuration-space variables of the PKPM
 * solver (p=1 serendipity): the nodal primitive variables (the flow M/rho and the thermal
 * speed; nodal quotients, the weak division at p=1), the closures with the velocity-space
 * edges, the in-cell gradients and grad(ln n), the vertex positivity limiter, the constraints,
 * the transversality projection, the explicit source, the initial-condition normalization and
 * the diagnostics.
 *
 * @param inp Input parameters defined in gkyl_dg_calc_pkpm_vars_inp struct.
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_vars *gkyl_dg_calc_pkpm_vars_inew(
  const struct gkyl_dg_calc_pkpm_vars_inp *inp
);

/**
 * Create new updater on NV-GPU.
 */
struct gkyl_dg_calc_pkpm_vars *gkyl_dg_calc_pkpm_vars_cu_dev_inew(
  const struct gkyl_dg_calc_pkpm_vars_inp *inp
);

/**
 * Compute the in-cell gradient of a three-component frame variable (u, or [rho, 0, 0] and [xi, ...]),
 * with the share blend of each face jump lifted into the volume: the p=1 expansions of d_i q_j for
 * every configuration direction i and component j, stored per cell as [direction][component][mode]
 * (cdim * 3 * num_basis entries).
 *
 * The remap face flux maps the frame jumps exactly, so its share (1 - blend) of a jump belongs to
 * the face. The label-paired flux does not, so its share is lifted into the volume:
 *   <psi_k, d_i q_j> = dx1 [psi_k(+1) qhat_up - psi_k(-1) qhat_lo - <d_xi psi_k, q>],
 *   qhat_up = q(+1) + blend_up/2 [q]_up,  qhat_lo = q(-1) - blend_lo/2 [q]_lo,
 * which is the plain in-cell gradient at blend = 0 and the weak (surface-averaged) gradient at
 * blend = 1. It needs q in the neighbors of every cell of conf_range (the ghost cells). It is the
 * gradient the kinetic sources are formed from. p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range (interior cells)
 * @param blend      Input weight of the label-paired flux of each lower face per direction
 *                   (gkyl_dg_pkpm_remap_flux_surf_advance), including the upper ghost-owned faces
 * @param q          Input three-component field, defined on the extended range
 * @param grad       Output gradient
 */
void gkyl_dg_calc_pkpm_vars_in_cell_grad(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *blend, const struct gkyl_array *q, struct gkyl_array *grad
);

/**
 * Closures c_0, alpha_0 = a_P/xi, c and a_perp = I_perp a of the normalized kinetic update of
 * the eight blocks [F, G, M0_perp, N_perp], with the characteristic velocity-space edges: per
 * Gauss-Legendre node, the constraints M0(F) = 1, M1(F) = 0, M2(F) + 2 M0(G) = 3 and
 * I_perp int M0_perp dw = 0 of the complete right-hand side (configuration pass + closure terms +
 * the HLL edge fluxes against an empty state outside the velocity domain, constructed at the
 * w-face nodes as on the interior w faces and projected onto p=1), a 3x3 solve after the
 * perpendicular rows are eliminated. The edge fluxes of all eight
 * blocks are written in the stored face format: the lower one into the lowest velocity cell's
 * slot of vel_flux_surf, the upper one into the slot of the velocity ghost cell above the
 * uppermost cell (the array layout must have that ghost cell). The time-step rate of all the
 * velocity-space terms of the column (the w fluxes with the new closures and the closure sources,
 * bounded with absolute coefficients at the larger |w| edge: no w of the column exceeds it) is
 * added to cflrate in every cell of the column. The edge speeds use the closures of the previous stage
 * (the output arrays on input). At a node with M0(F) <= 0 or a non-positive determinant all
 * closures and edge fluxes of the cell are zero; the number of such nodes is returned. p=1
 * serendipity only.
 *
 * @param up            Updater for computing pkpm variables
 * @param conf_range    Configuration space range
 * @param phase_range   Phase-space range (its velocity extent locates the edge cells)
 * @param phase_grid    Phase-space grid (velocity extent and cell size)
 * @param norm_moms     Raw moments [M0, M2, M0_G, M1, M3, M1_G, int M0_perp dw] of the normalized
 *                      distributions
 * @param rhs_moms      The same moments of the configuration-pass right-hand side
 * @param bvar          Magnetic field unit vector
 * @param wnode         Raw w-flux values at the w-face nodes (gkyl_dg_pkpm_harm_vol_coef)
 * @param theta         T_perp/T
 * @param fin           Normalized distributions [F, G, M0_perp, N_perp]
 * @param c_0           Output normalization closure (input: previous stage)
 * @param alpha_0       Output pressure force a_P/xi (input: previous stage)
 * @param c_h           Output heating rate (input: previous stage)
 * @param a_perp        Output perpendicular force [x, y, z]
 * @param vel_flux_surf Stored velocity face fluxes: the slots of the lowest cells (lower
 *                      edge) and of the velocity ghost cells above the uppermost cells
 *                      (upper edge) are written; phase_range must be laid out with that
 *                      ghost cell
 * @param cflrate       Phase-space time-step rates, incremented in every cell of each column
 * @return Number of nodes where the closures are undefined
 */
int gkyl_dg_calc_pkpm_vars_closure(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0,
  struct gkyl_array *c_h, struct gkyl_array *a_perp, struct gkyl_array *vel_flux_surf,
  struct gkyl_array *cflrate
);

/**
 * Compute the nodal primitive variables of the fluid state
 * [rho, rho ux, rho uy, rho uz, E]: the flow u = M/rho and [xi, xi^2, p] with
 * xi^2 = (2 E/rho - |u|^2)/3, p = rho xi^2, by the nodal quotient at the 2^cdim
 * Gauss-Legendre nodes interpolated to p=1. For p=1 the nodal quotient is the weak
 * division: int phi rho u = int phi M for every p=1 phi (the 2-point rule integrates
 * rho u phi exactly). Nodes with rho <= 0 get zero, nodes with xi^2 <= 0 zero [xi, xi^2, p];
 * both are counted. p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range
 * @param fluid      Input fluid state
 * @param u          Output flow velocity (3 expansions per cell)
 * @param xi         Output [xi, xi^2, p] (3 expansions per cell)
 * @return Number of inadmissible nodes in this call
 */
int gkyl_dg_calc_pkpm_vars_prim(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
);

/**
 * Cell integrals of the diagnostics of the [rho, M, E] state: [rho, rho ux, rho uy,
 * rho uz, rho ux^2, rho uy^2, rho uz^2, p_par, p_perp] (accumulated). p=1 serendipity only.
 *
 * @param up               Updater for computing pkpm variables
 * @param conf_range       Configuration space range
 * @param fluid            Input fluid state [rho, rho ux, rho uy, rho uz, E]
 * @param vlasov_pkpm_moms Input kinetic coupling moments [p_par, p_perp, q_par, q_perp]
 * @param u                Input flow velocity
 * @param int_pkpm_vars    Output cell integrals (9 per cell)
 */
void gkyl_dg_calc_pkpm_vars_int(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
);

/**
 * The constraints of the normalized distributions, [M0(F), M1(F), (M2(F) + 2 M0(G))/3,
 * int M0_perp dw (x, y, z)], from the raw normalized moments [M0, M2, M0_G, M1, M3, M1_G,
 * int M0_perp dw (x, y, z), ...]: they are [1, 0, 1, 0, 0, 0] in every mode of every cell when
 * the closures hold.
 *
 * @param up          Updater (not used)
 * @param conf_range  Configuration space range
 * @param num_basis   Number of configuration basis functions
 * @param norm_moms   Input normalized moments (at least 9 expansions per cell)
 * @param constraints Output constraints (6 expansions per cell)
 */
void gkyl_dg_calc_pkpm_vars_constraints(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, int num_basis,
  const struct gkyl_array *norm_moms, struct gkyl_array *constraints
);

/**
 * Transversality projection of the first-harmonic blocks: M0_perp and N_perp of every phase cell
 * of phase_range are projected onto the plane perpendicular to b at the 2^cdim interior
 * Gauss-Legendre nodes of configuration space, for every Legendre mode in w, in place. p=1
 * serendipity only.
 *
 * @param up          Updater for computing pkpm variables
 * @param conf_range  Configuration-space range indexing bvar
 * @param phase_range Phase-space range of the cells to project
 * @param bvar        Magnetic field unit vector (on the cells of phase_range)
 * @param f           Distributions [F, G, M0_perp, N_perp], projected in place
 */
void gkyl_dg_calc_pkpm_vars_transverse(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, struct gkyl_array *f
);

/**
 * Conditioning of the closure solves of the last gkyl_dg_calc_pkpm_vars_closure call: the
 * smallest determinant of the per-node 3x3 system relative to its edge-free value (1 when
 * the velocity-domain edges carry no trace; small or negative when they do, i.e. the
 * velocity extents are too small for the distribution).
 *
 * @param up Updater.
 * @return Minimum det/det0 over the nodes of the last call.
 */
double gkyl_dg_calc_pkpm_vars_closure_cond_min(const struct gkyl_dg_calc_pkpm_vars *up);

/**
 * Compute grad(ln n) = grad(rho)/rho, the configuration-space factor of the source
 * -v_par b.grad(ln n) of the normalized kinetic equations, by the nodal quotient of the
 * gradient of rho (computed with gkyl_dg_calc_pkpm_vars_in_cell_grad on the array
 * [rho, 0, 0]) at the 2^cdim Gauss-Legendre nodes. At a node where rho <= 0 the quotient is
 * set to zero; the number of such nodes is returned. p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range
 * @param rho        Input mass density
 * @param grad_rho   Input gradient of [rho, 0, 0] (cdim x 3 expansions per cell)
 * @param grad_ln_n  Output grad(ln n) (cdim expansions per cell)
 * @return Number of nodes where rho <= 0 in this call
 */
int gkyl_dg_calc_pkpm_vars_grad_ln_n(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
);

/**
 * Positivity limiter of the fluid state [rho, M, E] (Zhang-Shu scaling, exact for the
 * p=1 serendipity basis whose extreme values over a cell are its 2^cdim vertex values):
 * in every cell the density slopes are scaled so that the vertex density stays above
 * eps times the cell average, then all slopes are scaled jointly so that the vertex
 * internal energy E - |M|^2/(2 rho) stays above eps times that of the cell average. The
 * state is then admissible at the vertices, hence by convexity at the Gauss-Legendre nodes,
 * where the thermal speed is formed (gkyl_dg_calc_pkpm_vars_prim).
 * The cell averages are untouched. A cell whose average is itself inadmissible is left as is.
 * The per-cell record `count` (3 components: density scalings, energy scalings,
 * inadmissible averages) is incremented where the limiter acted, and the numbers of cells
 * of each kind in this call are returned through num_rho, num_e and num_avg.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range
 * @param eps        Relative floor of the vertex density and internal energy
 * @param fluid      In/out fluid state [rho, rho ux, rho uy, rho uz, E]
 * @param count      In/out per-cell record of the limiter actions (3 components)
 * @param num_rho    Output: cells whose density slopes were scaled in this call
 * @param num_e      Output: cells whose slopes were scaled for the internal energy
 * @param num_avg    Output: cells with an inadmissible average (left as is)
 */
void gkyl_dg_calc_pkpm_vars_positivity(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, double eps,
  struct gkyl_array *fluid, struct gkyl_array *count, int *num_rho, int *num_e, int *num_avg
);

/**
 * Rescale the distributions [F, G] to the constraints M0(F) = 1 and M2(F) + 2 M0(G) = 3
 * (initial conditions): fbar = [F/M0, G s/M0] with s = (3 - M2/M0)/(2 M0_G/M0), the factors
 * formed from the moments of f at the 2^cdim Gauss-Legendre nodes of configuration space. At
 * a node where M0 <= 0 or M0_G <= 0 both factors are set to zero and the node is counted. The
 * other blocks of fbar are not written. p=1 serendipity only.
 *
 * @param up          Updater for computing pkpm variables
 * @param conf_range  Configuration space range for indexing norm_moms
 * @param phase_range Phase space range to operate on
 * @param norm_moms   Input raw moments [M0, M2, M0_G, ...] of f (mass 1)
 * @param f           Input distributions [F, G, ...]
 * @param fbar        Output rescaled distributions [F, G]
 * @return Number of nodes where M0 <= 0 or M0_G <= 0 in this call
 */
int gkyl_dg_calc_pkpm_vars_normalize_dist(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *norm_moms,
  const struct gkyl_array *f, struct gkyl_array *fbar
);

/**
 * Construct the diagnostic arrays:
 * pkpm_prim_io = [P_xx, P_xy, P_xz, P_yy, P_yz, P_zz, q_x, q_y, q_z] (9 components), the
 * pressure tensor and the heat flux of the fluid equations,
 *   P_ij = (p_par - p_perp) b_i b_j + p_perp delta_ij + s_i b_j + b_i s_j,
 *   q_i = b_i (q_par + 2 q_perp)/2 + r_i, and
 * pkpm_vars_io = [u_x, u_y, u_z, T_perp/m, a_P, div b, bb:grad u, T_perp/m div b, xi, c]
 * (10 components), with the products projected onto p=1. p=1 serendipity only.
 *
 * @param up               Updater for computing pkpm variables
 * @param conf_range       Configuration space range
 * @param vlasov_pkpm_moms Input kinetic coupling moments [p_par, p_perp, q_par, q_perp,
 *                         s (x, y, z), r (x, y, z)]
 * @param bvar             Input magnetic field unit vector [bx, by, bz]
 * @param grad_b           Input weak gradient of b (gkyl_dg_calc_pkpm_bvar_grad)
 * @param u                Input flow velocity [ux, uy, uz]
 * @param grad_u           Input gradient of u (gkyl_dg_calc_pkpm_vars_in_cell_grad)
 * @param theta            Input T_perp/T
 * @param alpha_0          Input pressure force a_P/xi
 * @param xi               Input thermal speed [xi, xi^2, p]
 * @param c_h              Input heating rate c
 * @param pkpm_prim_io     Output pressure tensor and heat flux (9 components)
 * @param pkpm_vars_io     Output pkpm variables (10 components)
 */
void gkyl_dg_calc_pkpm_vars_diag(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *bvar,
  const struct gkyl_array *grad_b, const struct gkyl_array *u, const struct gkyl_array *grad_u,
  const struct gkyl_array *theta, const struct gkyl_array *alpha_0, const struct gkyl_array *xi,
  const struct gkyl_array *c_h, struct gkyl_array *pkpm_prim_io, struct gkyl_array *pkpm_vars_io
);

/**
 * Explicit electromagnetic source of the fluid state: rhs is incremented by
 * [0, q/m (rho E + M x B), q/m M.E]. The CFL rate of each cell is raised to the rate of the
 * source where that is larger: the plasma plus the cyclotron frequency, or the rate of the
 * force, from the cell averages (the source kernels of the DG Euler fluid). p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range
 * @param qmem       Input acceleration and q/m times the magnetic field
 *                   [ax, ay, az, q/m Bx, q/m By, q/m Bz, ...]
 * @param omega_p_sq Input array (1 component): plasma frequency squared of each cell, summed
 *                   over the species
 * @param fluid      Input fluid state [rho, rho ux, rho uy, rho uz, E]
 * @param cflrate    CFL rate of the fluid (1 component), raised to the rate of the source
 * @param rhs        Output right-hand side of the fluid state, incremented
 */
void gkyl_dg_calc_pkpm_vars_source(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *omega_p_sq,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Delete pointer to updater to compute pkpm variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_vars_release(struct gkyl_dg_calc_pkpm_vars *up);

/**
 * Host-side wrappers for pkpm vars operations on device
 */

int gkyl_dg_calc_pkpm_vars_closure_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0,
  struct gkyl_array *c_h, struct gkyl_array *a_perp, struct gkyl_array *vel_flux_surf,
  struct gkyl_array *cflrate
);
void gkyl_dg_calc_pkpm_vars_transverse_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, struct gkyl_array *f
);
int gkyl_dg_calc_pkpm_vars_prim_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
);
void gkyl_dg_calc_pkpm_vars_int_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
);
void gkyl_dg_calc_pkpm_vars_constraints_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, int num_basis,
  const struct gkyl_array *norm_moms, struct gkyl_array *constraints
);
int gkyl_dg_calc_pkpm_vars_grad_ln_n_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
);
void gkyl_dg_calc_pkpm_vars_positivity_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, double eps,
  struct gkyl_array *fluid, struct gkyl_array *count, int *num_rho, int *num_e, int *num_avg
);
int gkyl_dg_calc_pkpm_vars_normalize_dist_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *norm_moms,
  const struct gkyl_array *f, struct gkyl_array *fbar
);
void gkyl_dg_calc_pkpm_vars_diag_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *bvar,
  const struct gkyl_array *grad_b, const struct gkyl_array *u, const struct gkyl_array *grad_u,
  const struct gkyl_array *theta, const struct gkyl_array *alpha_0, const struct gkyl_array *xi,
  const struct gkyl_array *c_h, struct gkyl_array *pkpm_prim_io, struct gkyl_array *pkpm_vars_io
);

void gkyl_dg_calc_pkpm_vars_in_cell_grad_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *blend, const struct gkyl_array *q, struct gkyl_array *grad
);

void gkyl_dg_calc_pkpm_vars_source_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *omega_p_sq,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
);
