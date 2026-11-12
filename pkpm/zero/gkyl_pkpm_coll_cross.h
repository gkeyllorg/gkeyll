#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Object type
typedef struct gkyl_pkpm_coll_cross gkyl_pkpm_coll_cross;

// A species in the cross-species exchange (gkyl_pkpm_coll_cross_advance).
struct gkyl_pkpm_coll_cross_species {
  double mass; // species mass
  const struct gkyl_array *prim; // node values [n, u (3), vtsq] before the step
  struct gkyl_array *fluid; // fluid state [rho, M (3), E] times the Jacobian, updated in place
  struct gkyl_array *temp_ratio; // output: temperature before the step over that after it
  int num_partners; // number of species it collides with
  int partner[GKYL_MAX_SPECIES]; // positions of its partners in the list of species
  const struct gkyl_array *coef[GKYL_MAX_SPECIES]; // exchange coefficient K with each partner
};

/**
 * Create an updater for the cross-species collisions of the PKPM collision operators (LBO and
 * BGK) at the 2^cdim interior Gauss-Legendre nodes of configuration space.
 *
 * Species s colliding with species r at the frequency nu_sr exchanges momentum and energy as
 * in the Greene closure with delta (1 + beta) = 2 and three degrees of freedom, vtsq = T/m:
 *
 *   R_s = K (u_r - u_s),  K = 2 A_s A_r / (A_s + A_r),  A_s = rho_s nu_sr,  A_r = rho_r nu_rs,
 *   Q_s = u_s . R_s + K [3 (T_r - T_s) + m_r |u_r - u_s|^2] / (m_s + m_r),
 *
 * so that R_s + R_r = 0 and Q_s + Q_r = 0 at every node. In terms of a cross flow and thermal
 * speed, R_s = rho_s nu_sr (u_sr - u_s) and Q_s - u_s . R_s = 3 rho_s nu_sr (vtsq_sr - vtsq_s),
 *
 *   u_sr = u_s + f_s (u_r - u_s),  f_s = 2 A_r / (A_s + A_r),
 *   vtsq_sr = vtsq_s + f_s [(m_r vtsq_r - m_s vtsq_s) + m_r |u_r - u_s|^2 / 3] / (m_s + m_r).
 *
 * The exchange is in the fluid states, and it is independent of the collision operator. In the
 * normalized variables the kinetic operator has no density, momentum or energy moments: the
 * drift is the frame acceleration and the heating the change of the normalization, which leaves
 * the self operator at the rate
 *
 *   nu_sr [1 + heat_fac (vtsq_sr - vtsq_s) / vtsq_s]
 *
 * on every block. heat_fac = 1 is the LBO, exactly (the operator at the rate
 * nu_sr vtsq_sr/vtsq_s). heat_fac = 2 is the BGK operator with the unit target at the rate
 * nu_sr + d ln(n T_s)/dt, which relaxes every pressure component at nu_sr. The two parts meet
 * only through this rate, so the exchange is advanced on its own (gkyl_pkpm_coll_cross_advance).
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param heat_fac Factor of the heating in the kinetic rate: 1 (LBO) or 2 (BGK).
 * @param use_gpu Whether to run on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_pkpm_coll_cross *gkyl_pkpm_coll_cross_new(
  const struct gkyl_basis *conf_basis, double heat_fac, bool use_gpu
);

/**
 * Compute the node values [n, ux, uy, uz, vtsq] of a species from its fluid state.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param mass Species mass.
 * @param jacob Position-map Jacobian (first component, constant in a cell).
 * @param fluid Fluid state [rho, M (3), E] times the Jacobian (p=1 expansions).
 * @param prim Output node values [n, ux, uy, uz, vtsq].
 */
void gkyl_pkpm_coll_cross_prim(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double mass,
  const struct gkyl_array *jacob, const struct gkyl_array *fluid, struct gkyl_array *prim
);

/**
 * Compute the exchange coefficient K = 2 A_s A_r / (A_s + A_r) of a pair at the nodes. It is
 * the same for both members of the pair (zero where both frequencies vanish).
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param mass_self Mass of the species.
 * @param mass_other Mass of the partner.
 * @param nu_self Frequency nu_sr of the species colliding with the partner, at the nodes.
 * @param nu_other Frequency nu_rs of the partner colliding with the species, at the nodes.
 * @param prim_self Node values [n, u (3), vtsq] of the species.
 * @param prim_other Node values [n, u (3), vtsq] of the partner.
 * @param coef Output exchange coefficient at the nodes.
 */
void gkyl_pkpm_coll_cross_coef(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double mass_self,
  double mass_other, const struct gkyl_array *nu_self, const struct gkyl_array *nu_other,
  const struct gkyl_array *prim_self, const struct gkyl_array *prim_other, struct gkyl_array *coef
);

/**
 * Accumulate the rate that one partner adds to the kinetic operator of a species at the nodes,
 * nu_sr [1 + heat_fac (vtsq_sr - vtsq_s) / vtsq_s].
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param mass_self Mass of the species.
 * @param mass_other Mass of the partner.
 * @param nu_self Frequency nu_sr of the species colliding with the partner, at the nodes.
 * @param nu_other Frequency nu_rs of the partner colliding with the species, at the nodes.
 * @param prim_self Node values [n, u (3), vtsq] of the species.
 * @param prim_other Node values [n, u (3), vtsq] of the partner.
 * @param nu_kin Node values of the rate of the kinetic operator, incremented.
 */
void gkyl_pkpm_coll_cross_kinetic_rate(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double mass_self,
  double mass_other, const struct gkyl_array *nu_self, const struct gkyl_array *nu_other,
  const struct gkyl_array *prim_self, const struct gkyl_array *prim_other, struct gkyl_array *nu_kin
);

/**
 * Floor the rate of the kinetic operator (self plus every partner) at zero at the nodes.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param nu_kin Node values of the rate of the kinetic operator.
 */
void gkyl_pkpm_coll_cross_floor(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range,
  struct gkyl_array *nu_kin
);

/**
 * Backward-Euler step of the momentum and energy exchange of every species with its partners,
 * at the nodes, with the exchange coefficients of the state before the step:
 *
 *   rho_s (u_s' - u_s) = dt sum_r K_sr (u_r' - u_s')
 *   (3/2) n_s (T_s' - T_s) = dt sum_r K_sr [3 (T_r' - T_s') + m_r |u_r' - u_s'|^2] / (m_s + m_r)
 *                            + rho_s |u_s' - u_s|^2 / 2,
 *
 * the backward-Euler step of the momentum and of the total energy of each species. The total
 * momentum and the total energy of the species are conserved at every node. Each flow after the
 * step is a weighted average of the flows before it and each temperature a weighted average of
 * the temperatures before it plus a non-negative heating, so the temperatures stay positive.
 * The fluid states are updated in place; the densities do not change.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param dt Time step.
 * @param jacob Position-map Jacobian (first component, constant in a cell).
 * @param num_species Number of species in the exchange.
 * @param species The species: masses, node values before the step, partners and exchange
 *   coefficients (coef of the pair (s, r) must equal that of (r, s)), fluid states and outputs.
 */
void gkyl_pkpm_coll_cross_advance(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double dt,
  const struct gkyl_array *jacob, int num_species,
  const struct gkyl_pkpm_coll_cross_species species[]
);

/**
 * Shape factor of the implicit BGK step of the distributions of a species at the nodes,
 *
 *   S = min(1, temp_ratio / (1 + nu_tot dt)),
 *
 * with temp_ratio the temperature before the exchange over that after it: every pressure
 * component, the product of the temperature and a normalized moment, then takes the
 * backward-Euler step at the frequency nu_tot.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param dt Time step.
 * @param nu_tot Collision frequency (self plus every partner) at the nodes.
 * @param temp_ratio Temperature before the exchange over that after it, at the nodes.
 * @param shape Output shape factor at the nodes.
 */
void gkyl_pkpm_coll_cross_shape(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double dt,
  const struct gkyl_array *nu_tot, const struct gkyl_array *temp_ratio, struct gkyl_array *shape
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_coll_cross_release(struct gkyl_pkpm_coll_cross *up);
