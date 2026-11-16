#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_em_coupling gkyl_dg_calc_pkpm_em_coupling;

// Input packaged as a struct
struct gkyl_dg_calc_pkpm_em_coupling_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_range *mem_range; // Configuration-space range that sets the size of the batched solve.
  int num_species; // Number of fluid species in the update.
  const double *qbym; // Charge/mass of each fluid species (num_species values).
  double epsilon0; // Permittivity of free space.
  bool pkpm_field_static; // bool to determine if the EM fields are static (no dE/dt update).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create new updater for the locally implicit (time-centered) fluid-Maxwell source update
 * of the p=1 PKPM system, solved nodally at the 2^cdim Gauss-Legendre nodes of each cell
 * (Wang et al. JCP 2020 in the unknowns [Mbar_s, epsilon0 Ebar]):
 *
 *   Mbar_s - dt/2 q_s/m_s (rho_s Ebar + Mbar_s x B^n) = M_s^n + dt/2 (q_s/m_s rho_s E_ext + rho_s a_s)
 *   epsilon0 Ebar + dt/2 sum_s q_s/m_s Mbar_s = epsilon0 E^n - dt/2 J_app
 *
 * then M^{n+1} = 2 Mbar - M^n, E^{n+1} = 2 Ebar - E^n. The products rho E and M x B of two
 * p=1 fields are integrated exactly by the 2-point rule, so the nodal solve equals the
 * modal solve with exact weak products (no aliasing), and sum_s |M_s|^2/(2 rho_s) +
 * epsilon0 |E|^2/2 is conserved node by node. rho^n and B^n are unchanged by the sources.
 *
 * @param inp Input parameters defined in gkyl_dg_calc_pkpm_em_coupling_inp struct.
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_em_coupling *gkyl_dg_calc_pkpm_em_coupling_inew(const struct gkyl_dg_calc_pkpm_em_coupling_inp *inp);

/**
 * Create new updater for the implicit fluid-Maxwell source update on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_pkpm_em_coupling *gkyl_dg_calc_pkpm_em_coupling_cu_dev_inew(const struct gkyl_dg_calc_pkpm_em_coupling_inp *inp);

/**
 * Advance the fluid state of every species and the electric field by the implicit source
 * update (in place): the momentum and E by the nodal implicit solve, the fluid energy by
 * the change of the bulk kinetic energy (the work of the field); the density is unchanged.
 *
 * @param up                     Updater for computing fluid-EM coupling in the PKPM system.
 * @param dt                     Time step.
 * @param app_accel[num_species] Input array of applied accelerations (external forces) on each species
 * @param ext_em                 Input array of external electromagnetic fields
 * @param app_current            Input array of applied currents
 * @param fluid[num_species]     Fluid state [rho, rho ux, rho uy, rho uz, E] of each species: old time in, new time out
 * @param em                     EM fields: electric field at the old time in, new time out
 */
void gkyl_dg_calc_pkpm_em_coupling_advance(
  struct gkyl_dg_calc_pkpm_em_coupling *up, double dt,
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *ext_em,
  const struct gkyl_array *app_current, struct gkyl_array *fluid[GKYL_MAX_SPECIES], struct gkyl_array *em
);

/**
 * Delete the updater for the implicit fluid-Maxwell source update.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_em_coupling_release(struct gkyl_dg_calc_pkpm_em_coupling *up);

/**
 * Host-side wrapper for the implicit fluid-Maxwell source update on device
 */

void gkyl_dg_calc_pkpm_em_coupling_advance_cu(
  struct gkyl_dg_calc_pkpm_em_coupling *up, double dt,
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *ext_em,
  const struct gkyl_array *app_current, struct gkyl_array *fluid[GKYL_MAX_SPECIES], struct gkyl_array *em
);
