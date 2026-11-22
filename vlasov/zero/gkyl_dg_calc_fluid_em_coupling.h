#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_dg_calc_fluid_em_coupling gkyl_dg_calc_fluid_em_coupling;

// Input packaged as a struct
struct gkyl_dg_calc_fluid_em_coupling_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1).
  const struct gkyl_range
    *mem_range; // Configuration-space range that sets the size of the batched solve.
  int num_species; // Number of fluid species in the update.
  const double *qbym; // Charge/mass of each fluid species (num_species values).
  // Whether the state of each species carries an energy, [rho, rho u, E] (num_species values);
  // false for an isothermal fluid, [rho, rho u]. NULL: all species carry one.
  const bool *has_energy;
  double epsilon0; // Permittivity of free space.
  // True if the EM fields are static: the species currents do not enter the E rows (E then
  // changes only by the applied current handed to the advance method).
  bool is_field_static;
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create new updater for the locally implicit (time-centered) fluid-Maxwell source update
 * of p=1 DG fluids with state [rho, rho u, (E)], solved at each of the 2^cdim
 * Gauss-Legendre nodes of each cell (one linear system per node; Wang et al. JCP 2020 in the
 * unknowns [Mbar_s, epsilon0 Ebar]), with a resistive layer dE/dt = -sigma E (sigma >= 0):
 *
 *   Mbar_s - dt/2 q_s/m_s (rho_s Ebar + Mbar_s x B^n) = M_s^n + dt/2 (q_s/m_s rho_s E_ext + rho_s a_s)
 *   epsilon0 Ebar + dt phi2 sum_s q_s/m_s Mbar_s = epsilon0 phi E^n - dt phi2 J_app
 *
 * then M^{n+1} = 2 Mbar - M^n, E^{n+1} = (phi/phi2) Ebar - (phi^2/phi2 - g2) E^n, with
 * x = sigma dt, g2 = exp(-x), phi = (1 - g2)/x, phi2 = (1 - phi)/x. E^{n+1} is the exact solution
 * of dE/dt = -sigma E + S with the current S held at its midpoint value, and Ebar its exact time
 * average; at x = 0 (phi = 1, phi2 = 1/2) this is the implicit midpoint rule. The products rho E
 * and M x B of two p=1 fields are integrated exactly by the 2-point rule, so without a layer the
 * nodal solve equals the modal solve with exact weak products (no aliasing). Node by node,
 * sum_s |M_s|^2/(2 rho_s) + epsilon0 |E|^2/2 changes by -epsilon0 sigma int |E|^2 dt -
 * dt J_app.Ebar. rho^n and B^n are unchanged by the sources.
 *
 * @param inp Input parameters defined in gkyl_dg_calc_fluid_em_coupling_inp struct.
 * @return New updater pointer.
 */
struct gkyl_dg_calc_fluid_em_coupling *gkyl_dg_calc_fluid_em_coupling_inew(
  const struct gkyl_dg_calc_fluid_em_coupling_inp *inp
);

/**
 * Create new updater for the implicit fluid-Maxwell source update on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_fluid_em_coupling *gkyl_dg_calc_fluid_em_coupling_cu_dev_inew(
  const struct gkyl_dg_calc_fluid_em_coupling_inp *inp
);

/**
 * Advance the fluid state of every species and the electric field by the implicit source
 * update (in place): the momentum and E by the nodal implicit solve, the fluid energy by
 * the change of the bulk kinetic energy (the work of the field); the density is unchanged.
 *
 * @param up                     Updater for computing fluid-EM coupling.
 * @param dt                     Time step.
 * @param app_accel[num_species] Input array of applied accelerations (external forces) on each species
 * @param ext_em                 Input array of external electromagnetic fields
 * @param app_current            Input array of applied currents
 * @param sigma                  Input array of the resistive layer, the damping rate of E at the 2^cdim Gauss-Legendre nodes of each cell (NULL: none)
 * @param fluid[num_species]     Fluid state [rho, rho ux, rho uy, rho uz, E] of each species (no E without energy): old time in, new time out
 * @param em                     EM fields: electric field at the old time in, new time out
 */
void gkyl_dg_calc_fluid_em_coupling_advance(
  struct gkyl_dg_calc_fluid_em_coupling *up, double dt,
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *ext_em,
  const struct gkyl_array *app_current, const struct gkyl_array *sigma,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], struct gkyl_array *em
);

/**
 * Delete the updater for the implicit fluid-Maxwell source update.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_fluid_em_coupling_release(struct gkyl_dg_calc_fluid_em_coupling *up);

/**
 * Host-side wrapper for the implicit fluid-Maxwell source update on device
 */

void gkyl_dg_calc_fluid_em_coupling_advance_cu(
  struct gkyl_dg_calc_fluid_em_coupling *up, double dt,
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *ext_em,
  const struct gkyl_array *app_current, const struct gkyl_array *sigma,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], struct gkyl_array *em
);
