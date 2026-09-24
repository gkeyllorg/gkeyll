#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_ambi_bolt_potential gkyl_ambi_bolt_potential;

/**
 * Create new updater to compute the electrostatic potential assuming ambipolar
 * sheath particle fluxes and Boltzmann isothermal electrons.
 *
 * @param grid Cartesian grid dynamic field is defined on.
 * @param basis Basis object (configuration space).
 * @param mass_e Electron mass.
 * @param charge_e Electron charge.
 * @param temp_e Electron temperature.
 * @param use_gpu Boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_ambi_bolt_potential* gkyl_ambi_bolt_potential_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, double mass_e, double charge_e, double temp_e,
  bool use_gpu);

/**
 * Compute the ion density and electrostatic potential at the sheath entrance.
 * Below we use jac to refer to the fact that some of these quantities are
 * multiplied by the configuration space Jacobian (e.g. for field aligned
 * coordinates). In reality they are also multiplied by the phase space
 * Jacobian (approx B), but we focus on jac because need to divide by it to
 * get the potential (or multiply by 1/jac).
 *
 * @param Ambipolar, Boltzmann electron sheath potential updater. 
 * @param edge Lower (-1) or upper (1) boundary along the field line.
 * @param skin_r Global skin range along the field line.
 * @param ghost_r Corresponding global ghost range along the field line.
 * @param cmag Clebsch function in definition of the magnetic field.
 * @param jacobtot_inv Reciprocal of the phase-space and conf-space Jacobians (1/(J*B)).
 * @param gammai Ion particle flux at the sheath entrance times the
 *               conf-space Jacobian.
 * @param m0i Ion number density.
 * @param Jm0i Ion number density times the conf-space Jacobian.
 * @param sheath_vals Ion number density and potential at the sheath entrance.
 */
void
gkyl_ambi_bolt_potential_sheath_calc(struct gkyl_ambi_bolt_potential *up, enum gkyl_edge_loc edge,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_array *cmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *gammai, const struct gkyl_array *m0i, const struct gkyl_array *Jm0i,
  struct gkyl_array *sheath_vals);

/**
 * Compute the electrostatic potential in the domain as
 *  phi = phi_s + (T_e/e)*ln(n_i/n_is).
 *
 * @param up Ambipolar, Boltzmann electron sheath potential updater. 
 * @param local_r Local range.
 * @param extlocal_r Extended local range.
 * @param jacob_geo_inv Reciprocal of the configuration space Jacobian.
 * @param m0i Ion number density times the conf-space Jacobian.
 * @param sheath_vals Ion number density and potential at the sheath entrance.
 * @param phi electrostatic potential.
 */
void
gkyl_ambi_bolt_potential_phi_calc(struct gkyl_ambi_bolt_potential *up,
  const struct gkyl_range *local_r, const struct gkyl_range *extlocal_r,
  const struct gkyl_array *m0i, const struct gkyl_array *sheath_vals,
  struct gkyl_array *phi);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_ambi_bolt_potential_release(gkyl_ambi_bolt_potential *up);
