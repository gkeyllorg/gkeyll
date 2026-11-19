#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>

// Input packaged as a struct
struct gkyl_dg_pkpm_kinetic_vel_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range for indexing the configuration-space arrays.
  const struct gkyl_range *phase_range; // Phase-space range for indexing the stored fluxes.
  const struct gkyl_array *c_0; // Normalization closure c_0 (source c_0 U on every block).
  const struct gkyl_array *alpha_0; // Pressure force a_P/xi closing the first moment.
  const struct gkyl_array *c_h; // Heating rate closing the scalar second moment.
  const struct gkyl_array *a_perp; // Perpendicular force closing int M0_perp dw [x, y, z].
  const struct gkyl_array *vel_flux_surf; // Modal expansion of the fluxes on the v_par faces (each cell's lower face; the velocity ghost cell's slot is the upper domain edge).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new equation object for the velocity-space pass of the PKPM
 * kinetic update of the eight blocks [F, G, M0_perp, N_perp]: the volume terms of the closures
 * (forces (alpha_0 + c w) U, sources c_0 U and k c U with k = 0, 2, 1, 3 for F, G, M0_perp
 * and N_perp, and the a_perp couplings) and the surface term
 * scattering the stored v_par fluxes (the HLL edge fluxes of the closure at the velocity domain
 * edges). Run after the configuration-space pass and the closures, accumulating into the same
 * right-hand side.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_kinetic_vel_inp struct.
 * @return Pointer to the equation object.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_vel_inew(const struct gkyl_dg_pkpm_kinetic_vel_inp *inp);

/**
 * Create a new equation object for the velocity-space pass on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_vel_cu_dev_inew(
  const struct gkyl_dg_pkpm_kinetic_vel_inp *inp
);
