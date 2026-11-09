#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_vlasov_position_map.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_pkpm_kinetic_conf_auxfields {
  const struct gkyl_array *bvar; // Magnetic field unit vector [bx, by, bz].
  const struct gkyl_array *u; // Flow velocity [ux, uy, uz].
  const struct gkyl_array *xi; // Thermal speed [xi, xi^2, p] (gkyl_dg_calc_pkpm_vars_prim).
  const struct gkyl_array *
    conf_flux_surf; // Modal fluxes the two sides of each configuration-space face receive, stored on the face's owning cell as [dir][side][block][mode].
};

// Input packaged as a struct
struct gkyl_dg_pkpm_kinetic_conf_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range for indexing auxiliary fields.
  const struct gkyl_range *phase_range; // Phase-space range for indexing the stored fluxes.
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space position map (REQUIRED; the streaming divergence carries 1/J_d).
  struct gkyl_dg_pkpm_kinetic_conf_auxfields auxfields; // Auxiliary fields.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new equation object for the configuration-space pass (stage A) of the
 * PKPM kinetic update of the eight blocks [F, G, M0_perp, N_perp]: the streaming volume term of
 * every block and the surface term scattering the stored configuration-space fluxes of the eight
 * blocks. The known velocity-space terms of the pass (w forces, mirror and coupling terms,
 * sources) are the stored operator of gkyl_dg_pkpm_harm_vol; the closure terms and the v_par
 * surface term belong to the velocity-space pass (stage B).
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_kinetic_conf_inp struct.
 * @return Pointer to the equation object.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_conf_inew(const struct gkyl_dg_pkpm_kinetic_conf_inp *inp);

/**
 * Create a new equation object for the configuration-space pass on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_conf_cu_dev_inew(
  const struct gkyl_dg_pkpm_kinetic_conf_inp *inp
);

/**
 * Set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_kinetic_conf_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin
);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_kinetic_conf_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin
);
#endif
