#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>

// Input packaged as a struct
struct gkyl_dg_vlasov_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_range *conf_range; // Configuration-space range for use in indexing forces (EM fields, potentials, etc.).
  const struct gkyl_range *phase_range; // Range for indexing velocity-space flux.
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped.
  enum gkyl_model_id model_id; // enum to determine what type of Vlasov model (e.g., non-relativistic vs. relativistic).
  enum gkyl_hamil_id hamil_id; // enum for the Hamiltonian representation (sparse/dense velocity-space or phase-space expansion).
  bool has_E; // bool to determine whether we have electric fields (used for external forces too).
  bool has_phi; // bool to determine whether we have potentials (either electrostatic or gravitational).
  bool has_B; // bool to determine whether we have magnetic fields.
  bool has_rad; // bool to determine whether we have radiation drag forces.
  const struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object. Required: it also provides
                                                  // the velocity-space range used to index per-velocity-cell
                                                  // quantities (Jacobian, radiation drag).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object. Required: provides the
                                                  // (per-conf-cell constant) position-map Jacobian used to
                                                  // transform the streaming term to mapped coordinates.
  const struct gkyl_array *poisson_tensor_conf; // Configuration space poisson tensor used for nc poisson brackets.
  const struct gkyl_array *hamil; // Hamiltonian utilized to compute advection in configuration and velocity space. 
  const struct gkyl_array *qmem; // q/m*(E,B) electromagnetic fields (including external electromagnetic fields and forces).
  const struct gkyl_array *pot_tot; // (q/m*(phi + phi_ext) + m*phi_g, q/m*A_ext) total potentials. 
  const struct gkyl_array *rad; // Radiation drag forces, F_rad(v). 
  const struct gkyl_array *conf_flux_surf; // Modal expansion of fluxes at configuration space surfaces. 
  const struct gkyl_array *vel_flux_surf; // Modal expansion of fluxes at velocity space surfaces. 
  const struct gkyl_array *f_no_J; // Distribution function without velocity-space Jacobian. Used by magnetic field updates.
  bool use_lo; // bool to determine if using high-order kernels for non-canonical Hamiltonian models.  
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new Vlasov equation object. 
 *
 * @param inp Input parameters defined in gkyl_dg_vlasov_inp struct.
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_inew(const struct gkyl_dg_vlasov_inp *inp);

/**
 * Create a new Vlasov equation object that lives on NV-GPU
 *
 * @param inp Input parameters defined in gkyl_dg_vlasov_inp struct.
 * @return Pointer to Vlasov equation object on device. 
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_cu_dev_inew(const struct gkyl_dg_vlasov_inp *inp);
