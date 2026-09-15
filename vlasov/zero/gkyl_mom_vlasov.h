#pragma once

#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_vlasov_velocity_map.h>

// Input packaged as a struct
struct gkyl_mom_vlasov_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  enum gkyl_model_id model_id; // enum to determine what type of Vlasov model (e.g., non-relativistic vs. relativistic).
  enum gkyl_hamil_id hamil_id; // enum for the Hamiltonian representation (sparse/dense velocity-space or phase-space expansion).
  const struct gkyl_range *vel_range; // Range for indexing velocity-space Jacobian and velocity map.
  const struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object. NULL => uniform velocity grid.
  const struct gkyl_range *hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_array *hamil; // Hamiltonian utilized to compute certain moments (such as energy or dH/dv moment). 
  double v_thresh; // Threshold velocity for moments over a subset of the domain such as M0_UPPER and M0_LOWER. 
  double f_thresh; // Threshold f for accumulating moments over a subset of the domain such as M0_UPPER and M0_LOWER. 
  enum gkyl_distribution_moments mom_type; // Name of moment to compute. See gkyl_mom_type.h. 
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new Vlasov moment type object. 
 *
 * @param inp Input parameters defined in gkyl_mom_vlasov_inp struct.
 * @return Pointer to Vlasov moment type object.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_inew(const struct gkyl_mom_vlasov_inp *inp);

/**
 * Create a new Vlasov moment type object that lives on 
 * NV-GPU. See new() method for documentation.
 *
 * @param inp Input parameters defined in gkyl_mom_vlasov_inp struct.
 * @return Pointer to Vlasov moment type object on device. 
 */
struct gkyl_mom_type* gkyl_mom_vlasov_cu_dev_inew(const struct gkyl_mom_vlasov_inp *inp);

/**
 * Create a new Vlasov integrated moment type object. 
 *
 * @param inp Input parameters defined in gkyl_mom_vlasov_inp struct.
 * @return Pointer to Vlasov moment type object.
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_inew(const struct gkyl_mom_vlasov_inp *inp);

/**
 * Create a new Vlasov integrated moment type object that lives on 
 * NV-GPU. See new() method for documentation.
 *
 * @param inp Input parameters defined in gkyl_mom_vlasov_inp struct.
 * @return Pointer to Vlasov moment type object on device. 
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_cu_dev_inew(const struct gkyl_mom_vlasov_inp *inp);
