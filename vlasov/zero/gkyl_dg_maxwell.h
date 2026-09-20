#pragma once

#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_range.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>


// Input packaged as a struct
struct gkyl_dg_maxwell_inp {
  const struct gkyl_basis *cbasis; // Configuration-space basis functions.
  const struct gkyl_range *crange; // Configuration-space range for use in indexing forces (EM fields, potentials, etc.).
  const struct gkyl_array *jacob_pos; // Configuration-space (position-map) Jacobian, per conf cell (required).
                                      // Callers without a nonuniform map pass an identity position map's Jacobian.
  enum gkyl_field_id field_id; // enum to determine what type of maxwell model (e.g., non-relativistic vs. relativistic).
  double lightSpeed; // Speed of light
  double elcErrorSpeedFactor; // Factor multiplying lightSpeed for div E correction
  double mgnErrorSpeedFactor; // Factor multiplying lightSpeed for div B correction
  const struct gkyl_array *conf_flux_surf; // Modal expansion of fluxes at configuration space surfaces. 
  const struct gkyl_surf_and_vol_node_arrays *lapse; // nodal expansion of lapse
  const struct gkyl_surf_and_vol_node_arrays *shift; // nodal expansion of shift
  const struct gkyl_surf_and_vol_node_arrays *h_ij; // nodal expansion of spatial
  const struct gkyl_surf_and_vol_node_arrays *h_ij_inv; // nodal expansion of inverse spatial metric
  const struct gkyl_surf_and_vol_node_arrays *det_h; // nodal expansion of spatial metric determinant
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new Maxwell equation object.
 *
 * @param inp Input parameters defined in gkyl_dg_maxwell_inp struct.
 * @return Pointer to maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_inew(const struct gkyl_dg_maxwell_inp *inp);

/*
 * Create a new Maxwell equation object that lives on NV-GPU.
 *
 * @param inp Input parameters defined in gkyl_dg_maxwell_inp struct.
 * @return Pointer to maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_cu_dev_inew(const struct gkyl_dg_maxwell_inp *inp);
