#pragma once

#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type.
typedef struct gkyl_dg_gr_maxwell_lorentz_conf gkyl_dg_gr_maxwell_lorentz_conf;

// Input packaged as a struct.
struct gkyl_dg_gr_maxwell_lorentz_conf_inp {
  const struct gkyl_vlasov_position_map *pos_map; // REQUIRED: provides the cell-constant jacob_pos.
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  int vdim; // Velocity-space dimensions.
  double qbym; // Charge over mass.
  double chi; // Electric-field divergence correction speed.
  double gamma; // Magnetic-field divergence correction speed.
  double K_phi; // Electric-field cleaning damping coefficient.
  double K_psi; // Magnetic-field cleaning damping coefficient.
  bool use_gpu; // Bool to determine if on GPU.
};

/**
 * Create a new updater to compute local Lorentz-force fields from GR-Maxwell
 * fields and vierbein geometry.
 *
 * @param inp Input parameters defined in gkyl_dg_gr_maxwell_lorentz_conf_inp.
 * @return Pointer to GR-Maxwell Lorentz-force updater.
 */
struct gkyl_dg_gr_maxwell_lorentz_conf*
gkyl_dg_gr_maxwell_lorentz_conf_inew(const struct gkyl_dg_gr_maxwell_lorentz_conf_inp *inp);

/**
 * Create a new updater to compute local Lorentz-force fields on NV-GPU. See
 * new() method for documentation.
 */
struct gkyl_dg_gr_maxwell_lorentz_conf*
gkyl_dg_gr_maxwell_lorentz_conf_cu_dev_inew(const struct gkyl_dg_gr_maxwell_lorentz_conf_inp *inp);

/**
 * Add q/m times the local Lorentz-force fields to qmem.
 *
 * qmem may already contain applied accelerations or external EM fields; this
 * updater accumulates the transformed GR-Maxwell contribution into it.
 *
 * @param up GR-Maxwell Lorentz-force updater.
 * @param conf_range Configuration-space range for indexing fields and geometry.
 * @param lapse Lapse at volume nodes.
 * @param shift Shift at volume nodes.
 * @param h_ij Covariant spatial metric at volume nodes.
 * @param h_ij_inv Contravariant spatial metric at volume nodes.
 * @param det_h Spatial metric Jacobian at volume nodes.
 * @param vierb_cov Covariant vierbein at volume nodes.
 * @param vierb_con Contravariant vierbein at volume nodes.
 * @param field_con Input contravariant (JD,JB) field components.
 * @param qmem Output q/m*(E,B) force memory.
 */
void
gkyl_dg_gr_maxwell_lorentz_conf_advance(struct gkyl_dg_gr_maxwell_lorentz_conf *up, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_surf_and_vol_node_arrays *vierb_cov, 
  const struct gkyl_surf_and_vol_node_arrays *vierb_con, const struct gkyl_array *field_con, struct gkyl_array *qmem);

/**
 * Host-side wrapper for local GR-Maxwell Lorentz-force fields on device.
 */
void
gkyl_dg_gr_maxwell_lorentz_conf_advance_cu(struct gkyl_dg_gr_maxwell_lorentz_conf *up, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_surf_and_vol_node_arrays *vierb_cov, 
  const struct gkyl_surf_and_vol_node_arrays *vierb_con, const struct gkyl_array *field_con, struct gkyl_array *qmem);

void
gkyl_dg_gr_maxwell_lorentz_conf_release(struct gkyl_dg_gr_maxwell_lorentz_conf* up);
