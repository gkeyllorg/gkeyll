#pragma once

#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_range.h>

// Object type.
typedef struct gkyl_dg_gr_maxwell_current_deposition gkyl_dg_gr_maxwell_current_deposition;

// Input packaged as a struct.
struct gkyl_dg_gr_maxwell_current_deposition_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  int vdim; // Velocity-space dimensions.
  bool use_gpu; // Bool to determine if on GPU.
};

/**
 * Create a new updater to deposit GR-Maxwell current from triad moments.
 *
 * @param inp Input parameters.
 * @return Pointer to GR-Maxwell current deposition updater.
 */
struct gkyl_dg_gr_maxwell_current_deposition*
gkyl_dg_gr_maxwell_current_deposition_inew(const struct gkyl_dg_gr_maxwell_current_deposition_inp *inp);

/**
 * Create a new updater to deposit GR-Maxwell current from triad moments on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_gr_maxwell_current_deposition*
gkyl_dg_gr_maxwell_current_deposition_cu_dev_inew(
  const struct gkyl_dg_gr_maxwell_current_deposition_inp *inp);

/**
 * Accumulate current contribution into the Maxwell RHS.
 *
 * @param up GR-Maxwell current deposition updater.
 * @param conf_range Configuration-space range for indexing moments and geometry.
 * @param q_over_eps0 Signed charge over epsilon0; the kernel forms rho*beta - alpha*e^i_a*Jhat^a.
 * @param lapse Lapse at volume nodes.
 * @param shift Shift at volume nodes.
 * @param vierb_con Contravariant vierbein at volume nodes.
 * @param m0 Scalar density moment.
 * @param m1i Flat-frame momentum/current moment.
 * @param rhs Output Maxwell RHS; current is accumulated into the first 3 components.
 */
void
gkyl_dg_gr_maxwell_current_deposition_advance( struct gkyl_dg_gr_maxwell_current_deposition *up,
  const struct gkyl_range *conf_range, double q_over_eps0, const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift, const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *m0, const struct gkyl_array *m1i, struct gkyl_array *rhs);

/**
 * Host-side wrapper for GR-Maxwell current deposition on device.
 */
void
gkyl_dg_gr_maxwell_current_deposition_advance_cu(struct gkyl_dg_gr_maxwell_current_deposition *up,
  const struct gkyl_range *conf_range, double q_over_eps0, const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift, const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *m0, const struct gkyl_array *m1i, struct gkyl_array *rhs);

void
gkyl_dg_gr_maxwell_current_deposition_release(
  struct gkyl_dg_gr_maxwell_current_deposition *up);
