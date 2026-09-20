#pragma once

#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>

// Object type.
typedef struct gkyl_dg_gr_maxwell_geom_source gkyl_dg_gr_maxwell_geom_source;

// Input packaged as a struct.
struct gkyl_dg_gr_maxwell_geom_source_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  enum gkyl_field_id field_id; // Maxwell model.
  double chi; // Electric-field divergence correction speed.
  double gamma; // Magnetic-field divergence correction speed.
  double K_phi; // Geometric source coefficient for electric cleaning field.
  double K_psi; // Geometric source coefficient for magnetic cleaning field.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater to compute GR Maxwell geometric source terms.
 *
 * This updater is intentionally separate from the GR Maxwell volume kernels:
 * it consumes conserved fields Jc*(D,B,phi,psi) directly and computes an RHS
 * increment in the same conserved layout.
 *
 * @param inp Input parameters.
 * @return New geometric source updater.
 */
struct gkyl_dg_gr_maxwell_geom_source*
gkyl_dg_gr_maxwell_geom_source_inew(const struct gkyl_dg_gr_maxwell_geom_source_inp *inp);

/**
 * Create new updater to compute GR Maxwell geometric source terms on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_gr_maxwell_geom_source*
gkyl_dg_gr_maxwell_geom_source_cu_dev_inew(const struct gkyl_dg_gr_maxwell_geom_source_inp *inp);

/**
 * Compute geometric source contribution to the EM RHS.
 *
 * @param up Geometric source updater.
 * @param conf_range Configuration-space range for indexing fields.
 * @param geom_factor_con Contravariant geometric source factors at volume nodes.
 * @param field_con Input conserved fields Jc*(D,B,phi,psi).
 * @param rhs Output RHS increment in conserved field layout.
 */
void
gkyl_dg_gr_maxwell_geom_source_advance(struct gkyl_dg_gr_maxwell_geom_source *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *geom_factor_con,
  const struct gkyl_array *field_con, struct gkyl_array *rhs);

/**
 * Host-side wrapper for computing GR Maxwell geometric source terms on device.
 */
void
gkyl_dg_gr_maxwell_geom_source_advance_cu(struct gkyl_dg_gr_maxwell_geom_source *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *geom_factor_con,
  const struct gkyl_array *field_con, struct gkyl_array *rhs);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void
gkyl_dg_gr_maxwell_geom_source_release(struct gkyl_dg_gr_maxwell_geom_source *up);
