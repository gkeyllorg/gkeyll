#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>

// Auxiliary fields for the passive gyrokinetic DG equation.
struct gkyl_dg_gyrokinetic_passive_auxfields {
  const struct gkyl_array *flux_surf; // Surface expansion of the passive phase-space flux.
  const struct gkyl_array *speeds;    // Conf-space passive advection speeds (all directions packed).
};

/**
 * Create a new passive gyrokinetic equation object.
 *
 * The passive equation advects f in conf-space only, using a pre-computed
 * surface expansion of the passive phase-space flux (flux_surf). No phi or
 * electromagnetic fields are involved.
 *
 * @param cbasis Configuration space basis functions.
 * @param pbasis Phase-space basis functions.
 * @param conf_range Configuration space range (for indexing aux fields).
 * @param phase_range Phase space range (for indexing flux_surf).
 * @param charge Species charge.
 * @param mass Species mass.
 * @param gk_geom Geometry struct.
 * @param vel_map Velocity space mapping object.
 * @param use_gpu Boolean to determine if equation object is on device.
 * @return Pointer to passive gyrokinetic equation object.
 */
struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_passive_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map,
  bool use_gpu);

/**
 * Set the auxiliary fields for the passive gyrokinetic equation.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gyrokinetic_passive_set_auxfields(const struct gkyl_dg_eqn *eqn,
  struct gkyl_dg_gyrokinetic_passive_auxfields auxin);
