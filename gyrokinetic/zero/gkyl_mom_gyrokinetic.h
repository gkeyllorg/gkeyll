#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

/**
 * Create new gyrokinetic moment type object. Valid 'mom' strings are "M0",
 * "M1", "M2", "M2par", "M2perp", "M3par", "M3perp", "ThreeMoments",
 * "FourMoments", "HamiltonianMoments".
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param conf_range Configuration-space range.
 * @param mass Mass of species.
 * @param charge Charge of species.
 * @param vel_map Velocity space mapping object.
 * @param gk_geom Geometry object.
 * @param phi Electrostatic potential (for Hamiltonian moment).
 * @param mom_type Name of moment to compute.
 * @param use_gpu bool to determine if on GPU.
 */
struct gkyl_mom_type* 
gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  double mass, double charge, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool use_gpu);

/**
 * Create new integrated gyrokinetic moment type object.
 * Valid 'mom' strings are "ThreeMoments", "FourMoments", "HamiltonianMoments".
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param conf_range Configuration-space range.
 * @param mass Mass of species.
 * @param vel_map Velocity space mapping object.
 * @param gk_geom Geometry object.
 * @param phi Electrostatic potential (for Hamiltonian moment).
 * @param mom_type Name of moment to compute.
 * @param use_gpu bool to determine if on GPU.
 */
struct gkyl_mom_type* 
gkyl_int_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  double mass, double charge, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool use_gpu);
