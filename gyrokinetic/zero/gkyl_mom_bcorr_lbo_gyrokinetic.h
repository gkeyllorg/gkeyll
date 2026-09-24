#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_array.h>
#include <gkyl_velocity_map.h>

/**
 * Create new LBO Gyrokinetic boundary correction moment type object.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vBoundary Values at the edges of velocity space.
 * @param mass Mass of species
 * @param vel_map Velocity space mapping object.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_mom_bcorr_lbo_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis,  double mass,
  const struct gkyl_velocity_map *vel_map, bool use_gpu);
