// Private header for use in the Gyrokinetic app: do not include in
// user-facing header files! Context (and device helpers) for the
// computational-to-physical coordinate mappings passed to proj_on_basis.
#pragma once

#include <gkyl_position_map.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>

// Context for c2p function passed to proj_on_basis.
struct gk_proj_on_basis_c2p_func_ctx {
  int cdim, vdim;
  struct gkyl_position_map *pos_map;
  struct gkyl_velocity_map *vel_map;
};

#ifdef GKYL_HAVE_CUDA

/**
 * Return the device address of the phase-space c2p function (position map in
 * configuration space, velocity map in velocity space), for projecting with
 * proj_on_basis on the GPU. The context passed along with it must be a
 * GPU-resident gk_proj_on_basis_c2p_func_ctx whose pos_map/vel_map members
 * are the maps' device (on_dev) objects. Defined in
 * gk_species_projection_cu.cu.
 */
proj_on_basis_c2p_t gk_species_projection_c2p_phase_func_cu_dev_ptr(void);

#endif
