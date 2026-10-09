#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>
#include <gkyl_gk_bc_type.h>

// Forward declaration of msgpack metadata (defined in gkyl_array_rio.h).
struct gkyl_msgpack_data;

// Object type
typedef struct gkyl_bc_sheath_gyrokinetic gkyl_bc_sheath_gyrokinetic;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param vel_map Velocity space mapping object.
 * @param cdim Configuration space dimensions.
 * @param q2Dm charge-to-mass ratio times 2.
 * @param type Type of gyrokinetic BC to apply.
 * @param surrogate_model_path Path to the .kann surrogate model file (NULL if not using surrogate).
 * @param phase_grid Phase-space grid, used to build the vcutsq diagnostic grid.
 * @param phase_local Local phase-space range.
 * @param poly_order Polynomial order of the basis functions (used to build the vcutsq diagnostic grid).
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_sheath_gyrokinetic *gkyl_bc_sheath_gyrokinetic_new(
  int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *basis, const struct gkyl_range *skin_r,
  const struct gkyl_range *ghost_r, const struct gkyl_velocity_map *vel_map, int cdim, double q2Dm,
  enum gkyl_gyrokinetic_bc_type type, const char *surrogate_model_path,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_range *phase_local,
  unsigned int poly_order, bool use_gpu
);

/**
 * Compute the cutting velicity and apply the sheath BC with the bc_sheath_gyrokinetic object.
 *
 * @param up BC updater.
 * @param phi Electrostatic potential at the magnetic presheath entrance (simulation boundary).
 * @param phi_wall Electrostatic potential at the wall.
 * @param density Density at the magnetic presheath entrance.
 * @param temperature Temperature at the magnetic presheath entrance.
 * @param bmag Magnetic field strength at the magnetic presheath entrance.
 * @param bimpact_angle Angle of the magnetic field with respect to the wall at the magnetic presheath entrance.
 * @param distf Distribution function array to apply BC to.
 * @param conf_r Configuration space range (to index phi).
 */
void gkyl_bc_sheath_gyrokinetic_advance(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *density,
  const struct gkyl_array *temperature, const struct gkyl_array *bmag,
  const struct gkyl_array *bimpact_angle, struct gkyl_array *distf, const struct gkyl_range *conf_r
);

/**
 * Set the vcut_fact array used in the sheath BC. 
 * This is used to implement a mu-dependent vcut in the sheath BC, 
 * where vcut_fact is the mu dependent multiplying factor in the expression for vcut.
 * 
 * @param up BC updater.
 * @param vcut_fact The vcut_fact array to use in the sheath BC.
 */
void gkyl_bc_sheath_gyrokinetic_set_vcutsq(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *vcut_fact
);

/**
 * Write out the vcutsq array.
 *
 * @param up BC updater.
 * @param meta Metadata to write (may be NULL).
 * @param fname Name of output file (including .gkyl extension).
 */
void gkyl_bc_sheath_gyrokinetic_write_vcutsq(
  struct gkyl_bc_sheath_gyrokinetic *up, struct gkyl_msgpack_data *meta, const char *fname
);

/**
 * Free memory associated with bc_sheath_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up);