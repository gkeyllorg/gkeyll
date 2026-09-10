#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>
#include <gkyl_util.h>

#include <math.h>

typedef void (*mapc2p_t)(double t, const double *zc, double *vp, void *ctx);

// Velocity space mappings.
struct gkyl_mapc2p_inp {
  mapc2p_t mapping; // univariate mapping vp[0](zc[0]), vp[1](zc[1]), etc.
  void *ctx;  // Context for mapping.
};

// Object type.
typedef struct gkyl_velocity_map gkyl_velocity_map;

struct gkyl_velocity_map {
  bool is_identity; // =True if comp coords = phys coords.
  struct gkyl_rect_grid grid; // Phase space grid.
  struct gkyl_rect_grid grid_vel; // Velocity space grid.
  struct gkyl_range local, local_ext; // Local & extended local phase-space range.
  struct gkyl_range local_vel, local_ext_vel; // Local & extended local velocity-space range.
  struct gkyl_array *jacobvel; // Velocity space Jacobian.
  struct gkyl_basis *vmap_basis;  // Basis for velocity mapping.
  struct gkyl_array *vmap; // Velocity mapping in each velocity direction.
  struct gkyl_array *vmap_prime; // Derivative of the velocity mappings.
  struct gkyl_array *vmap_sq; // Velocity mapping in each velocity direction squared.
  struct gkyl_velocity_map *on_dev; // Device copy of itself.
  double vbounds[2*GKYL_MAX_VDIM]; // Velocity at the boundaries.
  // For internal/private use only:
  struct gkyl_array *vmap_ho; // Host copy of vmap.
  struct gkyl_basis vmap_basis_ho;  // Host basis for velocity mapping.
  uint32_t flags;
  struct gkyl_ref_count ref_count;
};

/**
 * Create a new velocity map object, containing the things needed
 * to use a (univariate) velocity map (e.g. to use nonununiform velocity
 * grids).
 *
 * @param mapc2p_in Comp. to phys. mapping input object (see definition above).
 * @param grid Phase space grid.
 * @param grid_vel Velocity grid.
 * @param local Local phase range.
 * @param local_ext Local extended phase range.
 * @param local_vel Local velocity range.
 * @param local_ext_vel Local extended velocity range.
 * @param use_gpu Whether to create a device copy of this new object.
 * @return New velocity map object.
 */
struct gkyl_velocity_map* gkyl_velocity_map_new(struct gkyl_mapc2p_inp mapc2p_in,
  struct gkyl_rect_grid grid, struct gkyl_rect_grid grid_vel,
  struct gkyl_range local, struct gkyl_range local_ext,
  struct gkyl_range local_vel, struct gkyl_range local_ext_vel, bool use_gpu);

/**
 * Write the velocity map and its jacobian to file.
 *
 * @param gvm Velocity map object.
 * @param species_comm Species communicator.
 * @param app_name Name of the app.
 * @param species_name Name of the species.
 */
void
gkyl_velocity_map_write(const struct gkyl_velocity_map* gvm, struct gkyl_comm* species_comm,
  const char* app_name, const char* species_name);

/**
 * Evaluate the velocity mappings at the v-space boundary to get the
 * boundary velocity values. NOTE: done on the host.
 *
 * @param gvm Velocity map object.
 * @param vbounds Host array of v-space boundary values (size 2*GKYL_MAX_VDIM).
 */
void
gkyl_velocity_map_get_boundary_values(const struct gkyl_velocity_map* gvm, double *vbounds);

/**
 * Reduce (i.e. get the min or max) the cell length in physical
 * velocity units over the whole velocity grid.
 *
 * @param gvm Velocity map object.
 * @param op GKYL_MIN or GKYL_MAX.
 * @param dv_m Min/Max velocity cell length on the grid.
 */
void
gkyl_velocity_map_reduce_dv(const struct gkyl_velocity_map* gvm, enum gkyl_array_op op, double *dv_m);

/**
 * Reduce (i.e. get the min or max) the cell length in physical
 * velocity units over a specified velocity range.
 *
 * @param gvm Velocity map object.
 * @param op GKYL_MIN or GKYL_MAX.
 * @param range_vel Velocity range to get min/max dv from.
 * @param dv_m Min/Max velocity cell length on the grid.
 */
void
gkyl_velocity_map_reduce_dv_range(const struct gkyl_velocity_map* gvm, enum gkyl_array_op op,
  double *dv_m, struct gkyl_range range_vel);

/**
 * Evaluate the velocity mapping at a specific computational (velocity) coordinate.
 * NOTE: done on the host.
 *
 * @param gvm Velocity map object.
 * @param zc Computational velocity coordinates.
 * @param vp Resulting physical velocity coordinates.
 */
void
gkyl_velocity_map_eval_c2p(const struct gkyl_velocity_map* gvm, const double *zc, double *vp);

/**
 * Evaluate the velocity mapping at a specific computational (velocity)
 * coordinate. Unlike gkyl_velocity_map_eval_c2p, this variant uses the vmap
 * and vmap_basis members (not their _ho host copies), so it is callable
 * inside device kernels when given the device object (gvm->on_dev). On GPU
 * builds do not call it on the host with the host object, whose vmap data
 * lives on the device.
 *
 * @param gvm Velocity map object (on_dev object inside kernels).
 * @param zc Computational velocity coordinates.
 * @param vp Resulting physical velocity coordinates.
 */
GKYL_CU_DH static inline void
gkyl_velocity_map_eval_c2p_dev(const struct gkyl_velocity_map* gvm, const double *zc, double *vp)
{
  // Find the index of the cell containing zc.
  int idx_zc[GKYL_MAX_VDIM];
  for (int d=0; d<gvm->local_ext_vel.ndim; d++) {
    int idx = gvm->local_ext_vel.lower[d] + (int) floor((zc[d] - (gvm->grid_vel.lower[d]) )/gvm->grid_vel.dx[d]);
    // Bound idx to the range in the grid. If it falls outside of that is due
    // to floating point arithmetic, or due to an error in the code.
    idx = GKYL_MIN2(idx, gvm->local_ext_vel.upper[d]);
    idx = GKYL_MAX2(idx, gvm->local_ext_vel.lower[d]);
    idx_zc[d] = idx;
  }

  // Fetch DG coefficients of the velocity map in idx_zc.
  long lidx_zc = gkyl_range_idx(&gvm->local_ext_vel, idx_zc);
  const double *vmap_c = (const double *) gkyl_array_cfetch(gvm->vmap, lidx_zc);

  double zc_cc[GKYL_MAX_VDIM];
  gkyl_rect_grid_cell_center(&gvm->grid_vel, idx_zc, zc_cc);

  for (int d=0; d<gvm->local_ext_vel.ndim; d++) {
    // Convert computational to logical coord.
    double zlog[] = {(zc[d] - zc_cc[d]) / (0.5*gvm->grid_vel.dx[d])};
    // Evaluate vmap expansion at logical coord.
    vp[d] = gvm->vmap_basis->eval_expand(zlog, &vmap_c[d*gvm->vmap_basis->num_basis]);
  }
}

/**
 * Indicate if this velocity map object is allocated on the GPU.
 *
 * @param gvm Velocity map object.
 * @return
 */
bool gkyl_velocity_map_is_cu_dev(const struct gkyl_velocity_map* gvm);

/**
 * Create a new pointer to the velocity map object.
 * Release it with gkyl_velocity_map_release.
 *
 * @param New pointer to the velocity map object.
 */
struct gkyl_velocity_map* gkyl_velocity_map_acquire(const struct gkyl_velocity_map* gvm);

/**
 * Release pointer to (and eventually memory associated with)
 * the velocity map object. 
 *
 * @param Velocity map object.
 */
void gkyl_velocity_map_release(const struct gkyl_velocity_map *gvm);
