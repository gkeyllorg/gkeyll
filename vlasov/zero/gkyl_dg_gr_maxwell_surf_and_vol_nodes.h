#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_evalf_def.h>

// Object type
typedef struct gkyl_surf_and_vol_node_arrays gkyl_surf_and_vol_node_arrays;
typedef struct gkyl_dg_gr_maxwell_surf_and_vol_nodes gkyl_dg_gr_maxwell_surf_and_vol_nodes;

typedef void (*dg_gr_maxwell_nodes_c2p_t)(const double *xcomp, double *xphys, void *ctx);

struct gkyl_surf_and_vol_node_arrays {
  struct gkyl_array *nodal_arr_vol; // array evaluated at local nodal locations (vol)
  struct gkyl_array *nodal_arr_surf_x; // array evaluated at local nodal locations (surf-x)
  struct gkyl_array *nodal_arr_surf_y; // array evaluated at local nodal locations (surf-y)
  struct gkyl_array *nodal_arr_surf_z; // array evaluated at local nodal locations (surf-z)
  struct gkyl_array
    *nodal_arr_vol_host; // Host side copy, array evaluated at local nodal locations (vol)
  struct gkyl_array
    *nodal_arr_surf_x_host; // Host side copy, array evaluated at local nodal locations (surf-x)
  struct gkyl_array
    *nodal_arr_surf_y_host; // Host side copy, array evaluated at local nodal locations (surf-y)
  struct gkyl_array
    *nodal_arr_surf_z_host; // Host side copy, array evaluated at local nodal locations (surf-z)
  bool use_gpu; // Boolean for using GPU
  int ndim; // Number of dimensions

  uint32_t flags;
  struct gkyl_ref_count ref_count;
  struct gkyl_surf_and_vol_node_arrays *on_dev; // Device side copy
};

struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp {
  int polyorder; // Polyorder
  const struct gkyl_rect_grid *grid; // Grid to project on
  const struct gkyl_basis *basis; // Basis functions
  int num_ret_vals; // Number of values 'eval' sets
  evalf_t eval; // Function to project.
  void *ctx; // Context for function evaluation. Can be NULL.

  dg_gr_maxwell_nodes_c2p_t c2p_func; // Function that transforms a set of ndim
    // computational coordinates to physical ones.
  void *c2p_func_ctx; // Context for c2p_func.
};

/**
 * Create new updater to compute function on nodes and calculate its
 * expansion on basis functions. Free using gkyl_dg_gr_maxwell_surf_and_vol_nodes_release
 * method.
 * 
 * @param inp Input parameters
 * @return New updater pointer.
 */
gkyl_dg_gr_maxwell_surf_and_vol_nodes *gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew(
  const struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp *inp
);

/**
 * Create new updater to compute function on volume and surface nodes. 
 * Free using gkyl_dg_gr_maxwell_surf_and_vol_nodes_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions to project on
 * @param num_ret_vals Number of values 'eval' sets
 * @param polyorder Polyorder of the spatial discretization used for determining number of nodes.
 * @param eval Function to project.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_dg_gr_maxwell_surf_and_vol_nodes *gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis, int num_ret_vals,
  int polyorder, evalf_t eval, void *ctx
);

/**
 * Create new stucture for storing all the arrays of values at various nodal and modal points
 * for both the surface and volume and on both host and gpu (if availible)
 *
 * @param info Projection information such as number of nodes requried to project the geometry
 * @param volume Volume to build the arrays over
 * @param use_gpu Boolean for using GPU
 * @return New updater pointer.
 */
gkyl_surf_and_vol_node_arrays *gkyl_surf_and_vol_node_arrays_new(
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *info, long volume, bool use_gpu
);

/**
 * Compute evaluation on nodes and corresponding expansion
 * coefficients. The update_rng MUST be a sub-range of the range on
 * which the array is defined. That is, it must be either the same
 * range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Eval on nodes updater to run
 * @param tm Time at which eval must be computed
 * @param update_rng Range on which to run eval.
 * @param out Output arrays at volume and nodal points
 */
void gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(
  const gkyl_dg_gr_maxwell_surf_and_vol_nodes *up, double tm, const struct gkyl_range *update_rng,
  struct gkyl_surf_and_vol_node_arrays *out
);

/**
 * Get the coordinates of a given node.
 *
 * @param up Project on basis updater.
 * @param node Index indicate the desired node.
 * @return Node coordinates.
 */
double *gkyl_dg_gr_maxwell_surf_and_vol_nodes_fetch_node(
  const gkyl_dg_gr_maxwell_surf_and_vol_nodes *up, long node
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(gkyl_dg_gr_maxwell_surf_and_vol_nodes *up);

/**
 * Copy surf and volume node array information to device.
 *
 * @param up Updater to delete.
 * @param ndim Number of dimensions.
 * @return New object on device.
 */
struct gkyl_surf_and_vol_node_arrays *gkyl_surf_and_vol_node_copy_to_device(
  struct gkyl_surf_and_vol_node_arrays *vol_surf_nodes, int ndim
);

/**
 * Cehck is strcut is a on cu dev.
 *
 * @param up Object to delete.
 * @return true if on cu dev, false otherwise.
 */
bool gkyl_surf_and_vol_node_arrays_is_cu_dev(const struct gkyl_surf_and_vol_node_arrays *up);

/**
 * Acquire array objects.
 *
 * @param up Object to acquire arrays from.
 * @return Acquired object.
 */
struct gkyl_surf_and_vol_node_arrays *gkyl_surf_and_vol_node_arrays_acquire(
  const struct gkyl_surf_and_vol_node_arrays *up
);

/**
 * Decrement reference count for array objects.
 *
 * @param up Object to decrement reference count too.
 */
void gkyl_surf_and_vol_node_arrays_release(const struct gkyl_surf_and_vol_node_arrays *up);
