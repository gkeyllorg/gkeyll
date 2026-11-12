#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes_priv.h>
#include <gkyl_range.h>

struct gkyl_dg_gr_maxwell_surf_and_vol_nodes {
  struct gkyl_rect_grid grid;
  int num_ret_vals; // number of values returned by eval function
  evalf_t eval; // function to project
  void *ctx; // evaluation context
  int num_nodes_vol; // number of basis functions
  int num_nodes_surf; // number of basis functions

  // Labatto nodes
  struct gkyl_array *nodes; // local nodal coordinates

  // volume and surface specific nodes
  struct gkyl_array *vol_nodes; // local nodal (volume) coordinates
  struct gkyl_array *surf_nodes_x; // local nodal (surface - x) coordinates (at -1 edge in x)
  struct gkyl_array *surf_nodes_y; // local nodal (surface - y) coordinates (at -1 edge in y)
  struct gkyl_array *surf_nodes_z; // local nodal (surface - z) coordinates (at -1 edge in z)

  dg_gr_maxwell_nodes_c2p_t c2p; // Function transformation comp to phys coords.
  void *c2p_ctx; // Context for the c2p mapping.
};

void
gkyl_surf_and_vol_node_arrays_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_surf_and_vol_node_arrays *up =
    container_of(ref, struct gkyl_surf_and_vol_node_arrays, ref_count);
  if (gkyl_surf_and_vol_node_arrays_is_cu_dev(up)) {
    gkyl_cu_free(up->on_dev);
  } else {
    gkyl_array_release(up->nodal_arr_vol);
    gkyl_array_release(up->nodal_arr_surf_x);
    if (up->ndim > 1) {
      gkyl_array_release(up->nodal_arr_surf_y);
    }
    if (up->ndim > 2) {
      gkyl_array_release(up->nodal_arr_surf_z);
    }
    gkyl_array_release(up->nodal_arr_vol_host);
    gkyl_array_release(up->nodal_arr_surf_x_host);
    if (up->ndim > 1) {
      gkyl_array_release(up->nodal_arr_surf_y_host);
    }
    if (up->ndim > 2) {
      gkyl_array_release(up->nodal_arr_surf_z_host);
    }
  }
  gkyl_free(up);
}

// Identity comp to phys coord mapping, for when user doesn't provide a map.
static inline void
c2p_identity(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_rect_grid *grid = ctx;
  int ndim = grid->ndim;
  for (int d = 0; d < ndim; d++) {
    xphys[d] = xcomp[d];
  }
}

struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *
gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis, int num_ret_vals,
  int polyorder, evalf_t eval, void *ctx
)
{
  return gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew(
    &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp){
      .polyorder = polyorder,
      .grid = grid,
      .basis = basis,
      .num_ret_vals = num_ret_vals,
      .eval = eval,
      .ctx = ctx,
      .c2p_func = 0,
      .c2p_func_ctx = NULL,
    }
  );
}

struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *
gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew(
  const struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp *inp
)
{
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *up =
    gkyl_malloc(sizeof(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes));

  up->grid = *inp->grid;
  up->num_ret_vals = inp->num_ret_vals;
  up->eval = inp->eval;
  up->ctx = inp->ctx;

  // Number of surface and volume nodes
  int ndim = inp->grid->ndim;
  int polyorder = inp->polyorder;
  int num_nodes_vol_per_ndim[3] = {
    (polyorder + 1), (polyorder + 1) * (polyorder + 1),
    (polyorder + 1) * (polyorder + 1) * (polyorder + 1)
  };
  up->num_nodes_vol = num_nodes_vol_per_ndim[ndim - 1];
  int num_nodes_surf_per_ndim[3] = {1, (polyorder + 1), (polyorder + 1) * (polyorder + 1)};
  int num_nodes_surf = num_nodes_surf_per_ndim[ndim - 1];
  up->num_nodes_surf = num_nodes_surf;

  // initialize nodes in local coordinates
  up->nodes = gkyl_array_new(GKYL_DOUBLE, ndim, up->num_nodes_vol);
  up->vol_nodes = gkyl_array_new(GKYL_DOUBLE, ndim, up->num_nodes_vol);
  up->surf_nodes_x = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_surf);
  up->surf_nodes_y = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_surf);
  up->surf_nodes_z = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_surf);

  // Gauss-Labatto Nodes (at corners)
  inp->basis->node_list(gkyl_array_fetch(up->nodes, 0));

  // Gauss-Legandre Nodes for surface and volumes
  vol_node_list(up->vol_nodes, ndim, polyorder);
  surf_node_list(up->surf_nodes_x, ndim, 0, polyorder);
  if (ndim > 1) {
    surf_node_list(up->surf_nodes_y, ndim, 1, polyorder);
  }
  if (ndim > 2) {
    surf_node_list(up->surf_nodes_z, ndim, 2, polyorder);
  }

  if (inp->c2p_func == 0) {
    up->c2p = c2p_identity;
    up->c2p_ctx = &up->grid; // Use grid as the context since all we need is ndim.
  } else {
    up->c2p = inp->c2p_func;
    up->c2p_ctx = inp->c2p_func_ctx;
  }

  return up;
}

struct gkyl_surf_and_vol_node_arrays *
gkyl_surf_and_vol_node_copy_to_device(struct gkyl_surf_and_vol_node_arrays *vol_surf_nodes, int ndim)
{
  struct gkyl_surf_and_vol_node_arrays *up =
    (struct gkyl_surf_and_vol_node_arrays *)gkyl_malloc(sizeof(struct gkyl_surf_and_vol_node_arrays)
    );

  up->use_gpu = true;
  up->ndim = ndim;

  // Copy host-side geometry onto device
  gkyl_array_copy(vol_surf_nodes->nodal_arr_vol, vol_surf_nodes->nodal_arr_vol_host);
  gkyl_array_copy(vol_surf_nodes->nodal_arr_surf_x, vol_surf_nodes->nodal_arr_surf_x_host);
  if (ndim > 1) {
    gkyl_array_copy(vol_surf_nodes->nodal_arr_surf_y, vol_surf_nodes->nodal_arr_surf_y_host);
  }
  if (ndim > 2) {
    gkyl_array_copy(vol_surf_nodes->nodal_arr_surf_z, vol_surf_nodes->nodal_arr_surf_z_host);
  }

  // indirection, pointing to device pointers for copying of device structure
  up->nodal_arr_vol = vol_surf_nodes->nodal_arr_vol->on_dev;
  up->nodal_arr_surf_x = vol_surf_nodes->nodal_arr_surf_x->on_dev;
  if (ndim > 1) {
    up->nodal_arr_surf_y = vol_surf_nodes->nodal_arr_surf_y->on_dev;
  }
  if (ndim > 2) {
    up->nodal_arr_surf_z = vol_surf_nodes->nodal_arr_surf_z->on_dev;
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_surf_and_vol_node_arrays_free);
  up->on_dev = up; // CPU eqn obj points to itself

  // copy the host struct to device struct
  struct gkyl_surf_and_vol_node_arrays *up_cu = (struct gkyl_surf_and_vol_node_arrays *)
    gkyl_cu_malloc(sizeof(struct gkyl_surf_and_vol_node_arrays));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_surf_and_vol_node_arrays), GKYL_CU_MEMCPY_H2D);

  up->on_dev = up_cu; // set the on_dev pointer to the device struct

  return up;
}

struct gkyl_surf_and_vol_node_arrays *
gkyl_surf_and_vol_node_arrays_new(
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *info, long volume, bool use_gpu
)
{
  struct gkyl_surf_and_vol_node_arrays *up =
    gkyl_malloc(sizeof(struct gkyl_surf_and_vol_node_arrays));

  // Allocate the arrays
  up->nodal_arr_vol = mkarr(use_gpu, info->num_ret_vals * info->num_nodes_vol, volume);
  up->nodal_arr_surf_x = mkarr(use_gpu, info->num_ret_vals * info->num_nodes_surf, volume);
  if (info->grid.ndim > 1) {
    up->nodal_arr_surf_y = mkarr(use_gpu, info->num_ret_vals * info->num_nodes_surf, volume);
  }
  if (info->grid.ndim > 2) {
    up->nodal_arr_surf_z = mkarr(use_gpu, info->num_ret_vals * info->num_nodes_surf, volume);
  }

  // Allocate host side arrays
  if (use_gpu) {
    up->nodal_arr_vol_host = mkarr(false, info->num_ret_vals * info->num_nodes_vol, volume);
    up->nodal_arr_surf_x_host = mkarr(false, info->num_ret_vals * info->num_nodes_surf, volume);
    if (info->grid.ndim > 1) {
      up->nodal_arr_surf_y_host = mkarr(false, info->num_ret_vals * info->num_nodes_surf, volume);
    }
    if (info->grid.ndim > 2) {
      up->nodal_arr_surf_z_host = mkarr(false, info->num_ret_vals * info->num_nodes_surf, volume);
    }
  } else {
    up->nodal_arr_vol_host = gkyl_array_acquire(up->nodal_arr_vol);
    up->nodal_arr_surf_x_host = gkyl_array_acquire(up->nodal_arr_surf_x);
    if (info->grid.ndim > 1) {
      up->nodal_arr_surf_y_host = gkyl_array_acquire(up->nodal_arr_surf_y);
    }
    if (info->grid.ndim > 2) {
      up->nodal_arr_surf_z_host = gkyl_array_acquire(up->nodal_arr_surf_z);
    }
  }

  up->use_gpu = use_gpu;
  up->ndim = info->grid.ndim;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_surf_and_vol_node_arrays_free);
  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

static inline void
log_to_comp(
  int ndim, const double *eta, const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xc,
  double *GKYL_RESTRICT xout
)
{
  for (int d = 0; d < ndim; ++d) {
    xout[d] = 0.5 * dx[d] * eta[d] + xc[d];
  }
}

static inline void
copy_double_arr(int n, const double *GKYL_RESTRICT inp, double *GKYL_RESTRICT out)
{
  for (int i = 0; i < n; ++i) {
    out[i] = inp[i];
  }
}

double *
gkyl_dg_gr_maxwell_surf_and_vol_nodes_fetch_node(
  const struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *up, long node
)
{
  return gkyl_array_fetch(up->nodes, node);
}

static void
gkyl_copy_to_nodal_array(
  const struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *up, const struct gkyl_array *fun_at_nodes,
  int num_nodes, double *f
)
{
  const double *fao = gkyl_array_cfetch(fun_at_nodes, 0); // pointer to values at nodes

  int num_ret_vals = up->num_ret_vals;
  for (int i = 0; i < num_ret_vals; ++i) {
    // Note this is swapped from eval_on_nodes and has the memory pattern of:
    // c0[0], c0[1], ... c1[0], c1[1], .... where c0 is the component and 0...num_nodes are the nodes.
    for (int k = 0; k < num_nodes; ++k) {
      f[k + num_nodes * i] = fao[num_ret_vals * k + i];
    }
  }
}

void
gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(
  const struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *up, double tm,
  const struct gkyl_range *update_range, struct gkyl_surf_and_vol_node_arrays *surf_vol_nodal_arrays
)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  int num_nodes_vol = up->num_nodes_vol;
  int num_nodes_surf = up->num_nodes_surf;
  struct gkyl_array *fun_at_nodes_vol = surf_vol_nodal_arrays->nodal_arr_vol_host;
  struct gkyl_array *fun_at_nodes_surf_x = surf_vol_nodal_arrays->nodal_arr_surf_x_host;
  struct gkyl_array *fun_at_nodes_surf_y = surf_vol_nodal_arrays->nodal_arr_surf_y_host;
  struct gkyl_array *fun_at_nodes_surf_z = surf_vol_nodal_arrays->nodal_arr_surf_z_host;

  // Create a cell per node:
  struct gkyl_array *vol_func = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, num_nodes_vol);
  struct gkyl_array *surf_func = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, num_nodes_surf);

  // Volume
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i = 0; i < num_nodes_vol; ++i) {
      log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->vol_nodes, i), up->grid.dx, xc, xmu);
      up->c2p(xmu, xmu, up->c2p_ctx);
      up->eval(tm, xmu, gkyl_array_fetch(vol_func, i), up->ctx);
    }

    // Copy local array to the nodal array
    long lidx = gkyl_range_idx(update_range, iter.idx);
    gkyl_copy_to_nodal_array(up, vol_func, num_nodes_vol, gkyl_array_fetch(fun_at_nodes_vol, lidx));
  }

  // x-surface
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i = 0; i < num_nodes_surf; ++i) {
      log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->surf_nodes_x, i), up->grid.dx, xc, xmu);
      up->c2p(xmu, xmu, up->c2p_ctx);
      up->eval(tm, xmu, gkyl_array_fetch(surf_func, i), up->ctx);
    }

    // Copy local array to the nodal array
    long lidx = gkyl_range_idx(update_range, iter.idx);
    gkyl_copy_to_nodal_array(
      up, surf_func, num_nodes_surf, gkyl_array_fetch(fun_at_nodes_surf_x, lidx)
    );
  }

  // y-surface
  if (up->grid.ndim > 1) {
    gkyl_range_iter_init(&iter, update_range);
    while (gkyl_range_iter_next(&iter)) {
      gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

      for (int i = 0; i < num_nodes_surf; ++i) {
        log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->surf_nodes_y, i), up->grid.dx, xc, xmu);
        up->c2p(xmu, xmu, up->c2p_ctx);
        up->eval(tm, xmu, gkyl_array_fetch(surf_func, i), up->ctx);
      }

      // Copy local array to the nodal array
      long lidx = gkyl_range_idx(update_range, iter.idx);
      gkyl_copy_to_nodal_array(
        up, surf_func, num_nodes_surf, gkyl_array_fetch(fun_at_nodes_surf_y, lidx)
      );
    }
  }

  // z-surface
  if (up->grid.ndim > 2) {
    gkyl_range_iter_init(&iter, update_range);
    while (gkyl_range_iter_next(&iter)) {
      gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

      for (int i = 0; i < num_nodes_surf; ++i) {
        log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->surf_nodes_z, i), up->grid.dx, xc, xmu);
        up->c2p(xmu, xmu, up->c2p_ctx);
        up->eval(tm, xmu, gkyl_array_fetch(surf_func, i), up->ctx);
      }

      // Copy local array to the nodal array
      long lidx = gkyl_range_idx(update_range, iter.idx);
      gkyl_copy_to_nodal_array(
        up, surf_func, num_nodes_surf, gkyl_array_fetch(fun_at_nodes_surf_z, lidx)
      );
    }
  }

  // Copy to host if using gpu
  if (surf_vol_nodal_arrays->use_gpu) {
    gkyl_array_copy(surf_vol_nodal_arrays->nodal_arr_vol, surf_vol_nodal_arrays->nodal_arr_vol_host);
    gkyl_array_copy(
      surf_vol_nodal_arrays->nodal_arr_surf_x, surf_vol_nodal_arrays->nodal_arr_surf_x_host
    );
    if (up->grid.ndim > 1) {
      gkyl_array_copy(
        surf_vol_nodal_arrays->nodal_arr_surf_y, surf_vol_nodal_arrays->nodal_arr_surf_y_host
      );
    }
    if (up->grid.ndim > 2) {
      gkyl_array_copy(
        surf_vol_nodal_arrays->nodal_arr_surf_z, surf_vol_nodal_arrays->nodal_arr_surf_z_host
      );
    }
  }

  // Release functions per node
  gkyl_array_release(vol_func);
  gkyl_array_release(surf_func);
}

void
gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *up)
{
  gkyl_array_release(up->nodes);
  gkyl_array_release(up->vol_nodes);
  gkyl_array_release(up->surf_nodes_x);
  gkyl_array_release(up->surf_nodes_y);
  gkyl_array_release(up->surf_nodes_z);
  gkyl_free(up);
}

bool
gkyl_surf_and_vol_node_arrays_is_cu_dev(const struct gkyl_surf_and_vol_node_arrays *vol_surf_nodes)
{
  return GKYL_IS_CU_ALLOC(vol_surf_nodes->flags);
}

struct gkyl_surf_and_vol_node_arrays *
gkyl_surf_and_vol_node_arrays_acquire(const struct gkyl_surf_and_vol_node_arrays *vol_surf_nodes)
{
  gkyl_ref_count_inc(&vol_surf_nodes->ref_count);
  return (struct gkyl_surf_and_vol_node_arrays *)vol_surf_nodes;
}

void
gkyl_surf_and_vol_node_arrays_release(const struct gkyl_surf_and_vol_node_arrays *vol_surf_nodes)
{
  gkyl_ref_count_dec(&vol_surf_nodes->ref_count);
}
