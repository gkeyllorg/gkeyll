#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_gr_maxwell_geom_source.h>
#include <gkyl_dg_gr_maxwell_geom_source_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>

gkyl_dg_gr_maxwell_geom_source*
gkyl_dg_gr_maxwell_geom_source_inew(const struct gkyl_dg_gr_maxwell_geom_source_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_gr_maxwell_geom_source_cu_dev_inew(inp);
  }
#endif

  struct gkyl_dg_gr_maxwell_geom_source *up = gkyl_malloc(sizeof(*up));

  up->cdim = inp->conf_basis->ndim;
  up->use_gpu = inp->use_gpu;
  up->conf_grid = *inp->conf_grid;
  up->gr_maxwell_data.chi = inp->chi;
  up->gr_maxwell_data.gamma = inp->gamma;
  up->gr_maxwell_data.K_phi = inp->K_phi;
  up->gr_maxwell_data.K_psi = inp->K_psi;

  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;

  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->geom_source = ser_geom_source_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      up->geom_source = ten_geom_source_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }

  // ensure non-NULL pointers
  assert(up->geom_source);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_gr_maxwell_geom_source_advance(struct gkyl_dg_gr_maxwell_geom_source *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *geom_factor_con,
  const struct gkyl_array *field_con, struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_gr_maxwell_geom_source_advance_cu(up,
      conf_range, geom_factor_con, field_con, rhs);
  }
#endif

  assert(up->geom_source);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->conf_grid, iter.idx, xc);

    const double *geom_factor_con_d = gkyl_array_cfetch(geom_factor_con->nodal_arr_vol, cidx);
    const double *field_con_d = gkyl_array_cfetch(field_con, cidx);
    double *rhs_d = gkyl_array_fetch(rhs, cidx);

    up->geom_source(&up->gr_maxwell_data, xc, up->conf_grid.dx,
      geom_factor_con_d, field_con_d, rhs_d);
  }
}

void
gkyl_dg_gr_maxwell_geom_source_release(struct gkyl_dg_gr_maxwell_geom_source *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  gkyl_free(up);
}
