#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_calc_radiation.h>
#include <gkyl_dg_vlasov_calc_radiation_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_vlasov_calc_radiation(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_vlasov_radiation_id radiation_id,
  const struct gkyl_vlasov_velocity_map *vel_map, double t_cool, double p0, struct gkyl_array *rad,
  bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_dg_vlasov_calc_radiation_cu(
      vel_grid, vel_basis, vel_range, radiation_id, vel_map, t_cool, p0, rad
    );
    return;
  }
#endif

  // The velocity map is required: the kernels always evaluate the velocity
  // coordinate from the stored map (identity map for uniform grids).
  assert(vel_map);
  const struct gkyl_array *vmap = vel_map->vmap;

  int vdim = vel_basis->ndim;
  int poly_order = vel_basis->poly_order;
  calc_radiation_t calc_radiation;
  switch (vel_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      calc_radiation = choose_ser_radiation_kern(radiation_id, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      calc_radiation = choose_tensor_radiation_kern(radiation_id, vdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, vel_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(vel_grid, iter.idx, xc);
    long vidx = gkyl_range_idx(vel_range, iter.idx);

    double *rad_d = gkyl_array_fetch(rad, vidx);
    calc_radiation(xc, vel_grid->dx, vmap ? gkyl_array_cfetch(vmap, vidx) : 0, t_cool, p0, rad_d);
  }
}