#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov_calc_hamil_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_vlasov_calc_hamil(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_model_id model_id,
  const struct gkyl_vlasov_velocity_map *vel_map, struct gkyl_array *hamil,
  struct gkyl_array *hamil_inv, bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_dg_vlasov_calc_hamil_cu(
      vel_grid, vel_basis, vel_range, model_id, vel_map, hamil, hamil_inv
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
  calc_hamil_t calc_hamil;
  switch (vel_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      calc_hamil = choose_ser_hamil_kern(model_id, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      calc_hamil = choose_tensor_hamil_kern(model_id, vdim, poly_order);
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

    double *hamil_d = gkyl_array_fetch(hamil, vidx);
    double *hamil_inv_d = gkyl_array_fetch(hamil_inv, vidx);
    calc_hamil(xc, vel_grid->dx, vmap ? gkyl_array_cfetch(vmap, vidx) : 0, hamil_d, hamil_inv_d);
  }
}