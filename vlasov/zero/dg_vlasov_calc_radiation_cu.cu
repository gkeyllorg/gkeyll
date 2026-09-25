/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_calc_radiation.h>
#include <gkyl_dg_vlasov_calc_radiation_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_vlasov_calc_radiation_cu_kernel(
  struct gkyl_rect_grid vel_grid, struct gkyl_basis vel_basis, struct gkyl_range vel_range,
  enum gkyl_vlasov_radiation_id radiation_id, const struct gkyl_array *vmap, double t_cool,
  double p0, struct gkyl_array *rad
)
{
  int vdim = vel_basis.ndim;
  int poly_order = vel_basis.poly_order;
  calc_radiation_t calc_radiation;
  switch (vel_basis.b_type) {
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

  int idx[GKYL_MAX_DIM];
  // Cell center array
  double xc[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < vel_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&vel_range, linc1, idx);
    gkyl_rect_grid_cell_center(&vel_grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&vel_range, idx);

    double *rad_d = (double *)gkyl_array_fetch(rad, loc);
    calc_radiation(
      xc, vel_grid.dx, vmap ? (const double *)gkyl_array_cfetch(vmap, loc) : 0, t_cool, p0, rad_d
    );
  }
}

// Host-side wrapper for initialization of radiation drag force.
void
gkyl_dg_vlasov_calc_radiation_cu(
  const struct gkyl_rect_grid *vel_grid, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *vel_range, enum gkyl_vlasov_radiation_id radiation_id,
  const struct gkyl_vlasov_velocity_map *vel_map, double t_cool, double p0, struct gkyl_array *rad
)
{
  int nblocks = vel_range->nblocks;
  int nthreads = vel_range->nthreads;
  // The velocity map is required: the kernels always evaluate the velocity
  // coordinate from the stored map (identity map for uniform grids).
  assert(vel_map);
  gkyl_dg_vlasov_calc_radiation_cu_kernel<<<nblocks, nthreads>>>(
    *vel_grid, *vel_basis, *vel_range, radiation_id, vel_map->vmap->on_dev, t_cool, p0, rad->on_dev
  );
}