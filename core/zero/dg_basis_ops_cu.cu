/* -*- c++ -*- */

extern "C" {
#include <gkyl_util.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_dg_basis_ops_priv.h>
#include <gkyl_range.h>
#include <float.h>
}

__global__ void
gkyl_dg_basis_ops_eval_array_at_coord_comp_cu_ker_eval(const struct gkyl_array *arr, const double *coord,
  const struct gkyl_basis *basis, struct gkyl_rect_grid grid, struct gkyl_range rng,
  double *out)
{
  int coord_idx[GKYL_MAX_DIM];
  gkyl_rect_grid_coord_idx(&grid, coord, coord_idx);

  if (gkyl_range_contains_idx(&rng, coord_idx)) {
    double xc[GKYL_MAX_DIM], coord_log[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, coord_idx, xc);
    for (int d=0; d<grid.ndim; d++)
      coord_log[d] = (2.0/grid.dx[d])*(coord[d]-xc[d]);
  
    long linidx = gkyl_range_idx(&rng, coord_idx);
    const double *arr_c = (const double*) gkyl_array_cfetch(arr, linidx);
    out[0] = basis->eval_expand(coord_log, arr_c);
  }
  else
    out[0] = -DBL_MAX;
}

__global__ void
gkyl_dg_basis_ops_eval_array_at_coord_comp_cu_ker_none(const struct gkyl_array *arr, const double *coord,
  const struct gkyl_basis *basis, struct gkyl_rect_grid grid, struct gkyl_range rng,
  double *out)
{
  out[0] = -DBL_MAX;
}

void
gkyl_dg_basis_ops_eval_array_at_coord_comp_cu(const struct gkyl_array *arr, const double *coord,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, const struct gkyl_range *rng,
  double *out)
{
  gkyl_dg_basis_ops_eval_array_at_coord_comp_cu_ker_eval<<<1,1>>>(arr->on_dev,
    coord, basis, *grid, *rng, out);
}
