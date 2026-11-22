/* -*- c++ -*- */

extern "C" {
#include <gkyl_dg_calc_resistive_layer.h>
#include <gkyl_dg_calc_resistive_layer_priv.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_dg_calc_resistive_layer_advance_cu_ker(
  struct gkyl_range range, int num_basis, int tot_quad,
  const struct gkyl_array *GKYL_RESTRICT weights,
  const struct gkyl_array *GKYL_RESTRICT basis_at_ords,
  const struct gkyl_array *GKYL_RESTRICT sigma, const struct gkyl_array *GKYL_RESTRICT em,
  struct gkyl_array *GKYL_RESTRICT rhs
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long loc = gkyl_range_idx(&range, idx);
    dg_calc_resistive_layer_cell(
      num_basis, tot_quad, (const double *)weights->data, (const double *)basis_at_ords->data,
      (const double *)gkyl_array_cfetch(sigma, loc), (const double *)gkyl_array_cfetch(em, loc),
      (double *)gkyl_array_fetch(rhs, loc)
    );
  }
}

void
gkyl_dg_calc_resistive_layer_advance_cu(
  const struct gkyl_dg_calc_resistive_layer *up, const struct gkyl_range *range,
  const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  int nblocks = range->nblocks, nthreads = range->nthreads;
  gkyl_dg_calc_resistive_layer_advance_cu_ker<<<nblocks, nthreads>>>(
    *range, up->num_basis, up->tot_quad, up->weights->on_dev, up->basis_at_ords->on_dev,
    up->sigma->on_dev, em->on_dev, rhs->on_dev
  );
}
