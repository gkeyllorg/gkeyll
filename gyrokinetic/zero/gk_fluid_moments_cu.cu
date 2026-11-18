extern "C" {
#include <gkyl_gk_fluid_moments_priv.h>
}

__global__ static void
gk_fluid_moments_final_integrate_kernel(
  int nb, int nm, double rmass_fac, struct gkyl_range range, const struct gkyl_array *result,
  struct gkyl_array *out
)
{
  int idx[GKYL_MAX_DIM];
  for (long cell = threadIdx.x + blockIdx.x * blockDim.x; cell < range.volume;
       cell += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&range, cell, idx);
    long loc = gkyl_range_idx(&range, idx);
    const double *dg = (const double *)gkyl_array_cfetch(result, loc);
    double *integral = (double *)gkyl_array_fetch(out, loc);
    for (int c = 0; c < nm; ++c) {
      integral[c] = rmass_fac * dg[c * nb];
    }
  }
}

void
gkyl_gk_fluid_moments_final_integrate_cu(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range, struct gkyl_array *out
)
{
  gk_fluid_moments_final_integrate_kernel<<<range->nblocks, range->nthreads>>>(
    up->basis.num_basis, up->num_mom_out, up->volnorm_fac, *range, up->result->on_dev, out->on_dev
  );
}
