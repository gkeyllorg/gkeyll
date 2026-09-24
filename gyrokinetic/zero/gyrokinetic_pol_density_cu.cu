/* -*- c++ -*- */

extern "C" {
#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_gyrokinetic_pol_density_priv.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_gk_pol_den_set_cu_ker_ptrs(struct gkyl_gyrokinetic_pol_density_kernels *kernels,
  struct gkyl_basis cbasis)
{
  int pdim = cbasis.ndim;
  enum gkyl_basis_type b_type = cbasis.b_type;
  int poly_order = cbasis.poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->pol_den = gk_pol_density_kern_list_ser[pdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gk_pol_den_choose_kernel_cu(struct gkyl_gyrokinetic_pol_density_kernels *kernels,
  struct gkyl_basis cbasis)
{
  gkyl_gk_pol_den_set_cu_ker_ptrs<<<1,1>>>(kernels, cbasis);
}

__global__ static void
gkyl_gyrokinetic_pol_density_advance_cu_ker(
  struct gkyl_gyrokinetic_pol_density_kernels *kers,
  const struct gkyl_rect_grid grid, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT pol_weight, const struct gkyl_array* GKYL_RESTRICT phi,
  struct gkyl_array* GKYL_RESTRICT npol)
{
  int cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long linidx = gkyl_range_idx(&conf_range, cidx);
    const double *pol_weight_d = (const double*) gkyl_array_cfetch(pol_weight, linidx);
    const double *phi_d = (const double*) gkyl_array_cfetch(phi, linidx);
    double *npol_d = (double*) gkyl_array_fetch(npol, linidx);

    // Compute the polarization density.
    kers->pol_den(grid.dx, pol_weight_d, phi_d, npol_d);
  }
}

void
gkyl_gyrokinetic_pol_density_advance_cu(gkyl_gyrokinetic_pol_density* up,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT pol_weight,
  const struct gkyl_array *GKYL_RESTRICT phi, struct gkyl_array *GKYL_RESTRICT npol)
{
  int nblocks = conf_rng->nblocks, nthreads = conf_rng->nthreads;

  gkyl_gyrokinetic_pol_density_advance_cu_ker<<<nblocks, nthreads>>>
    (up->kernels, up->grid, *conf_rng, pol_weight->on_dev, phi->on_dev, npol->on_dev);
}
