extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_gyrokinetic_omega_h_priv.h>
}

/**
 * Reduce sampled Jn*G^2 across cells, with one nonnegative atomic maximum per block.
 *
 * @param range Interior range indexing both density and cached weights.
 * @param basis_at_nodes Device basis evaluations at all samples.
 * @param weight_sq Device cached G^2 values.
 * @param density Device modal Jn coefficients.
 * @param max_out Device scalar initialized to zero before this kernel is launched.
 */
__global__ static void
gyrokinetic_omega_h_max_kernel(
  struct gkyl_range range, const struct gkyl_array *basis_at_nodes,
  const struct gkyl_array *weight_sq, const struct gkyl_array *density, double *max_out
)
{
  __shared__ double block_max[GKYL_DEFAULT_NUM_THREADS];
  double max_value = 0.0;
  for (long cell = threadIdx.x + blockIdx.x * blockDim.x; cell < range.volume;
       cell += blockDim.x * gridDim.x) {
    int idx[GKYL_MAX_DIM];
    gkyl_sub_range_inv_idx(&range, cell, idx);
    long loc = gkyl_range_idx(&range, idx);
    max_value = fmax(
      max_value, gyrokinetic_omega_h_cell_max(
                   basis_at_nodes, (const double *)gkyl_array_cfetch(weight_sq, loc),
                   (const double *)gkyl_array_cfetch(density, loc)
                 )
    );
  }
  block_max[threadIdx.x] = max_value;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (threadIdx.x < stride) {
      block_max[threadIdx.x] = fmax(block_max[threadIdx.x], block_max[threadIdx.x + stride]);
    }
    __syncthreads();
  }
  // Nonnegative IEEE doubles have the same ordering as their unsigned bit patterns.
  if (threadIdx.x == 0) {
    atomicMax((unsigned long long *)max_out, (unsigned long long)__double_as_longlong(block_max[0]));
  }
}

double
gkyl_gyrokinetic_omega_h_advance_cu(
  struct gkyl_gyrokinetic_omega_h *up, const struct gkyl_array *density
)
{
  gkyl_cu_memset(up->max_dev, 0, sizeof(double));
  if (up->range.volume > 0) {
    int nblocks =
      GKYL_MIN2(256, (up->range.volume + GKYL_DEFAULT_NUM_THREADS - 1) / GKYL_DEFAULT_NUM_THREADS);
    gyrokinetic_omega_h_max_kernel<<<nblocks, GKYL_DEFAULT_NUM_THREADS>>>(
      up->range, up->basis_at_nodes->on_dev, up->weight_sq->on_dev, density->on_dev, up->max_dev
    );
  }
  double max_value;
  gkyl_cu_memcpy(&max_value, up->max_dev, sizeof(double), GKYL_CU_MEMCPY_D2H);
  return max_value;
}
