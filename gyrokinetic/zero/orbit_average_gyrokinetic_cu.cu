/* -*- c++ -*- */
#include <float.h>
#include <math.h>
#include <stdint.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_orbit_average_gyrokinetic_priv.h>
}

enum { GKYL_ORBIT_CU_THREADS = 128 };

// The host descriptor owns device buffers. Only the three final diagnostics
// cross to the host during an update; the distribution stays on the device.
struct gk_orbit_cu {
  struct gkyl_orbit_average_gyrokinetic_block *blocks;
  long *cell_indices;
  int *cell_blocks, *reduced_blocks;
  double *q, *dual, *eval, *coeff, *correction, *blend;
  double *block_invalid, *block_limited, *status;
};

// Avoid arithmetic on infinity in the host wrapper, whose flags may include
// -ffast-math. The public interface explicitly supports the projection limit.
static double
gk_orbit_beta_cu(double dt, double frequency)
{
  volatile union {
    double value;
    uint64_t bits;
  } rate = {frequency};
  if (rate.bits == UINT64_C(0x7ff0000000000000) || (dt > 1.0 && frequency > DBL_MAX / dt)) {
    return 1.0;
  }
  const double product = dt * frequency;
  return product >= 1.0 ? 1.0 - 1.0 / (1.0 + product) : product / (1.0 + product);
}

// The bit test remains valid with the application's fast-math compiler flags.
__device__ static bool
gk_orbit_isfinite_cu(double value)
{
  return (__double_as_longlong(value) & 0x7ff0000000000000LL) != 0x7ff0000000000000LL;
}

__global__ static void
gk_orbit_reduce_cu_kernel(
  int num_basis, long num_reduced, const struct gkyl_orbit_average_gyrokinetic_block *blocks,
  const long *cell_indices, const int *reduced_blocks, const double *dual,
  const struct gkyl_array *fin, double *coeff
)
{
  __shared__ double buffer[GKYL_ORBIT_CU_THREADS];
  for (long reduced = blockIdx.x; reduced < num_reduced; reduced += gridDim.x) {
    const struct gkyl_orbit_average_gyrokinetic_block block = blocks[reduced_blocks[reduced]];
    const int column = reduced - block.reduce_offset;
    const double *dual_column = dual + block.basis_offset + (long)column * block.num_modes;
    double value = 0.0;
    for (int mode = threadIdx.x; mode < block.num_modes; mode += blockDim.x) {
      const long cell = cell_indices[block.cell_offset + mode / num_basis];
      const double *fin_d = (const double *)gkyl_array_cfetch(fin, cell);
      value += dual_column[mode] * fin_d[mode % num_basis];
    }
    buffer[threadIdx.x] = value;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
      if (threadIdx.x < stride) {
        buffer[threadIdx.x] += buffer[threadIdx.x + stride];
      }
      __syncthreads();
    }
    if (threadIdx.x == 0) {
      coeff[reduced] = buffer[0];
    }
    __syncthreads();
  }
}

__global__ static void
gk_orbit_correction_cu_kernel(
  int num_basis, long num_modes, const struct gkyl_orbit_average_gyrokinetic_block *blocks,
  const long *cell_indices, const int *cell_blocks, const double *q, const double *coeff,
  const struct gkyl_array *fin, double *correction
)
{
  for (long mode = threadIdx.x + (long)blockIdx.x * blockDim.x; mode < num_modes;
       mode += (long)blockDim.x * gridDim.x) {
    const long cell = mode / num_basis;
    const struct gkyl_orbit_average_gyrokinetic_block block = blocks[cell_blocks[cell]];
    const int block_mode = mode - block.cell_offset * num_basis;
    double value = 0.0;
    for (int column = 0; column < block.rank; ++column) {
      value += q[block.basis_offset + (long)column * block.num_modes + block_mode] *
               coeff[block.reduce_offset + column];
    }
    const double *fin_d = (const double *)gkyl_array_cfetch(fin, cell_indices[cell]);
    correction[mode] = value - fin_d[mode % num_basis];
  }
}

__global__ static void
gk_orbit_limit_cu_kernel(
  int num_basis, int num_pos, int num_blocks, bool positivity, double beta,
  const struct gkyl_orbit_average_gyrokinetic_block *blocks, const long *cell_indices,
  const double *eval, const double *correction, const struct gkyl_array *fin, double *blend,
  double *block_invalid, double *block_limited
)
{
  __shared__ double bounds[GKYL_ORBIT_CU_THREADS];
  __shared__ int invalid[GKYL_ORBIT_CU_THREADS];
  for (int block_id = blockIdx.x; block_id < num_blocks; block_id += gridDim.x) {
    const struct gkyl_orbit_average_gyrokinetic_block block = blocks[block_id];
    double bound = INFINITY;
    int bad = 0;
    // Check modal data even when positivity limiting is disabled. Nonfinite
    // data must not permit a partial update of otherwise valid blocks.
    for (int mode = threadIdx.x; mode < block.num_modes; mode += blockDim.x) {
      const long cell = block.cell_offset + mode / num_basis;
      const double *fin_d = (const double *)gkyl_array_cfetch(fin, cell_indices[cell]);
      if (!gk_orbit_isfinite_cu(fin_d[mode % num_basis]) ||
          !gk_orbit_isfinite_cu(correction[cell * num_basis + mode % num_basis])) {
        bad = 1;
      }
    }
    if (positivity) {
      for (long node = threadIdx.x; node < (long)block.num_cells * num_pos; node += blockDim.x) {
        const long cell = block.cell_offset + node / num_pos;
        const double *fin_d = (const double *)gkyl_array_cfetch(fin, cell_indices[cell]);
        const double *delta_d = correction + cell * num_basis;
        const double *eval_d = eval + (node % num_pos) * num_basis;
        double value = 0.0, delta = 0.0;
        for (int mode = 0; mode < num_basis; ++mode) {
          value += eval_d[mode] * fin_d[mode];
          delta += eval_d[mode] * delta_d[mode];
        }
        if (!gk_orbit_isfinite_cu(value) || !gk_orbit_isfinite_cu(delta) || value < 0.0) {
          bad = 1;
        }
        if (delta < 0.0) {
          bound = fmin(bound, value / -delta);
        }
      }
    }
    bounds[threadIdx.x] = bound;
    invalid[threadIdx.x] = bad;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
      if (threadIdx.x < stride) {
        bounds[threadIdx.x] = fmin(bounds[threadIdx.x], bounds[threadIdx.x + stride]);
        invalid[threadIdx.x] |= invalid[threadIdx.x + stride];
      }
      __syncthreads();
    }
    if (threadIdx.x == 0) {
      const bool limited = bounds[0] < beta;
      blend[block_id] = limited ? fmax(0.0, 0.99 * bounds[0]) : beta;
      block_invalid[block_id] = invalid[0];
      block_limited[block_id] = limited;
    }
    __syncthreads();
  }
}

__global__ static void
gk_orbit_status_cu_kernel(
  int num_blocks, double beta, const double *blend, const double *block_invalid,
  const double *block_limited, double *status
)
{
  // Only scalar diagnostics are reduced here, after the parallel node checks.
  double invalid = 0.0, min_blend = beta, num_limited = 0.0;
  for (int block_id = 0; block_id < num_blocks; ++block_id) {
    invalid = fmax(invalid, block_invalid[block_id]);
    min_blend = fmin(min_blend, blend[block_id]);
    num_limited += block_limited[block_id];
  }
  status[0] = invalid;
  status[1] = min_blend;
  status[2] = num_limited;
}

__global__ static void
gk_orbit_update_cu_kernel(
  int num_basis, long num_modes, const long *cell_indices, const int *cell_blocks,
  const double *correction, const double *blend, struct gkyl_array *f
)
{
  for (long mode = threadIdx.x + (long)blockIdx.x * blockDim.x; mode < num_modes;
       mode += (long)blockDim.x * gridDim.x) {
    const long cell = mode / num_basis;
    const double beta = blend[cell_blocks[cell]];
    double *f_d = (double *)gkyl_array_fetch(f, cell_indices[cell]);
    f_d[mode % num_basis] += beta * correction[mode];
  }
}

void
gkyl_orbit_average_gyrokinetic_cu_init(gkyl_orbit_average_gyrokinetic *up)
{
  struct gk_orbit_cu *gpu = (struct gk_orbit_cu *)gkyl_calloc(1, sizeof(*gpu));
  up->gpu = gpu;
  if (up->num_cells == 0) {
    return;
  }

  const size_t block_bytes = up->num_blocks * sizeof(*up->blocks);
  gpu->blocks = (struct gkyl_orbit_average_gyrokinetic_block *)gkyl_cu_malloc(block_bytes);
  gkyl_cu_memcpy(gpu->blocks, up->blocks, block_bytes, GKYL_CU_MEMCPY_H2D);
  gpu->cell_indices = (long *)gkyl_cu_malloc(up->num_cells * sizeof(long));
  gkyl_cu_memcpy(
    gpu->cell_indices, up->cell_indices, up->num_cells * sizeof(long), GKYL_CU_MEMCPY_H2D
  );

  int *cell_blocks = (int *)gkyl_malloc(up->num_cells * sizeof(int));
  int *reduced_blocks = (int *)gkyl_malloc(up->num_reduced * sizeof(int));
  for (int block_id = 0; block_id < up->num_blocks; ++block_id) {
    const struct gkyl_orbit_average_gyrokinetic_block *block = up->blocks + block_id;
    for (int cell = 0; cell < block->num_cells; ++cell) {
      cell_blocks[block->cell_offset + cell] = block_id;
    }
    for (int column = 0; column < block->rank; ++column) {
      reduced_blocks[block->reduce_offset + column] = block_id;
    }
  }
  gpu->cell_blocks = (int *)gkyl_cu_malloc(up->num_cells * sizeof(int));
  gpu->reduced_blocks = (int *)gkyl_cu_malloc(up->num_reduced * sizeof(int));
  gkyl_cu_memcpy(gpu->cell_blocks, cell_blocks, up->num_cells * sizeof(int), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(
    gpu->reduced_blocks, reduced_blocks, up->num_reduced * sizeof(int), GKYL_CU_MEMCPY_H2D
  );
  gkyl_free(cell_blocks);
  gkyl_free(reduced_blocks);

  gpu->q = (double *)gkyl_cu_malloc(up->basis_size * sizeof(double));
  gpu->dual = (double *)gkyl_cu_malloc(up->basis_size * sizeof(double));
  gkyl_cu_memcpy(gpu->q, up->q, up->basis_size * sizeof(double), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(gpu->dual, up->dual, up->basis_size * sizeof(double), GKYL_CU_MEMCPY_H2D);
  const size_t eval_bytes = up->num_pos * up->num_basis * sizeof(double);
  gpu->eval = (double *)gkyl_cu_malloc(eval_bytes);
  gkyl_cu_memcpy(gpu->eval, up->eval, eval_bytes, GKYL_CU_MEMCPY_H2D);
  gpu->coeff = (double *)gkyl_cu_malloc(up->num_reduced * sizeof(double));
  gpu->correction = (double *)gkyl_cu_malloc(up->num_cells * up->num_basis * sizeof(double));
  gpu->blend = (double *)gkyl_cu_malloc(up->num_blocks * sizeof(double));
  gpu->block_invalid = (double *)gkyl_cu_malloc(up->num_blocks * sizeof(double));
  gpu->block_limited = (double *)gkyl_cu_malloc(up->num_blocks * sizeof(double));
  gpu->status = (double *)gkyl_cu_malloc(3 * sizeof(double));
}

struct gkyl_orbit_average_gyrokinetic_result
gkyl_orbit_average_gyrokinetic_advance_cu(
  gkyl_orbit_average_gyrokinetic *up, double dt, double frequency, struct gkyl_array *f
)
{
  struct gkyl_orbit_average_gyrokinetic_result result = {true, 0.0, 0};
  if (dt == 0.0 || frequency == 0.0 || up->num_cells == 0) {
    return result;
  }
  const double beta = gk_orbit_beta_cu(dt, frequency);
  const struct gk_orbit_cu *gpu = (const struct gk_orbit_cu *)up->gpu;
  const long num_modes = up->num_cells * up->num_basis;
  const int mode_blocks = (int)fmin(65535.0, (double)((num_modes - 1) / GKYL_ORBIT_CU_THREADS + 1));
  const int reduced_blocks = (int)fmin(65535.0, (double)up->num_reduced);
  const int limit_blocks = up->num_blocks < 65535 ? up->num_blocks : 65535;
  gk_orbit_reduce_cu_kernel<<<reduced_blocks, GKYL_ORBIT_CU_THREADS>>>(
    up->num_basis, up->num_reduced, gpu->blocks, gpu->cell_indices, gpu->reduced_blocks, gpu->dual,
    f->on_dev, gpu->coeff
  );
  gk_orbit_correction_cu_kernel<<<mode_blocks, GKYL_ORBIT_CU_THREADS>>>(
    up->num_basis, num_modes, gpu->blocks, gpu->cell_indices, gpu->cell_blocks, gpu->q, gpu->coeff,
    f->on_dev, gpu->correction
  );
  gk_orbit_limit_cu_kernel<<<limit_blocks, GKYL_ORBIT_CU_THREADS>>>(
    up->num_basis, up->num_pos, up->num_blocks, up->positivity, beta, gpu->blocks,
    gpu->cell_indices, gpu->eval, gpu->correction, f->on_dev, gpu->blend, gpu->block_invalid,
    gpu->block_limited
  );
  gk_orbit_status_cu_kernel<<<1, 1>>>(
    up->num_blocks, beta, gpu->blend, gpu->block_invalid, gpu->block_limited, gpu->status
  );
  double status[3];
  gkyl_cu_memcpy(status, gpu->status, sizeof(status), GKYL_CU_MEMCPY_D2H);
  result.success = status[0] == 0.0;
  result.min_blend = status[1];
  result.num_limited_blocks = (long)status[2];
  if (result.success) {
    gk_orbit_update_cu_kernel<<<mode_blocks, GKYL_ORBIT_CU_THREADS>>>(
      up->num_basis, num_modes, gpu->cell_indices, gpu->cell_blocks, gpu->correction, gpu->blend,
      f->on_dev
    );
  }
  return result;
}

void
gkyl_orbit_average_gyrokinetic_cu_release(gkyl_orbit_average_gyrokinetic *up)
{
  struct gk_orbit_cu *gpu = (struct gk_orbit_cu *)up->gpu;
  if (gpu == 0) {
    return;
  }
  if (up->num_cells > 0) {
    gkyl_cu_free(gpu->blocks);
    gkyl_cu_free(gpu->cell_indices);
    gkyl_cu_free(gpu->cell_blocks);
    gkyl_cu_free(gpu->reduced_blocks);
    gkyl_cu_free(gpu->q);
    gkyl_cu_free(gpu->dual);
    gkyl_cu_free(gpu->eval);
    gkyl_cu_free(gpu->coeff);
    gkyl_cu_free(gpu->correction);
    gkyl_cu_free(gpu->blend);
    gkyl_cu_free(gpu->block_invalid);
    gkyl_cu_free(gpu->block_limited);
    gkyl_cu_free(gpu->status);
  }
  gkyl_free(gpu);
  up->gpu = 0;
}
