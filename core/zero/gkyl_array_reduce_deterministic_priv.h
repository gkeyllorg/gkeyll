#pragma once

// Fixed-order GPU sum: reduce ordered input cells to block partials, then
// combine those partials in a second kernel. USE_RANGE selects the index map.

template <bool USE_RANGE, unsigned int BLOCKSIZE>
__global__ static void
array_sum_partial_cu_ker(
  const struct gkyl_array *inp, struct gkyl_range range, long count, size_t comp,
  double *partial_dev
)
{
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;

  // Visit input cells in a fixed stride; range mode maps each cell to its array index.
  double local_sum = 0.0;
  for (long linc = (long)blockIdx.x * BLOCKSIZE + threadIdx.x; linc < count;
       linc += (long)gridDim.x * BLOCKSIZE) {
    long index = linc;
    if (USE_RANGE) {
      int idx[GKYL_MAX_DIM];
      gkyl_sub_range_inv_idx(&range, linc, idx);
      index = gkyl_range_idx(&range, idx);
    }
    local_sum += ((const double *)inp->data)[index * inp->ncomp + comp];
  }
  // Store one ordered partial sum per block.
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    partial_dev[blockIdx.x] = sum;
  }
}

template <unsigned int BLOCKSIZE>
__global__ static void
array_sum_finish_cu_ker(const double *partial_dev, int num_blocks, double *out_dev)
{
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;

  // Combine block partials in a fixed stride before the final block reduction.
  double local_sum = 0.0;
  for (int block = threadIdx.x; block < num_blocks; block += BLOCKSIZE) {
    local_sum += partial_dev[block];
  }
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    *out_dev = sum;
  }
}

template <bool USE_RANGE>
static void
array_sum_deterministic_cu(
  double *out_dev, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  // An empty input still writes a zero result for every component.
  long count = USE_RANGE ? range->volume : inp->size;
  if (count == 0) {
    const int nthreads = GKYL_DEFAULT_NUM_THREADS;
    array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
      out_dev, inp->ncomp, 0.0
    );
    cudaDeviceSynchronize();
    return;
  }

  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  const int max_blocks = 65535;
  long required_blocks = count / nthreads + (count % nthreads != 0);
  int num_blocks = (int)(required_blocks < max_blocks ? required_blocks : max_blocks);
  double *partial_dev = (double *)gkyl_cu_malloc((size_t)num_blocks * sizeof(double));
  struct gkyl_range input_range = {};
  if (USE_RANGE) {
    input_range = *range;
  }
  // Reuse the partial buffer while reducing components one at a time.
  for (size_t comp = 0; comp < inp->ncomp; ++comp) {
    array_sum_partial_cu_ker<USE_RANGE, nthreads>
      <<<num_blocks, nthreads>>>(inp->on_dev, input_range, count, comp, partial_dev);
    array_sum_finish_cu_ker<nthreads><<<1, nthreads>>>(partial_dev, num_blocks, out_dev + comp);
  }
  // The output may be mapped host memory, and the partials must outlive the launches.
  cudaDeviceSynchronize();
  gkyl_cu_free(partial_dev);
}
