/* -*- c++ -*- */
extern "C" {
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_reduce.h>
#include <assert.h>
#include <float.h>
}

// CUDA kernel for LESS_THAN_THRESHOLD masks
__global__ void
gkyl_dg_array_mask_less_than_kernel(struct gkyl_range mask_rng,
  const struct gkyl_array *arr_to_mask, struct gkyl_array *mask, double f_threshold)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < mask_rng.volume; tid += blockDim.x * gridDim.x) {

    gkyl_sub_range_inv_idx(&mask_rng, tid, idx);

    if (gkyl_range_contains_idx(&mask_rng, idx)) {
      long linidx = gkyl_range_idx(&mask_rng, idx);
      const double *distf_c = (const double *)gkyl_array_cfetch(arr_to_mask, linidx);
      double *mask_c = (double *)gkyl_array_fetch(mask, linidx);

      double abs_val = fabs(distf_c[0]);
      *mask_c = (abs_val < f_threshold) ? 1.0 : -1.0;
    }
  }
}

// CUDA kernel for GREATER_THAN_THRESHOLD masks
__global__ void
gkyl_dg_array_mask_greater_than_kernel(struct gkyl_range mask_rng,
  const struct gkyl_array *arr_to_mask, struct gkyl_array *mask, double f_threshold)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < mask_rng.volume; tid += blockDim.x * gridDim.x) {

    gkyl_sub_range_inv_idx(&mask_rng, tid, idx);

    if (gkyl_range_contains_idx(&mask_rng, idx)) {
      long linidx = gkyl_range_idx(&mask_rng, idx);
      const double *distf_c = (const double *)gkyl_array_cfetch(arr_to_mask, linidx);
      double *mask_c = (double *)gkyl_array_fetch(mask, linidx);

      double abs_val = fabs(distf_c[0]);
      *mask_c = (abs_val > f_threshold) ? 1.0 : -1.0;
    }
  }
}

// CUDA kernel to find max value in velocity space for each configuration cell
__global__ void
gkyl_dg_array_mask_find_local_max_kernel(struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  struct gkyl_range mask_rng, const struct gkyl_array *arr_to_mask,
  struct gkyl_array *local_max_arr)
{
  int conf_idx[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < conf_rng.volume; tid += blockDim.x * gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_rng, tid, conf_idx);

    if (gkyl_range_contains_idx(&conf_rng, conf_idx)) {
      double local_max = -DBL_MAX;

      // Iterate over velocity space for this configuration cell
      for (unsigned long vel_tid = 0; vel_tid < vel_rng.volume; vel_tid++) {
        int vel_idx[GKYL_MAX_DIM];
        gkyl_sub_range_inv_idx(&vel_rng, vel_tid, vel_idx);

        // Build phase space index
        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < conf_rng.ndim; d++) {
          pidx[d] = conf_idx[d];
        }
        for (int d = 0; d < vel_rng.ndim; d++) {
          pidx[conf_rng.ndim + d] = vel_idx[d];
        }

        long linidx_phase = gkyl_range_idx(&mask_rng, pidx);
        const double *arr_c = (const double *)gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double abs_val = fabs(arr_c[0]);

        if (abs_val > local_max) {
          local_max = abs_val;
        }
      }

      long conf_linidx = gkyl_range_idx(&conf_rng, conf_idx);
      double *local_max_c = (double *)gkyl_array_fetch(local_max_arr, conf_linidx);
      local_max_c[0] = local_max;
    }
  }
}

// CUDA kernel to apply spatial fractional mask (LESS_THAN)
__global__ void
gkyl_dg_array_mask_spatial_frac_less_than_kernel(struct gkyl_range conf_rng,
  struct gkyl_range vel_rng, struct gkyl_range mask_rng, const struct gkyl_array *arr_to_mask,
  struct gkyl_array *mask, const struct gkyl_array *local_max_arr, double frac_threshold)
{
  int conf_idx[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < conf_rng.volume; tid += blockDim.x * gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_rng, tid, conf_idx);

    if (gkyl_range_contains_idx(&conf_rng, conf_idx)) {
      long conf_linidx = gkyl_range_idx(&conf_rng, conf_idx);
      const double *local_max_c = (const double *)gkyl_array_cfetch(local_max_arr, conf_linidx);
      double local_threshold = frac_threshold * local_max_c[0];

      // Apply mask to all velocity cells in this configuration cell
      for (unsigned long vel_tid = 0; vel_tid < vel_rng.volume; vel_tid++) {
        int vel_idx[GKYL_MAX_DIM];
        gkyl_sub_range_inv_idx(&vel_rng, vel_tid, vel_idx);

        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < conf_rng.ndim; d++) {
          pidx[d] = conf_idx[d];
        }
        for (int d = 0; d < vel_rng.ndim; d++) {
          pidx[conf_rng.ndim + d] = vel_idx[d];
        }

        long linidx_phase = gkyl_range_idx(&mask_rng, pidx);
        const double *arr_c = (const double *)gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double *mask_c = (double *)gkyl_array_fetch(mask, linidx_phase);

        double abs_val = fabs(arr_c[0]);
        *mask_c = (abs_val < local_threshold) ? 1.0 : -1.0;
      }
    }
  }
}

// CUDA kernel to apply spatial fractional mask (GREATER_THAN)
__global__ void
gkyl_dg_array_mask_spatial_frac_greater_than_kernel(struct gkyl_range conf_rng,
  struct gkyl_range vel_rng, struct gkyl_range mask_rng, const struct gkyl_array *arr_to_mask,
  struct gkyl_array *mask, const struct gkyl_array *local_max_arr, double frac_threshold)
{
  int conf_idx[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < conf_rng.volume; tid += blockDim.x * gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_rng, tid, conf_idx);

    if (gkyl_range_contains_idx(&conf_rng, conf_idx)) {
      long conf_linidx = gkyl_range_idx(&conf_rng, conf_idx);
      const double *local_max_c = (const double *)gkyl_array_cfetch(local_max_arr, conf_linidx);
      double local_threshold = frac_threshold * local_max_c[0];

      // Apply mask to all velocity cells in this configuration cell
      for (unsigned long vel_tid = 0; vel_tid < vel_rng.volume; vel_tid++) {
        int vel_idx[GKYL_MAX_DIM];
        gkyl_sub_range_inv_idx(&vel_rng, vel_tid, vel_idx);

        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < conf_rng.ndim; d++) {
          pidx[d] = conf_idx[d];
        }
        for (int d = 0; d < vel_rng.ndim; d++) {
          pidx[conf_rng.ndim + d] = vel_idx[d];
        }

        long linidx_phase = gkyl_range_idx(&mask_rng, pidx);
        const double *arr_c = (const double *)gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double *mask_c = (double *)gkyl_array_fetch(mask, linidx_phase);

        double abs_val = fabs(arr_c[0]);
        *mask_c = (abs_val > local_threshold) ? 1.0 : -1.0;
      }
    }
  }
}

// Static host functions that launch the appropriate CUDA kernels.
// These are assigned to the advance_func_cu function pointer at init time.

static void
advance_cu_none(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask)
{
  // Do nothing for NONE type
}

static void
advance_cu_less_than(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask)
{
  int nblocks = arr_to_mask->nblocks;
  int nthreads = arr_to_mask->nthreads;
  gkyl_dg_array_mask_less_than_kernel<<<nblocks, nthreads>>>(
    *mask->mask_rng, arr_to_mask->on_dev, mask->mask_arr->on_dev, mask->threshold);
}

static void
advance_cu_greater_than(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask)
{
  int nblocks = arr_to_mask->nblocks;
  int nthreads = arr_to_mask->nthreads;
  gkyl_dg_array_mask_greater_than_kernel<<<nblocks, nthreads>>>(
    *mask->mask_rng, arr_to_mask->on_dev, mask->mask_arr->on_dev, mask->threshold);
}

static void
advance_cu_less_than_frac_conf(struct gkyl_dg_array_mask *mask,
  const struct gkyl_array *arr_to_mask)
{
  int nthreads = mask->mask_arr->nthreads;
  int conf_nblocks = mask->conf_rng->nblocks;

  // Phase 1: Find max in velocity space for each configuration cell
  gkyl_dg_array_mask_find_local_max_kernel<<<conf_nblocks, nthreads>>>(
    *mask->conf_rng, *mask->vel_rng, *mask->mask_rng, arr_to_mask->on_dev,
    mask->local_max_arr->on_dev);

  // Phase 2: Apply mask based on local thresholds
  gkyl_dg_array_mask_spatial_frac_less_than_kernel<<<conf_nblocks, nthreads>>>(
    *mask->conf_rng, *mask->vel_rng, *mask->mask_rng, arr_to_mask->on_dev,
    mask->mask_arr->on_dev, mask->local_max_arr->on_dev, mask->threshold);
}

static void
advance_cu_greater_than_frac_conf(struct gkyl_dg_array_mask *mask,
  const struct gkyl_array *arr_to_mask)
{
  int nthreads = mask->mask_arr->nthreads;
  int conf_nblocks = mask->conf_rng->nblocks;

  // Phase 1: Find max in velocity space for each configuration cell
  gkyl_dg_array_mask_find_local_max_kernel<<<conf_nblocks, nthreads>>>(
    *mask->conf_rng, *mask->vel_rng, *mask->mask_rng, arr_to_mask->on_dev,
    mask->local_max_arr->on_dev);

  // Phase 2: Apply mask based on local thresholds
  gkyl_dg_array_mask_spatial_frac_greater_than_kernel<<<conf_nblocks, nthreads>>>(
    *mask->conf_rng, *mask->vel_rng, *mask->mask_rng, arr_to_mask->on_dev,
    mask->mask_arr->on_dev, mask->local_max_arr->on_dev, mask->threshold);
}

// Host function to launch CUDA kernel.
void
gkyl_dg_array_mask_advance_cu(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask)
{
  mask->advance_func_cu(mask, arr_to_mask);
}

__global__ void
advance_threshold_none(struct gkyl_dg_array_mask *mask, const double global_max)
{
}

__global__ void
advance_threshold_frac_kernel(struct gkyl_dg_array_mask *mask, const double global_max)
{
  mask->threshold = mask->frac_threshold * global_max;
}

static void
advance_threshold_frac(struct gkyl_dg_array_mask *mask, const double global_max)
{
  advance_threshold_frac_kernel<<<1, 1>>>(mask->on_dev, global_max);
}

void
gkyl_dg_array_mask_advance_threshold_cu(struct gkyl_dg_array_mask *mask, const double global_max)
{
  mask->advance_threshold_func(mask, global_max);
}

struct dg_array_mask_idx {
  int idx[GKYL_MAX_DIM]; // Index.
};

__global__ static void
gkyl_dg_array_mask_eval_idx_kernel(struct gkyl_dg_array_mask *mask,
  struct dg_array_mask_idx idx_struct, bool *val)
{
  val[0] = mask->eval_idx_func(mask, idx_struct.idx);
}

void
gkyl_dg_array_mask_eval_idx_cu(struct gkyl_dg_array_mask *mask, const int *idx, bool *val)
{
  struct dg_array_mask_idx idx_struct = {};
  for (int d = 0; d < mask->mask_rng_ndim; d++) {
    idx_struct.idx[d] = idx[d];
  }

  gkyl_dg_array_mask_eval_idx_kernel<<<1, 1>>>(mask->on_dev, idx_struct, val);
}

__global__ static void
gkyl_dg_array_mask_set_dev_func_ptr(struct gkyl_dg_array_mask *mask, enum gkyl_dg_array_mask_types)
{
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE)
    mask->eval_idx_func = eval_idx_ker_disabled;
  else
    mask->eval_idx_func = eval_idx_ker_enabled;
}

struct gkyl_dg_array_mask*
gkyl_dg_array_mask_cu_dev_new(struct gkyl_dg_array_mask *mask_ho)
{
  struct gkyl_dg_array_mask *mask = (struct gkyl_dg_array_mask *)gkyl_malloc(sizeof(*mask));

  mask->type = mask_ho->type;
  mask->threshold = mask_ho->threshold;
  mask->frac_threshold = mask_ho->frac_threshold;
  mask->use_gpu = mask_ho->use_gpu;
  mask->conf_rng = mask_ho->conf_rng;
  mask->conf_rng_ext = mask_ho->conf_rng_ext;
  mask->vel_rng = mask_ho->vel_rng;
  mask->mask_rng = mask_ho->mask_rng;
  mask->mask_rng_ext = mask_ho->mask_rng_ext;
  mask->mask_rng_ndim = mask_ho->mask_rng_ndim;
  mask->mask_arr = 0;
  mask->local_max_arr = 0;

  // Copy host function pointers (used for CPU-side dispatch)
  mask->advance_func = mask_ho->advance_func;
  mask->advance_threshold_func = mask_ho->advance_threshold_func;
  mask->scale_by_cell_func = mask_ho->scale_by_cell_func;

  mask->flags = 0;
  GKYL_SET_CU_ALLOC(mask->flags);
  mask->ref_count = gkyl_ref_count_init(gkyl_dg_array_mask_free);

  mask->advance_func_cu = advance_cu_none;

  // Set GPU advance function pointer based on mask type
  switch (mask->type) {
    case GKYL_DG_ARRAY_MASK_NONE:
      mask->advance_func_cu = advance_cu_none;
      break;
    case GKYL_DG_ARRAY_MASK_C0_LESS:
      mask->advance_func_cu = advance_cu_less_than;
      break;
    case GKYL_DG_ARRAY_MASK_C0_GREATER:
      mask->advance_func_cu = advance_cu_greater_than;
      break;
    case GKYL_DG_ARRAY_MASK_C0_LESS_FRAC:
      mask->advance_func_cu = advance_cu_less_than;
      break;
    case GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC:
      mask->advance_func_cu = advance_cu_greater_than;
      break;
    case GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF:
      mask->advance_func_cu = advance_cu_less_than_frac_conf;
      break;
    case GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF:
      mask->advance_func_cu = advance_cu_greater_than_frac_conf;
      break;
    default:
      break;
  }

  if (mask->type == GKYL_DG_ARRAY_MASK_NONE)
    mask->eval_idx_func = eval_idx_ker_disabled;
  else
    mask->eval_idx_func = eval_idx_ker_enabled;

  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    // Initialize the device object.
    struct gkyl_dg_array_mask *mask_cu =
      (struct gkyl_dg_array_mask *)gkyl_cu_malloc(sizeof(*mask_cu));
    gkyl_cu_memcpy(mask_cu, mask, sizeof(struct gkyl_dg_array_mask), GKYL_CU_MEMCPY_H2D);
    mask->on_dev = mask_cu;
  }
  else {
    struct gkyl_array *mask_array = gkyl_array_cu_dev_new(GKYL_DOUBLE, mask_ho->mask_arr->ncomp,
      mask_ho->mask_arr->size);
    gkyl_array_copy(mask_array, mask_ho->mask_arr);
    mask->mask_arr = mask_array->on_dev;
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF) {
      mask->local_max_arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, mask->conf_rng_ext->volume);
    }
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC) {}

    // Initialize the device object.
    struct gkyl_dg_array_mask *mask_cu =
      (struct gkyl_dg_array_mask *)gkyl_cu_malloc(sizeof(*mask_cu));
    gkyl_cu_memcpy(mask_cu, mask, sizeof(struct gkyl_dg_array_mask), GKYL_CU_MEMCPY_H2D);
    mask->on_dev = mask_cu;

    // The returned object should store host pointer to gkyl_array.
    mask->mask_arr = mask_array;
  }

  gkyl_dg_array_mask_set_dev_func_ptr<<<1, 1>>>(mask->on_dev, mask->type);

  // For NONE type, don't allocate mask array
  return mask;
}
