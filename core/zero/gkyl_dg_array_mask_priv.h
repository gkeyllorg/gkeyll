#pragma once

// Private header for dg_array_mask object, not for direct use in user code.

#include <gkyl_dg_array_mask.h>
#include <gkyl_ref_count.h>

/**
 * Skip cell object definition.
 */
struct gkyl_dg_array_mask {
  enum gkyl_dg_array_mask_types type;     // Type of mask operation.
  bool default_value;                     // Default value for mask (true/false) if no masking is applied. Defaults to false (-1.0).
  struct gkyl_array *mask_arr;            // Mask array (1.0 is true, -1.0 is false).
  double threshold;                       // Threshold for marking cells as masked. Scaled absolute value for *_THRESHOLD types, fraction for *_FRAC_THRESHOLD types.
  double frac_threshold;                  // Fractional threshold, which is 0-1 for global fraction mask.
  const struct gkyl_range *mask_rng;      // Pointer to range over which mask is applied (phase_rng for kinetic, conf_rng for fluid).
  int mask_rng_ndim;                      // Dimensionality of the mask's range.
  const struct gkyl_range *mask_rng_ext;  // Pointer to extended range for mask allocation.
  const struct gkyl_range *conf_rng;      // Configuration-space range.
  const struct gkyl_range *conf_rng_ext;  // Extended configuration-space range.
  const struct gkyl_range *vel_rng;       // Velocity-space range.
  bool use_gpu;                           // Flag indicating GPU usage.

  struct gkyl_array *local_max_arr; // Pre-allocated config-space array for spatial fractional masks.

  // Function pointer for advance method (CPU), set at init time based on mask type.
  void (*advance_func)(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask);
  void (*advance_threshold_func)(struct gkyl_dg_array_mask *mask, const double global_max);

  // Function pointer for advance method (GPU), set at init time based on mask type.
  void (*advance_func_cu)(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask);

  // Function that evaluates the mask.
  bool (*eval_idx_func)(struct gkyl_dg_array_mask *mask, const int *idx);

  // Function pointer for scale_by_cell method, set at init time based on mask type.
  void (*scale_by_cell_func)(struct gkyl_dg_array_mask *mask,
    const struct gkyl_array *arr_to_multiply);

  uint32_t flags;
  struct gkyl_dg_array_mask *on_dev; // Pointer to device object.

  struct gkyl_ref_count ref_count; // Reference counter.
};

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref);

GKYL_CU_DH
static bool
eval_idx_ker_disabled(struct gkyl_dg_array_mask *mask, const int *idx)
{
  return mask->default_value;
}

GKYL_CU_DH
static bool
eval_idx_ker_enabled(struct gkyl_dg_array_mask *mask, const int *idx)
{
  long linidx = gkyl_range_idx(mask->mask_rng, idx);
  const double *mask_c = (const double *)gkyl_array_cfetch(mask->mask_arr, linidx);
  return *mask_c > 0; // Returns true if the mask is true.
}

/**
 * Function to evaluate the mask at an index, to be used inside kernels.
 *
 * @param mask Mask object.
 * @param idx Multi-dimensional index array.
 * @return Value of the mask at the given index.
 */
GKYL_CU_DH
static inline bool
gkyl_dg_array_mask_eval_idx_ker(struct gkyl_dg_array_mask *mask, const int *idx)
{
  return mask->eval_idx_func(mask, idx);
}

#ifdef GKYL_HAVE_CUDA

/**
 * Create a new dg_array_mask object on CUDA device.
 *
 * @param mask_ho Host-side dg_array_mask object.
 * @return New dg_array_mask object on device.
 */
struct gkyl_dg_array_mask* gkyl_dg_array_mask_cu_dev_new(struct gkyl_dg_array_mask *mask_ho);

/**
 * CUDA device function to update dg_array_mask on GPU.
 *
 * @param mask dg_array_mask object.
 * @param arr_to_mask Array to mask.
 */
void gkyl_dg_array_mask_advance_cu(struct gkyl_dg_array_mask *mask,
  const struct gkyl_array *arr_to_mask);


/**
 * CUDA device function to update the mask's threshold on GPU.
 *
 * @param mask dg_array_mask object.
 * @param global_max Global maximum value used for fractional thresholding.
 */
void
gkyl_dg_array_mask_advance_threshold_cu(struct gkyl_dg_array_mask *mask,
  const double global_max);

/**
 * CUDA device function to evaluate the mask at an index.
 *
 * @param mask dg_array_mask object.
 * @param idx Index to evaluate the mask at.
 * @param val Value of the mask at given index.
 */
void gkyl_dg_array_mask_eval_idx_cu(struct gkyl_dg_array_mask *mask, const int *idx, bool *val);
#endif
