#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object definition
typedef struct gkyl_dg_array_mask gkyl_dg_array_mask;

enum gkyl_dg_array_mask_types {
  GKYL_DG_ARRAY_MASK_NONE = 0,             // No mask applied.
  GKYL_DG_ARRAY_MASK_C0_LESS,              // Mask applied based on the 0th component of the array.
  GKYL_DG_ARRAY_MASK_C0_LESS_FRAC,         // Applies the mask based on the global maximum 0th component of the input array.
  GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF,    // Applies the mask based on the local maximum 0th component of the input array, which is spatially dependent.
  GKYL_DG_ARRAY_MASK_C0_GREATER,           // Mask applied based on the 0th component of the array.
  GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC,      // Applies the mask based on the global maximum 0th component of the input array.
  GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF, // Applies the mask based on the local maximum 0th component of the input array, which is spatially dependent.
};

// Input structure for creating a mask object.
struct gkyl_dg_array_mask_inp {
  enum gkyl_dg_array_mask_types type;     // Type of mask.
  bool default_value;                     // Default value for mask (true/false) if no masking is applied.
  double threshold;                       // Threshold for marking cells as masked. Absolute value for *_THRESHOLD types, fraction (0-1) for *_FRAC_THRESHOLD types.
  const struct gkyl_range *phase_rng;     // Phase-space range.
  const struct gkyl_range *phase_rng_ext; // Extended phase-space range.
  const struct gkyl_range *conf_rng;      // Configuration-space range.
  const struct gkyl_range *conf_rng_ext;  // Extended configuration-space range.
  const struct gkyl_range *vel_rng;       // Velocity-space range.
  bool use_gpu;                           // Flag indicating GPU usage.
};

/**
 * Create a new mask object.
 * A mask is an object that applies a conditional to an input array.
 * The conditional is specified by the type parameter in the input structure.
 * The mask itself is 1 for true and -1 for false.
 *
 * @param mask_inp Input parameters for mask.
 * @return New mask object.
 */
struct gkyl_dg_array_mask*
gkyl_dg_array_mask_new(struct gkyl_dg_array_mask_inp mask_inp);

/**
 * Compute the appropriate threshold for fractional threshold masks.
 *
 * @param mask Mask updater.
 * @param global_max Global maximum value from input array.
 */
void
gkyl_dg_array_mask_advance_threshold(struct gkyl_dg_array_mask *mask,
  const double global_max);

/**
 * Compute the appropriate mask given the input array.
 *
 * @param mask Mask updater.
 * @param arr_in Array which the mask is applied on top of.
 */
void
gkyl_dg_array_mask_advance(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_in);

/**
 * Evaluate if the conditional mask is true at a given multi-dimensional index.
 *
 * Note that this method is intended for use outside of a loop over the grid or
 * outside a device kernel. If you need to do this inside a grid loop/device
 * kernel, you must use gkyl_dg_array_mask_eval_idx_ker defined in the private
 * header.
 *
 * @param mask Mask object.
 * @param idx Multi-dimensional index array.
 * @param val Value of the mask at the given index.
 */
void gkyl_dg_array_mask_eval_idx(struct gkyl_dg_array_mask *mask, const int *idx, bool *val);

/**
 * Scale the dg_array_mask by an array, modifying the mask.
 * mask = mask * arr_in
 *
 * @param mask Mask object.
 * @param arr_in Array to be multiplied by the mask.
 */
void
gkyl_dg_array_mask_scale_by_cell(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_in);

/**
 * Acquire a reference to the mask object.
 */
struct gkyl_dg_array_mask*
gkyl_dg_array_mask_acquire(struct gkyl_dg_array_mask *mask);

/**
 * Get the underlying mask array.
 *
 * @param mask Mask object.
 * @return Pointer to the mask array.
 */
const struct gkyl_array*
gkyl_dg_array_mask_get_mask(const struct gkyl_dg_array_mask *mask);

/**
 * Get device pointer to mask object.
 *
 * @param mask Mask object.
 */
struct gkyl_dg_array_mask*
gkyl_dg_array_mask_get_dev_ptr(struct gkyl_dg_array_mask *mask);

/**
 * Release memory associated with mask object.
 *
 * @param mask Mask object to release.
 */
void
gkyl_dg_array_mask_release(struct gkyl_dg_array_mask *mask);
