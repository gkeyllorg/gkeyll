#pragma once

#include <gkyl_tensor_field.h>


/**
 * Lowers the index of tensor using the associated metric in place (in ten)
 * 
 * @param metric Metric written as a tensor field
 * @param lowered_idx Index to lower of ten
 * @param ten Tensor field
 */
void gkyl_tensor_field_lower_idx_in_place(struct gkyl_tensor_field *metric, int lowered_idx, 
  struct gkyl_tensor_field *ten, struct gkyl_tensor_field *mem);

/**
 * Raises the index of tensor using the associated metric in place (in ten)
 * 
 * @param metric Metric written as a tensor field
 * @param raised_idx Index to raise of ten
 * @param ten Tensor field
 */
void gkyl_tensor_field_raise_idx_in_place(struct gkyl_tensor_field *metric, int raised_idx,
  struct gkyl_tensor_field *ten, struct gkyl_tensor_field *mem);

/**
 * Lowers the index of tensor using the associated metric in place (in ten)
 * 
 * @param metric Metric written as a tensor field
 * @param lowered_idx Index to lower of ten
 * @param ten Tensor field
 * @param ten_out (output) Tensor field
 */
void gkyl_tensor_field_lower_idx_set(const struct gkyl_tensor_field *metric, int lowered_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *ten_out);

/**
 * Raises the index of tensor using the associated metric in place (in ten)
 * 
 * @param metric Metric written as a tensor field
 * @param raised_idx Index to raise of ten
 * @param ten Tensor field
 * @param ten_out (output) Tensor field
 */
void gkyl_tensor_field_raise_idx_set(const struct gkyl_tensor_field *metric, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *ten_out);

/**
 * Host-side wrappers for array operations
 */
void
tensor_field_raise_or_lower_idx_set_cu(const struct gkyl_tensor_field *met, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *tensor_out);