#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

/**
 * Compute df = d^order/dx_dir^order f, where f and df are DG fields
 * expanded in the given basis, using a local derivative (simply
 * differentiate the DG expansion in every cell).
 *
 * c_oop and c_iop select which component of a vector DG field to use.
 *
 * @param basis Basis functions used in the expansion.
 * @param dir   Direction in which to differentiate (0-based).
 * @param order Order of the derivative (1,2).
 * @param dx    Cell length in direction of the derivative.
 * @param c_oop Component of the output field to write into.
 * @param out   Output DG field.
 * @param c_iop Component of the input field to read from.
 * @param inp   Input DG field.
 */
void gkyl_dg_differentiate_op_local(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp);

/**
 * Same as gkyl_dg_differentiate_op_local, but applied only over the
 * specified range (a sub-range of the full index range of the arrays).
 *
 * @param basis Basis functions used in the expansion.
 * @param dir   Direction in which to differentiate (0-based).
 * @param order Order of the derivative (1,2).
 * @param dx    Cell length in direction of the derivative.
 * @param c_oop Component of the output field to write into.
 * @param out   Output DG field.
 * @param c_iop Component of the input field to read from.
 * @param inp   Input DG field.
 * @param range Range over which to apply the operator.
 */
void gkyl_dg_differentiate_op_local_range(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp,
  const struct gkyl_range *range);
