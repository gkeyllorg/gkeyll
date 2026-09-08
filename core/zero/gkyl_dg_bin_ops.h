#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

// Type for storing preallocating memory needed in various batch
// operations
typedef struct gkyl_dg_bin_op_mem gkyl_dg_bin_op_mem;

/**
 * Allocate memory for use in bin op (division operator). Free using
 * release method.
 *
 * @param nbatch Batch size
 * @param neqn Number of equations in each batch 
 */
gkyl_dg_bin_op_mem *gkyl_dg_bin_op_mem_new(size_t nbatch, size_t neqn);
// Same as above, except for GPUs
gkyl_dg_bin_op_mem *gkyl_dg_bin_op_mem_cu_dev_new(size_t nbatch, size_t neqn);

/**
 * Release memory needed in the bin ops.
 *
 * @param mem Memory to release
 */
void gkyl_dg_bin_op_mem_release(gkyl_dg_bin_op_mem *mem);

/**
 * Compute out = lop*rop. The c_oop, c_lop and c_rop are the
 * components into the DG fields to multiply (in case the field is a
 * vector field). For scalar fields c_oop = c_rop = c_lop = 0, for
 * example.
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 */
void gkyl_dg_mul_op(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out, int c_lop,
  const struct gkyl_array *lop, int c_rop, const struct gkyl_array *rop);

/**
 * Same as gkyl_dg_mul_op, except operator is applied only on
 * specified range (sub-range of range containing the DG fields).
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_mul_op_range(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_lop, const struct gkyl_array *lop, int c_rop, const struct gkyl_array *rop,
  const struct gkyl_range *range);

/**
 * Compute pout = cop*pop on specified range (sub-range of range
 * containing the DG fields), where pout and pop are phase-space
 * operands, and cop is a conf-space operand.
 *
 * @param cbasis Configuration space basis functions used in expansions.
 * @param pbasis Phase space basis functions used in expansions.
 * @param pout Output phase-space DG field.
 * @param cop Conf-space operand DG field.
 * @param pop Phase-space operand DG field.
 * @param crange Conf-space range to apply multiplication operator.
 * @param prange Phase-space range to apply multiplication operator.
 */
void gkyl_dg_mul_conf_phase_op_range(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array *pout, const struct gkyl_array *cop,
  const struct gkyl_array *pop, const struct gkyl_range *crange, const struct gkyl_range *prange);

/**
 * Compute pout += a*cop*pop on specified range (sub-range of range
 * containing the DG fields), where pout and pop are phase-space
 * operands, and cop is a conf-space operand. 
 *
 * @param cbasis Configuration space basis functions used in expansions.
 * @param pbasis Phase space basis functions used in expansions.
 * @param pout Output phase-space DG field.
 * @param a Factor to multiply accumulated output by (e.g., a = -1.0 for decrementing)
 * @param cop Conf-space operand DG field.
 * @param pop Phase-space operand DG field.
 * @param crange Conf-space range to apply multiplication operator.
 * @param prange Phase-space range to apply multiplication operator.
 */
void gkyl_dg_mul_conf_phase_op_accumulate_range(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array *pout, double a, const struct gkyl_array *cop,
  const struct gkyl_array *pop, const struct gkyl_range *crange, const struct gkyl_range *prange);

/**
 * Compute out = lop . rop, where lop and rop are vector fields.
 * If these vector fields only have one component the result is the
 * same as using dg_mul_op with c_oop=c_lop=c_rop=0.
 *
 * @param basis Basis functions used in expansions
 * @param out Output DG scalar field.
 * @param lop Left operand DG vector field.
 * @param rop Right operand DG vector field.
 */
void gkyl_dg_dot_product_op(const struct gkyl_basis *basis, struct gkyl_array *out,
  const struct gkyl_array *lop, const struct gkyl_array *rop);

/**
 * Same as gkyl_dg_dot_product_op, except operator is applied only on
 * specified range (sub-range of range containing the DG fields).
 *
 * @param basis Basis functions used in expansions.
 * @param out Output DG scalar field.
 * @param lop Left operand DG vector field.
 * @param rop Right operand DG vector field.
 * @param range Range to apply dot product operator.
 */
void gkyl_dg_dot_product_op_range(const struct gkyl_basis *basis, struct gkyl_array *out,
  const struct gkyl_array *lop, const struct gkyl_array *rop, const struct gkyl_range *range);

/**
 * Compute out = lop/rop. The c_oop, c_lop and c_rop are the
 * components into the DG fields to divide (in case the field is a
 * vector field). For scalar fields c_oop = c_rop = c_lop = 0, for
 * example.
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 */
void gkyl_dg_div_op(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis, int c_oop,
  struct gkyl_array *out, int c_lop, const struct gkyl_array *lop, int c_rop,
  const struct gkyl_array *rop);

/**
 * Same as gkyl_dg_div_op, except operator is applied only on
 * specified range (sub-range of range containing the DG fields).
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_div_op_range(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis, int c_oop,
  struct gkyl_array *out, int c_lop, const struct gkyl_array *lop, int c_rop,
  const struct gkyl_array *rop, const struct gkyl_range *range);

/**
 * Compute out = 1/iop. The c_oop and c_iop are the
 * components into the DG fields to use (in case the fields are
 * vector fields). For scalar fields c_oop = c_iop = 0, for
 * example.
 *
 * @param basis Basis functions used in expansions.
 * @param c_oop Component of output field in which to store product.
 * @param out Output DG field.
 * @param c_iop Component of input operand.
 * @param iop Input operand DG field.
 */
void gkyl_dg_inv_op(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out, int c_iop,
  const struct gkyl_array *iop);

/**
 * Compute out = 1/iop on specified range. The c_oop and c_iop are
 * the components into the DG fields to use (in case the fields are
 * vector fields). For scalar fields c_oop = c_iop = 0, for
 * example.
 *
 * @param basis Basis functions used in expansions.
 * @param c_oop Component of output field in which to store product.
 * @param out Output DG field.
 * @param c_iop Component of input operand.
 * @param iop Input operand DG field.
 */
void gkyl_dg_inv_op_range(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop, const struct gkyl_range *range);

/**
 * Compute the cell-average of input array iop and store it in out
 * array.
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field 
 * @param out Output DG field
 * @param c_iop Component of input DG field
 * @param iop Input DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_calc_average_range(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop, struct gkyl_range range);

/**
 * Compute the mean L2 norm of input array iop and store it in out
 * array.
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field 
 * @param out Output DG field
 * @param c_iop Component of input DG field
 * @param iop Input DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_calc_l2_range(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop, struct gkyl_range range);

/**
 * Host-side wrappers for dg_bin_op operations
 */
void gkyl_dg_mul_op_cu(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out, int c_lop,
  const struct gkyl_array *lop, int c_rop, const struct gkyl_array *rop);

void gkyl_dg_mul_op_range_cu(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_lop, const struct gkyl_array *lop, int c_rop, const struct gkyl_array *rop,
  const struct gkyl_range *range);

void gkyl_dg_mul_conf_phase_op_range_cu(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array *pout, const struct gkyl_array *cop,
  const struct gkyl_array *pop, const struct gkyl_range *crange, const struct gkyl_range *prange);

void gkyl_dg_mul_conf_phase_op_accumulate_range_cu(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array *pout, double a, const struct gkyl_array *cop,
  const struct gkyl_array *pop, const struct gkyl_range *crange, const struct gkyl_range *prange);

void gkyl_dg_dot_product_op_cu(const struct gkyl_basis *basis, struct gkyl_array *out,
  const struct gkyl_array *lop, const struct gkyl_array *rop);

void gkyl_dg_dot_product_op_range_cu(const struct gkyl_basis *basis, struct gkyl_array *out,
  const struct gkyl_array *lop, const struct gkyl_array *rop, const struct gkyl_range *range);

void gkyl_dg_div_op_cu(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis, int c_oop,
  struct gkyl_array *out, int c_lop, const struct gkyl_array *lop, int c_rop,
  const struct gkyl_array *rop);

void gkyl_dg_div_op_range_cu(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis, int c_oop,
  struct gkyl_array *out, int c_lop, const struct gkyl_array *lop, int c_rop,
  const struct gkyl_array *rop, const struct gkyl_range *range);

void gkyl_dg_inv_op_cu(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out, int c_iop,
  const struct gkyl_array *iop);

void gkyl_dg_inv_op_range_cu(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop, const struct gkyl_range *range);
