#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_dg_bin_ops.h>

// Object type
typedef struct gkyl_deflated_dg_bin_ops gkyl_deflated_dg_bin_ops;

/**
 * Create new updater to divide or multiply 2 DG fields
 * along surfaces that are constant in the last coordinate.
 * The surfaces are lines (1d) if the initial fields are 2d
 * and planes (2d) if the initial fields are 3d).
 *
 * The advance method deflate sthe input fields onto each
 * surface, multiply or divide them at each surface,
 * and then inflates the result to provide a field with
 * the same dimension as the input fields.
 *
 * Free using gkyl_deflated_dg_bin_ops_release method.
 *
 *
 * @param grid Grid object
 * @param basis_on_dev Basis functions of the DG field stored on the gpu.
 * @param basis Basis functions of the DG field.
 * @param local range on which the operation should take place
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_deflated_dg_bin_ops* gkyl_deflated_dg_bin_ops_new(struct gkyl_rect_grid grid, 
  struct gkyl_basis *basis_on_dev, struct gkyl_basis basis, struct gkyl_range local, bool use_gpu);

/**
 * Multiply the two input fields on surfaces constant in the last dimension
 * Compute out = lop*rop. The c_oop, c_lop and c_rop are the
 * components into the DG fields to multiply (in case the field is a
 * vector field). For scalar fields c_oop = c_rop = c_lop = 0, for
 * example.
 *
 * @param up deflated_dg_bin_ops updater
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 *
 */
void gkyl_deflated_dg_bin_ops_mul(struct gkyl_deflated_dg_bin_ops* up, int c_oop, 
  struct gkyl_array *out, int c_lop, struct gkyl_array *lop, int c_rop, struct gkyl_array* rop);

/**
* Divide the two input fields on surfaces constant in the last dimension
* Compute out = lop/rop. The c_oop, c_lop and c_rop are the
* components into the DG fields to multiply (in case the field is a
* vector field). For scalar fields c_oop = c_rop = c_lop = 0, for
* example.
*
* @param up deflated_dg_bin_ops updater
* @param c_oop Component of output field in which to store product
* @param out Output DG field
* @param c_lop Component of left operand to use in product
* @param lop Left operand DG field
* @param c_rop Component of right operand to use in product
* @param rop Right operand DG field
*
*/
void gkyl_deflated_dg_bin_ops_div(struct gkyl_deflated_dg_bin_ops* up, int c_oop, 
  struct gkyl_array *out, int c_lop, struct gkyl_array *lop, int c_rop, struct gkyl_array* rop);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflated_dg_bin_ops_release(struct gkyl_deflated_dg_bin_ops* up);
