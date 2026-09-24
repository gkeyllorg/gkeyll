#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_calc_priv.h>
#include <gkyl_prim_lbo_kernels.h>
#include <gkyl_mat.h>
#include <assert.h>
#include <gkyl_prim_lbo_gyrokinetic.h>

// "derived" class constructors
struct gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_rng, bool use_gpu)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_gyrokinetic_new(cbasis, pbasis, use_gpu);
  struct gkyl_prim_lbo_calc *calc = gkyl_prim_lbo_calc_new(grid, prim, use_gpu);
  // Since calc now has pointer to specific type, decrease reference counter of type
  // so that eventual gkyl_prim_lbo_calc_release method on calculator deallocates specific type data
  gkyl_prim_lbo_type_release(prim);
  return calc;
}