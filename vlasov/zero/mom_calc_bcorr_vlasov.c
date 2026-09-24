#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_calc_bcorr_priv.h>
#include <gkyl_util.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>

// "derived" class constructors
struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, bool use_gpu)
{
  struct gkyl_mom_type *bcorr_type; // LBO boundary corrections moment type
  bcorr_type = gkyl_mom_bcorr_lbo_vlasov_new(cbasis, pbasis, vBoundary, use_gpu);  
  struct gkyl_mom_calc_bcorr* calc = gkyl_mom_calc_bcorr_new(grid, bcorr_type, use_gpu);
  // Since calc now has pointer to specific type, decrease reference counter of type
  // so that eventual gkyl_mom_calc_bcorr_release method on calculator deallocates specific type data
  gkyl_mom_type_release(bcorr_type);
  return calc;
}