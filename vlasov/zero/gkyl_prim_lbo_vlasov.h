#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>

/**
 * Create a new Vlasov primitive moment object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov primitive moment object
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, bool use_gpu);

/**
 * Create a new Vlasov primitive type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);
