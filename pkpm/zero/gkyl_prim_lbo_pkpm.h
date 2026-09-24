#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>

/**
 * Create a new PKPM primitive moment object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration-space range
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov (with fluid coupling) primitive moment object
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_pkpm_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  bool use_gpu);

/**
 * Create a new Vlasov PKPM primitive moment object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_pkpm_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);
