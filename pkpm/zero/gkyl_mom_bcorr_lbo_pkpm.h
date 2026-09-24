#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new lbo Vlasov boundary correction moment type object
 * for the parallel-kinetic-perpendicular-moment (pkpm) model 
 * Note that pkpm model is in local rest frame so corrections only apply to energy.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vBoundary Values at the edges of velocity space.
 * @param mass Mass of species
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_mom_bcorr_lbo_pkpm_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const double* vBoundary, 
  double mass, bool use_gpu);

/**
 * Create new LBO Vlasov boundary correction moment type object
 * for the parallel-kinetic-perpendicular-moment (pkpm) model on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* 
gkyl_mom_bcorr_lbo_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const double* vBoundary, 
  double mass);
