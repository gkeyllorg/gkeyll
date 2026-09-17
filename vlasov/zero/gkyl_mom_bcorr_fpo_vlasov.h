#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_bcorr_fpo_vlasov_auxfields { 
  const struct gkyl_array *D; // (tensor) diffusion coefficient
};

/**
 * Create new fpo Vlasov boundary correction moment type object.
 * Fokker-Planck operator requires two sets of corrects for momentum and energy.
 * First set of corrections is identical to the LBO corrections.
 * Other set of corrections involves the integration of the diffusion tensor 
 * at the edge of velocity space.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param phase_range Phase-space range for indexing diffusion tensor
 * @param vBoundary Values at the edges of velocity space.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_mom_bcorr_fpo_vlasov_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, 
  const double* vBoundary, bool use_gpu);

/**
 * Create new fpo Vlasov boundary correction moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* 
gkyl_mom_bcorr_fpo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, 
  const double* vBoundary);

/**
 * Set the auxiliary fields needed in computing momentum and energy boundary correction moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_bcorr_fpo_vlasov_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_bcorr_fpo_vlasov_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields needed in computing 
 * momentum and energy boundary correction moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_bcorr_fpo_vlasov_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_bcorr_fpo_vlasov_auxfields auxin);

#endif
