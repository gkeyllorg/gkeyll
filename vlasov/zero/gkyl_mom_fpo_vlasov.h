#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_fpo_vlasov_auxfields { 
  const struct gkyl_array *a; // (vector) drag coefficient
  const struct gkyl_array *D; // (tensor) diffusion coefficient
};

/**
 * Create new moment type for computing volume momentum and energy corrections for
 * the Fokker-Planck collision operator. Computes a 4-component moment with components
 * (a + div(D)), a . v + div(D . v) 
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param phase_range Range for indexing drag (a) and diffusion (D) coefficients
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_fpo_vlasov_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, bool use_gpu);

/**
 * Create new Vlasov moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_fpo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range);

/**
 * Set the auxiliary fields needed in computing momentum and energy correction moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_fpo_vlasov_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_fpo_vlasov_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields needed in computing 
 * momentum and energy correction moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_fpo_vlasov_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_fpo_vlasov_auxfields auxin);
#endif

