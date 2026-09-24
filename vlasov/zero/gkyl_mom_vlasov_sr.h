#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_vlasov_sr_auxfields { 
  const struct gkyl_array *gamma; // gamma = sqrt(1 + p^2)
};

/**
 * Create new special relativistic Vlasov moment type object. 
 * Valid mom_types are M0, M1, M2, M3
 * Note: M2 is the integral(gamma*f) velocity moment and M3 the integral(p*f) velocity moment.
 * Ni = (M0, M1i) (1+vdim components).
 * Tij = stress-energy tensor (M2, M3i (vdim components), Stress tensor (vdim*(vdim+1))/2 components))
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param mom_type Name of moment to compute.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  const struct gkyl_range* vel_range, enum gkyl_distribution_moments mom_type, bool use_gpu);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object. Integrates (M0, M2, M3i) (vdim + 2) components.
 * Note the different order from non-relativistic since we compute the
 * integrated energy flux instead of the integrated mass flux. 
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param mom_type Name of moment to compute.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  const struct gkyl_range* vel_range, enum gkyl_distribution_moments mom_type, bool use_gpu);

/**
 * Set the auxiliary fields needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);
