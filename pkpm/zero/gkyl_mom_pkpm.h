#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

/**
 * Create new Vlasov parallel-kinetic-perpendicular-moment (pkpm) moment type object.
 * Can either compute moments needed for update (rho, p_par, p_perp) or diagnostic moments
 * (rho, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_range Range for indexing magnetic field unit vector and unit tensor
 * @param mass Mass of species (pkpm moments are mass weighted)
 * @param diag bool to determine if computing diagnostic moments 
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_pkpm_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, double mass, bool diag, bool use_gpu);

/**
 * Create new Vlasov moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, double mass, bool diag);
