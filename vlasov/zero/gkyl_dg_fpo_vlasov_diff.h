#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_fpo_vlasov_diff_auxfields {
  const struct gkyl_array *diff_coeff;
  const struct gkyl_array *diff_coeff_surf;
};

/**
 * Create a new FPO diffusion equation object.
 *
 * @param pbasis Phase-space basis functions
 * @param phase_range phase space range for use in indexing diffusion flux term
 * @return Pointer to fpo equation object
 */
struct gkyl_dg_eqn* gkyl_dg_fpo_vlasov_diff_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, bool use_gpu);

/**
 * Create a new FPO diffusion equation object that lives on NV-GPU
 *
 * @param pbasis Phase-space basis functions
 * @param phase_range phase space range for use in indexing diffusion flux term
 * @return Pointer to fpo equation object
 */
struct gkyl_dg_eqn* gkyl_dg_fpo_vlasov_diff_cu_dev_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range);

/**
 * Set auxiliary fields needed in updating the diffusion flux term (D = grad(grad(g)), g solved for externally).
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_fpo_vlasov_diff_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_fpo_vlasov_diff_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating diff flux term (D = grad(grad(g)), g solved for externally).
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_fpo_vlasov_diff_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_fpo_vlasov_diff_auxfields auxin);

#endif
