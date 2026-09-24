#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_pkpm_diff_auxfields { 
  const struct gkyl_array *nuSum;
  const struct gkyl_array *nuPrimMomsSum;
};

/**
 * Create a new LBO diffusion term equation object for Vlasov equation in 
 * parallel-kinetic-perpendicular-moment (pkpm) model.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @param pgrid Phase-space grid object (used for robustness checking of vth^2 compared to vmax^2).
 * @param use_gpu Bool to determine if equation object is on host or device
 * @return Pointer to LBO diffusion term equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_pkpm_diff_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  const struct gkyl_rect_grid *pgrid, bool use_gpu);

struct gkyl_dg_eqn* gkyl_dg_lbo_pkpm_diff_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  const struct gkyl_rect_grid *pgrid);

/**
 * Set auxiliary fields needed in updating the diffusion flux term, 
 * nu*vt^2 (collision frequency * thermal velocity squared = nu*T/m).
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_pkpm_diff_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_pkpm_diff_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating the diffusion flux term, 
 * nu*vt^2 (collision frequency * thermal velocity squared = nu*T/m).
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_pkpm_diff_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_pkpm_diff_auxfields auxin);

#endif
