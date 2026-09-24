#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_euler_pkpm_auxfields { 
  const struct gkyl_array *vlasov_pkpm_moms;
  const struct gkyl_array *pkpm_prim;
  const struct gkyl_array *pkpm_prim_surf;
  const struct gkyl_array *pkpm_p_ij;
  const struct gkyl_array *pkpm_lax;
  const struct gkyl_array *pkpm_penalization;
};

/**
 * Create a new Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing auxiliary variables
 * @return Pointer to Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model
 */
struct gkyl_dg_eqn* gkyl_dg_euler_pkpm_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* conf_range, bool use_gpu);

/**
 * Create new Euler equation object arallel-kinetic-perpendicular-moment (pkpm) model the lives on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_euler_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating Euler equation for parallel-kinetic-perpendicular-moment (pkpm) model.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_pkpm_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating Euler equation for parallel-kinetic-perpendicular-moment (pkpm) model.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_pkpm_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin);

#endif
