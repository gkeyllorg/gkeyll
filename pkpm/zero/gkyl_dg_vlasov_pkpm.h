#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_pkpm_auxfields { 
  const struct gkyl_array *bvar;
  const struct gkyl_array *bvar_surf;
  const struct gkyl_array *pkpm_prim;
  const struct gkyl_array *pkpm_prim_surf;
  const struct gkyl_array *max_b;
  const struct gkyl_array *pkpm_lax;
  const struct gkyl_array *div_b;
  const struct gkyl_array *pkpm_accel_vars;
  const struct gkyl_array *g_dist_source;
};

/**
 * Create a new Vlasov equation object for parallel-kinetic-perpendicular-moment (pkpm) model.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing configuration space fields
 * @param phase_range Phase space range for use in indexing phase space fields
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov equation object for parallel-kinetic-perpendicular-moment (pkpm) model.
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_pkpm_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* phase_range, bool use_gpu);

struct gkyl_dg_eqn* gkyl_dg_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* phase_range);

/**
 * Set the auxiliary fields 
 * (e.g. uvar, the bulk flow, pvar, the pressure tensor, and bvar, the magnetic field unit vector b and tensor bb)
 * needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_pkpm_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_pkpm_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields 
 * (e.g. uvar, the bulk flow, pvar, the pressure tensor, and bvar, the magnetic field unit vector b and tensor bb)
 * needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_pkpm_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_pkpm_auxfields auxin);


#endif
