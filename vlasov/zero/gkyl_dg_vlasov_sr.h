#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_sr_auxfields { 
  const struct gkyl_array *qmem; // q/m * EM
  const struct gkyl_array *gamma; // gamma = sqrt(1 + p^2), particle Lorentz factor
};

/**
 * Create a new special relativistic Vlasov equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param vel_range Velocity space range for use in indexing p/gamma (velocity)
 * @param field_id enum to determine what type of EM fields
 * (special relativistic Vlasov-Maxwell vs. special relativistic neutrals)
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to special relativistic Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  enum gkyl_field_id field_id, bool use_gpu);

/**
 * Create a new special relativistic Vlasov equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param vel_range Velocity space range for use in indexing p/gamma (velocity)
 * @param field_id enum to determine what type of EM fields 
 * (special relativistic Vlasov-Maxwell vs. special relativistic neutrals)
 * @return Pointer to special relativistic Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  enum gkyl_field_id field_id);

/**
 * Set the auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_sr_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_sr_auxfields auxin);


#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_sr_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_sr_auxfields auxin);

#endif
