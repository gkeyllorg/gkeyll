#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_advection_auxfields { 
  const struct gkyl_array *u_i;
};

/**
 * Create a new advection equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing advection velocity
 * @return Pointer to advection equation object
 */
struct gkyl_dg_eqn* gkyl_dg_advection_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu);

/**
 * Create a new advection equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing advection velocity
 * @return Pointer to advection equation object
 */
struct gkyl_dg_eqn* gkyl_dg_advection_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. advection velocity u) needed in updating advection equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_advection_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields (e.g. advection velocity u) needed in updating advection equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_advection_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin);

#endif
