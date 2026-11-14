#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_diffusion_fluid_auxfields {
  const struct gkyl_array *D;
};

/**
 * Create a new fluid diffusion equation object.
 *
 * @param basis Basis functions of the equation system.
 * @param num_equations Number of scalar fluid equations.
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param use_gpu Whether to run on host or device.
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn *gkyl_dg_diffusion_fluid_new(
  const struct gkyl_basis *basis, int num_equations, const bool *diff_in_dir, int diff_order,
  bool use_gpu
);

/**
 * Set the auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_dg_diffusion_fluid_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_diffusion_fluid_auxfields auxin
);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_dg_diffusion_fluid_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_diffusion_fluid_auxfields auxin
);

#endif
