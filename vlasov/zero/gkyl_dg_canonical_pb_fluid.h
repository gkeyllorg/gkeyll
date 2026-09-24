#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_wv_eqn.h>

// Struct containing the pointers to auxiliary fields.
// Potential is in *Canonical* cordinates.
struct gkyl_dg_canonical_pb_fluid_auxfields { 
  const struct gkyl_array *phi;
  const struct gkyl_array *alpha_surf;
  const struct gkyl_array *sgn_alpha_surf;
  const struct gkyl_array *const_sgn_alpha;
};

/**
 * Create a new canonical-pb equation object for fluid systems object
 * Utilized for systems such as incompressible Euler, Hasegawa-Mima,
 * or Hasegawa-Wakatani that can be described by a canonical Poisson Bracket.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for indexing potential and surface alpha.
 * @param wv_eqn Wave equation object which contains information about specific fluid system.
 *               For example, how many equations are we solving (1 for Hasegawa-Mima, 2 for Hasegawa-Wakatani).
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to special canonical-pb equation object
 */
struct gkyl_dg_eqn* gkyl_dg_canonical_pb_fluid_new(const struct gkyl_basis* cbasis,
  const struct gkyl_range* conf_range, const struct gkyl_wv_eqn *wv_eqn, bool use_gpu);

/**
 * Create a new canonical-pb equation object for fluid systems object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_canonical_pb_fluid_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_range* conf_range, const struct gkyl_wv_eqn *wv_eqn);

/**
 * Set the auxiliary fields
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_canonical_pb_fluid_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_fluid_auxfields auxin);


#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_canonical_pb_fluid_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_fluid_auxfields auxin);

#endif
