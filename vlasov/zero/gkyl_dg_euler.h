#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_euler_auxfields { 
  const struct gkyl_array *u;
  const struct gkyl_array *p;
  const struct gkyl_array *u_surf;
  const struct gkyl_array *p_surf;
};

/**
 * Create a new Euler equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing auxiliary variables
 * @param wv_eqn Wave equation object which contains information and functions for the specific fluid equation
 * @param geom Wave geometry object for computing fluctuations local to surfaces
 * @param use_gpu Boolean to determine whether equation object is on host or device
 * @return Pointer to euler equation object
 */
struct gkyl_dg_eqn* gkyl_dg_euler_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range,
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, bool use_gpu);

/**
 * Create new Euler equation object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_euler_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range,
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom);

/**
 * Set the auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating euler equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating euler equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_auxfields auxin);

#endif
