#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Auxiliary arrays of the Euler equations as passed through the fluid updater (aux_inp of
// gkyl_dg_updater_fluid_inp).
struct gkyl_dg_euler_auxfields {
  const struct gkyl_array *u;
  const struct gkyl_array *p;
};

// Input packaged as a struct
struct gkyl_dg_euler_inp {
  const struct gkyl_basis *cbasis; // Configuration-space basis functions.
  const struct gkyl_range
    *conf_range; // Configuration-space range for indexing u, p and the position-map Jacobian.
  const struct gkyl_wv_eqn *
    wv_eqn; // Wave equation object (Euler or isothermal Euler): its Riemann solver is the interface flux.
  const struct gkyl_wave_geom *geom; // Wave geometry for the frame of each interface.
  // Configuration-space position map (REQUIRED). On a non-identity map the state is the physical
  // one and each cell has its physical size J_d dx_d (the C^0 linear map has a constant J_d per
  // cell), so the update is the rate of the physical state.
  const struct gkyl_vlasov_position_map *pos_map;
  const struct gkyl_array *u; // Flow velocity [ux, uy, uz] (REQUIRED).
  const struct gkyl_array *p; // Pressure (REQUIRED).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new Euler equation object: Euler, [rho, rho ux, rho uy, rho uz, energy], or isothermal
 * Euler, [rho, rho ux, rho uy, rho uz], by the type of the wave equation object. Its boundary
 * surface term is the absorbing wall: the flux of the particles of the local Maxwellian of the
 * boundary cell that move toward the wall (half-range moments; none come back) on the wall face,
 * and the Riemann flux on the interior face. The updater applies it on the edges it is asked to
 * (the zero-flux flags of hyper_dg).
 *
 * @param inp Input parameters defined in gkyl_dg_euler_inp struct.
 * @return Pointer to euler equation object
 */
struct gkyl_dg_eqn *gkyl_dg_euler_inew(const struct gkyl_dg_euler_inp *inp);

/**
 * Create new Euler equation object on NV-GPU: see inew() method above for documentation.
 */
struct gkyl_dg_eqn *gkyl_dg_euler_cu_dev_inew(const struct gkyl_dg_euler_inp *inp);

/**
 * The flux through the absorbing wall at one edge, integrated over the wall face of each
 * boundary cell (with its physical area on a position map): the same flux as the boundary
 * surface term of the equation object. Component k of a cell is the integral of the flux of
 * component k of the state along the coordinate direction (negative when it leaves through a
 * lower wall, positive through an upper wall).
 *
 * @param eqn Euler equation object.
 * @param grid Configuration-space grid.
 * @param dir Direction of the wall.
 * @param edge -1 for the lower wall, +1 for the upper wall.
 * @param skin_range Range of the boundary cells of the edge.
 * @param fluid Physical fluid state.
 * @param flux_int Output: integrated wall flux of each boundary cell (num_equations components).
 */
void gkyl_dg_euler_wall_flux(
  const struct gkyl_dg_eqn *eqn, const struct gkyl_rect_grid *grid, int dir, int edge,
  const struct gkyl_range *skin_range, const struct gkyl_array *fluid, struct gkyl_array *flux_int
);

/**
 * Device version of gkyl_dg_euler_wall_flux.
 */
void gkyl_dg_euler_wall_flux_cu(
  const struct gkyl_dg_eqn *eqn, const struct gkyl_rect_grid *grid, int dir, int edge,
  const struct gkyl_range *skin_range, const struct gkyl_array *fluid, struct gkyl_array *flux_int
);
