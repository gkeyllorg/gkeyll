#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>
#include <gkyl_util.h>

// Object type
typedef struct gkyl_dg_calc_gk_rad_vars gkyl_dg_calc_gk_rad_vars;

// Struct to store drag coefficients for each
// cross-species collisions and for each density.
struct gkyl_gk_rad_drag {
  struct gkyl_array *arr; // Array with drag coeff.
  int num_dens; // Number of densities drag coeffs computed for.
  struct gkyl_gk_rad_drag *on_dev; // on_dev.arr has on_dev pointers.
  struct gkyl_gk_rad_drag *data; // data.arr has host pointer with array data and ->on_dev on the GPU (use_gpu=true).
};

/**
 * Allocate a new struct storing a drag coefficient for each cross-collision and for each density.
 *
 * @param num_collisions Number of radiating collisions.
 * @param num_densities Number of densities fits are calculate for, for each species.
 * @param num_comp Number of DG components/basis functions
 * @param sz Number of cells/size
 * @param use_gpu Whether to store drag coefficient on the gpu.
 */
struct gkyl_gk_rad_drag* gkyl_dg_calc_gk_rad_vars_drag_new(int num_collisions,
  const int *num_densities, int ncomp, long sz, bool use_gpu);

/**
 * Free memory associated with a gkyl_gk_rad_drag struct.
 *
 * @param drag Struct with drag coefficient for each cross-collision and density.
 * @param num_collision Number of radiating collisions.
 * @param use_gpu Whether data was stored on the GPU.
 */
void gkyl_dg_calc_gk_rad_vars_drag_release(struct gkyl_gk_rad_drag *drag,
  int num_collisions, bool use_gpu);

/**
 * Create new updater to compute the drag coefficients needed for 
 * radiation in gyrokinetic equations. Method computes:
 *    1. nu = nu(v) for both vparallel and mu updates, 
 *       where |v| = sqrt(v_par^2 + 2 mu B/m)
 *    Note that through the spatial variation of B = B(x,z), 
 *    both these drag coefficients depend on phase space, but a reduced (x,z,vpar,mu) phase space
 *    2. nvnu = sum_s n_{i_s} nu_s(v)
 * 
 * Both of these methods also return both the surface expansions 
 * (evaluated at constant vparallel and constant mu) and volume expansions
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param charge Species charge
 * @param mass Species mass
 * @param gk_geom Geometry struct
 * @param vel_map Velocity space mapping object.
 * @return New updater pointer.
 */
struct gkyl_dg_calc_gk_rad_vars* 
gkyl_dg_calc_gk_rad_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, double charge,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, bool use_gpu);

/**
 * Compute drag coefficients needed for radiation in gyrokinetic equations
 * 
 * @param up Updater for computing gyrokinetic radiation variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param a, alpha, beta, gamma, v0 Input fitting parameters for a given collision type
 * @param vnu_surf Output surface expansion of vpar drag coefficient
 * @param vnu Output volume expansion of vpar drag coefficient
 * @param vsqnu_surf Output surface expansion of mu drag coefficient
 * @param vsqnu Output volume expansion of mu drag coefficient
 */
void gkyl_dg_calc_gk_rad_vars_nu_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  double a, double alpha, double beta, double gamma, double v0, 
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu);

/**
 * Compute sum_s n_{i_s} nu_s(v) total drag coefficient for drag due to radiation in gyrokinetic equations
 * Sum is handled by repeated calling of this updater for each ion density and set of drag coefficients
 * 
 * @param up Updater for computing gyrokinetic radiation variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param vnu_surf Input surface expansion of vpar drag coefficient for a given collision producing radiation drag
 * @param vnu Input volume expansion of vpar drag coefficient for a given collision producing radiation drag
 * @param vsqnu_surf Input surface expansion of mu drag coefficient for a given collision producing radiation drag
 * @param vsqnu Input volume expansion of mu drag coefficient for a given collision producing radiation drag
 * @param n_elc_rad Array of electron densities where the vnuXXX arrays lie
 * @param n_elc Input volume expansion of electron density
 * @param nI Input volume expansion of ion density for a given collision producing radiation drag
 * @param nvnu_surf Output surface expansion of vpar component of sum_s n_{i_s} nu_s(v)
 * @param nvnu Output volume expansion of vpar component of sum_s n_{i_s} nu_s(v)
 * @param nvsqnu_surf Output surface expansion of mu component of sum_s n_{i_s} nu_s(v)
 * @param nvsqnu Output volume expansion of mu drag component of sum_s n_{i_s} nu_s(v)
 * @param vtsq_min_normalized Minimum vtsq for each fit divided by configuration space normalization
 * @param vtsq vtsq
 */
void gkyl_dg_calc_gk_rad_vars_nI_nu_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_gk_rad_drag* vnu_surf, const struct gkyl_gk_rad_drag* vnu, 
  const struct gkyl_gk_rad_drag* vsqnu_surf, const struct gkyl_gk_rad_drag* vsqnu,
  const struct gkyl_array* n_elc_rad, const struct gkyl_array* n_elc,
  const struct gkyl_array* nI, struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu,
  struct gkyl_array* vtsq_min_normalized, struct gkyl_array* vtsq);

/**
 * Delete pointer to updater to compute gyrokinetic variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_gk_rad_vars_release(struct gkyl_dg_calc_gk_rad_vars *up);
