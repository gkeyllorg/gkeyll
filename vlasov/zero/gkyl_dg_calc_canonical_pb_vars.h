#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_canonical_pb_vars gkyl_dg_calc_canonical_pb_vars;

/**
 * Create new updater to compute canonical_pb variables needed in 
 * updates for general geometry. Methods compute:
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_canonical_pb_vars* 
gkyl_dg_calc_canonical_pb_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,  bool use_gpu);

/**
 * Create new updater to compute canonical_pb general geometry variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_canonical_pb_vars* 
gkyl_dg_calc_canonical_pb_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis);


/**
 * Host-side wrappers for canonical_pb general geometry variable operations on device
 */

/**
 * Convert the contravaraint components to covariant components of the momentum
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing general geometry canonical_pb variables 
 * @param conf_range Configuration space range 
 * @param h_ij Covariant metric tensor expansion coefficients
 * @param V_drift Drift velocity moment expansion (contravariant components)
 * @param M1i Drift velocity times density moment expansion (contravariant components)
 * @param V_drift_cov Drift velocity moment expansion (covariant components)
 * @param M1i_cov Drift velocity times density moment expansion (covariant components)
 */
void gkyl_canonical_pb_contra_to_covariant_m1i(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij, const struct gkyl_array *V_drift, const struct gkyl_array *M1i, struct gkyl_array *V_drift_cov, 
 struct gkyl_array *M1i_cov);

 /**
 * Host-side wrappers for canonical_pb pressure operations on device
 */
void gkyl_canonical_pb_contra_to_covariant_m1i_cu(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij, const struct gkyl_array *V_drift, const struct gkyl_array *M1i, struct gkyl_array *V_drift_cov, 
 struct gkyl_array *M1i_cov);

/**
 * Compute the pressure moment from the energy and velocity moments
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing general geometry canonical_pb variables 
 * @param conf_range Configuration space range 
 * @param h_ij_inv Inverse of the metric tensor expansion
 * @param MEnergy MEnergy moment expansion: int(Hfd^3p)
 * @param V_drift Drift velocity moment expansion (contravariant components)
 * @param M1i Drift velocity times density moment expansion (contravariant components)
 * @param pressure Output, scalar pressure expansion
 */
void gkyl_canonical_pb_pressure(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij_inv, 
 const struct gkyl_array *MEnergy, const struct gkyl_array *V_drift, const struct gkyl_array *M1i,
 struct gkyl_array *pressure);

 /**
 * Host-side wrappers for canonical_pb pressure operations on device
 */
void gkyl_canonical_pb_pressure_cu(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij_inv, 
 const struct gkyl_array *MEnergy, const struct gkyl_array *V_drift, const struct gkyl_array *M1i,
 struct gkyl_array *pressure);

/**
 * Delete pointer to updater to compute canonical_pb general geometry variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_canonical_pb_vars_release(struct gkyl_dg_calc_canonical_pb_vars *up);