
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_dg_vlasov_pressure gkyl_dg_vlasov_pressure;

/**
 * New updater to compute pressure quantities from
 * calculated velocity moments.
 *
 * @param conf_basis Configuration-space basis functions
 * @param conf_range Configuration-space range
 * @param conf_range_ext Extended configuration-space range
 * @param vdim Number of velocity dimensions
 * @param mass Species mass
 * @param use_gpu bool for gpu useage
 * @return New updater pointer
 */
struct gkyl_dg_vlasov_pressure * gkyl_dg_vlasov_pressure_new(
  const struct gkyl_basis *conf_basis,
  const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext,
  int vdim, double mass, bool use_gpu
);

/**
 * Compute the scalar pressure from the velocity moments:
 * P = (m/vdim) * (M2 - V_drift dot M1i)
 *
 * @param up Pressure updater
 * @param conf_local Configuration-space range
 * @param M2 Scalar second moment
 * @param M1i First moment vector
 * @param V_drift Drift velocity vector
 * @param pressure Output scalar pressure
 */
void gkyl_dg_vlasov_pressure_scalar(
  struct gkyl_dg_vlasov_pressure *up,
  const struct gkyl_range *conf_local,
  const struct gkyl_array *M2,
  const struct gkyl_array *M1i,
  const struct gkyl_array *V_drift,
  struct gkyl_array *pressure
);

/**
 * Compute the pressure tensor from the velocity moments:
 * P_ij = m * (M2ij - V_i*M1j)
 *
 * @param up Pressure updater
 * @param conf_local Configuration-space range
 * @param M2ij Second moment tensor
 * @param M1i First moment vector
 * @param V_drift Drift velocity vector
 * @param pressure_tensor Output pressure tensor
 */
void gkyl_dg_vlasov_pressure_tensor(
  struct gkyl_dg_vlasov_pressure *up,
  const struct gkyl_range *conf_local,
  const struct gkyl_array *M2ij,
  const struct gkyl_array *M1i,
  const struct gkyl_array *V_drift,
  struct gkyl_array *pressure_tensor
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_vlasov_pressure_release(struct gkyl_dg_vlasov_pressure *up);
