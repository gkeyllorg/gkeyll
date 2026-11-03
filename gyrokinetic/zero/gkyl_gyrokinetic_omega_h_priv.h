#pragma once

#include <gkyl_gyrokinetic_omega_h.h>

struct gkyl_gyrokinetic_omega_h {
  struct gkyl_range range;
  struct gkyl_array *basis_at_nodes; // Basis values, one sample per array element.
  struct gkyl_array *weight_sq; // G^2 at each sample, one configuration cell per element.
  double *max_dev; // Scalar device reduction result; NULL on CPU.
  bool use_gpu;
};

/**
 * Evaluate the local product at each cached point and return its maximum.
 *
 * @param basis_at_nodes Basis evaluations, including the cell center.
 * @param weight_sq Cached G^2 at the evaluation points of this cell.
 * @param density Modal Jn coefficients for this cell.
 * @return Nonnegative maximum sampled Jn*G^2 in this cell.
 */
GKYL_CU_DH static inline double
gyrokinetic_omega_h_cell_max(
  const struct gkyl_array *basis_at_nodes, const double *weight_sq, const double *density
)
{
  double max_value = 0.0;
  for (size_t node = 0; node < basis_at_nodes->size; ++node) {
    const double *basis_val = (const double *)gkyl_array_cfetch(basis_at_nodes, node);
    double density_nodal = 0.0;
    for (size_t mode = 0; mode < basis_at_nodes->ncomp; ++mode) {
      density_nodal += density[mode] * basis_val[mode];
    }
    max_value = fmax(max_value, fmax(0.0, density_nodal) * weight_sq[node]);
  }
  return max_value;
}

#ifdef GKYL_HAVE_CUDA
/**
 * Evaluate and reduce the local product on the GPU, copying only the result back.
 *
 * @param up Host updater containing device arrays and reduction storage.
 * @param density Device array of modal Jn coefficients.
 * @return Local maximum on the host.
 */
double gkyl_gyrokinetic_omega_h_advance_cu(
  struct gkyl_gyrokinetic_omega_h *up, const struct gkyl_array *density
);
#endif
