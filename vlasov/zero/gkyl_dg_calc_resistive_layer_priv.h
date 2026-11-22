// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_dg_calc_resistive_layer.h>
#include <gkyl_util.h>

struct gkyl_dg_calc_resistive_layer {
  int num_basis; // number of configuration-space basis functions
  int tot_quad; // number of quadrature points per cell
  struct gkyl_array *weights; // quadrature weights
  struct gkyl_array *basis_at_ords; // basis functions at the quadrature points
  struct gkyl_array *sigma; // sigma at the quadrature points of each cell
  bool use_gpu;
};

// Accumulate -(sigma E) onto the three electric field expansions of rhs in one cell.
GKYL_CU_DH static inline void
dg_calc_resistive_layer_cell(
  int num_basis, int tot_quad, const double *GKYL_RESTRICT weights,
  const double *GKYL_RESTRICT basis_at_ords, const double *GKYL_RESTRICT sigma,
  const double *GKYL_RESTRICT em, double *GKYL_RESTRICT rhs
)
{
  for (int c = 0; c < 3; ++c) {
    for (int q = 0; q < tot_quad; ++q) {
      const double *b_ord = &basis_at_ords[q * num_basis];
      double E_q = 0.0;
      for (int k = 0; k < num_basis; ++k) {
        E_q += em[c * num_basis + k] * b_ord[k];
      }
      double f_q = weights[q] * sigma[q] * E_q;
      for (int k = 0; k < num_basis; ++k) {
        rhs[c * num_basis + k] -= f_q * b_ord[k];
      }
    }
  }
}

#ifdef GKYL_HAVE_CUDA
void gkyl_dg_calc_resistive_layer_advance_cu(
  const struct gkyl_dg_calc_resistive_layer *up, const struct gkyl_range *range,
  const struct gkyl_array *em, struct gkyl_array *rhs
);
#endif
