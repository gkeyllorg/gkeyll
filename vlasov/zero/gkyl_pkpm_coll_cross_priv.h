// Private header: not for direct use
#pragma once

#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_util.h>

// components of the fluid state the exchange changes, [M (3), E]
#define PKPM_COLL_CROSS_NUM_INCR (PKPM_FLUID_NUM_COMP - PKPM_FLUID_MOM)

struct gkyl_pkpm_coll_cross {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double heat_fac; // kinetic rate nu_sr (1 + heat_fac (vtsq_sr - vtsq_s) / vtsq_s)
  dg_conf_nodal_kern_t conf_to_nodal; // expansion to node values
  dg_conf_nodal_kern_t conf_to_modal; // node values to expansion
};

// A = L U in place, without pivoting (for the diagonally dominant matrices of the exchange).
GKYL_CU_DH static inline void
pkpm_coll_cross_factor(int num, double mat[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES])
{
  for (int k = 0; k < num; ++k) {
    for (int i = k + 1; i < num; ++i) {
      mat[i][k] /= mat[k][k];
      for (int j = k + 1; j < num; ++j) {
        mat[i][j] -= mat[i][k] * mat[k][j];
      }
    }
  }
}

// Solve L U x = rhs in place with the factors of pkpm_coll_cross_factor.
GKYL_CU_DH static inline void
pkpm_coll_cross_solve(
  int num, const double mat[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES], double rhs[GKYL_MAX_SPECIES]
)
{
  for (int i = 1; i < num; ++i) {
    for (int j = 0; j < i; ++j) {
      rhs[i] -= mat[i][j] * rhs[j];
    }
  }
  for (int i = num - 1; i >= 0; --i) {
    for (int j = i + 1; j < num; ++j) {
      rhs[i] -= mat[i][j] * rhs[j];
    }
    rhs[i] /= mat[i][i];
  }
}
