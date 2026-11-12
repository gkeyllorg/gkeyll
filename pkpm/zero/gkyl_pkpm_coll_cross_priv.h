// Private header: not for direct use
#pragma once

#include <gkyl_dg_conf_nodal_kernels.h>
#include <gkyl_util.h>

#define PKPM_COLL_CROSS_MAX_NODES 8
#define PKPM_COLL_CROSS_NUM_FLUID 5 // blocks of the fluid state [rho, M (3), E]
#define PKPM_COLL_CROSS_NUM_PRIM 5 // node fields [n, u (3), vtsq]
#define PKPM_COLL_CROSS_NUM_INCR 4 // blocks of the fluid state the exchange changes, [M (3), E]

typedef void (*pkpm_coll_cross_conf_map_t)(const double *in, double *GKYL_RESTRICT out);

GKYL_CU_D static const pkpm_coll_cross_conf_map_t pkpm_coll_cross_conf_to_nodal_list[] = {
  dg_conf_nodal_conf_to_nodal_1x_ser_p1, dg_conf_nodal_conf_to_nodal_2x_ser_p1,
  dg_conf_nodal_conf_to_nodal_3x_ser_p1
};

GKYL_CU_D static const pkpm_coll_cross_conf_map_t pkpm_coll_cross_nodal_to_conf_list[] = {
  dg_conf_nodal_conf_to_modal_1x_ser_p1, dg_conf_nodal_conf_to_modal_2x_ser_p1,
  dg_conf_nodal_conf_to_modal_3x_ser_p1
};

struct gkyl_pkpm_coll_cross {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double heat_fac; // kinetic rate nu_sr (1 + heat_fac (vtsq_sr - vtsq_s) / vtsq_s)
  pkpm_coll_cross_conf_map_t conf_to_nodal; // modal to nodal, for this cdim
  pkpm_coll_cross_conf_map_t nodal_to_conf; // nodal to modal, for this cdim
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
