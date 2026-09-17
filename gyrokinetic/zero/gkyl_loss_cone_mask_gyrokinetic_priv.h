// Private header: not for direct use
#pragma once

#include <float.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

GKYL_CU_DH static inline void nodal_coords(int ndim, int node, double *x)
{
  for (int d = 0; d < ndim; ++d) {
    int place_value = (int)pow(2.0, (double)(ndim - 1 - d));
    int digit = (node / place_value) % 2;
    x[d] = digit ? 1.0 : -1.0;
  }
}

struct gkyl_loss_cone_mask_gyrokinetic {
  int cdim; // Configuration-space dimension.
  int num_basis_conf; // Number of configuration-space basis functions.
  int num_nodes_conf; // Number of configuration-space cell nodes.
  bool use_gpu;

  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  double mass; // Species mass.
  double charge; // Species charge.
  enum gkyl_loss_cone_boundary_type lower_boundary; // Lower parallel-boundary behavior.
  enum gkyl_loss_cone_boundary_type upper_boundary; // Upper parallel-boundary behavior.

  struct gkyl_array *basis_at_nodes_conf; // Basis functions at configuration-space nodes.
};

GKYL_CU_DH static inline int conf_node_z_endpoint_index(int cdim, int conf_node, int zdim)
{
  double eta[GKYL_MAX_DIM] = {0.0};
  nodal_coords(cdim, conf_node, eta);
  return eta[zdim] > 0.0 ? 1 : 0;
}

GKYL_CU_DH static inline int conf_node_with_matching_perpendicular_coords(
  int cdim, int anchor_node, int zdim, int z_endpoint_index
)
{
  int num_nodes_conf = (int)pow(2.0, (double)cdim);

  double eta_anchor[GKYL_MAX_DIM] = {0.0};
  nodal_coords(cdim, anchor_node, eta_anchor);

  for (int cand = 0; cand < num_nodes_conf; ++cand) {
    if (conf_node_z_endpoint_index(cdim, cand, zdim) != z_endpoint_index) {
      continue;
    }

    double eta_cand[GKYL_MAX_DIM] = {0.0};
    nodal_coords(cdim, cand, eta_cand);

    bool same_transverse = true;
    for (int d = 0; d < cdim; ++d) {
      if (d == zdim) {
        continue;
      }
      if (eta_cand[d] != eta_anchor[d]) {
        same_transverse = false;
        break;
      }
    }

    if (same_transverse) {
      return cand;
    }
  }

  return anchor_node;
}

GKYL_CU_DH static inline double field_node_val(
  const struct gkyl_array *arr, const struct gkyl_array *basis_at_nodes, int num_basis, long linidx,
  int node
)
{
  const double *arr_d = (const double *)gkyl_array_cfetch(arr, linidx);
  const double *basis_d = (const double *)gkyl_array_cfetch(basis_at_nodes, node);

  double val = 0.0;
  for (int k = 0; k < num_basis; ++k) {
    val += arr_d[k] * basis_d[k];
  }

  return val;
}

GKYL_CU_DH static inline bool hamiltonian_below_barrier(
  double kinetic_energy, double magnetic_energy, double electric_energy, double barrier
)
{
  if (barrier == DBL_MAX) {
    return true;
  }

  double hamiltonian = kinetic_energy + magnetic_energy + electric_energy;
  double energy_scale =
    fabs(kinetic_energy) + fabs(magnetic_energy) + fabs(electric_energy) + fabs(barrier);
  double roundoff_tol = 32.0 * DBL_EPSILON * energy_scale;

  // Treat equality, including equality to within modal-evaluation roundoff,
  // as reaching the boundary rather than reflecting from it.
  return hamiltonian < barrier - roundoff_tol;
}

GKYL_CU_DH static inline void escape_barriers(
  int cdim, int num_basis_conf, const struct gkyl_range *conf_range,
  const struct gkyl_array *basis_at_nodes_conf, const struct gkyl_array *phi,
  const struct gkyl_array *bmag, const struct gkyl_array *phi_wall_lo,
  const struct gkyl_array *phi_wall_up, const int *base_idx, int z_cell, int anchor_conf_node,
  double mu, double charge, enum gkyl_loss_cone_boundary_type lower_boundary,
  enum gkyl_loss_cone_boundary_type upper_boundary, double *barrier_left, double *barrier_right
)
{
  int zdim = cdim - 1;

  int z_endpoint_index = conf_node_z_endpoint_index(cdim, anchor_conf_node, zdim);
  int z_upper_node = conf_node_with_matching_perpendicular_coords(cdim, anchor_conf_node, zdim, 1);
  int z_lower_node = conf_node_with_matching_perpendicular_coords(cdim, anchor_conf_node, zdim, 0);

  int scan_idx[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    scan_idx[d] = base_idx[d];
  }

  *barrier_left = -DBL_MAX;
  *barrier_right = -DBL_MAX;

  for (int iz = conf_range->lower[zdim]; iz <= conf_range->upper[zdim]; ++iz) {
    scan_idx[zdim] = iz;
    long linidx = gkyl_range_idx(conf_range, scan_idx);

    double phi_lower =
      field_node_val(phi, basis_at_nodes_conf, num_basis_conf, linidx, z_lower_node);
    double phi_upper =
      field_node_val(phi, basis_at_nodes_conf, num_basis_conf, linidx, z_upper_node);
    double bmag_lower =
      field_node_val(bmag, basis_at_nodes_conf, num_basis_conf, linidx, z_lower_node);
    double bmag_upper =
      field_node_val(bmag, basis_at_nodes_conf, num_basis_conf, linidx, z_upper_node);
    double u_lower = mu * bmag_lower + charge * phi_lower;
    double u_upper = mu * bmag_upper + charge * phi_upper;

    // A complete cell contributes both endpoint traces. In the anchor cell,
    // only the portion between the anchor node and the relevant wall belongs
    // to that directional path. Keeping both traces also handles discontinuous
    // input fields conservatively.
    if (iz < z_cell) {
      *barrier_left = GKYL_MAX2(*barrier_left, GKYL_MAX2(u_lower, u_upper));
    } else if (iz == z_cell) {
      *barrier_left = GKYL_MAX2(*barrier_left, u_lower);
      *barrier_right = GKYL_MAX2(*barrier_right, u_upper);
      if (z_endpoint_index == 1) {
        *barrier_left = GKYL_MAX2(*barrier_left, u_upper);
      } else {
        *barrier_right = GKYL_MAX2(*barrier_right, u_lower);
      }
    } else {
      *barrier_right = GKYL_MAX2(*barrier_right, GKYL_MAX2(u_lower, u_upper));
    }
  }

  if (lower_boundary == GKYL_LOSS_CONE_BC_CLOSED) {
    *barrier_left = DBL_MAX;
  } else if (lower_boundary == GKYL_LOSS_CONE_BC_SHEATH) {
    scan_idx[zdim] = conf_range->lower[zdim];
    long linidx = gkyl_range_idx(conf_range, scan_idx);
    double bmag_wall =
      field_node_val(bmag, basis_at_nodes_conf, num_basis_conf, linidx, z_lower_node);
    double phi_wall =
      field_node_val(phi_wall_lo, basis_at_nodes_conf, num_basis_conf, linidx, z_lower_node);
    *barrier_left = GKYL_MAX2(*barrier_left, mu * bmag_wall + charge * phi_wall);
  }

  if (upper_boundary == GKYL_LOSS_CONE_BC_CLOSED) {
    *barrier_right = DBL_MAX;
  } else if (upper_boundary == GKYL_LOSS_CONE_BC_SHEATH) {
    scan_idx[zdim] = conf_range->upper[zdim];
    long linidx = gkyl_range_idx(conf_range, scan_idx);
    double bmag_wall =
      field_node_val(bmag, basis_at_nodes_conf, num_basis_conf, linidx, z_upper_node);
    double phi_wall =
      field_node_val(phi_wall_up, basis_at_nodes_conf, num_basis_conf, linidx, z_upper_node);
    *barrier_right = GKYL_MAX2(*barrier_right, mu * bmag_wall + charge * phi_wall);
  }
}

#ifdef GKYL_HAVE_CUDA
#ifdef __cplusplus
extern "C" {
#endif

void gkyl_loss_cone_mask_gyrokinetic_advance_cu(
  gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out
);

#ifdef __cplusplus
}
#endif
#endif
