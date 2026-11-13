// Private header: not for direct use
#pragma once

#include <assert.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_pkpm_remap_flux_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_supermesh_1d.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

typedef void (*pkpm_remap_flux_trace_t)(
  double jdinv_l, double jdinv_r, const double *f_l, const double *f_r, double *GKYL_RESTRICT G_l,
  double *GKYL_RESTRICT G_r
);
typedef void (*pkpm_remap_flux_store_t)(
  const double *wpsi, const double *proj, double *GKYL_RESTRICT flux
);

// for use in kernel tables
typedef struct {
  pkpm_remap_flux_trace_t kernels[3];
} gkyl_dg_pkpm_remap_flux_trace_kern_list;
typedef struct {
  pkpm_remap_flux_store_t kernels[3];
} gkyl_dg_pkpm_remap_flux_store_kern_list;

// Block traces of each side per velocity cell (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_remap_flux_trace_kern_list ser_pkpm_remap_flux_trace_kernels[] =
  {
    {pkpm_remap_flux_trace_x_1x1v_ser_p1, NULL, NULL}, // 0
    {pkpm_remap_flux_trace_x_2x1v_ser_p1, pkpm_remap_flux_trace_y_2x1v_ser_p1, NULL}, // 1
    {pkpm_remap_flux_trace_x_3x1v_ser_p1, pkpm_remap_flux_trace_y_3x1v_ser_p1,
     pkpm_remap_flux_trace_z_3x1v_ser_p1} // 2
};

// Store of Legendre projections into surface modes (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_remap_flux_store_kern_list ser_pkpm_remap_flux_store_kernels[] =
  {
    {pkpm_remap_flux_store_x_1x1v_ser_p1, NULL, NULL}, // 0
    {pkpm_remap_flux_store_x_2x1v_ser_p1, pkpm_remap_flux_store_y_2x1v_ser_p1, NULL}, // 1
    {pkpm_remap_flux_store_x_3x1v_ser_p1, pkpm_remap_flux_store_y_3x1v_ser_p1,
     pkpm_remap_flux_store_z_3x1v_ser_p1} // 2
};

// Traces of the frame at one transverse node of a face, and the map of one side to the face frame.
struct pkpm_remap_flux_side {
  double rho, u[3], xi, b[3], theta; // traces at the node
  struct pkpm_harm_face_map map; // the map to the face frame
  double energy_g; // K_inv/r_inv^2, the factor of G in the block G frame-changed back to this side
  double cc_inv; // |c_perp_inv|^2
};

// The frame and the coefficients of the flux at one transverse node of a face.
struct pkpm_remap_flux_node {
  struct pkpm_remap_flux_side side[2];
  double rho_bar, u_bar[3], xi_bar, b_bar[3]; // the face frame
  double fc[PKPM_FC_NUM]; // face coefficients (pkpm_harm_face_coef)
  double kinks[3]; // sorted face-frame labels where the HLL weights have kinks
  int num_kinks;
};

// The cell a side's sweep is currently in, at one transverse node. Y and P are stored [mode][block].
struct pkpm_remap_flux_cell {
  int k; // velocity cell of the side (-1: outside its grid)
  double w_center, dw; // center and width of the cell
  double Y[3][8]; // block traces in the face frame, X = Y0 + Y1 eta + Y2 (3 eta^2 - 1)
  double P[3][8]; // accumulated Legendre projection of the face flux onto the cell, face frame
};

// The sum over a face of the flux weighted for the shares, per receiving side.
struct pkpm_remap_flux_share {
  double momentum_out[2], momentum_total[2]; // weight |w_s F|
  double energy_out[2], energy_total[2]; // weight w_s^2 F + 2 G
};

struct gkyl_dg_pkpm_remap_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int num_conf_basis; // Configuration-space basis size.
  int num_vel_cells; // Cells of the velocity grid.
  int num_surf_basis; // Surface modes per block, direction and side of the kinetic flux.
  int num_trans_modes; // Transverse modes of a face, the surface modes of the fluid flux.
  int num_nodes; // Transverse nodes of a face.
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map).
  pkpm_remap_flux_trace_t trace[3]; // Block traces in each direction.
  pkpm_remap_flux_store_t store[3]; // Store of projections in each direction.

  double *psi; // Transverse modes at the nodes, [node][mode].
  double *wpsi; // Node weight times the transverse modes, [node][mode].
  double *conf_trace[3]; // Conf basis at the face nodes per direction, [side][node][basis].
  double *conf_node; // Conf basis at the Gauss-Legendre nodes of a cell, [node][basis].
  double *vel_edges; // Edges of the velocity grid.

  // Workspace for one face.
  long *column; // Phase index of each velocity cell, [side][cell].
  double *block_traces; // Modal block traces, [side][cell][block][mode].
  double *face_flux; // Kinetic flux of the face, [side][cell][block][mode].
  struct gkyl_supermesh_1d_interval *mesh; // Supermesh of one node.
  struct pkpm_remap_flux_node *nodes; // Frame of every node.
};

// 3-point Gauss-Legendre on [-1, 1]: exact on every supermesh interval, the flux being at most cubic.
GKYL_CU_D static const double pkpm_remap_flux_gn3[3] = {
  -0.7745966692414834, 0.0, 0.7745966692414834
};
GKYL_CU_D static const double pkpm_remap_flux_gw3[3] = {
  0.5555555555555556, 0.8888888888888888, 0.5555555555555556
};
// Orthonormal Legendre modes at those nodes, [node][mode].
GKYL_CU_D static const double pkpm_remap_flux_leg3[3][3] = {
  {0.7071067811865475, -0.9486832980505138, 0.6324555320336759},
  {0.7071067811865475, 0.0, -0.7905694150420949},
  {0.7071067811865475, 0.9486832980505138, 0.6324555320336759}
};

// Coefficients of the quadratic A + B eta + C (3 eta^2 - 1) per orthonormal Legendre mode.
#define PKPM_REMAP_FLUX_LEG0 0.7071067811865475
#define PKPM_REMAP_FLUX_LEG1 1.224744871391589
#define PKPM_REMAP_FLUX_LEG2 0.7905694150420949

// Value at the node of a p=1 field from its modal coefficients.
GKYL_CU_DH static inline double
pkpm_remap_flux_trace_at(int num_basis, const double *table, const double *coeff)
{
  double val = 0.0;
  for (int m = 0; m < num_basis; ++m) {
    val += table[m] * coeff[m];
  }
  return val;
}

// Enter velocity cell k of a side at a node: the block traces in the face frame, and clear the
// accumulated projection.
GKYL_CU_DH static inline void
pkpm_remap_flux_cell_enter(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int side, int k, const double *psi,
  const struct pkpm_remap_flux_side *sd, struct pkpm_remap_flux_cell *cell
)
{
  int nv = up->num_vel_cells, nt = up->num_surf_basis, na = up->num_trans_modes;
  const double *G = &up->block_traces[(side * nv + k) * 8 * nt];
  double X[3][8];
  for (int c = 0; c < 8; ++c) {
    double q[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < na; ++a) {
      for (int m = 0; m < 3; ++m) {
        q[m] += psi[a] * G[c * nt + a * 3 + m];
      }
    }
    X[0][c] = PKPM_REMAP_FLUX_LEG0 * q[0];
    X[1][c] = PKPM_REMAP_FLUX_LEG1 * q[1];
    X[2][c] = PKPM_REMAP_FLUX_LEG2 * q[2];
  }
  for (int m = 0; m < 3; ++m) {
    pkpm_harm_frame_change(X[m], sd->map.K, sd->map.r, sd->map.c_perp, cell->Y[m]);
    for (int c = 0; c < 8; ++c) {
      cell->P[m][c] = 0.0;
    }
  }
  cell->k = k;
  cell->dw = up->vel_edges[k + 1] - up->vel_edges[k];
  cell->w_center = 0.5 * (up->vel_edges[k + 1] + up->vel_edges[k]);
}

// Leave the velocity cell of a side at a node: map the accumulated projection back to the side's
// frame and add it to the side's flux with the transverse weight.
GKYL_CU_DH static inline void
pkpm_remap_flux_cell_leave(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, const double *wpsi,
  const struct pkpm_remap_flux_side *sd, const struct pkpm_remap_flux_cell *cell, double *flux
)
{
  double Z[3][8];
  for (int m = 0; m < 3; ++m) {
    pkpm_harm_frame_change(cell->P[m], sd->map.K_inv, sd->map.r_inv, sd->map.c_perp_inv, Z[m]);
  }
  up->store[dir](wpsi, &Z[0][0], flux);
}

// One supermesh Gauss point: the HLL flux of the two sides' face-frame blocks at the label wb, with
// the flux projected per side, the face-frame moments and the shares accumulated.
GKYL_CU_DH static inline void
pkpm_remap_flux_point(
  int dir, double wb, double om, const struct pkpm_remap_flux_node *nd,
  struct pkpm_remap_flux_cell *cells, double *mom, struct pkpm_remap_flux_share *share
)
{
  double X[2][8], eta[2], w_side[2];
  for (int s = 0; s < 2; ++s) {
    const struct pkpm_remap_flux_cell *cell = &cells[s];
    w_side[s] = nd->side[s].map.r * wb + nd->side[s].map.c_par;
    if (cell->k < 0) {
      for (int c = 0; c < 8; ++c) {
        X[s][c] = 0.0;
      }
      continue;
    }
    eta[s] = 2.0 * (w_side[s] - cell->w_center) / cell->dw;
    double e2 = 3.0 * eta[s] * eta[s] - 1.0;
    for (int c = 0; c < 8; ++c) {
      X[s][c] = cell->Y[0][c] + cell->Y[1][c] * eta[s] + cell->Y[2][c] * e2;
    }
  }
  double phi_l[8], phi_r[8], bounds[2], flux[8];
  pkpm_harm_conf_flux(dir, X[0], X[1], wb, nd->fc, 1, phi_l, phi_r, bounds);
  for (int c = 0; c < 8; ++c) {
    flux[c] = phi_l[c] + phi_r[c];
  }

  mom[0] += om * flux[0];
  mom[1] += om * wb * flux[0];
  mom[2] += om * wb * wb * flux[0];
  mom[3] += om * flux[1];
  for (int k = 0; k < 3; ++k) {
    mom[4 + k] += om * flux[2 + k];
  }

  double mass = om * fabs(flux[0]);
  for (int s = 0; s < 2; ++s) {
    const struct pkpm_remap_flux_side *sd = &nd->side[s];
    struct pkpm_remap_flux_cell *cell = &cells[s];
    const double *ci = sd->map.c_perp_inv;
    double cM = ci[0] * flux[2] + ci[1] * flux[3] + ci[2] * flux[4];
    double y_f = sd->map.K_inv * flux[0];
    double y_g = sd->energy_g * (flux[1] - cM + 0.5 * sd->cc_inv * flux[0]);
    double w2 = w_side[s] * w_side[s];
    double momentum = mass * fabs(w_side[s]);
    double energy = om * (fabs(y_f) * w2 + 2.0 * fabs(y_g));
    share->momentum_total[s] += momentum;
    share->energy_total[s] += energy;
    if (cell->k < 0) {
      share->momentum_out[s] += momentum;
      share->energy_out[s] += energy;
      continue;
    }
    double wt = om * sd->map.r * 2.0 / cell->dw;
    double e = eta[s];
    double leg[3] = {
      PKPM_REMAP_FLUX_LEG0, PKPM_REMAP_FLUX_LEG1 * e, PKPM_REMAP_FLUX_LEG2 * (3.0 * e * e - 1.0)
    };
    for (int m = 0; m < 3; ++m) {
      for (int c = 0; c < 8; ++c) {
        cell->P[m][c] += wt * leg[m] * flux[c];
      }
    }
  }
}

// The largest share of the flux outside the receiving side's velocity grid.
GKYL_CU_DH static inline double
pkpm_remap_flux_blend(const struct pkpm_remap_flux_share *share)
{
  double momentum = 0.0, energy = 0.0;
  for (int s = 0; s < 2; ++s) {
    if (share->momentum_total[s] > 0.0) {
      momentum = fmax(momentum, share->momentum_out[s] / share->momentum_total[s]);
    }
    if (share->energy_total[s] > 0.0) {
      energy = fmax(energy, share->energy_out[s] / share->energy_total[s]);
    }
  }
  return fmax(momentum, energy);
}

// Thermal speed at one face of a node line of a cell, from the values of the thermal speed and the
// density at the line's two Gauss-Legendre nodes (lower, upper); sign is +1 for the upper face and
// -1 for the lower. tau = 1/xi^2 is extrapolated linearly to the face: exact where the pressure is
// uniform and the density linear. On the face beyond the hotter node the result is limited to the
// ratio of the pressure and density traces (p = rho xi^2 at the nodes, both linear along the line),
// and never below the hot node's own value. A line with a non-positive node value keeps its linear
// extrapolation.
GKYL_CU_DH static inline double
pkpm_remap_flux_xi_line(double xi_lo, double xi_up, double rho_lo, double rho_up, double sign)
{
  const double half_sqrt3 = 0.8660254037844386;
  if (!(xi_lo > 0.0 && xi_up > 0.0)) {
    return 0.5 * (xi_lo + xi_up) + sign * half_sqrt3 * (xi_up - xi_lo);
  }
  double tau_lo = 1.0 / (xi_lo * xi_lo), tau_up = 1.0 / (xi_up * xi_up);
  double tau_face = 0.5 * (tau_lo + tau_up) + sign * half_sqrt3 * (tau_up - tau_lo);
  double sign_hot = tau_lo < tau_up ? -1.0 : 1.0;
  if (sign == sign_hot) {
    double tau_hot = fmin(tau_lo, tau_up);
    double pres_lo = rho_lo / tau_lo, pres_up = rho_up / tau_up;
    double pres_face = 0.5 * (pres_lo + pres_up) + sign * half_sqrt3 * (pres_up - pres_lo);
    double rho_face = 0.5 * (rho_lo + rho_up) + sign * half_sqrt3 * (rho_up - rho_lo);
    double tau_min = (pres_face > 0.0 && rho_face > 0.0) ? fmin(tau_hot, rho_face / pres_face) :
                                                           tau_hot;
    tau_face = fmax(tau_face, tau_min);
  }
  return 1.0 / sqrt(tau_face);
}

// Thermal speed of a cell on one of its faces in direction dir (sign +1: upper, -1: lower), as the
// p=1 expansion xi_face whose trace on that face interpolates the face values of the cell's node
// lines along dir (pkpm_remap_flux_xi_line) across the face.
GKYL_CU_DH static inline void
pkpm_remap_flux_xi_face(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, double sign, const double *fluid,
  const double *xi, double *xi_face
)
{
  int nc = up->num_conf_basis, num_cell_nodes = 1 << up->cdim;
  double node_xi[8], node_rho[8], line[8];
  for (int q = 0; q < num_cell_nodes; ++q) {
    node_xi[q] = pkpm_remap_flux_trace_at(nc, &up->conf_node[q * nc], xi);
    node_rho[q] = pkpm_remap_flux_trace_at(nc, &up->conf_node[q * nc], fluid);
  }
  for (int q = 0; q < num_cell_nodes; ++q) {
    int q_up = q | (1 << dir);
    if (q_up != q) {
      line[q] =
        pkpm_remap_flux_xi_line(node_xi[q], node_xi[q_up], node_rho[q], node_rho[q_up], sign);
      line[q_up] = line[q];
    }
  }
  // the 2-point Gauss-Legendre rule (unit weights) is exact for the p=1 projection
  for (int m = 0; m < nc; ++m) {
    xi_face[m] = 0.0;
    for (int q = 0; q < num_cell_nodes; ++q) {
      xi_face[m] += up->conf_node[q * nc + m] * line[q];
    }
  }
}

// Fill the frame of one transverse node of a face; returns whether both sides admit the frame map.
GKYL_CU_DH static inline int
pkpm_remap_flux_node_traces(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, int node, const double *const *bvar,
  const double *const *fluid, const double *const *u, const double *const *xi,
  const double *const *theta, struct pkpm_remap_flux_node *nd
)
{
  int nc = up->num_conf_basis, nn = up->num_nodes, admissible = 1;
  for (int s = 0; s < 2; ++s) {
    const double *table = &up->conf_trace[dir][(s * nn + node) * nc];
    struct pkpm_remap_flux_side *sd = &nd->side[s];
    sd->rho = pkpm_remap_flux_trace_at(nc, table, fluid[s]);
    sd->xi = pkpm_remap_flux_trace_at(nc, table, xi[s]);
    sd->theta = pkpm_remap_flux_trace_at(nc, table, theta[s]);
    for (int k = 0; k < 3; ++k) {
      sd->u[k] = pkpm_remap_flux_trace_at(nc, table, &u[s][k * nc]);
      sd->b[k] = pkpm_remap_flux_trace_at(nc, table, &bvar[s][k * nc]);
    }
    admissible = admissible && sd->rho > 0.0 && sd->xi > 0.0;
  }
  const struct pkpm_remap_flux_side *sl = &nd->side[0], *sr = &nd->side[1];
  for (int k = 0; k < PKPM_FC_NUM; ++k) {
    nd->fc[k] = 0.0;
  }
  pkpm_harm_face_coef(
    dir, sl->u[dir], sr->u[dir], sl->xi, sr->xi, sl->b, sr->b, sl->theta, sr->theta, 1, nd->fc
  );
  return admissible;
}

// The face frame of a node and the map of each side to it, and the kinks of the HLL weights.
GKYL_CU_DH static inline void
pkpm_remap_flux_node_frame(struct pkpm_remap_flux_node *nd)
{
  const struct pkpm_remap_flux_side *sl = &nd->side[0], *sr = &nd->side[1];
  nd->rho_bar = 0.5 * (sl->rho + sr->rho);
  nd->xi_bar = 0.5 * (sl->xi + sr->xi);
  double bmag = 0.0;
  for (int k = 0; k < 3; ++k) {
    nd->u_bar[k] = 0.5 * (sl->u[k] + sr->u[k]);
    nd->b_bar[k] = 0.5 * (sl->b[k] + sr->b[k]);
    bmag += nd->b_bar[k] * nd->b_bar[k];
  }
  bmag = bmag > 0.0 ? 1.0 / sqrt(bmag) : 0.0;
  for (int k = 0; k < 3; ++k) {
    nd->b_bar[k] *= bmag;
  }
  for (int s = 0; s < 2; ++s) {
    struct pkpm_remap_flux_side *sd = &nd->side[s];
    pkpm_harm_face_map(
      sd->rho, sd->u, sd->xi, nd->rho_bar, nd->u_bar, nd->xi_bar, nd->b_bar, &sd->map
    );
    sd->energy_g = sd->map.K_inv / (sd->map.r_inv * sd->map.r_inv);
    const double *ci = sd->map.c_perp_inv;
    sd->cc_inv = ci[0] * ci[0] + ci[1] * ci[1] + ci[2] * ci[2];
  }
  double a0 = nd->fc[PKPM_FC_A0], a1 = nd->fc[PKPM_FC_A1], c = nd->fc[PKPM_FC_C];
  nd->num_kinks = 0;
  if (fabs(a1) > 0.0) {
    double z = -a0 / a1;
    if (c > 0.0) {
      double z1 = (-c - a0) / a1, z2 = (c - a0) / a1;
      nd->kinks[0] = fmin(z1, z2);
      nd->kinks[1] = z;
      nd->kinks[2] = fmax(z1, z2);
      nd->num_kinks = 3;
    } else {
      nd->kinks[0] = z;
      nd->num_kinks = 1;
    }
  }
}

// The sweep of one transverse node of a face over its supermesh: adds the remap flux of each side
// to the output and returns the face-frame moments in mom.
GKYL_CU_DH static inline void
pkpm_remap_flux_node_sweep(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, int node,
  const struct pkpm_remap_flux_node *nd, double *flux_out, int nflux, double *mom,
  struct pkpm_remap_flux_share *share
)
{
  int nv = up->num_vel_cells, na = up->num_trans_modes;
  const double *psi = &up->psi[node * na], *wpsi = &up->wpsi[node * na];
  double tol = 1e-13 * up->phase_grid.dx[up->cdim];
  int num_intervals = gkyl_supermesh_1d_merge(
    nv + 1, up->vel_edges, nd->side[0].map.r, nd->side[0].map.c_par, nv + 1, up->vel_edges,
    nd->side[1].map.r, nd->side[1].map.c_par, nd->num_kinks, nd->kinks, tol, up->mesh
  );
  for (int m = 0; m < 7; ++m) {
    mom[m] = 0.0;
  }
  struct pkpm_remap_flux_cell cells[2];
  cells[0].k = cells[1].k = -1;
  for (int q = 0; q < num_intervals; ++q) {
    const struct gkyl_supermesh_1d_interval *iv = &up->mesh[q];
    int k_side[2] = {iv->k_l, iv->k_r};
    for (int s = 0; s < 2; ++s) {
      if (k_side[s] == cells[s].k) {
        continue;
      }
      if (cells[s].k >= 0) {
        pkpm_remap_flux_cell_leave(
          up, dir, wpsi, &nd->side[s], &cells[s], &flux_out[(s * nv + cells[s].k) * nflux]
        );
      }
      if (k_side[s] >= 0) {
        pkpm_remap_flux_cell_enter(up, s, k_side[s], psi, &nd->side[s], &cells[s]);
      } else {
        cells[s].k = -1;
      }
    }
    double mid = 0.5 * (iv->lo + iv->up), half = 0.5 * (iv->up - iv->lo);
    for (int g = 0; g < 3; ++g) {
      pkpm_remap_flux_point(
        dir, mid + half * pkpm_remap_flux_gn3[g], half * pkpm_remap_flux_gw3[g], nd, cells, mom,
        share
      );
    }
  }
  for (int s = 0; s < 2; ++s) {
    if (cells[s].k >= 0) {
      pkpm_remap_flux_cell_leave(
        up, dir, wpsi, &nd->side[s], &cells[s], &flux_out[(s * nv + cells[s].k) * nflux]
      );
    }
  }
}

// Add blend times the label-paired flux of one transverse node to both sides' fluxes and return
// its fluid flux: the flux of the two sides paired at equal w, one velocity cell at a time, as
// the moment map of each one-sided contribution in that side's own frame.
GKYL_CU_DH static inline void
pkpm_remap_flux_node_label(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, int node, double blend,
  const struct pkpm_remap_flux_node *nd, double *flux_out, int nflux, double *fluid_flux
)
{
  int nv = up->num_vel_cells, nt = up->num_surf_basis, na = up->num_trans_modes;
  const double *psi = &up->psi[node * na];
  double blend_wpsi[4];
  for (int a = 0; a < na; ++a) {
    blend_wpsi[a] = blend * up->wpsi[node * na + a];
  }
  double moms[2][7];
  for (int s = 0; s < 2; ++s) {
    for (int m = 0; m < 7; ++m) {
      moms[s][m] = 0.0;
    }
  }
  for (int k = 0; k < nv; ++k) {
    double dw = up->vel_edges[k + 1] - up->vel_edges[k];
    double w_center = 0.5 * (up->vel_edges[k + 1] + up->vel_edges[k]);
    double q[2][8][3];
    for (int s = 0; s < 2; ++s) {
      const double *G = &up->block_traces[(s * nv + k) * 8 * nt];
      for (int c = 0; c < 8; ++c) {
        for (int m = 0; m < 3; ++m) {
          q[s][c][m] = 0.0;
          for (int a = 0; a < na; ++a) {
            q[s][c][m] += psi[a] * G[c * nt + a * 3 + m];
          }
        }
      }
    }
    double proj[3][8];
    for (int m = 0; m < 3; ++m) {
      for (int c = 0; c < 8; ++c) {
        proj[m][c] = 0.0;
      }
    }
    for (int j = 0; j < 3; ++j) {
      const double *leg = pkpm_remap_flux_leg3[j];
      double wq = w_center + 0.5 * pkpm_remap_flux_gn3[j] * dw;
      double X[2][8];
      for (int s = 0; s < 2; ++s) {
        for (int c = 0; c < 8; ++c) {
          X[s][c] = q[s][c][0] * leg[0] + q[s][c][1] * leg[1] + q[s][c][2] * leg[2];
        }
      }
      double phi_l[8], phi_r[8], bounds[2];
      pkpm_harm_conf_flux(dir, X[0], X[1], wq, nd->fc, 1, phi_l, phi_r, bounds);
      for (int m = 0; m < 3; ++m) {
        for (int c = 0; c < 8; ++c) {
          proj[m][c] += pkpm_remap_flux_gw3[j] * leg[m] * (phi_l[c] + phi_r[c]);
        }
      }
      double wgt = 0.5 * pkpm_remap_flux_gw3[j] * dw;
      const double *phi[2] = {phi_l, phi_r};
      for (int s = 0; s < 2; ++s) {
        moms[s][0] += wgt * phi[s][0];
        moms[s][1] += wgt * wq * phi[s][0];
        moms[s][2] += wgt * wq * wq * phi[s][0];
        moms[s][3] += wgt * phi[s][1];
        for (int a = 0; a < 3; ++a) {
          moms[s][4 + a] += wgt * phi[s][2 + a];
        }
      }
    }
    for (int s = 0; s < 2; ++s) {
      up->store[dir](blend_wpsi, &proj[0][0], &flux_out[(s * nv + k) * nflux]);
    }
  }
  for (int m = 0; m < 5; ++m) {
    fluid_flux[m] = 0.0;
  }
  for (int s = 0; s < 2; ++s) {
    const struct pkpm_remap_flux_side *sd = &nd->side[s];
    double side_flux[5];
    pkpm_harm_fluid_moment_map(moms[s], sd->rho, sd->u, sd->xi, sd->b, side_flux);
    for (int m = 0; m < 5; ++m) {
      fluid_flux[m] += side_flux[m];
    }
  }
}
