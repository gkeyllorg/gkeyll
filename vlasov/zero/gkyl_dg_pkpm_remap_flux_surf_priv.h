// Private header: not for direct use
#pragma once

#include <assert.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_pkpm_layout.h>
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

typedef void (*pkpm_remap_flux_interval_t)(
  const double *fc, double w_lo, double w_up, const double *side_l, const double *side_r,
  const double *cell_l, const double *cell_r, const double *Y_l, const double *Y_r,
  double *GKYL_RESTRICT P_l, double *GKYL_RESTRICT P_r, double *GKYL_RESTRICT mom,
  double *GKYL_RESTRICT share
);

typedef void (*pkpm_remap_flux_to_face_t)(
  const double *psi, const double *fwd, const double *G, double *GKYL_RESTRICT X,
  double *GKYL_RESTRICT Y
);
typedef void (*pkpm_remap_flux_from_face_t)(
  const double *wpsi, const double *bck, const double *P, double *GKYL_RESTRICT flux
);
typedef void (*pkpm_remap_flux_label_t)(
  const double *fc, double w_center, double dw, const double *X_l, const double *X_r,
  double *GKYL_RESTRICT proj, double *GKYL_RESTRICT mom_l, double *GKYL_RESTRICT mom_r,
  double *GKYL_RESTRICT speed
);

// for use in kernel tables
typedef struct {
  pkpm_remap_flux_trace_t kernels[3];
} gkyl_dg_pkpm_remap_flux_trace_kern_list;
typedef struct {
  pkpm_remap_flux_from_face_t kernels[3];
} gkyl_dg_pkpm_remap_flux_from_face_kern_list;
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

// Blocks of one velocity cell of a side at a transverse node, in the frame of the side and in
// the face frame (p=1 serendipity only; indexed [cdim-1])
GKYL_CU_D static const pkpm_remap_flux_to_face_t ser_pkpm_remap_flux_to_face_kernels[] = {
  pkpm_remap_flux_to_face_1x1v_ser_p1, pkpm_remap_flux_to_face_2x1v_ser_p1,
  pkpm_remap_flux_to_face_3x1v_ser_p1
};

// Face-frame projection of one velocity cell of a side at a transverse node back in the frame of
// the side, into the surface modes (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_remap_flux_from_face_kern_list
  ser_pkpm_remap_flux_from_face_kernels[] = {
    {pkpm_remap_flux_from_face_x_1x1v_ser_p1, NULL, NULL}, // 0
    {pkpm_remap_flux_from_face_x_2x1v_ser_p1, pkpm_remap_flux_from_face_y_2x1v_ser_p1, NULL}, // 1
    {pkpm_remap_flux_from_face_x_3x1v_ser_p1, pkpm_remap_flux_from_face_y_3x1v_ser_p1,
     pkpm_remap_flux_from_face_z_3x1v_ser_p1} // 2
};

// One interval of the velocity supermesh of a transverse node (any cdim; indexed [dir]): the flux
// at its 3 Gauss-Legendre points (exact, the flux being at most cubic there), projected onto the
// cell of each side, with its face-frame moments and its shares for the blend.
GKYL_CU_D static const pkpm_remap_flux_interval_t pkpm_remap_flux_interval_kernels[] = {
  pkpm_remap_flux_interval_x_1v_p2, pkpm_remap_flux_interval_y_1v_p2,
  pkpm_remap_flux_interval_z_1v_p2
};

// Label-paired flux of one velocity cell at a transverse node (any cdim; indexed [dir]): the flux
// of the two sides paired at equal w at the 3 Gauss-Legendre points of the cell, and its speed
// bounds there, for the rates.
GKYL_CU_D static const pkpm_remap_flux_label_t pkpm_remap_flux_label_kernels[] = {
  pkpm_remap_flux_label_x_1v_p2, pkpm_remap_flux_label_y_1v_p2, pkpm_remap_flux_label_z_1v_p2
};

// Legendre coefficients in w of the eight blocks of one velocity cell at a node, [mode][block].
#define PKPM_REMAP_FLUX_NUM_COEFF (PKPM_NUM_VEL_MODES * PKPM_DIST_NUM_BLOCKS)
// Largest number of transverse modes of a face (3x).
#define PKPM_REMAP_FLUX_MAX_TRANS_MODES (PKPM_MAX_CONF_BASIS / 2)

// Traces of the frame at one transverse node of a face, and the map of one side to the face frame.
struct pkpm_remap_flux_side {
  double rho, u[3], xi, b[3], theta; // traces at the node
  struct pkpm_harm_face_map map; // the map to the face frame
  double kin[PKPM_RF_SIDE_NUM]; // the part of the map the interval kernels read (PKPM_RF_SIDE_*)
  double fwd[PKPM_FRAME_NUM]; // factors of the change of frame to the face frame (PKPM_FRAME_*)
  double bck[PKPM_FRAME_NUM]; // factors of the change of frame back to the side
};

// The frame and the coefficients of the flux at one transverse node of a face.
struct pkpm_remap_flux_node {
  struct pkpm_remap_flux_side side[2];
  double rho_bar, u_bar[3], xi_bar, b_bar[3]; // the face frame
  double fc[PKPM_FC_NUM]; // face coefficients (pkpm_harm_face_coef)
  double kinks[3]; // sorted face-frame labels where the HLL weights have kinks
  int num_kinks;
};

struct gkyl_dg_pkpm_remap_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int num_conf_basis; // Configuration-space basis size.
  int num_vel_cells; // Cells of the velocity grid.
  int num_surf_basis; // Surface modes per block, direction and side of the kinetic flux.
  int num_trans_modes; // Transverse modes of a face, the surface modes of the fluid flux.
  int num_nodes; // Transverse nodes of a face.
  int num_cell_nodes; // Gauss-Legendre nodes of a cell.
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map).
  pkpm_remap_flux_trace_t trace[3]; // Block traces in each direction.
  pkpm_remap_flux_to_face_t to_face; // Blocks of a velocity cell at a node, in both frames.
  pkpm_remap_flux_interval_t interval[3]; // Supermesh interval in each direction.
  pkpm_remap_flux_from_face_t from_face[3]; // Face-frame projection back to a side, stored.
  pkpm_remap_flux_label_t label[3]; // Label-paired flux of a velocity cell in each direction.
  pkpm_remap_flux_store_t store[3]; // Store of projections in each direction.

  double *psi; // Transverse modes at the nodes, [node][mode].
  double *wpsi; // Node weight times the transverse modes, [node][mode].
  double *conf_trace[3]; // Conf basis at the face nodes per direction, [side][node][basis].
  double *vel_edges; // Edges of the velocity grid.
  double *vel_center; // Center of each velocity cell.
  double *vel_dw; // Width of each velocity cell.
  double *vel_rdw2; // 2/width of each velocity cell.

  // Gauss-Legendre nodes of a cell (those of gkyl_dg_conf_nodal).
  dg_conf_nodal_kern_t conf_to_nodal; // Expansion to node values.
  dg_conf_nodal_kern_t conf_to_modal; // Node values to expansion.
  int *line_upper; // Upper node of the line through a node along a direction, [dir][node].
  double node_to_face; // Distance from the midpoint of a line's nodes to a face over their spacing.

  // Workspace for one face.
  long *column; // Phase index of each velocity cell, [side][cell].
  double *block_traces; // Modal block traces, [side][cell][block][mode].
  double *remap_flux; // Remap flux of the face, [side][cell][block][mode].
  double *label_flux; // Label-paired flux of the face (both sides receive it), [cell][block][mode].
  // Largest speeds of the flux over the nodes of the face and the Gauss-Legendre points of each
  // velocity cell, [cell][fastest signal in either direction, largest spread of the two bounds].
  double *cell_speed;
  struct pkpm_remap_flux_node *nodes; // Frame of every node.
  // Workspace for one transverse node. A side has a slot per velocity cell (slot = cell + 1) and
  // slot 0 for the labels outside its grid, where it has no blocks: an interval of the supermesh
  // then reads and writes the slots of its two cells with no test.
  struct gkyl_supermesh_1d_interval *mesh; // Supermesh.
  double *node_X; // Blocks in the frame of the side, [side][cell][mode][block].
  double *node_Y; // Blocks in the face frame, [side][slot][mode][block] (slot 0: zero).
  double *node_P; // Projection of the face-frame flux, [side][slot][mode][block].
  double *node_cell; // Cell coordinate of each slot (PKPM_RF_CELL_*), [side][slot].
};

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

// The largest share of the flux outside the receiving side's velocity grid.
GKYL_CU_DH static inline double
pkpm_remap_flux_blend(const double *share)
{
  double momentum = 0.0, energy = 0.0;
  for (int s = 0; s < 2; ++s) {
    if (share[PKPM_RF_SHARE_MOM_TOT + s] > 0.0) {
      momentum =
        fmax(momentum, share[PKPM_RF_SHARE_MOM_OUT + s] / share[PKPM_RF_SHARE_MOM_TOT + s]);
    }
    if (share[PKPM_RF_SHARE_EN_TOT + s] > 0.0) {
      energy = fmax(energy, share[PKPM_RF_SHARE_EN_OUT + s] / share[PKPM_RF_SHARE_EN_TOT + s]);
    }
  }
  return fmax(momentum, energy);
}

// Thermal speed at one face of a node line of a cell, from the values of the thermal speed and the
// density at the line's two Gauss-Legendre nodes (lower, upper); sign is +1 for the upper face and
// -1 for the lower, and node_to_face is the distance from the nodes' midpoint to the face over the
// distance between the nodes. tau = 1/xi^2 is extrapolated linearly to the face: exact where the
// pressure is uniform and the density linear. On the face beyond the hotter node the result is
// limited to the ratio of the pressure and density traces (p = rho xi^2 at the nodes, both linear
// along the line), and never below the hot node's own value. A line with a non-positive node value
// keeps its linear extrapolation.
GKYL_CU_DH static inline double
pkpm_remap_flux_xi_line(
  double xi_lo, double xi_up, double rho_lo, double rho_up, double sign, double node_to_face
)
{
  // slope factor of the linear extrapolation from the two nodes to the face
  double slope = sign * node_to_face;
  if (!(xi_lo > 0.0 && xi_up > 0.0)) {
    return 0.5 * (xi_lo + xi_up) + slope * (xi_up - xi_lo);
  }
  double tau_lo = 1.0 / (xi_lo * xi_lo), tau_up = 1.0 / (xi_up * xi_up);
  double tau_face = 0.5 * (tau_lo + tau_up) + slope * (tau_up - tau_lo);
  double sign_hot = tau_lo < tau_up ? -1.0 : 1.0;
  if (sign == sign_hot) {
    double tau_hot = fmin(tau_lo, tau_up);
    double pres_lo = rho_lo / tau_lo, pres_up = rho_up / tau_up;
    double pres_face = 0.5 * (pres_lo + pres_up) + slope * (pres_up - pres_lo);
    double rho_face = 0.5 * (rho_lo + rho_up) + slope * (rho_up - rho_lo);
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
  int num_cell_nodes = up->num_cell_nodes;
  const int *line_upper = &up->line_upper[dir * num_cell_nodes];
  double node_xi[PKPM_MAX_CONF_BASIS], node_rho[PKPM_MAX_CONF_BASIS], line[PKPM_MAX_CONF_BASIS];
  up->conf_to_nodal(xi, node_xi);
  up->conf_to_nodal(&fluid[PKPM_FLUID_RHO * up->num_conf_basis], node_rho);
  for (int q = 0; q < num_cell_nodes; ++q) {
    int q_up = line_upper[q];
    if (q_up != q) {
      line[q] = pkpm_remap_flux_xi_line(
        node_xi[q], node_xi[q_up], node_rho[q], node_rho[q_up], sign, up->node_to_face
      );
      line[q_up] = line[q];
    }
  }
  up->conf_to_modal(line, xi_face);
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

// The face frame of a node and the map of each side to it, and the kinks of the HLL weights. A
// face that is not admissible (a non-positive density or thermal-speed trace at one of its nodes)
// takes the label-paired flux alone: its remap flux carries no weight and is formed with unit
// density and thermal speed on both sides, which keeps it finite.
GKYL_CU_DH static inline void
pkpm_remap_flux_node_frame(struct pkpm_remap_flux_node *nd, int admissible)
{
  double rho[2], xi[2];
  for (int s = 0; s < 2; ++s) {
    rho[s] = admissible ? nd->side[s].rho : 1.0;
    xi[s] = admissible ? nd->side[s].xi : 1.0;
  }
  const struct pkpm_remap_flux_side *sl = &nd->side[0], *sr = &nd->side[1];
  nd->rho_bar = 0.5 * (rho[0] + rho[1]);
  nd->xi_bar = 0.5 * (xi[0] + xi[1]);
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
      rho[s], sd->u, xi[s], nd->rho_bar, nd->u_bar, nd->xi_bar, nd->b_bar, &sd->map
    );
    pkpm_harm_frame_factors(sd->map.K, sd->map.r, sd->map.c_perp, sd->fwd);
    pkpm_harm_frame_factors(sd->map.K_inv, sd->map.r_inv, sd->map.c_perp_inv, sd->bck);
    sd->kin[PKPM_RF_SIDE_R] = sd->map.r;
    sd->kin[PKPM_RF_SIDE_C] = sd->map.c_par;
    sd->kin[PKPM_RF_SIDE_KINV] = sd->map.K_inv;
    sd->kin[PKPM_RF_SIDE_EG] = sd->bck[PKPM_FRAME_K2];
    sd->kin[PKPM_RF_SIDE_CC] = sd->bck[PKPM_FRAME_CC];
    for (int k = 0; k < 3; ++k) {
      sd->kin[PKPM_RF_SIDE_CPERP + k] = sd->map.c_perp_inv[k];
    }
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

// The remap flux of one transverse node of a face: every velocity cell of both sides goes to the
// face frame, the flux is integrated over each interval of the supermesh of the two mapped grids
// and projected on the cells of both sides, and every projection goes back to its side and into
// the remap flux of the face. Returns the face-frame moments of the flux in mom and adds to the
// shares of the face.
GKYL_CU_DH static inline void
pkpm_remap_flux_node_remap(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, int node,
  const struct pkpm_remap_flux_node *nd, int nflux, double *mom, double *share
)
{
  int nv = up->num_vel_cells, na = up->num_trans_modes, nslot = nv + 1;
  int ncoef = PKPM_REMAP_FLUX_NUM_COEFF;
  const double *psi = &up->psi[node * na], *wpsi = &up->wpsi[node * na];
  for (int s = 0; s < 2; ++s) {
    const struct pkpm_remap_flux_side *sd = &nd->side[s];
    for (int k = 0; k < nv; ++k) {
      long slot = s * nslot + k + 1;
      up->to_face(
        psi, sd->fwd, &up->block_traces[(s * nv + k) * nflux], &up->node_X[(s * nv + k) * ncoef],
        &up->node_Y[slot * ncoef]
      );
      // eta = 2 (r wbar + c - w_center)/dw at the face-frame label wbar
      double *cell = &up->node_cell[slot * PKPM_RF_CELL_NUM];
      cell[PKPM_RF_CELL_ETA0] = (sd->map.c_par - up->vel_center[k]) * up->vel_rdw2[k];
      cell[PKPM_RF_CELL_ETA1] = sd->map.r * up->vel_rdw2[k];
      for (int m = 0; m < ncoef; ++m) {
        up->node_P[slot * ncoef + m] = 0.0;
      }
    }
  }

  double tol = 1e-13 * up->phase_grid.dx[up->cdim];
  int num_intervals = gkyl_supermesh_1d_merge(
    nv + 1, up->vel_edges, nd->side[0].map.r, nd->side[0].map.c_par, nv + 1, up->vel_edges,
    nd->side[1].map.r, nd->side[1].map.c_par, nd->num_kinks, nd->kinks, tol, up->mesh
  );
  for (int m = 0; m < PKPM_RF_NUM_MOMS; ++m) {
    mom[m] = 0.0;
  }
  for (int q = 0; q < num_intervals; ++q) {
    const struct gkyl_supermesh_1d_interval *iv = &up->mesh[q];
    long slot_l = iv->k_l + 1, slot_r = nslot + iv->k_r + 1;
    up->interval[dir](
      nd->fc, iv->lo, iv->up, nd->side[0].kin, nd->side[1].kin,
      &up->node_cell[slot_l * PKPM_RF_CELL_NUM], &up->node_cell[slot_r * PKPM_RF_CELL_NUM],
      &up->node_Y[slot_l * ncoef], &up->node_Y[slot_r * ncoef], &up->node_P[slot_l * ncoef],
      &up->node_P[slot_r * ncoef], mom, share
    );
  }

  for (int s = 0; s < 2; ++s) {
    for (int k = 0; k < nv; ++k) {
      up->from_face[dir](
        wpsi, nd->side[s].bck, &up->node_P[(s * nslot + k + 1) * ncoef],
        &up->remap_flux[(s * nv + k) * nflux]
      );
    }
  }
}

// The label-paired flux of one transverse node of a face, into the label-paired flux of the
// face: the flux of the two sides paired at equal w, one velocity cell at a time (the blocks in
// the frame of each side are those pkpm_remap_flux_node_remap left for the node). Returns its
// fluid flux, the moment map of each one-sided contribution in the frame of that side, and raises
// the speeds of the cells of the face to those of this node.
GKYL_CU_DH static inline void
pkpm_remap_flux_node_label(
  const struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, int node,
  const struct pkpm_remap_flux_node *nd, int nflux, double *fluid_flux
)
{
  int nv = up->num_vel_cells, na = up->num_trans_modes, ncoef = PKPM_REMAP_FLUX_NUM_COEFF;
  const double *wpsi = &up->wpsi[node * na];
  double moms[2][PKPM_RF_NUM_MOMS];
  for (int s = 0; s < 2; ++s) {
    for (int m = 0; m < PKPM_RF_NUM_MOMS; ++m) {
      moms[s][m] = 0.0;
    }
  }
  for (int k = 0; k < nv; ++k) {
    double proj[PKPM_REMAP_FLUX_NUM_COEFF];
    up->label[dir](
      nd->fc, up->vel_center[k], up->vel_dw[k], &up->node_X[k * ncoef],
      &up->node_X[(nv + k) * ncoef], proj, moms[0], moms[1], &up->cell_speed[2 * k]
    );
    up->store[dir](wpsi, proj, &up->label_flux[k * nflux]);
  }
  for (int m = 0; m < PKPM_FLUID_NUM_COMP; ++m) {
    fluid_flux[m] = 0.0;
  }
  for (int s = 0; s < 2; ++s) {
    const struct pkpm_remap_flux_side *sd = &nd->side[s];
    double side_flux[PKPM_FLUID_NUM_COMP];
    pkpm_harm_fluid_moment_map(moms[s], sd->rho, sd->u, sd->xi, sd->b, side_flux);
    for (int m = 0; m < PKPM_FLUID_NUM_COMP; ++m) {
      fluid_flux[m] += side_flux[m];
    }
  }
}
