#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_dg_pkpm_remap_flux_surf_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>

gkyl_dg_pkpm_remap_flux_surf *
gkyl_dg_pkpm_remap_flux_surf_inew(const struct gkyl_dg_pkpm_remap_flux_surf_inp *inp)
{
  assert(!inp->use_gpu); // host only: there is no device implementation
  assert(inp->pos_map);
  assert(inp->conf_basis->poly_order == 1);
  assert(
    (inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY ||
     inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR)
  );
  assert(inp->num_face_nodes >= 3);

  struct gkyl_dg_pkpm_remap_flux_surf *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, n_face = inp->num_face_nodes, nc = inp->conf_basis->num_basis;
  assert(inp->phase_basis->ndim == cdim + 1);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->num_conf_basis = nc;
  up->num_vel_cells = inp->phase_grid->cells[cdim];
  up->num_surf_basis = inp->phase_basis->num_basis / 2;
  up->num_trans_modes = nc / 2;
  up->num_nodes = 1;
  for (int d = 1; d < cdim; ++d) {
    up->num_nodes *= n_face;
  }
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  up->to_face = ser_pkpm_remap_flux_to_face_kernels[cdim - 1];
  for (int d = 0; d < cdim; ++d) {
    up->trace[d] = ser_pkpm_remap_flux_trace_kernels[cdim - 1].kernels[d];
    up->interval[d] = pkpm_remap_flux_interval_kernels[d];
    up->from_face[d] = ser_pkpm_remap_flux_from_face_kernels[cdim - 1].kernels[d];
    up->label[d] = pkpm_remap_flux_label_kernels[d];
    up->store[d] = ser_pkpm_remap_flux_store_kernels[cdim - 1].kernels[d];
    assert(up->trace[d] && up->from_face[d] && up->store[d]);
  }

  int nv = up->num_vel_cells, nn = up->num_nodes, na = up->num_trans_modes;
  double wmin = inp->phase_grid->lower[cdim], dw = inp->phase_grid->dx[cdim];
  up->vel_edges = gkyl_malloc((nv + 1) * sizeof(double));
  for (int k = 0; k <= nv; ++k) {
    up->vel_edges[k] = wmin + k * dw;
  }
  up->vel_center = gkyl_malloc(nv * sizeof(double));
  up->vel_dw = gkyl_malloc(nv * sizeof(double));
  up->vel_rdw2 = gkyl_malloc(nv * sizeof(double));
  for (int k = 0; k < nv; ++k) {
    up->vel_center[k] = 0.5 * (up->vel_edges[k + 1] + up->vel_edges[k]);
    up->vel_dw[k] = up->vel_edges[k + 1] - up->vel_edges[k];
    up->vel_rdw2[k] = 2.0 / up->vel_dw[k];
  }

  // The kernels take the velocity modes of the basis to be the p=2 Legendre functions
  // [1, eta, 3 eta^2 - 1] times the factors pkpmRemapFluxVelNorm of their generator: checked
  // here at the ordinates of their 3-point Gauss-Legendre rule.
  struct gkyl_basis vel_basis;
  gkyl_cart_modal_serendip(&vel_basis, 1, 2);
  assert(vel_basis.num_basis == PKPM_NUM_VEL_MODES);
  double eta_edge = 1.0, modes_edge[PKPM_NUM_VEL_MODES];
  vel_basis.eval(&eta_edge, modes_edge);
  double vel_norm[PKPM_NUM_VEL_MODES] = {modes_edge[0], modes_edge[1], 0.5 * modes_edge[2]};
  for (int g = 0; g < 3; ++g) {
    double eta = gkyl_gauss_ordinates[3][g], vel_modes[PKPM_NUM_VEL_MODES];
    vel_basis.eval(&eta, vel_modes);
    double shape[PKPM_NUM_VEL_MODES] = {1.0, eta, 3.0 * eta * eta - 1.0};
    for (int m = 0; m < PKPM_NUM_VEL_MODES; ++m) {
      assert(fabs(vel_modes[m] - vel_norm[m] * shape[m]) < 1e-14);
    }
  }

  // Gauss-Legendre rule of the transverse directions of a face, in row-major order, and the
  // transverse modes (the p=1 basis in those directions) and the configuration basis of each
  // side at its nodes. In 1x the face is a point: one node of unit weight and a unit mode.
  int num_trans = cdim - 1;
  double ordinates1[n_face], weights1[n_face];
  if (n_face <= gkyl_gauss_max) {
    memcpy(ordinates1, gkyl_gauss_ordinates[n_face], sizeof(double[n_face]));
    memcpy(weights1, gkyl_gauss_weights[n_face], sizeof(double[n_face]));
  } else {
    gkyl_gauleg(-1.0, 1.0, ordinates1, weights1, n_face);
  }
  struct gkyl_basis trans_basis;
  if (num_trans > 0) {
    gkyl_cart_modal_serendip(&trans_basis, num_trans, 1);
    assert(trans_basis.num_basis == na);
  }
  up->psi = gkyl_malloc(nn * na * sizeof(double));
  up->wpsi = gkyl_malloc(nn * na * sizeof(double));
  for (int d = 0; d < 3; ++d) {
    up->conf_trace[d] = 0;
  }
  for (int d = 0; d < cdim; ++d) {
    up->conf_trace[d] = gkyl_malloc(2 * nn * nc * sizeof(double));
  }
  int face_lower[GKYL_MAX_CDIM] = {0}, face_upper[GKYL_MAX_CDIM] = {0};
  for (int t = 0; t < num_trans; ++t) {
    face_upper[t] = n_face - 1;
  }
  struct gkyl_range face_range;
  gkyl_range_init(&face_range, num_trans > 0 ? num_trans : 1, face_lower, face_upper);
  assert(face_range.volume == nn);
  struct gkyl_range_iter face_iter;
  gkyl_range_iter_init(&face_iter, &face_range);
  while (gkyl_range_iter_next(&face_iter)) {
    long i = gkyl_range_idx(&face_range, face_iter.idx);
    double z_trans[GKYL_MAX_CDIM], weight = 1.0;
    for (int t = 0; t < num_trans; ++t) {
      z_trans[t] = ordinates1[face_iter.idx[t]];
      weight *= weights1[face_iter.idx[t]];
    }
    if (num_trans > 0) {
      trans_basis.eval(z_trans, &up->psi[i * na]);
    } else {
      up->psi[i * na] = 1.0;
    }
    for (int a = 0; a < na; ++a) {
      up->wpsi[i * na + a] = weight * up->psi[i * na + a];
    }
    for (int d = 0; d < cdim; ++d) {
      for (int s = 0; s < 2; ++s) {
        double z[GKYL_MAX_CDIM];
        z[d] = s == 0 ? 1.0 : -1.0; // the left cell at its upper edge, the right cell at its lower
        for (int t = 0, dt = 0; dt < cdim; ++dt) {
          if (dt != d) {
            z[dt] = z_trans[t++];
          }
        }
        inp->conf_basis->eval(z, &up->conf_trace[d][(s * nn + i) * nc]);
      }
    }
  }

  // The Gauss-Legendre nodes of a cell (those of gkyl_dg_conf_nodal, row-major): the conversions
  // to and from the node values, for every direction the upper node of the line through each node
  // along it, and the factor of the linear extrapolation from the two nodes of a line to a face.
  up->conf_to_nodal = dg_conf_nodal_conf_to_nodal_kernels.kernels[cdim - 1];
  up->conf_to_modal = dg_conf_nodal_conf_to_modal_kernels.kernels[cdim - 1];
  up->node_to_face = 0.5 / gkyl_gauss_ordinates[2][1];
  int node_lower[GKYL_MAX_CDIM], node_upper[GKYL_MAX_CDIM];
  for (int d = 0; d < cdim; ++d) {
    node_lower[d] = 0;
    node_upper[d] = 1;
  }
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, cdim, node_lower, node_upper);
  up->num_cell_nodes = node_range.volume;
  assert(up->num_cell_nodes == nc);
  up->line_upper = gkyl_malloc(cdim * up->num_cell_nodes * sizeof(int));
  struct gkyl_range_iter node_iter;
  gkyl_range_iter_init(&node_iter, &node_range);
  while (gkyl_range_iter_next(&node_iter)) {
    long q = gkyl_range_idx(&node_range, node_iter.idx);
    for (int d = 0; d < cdim; ++d) {
      int idx_up[GKYL_MAX_CDIM];
      gkyl_copy_int_arr(cdim, node_iter.idx, idx_up);
      idx_up[d] = 1;
      up->line_upper[d * up->num_cell_nodes + q] = gkyl_range_idx(&node_range, idx_up);
    }
  }

  int nflux = PKPM_DIST_NUM_BLOCKS * up->num_surf_basis, nslot = nv + 1;
  int ncoef = PKPM_REMAP_FLUX_NUM_COEFF;
  up->column = gkyl_malloc(2 * nv * sizeof(long));
  up->block_traces = gkyl_malloc(2 * nv * nflux * sizeof(double));
  up->remap_flux = gkyl_malloc(2 * nv * nflux * sizeof(double));
  up->label_flux = gkyl_malloc(nv * nflux * sizeof(double));
  up->cell_speed = gkyl_malloc(nv * 2 * sizeof(double));
  up->nodes = gkyl_malloc(nn * sizeof(struct pkpm_remap_flux_node));
  up->mesh = gkyl_malloc(
    gkyl_supermesh_1d_max_intervals(nv + 1, nv + 1, 3) * sizeof(struct gkyl_supermesh_1d_interval)
  );
  up->node_X = gkyl_malloc(2 * nv * ncoef * sizeof(double));
  up->node_Y = gkyl_malloc(2 * nslot * ncoef * sizeof(double));
  up->node_P = gkyl_malloc(2 * nslot * ncoef * sizeof(double));
  up->node_cell = gkyl_malloc(2 * nslot * PKPM_RF_CELL_NUM * sizeof(double));
  // Slot 0 of each side, the labels outside its grid: no blocks, a projection that receives
  // zero, and every share of the flux counted as outside. The slots of the cells are inside.
  for (int s = 0; s < 2; ++s) {
    for (int m = 0; m < ncoef; ++m) {
      up->node_Y[s * nslot * ncoef + m] = 0.0;
      up->node_P[s * nslot * ncoef + m] = 0.0;
    }
    for (int k = 0; k < nslot; ++k) {
      double *cell = &up->node_cell[(s * nslot + k) * PKPM_RF_CELL_NUM];
      cell[PKPM_RF_CELL_ETA0] = 0.0;
      cell[PKPM_RF_CELL_ETA1] = 0.0;
      cell[PKPM_RF_CELL_OUT] = k == 0 ? 1.0 : 0.0;
    }
  }

  return up;
}

// The face between the configuration cells idx_l (lower) and idx_r (upper) in direction dir. The
// fluxes, the blend and the rates of the face are stored on idx_r. The rates take the smaller of
// the two cells' position-map Jacobians. Every face runs the same operations: the remap flux and
// the label-paired flux of every node, combined with the blend weight of the face.
static void
remap_flux_face(
  struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, const int *idx_l, const int *idx_r,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar,
  const struct gkyl_array *fluid, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *fluid_flux_surf, struct gkyl_array *rate_surf_f,
  struct gkyl_array *rate_surf_fluid, struct gkyl_array *blend
)
{
  int cdim = up->cdim, nv = up->num_vel_cells, nk = up->num_surf_basis;
  int nflux = PKPM_DIST_NUM_BLOCKS * nk, na = up->num_trans_modes, nn = up->num_nodes;
  int num_fluid = PKPM_FLUID_NUM_COMP * na;
  const int *idx[2] = {idx_l, idx_r};
  long cloc[2];
  double jdinv[2];
  const double *bvar_s[2], *fluid_s[2], *u_s[2], *xi_s[2], *theta_s[2];
  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    rem_dir[d] = 1;
  }
  for (int s = 0; s < 2; ++s) {
    cloc[s] = gkyl_range_idx(conf_range_ext, idx[s]);
    jdinv[s] = 1.0 / ((const double *)gkyl_array_cfetch(up->jacob_pos, cloc[s]))[2 * dir];
    bvar_s[s] = gkyl_array_cfetch(bvar, cloc[s]);
    fluid_s[s] = gkyl_array_cfetch(fluid, cloc[s]);
    u_s[s] = gkyl_array_cfetch(u, cloc[s]);
    xi_s[s] = gkyl_array_cfetch(xi, cloc[s]);
    theta_s[s] = gkyl_array_cfetch(theta, cloc[s]);

    int col_idx[GKYL_MAX_DIM] = {0};
    for (int d = 0; d < cdim; ++d) {
      col_idx[d] = idx[s][d];
    }
    struct gkyl_range column;
    gkyl_range_deflate(&column, phase_range_ext, rem_dir, col_idx);
    for (int k = 0; k < nv; ++k) {
      int jv = phase_range->lower[cdim] + k;
      up->column[s * nv + k] = gkyl_range_idx(&column, &jv);
    }
  }

  // 1. The modal block traces of both velocity columns, once per face.
  for (int k = 0; k < nv; ++k) {
    up->trace[dir](
      jdinv[0], jdinv[1], gkyl_array_cfetch(fin, up->column[k]),
      gkyl_array_cfetch(fin, up->column[nv + k]), &up->block_traces[k * nflux],
      &up->block_traces[(nv + k) * nflux]
    );
  }

  // 2. The traces of the frame and the face coefficients of every node. The thermal speed of each
  // side on this face comes from 1/xi^2 along the node lines of its cell.
  double xi_face[2][PKPM_MAX_CONF_BASIS];
  const double *xi_face_s[2] = {xi_face[0], xi_face[1]};
  for (int s = 0; s < 2; ++s) {
    pkpm_remap_flux_xi_face(up, dir, s == 0 ? 1.0 : -1.0, fluid_s[s], xi_s[s], xi_face[s]);
  }
  int admissible = 1;
  for (int i = 0; i < nn; ++i) {
    admissible = pkpm_remap_flux_node_traces(
                   up, dir, i, bvar_s, fluid_s, u_s, xi_face_s, theta_s, &up->nodes[i]
                 ) &&
                 admissible;
  }

  // 3. The remap flux and the label-paired flux of every node, kinetic and fluid, the shares of
  // the remap flux summed over the nodes, and the largest speeds of the flux in every velocity
  // cell.
  double fluid_remap[PKPM_FLUID_NUM_COMP * PKPM_REMAP_FLUX_MAX_TRANS_MODES] = {0.0};
  double fluid_label[PKPM_FLUID_NUM_COMP * PKPM_REMAP_FLUX_MAX_TRANS_MODES] = {0.0};
  double share[PKPM_RF_SHARE_NUM] = {0.0};
  for (int m = 0; m < 2 * nv * nflux; ++m) {
    up->remap_flux[m] = 0.0;
  }
  for (int m = 0; m < nv * nflux; ++m) {
    up->label_flux[m] = 0.0;
  }
  for (int m = 0; m < 2 * nv; ++m) {
    up->cell_speed[m] = 0.0;
  }
  for (int i = 0; i < nn; ++i) {
    struct pkpm_remap_flux_node *nd = &up->nodes[i];
    const double *wpsi = &up->wpsi[i * na];
    double mom[PKPM_RF_NUM_MOMS], node_flux[PKPM_FLUID_NUM_COMP];
    pkpm_remap_flux_node_frame(nd, admissible);
    pkpm_remap_flux_node_remap(up, dir, i, nd, nflux, mom, share);
    pkpm_harm_fluid_moment_map(mom, nd->rho_bar, nd->u_bar, nd->xi_bar, nd->b_bar, node_flux);
    for (int m = 0; m < PKPM_FLUID_NUM_COMP; ++m) {
      for (int a = 0; a < na; ++a) {
        fluid_remap[m * na + a] += wpsi[a] * node_flux[m];
      }
    }
    pkpm_remap_flux_node_label(up, dir, i, nd, nflux, node_flux);
    for (int m = 0; m < PKPM_FLUID_NUM_COMP; ++m) {
      for (int a = 0; a < na; ++a) {
        fluid_label[m * na + a] += wpsi[a] * node_flux[m];
      }
    }
  }

  // 4. The rates of the face, from the largest speeds of every velocity cell and the narrower of
  // the two cells.
  double rate = 1.5 * (2.0 / up->phase_grid.dx[dir]) * fmax(jdinv[0], jdinv[1]);
  double alpha_fluid = 0.0;
  for (int k = 0; k < nv; ++k) {
    ((double *)gkyl_array_fetch(rate_surf_f, up->column[nv + k]))[dir] =
      rate * up->cell_speed[2 * k];
    alpha_fluid = fmax(alpha_fluid, up->cell_speed[2 * k + 1]);
  }
  ((double *)gkyl_array_fetch(rate_surf_fluid, cloc[1]))[dir] = rate * alpha_fluid;

  // 5. The blend: the largest share of the remap flux outside the receiving grid, and 1 on a face
  // that is not admissible. The flux is (1 - blend) remap + blend label-paired; both sides'
  // kinetic fluxes go on the owning cell of the face.
  double face_blend = admissible ? pkpm_remap_flux_blend(share) : 1.0;
  double remap_weight = 1.0 - face_blend;
  double *fluid_flux = &((double *)gkyl_array_fetch(fluid_flux_surf, cloc[1]))[dir * num_fluid];
  for (int m = 0; m < num_fluid; ++m) {
    fluid_flux[m] = remap_weight * fluid_remap[m] + face_blend * fluid_label[m];
  }
  for (int k = 0; k < nv; ++k) {
    double *out =
      &((double *)gkyl_array_fetch(conf_flux_surf, up->column[nv + k]))[2 * dir * nflux];
    const double *label = &up->label_flux[k * nflux];
    for (int s = 0; s < 2; ++s) {
      const double *remap = &up->remap_flux[(s * nv + k) * nflux];
      for (int m = 0; m < nflux; ++m) {
        out[s * nflux + m] = remap_weight * remap[m] + face_blend * label[m];
      }
    }
  }
  ((double *)gkyl_array_fetch(blend, cloc[1]))[dir] = face_blend;
}

void
gkyl_dg_pkpm_remap_flux_surf_advance(
  struct gkyl_dg_pkpm_remap_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar,
  const struct gkyl_array *fluid, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *fluid_flux_surf, struct gkyl_array *rate_surf_f,
  struct gkyl_array *rate_surf_fluid, struct gkyl_array *cflrate_f,
  struct gkyl_array *cflrate_fluid, struct gkyl_array *blend
)
{
  int cdim = up->cdim, pdim = cdim + 1;
  int idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM], idx_up[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    for (int dir = 0; dir < cdim; ++dir) {
      // Lower face: the lower neighbor (interior or ghost) is the l side.
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      remap_flux_face(
        up, dir, idx_l, iter.idx, conf_range_ext, phase_range, phase_range_ext, bvar, fluid, u, xi,
        theta, fin, conf_flux_surf, fluid_flux_surf, rate_surf_f, rate_surf_fluid, blend
      );
      // Upper end of the range: the face owned by the upper neighbor, with the current cell as
      // the l side.
      if (iter.idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        remap_flux_face(
          up, dir, iter.idx, idx_r, conf_range_ext, phase_range, phase_range_ext, bvar, fluid, u,
          xi, theta, fin, conf_flux_surf, fluid_flux_surf, rate_surf_f, rate_surf_fluid, blend
        );
      }
    }
  }

  // A cell takes, in each direction, the larger of the rates of its two faces (its own lower
  // face and that of its upper neighbor); the rates of the directions add.
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(phase_range_ext, iter.idx);
    const double *rate_lo = gkyl_array_cfetch(rate_surf_f, loc);
    double *cflrate_d = gkyl_array_fetch(cflrate_f, loc);
    for (int dir = 0; dir < cdim; ++dir) {
      gkyl_copy_int_arr(pdim, iter.idx, idx_up);
      idx_up[dir] = idx_up[dir] + 1;
      const double *rate_up =
        gkyl_array_cfetch(rate_surf_f, gkyl_range_idx(phase_range_ext, idx_up));
      cflrate_d[0] += fmax(rate_lo[dir], rate_up[dir]);
    }
  }
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range_ext, iter.idx);
    const double *rate_lo = gkyl_array_cfetch(rate_surf_fluid, loc);
    double *cflrate_d = gkyl_array_fetch(cflrate_fluid, loc);
    for (int dir = 0; dir < cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idx_up);
      idx_up[dir] = idx_up[dir] + 1;
      const double *rate_up =
        gkyl_array_cfetch(rate_surf_fluid, gkyl_range_idx(conf_range_ext, idx_up));
      cflrate_d[0] += fmax(rate_lo[dir], rate_up[dir]);
    }
  }
}

int
gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_remap_flux_surf *up)
{
  return up->num_surf_basis;
}

int
gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(const struct gkyl_dg_pkpm_remap_flux_surf *up)
{
  return up->num_trans_modes;
}

void
gkyl_dg_pkpm_remap_flux_surf_release(struct gkyl_dg_pkpm_remap_flux_surf *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  gkyl_free(up->vel_edges);
  gkyl_free(up->vel_center);
  gkyl_free(up->vel_dw);
  gkyl_free(up->vel_rdw2);
  gkyl_free(up->psi);
  gkyl_free(up->wpsi);
  for (int d = 0; d < up->cdim; ++d) {
    gkyl_free(up->conf_trace[d]);
  }
  gkyl_free(up->line_upper);
  gkyl_free(up->column);
  gkyl_free(up->block_traces);
  gkyl_free(up->remap_flux);
  gkyl_free(up->label_flux);
  gkyl_free(up->cell_speed);
  gkyl_free(up->nodes);
  gkyl_free(up->mesh);
  gkyl_free(up->node_X);
  gkyl_free(up->node_Y);
  gkyl_free(up->node_P);
  gkyl_free(up->node_cell);
  gkyl_free(up);
}
