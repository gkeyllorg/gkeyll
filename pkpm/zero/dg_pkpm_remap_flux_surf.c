#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_dg_pkpm_remap_flux_surf_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>

gkyl_dg_pkpm_remap_flux_surf *
gkyl_dg_pkpm_remap_flux_surf_inew(const struct gkyl_dg_pkpm_remap_flux_surf_inp *inp)
{
  assert(!inp->use_gpu); // the GPU port follows the CPU verification
  assert(inp->pos_map);
  assert(inp->conf_basis->poly_order == 1);
  assert(inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
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

  for (int d = 0; d < cdim; ++d) {
    up->trace[d] = ser_pkpm_remap_flux_trace_kernels[cdim - 1].kernels[d];
    up->store[d] = ser_pkpm_remap_flux_store_kernels[cdim - 1].kernels[d];
    assert(up->trace[d] && up->store[d]);
  }

  int nv = up->num_vel_cells, nn = up->num_nodes, na = up->num_trans_modes;
  double wmin = inp->phase_grid->lower[cdim], dw = inp->phase_grid->dx[cdim];
  up->vel_edges = gkyl_malloc((nv + 1) * sizeof(double));
  for (int k = 0; k <= nv; ++k) {
    up->vel_edges[k] = wmin + k * dw;
  }

  // Gauss-Legendre rule of the transverse directions of a face (1x: one node with unit weight)
  double xg[n_face], wg[n_face];
  gkyl_gauleg(-1.0, 1.0, xg, wg, n_face);
  int num_trans = cdim - 1;
  up->psi = gkyl_malloc(nn * na * sizeof(double));
  up->wpsi = gkyl_malloc(nn * na * sizeof(double));
  for (int d = 0; d < 3; ++d) {
    up->conf_trace[d] = 0;
  }
  for (int d = 0; d < cdim; ++d) {
    up->conf_trace[d] = gkyl_malloc(2 * nn * nc * sizeof(double));
  }
  for (int i = 0; i < nn; ++i) {
    // the transverse Gauss indices of node i, the first transverse direction varying slowest
    int ig[2] = {0, 0};
    if (num_trans == 1) {
      ig[0] = i;
    } else if (num_trans == 2) {
      ig[0] = i / n_face;
      ig[1] = i % n_face;
    }
    double weight = 1.0;
    for (int t = 0; t < num_trans; ++t) {
      weight *= wg[ig[t]];
    }
    for (int a = 0; a < na; ++a) {
      double val = 1.0;
      for (int t = 0; t < num_trans; ++t) {
        val *= (a >> t) & 1 ? 1.224744871391589 * xg[ig[t]] : 0.7071067811865475;
      }
      up->psi[i * na + a] = val;
      up->wpsi[i * na + a] = weight * val;
    }
    for (int d = 0; d < cdim; ++d) {
      for (int s = 0; s < 2; ++s) {
        double z[GKYL_MAX_CDIM];
        z[d] = s == 0 ? 1.0 : -1.0; // the left cell at its upper edge, the right cell at its lower
        for (int t = 0, dt = 0; dt < cdim; ++dt) {
          if (dt != d) {
            z[dt] = xg[ig[t++]];
          }
        }
        inp->conf_basis->eval(z, &up->conf_trace[d][(s * nn + i) * nc]);
      }
    }
  }

  int nflux = 8 * up->num_surf_basis;
  up->column = gkyl_malloc(2 * nv * sizeof(long));
  up->block_traces = gkyl_malloc(2 * nv * nflux * sizeof(double));
  up->face_flux = gkyl_malloc(2 * nv * nflux * sizeof(double));
  up->mesh = gkyl_malloc(
    gkyl_supermesh_1d_max_intervals(nv + 1, nv + 1, 3) * sizeof(struct gkyl_supermesh_1d_interval)
  );
  up->nodes = gkyl_malloc(nn * sizeof(struct pkpm_remap_flux_node));

  return up;
}

// The face between the configuration cells idx_l (lower) and idx_r (upper) in direction dir. The
// fluxes are stored on idx_r; the rates go to idx_r if owner_r and to idx_l otherwise.
static void
remap_flux_face(
  struct gkyl_dg_pkpm_remap_flux_surf *up, int dir, const int *idx_l, const int *idx_r, int owner_r,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar,
  const struct gkyl_array *fluid, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *fluid_flux_surf, struct gkyl_array *cflrate_f,
  struct gkyl_array *cflrate_fluid, struct gkyl_array *blend
)
{
  int cdim = up->cdim, nv = up->num_vel_cells, nk = up->num_surf_basis, nflux = 8 * nk;
  int na = up->num_trans_modes, nn = up->num_nodes;
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

  // 2. The traces of the frame and the face coefficients of every node.
  int admissible = 1;
  for (int i = 0; i < nn; ++i) {
    admissible =
      pkpm_remap_flux_node_traces(up, dir, i, bvar_s, fluid_s, u_s, xi_s, theta_s, &up->nodes[i]) &&
      admissible;
  }

  // 3. The rates of the owner cell, from the speed bounds of every node and velocity node.
  int owner = owner_r ? 1 : 0;
  double rate = 1.5 * (2.0 / up->phase_grid.dx[dir]) * jdinv[owner];
  double alpha_fluid = 0.0;
  for (int k = 0; k < nv; ++k) {
    double dw = up->vel_edges[k + 1] - up->vel_edges[k];
    double w_center = 0.5 * (up->vel_edges[k + 1] + up->vel_edges[k]);
    double alpha_kin = 0.0;
    for (int i = 0; i < nn; ++i) {
      for (int j = 0; j < 3; ++j) {
        double bounds[2];
        pkpm_harm_conf_bounds(
          w_center + 0.5 * pkpm_remap_flux_gn3[j] * dw, up->nodes[i].fc, 1, bounds
        );
        alpha_kin = fmax(alpha_kin, fmax(-bounds[0], bounds[1]));
        alpha_fluid = fmax(alpha_fluid, bounds[1] - bounds[0]);
      }
    }
    ((double *)gkyl_array_fetch(cflrate_f, up->column[owner * nv + k]))[0] += rate * alpha_kin;
  }
  ((double *)gkyl_array_fetch(cflrate_fluid, cloc[owner]))[0] += rate * alpha_fluid;

  // 4. The remap flux of every node; the blend from the shares summed over the nodes.
  double *fluid_flux = &((double *)gkyl_array_fetch(fluid_flux_surf, cloc[1]))[dir * 5 * na];
  for (int m = 0; m < 5 * na; ++m) {
    fluid_flux[m] = 0.0;
  }
  for (int m = 0; m < 2 * nv * nflux; ++m) {
    up->face_flux[m] = 0.0;
  }
  double face_blend = 1.0;
  if (admissible) {
    struct pkpm_remap_flux_share share = {{0.0}};
    for (int i = 0; i < nn; ++i) {
      struct pkpm_remap_flux_node *nd = &up->nodes[i];
      pkpm_remap_flux_node_frame(nd);
      double mom[7], node_flux[5];
      pkpm_remap_flux_node_sweep(up, dir, i, nd, up->face_flux, nflux, mom, &share);
      pkpm_harm_fluid_moment_map(mom, nd->rho_bar, nd->u_bar, nd->xi_bar, nd->b_bar, node_flux);
      for (int m = 0; m < 5; ++m) {
        for (int a = 0; a < na; ++a) {
          fluid_flux[m * na + a] += up->wpsi[i * na + a] * node_flux[m];
        }
      }
    }
    face_blend = pkpm_remap_flux_blend(&share);
  }

  // 5. The blend: scale the remap flux and add the label-paired flux of the same nodes.
  if (face_blend > 0.0) {
    double remap_weight = 1.0 - face_blend;
    for (int m = 0; m < 2 * nv * nflux; ++m) {
      up->face_flux[m] *= remap_weight;
    }
    for (int m = 0; m < 5 * na; ++m) {
      fluid_flux[m] *= remap_weight;
    }
    for (int i = 0; i < nn; ++i) {
      double node_flux[5];
      pkpm_remap_flux_node_label(
        up, dir, i, face_blend, &up->nodes[i], up->face_flux, nflux, node_flux
      );
      for (int m = 0; m < 5; ++m) {
        for (int a = 0; a < na; ++a) {
          fluid_flux[m * na + a] += face_blend * up->wpsi[i * na + a] * node_flux[m];
        }
      }
    }
  }

  // 6. Both sides' kinetic fluxes on the owning cell of the face.
  for (int k = 0; k < nv; ++k) {
    double *out =
      &((double *)gkyl_array_fetch(conf_flux_surf, up->column[nv + k]))[2 * dir * nflux];
    for (int m = 0; m < nflux; ++m) {
      out[m] = up->face_flux[k * nflux + m];
      out[nflux + m] = up->face_flux[(nv + k) * nflux + m];
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
  struct gkyl_array *fluid_flux_surf, struct gkyl_array *cflrate_f,
  struct gkyl_array *cflrate_fluid, struct gkyl_array *blend
)
{
  int cdim = up->cdim;
  int idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    for (int dir = 0; dir < cdim; ++dir) {
      // Lower face: the left neighbor (interior or ghost) is the l side.
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      remap_flux_face(
        up, dir, idx_l, iter.idx, 1, conf_range_ext, phase_range, phase_range_ext, bvar, fluid, u,
        xi, theta, fin, conf_flux_surf, fluid_flux_surf, cflrate_f, cflrate_fluid, blend
      );
      // Upper domain boundary: the face owned by the upper ghost cell, with the current cell as
      // the l side; its speeds enter the rate of the current cell.
      if (iter.idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        remap_flux_face(
          up, dir, iter.idx, idx_r, 0, conf_range_ext, phase_range, phase_range_ext, bvar, fluid, u,
          xi, theta, fin, conf_flux_surf, fluid_flux_surf, cflrate_f, cflrate_fluid, blend
        );
      }
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
  gkyl_free(up->psi);
  gkyl_free(up->wpsi);
  for (int d = 0; d < up->cdim; ++d) {
    gkyl_free(up->conf_trace[d]);
  }
  gkyl_free(up->column);
  gkyl_free(up->block_traces);
  gkyl_free(up->face_flux);
  gkyl_free(up->mesh);
  gkyl_free(up->nodes);
  gkyl_free(up);
}
