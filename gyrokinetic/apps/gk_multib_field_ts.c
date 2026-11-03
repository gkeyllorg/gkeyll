#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_conn.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_translate_dim.h>
#include <gkyl_fem_proj.h>
#include <gkyl_array_reduce.h>

// Each rank holds its block-global boundary planes. Radial gathers replicate the
// connected plane, just as the perpendicular Poisson solve replicates its RHS.
struct gk_multib_field_ts_block {
  struct gkyl_basis basis;
  struct gkyl_basis *basis_on_dev;
  struct gkyl_range local, ext, connected, connected_ext, block_in_connected;
  struct gkyl_range packed_ghost;
  struct gkyl_range solve[2], neighbor[2][2], boundary_skin[2][2];
  int ghost_owner[2];
  bool shifted[2], check_boundary[2][2];
  double boundary_value[2][2];
  double *boundary_minmax;
  struct gkyl_array *ghost_buffer, *bc_buffer, *weight;
  struct gkyl_array *boundary[2], *boundary_buffer[2][2];
  struct gkyl_bc_basic_gyrokinetic *physical[2], *copy_boundary[2][2];
  struct gkyl_deflate_zsurf *deflate[2];
  struct gkyl_translate_dim *inflate;
  struct gkyl_fem_proj *projector[2];
  struct gkyl_fem_parproj *smooth, *smooth_rho;
};

struct gk_multib_field_ts {
  int num_blocks;
  bool project, use_gpu;
  struct gk_multib_field_ts_block *blocks;
  struct gkyl_multib_comm_conn **send, **recv;
  struct gkyl_array **local[2], **connected[2], **fem[2];
};

static bool
gk_multib_field_ts_shifted(const struct gkyl_gyrokinetic_multib_app *app, int bid, int edge)
{
  return app->block_topo->conn[bid].connections[2][edge].edge != GKYL_PHYSICAL;
}

static double
gk_multib_field_ts_radial_bc(const struct gkyl_gyrokinetic_multib *inp, int bid, int edge)
{
  for (int i = 0; i < inp->field.num_physical_bcs; ++i) {
    const struct gkyl_gyrokinetic_bc *bc = &inp->field.bcs[i];
    if (bc->bidx == bid && bc->dir == 0 && bc->edge == edge) {
      if (bc->type != GKYL_BC_GK_FIELD_DIRICHLET) {
        gkyl_exit("Multiblock TS projection requires constant radial Dirichlet boundaries.");
      }
      return bc->value[0];
    }
  }
  gkyl_exit("Multiblock TS projection is missing a physical radial boundary condition.");
  return 0.0;
}

// Only real physical radial faces and user biases constrain the parallel solve.
// The child app's artificial Dirichlet values at block interfaces must not enter.
static void
gk_multib_field_ts_new_smoother(
  const struct gkyl_gyrokinetic_multib *inp, struct gkyl_gyrokinetic_app *sbapp, int bid,
  struct gk_multib_field_ts_block *block
)
{
  int nuser = inp->field.bias_line_list ? inp->field.bias_line_list->num_bias_line : 0;
  struct gkyl_poisson_bias_line *lines =
    gkyl_malloc(GKYL_MAX2(1, 2 * inp->field.num_physical_bcs + nuser) * sizeof(*lines));
  int count = 0;
  for (int i = 0; i < inp->field.num_physical_bcs; ++i) {
    const struct gkyl_gyrokinetic_bc *bc = &inp->field.bcs[i];
    if (bc->bidx != bid || bc->dir != 0 || bc->type != GKYL_BC_GK_FIELD_DIRICHLET) {
      continue;
    }
    for (int edge = 0; edge < 2; ++edge) {
      lines[count++] = (struct gkyl_poisson_bias_line){
        .perp_dirs = {0, 2},
        .perp_coords =
          {bc->edge == GKYL_LOWER_EDGE ? sbapp->grid.lower[0] : sbapp->grid.upper[0],
           edge == 0 ? sbapp->grid.lower[2] : sbapp->grid.upper[2]},
        .val = bc->value[0],
      };
    }
  }
  for (int i = 0; i < nuser; ++i) {
    lines[count++] = inp->field.bias_line_list->bl[i];
  }
  struct gkyl_poisson_bias_line_list biases = {.num_bias_line = count, .bl = lines};
  block->smooth = gkyl_fem_parproj_new(
    &sbapp->global, &sbapp->grid, &sbapp->basis, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, &biases, NULL,
    NULL, sbapp->use_gpu
  );
  block->smooth_rho = gkyl_fem_parproj_new(
    &sbapp->global, &sbapp->grid, &sbapp->basis, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, NULL, NULL, NULL,
    sbapp->use_gpu
  );
  gkyl_free(lines);
}

static void
gk_multib_field_ts_new_projection(
  const struct gkyl_gyrokinetic_multib *inp, struct gkyl_gyrokinetic_multib_app *app,
  struct gk_multib_field_ts *ts, int bI, const int *ids, int count, int own, const int *offsets
)
{
  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[bI];
  struct gk_multib_field_ts_block *block = &ts->blocks[bI];
  // Unit index-space grid: cell measures in weight recover the computational
  // dx dy inner product even when radial blocks have different cell widths.
  struct gkyl_rect_grid grid;
  int cells[] = {offsets[count], sbapp->grid.cells[1]};
  gkyl_rect_grid_init(&grid, 2, (double[]){0, 0}, (double[]){cells[0], cells[1]}, cells);
  gkyl_create_grid_ranges(&grid, (int[]){1, 1}, &block->connected_ext, &block->connected);
  gkyl_sub_range_init(
    &block->block_in_connected, &block->connected_ext, (int[]){offsets[own] + 1, 1},
    (int[]){offsets[own + 1], cells[1]}
  );
  block->weight = mkarr(app->use_gpu, 4, block->connected_ext.volume);
  for (int i = 0; i < count; ++i) {
    const struct gkyl_gk_block_geom_info *geom =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, ids[i]);
    if (geom->cells[1] != cells[1] ||
        fabs(geom->lower[1] - sbapp->grid.lower[1]) > 1e-12 * sbapp->grid.dx[1] ||
        fabs(geom->upper[1] - sbapp->grid.upper[1]) > 1e-12 * sbapp->grid.dx[1]) {
      gkyl_exit("Multiblock TS projection requires matching y grids across radial blocks.");
    }
    struct gkyl_range range;
    gkyl_sub_range_init(
      &range, &block->connected_ext, (int[]){offsets[i] + 1, 1}, (int[]){offsets[i + 1], cells[1]}
    );
    double dx = (geom->upper[0] - geom->lower[0]) / geom->cells[0];
    if (dx <= 0) {
      gkyl_exit("Multiblock TS projection requires positive radial cell widths.");
    }
    gkyl_array_shiftc_range(block->weight, 2 * dx / sbapp->grid.dx[0], 0, &range);
  }
  for (int edge = 0; edge < 2; ++edge) {
    ts->local[edge][bI] = mkarr(app->use_gpu, 4, block->ext.volume);
    ts->connected[edge][bI] = mkarr(app->use_gpu, 4, block->connected_ext.volume);
    ts->fem[edge][bI] = mkarr(app->use_gpu, 4, block->connected_ext.volume);
    block->deflate[edge] =
      gkyl_deflate_zsurf_new(&sbapp->basis, &block->basis, 1 - edge, app->use_gpu);
    if (!block->shifted[edge]) {
      continue;
    }
    int first = own, last = own;
    while (first > 0 && gk_multib_field_ts_shifted(app, ids[first - 1], edge)) {
      --first;
    }
    while (last + 1 < count && gk_multib_field_ts_shifted(app, ids[last + 1], edge)) {
      ++last;
    }
    gkyl_sub_range_init(
      &block->solve[edge], &block->connected_ext, (int[]){offsets[first] + 1, 1},
      (int[]){offsets[last + 1], cells[1]}
    );
    struct gkyl_fem_proj_bc bcs = {
      .lo_type = {GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_PERIODIC},
      .up_type = {GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_PERIODIC},
    };
    block->boundary[edge] = mkarr(app->use_gpu, 4, block->connected_ext.volume);
    for (int side = 0; side < 2; ++side) {
      int end = side ? last : first;
      enum gkyl_fem_proj_bc_type *type = side ? &bcs.up_type[0] : &bcs.lo_type[0];
      double *value = side ? &bcs.up_value[0] : &bcs.lo_value[0];
      if (side ? last == count - 1 : first == 0) {
        *value = gk_multib_field_ts_radial_bc(inp, ids[end], side);
      } else {
        *type = GKYL_FEM_PROJ_DIRICHLET_VARYING;
        int skin = side ? offsets[last + 1] : offsets[first] + 1;
        int neighbor = skin + (side ? 1 : -1);
        gkyl_sub_range_init(
          &block->boundary_skin[edge][side], &block->connected_ext, (int[]){skin, 1},
          (int[]){skin, cells[1]}
        );
        gkyl_sub_range_init(
          &block->neighbor[edge][side], &block->connected_ext, (int[]){neighbor, 1},
          (int[]){neighbor, cells[1]}
        );
        block->boundary_buffer[edge][side] = mkarr(app->use_gpu, 4, cells[1]);
        block->copy_boundary[edge][side] = gkyl_bc_basic_gyrokinetic_new(
          0, side ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
          block->basis_on_dev, &block->neighbor[edge][side], &block->boundary_skin[edge][side], 4,
          2, app->use_gpu
        );
      }
    }
    int nuser = inp->field.bias_line_list ? inp->field.bias_line_list->num_bias_line : 0;
    struct gkyl_fem_proj_bias_plane *planes = gkyl_malloc((nuser + 1) * sizeof(*planes));
    int nplanes = 0;
    for (int i = 0; i < nuser; ++i) {
      const struct gkyl_poisson_bias_line *line = &inp->field.bias_line_list->bl[i];
      if (line->perp_dirs[0] != 0 || line->perp_dirs[1] != 2) {
        gkyl_exit("Multiblock TS projection supports bias lines parallel to y only.");
      }
      double z = edge ? sbapp->grid.upper[2] : sbapp->grid.lower[2];
      if (fabs(line->perp_coords[1] - z) > 1e-9 * sbapp->grid.dx[2]) {
        continue;
      }
      for (int j = first; j <= last; ++j) {
        const struct gkyl_gk_block_geom_info *geom =
          gkyl_gk_block_geom_get_block(app->gk_block_geom, ids[j]);
        double position = (line->perp_coords[0] - geom->lower[0]) * geom->cells[0] /
                          (geom->upper[0] - geom->lower[0]);
        if (position < -1e-9 || position > geom->cells[0] + 1e-9) {
          continue;
        }
        if (!isfinite(position) || !isfinite(line->val) ||
            fabs(position - round(position)) > 1e-9) {
          gkyl_exit("Multiblock TS projection requires finite, grid-aligned bias lines.");
        }
        int face = offsets[j] + (int)round(position);
        // An interface bias must agree with the unmodified neighbor. Keep the
        // varying constraint and check that agreement when preparing the RHS.
        if ((face == offsets[first] && first > 0) ||
            (face == offsets[last + 1] && last + 1 < count)) {
          int side = face == offsets[first] ? 0 : 1;
          if (block->check_boundary[edge][side] && block->boundary_value[edge][side] != line->val) {
            gkyl_exit("Multiblock TS projection: conflicting interface biases.");
          }
          block->check_boundary[edge][side] = true;
          block->boundary_value[edge][side] = line->val;
          if (!block->boundary_minmax) {
            block->boundary_minmax = app->use_gpu ? gkyl_cu_malloc(8 * sizeof(double)) :
                                                    gkyl_malloc(8 * sizeof(double));
          }
          break;
        }
        planes[nplanes++] =
          (struct gkyl_fem_proj_bias_plane){.dir = 0, .loc = face, .val = line->val};
        break;
      }
    }
    struct gkyl_fem_proj_bias_plane_list biases = {.num_bias_plane = nplanes, .bp = planes};
    block->projector[edge] = gkyl_fem_proj_inew(&(struct gkyl_fem_proj_inp){
      .grid = &grid,
      .solve_range = &block->solve[edge],
      .basis = &block->basis,
      .bcs = bcs,
      .weight_left = block->weight,
      .weight_right = block->weight,
      .bias_plane_list = &biases,
      .use_gpu = app->use_gpu,
    });
    gkyl_free(planes);
  }
  block->inflate =
    gkyl_translate_dim_new(2, block->basis, 3, sbapp->basis, 2, GKYL_UPPER_EDGE, app->use_gpu);
}

static void
gk_multib_field_ts_new_connections(
  struct gkyl_gyrokinetic_multib_app *app, struct gk_multib_field_ts *ts, int bI, const int *ids,
  int count, const int *offsets
)
{
  struct gk_multib_field_ts_block *block = &ts->blocks[bI];
  int rank, nranks;
  gkyl_comm_get_rank(app->singleb_apps[bI]->comm, &rank);
  gkyl_comm_get_size(app->comm, &nranks);
  int *ranks = gkyl_malloc(nranks * sizeof(int));
  int nsend = 0;
  if (rank == 0) {
    for (int i = 0; i < count; ++i) {
      nsend += app->decomp[ids[i]]->ndecomp;
    }
  }
  struct gkyl_comm_conn *send = gkyl_malloc(GKYL_MAX2(1, nsend) * sizeof(*send));
  struct gkyl_comm_conn *recv = gkyl_malloc(count * sizeof(*recv));
  int next = 0;
  for (int i = 0; i < count; ++i) {
    gkyl_rrobin_decomp_getranks(app->round_robin, ids[i], ranks);
    recv[i] =
      (struct gkyl_comm_conn){.sr = GKYL_COMM_CONN_RECV, .block_id = ids[i], .rank = ranks[0]};
    gkyl_sub_range_init(
      &recv[i].range, &block->connected_ext, (int[]){offsets[i] + 1, 1},
      (int[]){offsets[i + 1], block->local.upper[1]}
    );
    if (rank == 0) {
      for (int r = 0; r < app->decomp[ids[i]]->ndecomp; ++r) {
        send[next++] = (struct gkyl_comm_conn){
          .sr = GKYL_COMM_CONN_SEND,
          .block_id = ids[i],
          .rank = ranks[r],
          .range = block->local,
        };
      }
    }
  }
  ts->send[bI] = gkyl_multib_comm_conn_new(nsend, send);
  ts->recv[bI] = gkyl_multib_comm_conn_new(count, recv);
  gkyl_multib_comm_conn_sort(ts->send[bI]);
  gkyl_multib_comm_conn_sort(ts->recv[bI]);
  gkyl_free(send);
  gkyl_free(recv);
  gkyl_free(ranks);
}

void
gk_multib_field_ts_new(
  const struct gkyl_gyrokinetic_multib *inp, struct gkyl_gyrokinetic_multib_app *app,
  struct gk_multib_field *mbf
)
{
  mbf->ts = NULL;
  if (mbf->cdim != 3) {
    return;
  }
  struct gk_multib_field_ts *ts = gkyl_calloc(1, sizeof(*ts));
  mbf->ts = ts;
  ts->num_blocks = app->num_local_blocks;
  ts->use_gpu = app->use_gpu;
  bool any_ts = false;
  for (int bid = 0; bid < app->block_topo->num_blocks; ++bid) {
    for (int edge = 0; edge < 2; ++edge) {
      any_ts |= gk_multib_field_ts_shifted(app, bid, edge);
    }
  }
  ts->project = !inp->field.ts_field_dg && any_ts;
  ts->blocks = gkyl_calloc(ts->num_blocks, sizeof(*ts->blocks));
  ts->send = gkyl_calloc(ts->num_blocks, sizeof(*ts->send));
  ts->recv = gkyl_calloc(ts->num_blocks, sizeof(*ts->recv));
  for (int edge = 0; edge < 2; ++edge) {
    ts->local[edge] = gkyl_calloc(ts->num_blocks, sizeof(struct gkyl_array *));
    ts->connected[edge] = gkyl_calloc(ts->num_blocks, sizeof(struct gkyl_array *));
    ts->fem[edge] = gkyl_calloc(ts->num_blocks, sizeof(struct gkyl_array *));
  }
  for (int bI = 0; bI < ts->num_blocks; ++bI) {
    int bid = app->local_blocks[bI];
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[bI];
    struct gk_multib_field_ts_block *block = &ts->blocks[bI];
    if (sbapp->gk_geom->has_LCFS) {
      gkyl_exit("Multiblock 3D field boundaries require core/SOL interfaces on block faces.");
    }
    for (int d = 0; d < 2; ++d) {
      if (sbapp->local.lower[d] != sbapp->global.lower[d] ||
          sbapp->local.upper[d] != sbapp->global.upper[d]) {
        gkyl_exit("Multiblock 3D field boundaries support decomposition only along z.");
      }
    }
    gkyl_range_init(
      &block->packed_ghost, 3, (int[]){1, 1, 1},
      (int[]){sbapp->grid.cells[0], sbapp->grid.cells[1], 1}
    );
    block->ghost_buffer = mkarr(app->use_gpu, sbapp->basis.num_basis, block->packed_ghost.volume);
    block->bc_buffer = mkarr(app->use_gpu, sbapp->basis.num_basis, block->packed_ghost.volume);
    for (int edge = 0; edge < 2; ++edge) {
      block->shifted[edge] = gk_multib_field_ts_shifted(app, bid, edge);
      int z = edge ? sbapp->global.upper[2] : sbapp->global.lower[2];
      block->ghost_owner[edge] = -1;
      for (int r = 0; r < sbapp->decomp->ndecomp; ++r) {
        const struct gkyl_range *range = &sbapp->decomp->ranges[r];
        if (range->lower[2] <= z && z <= range->upper[2]) {
          block->ghost_owner[edge] = r;
        }
      }
      block->physical[edge] = gkyl_bc_basic_gyrokinetic_new(
        2, edge, GKYL_BC_GK_FIELD_BOUNDARY_VALUE, sbapp->basis_on_dev,
        edge ? &sbapp->global_upper_skin[2] : &sbapp->global_lower_skin[2],
        edge ? &sbapp->global_upper_ghost[2] : &sbapp->global_lower_ghost[2],
        sbapp->basis.num_basis, 3, app->use_gpu
      );
    }
    gk_multib_field_ts_new_smoother(inp, sbapp, bid, block);
    if (!ts->project) {
      continue;
    }
    if (sbapp->poly_order != 1 || sbapp->basis.b_type != GKYL_BASIS_MODAL_SERENDIPITY ||
        mbf->gkfield_id != GKYL_GK_FIELD_ES || sbapp->field->use_flr) {
      gkyl_exit("Multiblock TS projection requires electrostatic p1 serendipity without FLR.");
    }
    if (sbapp->field->poisson_bcs.lo_type[1] != GKYL_POISSON_PERIODIC ||
        sbapp->field->poisson_bcs.up_type[1] != GKYL_POISSON_PERIODIC) {
      gkyl_exit("Multiblock TS projection requires periodic y field boundaries.");
    }
    gkyl_cart_modal_serendip(&block->basis, 2, 1);
    block->basis_on_dev = app->use_gpu ? gkyl_cart_modal_serendip_cu_dev_new(2, 1) : &block->basis;
    struct gkyl_rect_grid grid;
    gkyl_rect_grid_init(&grid, 2, sbapp->grid.lower, sbapp->grid.upper, sbapp->grid.cells);
    gkyl_create_grid_ranges(&grid, (int[]){1, 1}, &block->ext, &block->local);
    int count = gkyl_multib_conn_get_num_connected(app->block_topo, bid, 0, 0, GKYL_CONN_ALL);
    int ids[count], offsets[count + 1], own = -1;
    gkyl_multib_conn_get_connection(app->block_topo, bid, 0, 0, GKYL_CONN_ALL, ids);
    offsets[0] = 0;
    for (int i = 0; i < count; ++i) {
      offsets[i + 1] = offsets[i] + gkyl_range_shape(&app->decomp[ids[i]]->parent_range, 0);
      if (ids[i] == bid) {
        own = i;
      }
      const struct gkyl_block_connections *conn = &app->block_topo->conn[ids[i]];
      for (int side = 0; side < 2; ++side) {
        if (conn->connections[0][side].edge != GKYL_PHYSICAL &&
            (conn->connections[0][side].dir != 0 ||
             conn->connections[0][side].edge != (side ? GKYL_LOWER_POSITIVE : GKYL_UPPER_POSITIVE)
            )) {
          gkyl_exit("Multiblock TS projection requires aligned radial block connections.");
        }
      }
    }
    assert(own >= 0);
    gk_multib_field_ts_new_projection(inp, app, ts, bI, ids, count, own, offsets);
    gk_multib_field_ts_new_connections(app, ts, bI, ids, count, offsets);
  }
}

void
gk_multib_field_ts_gather(
  struct gkyl_gyrokinetic_multib_app *app, struct gk_multib_field *mbf, struct gkyl_array **local,
  struct gkyl_array **global
)
{
  for (int bI = 0; bI < app->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[bI];
    struct gk_multib_field_ts_block *block = &mbf->ts->blocks[bI];
    gkyl_comm_array_allgather(sbapp->comm, &sbapp->local, &sbapp->global, local[bI], global[bI]);
    int rank;
    gkyl_comm_get_rank(sbapp->comm, &rank);
    for (int edge = 0; edge < 2; ++edge) {
      if (rank == block->ghost_owner[edge]) {
        gkyl_array_copy_range_to_range(
          block->ghost_buffer, local[bI], &block->packed_ghost,
          edge ? &sbapp->local_upper_ghost[2] : &sbapp->local_lower_ghost[2]
        );
      }
      gkyl_comm_array_bcast(
        sbapp->comm, block->ghost_buffer, block->ghost_buffer, block->ghost_owner[edge]
      );
      gkyl_array_copy_range_to_range(
        global[bI], block->ghost_buffer,
        edge ? &sbapp->global_upper_ghost[2] : &sbapp->global_lower_ghost[2], &block->packed_ghost
      );
    }
  }
}

void
gk_multib_field_ts_fill(
  struct gkyl_gyrokinetic_multib_app *app, struct gk_multib_field *mbf, int bI,
  struct gkyl_array *global
)
{
  struct gk_multib_field_ts_block *block = &mbf->ts->blocks[bI];
  struct gk_field *field = app->singleb_apps[bI]->field;
  for (int edge = 0; edge < 2; ++edge) {
    if (block->shifted[edge]) {
      gkyl_bc_twistshift_advance(edge ? field->bc_ts_up : field->bc_ts_lo, global, global);
    } else {
      gkyl_bc_basic_gyrokinetic_advance(block->physical[edge], block->bc_buffer, global);
    }
  }
}

void
gk_multib_field_ts_project(
  struct gkyl_gyrokinetic_multib_app *app, struct gk_multib_field *mbf, struct gkyl_array **global
)
{
  struct gk_multib_field_ts *ts = mbf->ts;
  if (!ts->project) {
    return;
  }
  for (int edge = 0; edge < 2; ++edge) {
    for (int bI = 0; bI < ts->num_blocks; ++bI) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[bI];
      struct gk_multib_field_ts_block *block = &ts->blocks[bI];
      int z = edge ? sbapp->global.upper[2] + 1 : sbapp->global.lower[2] - 1;
      gkyl_deflate_zsurf_advance(
        block->deflate[edge], z, &sbapp->global_ext, &block->local, global[bI], ts->local[edge][bI],
        1
      );
    }
    gkyl_multib_comm_conn_array_transfer(
      app->comm, ts->num_blocks, app->local_blocks, ts->send, ts->recv, ts->local[edge],
      ts->connected[edge]
    );
    for (int bI = 0; bI < ts->num_blocks; ++bI) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[bI];
      struct gk_multib_field_ts_block *block = &ts->blocks[bI];
      if (!block->projector[edge]) {
        continue;
      }
      gkyl_array_copy(block->boundary[edge], ts->connected[edge][bI]);
      for (int side = 0; side < 2; ++side) {
        if (block->copy_boundary[edge][side]) {
          gkyl_bc_basic_gyrokinetic_advance(
            block->copy_boundary[edge][side], block->boundary_buffer[edge][side],
            block->boundary[edge]
          );
        }
      }
      for (int side = 0; side < 2; ++side) {
        if (block->check_boundary[edge][side]) {
          gkyl_array_reduce_range(
            block->boundary_minmax, block->boundary[edge], GKYL_MIN,
            &block->boundary_skin[edge][side]
          );
          gkyl_array_reduce_range(
            block->boundary_minmax + 4, block->boundary[edge], GKYL_MAX,
            &block->boundary_skin[edge][side]
          );
          double values[8];
          if (app->use_gpu) {
            gkyl_cu_memcpy(values, block->boundary_minmax, sizeof(values), GKYL_CU_MEMCPY_D2H);
          } else {
            memcpy(values, block->boundary_minmax, sizeof(values));
          }
          double value = block->boundary_value[edge][side];
          for (int k = 0; k < 8; ++k) {
            double expected = k % 4 == 0 ? 2 * value : 0;
            if (fabs(values[k] - expected) > 1e-10 * fmax(1.0, fabs(value))) {
              gkyl_exit("Multiblock TS projection: bias conflicts with the non-TS interface trace."
              );
            }
          }
        }
      }
      gkyl_fem_proj_set_rhs(block->projector[edge], ts->connected[edge][bI], block->boundary[edge]);
      gkyl_fem_proj_solve(block->projector[edge], ts->fem[edge][bI]);
      gkyl_array_copy_range_to_range(
        ts->local[edge][bI], ts->fem[edge][bI], &block->local, &block->block_in_connected
      );
      gkyl_translate_dim_advance(
        block->inflate, &block->local,
        edge ? &sbapp->global_upper_ghost[2] : &sbapp->global_lower_ghost[2], ts->local[edge][bI],
        1, global[bI]
      );
    }
  }
}

void
gk_multib_field_ts_smooth(
  struct gk_multib_field *mbf, int bI, struct gkyl_array *global, struct gkyl_array *smooth,
  bool is_potential
)
{
  struct gkyl_fem_parproj *proj = is_potential ? mbf->ts->blocks[bI].smooth :
                                                 mbf->ts->blocks[bI].smooth_rho;
  gkyl_fem_parproj_set_rhs(proj, global, global);
  gkyl_fem_parproj_solve(proj, smooth);
}

void
gk_multib_field_ts_release(struct gk_multib_field *mbf)
{
  struct gk_multib_field_ts *ts = mbf->ts;
  if (!ts) {
    return;
  }
  for (int bI = 0; bI < ts->num_blocks; ++bI) {
    struct gk_multib_field_ts_block *block = &ts->blocks[bI];
    gkyl_array_release(block->ghost_buffer);
    gkyl_array_release(block->bc_buffer);
    gkyl_fem_parproj_release(block->smooth);
    gkyl_fem_parproj_release(block->smooth_rho);
    for (int edge = 0; edge < 2; ++edge) {
      gkyl_bc_basic_gyrokinetic_release(block->physical[edge]);
    }
    if (!ts->project) {
      continue;
    }
    gkyl_multib_comm_conn_release(ts->send[bI]);
    gkyl_multib_comm_conn_release(ts->recv[bI]);
    if (block->boundary_minmax) {
      if (ts->use_gpu) {
        gkyl_cu_free(block->boundary_minmax);
      } else {
        gkyl_free(block->boundary_minmax);
      }
    }
    gkyl_array_release(block->weight);
    gkyl_translate_dim_release(block->inflate);
    for (int edge = 0; edge < 2; ++edge) {
      gkyl_array_release(ts->local[edge][bI]);
      gkyl_array_release(ts->connected[edge][bI]);
      gkyl_array_release(ts->fem[edge][bI]);
      gkyl_deflate_zsurf_release(block->deflate[edge]);
      if (!block->projector[edge]) {
        continue;
      }
      gkyl_fem_proj_release(block->projector[edge]);
      gkyl_array_release(block->boundary[edge]);
      for (int side = 0; side < 2; ++side) {
        if (block->copy_boundary[edge][side]) {
          gkyl_bc_basic_gyrokinetic_release(block->copy_boundary[edge][side]);
          gkyl_array_release(block->boundary_buffer[edge][side]);
        }
      }
    }
    if (block->basis_on_dev != &block->basis) {
      gkyl_cart_modal_basis_release_cu(block->basis_on_dev);
    }
  }
  for (int edge = 0; edge < 2; ++edge) {
    gkyl_free(ts->local[edge]);
    gkyl_free(ts->connected[edge]);
    gkyl_free(ts->fem[edge]);
  }
  gkyl_free(ts->send);
  gkyl_free(ts->recv);
  gkyl_free(ts->blocks);
  gkyl_free(ts);
}
