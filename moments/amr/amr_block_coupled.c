#include <gkyl_amr_block_coupled_priv.h>

void
five_moment_wall_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

void
ten_moment_wall_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 10; i++) {
    if (i == 1 || i == 5 || i == 6) {
      ghost[i] = -skin[i];
    }
    else {
      ghost[i] = skin[i];
    }
  }
}

void
maxwell_wall_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 8; i++) {
    if (i == 1 || i == 2 || i == 3 || i == 6) {
      ghost[i] = -skin[i];
    }
    else {
      ghost[i] = skin[i];
    }
  }
}

void
five_moment_copy_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++) {
    ghost[i] = skin[i];
  }
}

void
ten_moment_copy_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 10; i++) {
    ghost[i] = skin[i];
  }
}

void
maxwell_copy_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 8; i++) {
    ghost[i] = skin[i];
  }
}

void
five_moment_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[9];
  for (int i = 0; i < 9; i++) {
    nghost[i] = 2;
  }

  bool wall_x = bdata->wall_x;
  bool wall_y = bdata->wall_y;

  bool copy_x = bdata->copy_x;
  bool copy_y = bdata->copy_y;

  for (int d = 0; d < 2; d++) {
    bdata->lower_bc_elc[d] = bdata->upper_bc_elc[d] = 0;
    bdata->lower_bc_ion[d] = bdata->upper_bc_ion[d] = 0;
    bdata->lower_bc_maxwell[d] = bdata->upper_bc_maxwell[d] = 0;

    if ((d == 0 && wall_x) || (d == 1 && wall_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }
    }
    else if ((d == 0 && copy_x) || (d == 1 && copy_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }
    }
  }

  skin_ghost_ranges_init_block(&bdata->skin_ghost, &bdata->ext_range, nghost);
  long buff_sz = 0;

  for (int d = 0; d < 2; d++) {
    long vol = bdata->skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol) {
      buff_sz = vol;
    }
  }

  bdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  bdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  bdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
five_moment_nested_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[25];
  for (int i = 0; i < 25; i++) {
    nghost[i] = 2;
  }

  bool wall_x = bdata->wall_x;
  bool wall_y = bdata->wall_y;

  bool copy_x = bdata->copy_x;
  bool copy_y = bdata->copy_y;

  for (int d = 0; d < 2; d++) {
    bdata->lower_bc_elc[d] = bdata->upper_bc_elc[d] = 0;
    bdata->lower_bc_ion[d] = bdata->upper_bc_ion[d] = 0;
    bdata->lower_bc_maxwell[d] = bdata->upper_bc_maxwell[d] = 0;

    if ((d == 0 && wall_x) || (d == 1 && wall_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_wall_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }
    }
    else if ((d == 0 && copy_x) || (d == 1 && copy_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          five_moment_copy_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }
    }
  }

  skin_ghost_ranges_init_block(&bdata->skin_ghost, &bdata->ext_range, nghost);
  long buff_sz = 0;

  for (int d = 0; d < 2; d++) {
    long vol = bdata->skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol) {
      buff_sz = vol;
    }
  }

  bdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  bdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  bdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
ten_moment_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[9];
  for (int i = 0; i < 9; i++) {
    nghost[i] = 2;
  }

  bool wall_x = bdata->wall_x;
  bool wall_y = bdata->wall_y;

  bool copy_x = bdata->copy_x;
  bool copy_y = bdata->copy_y;

  for (int d = 0; d < 2; d++) {
    bdata->lower_bc_elc[d] = bdata->upper_bc_elc[d] = 0;
    bdata->lower_bc_ion[d] = bdata->upper_bc_ion[d] = 0;
    bdata->lower_bc_maxwell[d] = bdata->upper_bc_maxwell[d] = 0;

    if ((d == 0 && wall_x) || (d == 1 && wall_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }
    }
    else if ((d == 0 && copy_x) || (d == 1 && copy_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }
    }
  }

  skin_ghost_ranges_init_block(&bdata->skin_ghost, &bdata->ext_range, nghost);
  long buff_sz = 0;

  for (int d = 0; d < 2; d++) {
    long vol = bdata->skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol) {
      buff_sz = vol;
    }
  }

  bdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  bdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  bdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
ten_moment_nested_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[25];
  for (int i = 0; i < 25; i++) {
    nghost[i] = 2;
  }

  bool wall_x = bdata->wall_x;
  bool wall_y = bdata->wall_y;

  bool copy_x = bdata->copy_x;
  bool copy_y = bdata->copy_y;

  for (int d = 0; d < 2; d++) {
    bdata->lower_bc_elc[d] = bdata->upper_bc_elc[d] = 0;
    bdata->lower_bc_ion[d] = bdata->upper_bc_ion[d] = 0;
    bdata->lower_bc_maxwell[d] = bdata->upper_bc_maxwell[d] = 0;

    if ((d == 0 && wall_x) || (d == 1 && wall_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_wall_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_wall_bc, 0);
      }
    }
    else if ((d == 0 && copy_x) || (d == 1 && copy_y)) {
      if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
        bdata->lower_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->lower_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->lower_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }

      if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
        bdata->upper_bc_elc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_elc, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->upper_bc_ion[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler_ion, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          ten_moment_copy_bc, 0);
        bdata->upper_bc_maxwell[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->maxwell, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
          maxwell_copy_bc, 0);
      }
    }
  }

  skin_ghost_ranges_init_block(&bdata->skin_ghost, &bdata->ext_range, nghost);
  long buff_sz = 0;

  for (int d = 0; d < 2; d++) {
    long vol = bdata->skin_ghost.lower_skin[d].volume;
    
    if (buff_sz <= vol) {
      buff_sz = vol;
    }
  }

  bdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  bdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  bdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
five_moment_block_bc_updaters_release(struct five_moment_block_data* bdata)
{
  for (int d = 0; d < 2; d++) {
    if (bdata->lower_bc_elc[d]) {
      gkyl_wv_apply_bc_release(bdata->lower_bc_elc[d]);
    }
    if (bdata->lower_bc_ion[d]) {
      gkyl_wv_apply_bc_release(bdata->lower_bc_ion[d]);
    }
    if (bdata->lower_bc_maxwell[d]) {
      gkyl_wv_apply_bc_release(bdata->lower_bc_maxwell[d]);
    }

    if (bdata->upper_bc_elc[d]) {
      gkyl_wv_apply_bc_release(bdata->upper_bc_elc[d]);
    }
    if (bdata->upper_bc_ion[d]) {
      gkyl_wv_apply_bc_release(bdata->upper_bc_ion[d]);
    }
    if (bdata->upper_bc_maxwell[d]) {
      gkyl_wv_apply_bc_release(bdata->upper_bc_maxwell[d]);
    }
  }

  gkyl_array_release(bdata->bc_buffer_elc);
  gkyl_array_release(bdata->bc_buffer_ion);
  gkyl_array_release(bdata->bc_buffer_maxwell);
}

void
five_moment_block_bc_updaters_apply(const struct five_moment_block_data* bdata, double tm,
  struct gkyl_array* fld_elc, struct gkyl_array* fld_ion, struct gkyl_array* fld_maxwell)
{
  for (int d = 0; d < 2; d++) {
    if (bdata->lower_bc_elc[d]) {
      gkyl_wv_apply_bc_advance(bdata->lower_bc_elc[d], tm, &bdata->range, fld_elc);
    }
    if (bdata->lower_bc_ion[d]) {
      gkyl_wv_apply_bc_advance(bdata->lower_bc_ion[d], tm, &bdata->range, fld_ion);
    }
    if (bdata->lower_bc_maxwell[d]) {
      gkyl_wv_apply_bc_advance(bdata->lower_bc_maxwell[d], tm, &bdata->range, fld_maxwell);
    }

    if (bdata->upper_bc_elc[d]) {
      gkyl_wv_apply_bc_advance(bdata->upper_bc_elc[d], tm, &bdata->range, fld_elc);
    }
    if (bdata->upper_bc_ion[d]) {
      gkyl_wv_apply_bc_advance(bdata->upper_bc_ion[d], tm, &bdata->range, fld_ion);
    }
    if (bdata->upper_bc_maxwell[d]) {
      gkyl_wv_apply_bc_advance(bdata->upper_bc_maxwell[d], tm, &bdata->range, fld_maxwell);
    }
  }
}

void
block_coupled_ll_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

  double ref_factor_inv = ((double)bdata[i].skin_ghost.lower_skin[d].volume / (double)bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_ll_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

  int ref_factor = (int)(bdata[i].skin_ghost.lower_skin[d].volume / bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (ref_factor * count), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (ref_factor * count), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (ref_factor * count), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_lu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

  double ref_factor_inv = ((double)bdata[i].skin_ghost.lower_skin[d].volume / (double)bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_lu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

  int ref_factor = (int)(bdata[i].skin_ghost.lower_skin[d].volume / bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (ref_factor * count), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (ref_factor * count), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (ref_factor * count), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_ul_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

  double ref_factor_inv = ((double)bdata[i].skin_ghost.upper_skin[d].volume / (double)bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);
    
    if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_ul_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

  int ref_factor = (int)(bdata[i].skin_ghost.upper_skin[d].volume / bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (ref_factor * count), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (ref_factor * count), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (ref_factor * count), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_uu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

  double ref_factor_inv = ((double)bdata[i].skin_ghost.upper_skin[d].volume / (double)bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((int)(ref_factor_inv * count)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
block_coupled_uu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

  int ref_factor = (int)(bdata[i].skin_ghost.upper_skin[d].volume / bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);
    
    if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * (ref_factor * count), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * (ref_factor * count), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * (ref_factor * count), fld_maxwell[tbid]->esznc);
      count += 1;
    }
    else {
      memcpy(gkyl_array_fetch(fld_elc[tbid], start),
        ((char*) bc_buffer_elc->data) + fld_elc[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_elc[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_ion[tbid], start),
        ((char*) bc_buffer_ion->data) + fld_ion[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_ion[tbid]->esznc);
      memcpy(gkyl_array_fetch(fld_maxwell[tbid], start),
        ((char*) bc_buffer_maxwell->data) + fld_maxwell[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld_maxwell[tbid]->esznc);
      count += 1;
    }
  }
}

void
five_moment_sync_blocks(const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[],
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  int num_blocks = btopo->num_blocks;
  int ndim = btopo->ndim;

  for (int i = 0; i < num_blocks; i++) {
    for (int d = 0; d < ndim; d++) {
      const struct gkyl_target_edge *te = btopo->conn[i].connections[d];

      if (te[0].edge != GKYL_PHYSICAL) {
        struct gkyl_array *bc_buffer_elc = bdata[i].bc_buffer_elc;
        struct gkyl_array *bc_buffer_ion = bdata[i].bc_buffer_ion;
        struct gkyl_array *bc_buffer_maxwell = bdata[i].bc_buffer_maxwell;

        gkyl_array_copy_to_buffer(bc_buffer_elc->data, fld_elc[i], &(bdata[i].skin_ghost.lower_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_ion->data, fld_ion[i], &(bdata[i].skin_ghost.lower_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_maxwell->data, fld_maxwell[i], &(bdata[i].skin_ghost.lower_skin[d]));

        int tbid = te[0].bid;
        int tdir = te[0].dir;

        if (te[0].edge == GKYL_LOWER_POSITIVE) {
          if (bdata[i].skin_ghost.lower_skin[d].volume == bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume > bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            block_coupled_ll_restriction_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume < bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            block_coupled_ll_projection_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
        }
        else if (te[0].edge == GKYL_UPPER_POSITIVE) {
          if (bdata[i].skin_ghost.lower_skin[d].volume == bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume > bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            block_coupled_lu_restriction_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume < bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            block_coupled_lu_projection_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
        }
      }

      if (te[1].edge != GKYL_PHYSICAL) {
        struct gkyl_array *bc_buffer_elc = bdata[i].bc_buffer_elc;
        struct gkyl_array *bc_buffer_ion = bdata[i].bc_buffer_ion;
        struct gkyl_array *bc_buffer_maxwell = bdata[i].bc_buffer_maxwell;

        gkyl_array_copy_to_buffer(bc_buffer_elc->data, fld_elc[i], &(bdata[i].skin_ghost.upper_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_ion->data, fld_ion[i], &(bdata[i].skin_ghost.upper_skin[d]));
        gkyl_array_copy_to_buffer(bc_buffer_maxwell->data, fld_maxwell[i], &(bdata[i].skin_ghost.upper_skin[d]));

        int tbid = te[1].bid;
        int tdir = te[1].dir;

        if (te[1].edge == GKYL_LOWER_POSITIVE) {
          if (bdata[i].skin_ghost.upper_skin[d].volume == bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume > bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            block_coupled_ul_restriction_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume < bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            block_coupled_ul_projection_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
        }
        else if (te[1].edge == GKYL_UPPER_POSITIVE) {
          if (bdata[i].skin_ghost.upper_skin[d].volume == bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
            gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume > bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            block_coupled_uu_restriction_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume < bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            block_coupled_uu_projection_op(tbid, tdir, i, d, bdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
              fld_elc, fld_ion, fld_maxwell);
          }
        }
      }
    }
  }
}

void
five_moment_block_data_write(const char* file_nm_elc, const char* file_nm_ion, const char* file_nm_maxwell, const struct five_moment_block_data* bdata)
{
  gkyl_grid_sub_array_write(&bdata->grid, &bdata->range, 0, bdata->f_elc[0], file_nm_elc);
  gkyl_grid_sub_array_write(&bdata->grid, &bdata->range, 0, bdata->f_ion[0], file_nm_ion);
  gkyl_grid_sub_array_write(&bdata->grid, &bdata->range, 0, bdata->f_maxwell[0], file_nm_maxwell);
}

double
five_moment_block_data_max_dt(const struct five_moment_block_data* bdata)
{
  double dt = DBL_MAX;

  for (int d = 0; d < 2; d++) {
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata->slvr_elc[d], &bdata->range, bdata->f_elc[0]));
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata->slvr_ion[d], &bdata->range, bdata->f_ion[0]));
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata->slvr_maxwell[d], &bdata->range, bdata->f_maxwell[0]));
  }

  return dt;
}

void
five_moment_update_block_job_func(void* ctx)
{
  struct five_moment_update_block_ctx *ub_ctx = ctx;
  const struct five_moment_block_data *bdata = ub_ctx->bdata;

  int d = ub_ctx->dir;
  double t_curr = ub_ctx->t_curr;
  double dt = ub_ctx->dt;

  ub_ctx->stat_elc = gkyl_wave_prop_advance(bdata->slvr_elc[d], t_curr, dt, &bdata->range, NULL, bdata->f_elc[d], bdata->f_elc[d + 1]);
  ub_ctx->stat_ion = gkyl_wave_prop_advance(bdata->slvr_ion[d], t_curr, dt, &bdata->range, NULL, bdata->f_ion[d], bdata->f_ion[d + 1]);
  ub_ctx->stat_maxwell = gkyl_wave_prop_advance(bdata->slvr_maxwell[d], t_curr, dt, &bdata->range, NULL, bdata->f_maxwell[d], bdata->f_maxwell[d + 1]);

  five_moment_block_bc_updaters_apply(bdata, t_curr, bdata->f_elc[d + 1], bdata->f_ion[d + 1], bdata->f_maxwell[d + 1]);
}

void
five_moment_update_block_job_func_source(void* ctx)
{
  struct five_moment_update_block_ctx *ub_ctx = ctx;
  const struct five_moment_block_data *bdata = ub_ctx->bdata;

  int d = ub_ctx->dir;
  double t_curr = ub_ctx->t_curr;
  double dt = ub_ctx->dt;
  int nstrang = ub_ctx->nstrang;

  struct gkyl_array *fluids[2];
  fluids[0] = bdata->f_elc[nstrang];
  fluids[1] = bdata->f_ion[nstrang];

  const struct gkyl_array *app_accel[2];
  app_accel[0] = bdata->app_accel_elc;
  app_accel[1] = bdata->app_accel_ion;

  const struct gkyl_array *rhs_source[2];
  rhs_source[0] = bdata->rhs_source_elc;
  rhs_source[1] = bdata->rhs_source_ion;

  const struct gkyl_array *nT_source[2];
  nT_source[0] = bdata->nT_source_elc;
  nT_source[1] = bdata->nT_source_ion;

  gkyl_moment_em_coupling_implicit_advance(bdata->src_slvr, t_curr, dt, &bdata->range, fluids, app_accel, rhs_source,
    bdata->f_maxwell[nstrang], bdata->app_current, bdata->ext_em, nT_source);

  five_moment_block_bc_updaters_apply(bdata, t_curr, bdata->f_elc[nstrang], bdata->f_ion[nstrang], bdata->f_maxwell[nstrang]);
}

struct gkyl_update_status
five_moment_update_all_blocks(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct five_moment_block_data bdata[], double t_curr, double dt)
{
  int num_blocks = btopo->num_blocks;
  int ndim = btopo->ndim;

  double dt_suggested = DBL_MAX;

  for (int d = 0; d < ndim; d++) {
    struct five_moment_update_block_ctx five_moment_block_ctx[num_blocks];

    for (int i = 0; i < num_blocks; i++) {
      five_moment_block_ctx[i] = (struct five_moment_update_block_ctx) {
        .bdata = &bdata[i],
        .t_curr = t_curr,
        .dir = d,
        .dt = dt,
        .bidx = i,
        .nstrang = 0,
      };
    }

#ifdef AMR_USETHREADS
    for (int i = 0; i < num_blocks; i++) {
      gkyl_job_pool_add_work(job_pool, five_moment_update_block_job_func, &five_moment_block_ctx[i]);
    }
    gkyl_job_pool_wait(job_pool);
#else
    for (int i = 0; i < num_blocks; i++) {
      five_moment_update_block_job_func(&five_moment_block_ctx[i]);
    }
#endif

    struct gkyl_array *fld_elc[num_blocks];
    struct gkyl_array *fld_ion[num_blocks];
    struct gkyl_array *fld_maxwell[num_blocks];

    for (int i = 0; i < num_blocks; i++) {
      if (five_moment_block_ctx[i].stat_elc.success == false || five_moment_block_ctx[i].stat_ion.success == false || five_moment_block_ctx[i].stat_maxwell.success == false) {
        dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_elc.dt_suggested);
        dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_ion.dt_suggested);
        dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_maxwell.dt_suggested);

        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = dt_suggested,
        };
      }

      dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_elc.dt_suggested);
      dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_ion.dt_suggested);
      dt_suggested = fmin(dt_suggested, five_moment_block_ctx[i].stat_maxwell.dt_suggested);

      fld_elc[i] = bdata[i].f_elc[d + 1];
      fld_ion[i] = bdata[i].f_ion[d + 1];
      fld_maxwell[i] = bdata[i].f_maxwell[d + 1];
    }

    five_moment_sync_blocks(btopo, bdata, fld_elc, fld_ion, fld_maxwell);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested,
  };
}

void
five_moment_update_all_blocks_source(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct five_moment_block_data bdata[], double t_curr, double dt, int nstrang)
{
  int num_blocks = btopo->num_blocks;

  struct five_moment_update_block_ctx five_moment_block_ctx[num_blocks];

  for (int i = 0; i < num_blocks; i++) {
    five_moment_block_ctx[i] = (struct five_moment_update_block_ctx) {
      .bdata = &bdata[i],
      .t_curr = t_curr,
      .dir = 0,
      .dt = dt,
      .bidx = i,
      .nstrang = nstrang,
    };
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_blocks; i++) {
    gkyl_job_pool_add_work(job_pool, five_moment_update_block_job_func_source, &five_moment_block_ctx[i]);
  }
  gkyl_job_pool_wait(job_pool);
#else
  for (int i = 0; i < num_blocks; i++) {
    five_moment_update_block_job_func_source(&five_moment_block_ctx[i]);
  }
#endif

  struct gkyl_array *fld_elc[num_blocks];
  struct gkyl_array *fld_ion[num_blocks];
  struct gkyl_array *fld_maxwell[num_blocks];

  for (int i = 0; i < num_blocks; i++) {
    fld_elc[i] = bdata[i].f_elc[nstrang];
    fld_ion[i] = bdata[i].f_ion[nstrang];
    fld_maxwell[i] = bdata[i].f_maxwell[nstrang];
  }

  five_moment_sync_blocks(btopo, bdata, fld_elc, fld_ion, fld_maxwell);
}

void
five_moment_init_job_func_block(void* ctx)
{
  struct five_moment_block_data *bdata = ctx;

  gkyl_fv_proj_advance(bdata->fv_proj_elc, 0.0, &bdata->ext_range, bdata->f_elc[0]);
  gkyl_fv_proj_advance(bdata->fv_proj_ion, 0.0, &bdata->ext_range, bdata->f_ion[0]);
  gkyl_fv_proj_advance(bdata->fv_proj_maxwell, 0.0, &bdata->ext_range, bdata->f_maxwell[0]);
}

void
five_moment_copy_job_func(void* ctx)
{
  struct five_moment_copy_job_ctx *j_ctx = ctx;

  gkyl_array_copy(j_ctx->out_elc, j_ctx->inp_elc);
  gkyl_array_copy(j_ctx->out_ion, j_ctx->inp_ion);
  gkyl_array_copy(j_ctx->out_maxwell, j_ctx->inp_maxwell);
}

struct gkyl_update_status
five_moment_update_block(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct five_moment_block_data bdata[], double t_curr, double dt0, struct sim_stats* stats)
{
  int num_blocks = btopo->num_blocks;
  double dt_suggested = DBL_MAX;

  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    BLOCK_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct five_moment_copy_job_ctx five_moment_copy_ctx[num_blocks];
  double dt = dt0;

  while (state != UPDATE_DONE) {
    if (state == PRE_UPDATE) {
      state = FIRST_COUPLING_UPDATE;

      for (int i = 0; i < num_blocks; i++) {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = bdata[i].f_elc[0],
          .inp_ion = bdata[i].f_ion[0],
          .inp_maxwell = bdata[i].f_maxwell[0],
          .out_elc = bdata[i].fdup_elc,
          .out_ion = bdata[i].fdup_ion,
          .out_maxwell = bdata[i].fdup_maxwell,
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_blocks; i++) {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_blocks; i++) {
        five_moment_copy_job_func(&five_moment_copy_ctx[i]);
      }
#endif
    }
    else if (state == FIRST_COUPLING_UPDATE) {
      state = BLOCK_UPDATE;

      five_moment_update_all_blocks_source(job_pool, btopo, bdata, t_curr, 0.5 * dt, 0);
    }
    else if (state == BLOCK_UPDATE) {
      state = SECOND_COUPLING_UPDATE;

      struct gkyl_update_status s = five_moment_update_all_blocks(job_pool, btopo, bdata, t_curr, dt);

      if (!s.success) {
        stats->nfail += 1;
        dt = s.dt_suggested;
        state = UPDATE_REDO;
      }
      else {
        dt_suggested = fmin(dt_suggested, s.dt_suggested);
      }
    }
    else if (state == SECOND_COUPLING_UPDATE) {
      state = POST_UPDATE;

      five_moment_update_all_blocks_source(job_pool, btopo, bdata, t_curr, 0.5 * dt, 2);
    }
    else if (state == POST_UPDATE) {
      state = UPDATE_DONE;

      for (int i = 0; i < num_blocks; i++) {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = bdata[i].f_elc[2],
          .inp_ion = bdata[i].f_ion[2],
          .inp_maxwell = bdata[i].f_maxwell[2],
          .out_elc = bdata[i].f_elc[0],
          .out_ion = bdata[i].f_ion[0],
          .out_maxwell = bdata[i].f_maxwell[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_blocks; i++) {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_blocks; i++) {
        five_moment_copy_job_func(&five_moment_copy_ctx[i]);
      }
#endif
    }
    else if (state == UPDATE_REDO) {
      state = PRE_UPDATE;

      for (int i = 0; i < num_blocks; i++) {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = bdata[i].fdup_elc,
          .inp_ion = bdata[i].fdup_ion,
          .inp_maxwell = bdata[i].fdup_maxwell,
          .out_elc = bdata[i].f_elc[0],
          .out_ion = bdata[i].f_ion[0],
          .out_maxwell = bdata[i].f_maxwell[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_blocks; i++) {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_blocks; i++) {
        five_moment_copy_job_func(&five_moment_copy_ctx[i]);
      }
#endif
    }
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

void
five_moment_write_sol_block(const char* fbase, int num_blocks, const struct five_moment_block_data bdata[])
{
  for (int i = 0; i < num_blocks; i++) {
    const char *fmt_elc = "%s_elc_b%d.gkyl";
    int sz_elc = snprintf(0, 0, fmt_elc, fbase, i);
    char file_nm_elc[sz_elc + 1];

    const char *fmt_ion = "%s_ion_b%d.gkyl";
    int sz_ion = snprintf(0, 0, fmt_ion, fbase, i);
    char file_nm_ion[sz_ion + 1];

    const char *fmt_field = "%s_field_b%d.gkyl";
    int sz_field = snprintf(0, 0, fmt_field, fbase, i);
    char file_nm_field[sz_field + 1];

    snprintf(file_nm_elc, sizeof file_nm_elc, fmt_elc, fbase, i);
    snprintf(file_nm_ion, sizeof file_nm_ion, fmt_ion, fbase, i);
    snprintf(file_nm_field, sizeof file_nm_field, fmt_field, fbase, i);

    five_moment_block_data_write(file_nm_elc, file_nm_ion, file_nm_field, &bdata[i]);
  }
}

double
five_moment_max_dt_block(int num_blocks, const struct five_moment_block_data bdata[])
{
  double dt = DBL_MAX;

  for (int i = 0; i < num_blocks; i++) {
    dt = fmin(dt, five_moment_block_data_max_dt(&bdata[i]));
  }

  return dt;
}
