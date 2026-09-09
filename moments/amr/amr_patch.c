#include <gkyl_amr_patch_priv.h>
#include <gkyl_wv_euler_mixture_priv.h>

void
skin_ghost_ranges_init_patch(struct skin_ghost_ranges_patch* sgr, const struct gkyl_range* parent, const int* ghost)
{
  gkyl_skin_ghost_ranges(&sgr->lower_skin[0], &sgr->lower_ghost[0], 0, GKYL_LOWER_EDGE, parent, ghost);
  gkyl_skin_ghost_ranges(&sgr->upper_skin[0], &sgr->upper_ghost[0], 0, GKYL_UPPER_EDGE, parent, ghost);
}

void
euler_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[3];
  for (int i = 0; i < 3; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc[0] = pdata->upper_bc[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      euler_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      euler_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
}

void
euler_nested_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[5];
  for (int i = 0; i < 5; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc[0] = pdata->upper_bc[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      euler_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      euler_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
}

void
gr_euler_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[3];
  for (int i = 0; i < 3; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc[0] = pdata->upper_bc[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      gr_euler_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      gr_euler_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 29, buff_sz);
}

void
gr_euler_nested_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[5];
  for (int i = 0; i < 5; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc[0] = pdata->upper_bc[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      gr_euler_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      gr_euler_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 29, buff_sz);
}

void
euler_mixture_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  int nghost[3];
  for (int i = 0; i < 3; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc[0] = pdata->upper_bc[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      euler_mixture_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      euler_mixture_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), buff_sz);
}

void
euler_mixture_nested_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  int nghost[5];
  for (int i = 0; i < 5; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc[0] = pdata->upper_bc[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      euler_mixture_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      euler_mixture_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), buff_sz);
}

void
euler_patch_bc_updaters_release(struct euler_patch_data* pdata)
{
  if (pdata->lower_bc[0]) {
    gkyl_wv_apply_bc_release(pdata->lower_bc[0]);
  }
  
  if (pdata->upper_bc[0]) {
    gkyl_wv_apply_bc_release(pdata->upper_bc[0]);
  }

  gkyl_array_release(pdata->bc_buffer);
}

void
euler_patch_bc_updaters_apply(const struct euler_patch_data* pdata, double tm, struct gkyl_array* fld)
{
  if (pdata->lower_bc[0]) {
    gkyl_wv_apply_bc_advance(pdata->lower_bc[0], tm, &pdata->range, fld);
  }

  if (pdata->upper_bc[0]) {
    gkyl_wv_apply_bc_advance(pdata->upper_bc[0], tm, &pdata->range, fld);
  }
}

void
patch_ll_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.lower_skin[d].volume / (double)pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_ll_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.lower_skin[d].volume / pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_lu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.lower_skin[d].volume / (double)pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_lu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.lower_skin[d].volume / pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_ul_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.upper_skin[d].volume / (double)pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);
    
    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_ul_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.upper_skin[d].volume / pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_uu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.upper_skin[d].volume / (double)pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
patch_uu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.upper_skin[d].volume / pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);
    
    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
    }
    else {
      memcpy(gkyl_array_fetch(fld[tbid], start),
        ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
      count += 1;
    }
  }
}

void
euler_sync_patches(const struct gkyl_block_topo* ptopo, const struct euler_patch_data pdata[], struct gkyl_array* fld[])
{
  int num_patches = ptopo->num_blocks;

  for (int i = 0; i < num_patches; i++) {
    const struct gkyl_target_edge *te = ptopo->conn[i].connections[0];

    if (te[0].edge != GKYL_PHYSICAL) {
      struct gkyl_array *bc_buffer = pdata[i].bc_buffer;
      
      gkyl_array_copy_to_buffer(bc_buffer->data, fld[i], &(pdata[i].skin_ghost.lower_skin[0]));

      int tbid = te[0].bid;
      int tdir = te[0].dir;

      if (te[0].edge == GKYL_LOWER_POSITIVE) {
        if (pdata[i].skin_ghost.lower_skin[0].volume == pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume > pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_ll_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume < pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_ll_projection_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
      }
      else if (te[0].edge == GKYL_UPPER_POSITIVE) {
        if (pdata[i].skin_ghost.lower_skin[0].volume == pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume > pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_lu_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume < pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_lu_projection_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
      }
    }

    if (te[1].edge != GKYL_PHYSICAL) {
      struct gkyl_array *bc_buffer = pdata[i].bc_buffer;

      gkyl_array_copy_to_buffer(bc_buffer->data, fld[i], &(pdata[i].skin_ghost.upper_skin[0]));

      int tbid = te[1].bid;
      int tdir = te[1].dir;

      if (te[1].edge == GKYL_LOWER_POSITIVE) {
        if (pdata[i].skin_ghost.upper_skin[0].volume == pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume > pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_ul_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume < pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_ul_projection_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
      }
      else if (te[1].edge == GKYL_UPPER_POSITIVE) {
        if (pdata[i].skin_ghost.upper_skin[0].volume == pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume > pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_uu_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume < pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_uu_projection_op(tbid, tdir, i, 0, pdata, bc_buffer, fld);
        }
      }
    }
  }
}

void
euler_patch_data_write(const char* file_nm, const struct euler_patch_data* pdata)
{
  gkyl_grid_sub_array_write(&pdata->grid, &pdata->range, 0, pdata->f[0], file_nm);
}

double
euler_patch_data_max_dt(const struct euler_patch_data* pdata)
{
  double dt = DBL_MAX;

  dt = fmin(dt, gkyl_wave_prop_max_dt(pdata->slvr[0], &pdata->range, pdata->f[0]));

  return dt;
}

void
euler_update_patch_job_func(void* ctx)
{
  struct euler_update_patch_ctx *up_ctx = ctx;
  const struct euler_patch_data *pdata = up_ctx->pdata;

  int d = up_ctx->dir;
  double t_curr = up_ctx->t_curr;
  double dt = up_ctx->dt;

  up_ctx->stat = gkyl_wave_prop_advance(pdata->slvr[d], t_curr, dt, &pdata->range, NULL, pdata->f[d], pdata->f[d + 1]);

  euler_patch_bc_updaters_apply(pdata, t_curr, pdata->f[d + 1]);
}

struct gkyl_update_status
euler_update_all_patches(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct euler_patch_data pdata[], double t_curr, double dt)
{
  int num_patches = ptopo->num_blocks;
  double dt_suggested = DBL_MAX;

  struct euler_update_patch_ctx euler_patch_ctx[num_patches];

  for (int i = 0; i < num_patches; i++) {
    euler_patch_ctx[i] = (struct euler_update_patch_ctx) {
      .pdata = &pdata[i],
      .t_curr = t_curr,
      .dir = 0,
      .dt = dt,
      .pidx = i,
    };
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_patches; i++) {
    gkyl_job_pool_add_work(job_pool, euler_update_patch_job_func, &euler_patch_ctx[i]);
  }
  gkyl_job_pool_wait(job_pool);
#else
  for (int i = 0; i < num_patches; i++) {
    euler_update_patch_job_func(&euler_patch_ctx[i]);
  }
#endif
  
  struct gkyl_array *fld[num_patches];

  for (int i = 0; i < num_patches; i++) {
    if (euler_patch_ctx[i].stat.success == false) {
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = euler_patch_ctx[i].stat.dt_suggested,
      };
    }

    dt_suggested = fmin(dt_suggested, euler_patch_ctx[i].stat.dt_suggested);
    fld[i] = pdata[i].f[1];
  }

  euler_sync_patches(ptopo, pdata, fld);

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested,
  };
}

void
euler_init_job_func_patch(void* ctx)
{
  struct euler_patch_data *pdata = ctx;

  gkyl_fv_proj_advance(pdata->fv_proj, 0.0, &pdata->ext_range, pdata->f[0]);
}

struct gkyl_update_status
euler_update_patch(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct euler_patch_data pdata[], double t_curr, double dt0, struct sim_stats* stats)
{
  int num_patches = ptopo->num_blocks;
  double dt_suggested = DBL_MAX;

  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FLUID_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct copy_job_ctx euler_copy_ctx[num_patches];
  double dt = dt0;

  while (state != UPDATE_DONE) {
    if (state == PRE_UPDATE) {
      state = FLUID_UPDATE;

      for (int i = 0; i < num_patches; i++) {
        euler_copy_ctx[i] = (struct copy_job_ctx) {
          .bidx = i,
          .inp = pdata[i].f[0],
          .out = pdata[i].fdup,
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_patches; i++) {
        gkyl_job_pool_add_work(job_pool, copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_patches; i++) {
        copy_job_func(&euler_copy_ctx[i]);
      }
#endif
    }
    else if (state == FLUID_UPDATE) {
      state = POST_UPDATE;

      struct gkyl_update_status s = euler_update_all_patches(job_pool, ptopo, pdata, t_curr, dt);

      if (!s.success) {
        stats->nfail += 1;
        dt = s.dt_suggested;
        state = UPDATE_REDO;
      }
      else {
        dt_suggested = fmin(dt_suggested, s.dt_suggested);
      }
    }
    else if (state == POST_UPDATE) {
      state = UPDATE_DONE;

      for (int i = 0; i < num_patches; i++) {
        euler_copy_ctx[i] = (struct copy_job_ctx) {
          .bidx = i,
          .inp = pdata[i].f[1],
          .out = pdata[i].f[0],
        };
      }

#ifdef AMR_USETHREADAS
      for (int i = 0; i < num_patches; i++) {
        gkyl_job_pool_add_work(job_pool, copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_patches; i++) {
        copy_job_func(&euler_copy_ctx[i]);
      }
#endif
    }
    else if (state == UPDATE_REDO) {
      state = PRE_UPDATE;

      for (int i = 0; i < num_patches; i++) {
        euler_copy_ctx[i] = (struct copy_job_ctx) {
          .bidx = i,
          .inp = pdata[i].fdup,
          .out = pdata[i].f[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_patches; i++) {
        gkyl_job_pool_add_work(job_pool, copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_patches; i++) {
        copy_job_func(&euler_copy_ctx[i]);
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
euler_write_sol_patch(const char* fbase, int num_patches, const struct euler_patch_data pdata[])
{
  for (int i = 0; i < num_patches; i++) {
    const char *fmt = "%s_p%d.gkyl";
    int sz = snprintf(0, 0, fmt, fbase, i);
    char file_nm[sz + 1];

    snprintf(file_nm, sizeof file_nm, fmt, fbase, i);
    euler_patch_data_write(file_nm, &pdata[i]);
  }
}

double
euler_max_dt_patch(int num_patches, const struct euler_patch_data pdata[])
{
  double dt = DBL_MAX;

  for (int i = 0; i < num_patches; i++) {
    dt = fmin(dt, euler_patch_data_max_dt(&pdata[i]));
  }

  return dt;
}

struct gkyl_block_topo*
create_patch_topo()
{
  struct gkyl_block_topo *ptopo = gkyl_block_topo_new(1, 3);

  ptopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
  };
  ptopo->conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
  };
  ptopo->conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
  };

  return ptopo;
}

struct gkyl_block_topo*
create_nested_patch_topo()
{
  struct gkyl_block_topo *ptopo = gkyl_block_topo_new(1, 5);

  ptopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
  };
  ptopo->conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
  };
  ptopo->conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
  };
  ptopo->conn[3] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
  };
  ptopo->conn[4] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
  };

  return ptopo;
}
