#include <gkyl_amr_block_coupled_priv.h>
#include <gkyl_amr_patch_coupled_priv.h>

void
five_moment_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[3];
  for (int i = 0; i < 3; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc_elc[0] = pdata->upper_bc_elc[0] = 0;
  pdata->lower_bc_ion[0] = pdata->upper_bc_ion[0] = 0;
  pdata->lower_bc_maxwell[0] = pdata->upper_bc_maxwell[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->lower_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->lower_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->upper_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->upper_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  pdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  pdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
five_moment_nested_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[5];
  for (int i = 0; i < 5; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc_elc[0] = pdata->upper_bc_elc[0] = 0;
  pdata->lower_bc_ion[0] = pdata->upper_bc_ion[0] = 0;
  pdata->lower_bc_maxwell[0] = pdata->upper_bc_maxwell[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->lower_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->lower_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->upper_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      five_moment_copy_bc, 0);
    pdata->upper_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  pdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
  pdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
ten_moment_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[3];
  for (int i = 0; i < 3; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc_elc[0] = pdata->upper_bc_elc[0] = 0;
  pdata->lower_bc_ion[0] = pdata->upper_bc_ion[0] = 0;
  pdata->lower_bc_maxwell[0] = pdata->upper_bc_maxwell[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->lower_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->lower_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->upper_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->upper_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  pdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  pdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
ten_moment_nested_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn)
{
  int nghost[5];
  for (int i = 0; i < 5; i++) {
    nghost[i] = 2;
  }

  pdata->lower_bc_elc[0] = pdata->upper_bc_elc[0] = 0;
  pdata->lower_bc_ion[0] = pdata->upper_bc_ion[0] = 0;
  pdata->lower_bc_maxwell[0] = pdata->upper_bc_maxwell[0] = 0;

  if (conn->connections[0][0].edge == GKYL_PHYSICAL) {
    pdata->lower_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->lower_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->lower_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_LOWER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  if (conn->connections[0][1].edge == GKYL_PHYSICAL) {
    pdata->upper_bc_elc[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_elc, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->upper_bc_ion[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->euler_ion, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      ten_moment_copy_bc, 0);
    pdata->upper_bc_maxwell[0] = gkyl_wv_apply_bc_new(&pdata->grid, pdata->maxwell, pdata->geom, 0, GKYL_UPPER_EDGE, nghost,
      maxwell_copy_bc, 0);
  }

  skin_ghost_ranges_init_patch(&pdata->skin_ghost, &pdata->ext_range, nghost);
  long buff_sz = 0;

  long vol = pdata->skin_ghost.lower_skin[0].volume;
  
  if (buff_sz <= vol) {
    buff_sz = vol;
  }

  pdata->bc_buffer_elc = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  pdata->bc_buffer_ion = gkyl_array_new(GKYL_DOUBLE, 10, buff_sz);
  pdata->bc_buffer_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, buff_sz);
}

void
five_moment_patch_bc_updaters_release(struct five_moment_patch_data* pdata)
{
  if (pdata->lower_bc_elc[0]) {
    gkyl_wv_apply_bc_release(pdata->lower_bc_elc[0]);
  }
  if (pdata->lower_bc_ion[0]) {
    gkyl_wv_apply_bc_release(pdata->lower_bc_ion[0]);
  }
  if (pdata->lower_bc_maxwell[0]) {
    gkyl_wv_apply_bc_release(pdata->lower_bc_maxwell[0]);
  }
  
  if (pdata->upper_bc_elc[0]) {
    gkyl_wv_apply_bc_release(pdata->upper_bc_elc[0]);
  }
  if (pdata->upper_bc_ion[0]) {
    gkyl_wv_apply_bc_release(pdata->upper_bc_ion[0]);
  }
  if (pdata->upper_bc_maxwell[0]) {
    gkyl_wv_apply_bc_release(pdata->upper_bc_maxwell[0]);
  }

  gkyl_array_release(pdata->bc_buffer_elc);
  gkyl_array_release(pdata->bc_buffer_ion);
  gkyl_array_release(pdata->bc_buffer_maxwell);
}

void
five_moment_patch_bc_updaters_apply(const struct five_moment_patch_data* pdata, double tm,
  struct gkyl_array* fld_elc, struct gkyl_array* fld_ion, struct gkyl_array* fld_maxwell)
{
  if (pdata->lower_bc_elc[0]) {
    gkyl_wv_apply_bc_advance(pdata->lower_bc_elc[0], tm, &pdata->range, fld_elc);
  }
  if (pdata->lower_bc_ion[0]) {
    gkyl_wv_apply_bc_advance(pdata->lower_bc_ion[0], tm, &pdata->range, fld_ion);
  }
  if (pdata->lower_bc_maxwell[0]) {
    gkyl_wv_apply_bc_advance(pdata->lower_bc_maxwell[0], tm, &pdata->range, fld_maxwell);
  }
  
  if (pdata->upper_bc_elc[0]) {
    gkyl_wv_apply_bc_advance(pdata->upper_bc_elc[0], tm, &pdata->range, fld_elc);
  }
  if (pdata->upper_bc_ion[0]) {
    gkyl_wv_apply_bc_advance(pdata->upper_bc_ion[0], tm, &pdata->range, fld_ion);
  }
  if (pdata->upper_bc_maxwell[0]) {
    gkyl_wv_apply_bc_advance(pdata->upper_bc_maxwell[0], tm, &pdata->range, fld_maxwell);
  }
}

void
patch_coupled_ll_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.lower_skin[d].volume / (double)pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_ll_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.lower_skin[d].volume / pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_lu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.lower_skin[d].volume / (double)pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_lu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.lower_skin[d].volume / pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_ul_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.upper_skin[d].volume / (double)pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);
    
    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_ul_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.upper_skin[d].volume / pdata[tbid].skin_ghost.lower_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_uu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  double ref_factor_inv = ((double)pdata[i].skin_ghost.upper_skin[d].volume / (double)pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
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
patch_coupled_uu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));

  int ref_factor = (int)(pdata[i].skin_ghost.upper_skin[d].volume / pdata[tbid].skin_ghost.upper_ghost[tdir].volume);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&(pdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);
    
    if ((pdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - pdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
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
five_moment_sync_patches(const struct gkyl_block_topo* ptopo, const struct five_moment_patch_data pdata[],
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[])
{
  int num_patches = ptopo->num_blocks;

  for (int i = 0; i < num_patches; i++) {
    const struct gkyl_target_edge *te = ptopo->conn[i].connections[0];

    if (te[0].edge != GKYL_PHYSICAL) {
      struct gkyl_array *bc_buffer_elc = pdata[i].bc_buffer_elc;
      struct gkyl_array *bc_buffer_ion = pdata[i].bc_buffer_ion;
      struct gkyl_array *bc_buffer_maxwell = pdata[i].bc_buffer_maxwell;

      gkyl_array_copy_to_buffer(bc_buffer_elc->data, fld_elc[i], &(pdata[i].skin_ghost.lower_skin[0]));
      gkyl_array_copy_to_buffer(bc_buffer_ion->data, fld_ion[i], &(pdata[i].skin_ghost.lower_skin[0]));
      gkyl_array_copy_to_buffer(bc_buffer_maxwell->data, fld_maxwell[i], &(pdata[i].skin_ghost.lower_skin[0]));

      int tbid = te[0].bid;
      int tdir = te[0].dir;

      if (te[0].edge == GKYL_LOWER_POSITIVE) {
        if (pdata[i].skin_ghost.lower_skin[0].volume == pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume > pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_coupled_ll_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume < pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_coupled_ll_projection_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
      }
      else if (te[0].edge == GKYL_UPPER_POSITIVE) {
        if (pdata[i].skin_ghost.lower_skin[0].volume == pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume > pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_coupled_lu_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
        else if (pdata[i].skin_ghost.lower_skin[0].volume < pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_coupled_lu_projection_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
      }
    }

    if (te[1].edge != GKYL_PHYSICAL) {
      struct gkyl_array *bc_buffer_elc = pdata[i].bc_buffer_elc;
      struct gkyl_array *bc_buffer_ion = pdata[i].bc_buffer_ion;
      struct gkyl_array *bc_buffer_maxwell = pdata[i].bc_buffer_maxwell;

      gkyl_array_copy_to_buffer(bc_buffer_elc->data, fld_elc[i], &(pdata[i].skin_ghost.upper_skin[0]));
      gkyl_array_copy_to_buffer(bc_buffer_ion->data, fld_ion[i], &(pdata[i].skin_ghost.upper_skin[0]));
      gkyl_array_copy_to_buffer(bc_buffer_maxwell->data, fld_maxwell[i], &(pdata[i].skin_ghost.upper_skin[0]));

      int tbid = te[1].bid;
      int tdir = te[1].dir;

      if (te[1].edge == GKYL_LOWER_POSITIVE) {
        if (pdata[i].skin_ghost.upper_skin[0].volume == pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(pdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume > pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_coupled_ul_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume < pdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
          patch_coupled_ul_projection_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
      }
      else if (te[1].edge == GKYL_UPPER_POSITIVE) {
      if (pdata[i].skin_ghost.upper_skin[0].volume == pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          gkyl_array_copy_from_buffer(fld_elc[tbid], bc_buffer_elc->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_ion[tbid], bc_buffer_ion->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
          gkyl_array_copy_from_buffer(fld_maxwell[tbid], bc_buffer_maxwell->data, &(pdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume > pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_coupled_uu_restriction_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
        else if (pdata[i].skin_ghost.upper_skin[0].volume < pdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
          patch_coupled_uu_projection_op(tbid, tdir, i, 0, pdata, bc_buffer_elc, bc_buffer_ion, bc_buffer_maxwell,
            fld_elc, fld_ion, fld_maxwell);
        }
      }
    }
  }
}

void
five_moment_patch_data_write(const char* file_nm_elc, const char* file_nm_ion, const char* file_nm_maxwell, const struct five_moment_patch_data* pdata)
{
  gkyl_grid_sub_array_write(&pdata->grid, &pdata->range, 0, pdata->f_elc[0], file_nm_elc);
  gkyl_grid_sub_array_write(&pdata->grid, &pdata->range, 0, pdata->f_ion[0], file_nm_ion);
  gkyl_grid_sub_array_write(&pdata->grid, &pdata->range, 0, pdata->f_maxwell[0], file_nm_maxwell);
}

double
five_moment_patch_data_max_dt(const struct five_moment_patch_data* pdata)
{
  double dt = DBL_MAX;

  dt = fmin(dt, gkyl_wave_prop_max_dt(pdata->slvr_elc[0], &pdata->range, pdata->f_elc[0]));
  dt = fmin(dt, gkyl_wave_prop_max_dt(pdata->slvr_ion[0], &pdata->range, pdata->f_ion[0]));
  dt = fmin(dt, gkyl_wave_prop_max_dt(pdata->slvr_maxwell[0], &pdata->range, pdata->f_maxwell[0]));

  return dt;
}

void
five_moment_update_patch_job_func(void* ctx)
{
  struct five_moment_update_patch_ctx *up_ctx = ctx;
  const struct five_moment_patch_data *pdata = up_ctx->pdata;

  int d = up_ctx->dir;
  double t_curr = up_ctx->t_curr;
  double dt = up_ctx->dt;

  up_ctx->stat_elc = gkyl_wave_prop_advance(pdata->slvr_elc[d], t_curr, dt, &pdata->range, NULL, pdata->f_elc[d], pdata->f_elc[d + 1]);
  up_ctx->stat_ion = gkyl_wave_prop_advance(pdata->slvr_ion[d], t_curr, dt, &pdata->range, NULL, pdata->f_ion[d], pdata->f_ion[d + 1]);
  up_ctx->stat_maxwell = gkyl_wave_prop_advance(pdata->slvr_maxwell[d], t_curr, dt, &pdata->range, NULL, pdata->f_maxwell[d], pdata->f_maxwell[d + 1]);

  five_moment_patch_bc_updaters_apply(pdata, t_curr, pdata->f_elc[d + 1], pdata->f_ion[d + 1], pdata->f_maxwell[d + 1]);
}

void
five_moment_update_patch_job_func_source(void* ctx)
{
  struct five_moment_update_patch_ctx *up_ctx = ctx;
  const struct five_moment_patch_data *pdata = up_ctx->pdata;

  int d = up_ctx->dir;
  double t_curr = up_ctx->t_curr;
  double dt = up_ctx->dt;
  int nstrang = up_ctx->nstrang;

  struct gkyl_array *fluids[2];
  fluids[0] = pdata->f_elc[nstrang];
  fluids[1] = pdata->f_ion[nstrang];

  const struct gkyl_array *app_accel[2];
  app_accel[0] = pdata->app_accel_elc;
  app_accel[1] = pdata->app_accel_ion;

  const struct gkyl_array *rhs_source[2];
  rhs_source[0] = pdata->rhs_source_elc;
  rhs_source[1] = pdata->rhs_source_ion;

  const struct gkyl_array *nT_source[2];
  nT_source[0] = pdata->nT_source_elc;
  nT_source[1] = pdata->nT_source_ion;
  
  gkyl_moment_em_coupling_implicit_advance(pdata->src_slvr, t_curr, dt, &pdata->range, fluids, app_accel, rhs_source,
    pdata->f_maxwell[nstrang], pdata->app_current, pdata->ext_em, nT_source);

  five_moment_patch_bc_updaters_apply(pdata, t_curr, pdata->f_elc[nstrang], pdata->f_ion[nstrang], pdata->f_maxwell[nstrang]);
}

struct gkyl_update_status
five_moment_update_all_patches(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct five_moment_patch_data pdata[], double t_curr, double dt)
{
  int num_patches = ptopo->num_blocks;
  double dt_suggested = DBL_MAX;

  struct five_moment_update_patch_ctx five_moment_patch_ctx[num_patches];

  for (int i = 0; i < num_patches; i++) {
    five_moment_patch_ctx[i] = (struct five_moment_update_patch_ctx) {
      .pdata = &pdata[i],
      .t_curr = t_curr,
      .dir = 0,
      .dt = dt,
      .pidx = i,
      .nstrang = 0,
    };
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_patches; i++) {
    gkyl_job_pool_add_work(job_pool, five_moment_update_patch_job_func, &five_moment_patch_ctx[i]);
  }
  gkyl_job_pool_wait(job_pool);
#else
  for (int i = 0; i < num_patches; i++) {
    five_moment_update_patch_job_func(&five_moment_patch_ctx[i]);
  }
#endif

  struct gkyl_array *fld_elc[num_patches];
  struct gkyl_array *fld_ion[num_patches];
  struct gkyl_array *fld_maxwell[num_patches];

  for (int i = 0; i < num_patches; i++) {
    if (five_moment_patch_ctx[i].stat_elc.success == false || five_moment_patch_ctx[i].stat_ion.success == false || five_moment_patch_ctx[i].stat_maxwell.success == false) {
      dt_suggested = fmin(dt_suggested, five_moment_patch_ctx[i].stat_elc.dt_suggested);
      dt_suggested = fmin(dt_suggested, five_moment_patch_ctx[i].stat_ion.dt_suggested);
      dt_suggested = fmin(dt_suggested, five_moment_patch_ctx[i].stat_maxwell.dt_suggested);

      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = dt_suggested,
      };
    }

    dt_suggested = fmin(dt_suggested, five_moment_patch_ctx[i].stat_elc.dt_suggested);
    dt_suggested = fmin(dt_suggested, five_moment_patch_ctx[i].stat_ion.dt_suggested);
    dt_suggested = fmin(dt_suggested, five_moment_patch_ctx[i].stat_maxwell.dt_suggested);

    fld_elc[i] = pdata[i].f_elc[1];
    fld_ion[i] = pdata[i].f_ion[1];
    fld_maxwell[i] = pdata[i].f_maxwell[1];
  }

  five_moment_sync_patches(ptopo, pdata, fld_elc, fld_ion, fld_maxwell);

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested,
  };
}

void
five_moment_update_all_patches_source(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct five_moment_patch_data pdata[], double t_curr, double dt, int nstrang)
{
  int num_patches = ptopo->num_blocks;

  struct five_moment_update_patch_ctx five_moment_patch_ctx[num_patches];

  for (int i = 0; i < num_patches; i++) {
    five_moment_patch_ctx[i] = (struct five_moment_update_patch_ctx) {
      .pdata = &pdata[i],
      .t_curr = t_curr,
      .dir = 0,
      .dt = dt,
      .pidx = i,
      .nstrang = nstrang,
    };
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_patches; i++) {
    gkyl_job_pool_add_work(job_pool, five_moment_update_patch_job_func_source, &five_moment_patch_ctx[i]);
  }
  gkyl_job_pool_wait(job_pool);
#else
  for (int i = 0; i < num_patches; i++) {
    five_moment_update_patch_job_func_source(&five_moment_patch_ctx[i]);
  }
#endif

  struct gkyl_array *fld_elc[num_patches];
  struct gkyl_array *fld_ion[num_patches];
  struct gkyl_array *fld_maxwell[num_patches];

  for (int i = 0; i < num_patches; i++) {
    fld_elc[i] = pdata[i].f_elc[nstrang];
    fld_ion[i] = pdata[i].f_ion[nstrang];
    fld_maxwell[i] = pdata[i].f_maxwell[nstrang];
  }

  five_moment_sync_patches(ptopo, pdata, fld_elc, fld_ion, fld_maxwell);
}

void
five_moment_init_job_func_patch(void* ctx)
{
  struct five_moment_patch_data *pdata = ctx;

  gkyl_fv_proj_advance(pdata->fv_proj_elc, 0.0, &pdata->ext_range, pdata->f_elc[0]);
  gkyl_fv_proj_advance(pdata->fv_proj_ion, 0.0, &pdata->ext_range, pdata->f_ion[0]);
  gkyl_fv_proj_advance(pdata->fv_proj_maxwell, 0.0, &pdata->ext_range, pdata->f_maxwell[0]);
}

struct gkyl_update_status
five_moment_update_patch(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct five_moment_patch_data pdata[], double t_curr, double dt0, struct sim_stats* stats)
{
  int num_patches = ptopo->num_blocks;
  double dt_suggested = DBL_MAX;

  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    PATCH_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct five_moment_copy_job_ctx five_moment_copy_ctx[num_patches];
  double dt = dt0;

  while (state != UPDATE_DONE) {
    if (state == PRE_UPDATE) {
      state = FIRST_COUPLING_UPDATE;

      for (int i = 0; i < num_patches; i++) {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = pdata[i].f_elc[0],
          .inp_ion = pdata[i].f_ion[0],
          .inp_maxwell = pdata[i].f_maxwell[0],
          .out_elc = pdata[i].fdup_elc,
          .out_ion = pdata[i].fdup_ion,
          .out_maxwell = pdata[i].fdup_maxwell,
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_patches; i++) {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_patches; i++) {
        five_moment_copy_job_func(&five_moment_copy_ctx[i]);
      }
#endif
    }
    else if (state == FIRST_COUPLING_UPDATE) {
      state = PATCH_UPDATE;

      five_moment_update_all_patches_source(job_pool, ptopo, pdata, t_curr, 0.5 * dt, 0);
    }
    else if (state == PATCH_UPDATE) {
      state = SECOND_COUPLING_UPDATE;

      struct gkyl_update_status s = five_moment_update_all_patches(job_pool, ptopo, pdata, t_curr, dt);

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

      five_moment_update_all_patches_source(job_pool, ptopo, pdata, t_curr, 0.5 * dt, 2);
    }
    else if (state == POST_UPDATE) {
      state = UPDATE_DONE;

      for (int i = 0; i < num_patches; i++) {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = pdata[i].f_elc[1],
          .inp_ion = pdata[i].f_ion[1],
          .inp_maxwell = pdata[i].f_maxwell[1],
          .out_elc = pdata[i].f_elc[0],
          .out_ion = pdata[i].f_ion[0],
          .out_maxwell = pdata[i].f_maxwell[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_patches; i++) {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_patches; i++) {
        five_moment_copy_job_func(&five_moment_copy_ctx[i]);
      }
#endif
    }
    else if (state == UPDATE_REDO) {
      state = PRE_UPDATE;

      for (int i = 0; i < num_patches; i++) {
        five_moment_copy_ctx[i] = (struct five_moment_copy_job_ctx) {
          .bidx = i,
          .inp_elc = pdata[i].fdup_elc,
          .inp_ion = pdata[i].fdup_ion,
          .inp_maxwell = pdata[i].fdup_maxwell,
          .out_elc = pdata[i].f_elc[0],
          .out_ion = pdata[i].f_ion[0],
          .out_maxwell = pdata[i].f_maxwell[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_patches; i++) {
        gkyl_job_pool_add_work(job_pool, five_moment_copy_job_func, &five_moment_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_patches; i++) {
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
five_moment_write_sol_patch(const char* fbase, int num_patches, const struct five_moment_patch_data pdata[])
{
  for (int i = 0; i < num_patches; i++) {
    const char *fmt_elc = "%s_elc_p%d.gkyl";
    int sz_elc = snprintf(0, 0, fmt_elc, fbase, i);
    char file_nm_elc[sz_elc + 1];

    const char *fmt_ion = "%s_ion_p%d.gkyl";
    int sz_ion = snprintf(0, 0, fmt_ion, fbase, i);
    char file_nm_ion[sz_ion + 1];

    const char *fmt_field = "%s_field_p%d.gkyl";
    int sz_field = snprintf(0, 0, fmt_field, fbase, i);
    char file_nm_field[sz_field + 1];

    snprintf(file_nm_elc, sizeof file_nm_elc, fmt_elc, fbase, i);
    snprintf(file_nm_ion, sizeof file_nm_ion, fmt_ion, fbase, i);
    snprintf(file_nm_field, sizeof file_nm_field, fmt_field, fbase, i);

    five_moment_patch_data_write(file_nm_elc, file_nm_ion, file_nm_field, &pdata[i]);
  }
}

double
five_moment_max_dt_patch(int num_patches, const struct five_moment_patch_data pdata[])
{
  double dt = DBL_MAX;

  for (int i = 0; i < num_patches; i++) {
    dt = fmin(dt, five_moment_patch_data_max_dt(&pdata[i]));
  }

  return dt;
}
