#include <gkyl_moment_priv.h>

// Compute the nc intergated values over the update_rgn, storing the
// result in the integ_q
void
calc_integ_quant(const struct gkyl_wv_eqn *eqn, double vol, const struct gkyl_array *q, const struct gkyl_wave_geom *geom,
  struct gkyl_range update_rng, double *integ_q)
{
  int nc = eqn->num_diag;
  double integ_out[nc];
  for (int i=0; i<nc; ++i) integ_q[i] = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_rng);
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(geom, iter.idx);
    const double *qcell = gkyl_array_cfetch(q, gkyl_range_idx(&update_rng, iter.idx));

    eqn->cons_to_diag(eqn, qcell, integ_out);
    for (int i=0; i<nc; ++i) integ_q[i] += vol*cg->kappa*integ_out[i];
  }
}

// Check if nan occurs in the array, returning true if they do and
// false otherwise
bool
check_for_nans(const struct gkyl_array *q, struct gkyl_range update_rng)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_rng);
  while (gkyl_range_iter_next(&iter)) {
    const double *qcell = gkyl_array_cfetch(q, gkyl_range_idx(&update_rng, iter.idx));
    for (int i=0; i<q->ncomp; ++i)
      if (isnan(qcell[i])) return true;
  }
  return false;
}

void
moment_apply_periodic_corner_sync_2d(const gkyl_moment_app *app, struct gkyl_array *f)
{
  long idx_src, idx_dest;
  double *out;
  const double *inp;
  if (app->ndim == 2) {
    // LL skin cell -> UU ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.lower[0],app->local.lower[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.upper[0]+1,app->local.upper[1]+1});

    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);

    // LU skin cell -> UL ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.lower[0],app->local.upper[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.upper[0]+1,app->local.lower[1]-1});

    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);

    // UL skin cell -> LU ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.upper[0],app->local.lower[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.lower[0]-1,app->local.upper[1]+1});

    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);

    // UU skin cell -> LL ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.upper[0],app->local.upper[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.lower[0]-1,app->local.lower[1]-1});

    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);
  }
}

// apply wedge BCs
void
moment_apply_wedge_bc(const gkyl_moment_app *app, double tcurr,
  const struct gkyl_range *update_rng, struct gkyl_array *bc_buffer,
  int dir, const struct gkyl_wv_apply_bc *lo, const struct gkyl_wv_apply_bc *up,
  struct gkyl_array *f)
{
  gkyl_wv_apply_bc_to_buff(lo, tcurr, update_rng, f, bc_buffer->data);
  gkyl_array_copy_from_buffer(f, bc_buffer->data, &(app->skin_ghost.upper_ghost[dir]));

  gkyl_wv_apply_bc_to_buff(up, tcurr, update_rng, f, bc_buffer->data);
  gkyl_array_copy_from_buffer(f, bc_buffer->data, &(app->skin_ghost.lower_ghost[dir]));
}

