#include <assert.h>
#include <gkyl_vlasov_priv.h>

void
vm_species_emission_init(struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit,
  int dir, enum gkyl_edge_loc edge, void *ctx)
{
  struct gkyl_bc_emission_ctx *params = ctx;
  emit->params = params;
  emit->num_species = params->num_species;
  emit->edge = edge;
  emit->dir = dir;
  emit->elastic = params->elastic;
  emit->t_bound = params->t_bound;
}

void
vm_species_emission_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_emitting_wall *emit)
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  int bdir = (emit->edge == GKYL_LOWER_EDGE) ? 2*emit->dir : 2*emit->dir+1;

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    ghost[cdim+d] = 0;
  }

  if (emit->edge == GKYL_LOWER_EDGE) {
    emit->write = gkyl_range_is_on_lower_edge(emit->dir, &s->lower_skin[emit->dir], &s->global);
  } else {
    emit->write = gkyl_range_is_on_upper_edge(emit->dir, &s->upper_skin[emit->dir], &s->global);
  }

  emit->emit_grid = &s->bflux.boundary_grid[bdir];
  emit->emit_buff_r = &s->bflux.flux_r[bdir];
  emit->emit_ghost_r = (emit->edge == GKYL_LOWER_EDGE) ? &s->lower_ghost[emit->dir] : &s->upper_ghost[emit->dir];
  emit->emit_skin_r = (emit->edge == GKYL_LOWER_EDGE) ? &s->lower_skin[emit->dir] : &s->upper_skin[emit->dir];
  emit->buffer = s->bc_buffer;
  emit->f_emit = mkarr(app->use_gpu, app->basis.num_basis, emit->emit_buff_r->volume);
  emit->f_emit_host = app->use_gpu? mkarr(false, emit->f_emit->ncomp, emit->f_emit->size)
	                          : gkyl_array_acquire(emit->f_emit);

  struct gkyl_array *proj_buffer = mkarr(false, app->basis.num_basis, emit->emit_buff_r->volume);

  // Initialize elastic component of emission
  if (emit->elastic) {
    emit->elastic_yield = mkarr(app->use_gpu, app->basis.num_basis, emit->emit_buff_r->volume);
    emit->elastic_update = gkyl_bc_emission_elastic_new(emit->params->elastic_model,
      emit->elastic_yield, emit->dir, emit->edge, cdim, vdim, s->info.mass, s->f->ncomp, emit->emit_grid,
      emit->emit_buff_r, app->poly_order, app->basis_on_dev.basis, &app->basis, proj_buffer,
      app->use_gpu);
  }

  // Initialize inelastic emission spectrums
  for (int i=0; i<emit->num_species; ++i) {
    emit->impact_species[i] = vm_find_species(app, emit->params->in_species[i]);
    emit->impact_grid[i] = &emit->impact_species[i]->bflux.boundary_grid[bdir];

    emit->flux_slvr[i] = gkyl_dg_updater_moment_new(emit->impact_grid[i], &app->confBasis,
      &app->basis, NULL, NULL, NULL, emit->impact_species[i]->model_id, 0, GKYL_F_MOMENT_M0M1M2, true,
      app->use_gpu);

    emit->impact_skin_r[i] = (emit->edge == GKYL_LOWER_EDGE) ? &emit->impact_species[i]->lower_skin[emit->dir] : &emit->impact_species[i]->upper_skin[emit->dir];
    emit->impact_ghost_r[i] = (emit->edge == GKYL_LOWER_EDGE) ? &emit->impact_species[i]->lower_ghost[emit->dir] : &emit->impact_species[i]->upper_ghost[emit->dir];
    emit->impact_buff_r[i] = &emit->impact_species[i]->bflux.flux_r[bdir];
    emit->impact_cbuff_r[i] = &emit->impact_species[i]->bflux.conf_r[bdir];

    emit->yield[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->impact_buff_r[i]->volume);
    emit->spectrum[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->emit_buff_r->volume);
    emit->weight[i] = mkarr(app->use_gpu, app->confBasis.num_basis,
      emit->impact_cbuff_r[i]->volume);
    emit->flux[i] = mkarr(app->use_gpu, app->confBasis.num_basis, emit->impact_cbuff_r[i]->volume);
    emit->bflux_arr[i] = emit->impact_species[i]->bflux.flux_arr[bdir];
    emit->k[i] = mkarr(app->use_gpu, app->confBasis.num_basis, emit->impact_cbuff_r[i]->volume);

    gkyl_bc_emission_flux_ranges(&emit->impact_normal_r[i], emit->dir + cdim,
      emit->impact_buff_r[i], ghost, emit->edge);
    
    emit->update[i] = gkyl_bc_emission_spectrum_new(emit->params->spectrum_model[i],
      emit->params->yield_model[i], emit->yield[i], emit->spectrum[i], emit->dir, emit->edge,
      cdim, vdim, emit->impact_species[i]->info.mass, s->info.mass, emit->impact_buff_r[i],
      emit->emit_buff_r, emit->impact_grid[i], emit->emit_grid, app->poly_order,
      &app->basis, proj_buffer, app->use_gpu);
  }
  gkyl_array_release(proj_buffer);
}

void
vm_species_emission_apply_bc(struct gkyl_vlasov_app *app, const struct vm_emitting_wall *emit,
  struct gkyl_array *fout, double tcurr)
{
  // Optional scaling of emission with time
  double t_scale = 1.0;
  if (tcurr < emit->t_bound)
    t_scale = sin(M_PI*tcurr/(2.0*emit->t_bound));

  gkyl_array_clear(emit->f_emit, 0.0); // Zero emitted distribution before beginning accumulate

  // Elastic emission contribution
  if (emit->elastic) {
    gkyl_bc_emission_elastic_advance(emit->elastic_update, emit->emit_skin_r, emit->buffer, fout,
      emit->f_emit, emit->elastic_yield, &app->basis);
  }
  // Inelastic emission contribution
  for (int i=0; i<emit->num_species; ++i) {
    int species_idx;
    species_idx = vm_find_species_idx(app, emit->impact_species[i]->info.name);
    gkyl_dg_updater_moment_advance(emit->flux_slvr[i], &emit->impact_normal_r[i],
      emit->impact_cbuff_r[i], emit->bflux_arr[i], emit->flux[i]);
    
    gkyl_bc_emission_spectrum_advance(emit->update[i], emit->impact_buff_r[i],
      emit->impact_cbuff_r[i], emit->emit_buff_r, emit->bflux_arr[i],
      emit->f_emit, emit->yield[i], emit->spectrum[i], emit->weight[i], emit->flux[i],
      emit->k[i]);
  }
  gkyl_array_set_range_to_range(fout, t_scale, emit->f_emit, emit->emit_ghost_r,
    emit->emit_buff_r);
}


// KB - The write function only works in 1x at the moment.
// It expects a single rank to own the whole emit range.
void
vm_species_emission_write(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_emitting_wall *emit, struct gkyl_msgpack_data *mt, int frame)
{
  const char *fmt = (emit->edge == GKYL_LOWER_EDGE) ? "%s-%s_bc_lo_%d.gkyl" : "%s-%s_bc_up_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  if (emit->write) {
    if (app->use_gpu)
      gkyl_array_copy(emit->f_emit_host, emit->f_emit);
    gkyl_grid_sub_array_write(emit->emit_grid, emit->emit_buff_r, mt, emit->f_emit_host, fileNm);
  }
}

void
vm_species_emission_release(const struct vm_emitting_wall *emit)
{
  gkyl_array_release(emit->f_emit_host);
  gkyl_array_release(emit->f_emit);
  if (emit->elastic) {
    gkyl_array_release(emit->elastic_yield);
    gkyl_bc_emission_elastic_release(emit->elastic_update);
  }
  for (int i=0; i<emit->num_species; ++i) {
    gkyl_array_release(emit->yield[i]);
    gkyl_array_release(emit->spectrum[i]);
    gkyl_array_release(emit->weight[i]);
    gkyl_array_release(emit->flux[i]);
    gkyl_array_release(emit->k[i]);
    gkyl_dg_updater_moment_release(emit->flux_slvr[i]);
    gkyl_bc_emission_spectrum_release(emit->update[i]);
  }
}
