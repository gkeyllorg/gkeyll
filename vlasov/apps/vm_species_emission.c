#include <assert.h>
#include <gkyl_vlasov_priv.h>

void
vm_species_emission_init(
  struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit, int dir, enum gkyl_edge_loc edge,
  void *ctx
)
{
  struct gkyl_bc_emission_ctx *params = ctx;
  emit->params = params;
  emit->num_species = params->num_species;
  emit->edge = edge;
  emit->dir = dir;
  emit->elastic = params->elastic;
  emit->t_bound = params->t_bound;
}

// Multiply a phase-space array on the emission buffer range (the ghost cells of the emitting
// edge) by the configuration-space Jacobian of those cells (the adjacent skin cell's).
static void
vm_species_emission_rescale_jacobpos(
  struct gkyl_vlasov_app *app, const struct vm_species *vms, const struct vm_emitting_wall *emit,
  struct gkyl_array *arr
)
{
  const struct gkyl_vlasov_position_map *vpm = vms->pos_map;
  if (vpm->is_identity) {
    return;
  }
  int cdim = app->cdim;
  struct gkyl_array *arr_host = app->use_gpu ? mkarr(false, arr->ncomp, arr->size) :
                                               gkyl_array_acquire(arr);
  if (app->use_gpu) {
    gkyl_array_copy(arr_host, arr);
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, emit->emit_buff_r);
  int idx_conf[GKYL_MAX_CDIM];
  while (gkyl_range_iter_next(&iter)) {
    for (int d = 0; d < cdim; ++d) {
      idx_conf[d] = iter.idx[d];
    }
    const double *jacob =
      gkyl_array_cfetch(vpm->jacob_pos_gauss_host, gkyl_range_idx(&vpm->local_ext_pos, idx_conf));
    double *a = gkyl_array_fetch(arr_host, gkyl_range_idx(emit->emit_buff_r, iter.idx));
    for (int k = 0; k < arr->ncomp; ++k) {
      a[k] *= jacob[0];
    }
  }
  if (app->use_gpu) {
    gkyl_array_copy(arr, arr_host);
  }
  gkyl_array_release(arr_host);
}

void
vm_species_emission_cross_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_emitting_wall *emit
)
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  int bdir = (emit->edge == GKYL_LOWER_EDGE) ? 2 * emit->dir : 2 * emit->dir + 1;

  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d = 0; d < vdim; ++d) {
    ghost[cdim + d] = 0;
  }

  if (emit->edge == GKYL_LOWER_EDGE) {
    emit->write = gkyl_range_is_on_lower_edge(emit->dir, &vms->lower_skin[emit->dir], &vms->global);
  } else {
    emit->write = gkyl_range_is_on_upper_edge(emit->dir, &vms->upper_skin[emit->dir], &vms->global);
  }

  emit->emit_grid = &vms->bflux.boundary_grid[bdir];
  emit->emit_buff_r = &vms->bflux.flux_r[bdir];
  emit->emit_ghost_r = (emit->edge == GKYL_LOWER_EDGE) ? &vms->lower_ghost[emit->dir] :
                                                         &vms->upper_ghost[emit->dir];
  emit->emit_skin_r = (emit->edge == GKYL_LOWER_EDGE) ? &vms->lower_skin[emit->dir] :
                                                        &vms->upper_skin[emit->dir];
  emit->buffer = vms->bc_buffer;
  emit->f_emit = mkarr(app->use_gpu, vms->basis.num_basis, emit->emit_buff_r->volume);
  emit->f_emit_host = app->use_gpu ? mkarr(false, emit->f_emit->ncomp, emit->f_emit->size) :
                                     gkyl_array_acquire(emit->f_emit);

  struct gkyl_array *proj_buffer = mkarr(false, vms->basis.num_basis, emit->emit_buff_r->volume);

  // The emission models are functions of the physical velocity, so the projections and yield
  // evaluations go through the species' coordinate maps (identity on uniform grids).
  emit->c2p_emit =
    (struct vm_proj_c2p_ctx){.cdim = cdim, .pos_map = vms->pos_map, .vel_map = vms->vel_map};

  // Initialize elastic component of emission
  if (emit->elastic) {
    // The elastic yield multiplies the stored distribution with a phase-space weak multiply;
    // only the serendipity and tensor phase bases have those kernels, so refuse the tensor p=1
    // hybrid basis with an actionable message instead of an assert in gkyl_dg_mul_op.
    if (vms->basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
      gkyl_exit(
        "vm_species_emission: elastic emission is not supported on the tensor p=1 (hybrid) basis; use serendipity or tensor p>=2."
      );
    }
    emit->elastic_yield = mkarr(app->use_gpu, vms->basis.num_basis, emit->emit_buff_r->volume);
    emit->elastic_update = gkyl_bc_emission_elastic_new(
      emit->params->elastic_model, emit->elastic_yield, emit->dir, emit->edge, cdim, vdim,
      vms->mass, vms->f->ncomp, emit->emit_grid, emit->emit_buff_r, app->poly_order,
      vms->basis_on_dev, &vms->basis, proj_buffer, vm_proj_c2p_phase, &emit->c2p_emit, app->use_gpu
    );
  }

  // Initialize inelastic emission spectrums
  for (int i = 0; i < emit->num_species; ++i) {
    emit->impact_species[i] = vm_find_species(app, emit->params->in_species[i]);
    // in_species must name an existing kinetic species.
    assert(emit->impact_species[i]);
    struct vm_species *imp = emit->impact_species[i];
    emit->impact_grid[i] = &imp->bflux.boundary_grid[bdir];

    // Boundary-flux moments of the impact species' distribution: the moment
    // type captures that species' basis/velocity map/Hamiltonian.
    struct gkyl_mom_vlasov_inp inp_mom = {
      .conf_basis = &app->basis,
      .phase_basis = &imp->basis,
      .vel_range = &imp->local_vel,
      .vel_map = imp->vel_map,
      .hamil_range = &imp->hamil_range,
      .hamil = imp->hamil,
      .model_id = imp->model_id,
      .hamil_id = imp->hamil_id,
      .mom_type = GKYL_F_MOMENT_M0M1M2,
      .use_gpu = app->use_gpu,
    };
    emit->mom_type[i] = gkyl_int_mom_vlasov_inew(&inp_mom);

    emit->flux_slvr[i] = gkyl_mom_calc_new(emit->impact_grid[i], emit->mom_type[i], app->use_gpu);

    emit->impact_skin_r[i] = (emit->edge == GKYL_LOWER_EDGE) ?
                               &emit->impact_species[i]->lower_skin[emit->dir] :
                               &emit->impact_species[i]->upper_skin[emit->dir];
    emit->impact_ghost_r[i] = (emit->edge == GKYL_LOWER_EDGE) ?
                                &emit->impact_species[i]->lower_ghost[emit->dir] :
                                &emit->impact_species[i]->upper_ghost[emit->dir];
    emit->impact_buff_r[i] = &emit->impact_species[i]->bflux.flux_r[bdir];
    emit->impact_cbuff_r[i] = &emit->impact_species[i]->bflux.conf_r[bdir];

    emit->yield[i] = mkarr(app->use_gpu, vms->basis.num_basis, emit->impact_buff_r[i]->volume);
    emit->spectrum[i] = mkarr(app->use_gpu, vms->basis.num_basis, emit->emit_buff_r->volume);
    emit->weight[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->impact_cbuff_r[i]->volume);
    emit->flux[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->impact_cbuff_r[i]->volume);
    emit->bflux_arr[i] = emit->impact_species[i]->bflux.flux_arr[bdir];
    emit->k[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->impact_cbuff_r[i]->volume);
    emit->weight_host[i] = app->use_gpu ?
                             mkarr(false, emit->weight[i]->ncomp, emit->weight[i]->size) :
                             gkyl_array_acquire(emit->weight[i]);
    emit->flux_host[i] = app->use_gpu ? mkarr(false, emit->flux[i]->ncomp, emit->flux[i]->size) :
                                        gkyl_array_acquire(emit->flux[i]);

    gkyl_bc_emission_flux_ranges(
      &emit->impact_normal_r[i], emit->dir + cdim, emit->impact_buff_r[i], ghost, emit->edge
    );

    emit->c2p_impact[i] =
      (struct vm_proj_c2p_ctx){.cdim = cdim, .pos_map = imp->pos_map, .vel_map = imp->vel_map};
    emit->update[i] = gkyl_bc_emission_spectrum_new(
      emit->params->spectrum_model[i], emit->params->yield_model[i], emit->yield[i],
      emit->spectrum[i], emit->dir, emit->edge, cdim, vdim, emit->impact_species[i]->mass,
      vms->mass, emit->impact_buff_r[i], emit->emit_buff_r, emit->impact_grid[i], emit->emit_grid,
      app->poly_order, &vms->basis, proj_buffer, vm_proj_c2p_phase, &emit->c2p_impact[i],
      vm_proj_c2p_phase, &emit->c2p_emit, app->use_gpu
    );

    // The solver evolves J_x J_v f, and the ghost cells carry the Jacobians of the adjacent
    // skin cell, so the projected (physical) spectrum is weighted by both Jacobians: the
    // emitted distribution then carries the flux the normalization asks for. (The elastic
    // part reflects the stored J_x J_v f and needs no weighting.)
    gkyl_vlasov_velocity_map_rescale_jacobvel(
      vms->vel_map, &app->basis, &vms->basis, emit->emit_buff_r, emit->spectrum[i],
      emit->spectrum[i]
    );
    vm_species_emission_rescale_jacobpos(app, vms, emit, emit->spectrum[i]);
  }
  gkyl_array_release(proj_buffer);

  emit->yield_diag = gkyl_dynvec_new(GKYL_DOUBLE, 2 * emit->num_species);
  emit->is_first_yield_write_call = true;
}

void
vm_species_emission_apply_bc(
  struct gkyl_vlasov_app *app, const struct vm_species *vms, const struct vm_emitting_wall *emit,
  struct gkyl_array *fout, double tcurr
)
{
  // Optional scaling of emission with time
  double t_scale = 1.0;
  if (tcurr < emit->t_bound) {
    t_scale = sin(M_PI * tcurr / (2.0 * emit->t_bound));
  }

  gkyl_array_clear(emit->f_emit, 0.0); // Zero emitted distribution before beginning accumulate

  // Elastic emission contribution
  if (emit->elastic) {
    gkyl_bc_emission_elastic_advance(
      emit->elastic_update, emit->emit_skin_r, emit->buffer, fout, emit->f_emit,
      emit->elastic_yield, &vms->basis
    );
  }
  // Inelastic emission contribution
  for (int i = 0; i < emit->num_species; ++i) {
    gkyl_mom_calc_advance(
      emit->flux_slvr[i], &emit->impact_normal_r[i], emit->impact_cbuff_r[i], emit->bflux_arr[i],
      emit->flux[i]
    );

    gkyl_bc_emission_spectrum_advance(
      emit->update[i], emit->impact_buff_r[i], emit->impact_cbuff_r[i], emit->emit_buff_r,
      emit->bflux_arr[i], emit->f_emit, emit->yield[i], emit->spectrum[i], emit->weight[i],
      emit->flux[i], emit->k[i]
    );
  }
  gkyl_array_set_range_to_range(fout, t_scale, emit->f_emit, emit->emit_ghost_r, emit->emit_buff_r);
}

// KB - The write function only works in 1x at the moment.
// It expects a single rank to own the whole emit range.
void
vm_species_emission_write(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_emitting_wall *emit,
  struct gkyl_msgpack_data *mt, int frame
)
{
  const char *fmt = (emit->edge == GKYL_LOWER_EDGE) ? "%s-%s_bc_lo_%d.gkyl" : "%s-%s_bc_up_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);

  if (emit->write) {
    if (app->use_gpu) {
      gkyl_array_copy(emit->f_emit_host, emit->f_emit);
    }
    gkyl_grid_sub_array_write(emit->emit_grid, emit->emit_buff_r, mt, emit->f_emit_host, fileNm);

    // Yield diagnostic: for each impact species the effective yield (the flux-weighted mean
    // of the yield over the incoming velocity cells, summed over the boundary cells) and the
    // incoming flux (summed over the boundary cells), as used by the last emission update.
    double vals[2 * GKYL_MAX_SPECIES];
    for (int i = 0; i < emit->num_species; ++i) {
      if (app->use_gpu) {
        gkyl_array_copy(emit->weight_host[i], emit->weight[i]);
        gkyl_array_copy(emit->flux_host[i], emit->flux[i]);
      }
      double w0 = 0.0, w1 = 0.0, flux = 0.0;
      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, emit->impact_cbuff_r[i]);
      while (gkyl_range_iter_next(&iter)) {
        long loc = gkyl_range_idx(emit->impact_cbuff_r[i], iter.idx);
        const double *w = gkyl_array_cfetch(emit->weight_host[i], loc);
        const double *fl = gkyl_array_cfetch(emit->flux_host[i], loc);
        w0 += w[0];
        w1 += w[1];
        flux += fl[0];
      }
      vals[2 * i] = w1 != 0.0 ? w0 / w1 : 0.0;
      vals[2 * i + 1] = flux;
    }
    gkyl_dynvec_append(emit->yield_diag, app->tcurr, vals);

    const char *fmt_yield = (emit->edge == GKYL_LOWER_EDGE) ? "%s-%s_bc_lo_yield.gkyl" :
                                                              "%s-%s_bc_up_yield.gkyl";
    int sz_yield = gkyl_calc_strlen(fmt_yield, app->name, vms->name);
    char fileNm_yield[sz_yield + 1]; // ensures no buffer overflow
    snprintf(fileNm_yield, sizeof fileNm_yield, fmt_yield, app->name, vms->name);
    if (emit->is_first_yield_write_call) {
      gkyl_dynvec_write(emit->yield_diag, fileNm_yield);
      emit->is_first_yield_write_call = false;
    } else {
      gkyl_dynvec_awrite(emit->yield_diag, fileNm_yield);
    }
    gkyl_dynvec_clear(emit->yield_diag);
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
  for (int i = 0; i < emit->num_species; ++i) {
    gkyl_array_release(emit->yield[i]);
    gkyl_array_release(emit->spectrum[i]);
    gkyl_array_release(emit->weight[i]);
    gkyl_array_release(emit->flux[i]);
    gkyl_array_release(emit->k[i]);
    gkyl_array_release(emit->weight_host[i]);
    gkyl_array_release(emit->flux_host[i]);
    gkyl_mom_type_release(emit->mom_type[i]);
    gkyl_mom_calc_release(emit->flux_slvr[i]);
    gkyl_bc_emission_spectrum_release(emit->update[i]);
  }
  gkyl_dynvec_release(emit->yield_diag);
}
