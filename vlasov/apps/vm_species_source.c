#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source *src)
{
  int vdim = app->vdim; 
  src->write_source = vms->info.source.write_source; // Optional flag to write out source and source moments. 
  src->evolve_source = vms->info.source.evolve_source; // Are the sources time-dependent?

  src->rescale_m0 = false; 
  src->calc_bflux = false;
  src->scale_factor = 1.0;
  if (vms->source_id == GKYL_BFLUX_SOURCE) {
    src->calc_bflux = true;
    assert(vms->info.source.source_length);
    assert(vms->info.source.source_species);
    src->source_length = vms->info.source.source_length;
    src->source_species = vm_find_species(app, vms->info.source.source_species);
    src->source_species_idx = vm_find_species_idx(app, vms->info.source.source_species);
    if (app->use_gpu) {
      src->scale_ptr = gkyl_cu_malloc((vdim+2)*sizeof(double));
    }
    else {
      src->scale_ptr = gkyl_malloc((vdim+2)*sizeof(double));
    }
  }
  else if (vms->source_id == GKYL_PROJ_ADAPT_DENSITY_SOURCE) {
    src->rescale_m0 = true; 
    struct gkyl_mom_vlasov_inp inp_mom = {
      .conf_basis = &app->basis,
      .phase_basis = &vms->basis,
      .vel_range = &vms->local_vel,
      .vel_map = vms->vel_map,
      .hamil_range = &vms->mom_hamil_range,
      .hamil = vms->mom_hamil,
      .model_id = vms->model_id,
      .hamil_id = vms->mom_hamil_id,
      .use_gpu = app->use_gpu,
    };
    src->num_cross_source = vms->info.source.num_cross_source; 
    for (int i=0; i<src->num_cross_source; i++) {
      src->adapt_source_species[i] = vm_find_species(app, vms->info.source.source_with[i]);
      src->adapt_source_species_idx[i] = vm_find_species_idx(app, vms->info.source.source_with[i]);
      // Threshold velocity for integration of moments over a subset of the domain. 
      // Also set threshold for whether we accumulate moment over subset of the domain. 
      inp_mom.v_thresh = vms->info.source.source_with_v_thresh[i];
      inp_mom.f_thresh = vms->info.source.source_with_f_thresh[i];
      struct gkyl_mom_type *m0_reduced;
      if (vms->info.source.source_with_upper_half[i]) {
        inp_mom.mom_type = GKYL_F_MOMENT_M0_UPPER; 
        m0_reduced = gkyl_mom_vlasov_inew(&inp_mom);    
      }
      else {
        inp_mom.mom_type = GKYL_F_MOMENT_M0_LOWER; 
        m0_reduced = gkyl_mom_vlasov_inew(&inp_mom);        
      }
      src->m0_reduced[i] = gkyl_mom_calc_new(&vms->grid, m0_reduced, app->use_gpu);
      gkyl_mom_type_release(m0_reduced); 

      src->scale_m0[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      src->scale_m0_host[i] = app->use_gpu ? mkarr(false, app->basis.num_basis, app->local_ext.volume)
                                  : gkyl_array_acquire(src->scale_m0[i]); 
      src->adapt_source[i] = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume); 
      src->adapt_proj_source[i] = vms->info.source.source_with_proj[i];
    }
    src->filter = false;
    if (vms->info.source.filter) {
      // Create Gaussian filter for use with adaptive density source. 
      // Always filter at least once, but optionally filter repeatedly. 
      src->filter = true;
      src->num_filters = vms->info.source.num_filters > 0 ? vms->info.source.num_filters : 1;
      struct gkyl_dg_gaussian_filter_inp inp = {
        .conf_grid = &app->grid,
        .conf_basis = &app->basis,
        .conf_range = &app->local,
        .conf_range_ext = &app->local_ext,
        .extend_filter = false,
        .use_gpu = app->use_gpu,
      };
      src->gauss_filter = gkyl_dg_gaussian_filter_inew(&inp);
    }    
  }

  // We need to ensure source has same shape as distribution function. 
  src->source = mkarr(app->use_gpu, vms->f->ncomp, vms->f->size);
  src->source_host = app->use_gpu ? mkarr(false, src->source->ncomp, src->source->size)
                                  : gkyl_array_acquire(src->source); 

  src->num_sources = vms->info.source.num_sources;
  for (int k=0; k<vms->info.source.num_sources; k++) {
    vm_species_projection_init(app, vms, vms->info.source.projection[k], &src->proj_source[k]);
  }

  // Allocate temporary variable for accumulating multiple sources if multiple sources present.
  if (src->num_sources > 1) {
    src->source_tmp = mkarr(app->use_gpu,  src->source->ncomp, src->source->size);
  }

  // Allocate data and updaters for diagnostic moments.
  src->num_diag_moments = vms->info.num_diag_moments;
  vms->src.moms = gkyl_malloc(sizeof(struct vm_species_moment[src->num_diag_moments]));
  for (int m=0; m<src->num_diag_moments; ++m) {
    vm_species_moment_init(app, vms, &vms->src.moms[m], vms->info.diag_moments[m], false);
  }

  // Allocate data and updaters for integrated moments.
  vm_species_moment_init(app, vms, &vms->src.integ_moms, GKYL_F_MOMENT_M0M1M2, true);
  if (app->use_gpu) {
    vms->src.red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } 
  // allocate dynamic-vector to store all-reduced integrated moments 
  vms->src.integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  vms->src.is_first_integ_write_call = true;
}

void
vm_species_source_calc(gkyl_vlasov_app *app, const struct vm_species *vms, 
  struct vm_source *src, double tm)
{
  if (vms->source_id) {
    if (src->rescale_m0) {
      for (int i=0; i<src->num_cross_source; i++) {
        vm_species_projection_calc(app, vms, &src->proj_source[src->adapt_proj_source[i]], 
          src->adapt_source[i], tm);
      } 
    }
    else if (src->num_sources > 1) {
      gkyl_array_clear(src->source, 0.0);
      for (int k=0; k<src->num_sources; k++) {
        vm_species_projection_calc(app, vms, &src->proj_source[k], src->source_tmp, tm);
        gkyl_array_accumulate(src->source, 1.0, src->source_tmp);
      }
    }
    else {
      vm_species_projection_calc(app, vms, &src->proj_source[0], src->source, tm);
    }
  }
}

void
vm_species_source_adapt_moms(gkyl_vlasov_app *app, const struct vm_species *vms, 
  struct vm_source *src, const struct gkyl_array *fin)
{
  if (vms->source_id == GKYL_PROJ_ADAPT_DENSITY_SOURCE) {  
    for (int i=0; i<src->num_cross_source; i++) {
      gkyl_mom_calc_advance(src->m0_reduced[i], &vms->local, &app->local, fin, src->scale_m0[i]);
      if (src->filter) {
        // Optionally filter repeatedly.
        for (int j=0; j<src->num_filters; j++) {
          // Synchronize ghost cells before filtering so ghost cells are included in filter. 
          gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, src->scale_m0[i]);
          gkyl_dg_gaussian_filter_advance(src->gauss_filter, &app->local, src->scale_m0[i]);
        }
      }
    }  
  }
}

void
vm_species_source_adapt(gkyl_vlasov_app *app, const struct vm_species *vms, 
  struct vm_source *src)
{
  int species_idx;
  species_idx = vm_find_species_idx(app, vms->info.name);  
  if (vms->source_id == GKYL_PROJ_ADAPT_DENSITY_SOURCE) {
    gkyl_array_clear(src->source, 0.0);
    for (int i=0; i<src->num_cross_source; i++) {
      // First compute the adaptive source from self-sourcing. 
      gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &vms->basis, src->source, 
        1.0, src->scale_m0[i], src->adapt_source[i], &app->local, &vms->local); 

      // Next compute the adaptive source from the cross species. 
      gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &vms->basis, src->source, 
        1.0, src->adapt_source_species[i]->src.scale_m0[i], src->adapt_source[i], &app->local, &vms->local); 
    } 
  }
}

// computes rhs of the boundary flux
void
vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  int species_idx;
  species_idx = vm_find_species_idx(app, vms->info.name);
  // use boundary fluxes to scale source profile
  if (src->calc_bflux) {
    src->scale_factor = 0.0;
    double z[app->basis.num_basis];
    double red_mom[1] = { 0.0 };

    for (int d=0; d<app->cdim; ++d) {
      gkyl_array_reduce(src->scale_ptr, src->source_species->bflux.mom_arr[2*d], GKYL_SUM);
      if (app->use_gpu) {
        gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
      }
      else {
        red_mom[0] = src->scale_ptr[0];
      }
      double red_mom_global[1] = { 0.0 };
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, red_mom, red_mom_global);
      src->scale_factor += red_mom_global[0];
      gkyl_array_reduce(src->scale_ptr, src->source_species->bflux.mom_arr[2*d+1], GKYL_SUM);
      if (app->use_gpu) {
        gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
      }
      else {
        red_mom[0] = src->scale_ptr[0];
      }
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, red_mom, red_mom_global);
      src->scale_factor += red_mom_global[0];
    }
    src->scale_factor = src->scale_factor/src->source_length;
  }
  gkyl_array_accumulate(rhs[species_idx], src->scale_factor, src->source);
}

void
vm_species_source_calc_integrated_mom(gkyl_vlasov_app* app, const struct vm_species *vms, 
  struct vm_source *src, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  int vdim = app->vdim;
  double avals[2+vdim], avals_global[2+vdim];
  
  vm_species_moment_calc(&src->integ_moms, vms->local, app->local, src->source);
  app->stat.n_mom += 1;

  // reduce to compute sum over whole domain, append to diagnostics
  if (app->use_gpu) {
    gkyl_array_reduce_range(src->red_integ_diag, src->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals, src->red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    gkyl_array_reduce_range(avals, src->integ_moms.marr_host, GKYL_SUM, &app->local);
  }

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
  gkyl_dynvec_append(src->integ_diag, tm, avals_global);  

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_species_source_write(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = vms->basis.id
    }
  );
  const char *fmt = "%s-%s_source_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, frame);
  
  // Calculate adaptive source before I/O if source is adaptive. 
  vm_species_source_adapt(app, vms, src); 

  // Divide out the velocity space Jacobian from source distribution if present
  // We do the division before I/O to increase the accuracy since we know
  // the velocity-space Jacobian at specific quadrature points. 
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, src->source, vms->f_no_J);
  // Also divide out the configuration-space Jacobian for physical output.
  gkyl_vlasov_position_map_divide_jacobpos(vms->pos_map, &vms->basis,
    &vms->local, vms->f_no_J, vms->f_no_J);

  // If we are on device, copy the source distribution function without the velocity-space
  // Jacobian to the host, otherwise just write out the f_no_J array. 
  if (app->use_gpu) {
    gkyl_array_copy(src->source_host, vms->f_no_J);
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, src->source_host, fileNm);
  }
  else {
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, vms->f_no_J, fileNm);
  }
    
  vlasov_array_meta_release(mt);  

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

void
vm_species_source_write_mom(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  for (int m=0; m<src->num_diag_moments; ++m) {
    struct timespec wtm = gkyl_wall_clock();
    vm_species_moment_calc(&src->moms[m], vms->local, app->local, src->source);
    app->stat.n_mom += 1;
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wtm);
      
    struct timespec wst = gkyl_wall_clock();
    if (app->use_gpu) {
      gkyl_array_copy(src->moms[m].marr_host, src->moms[m].marr);
    }

    const char *fmt = "%s-%s_source_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name,
      gkyl_distribution_moments_strs[vms->info.diag_moments[m]], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name,
      gkyl_distribution_moments_strs[vms->info.diag_moments[m]], frame);
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      src->moms[m].marr_host, fileNm);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;
  }

  // Write out the adaptive density 
  if (vms->source_id == GKYL_PROJ_ADAPT_DENSITY_SOURCE) {
    // Calculate adaptive source density before I/O if source is adaptive. 
    vm_species_source_adapt_moms(app, vms, src, vms->f); 

    const char *fmt_source_M0 = "%s-%s_source_M0_adapt_%d.gkyl";
    int sz_source_M0 = gkyl_calc_strlen(fmt_source_M0, app->name, vms->info.name, frame);
    char fileNm_source_M0[sz_source_M0+1]; // ensures no buffer overflow
    snprintf(fileNm_source_M0, sizeof fileNm_source_M0, fmt_source_M0, app->name, vms->info.name, frame);
    if (app->use_gpu) {
      gkyl_array_copy(src->scale_m0_host[0], src->scale_m0[0]);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->scale_m0_host[0], fileNm_source_M0);
    }
    else {
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->scale_m0[0], fileNm_source_M0);
    }
  }  

  vlasov_array_meta_release(mt); 

  app->stat.n_diag += 1;
}

void
vm_species_source_write_integrated_mom(gkyl_vlasov_app* app, 
  const struct vm_species *vms, struct vm_source *src)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s-source-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, "imom");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, "imom");
    
    if (src->is_first_integ_write_call) {
      gkyl_dynvec_write(src->integ_diag, fileNm);
      src->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(src->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(src->integ_diag);
  app->stat.n_diag_io += 1;
  
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src)
{
  gkyl_array_release(src->source);
  gkyl_array_release(src->source_host);

  if (src->calc_bflux) {
    if (app->use_gpu) {
      gkyl_cu_free(src->scale_ptr);
    } 
    else {
      gkyl_free(src->scale_ptr);
    }
  }

  if (src->rescale_m0) {
    for (int i=0; i<src->num_cross_source; i++) {
      gkyl_mom_calc_release(src->m0_reduced[i]);
      gkyl_array_release(src->scale_m0[i]);
      gkyl_array_release(src->scale_m0_host[i]);
      gkyl_array_release(src->adapt_source[i]);
    }
    gkyl_dg_gaussian_filter_release(src->gauss_filter);
  }

  for (int k=0; k<src->num_sources; k++) {
    vm_species_projection_release(app, &src->proj_source[k]);
  }

  if (src->num_sources > 1) {
    gkyl_array_release(src->source_tmp);
  }

  // Release moment data.
  for (int i=0; i<src->num_diag_moments; ++i) {
    vm_species_moment_release(app, &src->moms[i]);
  }
  gkyl_free(src->moms);
  vm_species_moment_release(app, &src->integ_moms); 
  if (app->use_gpu) {
    gkyl_cu_free(src->red_integ_diag);
  }
  gkyl_dynvec_release(src->integ_diag);  
}
