#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src)
{
  int vdim = app->vdim;  
  src->calc_bflux = false;
  if (s->source_id == GKYL_BFLUX_SOURCE) {
    src->calc_bflux = true;
    assert(s->info.source.source_length);
    assert(s->info.source.source_species);
    src->source_length = s->info.source.source_length;
    src->source_species = vm_find_species(app, s->info.source.source_species);
    src->source_species_idx = vm_find_species_idx(app, s->info.source.source_species);
    if (app->use_gpu) {
      src->scale_ptr = gkyl_cu_malloc((vdim+2)*sizeof(double));
    }
    else {
      src->scale_ptr = gkyl_malloc((vdim+2)*sizeof(double));
    }
  }

  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  src->source_host = src->source;
  if (app->use_gpu) {
    src->source_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
  }
  src->scale_factor = 1.0;
  
  src->write_source = s->info.source.write_source; // optional flag to write out source

  src->num_sources = s->info.source.num_sources;
  for (int k=0; k<s->info.source.num_sources; k++) {
    vm_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);
  }

  // Allocate temporary variable for accumulating multiple sources if multiple sources present
  if (src->num_sources > 1) {
    src->source_tmp = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  }

  // Allocate data and updaters for diagnostic moments.
  src->num_diag_moments = s->info.num_diag_moments;
  s->src.moms = gkyl_malloc(sizeof(struct vm_species_moment[src->num_diag_moments]));
  for (int m=0; m<src->num_diag_moments; ++m) {
    vm_species_moment_init(app, s, &s->src.moms[m], s->info.diag_moments[m], false);
  }

  // Allocate data and updaters for integrated moments.
  vm_species_moment_init(app, s, &s->src.integ_moms, GKYL_F_MOMENT_M0M1M2, true);
  if (app->use_gpu) {
    s->src.red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
    s->src.red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } 
  else {
    s->src.red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
    s->src.red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  s->src.integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->src.is_first_integ_write_call = true;
}

void
vm_species_source_calc(gkyl_vlasov_app *app, struct vm_species *s, 
  struct vm_source *src, double tm)
{
  if (s->source_id) {
    if (src->num_sources > 1) {
      gkyl_array_clear(src->source, 0.0);
      for (int k=0; k<src->num_sources; k++) {
        vm_species_projection_calc(app, s, &src->proj_source[k], src->source_tmp, tm);
        gkyl_array_accumulate(src->source, 1.0, src->source_tmp);
      }
    }
    else {
      vm_species_projection_calc(app, s, &src->proj_source[0], src->source, tm);
    }
  }
}

// computes rhs of the boundary flux
void
vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  int species_idx;
  species_idx = vm_find_species_idx(app, species->info.name);
  // use boundary fluxes to scale source profile
  if (src->calc_bflux) {
    src->scale_factor = 0.0;
    double z[app->confBasis.num_basis];
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
vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src)
{
  gkyl_array_release(src->source);

  if (app->use_gpu) {
    gkyl_array_release(src->source_host);
  }

  if (src->calc_bflux) {
    if (app->use_gpu) {
      gkyl_cu_free(src->scale_ptr);
    } 
    else {
      gkyl_free(src->scale_ptr);
    }
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
    gkyl_cu_free(src->red_integ_diag_global);
  }
  else {
    gkyl_free(src->red_integ_diag);
    gkyl_free(src->red_integ_diag_global);
  }  
  gkyl_dynvec_release(src->integ_diag);  
}
