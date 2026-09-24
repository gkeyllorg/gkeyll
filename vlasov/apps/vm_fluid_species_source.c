#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_fluid_species_source_init(struct gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct vm_fluid_source *src)
{
  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  
  src->source_host = src->source;
  if (app->use_gpu)
    src->source_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

    src->source_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp)    {
      .grid = &app->grid,
      .basis = &app->confBasis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = app->basis.poly_order+1,
      .num_ret_vals = 1,
      .eval = fluid_species->info.source.profile,
      .ctx = fluid_species->info.source.ctx
    }
  );
}

void
vm_fluid_species_source_calc(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm)
{
  if (fluid_species->source_id) {
    gkyl_proj_on_basis_advance(fluid_species->src.source_proj, tm, &app->local_ext, fluid_species->src.source_host);
    if (app->use_gpu) // note: source_host is same as source when not on GPUs
      gkyl_array_copy(fluid_species->src.source, fluid_species->src.source_host);
  }
}

// computes rhs of the boundary flux
void
vm_fluid_species_source_rhs(gkyl_vlasov_app *app, const struct vm_fluid_species *species,
  struct vm_fluid_source *src, const struct gkyl_array *fluid[], struct gkyl_array *rhs[])
{
  int species_idx;
  species_idx = vm_find_fluid_species_idx(app, species->info.name);
  
  gkyl_array_accumulate(rhs[species_idx], 1.0, src->source);
}

void
vm_fluid_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_fluid_source *src)
{
  gkyl_array_release(src->source);
  if (app->use_gpu) {
    gkyl_array_release(src->source_host);
  }
  gkyl_proj_on_basis_release(src->source_proj);
}
