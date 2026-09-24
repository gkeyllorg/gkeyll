#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_radiation_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_rad_drag *rad)
{
  int vdim = app->vdim;

  // allocate nu and initialize it
  rad->nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *nu_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  
  gkyl_proj_on_basis *proj_nu = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.radiation.nu, s->info.radiation.ctx_nu);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &app->local, nu_host);
  gkyl_proj_on_basis_release(proj_nu);
  gkyl_array_copy(rad->nu, nu_host);
  gkyl_array_release(nu_host);

  // Using the LBO operator, which expects a vdim+1 primitive moments array
  // Last array entry is T/m which is checked for positivity before applying drag operator. 
  rad->nu_rad_drag = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *nu_rad_drag_host = mkarr(false, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  
  gkyl_proj_on_basis *proj_rad_drag = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, vdim+1, s->info.radiation.nu_rad_drag, s->info.radiation.ctx_nu_rad_drag);
  gkyl_proj_on_basis_advance(proj_rad_drag, 0.0, &app->local, nu_rad_drag_host);
  gkyl_proj_on_basis_release(proj_rad_drag);
  gkyl_array_copy(rad->nu_rad_drag, nu_rad_drag_host);
  gkyl_array_release(nu_rad_drag_host);

  // Radiation operator uses the LBO kernels, so create auxiliary field struct for LBO
  struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = { .nuSum = rad->nu, .nuPrimMomsSum = rad->nu_rad_drag };
  rad->rad_slvr = gkyl_dg_updater_rad_vlasov_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, &drag_inp, app->use_gpu);
}

// updates the radiation terms in the rhs 
void
vm_species_radiation_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_rad_drag *rad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{    
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_rad_vlasov_advance(rad->rad_slvr, &species->local,
    fin, species->cflrate, rhs);
}

void 
vm_species_radiation_release(const struct gkyl_vlasov_app *app, const struct vm_rad_drag *rad)
{
  gkyl_array_release(rad->nu);
  gkyl_array_release(rad->nu_rad_drag);

  gkyl_dg_updater_rad_vlasov_release(rad->rad_slvr);
 }
