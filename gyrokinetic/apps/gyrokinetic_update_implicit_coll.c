#include <gkyl_gyrokinetic_priv.h>

// Take time-step using an implicit method for collisions.
// Use the actual timestep used to update explicit advection. 
void
gyrokinetic_update_implicit_coll(gkyl_gyrokinetic_app* app, double dt0)
{
  // Add contribution from implicit collision terms.
  double tcurr = app->tcurr, dt = dt0;

  struct gkyl_update_status st = { .success = true };

  const int ns = app->num_species;  
  struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  struct gkyl_array **bflux_in[ns];
  struct gkyl_array **bflux_out[ns];

  const int neuts = app->num_neut_species;  
  struct gkyl_array *fin_neut[neuts];
  struct gkyl_array *fout_neut[neuts];
  struct gkyl_array **bflux_in_neut[ns];
  struct gkyl_array **bflux_out_neut[ns];

  // Fetch input and output arrays.
  for (int i=0; i<ns; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f;
    fout[i] = gks->f1;
    bflux_in[i] = gks->bflux.f;
    bflux_out[i] = gks->bflux.f;
  }
  for (int i=0; i<neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f;
    fout_neut[i] = gkns->f1;
    bflux_in_neut[i] = gkns->bflux.f;
    bflux_out_neut[i] = gkns->bflux.f;
  }

  // Compute df/dt from implicit terms.
  gyrokinetic_rhs_implicit(app, tcurr, dt, fin, fout, bflux_out, fin_neut, fout_neut, bflux_out_neut, &st);

  // Apply boundary conditions and copy solution.
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);

  for (int i=0; i<ns; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_copy_range(gks, fin[i], fout[i], &gks->local_ext);
  }
  for (int i=0; i<neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_copy_range(gkns, fin_neut[i], fout_neut[i], &gkns->local_ext);
  };
}
