#include <gkyl_vlasov_priv.h>

// Take time-step using the RK3 method for the explicit advective comps.
// Use the actual timestep used to update 
void
vlasov_update_implicit_coll(gkyl_vlasov_app* app, double dt0)
{
  int ns = app->num_species;  
  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];

  for (int i=0; i<ns; ++i) {
    fin[i] = app->species[i].f;
    fout[i] = app->species[i].f1;
  }

  // compute necessary moments and boundary corrections for collisions
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      vm_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
  }
  
  // implicit BGK contributions
  for (int i=0; i<app->num_species; ++i) {
    app->species[i].bgk.implicit_step = true;
    app->species[i].bgk.dt_implicit = dt0;
    vm_species_rhs_implicit(app, &app->species[i], fin[i], fout[i], dt0);
  }

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dt0), 1.0, fin[i]);
    vm_species_apply_bc(app, &app->species[i], fout[i], app->tcurr);
  }
  
  for (int i=0; i<ns; ++i) {
    gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
  };
}