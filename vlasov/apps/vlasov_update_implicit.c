#include <gkyl_vlasov_priv.h>

// The implicit half of the operator-split step, taken with the actual time-step
// of the explicit SSP-RK3 half: implicit BGK collisions for each species, then
// the implicit fluid-EM coupling.

// Implicit BGK collisions in three phases over all species: moments, implicit
// RHS (which computes cross moments, so it needs all species' moments first),
// then BCs + copy-back into the solution.
static void
update_implicit_coll(gkyl_vlasov_app* app, double dt0)
{
  int num_species = app->num_species + app->num_fluid_species;
  const struct gkyl_array *fin[num_species];
  struct gkyl_array *fout[num_species];

  for (int i=0; i<num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    fin[i] = sp->dist ? sp->dist->f : 0;
    fout[i] = sp->dist ? sp->dist->f1 : 0;
  }

  for (int i=0; i<num_species; ++i)
    vlasov_species_calc_implicit_moms(app, &app->species[i], fin[i]);

  for (int i=0; i<num_species; ++i)
    vlasov_species_rhs_implicit(app, &app->species[i], fin[i], fout[i], dt0);

  for (int i=0; i<num_species; ++i)
    vlasov_species_finish_implicit_update(app, &app->species[i], fout[i], app->tcurr);
}

void
vlasov_update_implicit(gkyl_vlasov_app* app, double dt0)
{
  if (app->has_implicit_coll_scheme)
    update_implicit_coll(app, dt0);

  if (app->has_fluid_em_coupling)
    vm_fluid_em_coupling_update(app, app->fl_em, app->tcurr, dt0);
}
