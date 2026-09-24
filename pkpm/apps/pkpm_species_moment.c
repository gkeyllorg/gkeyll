#include <assert.h>
#include <gkyl_pkpm_priv.h>

// initialize species moment object
void
pkpm_species_moment_init(struct gkyl_pkpm_app *app, struct pkpm_species *s,
  struct pkpm_species_moment *sm, bool is_diag)
{
  sm->mcalc = gkyl_dg_updater_moment_pkpm_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local, &s->local_vel, s->info.mass, is_diag, app->use_gpu);    
  int num_mom = gkyl_dg_updater_moment_pkpm_num_mom(sm->mcalc);

  sm->marr = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
  sm->marr_host = sm->marr;
  if (app->use_gpu)
    sm->marr_host = mkarr(false, num_mom*app->confBasis.num_basis, app->local_ext.volume);
}

void
pkpm_species_moment_calc(const struct pkpm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_dg_updater_moment_pkpm_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

// release memory for moment data object
void
pkpm_species_moment_release(const struct gkyl_pkpm_app *app, const struct pkpm_species_moment *sm)
{
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_dg_updater_moment_pkpm_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}
