#include <gkyl_gyrokinetic_priv.h>

// Take time-step using an implicit method for collisions.
// Use the actual timestep used to update explicit advection.
void
gyrokinetic_update_implicit_coll(gkyl_gyrokinetic_app *app, double dt0)
{
  // Add contribution from implicit collision terms.
  double tcurr = app->tcurr, dt = dt0;

  struct gkyl_update_status st = {.success = true};

  const int ns = app->num_species;
  const int neuts = app->num_neut_species;
  struct gkyl_gyrokinetic_fdot_args *fdot_args = &app->fdot_args;
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***smom_in = fdot_args->smom_in;
  struct gkyl_array ***smom_in_neut = fdot_args->smom_in_neut;
  struct gkyl_array ***smom_out = fdot_args->smom_out;
  struct gkyl_array ***smom_out_neut = fdot_args->smom_out_neut;

  // Fetch input and output arrays.
  for (int i = 0; i < ns; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f;
    fout[i] = gks->f1;
    smom_in[i] = gks->step_moms.f;
    smom_out[i] = gks->step_moms.f;
  }
  for (int i = 0; i < neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f;
    fout_neut[i] = gkns->f1;
    smom_in_neut[i] = gkns->step_moms.f;
    smom_out_neut[i] = gkns->step_moms.f;
  }

  // Compute df/dt from implicit terms.
  gyrokinetic_rhs_implicit(app, tcurr, dt, fdot_args, &st);

  // Apply boundary conditions and copy solution.
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, smom_out, fout_neut);

  for (int i = 0; i < ns; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_copy_range(gks, gks->f, fout[i], &gks->local_ext);

    // Rescale species (e.g. fixed fraction, adiabatic).
    gk_species_scaling_apply(app, gks, &gks->sca, gks->f, smom_out);
  }
  for (int i = 0; i < neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_copy_range(gkns, gkns->f, fout_neut[i], &gkns->local_ext);
  };
}
