// Functions required for coupling with FLASH
//
#include <gkyl_dg_bin_ops.h>
#include <gkyl_flash_multib_priv.h>
#include <gkyl_flash_priv.h>


void
gkyl_flash_app_set_fluid_vel_from_fv(gkyl_flash_app *app, struct gkyl_array *fluid_vel_in)
{
  // Each species uses the same bulk velocity
  for (int i=0; i<app->num_species; ++i) {
    struct flash_species *s = &app->species[i];

    gkyl_flash_interp_vars_vel_to_dg(app->interp_vars, 
      &app->local_kinetic_ext, &app->local_fluid_ext, fluid_vel_in, s->fluid_rhou);

    // Ensure pkpm_moms are computed before scaling velocity by density
    flash_species_calc_pkpm_moms(app, &app->species[i], app->species[i].f);

    // Multiply by density because the PKPM updater looks for momentum density instead of velocity
    for (int d=0; d<3; ++d)
      gkyl_dg_mul_op(&app->confBasis, d, s->fluid_rhou, d, s->fluid_rhou, 0, s->pkpm_moms);
  }
}

void
gkyl_flash_multib_app_set_fluid_vel_from_fv(gkyl_flash_multib_app *app, struct gkyl_array *fluid_vel_in[])
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_flash_app_set_fluid_vel_from_fv(app->block_apps[b], fluid_vel_in[b]);

  // Sync momentum density across connected blocks once all blocks are set.
  flash_multib_sync_fluid_rhou(app);
}

void
gkyl_flash_app_set_init_dist_moms_from_fv(gkyl_flash_app *app, struct gkyl_array *init_dist_moms)
{
  for (int i=0; i<app->num_species; ++i) {
    if (!app->species[i].init_dist_moms)
      gkyl_exit("gkyl_flash_app_set_init_dist_moms: species must be created with init_dist_from_arrays = true.");

    gkyl_flash_interp_vars_init_dist_moms_to_dg(app->interp_vars,
      &app->local_kinetic_ext, &app->local_fluid_ext, init_dist_moms, app->species[i].init_dist_moms);
  }
}

void
gkyl_flash_app_init_maxwellian_from_moms_from_fv(gkyl_flash_app *app, struct gkyl_array *init_dist_moms)
{
  gkyl_flash_app_set_init_dist_moms_from_fv(app, init_dist_moms);
  flash_app_init_maxwellian_from_stored_moms(app);
}

void
gkyl_flash_app_get_kinetic_fluxes_to_fv(gkyl_flash_app *app, struct gkyl_array *coupling_fluxes)
{
  for (int i=0; i<app->num_species; ++i) {
    // Verify that array exists
    if (coupling_fluxes->elemsz == 0) {
      char msg[128];
      snprintf(msg, 128, 
        "Output coupling flux array is empty. Check that the proper number of arrays have been initialized. Exiting...");
      gkyl_exit(msg);
    }

    // Compute PKPM moments, [rho, p_par, p_perp].
    flash_species_calc_pkpm_moms(app, &app->species[i], app->species[i].f);

    // Compute coupling moments, self-consistent pressure
    // force (div(p_par b_hat)), and primitive variables
    flash_species_calc_pkpm_vars(app, &app->species[i], app->species[i].f);

    // Compute remaining PKPM variables: acceleration variables in kinetic equation and 
    // source distribution functions for Laguerre couplings.
    flash_species_calc_pkpm_update_vars(app, &app->species[i], app->species[i].f);

    // Compute fluxes
    flash_species_calc_coupling_fluxes(app, &app->species[i], app->species[i].f);

    // Interpolate to FV
    gkyl_flash_interp_vars_coupling_fluxes_to_fv(app->interp_vars, 
      &app->local_kinetic_ext, &app->local_fluid_ext, 
      app->species[i].flash_coupling_fluxes, coupling_fluxes);
  }
}

void
gkyl_flash_multib_app_get_kinetic_fluxes_to_fv(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *coupling_fluxes)
{
  gkyl_flash_app_get_kinetic_fluxes(app->block_apps[blockID], coupling_fluxes);
}
