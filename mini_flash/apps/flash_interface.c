// Functions required for coupling with FLASH
//
#include <gkyl_alloc.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_flash_multib_priv.h>
#include <gkyl_flash_priv.h>
#include <gkyl_flash_interface.h>

// Check that an FV array passed to/from the app has the expected number of
// components and lives on the extended local fluid range.
static void
flash_check_fv_array(
  const gkyl_flash_app *app, const struct gkyl_array *arr, int ncomp, const char *name
)
{
  if (arr->ncomp != ncomp || arr->size != app->local_fluid_ext.volume) {
    char msg[256];
    snprintf(
      msg, sizeof(msg),
      "%s: array has ncomp = %zu and size = %zu, but must have ncomp = %d "
      "and size = %ld (the extended local fluid range).",
      name, arr->ncomp, arr->size, ncomp, app->local_fluid_ext.volume
    );
    gkyl_exit(msg);
  }
}

void
gkyl_flash_app_set_fluid_vel_from_fv(gkyl_flash_app *app, struct gkyl_array *fluid_vel_in)
{
  flash_check_fv_array(app, fluid_vel_in, 3, "gkyl_flash_app_set_fluid_vel_from_fv");

  // Each species uses the same bulk velocity
  for (int i = 0; i < app->num_species; ++i) {
    struct flash_species *s = &app->species[i];

    gkyl_flash_interp_vars_vel_to_dg(
      app->interp_vars, &app->local_kinetic_ext, &app->local_fluid_ext, fluid_vel_in, s->fluid_rhou
    );

    // Ensure pkpm_moms are computed before scaling velocity by density
    flash_species_calc_pkpm_moms(app, &app->species[i], app->species[i].f);

    // Multiply by density because the PKPM updater looks for momentum density instead of velocity
    for (int d = 0; d < 3; ++d) {
      gkyl_dg_mul_op(&app->confBasis, d, s->fluid_rhou, d, s->fluid_rhou, 0, s->pkpm_moms);
    }
  }
}

void
gkyl_flash_multib_app_set_fluid_vel_from_fv_block(
  gkyl_flash_multib_app *app, int blockId_local, struct gkyl_array *fluid_vel_in
)
{
  gkyl_flash_app_set_fluid_vel_from_fv(app->block_apps[blockId_local], fluid_vel_in);
}

void
gkyl_flash_multib_app_sync_fluid_vel(gkyl_flash_multib_app *app)
{
  flash_multib_sync_fluid_rhou(app);
}

void
gkyl_flash_multib_app_set_fluid_vel_from_fv(
  gkyl_flash_multib_app *app, struct gkyl_array *fluid_vel_in[]
)
{
  for (int b = 0; b < app->num_local_blocks; ++b) {
    gkyl_flash_multib_app_set_fluid_vel_from_fv_block(app, b, fluid_vel_in[b]);
  }

  // Sync momentum density across connected blocks once all blocks are set.
  gkyl_flash_multib_app_sync_fluid_vel(app);
}

void
gkyl_flash_app_set_init_dist_moms_from_fv(gkyl_flash_app *app, struct gkyl_array *init_dist_moms)
{
  flash_check_fv_array(app, init_dist_moms, 2, "gkyl_flash_app_set_init_dist_moms_from_fv");

  for (int i = 0; i < app->num_species; ++i) {
    if (!app->species[i].init_dist_moms) {
      gkyl_exit(
        "gkyl_flash_app_set_init_dist_moms: species must be created with init_dist_from_arrays = true."
      );
    }

    gkyl_flash_interp_vars_init_dist_moms_to_dg(
      app->interp_vars, &app->local_kinetic_ext, &app->local_fluid_ext, init_dist_moms,
      app->species[i].init_dist_moms
    );
    app->species[i].has_init_dist_moms = true;
  }
}

void
gkyl_flash_app_init_maxwellian_from_moms_from_fv(
  gkyl_flash_app *app, struct gkyl_array *init_dist_moms
)
{
  gkyl_flash_app_set_init_dist_moms_from_fv(app, init_dist_moms);
  flash_app_init_maxwellian_from_stored_moms(app);
}

void
gkyl_flash_multib_app_set_init_dist_moms_from_fv_block(
  gkyl_flash_multib_app *app, int blockId_local, struct gkyl_array *init_dist_moms
)
{
  gkyl_flash_app_set_init_dist_moms_from_fv(app->block_apps[blockId_local], init_dist_moms);
}

void
gkyl_flash_multib_app_init_maxwellian_from_moms_from_fv(gkyl_flash_multib_app *app)
{
  // Verify that all species on each block app have their init_dist_moms arrays set.
  for (int b = 0; b < app->num_local_blocks; ++b) {
    for (int n = 0; n < app->block_apps[b]->num_species; ++n) {
      if (!app->block_apps[b]->species[n].has_init_dist_moms) {
        int my_rank;
        gkyl_comm_get_rank(app->comm, &my_rank);
        char msg[128];
        snprintf(
          msg, 128,
          "Trying to initialize a Maxwellian for block %d on rank %d with no init_dist_moms.\n", b,
          my_rank
        );
        gkyl_exit(msg);
      }
    }
  }

  flash_multib_app_init_maxwellian_from_stored_moms(app);
}

double
gkyl_flash_multib_app_compute_dt_max_block(gkyl_flash_multib_app *app, int blockId_local)
{
  return gkyl_flash_app_compute_dt_max(app->block_apps[blockId_local]);
}

void
gkyl_flash_multib_app_kinetic_update_block(gkyl_flash_multib_app *app, int blockId_local, double dt)
{
  struct gkyl_update_status st = flash_multib_kinetic_update_block(app, blockId_local, dt);

  // A smaller step on this block than on the others would desynchronize the blocks.
  // Written so that a NaN timestep (from NaN in the RHS) also fails.
  if (!(st.dt_actual >= dt)) {
    char msg[256];
    snprintf(
      msg, sizeof(msg),
      "gkyl_flash_multib_app_kinetic_update_block: local block %d (global %d) "
      "can only take dt = %.6e < requested dt = %.6e. dt must not exceed gkyl_flash_multib_app_compute_dt_max.",
      blockId_local, app->local_blocks[blockId_local], st.dt_actual, dt
    );
    gkyl_exit(msg);
  }
}

void
gkyl_flash_multib_app_kinetic_update_sync(gkyl_flash_multib_app *app, double dt)
{
  flash_multib_kinetic_update_sync(app, dt);
}

void
gkyl_flash_app_get_kinetic_fluxes_to_fv(gkyl_flash_app *app, struct gkyl_array *coupling_fluxes)
{
  flash_check_fv_array(app, coupling_fluxes, 3, "gkyl_flash_app_get_kinetic_fluxes_to_fv");

  for (int i = 0; i < app->num_species; ++i) {
    // Verify that array exists
    if (coupling_fluxes->elemsz == 0) {
      char msg[128];
      snprintf(
        msg, 128,
        "Output coupling flux array is empty. Check that the proper number of arrays have been initialized. Exiting..."
      );
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
    gkyl_flash_interp_vars_coupling_fluxes_to_fv(
      app->interp_vars, &app->local_kinetic_ext, &app->local_fluid_ext,
      app->species[i].flash_coupling_fluxes, coupling_fluxes
    );
  }
}

void
gkyl_flash_multib_app_get_kinetic_fluxes_to_fv_block(
  gkyl_flash_multib_app *app, int blockID, struct gkyl_array *coupling_fluxes
)
{
  gkyl_flash_app_get_kinetic_fluxes_to_fv(app->block_apps[blockID], coupling_fluxes);
}
