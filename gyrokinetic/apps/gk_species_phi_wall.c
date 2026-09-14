#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_phi_wall_advance_disabled(gkyl_gyrokinetic_app *app,
  const struct gk_species_wall_potential *wall, double tm)
{
}

static void
gk_species_phi_wall_advance_enabled(gkyl_gyrokinetic_app *app,
  const struct gk_species_wall_potential *wall, double tm)
{
  gkyl_eval_on_nodes_advance(wall->projector, tm, &app->local_ext, wall->phi_host);
  if (app->use_gpu)
    gkyl_array_copy(wall->phi, wall->phi_host);
}

void
gk_species_phi_wall_init(gkyl_gyrokinetic_app *app,
  const struct gkyl_gyrokinetic_bc *bc, struct gk_species_wall_potential *wall)
{
  *wall = (struct gk_species_wall_potential) {
    .advance_func = gk_species_phi_wall_advance_disabled,
  };
  if (bc->type != GKYL_BC_GK_SPECIES_SHEATH)
    return;

  wall->phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(wall->phi, 0.0);

  wall->phi_host = wall->phi;
  if (bc->aux_profile) {
    if (app->use_gpu)
      wall->phi_host = mkarr(false, wall->phi->ncomp, wall->phi->size);

    wall->projector = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
      1, bc->aux_profile, bc->aux_ctx);
    wall->advance_func = gk_species_phi_wall_advance_enabled;
    gk_species_phi_wall_advance(app, wall, 0.0);
  }
}

void
gk_species_phi_wall_advance(gkyl_gyrokinetic_app *app,
  const struct gk_species_wall_potential *wall, double tm)
{
  wall->advance_func(app, wall, tm);
}

void
gk_species_phi_wall_release(const gkyl_gyrokinetic_app *app,
  const struct gk_species_wall_potential *wall)
{
  if (!wall->phi)
    return;

  gkyl_array_release(wall->phi);
  if (wall->projector) {
    gkyl_eval_on_nodes_release(wall->projector);
    if (app->use_gpu)
      gkyl_array_release(wall->phi_host);
  }
}
