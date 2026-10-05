#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_collisionless_influx_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gk_collisionless_influx *influx,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
)
{
}

static void
gk_species_collisionless_influx_sheath_flux(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gk_collisionless_influx *influx,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
)
{
  if (influx->lower_sheath_flux) {
    gkyl_gk_sheath_conducting_flux_advance(
      influx->lower_sheath_flux, app->field->phi_smooth, species->gyro_phi, app->field->phi_wall_lo,
      fin, flux_surf
    );
  }
  if (influx->upper_sheath_flux) {
    gkyl_gk_sheath_conducting_flux_advance(
      influx->upper_sheath_flux, app->field->phi_smooth, species->gyro_phi, app->field->phi_wall_up,
      fin, flux_surf
    );
  }
}

void
gk_species_collisionless_influx_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_collisionless_influx *influx
)
{
  *influx = (struct gk_collisionless_influx){.advance = gk_species_collisionless_influx_disabled};

  int par_dir = app->cdim - 1;
  bool lower = gks->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_SHEATH_FLUX;
  bool upper = gks->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_SHEATH_FLUX;
  for (int dir = 0; dir < par_dir; ++dir) {
    assert(gks->lower_bc[dir].type != GKYL_BC_GK_SPECIES_SHEATH_FLUX);
    assert(gks->upper_bc[dir].type != GKYL_BC_GK_SPECIES_SHEATH_FLUX);
  }
  if (!(lower || upper)) {
    return;
  }

  int vdim = gks->local.ndim - app->cdim;
  assert(gks->collisionless.collisionless_id == GKYL_GK_COLLISIONLESS_ES);
  assert(gks->basis.poly_order == 1);
  assert(vdim == 2);

  struct gkyl_gk_sheath_conducting_flux_inp inp = {
    .dir = par_dir,
    .cdim = app->cdim,
    .phase_grid = &gks->grid,
    .phase_basis = &gks->basis,
    .phase_basis_cu = gks->basis_on_dev,
    .conf_range = &app->local,
    .conf_ext_range = &app->local_ext,
    .phase_ext_range = &gks->local_ext,
    .flux_op = gks->collisionless.surf_flux_op,
    .charge = gks->info.charge,
    .mass = gks->info.mass,
    .use_gpu = app->use_gpu,
  };
  if (lower) {
    inp.edge = GKYL_LOWER_EDGE;
    inp.skin_range = &gks->local_lower_skin[par_dir];
    influx->lower_sheath_flux = gkyl_gk_sheath_conducting_flux_inew(&inp);
  }
  if (upper) {
    inp.edge = GKYL_UPPER_EDGE;
    inp.skin_range = &gks->local_upper_skin[par_dir];
    influx->upper_sheath_flux = gkyl_gk_sheath_conducting_flux_inew(&inp);
  }
  influx->advance = gk_species_collisionless_influx_sheath_flux;
}

void
gk_species_collisionless_influx_advance(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gk_collisionless_influx *influx,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
)
{
  influx->advance(app, species, influx, fin, flux_surf);
}

void
gk_species_collisionless_influx_release(const struct gk_collisionless_influx *influx)
{
  gkyl_gk_sheath_conducting_flux_release(influx->lower_sheath_flux);
  gkyl_gk_sheath_conducting_flux_release(influx->upper_sheath_flux);
}
