#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_orbit_average_gyrokinetic.h>
#include <limits.h>

static int
gk_orbit_reduce(void *ctx, bool sum, long count, const double *input, double *output)
{
  struct gkyl_comm *comm = ctx;
  while (count > 0) {
    int chunk = count > INT_MAX ? INT_MAX : (int)count;
    int error =
      gkyl_comm_allreduce(comm, GKYL_DOUBLE, sum ? GKYL_SUM : GKYL_MIN, chunk, input, output);
    if (error) {
      return error;
    }
    count -= chunk;
    input += chunk;
    output += chunk;
  }
  return 0;
}

static enum gkyl_gk_trapped_passing_orbit_type
gk_orbit_boundary(enum gkyl_gyrokinetic_bc_type type)
{
  switch (type) {
    case GKYL_BC_GK_SPECIES_REFLECT:
    case GKYL_BC_GK_SPECIES_ZERO_FLUX:
      return GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL;
    case GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING:
    case GKYL_BC_GK_SPECIES_SHEATH_SURROGATE:
      return GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH;
    default:
      return GKYL_GK_TRAP_PASS_ORBIT_PASSING;
  }
}

void
gk_species_orbit_average_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_orbit_average *orbit
)
{
  if (gks->info.orbit_average.enabled && app->is_multib) {
    gkyl_exit("Orbit averaging requires the single-block app time integrator.");
  }
  *orbit = (struct gk_orbit_average){};
}

bool
gk_species_orbit_average_prepare(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  const struct gkyl_gyrokinetic_orbit_average *options = &gks->info.orbit_average;
  if (!options->enabled) {
    return true;
  }

  // The frozen projector describes parallel electrostatic orbits,
  // independently of the collisionless RHS. It cannot represent electromagnetic orbits,
  // transverse drifts, or periodic circulation by reflecting endpoints.
  enum gkyl_gk_collisionless_type type = gks->info.collisionless.type;
  if (gks->is_fluid || gks->info.is_static || app->cdim != 1 || gks->info.vdim != 2 ||
      app->poly_order != 1 || app->field->update_field || app->num_periodic_dir ||
      app->field->gkfield_id == GKYL_GK_FIELD_EM || gks->info.flr.type || gks->info.scaling.type ||
      (type != GKYL_GK_COLLISIONLESS_NONE && type != GKYL_GK_COLLISIONLESS_ES &&
       type != GKYL_GK_COLLISIONLESS_ES_NO_BY) ||
      (!options->frequency_is_infinite && options->frequency < 0.0) ||
      options->num_energy_nodes < 0 || options->num_energy_nodes == 1) {
    gkyl_gyrokinetic_app_cout(
      app, stderr,
      "Orbit averaging for %s requires a dynamic 1x2v P1 kinetic species, a frozen "
      "electrostatic field, nonperiodic boundaries, no FLR/scaling, and valid rate/knots.\n",
      gks->info.name
    );
    return false;
  }
  for (int edge = 0; edge < 2; ++edge) {
    const struct gkyl_gyrokinetic_bc *bc = edge ? &gks->upper_bc[0] : &gks->lower_bc[0];
    if (bc->type == GKYL_BC_GK_SPECIES_PERIODIC || bc->type == GKYL_BC_GK_SPECIES_TWISTSHIFT ||
        bc->aux_profile) {
      gkyl_gyrokinetic_app_cout(
        app, stderr,
        "Orbit averaging does not support periodic/twistshift or prescribed evolving wall potentials.\n"
      );
      return false;
    }
  }
  struct gk_orbit_average *orbit = &gks->orbit_average;
  if (orbit->up) {
    return true;
  }

  // These communication and geometry buffers exist only at frozen-cache setup.
  // Global arrays retain global_ext strides, just like the supplied global ranges.
  struct gkyl_array *bmag = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  struct gkyl_array *phi = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  struct gkyl_array *jacobgeo = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  struct gkyl_array *mask = mkarr(app->use_gpu, 1, gks->global_ext.volume);
  struct gkyl_array *wall_lo = 0, *wall_up = 0;
  enum gkyl_gk_trapped_passing_orbit_type lower = gk_orbit_boundary(gks->lower_bc[0].type);
  enum gkyl_gk_trapped_passing_orbit_type upper = gk_orbit_boundary(gks->upper_bc[0].type);
  int error = gkyl_comm_array_allgather(
    app->comm, &app->local, &app->global, app->gk_geom->geo_corn.bmag, bmag
  );
  error |=
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->field->phi_smooth, phi);
  error |= gkyl_comm_array_allgather(
    app->comm, &app->local, &app->global, app->gk_geom->geo_int.jacobgeo, jacobgeo
  );
  if (lower == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH) {
    wall_lo = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    error |= gkyl_comm_array_allgather(
      app->comm, &app->local, &app->global, gks->phi_wall_lo.phi, wall_lo
    );
  }
  if (upper == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH) {
    wall_up = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    error |= gkyl_comm_array_allgather(
      app->comm, &app->local, &app->global, gks->phi_wall_up.phi, wall_up
    );
  }
  enum gkyl_orbit_average_gyrokinetic_status status = GKYL_ORBIT_AVERAGE_INVALID_INPUT;
  if (!error) {
    struct gkyl_loss_cone_mask_gyrokinetic_inp mask_inp = {
      .conf_basis = &app->basis,
      .vel_map = gks->vel_map,
      .mass = gks->info.mass,
      .charge = gks->info.charge,
      .use_gpu = app->use_gpu,
      .lower_orbit = lower,
      .upper_orbit = upper,
    };
    gkyl_loss_cone_mask_gyrokinetic *mask_up = gkyl_loss_cone_mask_gyrokinetic_inew(&mask_inp);
    gkyl_loss_cone_mask_gyrokinetic_advance(
      mask_up, &gks->global, &app->global, bmag, phi, wall_lo, wall_up, mask
    );
    struct gkyl_orbit_average_gyrokinetic_inp inp = {
      .grid = &gks->grid,
      .conf_basis = &app->basis,
      .phase_basis = &gks->basis,
      .phase_range = &gks->global,
      .local_range = &gks->local,
      .reduce = gk_orbit_reduce,
      .reduce_ctx = gks->comm,
      .conf_range = &app->global,
      .vel_map = gks->vel_map,
      .bmag = bmag,
      .phi = phi,
      .jacobgeo = jacobgeo,
      .mask = mask,
      .phi_wall_lo = wall_lo,
      .phi_wall_up = wall_up,
      .lower_orbit = lower,
      .upper_orbit = upper,
      .mass = gks->info.mass,
      .charge = gks->info.charge,
      .num_energy_nodes = options->num_energy_nodes ? options->num_energy_nodes : 16,
      .separatrix_margin = options->separatrix_margin,
      .positivity = options->positivity,
      .frequency_is_infinite = options->frequency_is_infinite,
      .use_gpu = app->use_gpu,
    };
    orbit->up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
    gkyl_loss_cone_mask_gyrokinetic_release(mask_up);
  }
  gkyl_array_release(bmag);
  gkyl_array_release(phi);
  gkyl_array_release(jacobgeo);
  gkyl_array_release(mask);
  if (wall_lo) {
    gkyl_array_release(wall_lo);
  }
  if (wall_up) {
    gkyl_array_release(wall_up);
  }
  if (!orbit->up) {
    gkyl_gyrokinetic_app_cout(
      app, stderr, "Orbit-average setup failed for %s (status %d, communication %d).\n",
      gks->info.name, status, error
    );
    return false;
  }
  struct gkyl_orbit_average_gyrokinetic_info info =
    gkyl_orbit_average_gyrokinetic_get_info(orbit->up);
  if (!info.num_selected_cells) {
    gkyl_gyrokinetic_app_cout(
      app, stdout,
      "Orbit-average warning: no accepted closed cells for %s; relaxation is the identity.\n",
      gks->info.name
    );
  }
  return true;
}

void
gk_species_orbit_average_advance(gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt)
{
  if (!gks->info.orbit_average.enabled) {
    return;
  }
  struct gk_orbit_average *orbit = &gks->orbit_average;
  struct gkyl_orbit_average_gyrokinetic_result result = gkyl_orbit_average_gyrokinetic_advance(
    orbit->up, dt, gks->info.orbit_average.frequency, gks->f
  );
  if (!result.success) {
    gkyl_exit(
      "Orbit relaxation failed: check timestep/rate and positivity of the incoming DG state."
    );
  }
  orbit->min_blend = result.min_blend;
  orbit->num_limited_blocks = result.num_limited_blocks;
}

void
gyrokinetic_update_orbit_average(gkyl_gyrokinetic_app *app, double dt)
{
  bool active = false;
  for (int i = 0; i < app->num_species; ++i) {
    active |= app->species[i].info.orbit_average.enabled;
    gk_species_orbit_average_advance(app, &app->species[i], dt);
  }
  if (!active) {
    return;
  }
  struct gkyl_array *fin[app->num_species], **bflux[app->num_species];
  struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i = 0; i < app->num_species; ++i) {
    fin[i] = app->species[i].f;
    bflux[i] = app->species[i].bflux.f;
  }
  for (int i = 0; i < app->num_neut_species; ++i) {
    fin_neut[i] = app->neut_species[i].f;
  }
  gyrokinetic_calc_field_and_apply_bc(app, app->tcurr + dt, fin, bflux, fin_neut);
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    if (gks->info.orbit_average.enabled) {
      gk_species_calc_int_mom_dt(app, gks, dt, gks->fdot_mom_new);
    }
  }
}

struct gkyl_gyrokinetic_orbit_average_stat
gkyl_gyrokinetic_app_species_orbit_average_stat(
  const gkyl_gyrokinetic_app *app, const char *species_name
)
{
  struct gkyl_gyrokinetic_orbit_average_stat stat = {};
  const struct gk_species *gks = gk_find_species(app, species_name);
  if (!gks) {
    return stat;
  }
  stat.enabled = gks->info.orbit_average.enabled;
  if (gks->orbit_average.up) {
    stat.prepared = true;
    struct gkyl_orbit_average_gyrokinetic_info info =
      gkyl_orbit_average_gyrokinetic_get_info(gks->orbit_average.up);
    stat.num_selected_cells = info.num_selected_cells;
    stat.num_rejected_cells = info.num_rejected_cells;
    stat.num_modes = info.num_modes;
    stat.min_blend = gks->orbit_average.min_blend;
    stat.num_limited_blocks = gks->orbit_average.num_limited_blocks;
  }
  return stat;
}

void
gk_species_orbit_average_release(const struct gk_orbit_average *orbit)
{
  if (orbit->up) {
    gkyl_orbit_average_gyrokinetic_release(orbit->up);
  }
}

void
gk_species_orbit_average_reset(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_orbit_average *orbit,
  struct gkyl_gyrokinetic_orbit_average inp
)
{
  gk_species_orbit_average_release(orbit);
  gks->info.orbit_average = inp;
  gk_species_orbit_average_init(app, gks, orbit);

  // Keep the split while any species needs an implicit substep.
  bool use_split = app->has_implicit_coll_scheme;
  for (int i = 0; i < app->num_species; ++i) {
    use_split |= app->species[i].info.orbit_average.enabled;
  }
  app->update_func = use_split ? gyrokinetic_update_op_split : gyrokinetic_update_ssp_rk3;
}
