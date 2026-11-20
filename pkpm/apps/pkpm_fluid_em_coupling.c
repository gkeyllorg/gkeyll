#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_pkpm_priv.h>

// The solve of a source step. It advances the momentum of the fluids and the electric field, and
// the energy of the fluids by the change of the bulk kinetic energy (the work of the field): the
// density is unchanged, and so is the internal energy at the nodes of the solve. It acts on the
// physical fluids and field (without the position-map Jacobian).
static void
fluid_em_coupling_advance(
  struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em, double dt,
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES], struct gkyl_array *fluids[GKYL_MAX_SPECIES]
)
{
  gkyl_dg_calc_fluid_em_coupling_advance(
    pkpm_em->slvr, dt, app_accels, app->field->ext_em, app->field->curl_current, fluids,
    app->field->em_no_J
  );
}

// Domain mean of E_x of the physical field of the solve. The mean is taken of the field with
// the position-map Jacobian, whose cell averages sum to the integral over the physical domain;
// the evolved field is used as scratch for it (the update rescales into it at the end).
static double
fluid_em_coupling_ex_mean(struct gkyl_pkpm_app *app)
{
  struct pkpm_field *field = app->field;
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local, field->em_no_J, field->em
  );
  return pkpm_field_ghost_current_mean(app, field, 1.0, field->em);
}

// The solve with the ghost current (1x periodic): a uniform current along x whose strength keeps
// the average of E_x at its value at the start of the step, so that the average current of the
// fluids over the step drives no uniform field. The solve is linear in the current at a fixed
// initial state: two trial solves, without and with a unit ghost current, give the strength,
// and a third takes the step (the energy of the fluids is not linear in the current).
static void
fluid_em_coupling_advance_ghost(
  struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em, double dt,
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES], struct gkyl_array *fluids[GKYL_MAX_SPECIES]
)
{
  struct pkpm_field *field = app->field;
  int num_species = app->num_species;

  for (int i = 0; i < num_species; ++i) {
    gkyl_array_copy(pkpm_em->fluid_ghost[i], fluids[i]);
  }
  gkyl_array_copy(pkpm_em->em_ghost, field->em_no_J);
  double ex_mean_start = pkpm_field_ghost_current_mean(app, field, 1.0, field->em);

  double ghost_current[3] = {0.0, 1.0, 0.0}; // the trial strengths; the last is the solution
  double ex_mean_trial[2];
  double ghost_current_prev = 0.0;
  for (int trial = 0; trial < 3; ++trial) {
    if (trial == 2) {
      ghost_current[2] = (ex_mean_start - ex_mean_trial[0]) / (ex_mean_trial[1] - ex_mean_trial[0]);
    }
    gkyl_array_shiftc(field->curl_current, ghost_current[trial] - ghost_current_prev, 0);
    ghost_current_prev = ghost_current[trial];
    for (int i = 0; i < num_species; ++i) {
      gkyl_array_copy(fluids[i], pkpm_em->fluid_ghost[i]);
    }
    gkyl_array_copy(field->em_no_J, pkpm_em->em_ghost);
    fluid_em_coupling_advance(app, pkpm_em, dt, app_accels, fluids);
    if (trial < 2) {
      ex_mean_trial[trial] = fluid_em_coupling_ex_mean(app);
    }
  }
}

// initialize implicit fluid-EM coupling object for PKPM system
struct pkpm_fluid_em_coupling *
pkpm_fluid_em_coupling_init(struct gkyl_pkpm_app *app)
{
  struct pkpm_fluid_em_coupling *pkpm_em = gkyl_malloc(sizeof(struct pkpm_fluid_em_coupling));

  int num_species = app->num_species;
  double qbym[GKYL_MAX_SPECIES] = {0.0};
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    qbym[i] = s->charge / s->mass;
  }

  // Initialize solver
  pkpm_em->slvr = gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
    .conf_basis = &app->confBasis,
    .mem_range = &app->local,
    .num_species = num_species,
    .qbym = qbym,
    .epsilon0 = app->field->info.epsilon0,
    .is_field_static = app->field->info.is_static,
    .use_gpu = app->use_gpu,
  });

  // Ghost current: the fluids' average current over a source step is removed (a static field
  // takes no current).
  pkpm_em->advance_func = fluid_em_coupling_advance;
  pkpm_em->em_ghost = 0;
  if (app->field->use_ghost_current && !app->field->info.is_static) {
    pkpm_em->advance_func = fluid_em_coupling_advance_ghost;
    for (int i = 0; i < num_species; ++i) {
      const struct gkyl_array *fluid = app->species[i].fluid->fluid;
      pkpm_em->fluid_ghost[i] = mkarr(app->use_gpu, fluid->ncomp, fluid->size);
    }
    pkpm_em->em_ghost = mkarr(app->use_gpu, app->field->em->ncomp, app->field->em->size);
  }

  return pkpm_em;
}

void
pkpm_fluid_em_coupling_update(
  struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em, double tcurr, double dt
)
{
  int num_species = app->num_species;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];

  // The nodal implicit solve advances the fluid state [rho, M, E] (momentum by the solve,
  // energy by the work of the field) and the electric field in place. It works with the
  // PHYSICAL state: the evolved fluid and field carry the position-map Jacobian, which is
  // divided out before and restored after (identity map: copies).
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_vlasov_position_map_divide_jacobpos_conf(
      app->pos_map, &app->local, s->fluid->fluid->ncomp, s->fluid->fluid, s->fluid->fluid_no_J
    );
    fluids[i] = s->fluid->fluid_no_J;
  }
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local, app->field->em->ncomp, app->field->em, app->field->em_no_J
  );

  // Compute applied acceleration, external EM fields, and applied currents if present.
  // Note: these quantities use proj_on_basis so we copy to GPU every call if app->use_gpu = true.
  if (app->field->ext_em_evolve) {
    pkpm_field_calc_ext_em(app, app->field, tcurr);
  }
  if (app->field->app_current_evolve) {
    pkpm_field_calc_app_current(app, app->field, tcurr);
  }
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    if (s->fluid->app_accel_evolve) {
      pkpm_species_calc_app_accel(app, s, tcurr);
    }
    app_accels[i] = s->fluid->app_accel;
  }

  // The currents of the solve besides the species': the applied current and the effective
  // current of the curl. A static field takes neither (curl_current stays zero), so its E does
  // not change.
  struct pkpm_field *field = app->field;
  if (!field->info.is_static) {
    // Stormer-Verlet in the curl terms: E advances only here, with c^2 curl B of the field at
    // the start of the step (B is not changed by the sources) entering as the effective current
    // -eps0 c^2 curl B; the transport step then drifts B with this E.
    gkyl_array_clear(field->curl_rhs, 0.0);
    gkyl_hyper_dg_advance(field->slvr, &app->local, field->em, field->cflrate, field->curl_rhs);
    gkyl_vlasov_position_map_divide_jacobpos_conf(
      app->pos_map, &app->local, field->curl_rhs->ncomp, field->curl_rhs, field->curl_rhs_no_J
    );
    gkyl_array_set(field->curl_current, 1.0, field->app_current_no_J);
    gkyl_array_accumulate_offset(
      field->curl_current, -field->info.epsilon0, field->curl_rhs_no_J, 0
    );
  }

  pkpm_em->advance_func(app, pkpm_em, dt, app_accels, fluids);

  // Apply BCs to the updated fluids and field.
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_vlasov_position_map_rescale_jacobpos_conf(
      app->pos_map, &app->local, s->fluid->fluid_no_J, s->fluid->fluid
    );
    pkpm_fluid_species_apply_bc(app, s->fluid, s->fluid->fluid);
  }
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local, field->em_no_J, field->em
  );
  pkpm_field_apply_bc(app, field, field->em);
}

void
pkpm_fluid_em_coupling_limit(struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em)
{
  // The positivity limiter keeps the internal energy positive at the vertices of the fluids.
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    pkpm_species_apply_positivity(app, s->fluid, s->fluid->fluid);
    pkpm_fluid_species_apply_bc(app, s->fluid, s->fluid->fluid);
  }
}

void
pkpm_fluid_em_coupling_release(struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em)
{
  gkyl_dg_calc_fluid_em_coupling_release(pkpm_em->slvr);
  if (pkpm_em->em_ghost) {
    for (int i = 0; i < app->num_species; ++i) {
      gkyl_array_release(pkpm_em->fluid_ghost[i]);
    }
    gkyl_array_release(pkpm_em->em_ghost);
  }
  gkyl_free(pkpm_em);
}
