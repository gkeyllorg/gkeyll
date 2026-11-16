#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_pkpm_priv.h>

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
  pkpm_em->slvr = gkyl_dg_calc_pkpm_em_coupling_inew(&(struct gkyl_dg_calc_pkpm_em_coupling_inp){
    .conf_basis = &app->confBasis,
    .mem_range = &app->local,
    .num_species = num_species,
    .qbym = qbym,
    .epsilon0 = app->field->info.epsilon0,
    .pkpm_field_static = app->field->info.is_static,
    .use_gpu = app->use_gpu,
  });

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
    gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &app->local, s->fluid->fluid->ncomp, s->fluid->fluid, s->fluid->fluid_no_J);
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

  // The currents of the solve besides the species': the applied current, the ghost current and
  // the effective current of the curl. A static field takes none of them (curl_current stays
  // zero), so its E does not change.
  const struct gkyl_array *app_current = app->field->curl_current;
  if (!app->field->info.is_static) {
    struct pkpm_field *field = app->field;
    // The ghost current enters the implicit solve as an applied current -<J_x> (frozen over
    // the step, like the applied current itself); its domain mean is that of the J-weighted
    // (evolved) momenta, and it is used as a physical current like the rest of the solve.
    const struct gkyl_array *applied = field->app_current_no_J;
    if (field->use_ghost_current) {
      const struct gkyl_array *fluids_J[GKYL_MAX_SPECIES];
      for (int i = 0; i < num_species; ++i) {
        fluids_J[i] = app->species[i].fluid->fluid;
      }
      pkpm_field_calc_ghost_current(app, field, fluids_J);
      gkyl_array_set(field->app_current_tot, 1.0, field->app_current_no_J);
      gkyl_array_accumulate_range(field->app_current_tot, -1.0, field->ghost_current, &app->local);
      applied = field->app_current_tot;
    }
    // Stormer-Verlet in the curl terms: E advances only here, with c^2 curl B of the field at
    // the start of the half step (B is not changed by the sources) entering as the effective
    // current -eps0 c^2 curl B; the transport step then drifts B with this E.
    gkyl_array_clear(field->curl_rhs, 0.0);
    gkyl_hyper_dg_advance(field->slvr, &app->local, field->em, field->cflrate, field->curl_rhs);
    gkyl_vlasov_position_map_divide_jacobpos_conf(
      app->pos_map, &app->local, field->curl_rhs->ncomp, field->curl_rhs, field->curl_rhs_no_J
    );
    gkyl_array_set(field->curl_current, 1.0, applied);
    gkyl_array_accumulate_offset(field->curl_current, -field->info.epsilon0, field->curl_rhs_no_J, 0);
  }
  gkyl_dg_calc_pkpm_em_coupling_advance(
    pkpm_em->slvr, dt, app_accels, app->field->ext_em, app_current, fluids, app->field->em_no_J
  );

  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_vlasov_position_map_rescale_jacobpos_conf(app->pos_map, &app->local, s->fluid->fluid_no_J, s->fluid->fluid);
    pkpm_species_apply_positivity(app, s->fluid, s->fluid->fluid);
    pkpm_fluid_species_apply_bc(app, s->fluid, s->fluid->fluid);
  }
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local, app->field->em_no_J, app->field->em
  );
  pkpm_field_apply_bc(app, app->field, app->field->em);
}

void
pkpm_fluid_em_coupling_release(struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em)
{
  gkyl_dg_calc_pkpm_em_coupling_release(pkpm_em->slvr);
  gkyl_free(pkpm_em);
}
