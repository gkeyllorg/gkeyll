#include <assert.h>
#include <gkyl_pkpm_priv.h>

// initialize species moment object with the given mass factor (1.0 for the moments of
// the normalized distributions)
void
pkpm_species_moment_init_mass(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_species_moment *sm, bool is_diag,
  double mass
)
{
  sm->mcalc = gkyl_dg_updater_moment_pkpm_inew(&(struct gkyl_dg_updater_moment_pkpm_inp){
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .vel_range = &s->local_vel,
    .mass = mass,
    .is_diag = is_diag,
    .use_gpu = app->use_gpu,
  });
  int num_mom = gkyl_dg_updater_moment_pkpm_num_mom(sm->mcalc);

  sm->marr = mkarr(app->use_gpu, num_mom * app->confBasis.num_basis, app->local_ext.volume);
  sm->marr_host = sm->marr;
  if (app->use_gpu) {
    sm->marr_host = mkarr(false, num_mom * app->confBasis.num_basis, app->local_ext.volume);
  }
}

// initialize species moment object with the species mass
void
pkpm_species_moment_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_species_moment *sm, bool is_diag
)
{
  pkpm_species_moment_init_mass(app, s, sm, is_diag, s->info.mass);
}

void
pkpm_species_moment_calc(
  const struct pkpm_species_moment *sm, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin
)
{
  gkyl_dg_updater_moment_pkpm_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

// compute moments of an arbitrary phase-space array into a separate output array
// (for the moments of the configuration pass's right-hand side, keeping marr untouched)
void
pkpm_species_moment_calc_into(
  const struct pkpm_species_moment *sm, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin, struct gkyl_array *mout
)
{
  gkyl_dg_updater_moment_pkpm_advance(sm->mcalc, &phase_rng, &conf_rng, fin, mout);
}

// Diagnostic moments of the normalized state: the moments [rho, M1, p_par, p_perp, q_par,
// q_perp, r_parpar, r_parperp] of the normalized distributions [F, G] (mass-weighted)
// times n xi^k, k = [0, 1, 2, 2, 3, 3, 4, 4], n = rho/m from the fluid array (leading block) and xi the
// nodal thermal speed (weak products projected onto p=1). fluid_io is used as scratch
// (it is rebuilt by the diagnostics afterwards); pkpm_xi must be current.
void
pkpm_species_moment_diag_calc(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, const struct gkyl_range phase_rng,
  const struct gkyl_range conf_rng, const struct gkyl_array *fin, const struct gkyl_array *fluid
)
{
  int nc = app->confBasis.num_basis;
  int num_mom = gkyl_dg_updater_moment_pkpm_num_mom(s->pkpm_moms_diag.mcalc);
  static const int xi_pow[8] = {0, 1, 2, 2, 3, 3, 4, 4};
  // n = rho/m from the PHYSICAL density; the moments of the evolved J [F, G] carry one factor
  // of the position-map Jacobian, divided out at the end.
  gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &conf_rng, fluid->ncomp, fluid, s->fluid_no_J);
  gkyl_array_set_offset(s->pkpm_rho, 1.0 / s->info.mass, s->fluid_no_J, 0);
  gkyl_dg_updater_moment_pkpm_advance(s->pkpm_moms_diag.mcalc, &phase_rng, &conf_rng, fin, s->fluid_io);
  for (int c = 0; c < num_mom; ++c) {
    gkyl_dg_mul_op_range(
      &app->confBasis, c, s->pkpm_moms_diag.marr, 0, s->pkpm_rho, c, s->fluid_io, &conf_rng
    );
    for (int k = xi_pow[c]; k > 0; k -= 2) {
      gkyl_dg_mul_op_range(
        &app->confBasis, 0, s->pkpm_scratch, c, s->pkpm_moms_diag.marr, k > 1 ? 1 : 0, s->pkpm_xi,
        &conf_rng
      );
      gkyl_array_set_offset_range(s->pkpm_moms_diag.marr, 1.0, s->pkpm_scratch, c * nc, &conf_rng);
    }
  }
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &conf_rng, s->pkpm_moms_diag.marr->ncomp, s->pkpm_moms_diag.marr, s->pkpm_moms_diag.marr
  );
}

// release memory for moment data object
void
pkpm_species_moment_release(const struct gkyl_pkpm_app *app, const struct pkpm_species_moment *sm)
{
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }

  gkyl_dg_updater_moment_pkpm_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}
