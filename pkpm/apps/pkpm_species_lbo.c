#include <assert.h>
#include <gkyl_pkpm_priv.h>

// Primitive moments at the nodes.
static void
pkpm_species_lbo_moms_nodal(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  gkyl_pkpm_lbo_moms_advance(
    lbo->up_moms, &app->local, &species->local, species->pkpm_nmoms.marr, species->f_nodal,
    lbo->prim_moms
  );
}

// The operator at the nodes, accumulated onto the nodal right-hand side of the species.
static void
pkpm_species_lbo_kinetic_rhs_nodal(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  gkyl_dg_updater_lbo_pkpm_nodal_advance(
    lbo->up_coll, &species->local, species->f_nodal, species->cflrate_f, species->rhs_nodal
  );
}

// The operator at the configuration-space nodes: the collision frequency (the profile of the
// input, scaled by nuFrac; without one, the species' collision-frequency object computes it every
// stage), the updater for the primitive moments and the updater for the operator.
static void
pkpm_species_lbo_init_nodal(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
)
{
  int num_nodes = app->confBasis.num_basis;
  double nu_frac = s->info.collisions.nuFrac > 0.0 ? s->info.collisions.nuFrac : 1.0;
  lbo->nu_self = pkpm_species_coll_self_nu_new(app, s, nu_frac);
  // Without partners the self frequency is the frequency of the operator. With partners the
  // frequency of the operator is rebuilt from it every stage.
  if (s->info.collisions.num_cross_collisions) {
    lbo->nu = mkarr(false, num_nodes, app->local_ext.volume);
    gkyl_array_copy(lbo->nu, lbo->nu_self);
  } else {
    lbo->nu = gkyl_array_acquire(lbo->nu_self);
  }

  double vel_bounds[2] = {s->info.lower[0], s->info.upper[0]};
  lbo->up_moms = gkyl_pkpm_lbo_moms_new(&app->confBasis, &app->basis, vel_bounds, false);
  struct gkyl_dg_lbo_pkpm_nodal_auxfields nodal_inp = {.nu = lbo->nu, .prim_moms = lbo->prim_moms};
  lbo->up_coll = gkyl_dg_updater_lbo_pkpm_nodal_new(
    &s->grid, &app->confBasis, &app->basis, &app->local, &nodal_inp, false
  );
  lbo->kinetic_rhs = pkpm_species_lbo_kinetic_rhs_nodal;
  lbo->calc_moms = pkpm_species_lbo_moms_nodal;
}

void
pkpm_species_lbo_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
)
{
  assert(!app->use_gpu); // the operator at the nodes is not available on the GPU yet
  int vdim = app->vdim;

  // Primitive moments [u, vtsq] at the nodes.
  lbo->prim_moms = mkarr(false, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);

  pkpm_species_lbo_init_nodal(app, s, lbo);
}

void
pkpm_species_lbo_moms(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  struct timespec wst = gkyl_wall_clock();
  lbo->calc_moms(app, species, lbo, fin);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_lbo_rhs(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();

  // accumulate update due to collisions onto rhs
  lbo->kinetic_rhs(app, species, lbo, fin, rhs);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_lbo_release(const struct gkyl_pkpm_app *app, const struct pkpm_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->prim_moms);
  gkyl_array_release(lbo->nu);
  gkyl_array_release(lbo->nu_self);
  gkyl_pkpm_lbo_moms_release(lbo->up_moms);
  gkyl_dg_updater_lbo_pkpm_nodal_release(lbo->up_coll);
}
