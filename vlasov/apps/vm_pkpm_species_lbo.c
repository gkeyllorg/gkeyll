#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Primitive moments at the nodes.
static void
vm_pkpm_species_lbo_moms_nodal(
  gkyl_vlasov_app *app, const struct vlasov_species *species, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  gkyl_pkpm_lbo_moms_advance(
    lbo->up_moms, &app->local, &species->kinetic->local, species->kinetic->pkpm.nmoms.marr,
    species->kinetic->pkpm.f_nodal, lbo->prim_moms
  );
}

// The operator at the nodes, accumulated onto the nodal right-hand side of the species.
static void
vm_pkpm_species_lbo_kinetic_rhs_nodal(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  gkyl_dg_updater_lbo_pkpm_nodal_advance(
    lbo->up_coll, &kin->local, kin->pkpm.f_nodal, kin->cflrate, kin->pkpm.rhs_nodal
  );
}

// The operator at the configuration-space nodes: the collision frequency (the profile of the
// input, scaled by nu_frac; without one, the species' collision-frequency object computes it every
// stage), the updater for the primitive moments and the updater for the operator.
static void
vm_pkpm_species_lbo_init_nodal(
  struct gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_lbo_collisions *lbo
)
{
  int num_nodes = app->basis.num_basis;
  double nu_frac = s->kinetic->info.collisions.nu_frac > 0.0 ? s->kinetic->info.collisions.nu_frac :
                                                               1.0;
  lbo->nu_self = vm_pkpm_species_coll_self_nu_new(app, s, nu_frac);
  // Without partners the self frequency is the frequency of the operator. With partners the
  // frequency of the operator is rebuilt from it every stage.
  if (s->kinetic->info.collisions.num_cross_collisions) {
    lbo->nu = mkarr(false, num_nodes, app->local_ext.volume);
    gkyl_array_copy(lbo->nu, lbo->nu_self);
  } else {
    lbo->nu = gkyl_array_acquire(lbo->nu_self);
  }

  double vel_bounds[2] = {s->kinetic->info.lower[0], s->kinetic->info.upper[0]};
  lbo->up_moms = gkyl_pkpm_lbo_moms_inew(&(struct gkyl_pkpm_lbo_moms_inp){
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .vel_bounds = vel_bounds,
    .use_gpu = false,
  });
  lbo->up_coll = gkyl_dg_updater_lbo_pkpm_nodal_inew(&(struct gkyl_dg_updater_lbo_pkpm_nodal_inp){
    .phase_grid = &s->kinetic->grid,
    .conf_basis = &app->basis,
    .phase_basis = &s->kinetic->basis,
    .conf_range = &app->local,
    .nu = lbo->nu,
    .prim_moms = lbo->prim_moms,
    .use_gpu = false,
  });
  lbo->kinetic_rhs = vm_pkpm_species_lbo_kinetic_rhs_nodal;
  lbo->calc_moms = vm_pkpm_species_lbo_moms_nodal;
}

void
vm_pkpm_species_lbo_init(
  struct gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_lbo_collisions *lbo
)
{
  assert(!app->use_gpu); // the operator at the nodes is not available on the GPU yet
  int vdim = 1; // the normalized parallel velocity

  // Primitive moments [u, vtsq] at the nodes.
  lbo->prim_moms = mkarr(false, (vdim + 1) * app->basis.num_basis, app->local_ext.volume);

  vm_pkpm_species_lbo_init_nodal(app, s, lbo);
}

void
vm_pkpm_species_lbo_moms(
  gkyl_vlasov_app *app, const struct vlasov_species *species, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  struct timespec wst = gkyl_wall_clock();
  lbo->calc_moms(app, species, lbo, fin);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_pkpm_species_lbo_rhs(
  gkyl_vlasov_app *app, const struct vm_species *kin, struct vm_pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();

  // accumulate update due to collisions onto rhs
  lbo->kinetic_rhs(app, kin, lbo, fin, rhs);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_pkpm_species_lbo_release(
  const struct gkyl_vlasov_app *app, const struct vm_pkpm_lbo_collisions *lbo
)
{
  gkyl_array_release(lbo->prim_moms);
  gkyl_array_release(lbo->nu);
  gkyl_array_release(lbo->nu_self);
  gkyl_pkpm_lbo_moms_release(lbo->up_moms);
  gkyl_dg_updater_lbo_pkpm_nodal_release(lbo->up_coll);
}
