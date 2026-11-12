#include <assert.h>
#include <gkyl_pkpm_priv.h>

// Boundary corrections, primitive moments and collision frequency as p=1 expansions (GPU runs).
static void
pkpm_species_lbo_moms_modal(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  // construct boundary corrections
  gkyl_mom_calc_bcorr_advance(
    lbo->bcorr_calc, &species->local, &app->local, fin, lbo->boundary_corrections
  );

  // Compute primitive moments of the normalized distributions.
  gkyl_prim_lbo_calc_advance(
    lbo->coll_pcalc, &app->local, species->pkpm_nmoms.marr, lbo->boundary_corrections, lbo->self_nu,
    lbo->prim_moms
  );

  for (int d = 0; d < app->vdim; d++) {
    gkyl_dg_mul_op(&app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
  }
  gkyl_dg_mul_op(
    &app->confBasis, app->vdim, lbo->nu_prim_moms, app->vdim, lbo->prim_moms, 0, lbo->self_nu
  );
}

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

// The Galerkin drag and diffusion updaters (GPU runs).
static void
pkpm_species_lbo_kinetic_rhs_modal(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  gkyl_dg_updater_lbo_pkpm_advance(lbo->coll_slvr, &species->local, fin, species->cflrate_f, rhs);
}

// The expansion-based operator (GPU runs): the collision frequency of the input, the boundary
// corrections, the Galerkin solve for the primitive moments and the drag and diffusion updaters.
static void
pkpm_species_lbo_init_modal(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
)
{
  int vdim = app->vdim;
  double v_bounds[2 * GKYL_MAX_DIM];
  for (int d = 0; d < vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }

  // The kernels of the operator read the expansions of nu and of nu [u, vtsq].
  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_prim_moms =
    mkarr(app->use_gpu, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);

  // allocate nu and initialize it (the computed frequency is not available in this form)
  assert(s->info.collisions.self_nu);
  double nu_frac = s->info.collisions.nuFrac > 0.0 ? s->info.collisions.nuFrac : 1.0;
  lbo->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = pkpm_species_coll_profile_new(
    app, s->info.collisions.self_nu, s->info.collisions.ctx, nu_frac, false
  );
  gkyl_array_copy(lbo->self_nu, self_nu);
  gkyl_array_copy(lbo->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  lbo->boundary_corrections =
    mkarr(app->use_gpu, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);

  // edge of velocity space corrections to momentum and energy
  // Note: PKPM model still has a momentum correction to insure M1 = 0 from collisions
  // The collision operator acts on the normalized distributions on the w grid: the
  // boundary corrections and primitive moments are formed from the mass-1 normalized
  // moments (self collisions: u = 0 and vth^2 = 1 up to the corrections).
  lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_pkpm_new(
    &s->grid, &app->confBasis, &app->basis, v_bounds, 1.0, app->use_gpu
  );
  // primitive moment calculators
  lbo->coll_pcalc =
    gkyl_prim_lbo_pkpm_calc_new(&s->grid, &app->confBasis, &app->basis, &app->local, app->use_gpu);
  lbo->calc_moms = pkpm_species_lbo_moms_modal;

  // LBO updater
  struct gkyl_dg_lbo_pkpm_drag_auxfields drag_inp = {
    .nuSum = lbo->nu_sum,
    .nuPrimMomsSum = lbo->nu_prim_moms,
  };
  struct gkyl_dg_lbo_pkpm_diff_auxfields diff_inp = {
    .nuSum = lbo->nu_sum,
    .nuPrimMomsSum = lbo->nu_prim_moms,
  };
  lbo->coll_slvr = gkyl_dg_updater_lbo_pkpm_new(
    &s->grid, &app->confBasis, &app->basis, &app->local, &drag_inp, &diff_inp, app->use_gpu
  );
  lbo->kinetic_rhs = pkpm_species_lbo_kinetic_rhs_modal;
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
  int vdim = app->vdim;

  // Primitive moments [u, vtsq]: node values, or expansions on GPU runs.
  lbo->prim_moms =
    mkarr(app->use_gpu, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);

  // GPU runs take the expansion-based moments.
  lbo->use_nodal = !app->use_gpu;
  if (lbo->use_nodal) {
    pkpm_species_lbo_init_nodal(app, s, lbo);
  } else {
    pkpm_species_lbo_init_modal(app, s, lbo);
  }
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

  if (lbo->use_nodal) {
    gkyl_array_release(lbo->nu);
    gkyl_array_release(lbo->nu_self);
    gkyl_pkpm_lbo_moms_release(lbo->up_moms);
    gkyl_dg_updater_lbo_pkpm_nodal_release(lbo->up_coll);
  } else {
    gkyl_array_release(lbo->nu_sum);
    gkyl_array_release(lbo->nu_prim_moms);
    gkyl_dg_updater_lbo_pkpm_release(lbo->coll_slvr);
    gkyl_array_release(lbo->boundary_corrections);
    gkyl_array_release(lbo->self_nu);
    gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
    gkyl_prim_lbo_calc_release(lbo->coll_pcalc);
  }
}
