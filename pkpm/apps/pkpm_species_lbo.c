#include <assert.h>
#include <gkyl_pkpm_priv.h>

// Boundary corrections, primitive moments and collision frequency as p=1 expansions (GPU runs).
static void
pkpm_species_lbo_moms_modal(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  gkyl_array_set_range(lbo->m0, 1.0, species->pkpm_rho, &app->local);

  // construct boundary corrections
  gkyl_mom_calc_bcorr_advance(
    lbo->bcorr_calc, &species->local, &app->local, fin, lbo->boundary_corrections
  );

  // Compute primitive moments of the normalized distributions.
  gkyl_prim_lbo_calc_advance(
    lbo->coll_pcalc, &app->local, species->pkpm_nmoms.marr, lbo->boundary_corrections, lbo->self_nu,
    lbo->prim_moms
  );

  // PKPM moments already computed before this, so just fetch results
  if (lbo->normNu) {
    gkyl_array_clear(lbo->norm_nu, 0.0);
    gkyl_array_clear(lbo->self_nu, 0.0);
    // Get density information (and scale out mass factor in PKPM moments)
    gkyl_array_set_range(
      lbo->maxwellian_moms, 1.0 / species->info.mass, species->pkpm_rho, &app->local
    );
    gkyl_array_set_offset_range(
      lbo->maxwellian_moms, 1.0, lbo->prim_moms, app->confBasis.num_basis, &app->local
    );
    // physical thermal speed squared for the Spitzer frequency: xi^2 times the normalized one
    gkyl_dg_mul_op_range(
      &app->confBasis, 0, species->pkpm_scratch, app->vdim, lbo->prim_moms, 1, species->pkpm_xi,
      &app->local
    );
    gkyl_array_set_offset_range(
      lbo->maxwellian_moms, 1.0, species->pkpm_scratch, (1 + app->vdim) * app->confBasis.num_basis,
      &app->local
    );
    gkyl_spitzer_coll_freq_advance_normnu(
      lbo->spitzer_calc, &app->local, lbo->maxwellian_moms, 0., lbo->maxwellian_moms, 0., 1.0,
      lbo->norm_nu
    );
    gkyl_dg_mul_op(&app->confBasis, 0, lbo->self_nu, 0, lbo->nu_init, 0, lbo->norm_nu);
  }
  for (int d = 0; d < app->vdim; d++) {
    gkyl_dg_mul_op(&app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
  }
  gkyl_dg_mul_op(
    &app->confBasis, app->vdim, lbo->nu_prim_moms, app->vdim, lbo->prim_moms, 0, lbo->self_nu
  );
}

// Primitive moments at the nodes with the collision frequency of the input.
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

// Primitive moments at the nodes with the Spitzer scaling of the collision frequency,
// nu_ref n / (2 vtsq)^(3/2) with the physical vtsq = xi^2 times the normalized one.
static void
pkpm_species_lbo_moms_nodal_spitzer(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  gkyl_pkpm_lbo_moms_advance(
    lbo->up_moms, &app->local, &species->local, species->pkpm_nmoms.marr, species->f_nodal,
    lbo->prim_moms
  );

  int num_nodes = app->confBasis.num_basis;
  gkyl_dg_conf_nodal_conf_to_nodal(app->conf_nodal, &app->local, species->pkpm_rho, lbo->n_nodal);
  gkyl_array_scale_range(lbo->n_nodal, 1.0 / species->info.mass, &app->local);
  gkyl_dg_conf_nodal_conf_to_nodal(app->conf_nodal, &app->local, species->pkpm_xi, lbo->xi_nodal);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);
    const double *xi_sq = (const double *)gkyl_array_cfetch(lbo->xi_nodal, loc) + num_nodes;
    const double *vtsq_norm = (const double *)gkyl_array_cfetch(lbo->prim_moms, loc) + num_nodes;
    double *vtsq = gkyl_array_fetch(lbo->vtsq_nodal, loc);
    for (int q = 0; q < num_nodes; ++q) {
      vtsq[q] = xi_sq[q] * vtsq_norm[q];
    }
  }
  gkyl_pkpm_coll_nu_advance(
    lbo->up_nu, &app->local, lbo->nu_ref, lbo->n_nodal, lbo->vtsq_nodal, lbo->vtsq_nodal, lbo->nu
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

// The expansion-based operator (GPU runs): the collision frequency, the boundary corrections,
// the Galerkin solve for the primitive moments and the drag and diffusion updaters.
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

  // allocate nu and initialize it
  lbo->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = pkpm_species_coll_profile_new(
    app, s->info.collisions.self_nu, s->info.collisions.ctx, 1.0, false
  );
  gkyl_array_copy(lbo->self_nu, self_nu);
  gkyl_array_copy(lbo->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  if (lbo->normNu) {
    double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
    double eps0 = app->field->info.epsilon0 ? app->field->info.epsilon0 : 1.0;
    double hbar = s->info.collisions.hbar ? s->info.collisions.hbar : 1.0;
    lbo->spitzer_calc = gkyl_spitzer_coll_freq_new(
      &app->confBasis, app->poly_order + 1, nuFrac, eps0, hbar, app->use_gpu
    );
    // Create arrays for scaling collisionality by normalization factor
    // norm_nu is computed from Spitzer calc and is the normalization factor for the local
    // density and thermal velocity, norm_nu_sr = n/(vth_s^2 + vth_r^2)^(3/2)
    // nu_init is the inital collisionality profile, which must be stored so that at every time
    // time step the collisionality profile is properly scaled and the effects are not cumulative
    lbo->norm_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->nu_init = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_copy(lbo->nu_init, lbo->self_nu);
    lbo->maxwellian_moms =
      mkarr(app->use_gpu, (vdim + 2) * app->confBasis.num_basis, app->local_ext.volume);
  }
  lbo->boundary_corrections =
    mkarr(app->use_gpu, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);
  lbo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

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
// input, scaled by nuFrac, or its Spitzer scaling), the updater for the primitive moments and the
// updater for the operator.
static void
pkpm_species_lbo_init_nodal(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
)
{
  int num_nodes = app->confBasis.num_basis;
  double nu_frac = s->info.collisions.nuFrac > 0.0 ? s->info.collisions.nuFrac : 1.0;
  lbo->nu_ref = pkpm_species_coll_profile_new(
    app, s->info.collisions.self_nu, s->info.collisions.ctx, nu_frac, true
  );
  lbo->nu = mkarr(false, num_nodes, app->local_ext.volume);
  gkyl_array_copy(lbo->nu, lbo->nu_ref);
  lbo->up_nu = gkyl_pkpm_coll_nu_new(&app->confBasis, lbo->normNu, 0.0, false);

  double vel_bounds[2] = {s->info.lower[0], s->info.upper[0]};
  lbo->up_moms = gkyl_pkpm_lbo_moms_new(&app->confBasis, &app->basis, vel_bounds, false);
  struct gkyl_dg_lbo_pkpm_nodal_auxfields nodal_inp = {.nu = lbo->nu, .prim_moms = lbo->prim_moms};
  lbo->up_coll = gkyl_dg_updater_lbo_pkpm_nodal_new(
    &s->grid, &app->confBasis, &app->basis, &app->local, &nodal_inp, false
  );
  lbo->kinetic_rhs = pkpm_species_lbo_kinetic_rhs_nodal;

  if (lbo->normNu) {
    lbo->n_nodal = mkarr(false, num_nodes, app->local_ext.volume);
    lbo->xi_nodal = mkarr(false, s->pkpm_xi->ncomp, app->local_ext.volume);
    lbo->vtsq_nodal = mkarr(false, num_nodes, app->local_ext.volume);
    lbo->calc_moms = pkpm_species_lbo_moms_nodal_spitzer;
  } else {
    lbo->calc_moms = pkpm_species_lbo_moms_nodal;
  }
}

void
pkpm_species_lbo_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
)
{
  int vdim = app->vdim;
  lbo->num_cross_collisions = s->info.collisions.num_cross_collisions;
  lbo->normNu = s->info.collisions.normNu;
  lbo->spitzer_calc = 0;

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
pkpm_species_lbo_cross_init(
  struct gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_lbo_collisions *lbo
)
{
  int vdim = app->vdim;

  lbo->cross_calc = gkyl_prim_lbo_vlasov_cross_calc_new(
    &s->grid, &app->confBasis, &app->basis, &app->local, app->use_gpu
  );

  lbo->cross_nu_prim_moms =
    mkarr(app->use_gpu, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);

  lbo->greene_factor_mem = 0;
  if (app->use_gpu) {
    lbo->greene_factor_mem =
      gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
  } else {
    lbo->greene_factor_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
  }

  // Cross-species collisions are not supported in the thermal-velocity normalized
  // formulation yet (the other species' primitive moments would need the frame change).
  assert(lbo->num_cross_collisions == 0);
  // set pointers to species we cross-collide with
  for (int i = 0; i < lbo->num_cross_collisions; ++i) {
    lbo->collide_with[i] = pkpm_find_species(app, s->info.collisions.collide_with[i]);
    lbo->other_m[i] = lbo->collide_with[i]->info.mass;
    lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
    lbo->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_prim_moms[i] =
      mkarr(app->use_gpu, (vdim + 1) * app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_num[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_den[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_factor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    if (lbo->other_m[i] > s->info.mass) {
      gkyl_array_set(lbo->cross_nu[i], sqrt(2), lbo->self_nu);
      gkyl_array_set(lbo->other_nu[i], (s->info.mass) / (lbo->other_m[i]), lbo->self_nu);
    } else {
      gkyl_array_set(
        lbo->cross_nu[i], (lbo->other_m[i]) / (s->info.mass), lbo->collide_with[i]->lbo.self_nu
      );
      gkyl_array_set(lbo->other_nu[i], sqrt(2), lbo->collide_with[i]->lbo.self_nu);
    }

    gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

    lbo->other_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->other_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    lbo->self_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->self_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    gkyl_array_set(lbo->self_mnu[i], s->info.mass, lbo->cross_nu[i]);
    gkyl_array_set(lbo->other_mnu[i], lbo->other_m[i], lbo->other_nu[i]);
  }

  lbo->betaGreenep1 = 1.0;
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

// computes moments from cross-species collisions
void
pkpm_species_lbo_cross_moms(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct pkpm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  struct timespec wst = gkyl_wall_clock();

  wst = gkyl_wall_clock();
  for (int i = 0; i < lbo->num_cross_collisions; ++i) {
    gkyl_dg_mul_op_range(
      &app->confBasis, 0, lbo->self_mnu_m0[i], 0, lbo->self_mnu[i], 0, lbo->m0, &app->local
    );
    gkyl_dg_mul_op_range(
      &app->confBasis, 0, lbo->other_mnu_m0[i], 0, lbo->other_mnu[i], 0,
      lbo->collide_with[i]->lbo.m0, &app->local
    );

    gkyl_dg_mul_op_range(
      &app->confBasis, 0, lbo->greene_num[i], 0, lbo->other_mnu_m0[i], 0, lbo->m0, &app->local
    );

    gkyl_array_set(lbo->greene_den[i], 1.0, lbo->self_mnu_m0[i]);
    gkyl_array_accumulate(lbo->greene_den[i], 1.0, lbo->other_mnu_m0[i]);

    gkyl_dg_div_op_range(
      lbo->greene_factor_mem, &app->confBasis, 0, lbo->greene_factor[i], 0, lbo->greene_num[i], 0,
      lbo->greene_den[i], &app->local
    );
    gkyl_array_scale(lbo->greene_factor[i], 2 * lbo->betaGreenep1);

    gkyl_prim_lbo_cross_calc_advance(
      lbo->cross_calc, &app->local, lbo->greene_factor[i], species->info.mass, lbo->moms.marr,
      lbo->prim_moms, lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
      lbo->boundary_corrections, lbo->cross_nu[i], lbo->cross_prim_moms[i]
    );

    for (int d = 0; d < app->vdim; d++) {
      gkyl_dg_mul_op(
        &app->confBasis, d, lbo->cross_nu_prim_moms, d, lbo->cross_prim_moms[i], 0, lbo->cross_nu[i]
      );
    }
    gkyl_dg_mul_op(
      &app->confBasis, app->vdim, lbo->cross_nu_prim_moms, app->vdim, lbo->cross_prim_moms[i], 0,
      lbo->cross_nu[i]
    );

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, lbo->cross_nu_prim_moms);
  }
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
    gkyl_array_release(lbo->nu_ref);
    gkyl_array_release(lbo->nu);
    gkyl_pkpm_coll_nu_release(lbo->up_nu);
    gkyl_pkpm_lbo_moms_release(lbo->up_moms);
    gkyl_dg_updater_lbo_pkpm_nodal_release(lbo->up_coll);
    if (lbo->normNu) {
      gkyl_array_release(lbo->n_nodal);
      gkyl_array_release(lbo->xi_nodal);
      gkyl_array_release(lbo->vtsq_nodal);
    }
  } else {
    gkyl_array_release(lbo->nu_sum);
    gkyl_array_release(lbo->nu_prim_moms);
    gkyl_dg_updater_lbo_pkpm_release(lbo->coll_slvr);
    gkyl_array_release(lbo->boundary_corrections);
    gkyl_array_release(lbo->self_nu);
    gkyl_array_release(lbo->m0);
    gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
    gkyl_prim_lbo_calc_release(lbo->coll_pcalc);
    if (lbo->normNu) {
      gkyl_array_release(lbo->norm_nu);
      gkyl_array_release(lbo->nu_init);
      gkyl_spitzer_coll_freq_release(lbo->spitzer_calc);
      gkyl_array_release(lbo->maxwellian_moms);
    }
  }

  if (lbo->num_cross_collisions) {
    gkyl_dg_bin_op_mem_release(lbo->greene_factor_mem);
    gkyl_array_release(lbo->cross_nu_prim_moms);
    for (int i = 0; i < lbo->num_cross_collisions; ++i) {
      gkyl_array_release(lbo->cross_prim_moms[i]);
      gkyl_array_release(lbo->cross_nu[i]);
      gkyl_array_release(lbo->other_nu[i]);
      gkyl_array_release(lbo->self_mnu[i]);
      gkyl_array_release(lbo->self_mnu_m0[i]);
      gkyl_array_release(lbo->other_mnu[i]);
      gkyl_array_release(lbo->other_mnu_m0[i]);
      gkyl_array_release(lbo->greene_num[i]);
      gkyl_array_release(lbo->greene_den[i]);
      gkyl_array_release(lbo->greene_factor[i]);
    }
    gkyl_prim_lbo_cross_calc_release(lbo->cross_calc);
  }
}
