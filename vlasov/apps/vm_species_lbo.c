#include <assert.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_const.h>

static void
vmlbo_moms_disabled(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  // Empty method.
}

static void
vmlbo_moms_enabled(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute LTE moments (n, V_drift, T/m) if a computed (Spitzer) collision
  // frequency needs them: this species' own self/cross nu, or a partner's
  // computed cross nu, which reads this species' LTE moments (see cross_init).
  if (lbo->needs_lte_moms) {
    vm_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  }

  // Compute M0, M1i, M2 moments .
  vm_species_moment_calc(&lbo->moms, species->local, app->local, fin);

  // Construct boundary corrections.
  gkyl_mom_calc_bcorr_advance(
    lbo->bcorr_calc, &species->local, &app->local, fin, lbo->boundary_corrections
  );

  // Calculate nu_ss (and multibly moms and corrections by it for norm_nu).
  lbo->self_nu_func(app, species, lbo, fin);

  // Multiply M0, M1, M2 moments and boundary corrections by self nu.
  for (int d = 0; d < app->vdim + 2; d++) {
    gkyl_dg_mul_op(&app->basis, d, lbo->nu_moms, d, lbo->moms.marr, 0, lbo->self_nu);
  }
  for (int d = 0; d < app->vdim + 1; d++) {
    gkyl_dg_mul_op(
      &app->basis, d, lbo->nu_boundary_corrections, d, lbo->boundary_corrections, 0, lbo->self_nu
    );
  }

  // Construct primitive moments.
  gkyl_prim_lbo_calc_advance(
    lbo->coll_pcalc, &app->local, lbo->nu_moms, lbo->nu_boundary_corrections, lbo->self_nu,
    lbo->prim_moms
  );

  // Scale upar and vtSq by self nu.
  for (int d = 0; d < app->vdim + 1; d++) {
    gkyl_dg_mul_op(&app->basis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
  }

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

// Every stage rebuilds nu_sum = self_nu + sum_r cross_nu_sr: the self method
// (const or computed) resets it and each cross method accumulates onto it, so
// the four self/cross mode combinations compose correctly.
static void
vmlbo_self_nu_calc_constNu(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  // self_nu is static (projected at init); reset nu_sum from it so the cross
  // methods accumulate onto a fresh sum each stage.
  gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);
}

static void
vmlbo_self_nu_calc_normNu(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  // Calculate nu_ss(x,t) (LTE moments are staged by vmlbo_moms_enabled).
  gkyl_spitzer_coll_freq_advance_normnu(
    lbo->spitzer_calc, &app->local, species->lte.moms.marr, lbo->vtsq_min, species->lte.moms.marr,
    lbo->vtsq_min, lbo->norm_nu_fac_self, lbo->self_nu
  );

  gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);
}

static void
vmlbo_cross_nu_calc_constNu(
  gkyl_vlasov_app *app, const struct vm_species *s, struct vm_lbo_collisions *lbo, int coll_idx
)
{
  // cross_nu is static (projected at init); accumulate onto the freshly-reset nu_sum.
  gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[coll_idx]);
}

static void
vmlbo_cross_nu_calc_normNu(
  gkyl_vlasov_app *app, const struct vm_species *s, struct vm_lbo_collisions *lbo, int coll_idx
)
{
  // Calculate nu_sr(x,t).
  gkyl_spitzer_coll_freq_advance_normnu(
    lbo->spitzer_calc, &app->local, s->lte.moms.marr, lbo->vtsq_min,
    lbo->collide_with[coll_idx]->lte.moms.marr, lbo->collide_with[coll_idx]->lbo.vtsq_min,
    lbo->norm_nu_fac_cross[coll_idx], lbo->cross_nu[coll_idx]
  );

  gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[coll_idx]);
}

static void
vmlbo_alpha_E_constNu(
  gkyl_vlasov_app *app, const struct vm_species *s, struct vm_lbo_collisions *lbo, int coll_idx
)
{
  gkyl_array_clear(lbo->alpha_E, 0.0);
  gkyl_array_shiftc(lbo->alpha_E, lbo->alpha_E_fac[coll_idx], 0);
}

static void
vmlbo_alpha_E_normNu(
  gkyl_vlasov_app *app, const struct vm_species *s, struct vm_lbo_collisions *lbo, int coll_idx
)
{
  gkyl_dg_mul_op_range(
    &app->basis, 0, lbo->alpha_E, 0, lbo->cross_nu[coll_idx], 0, s->lte.moms.marr, &app->local
  );
  gkyl_array_scale_range(lbo->alpha_E, lbo->alpha_E_fac[coll_idx], &app->local);
}

static void
vmlbo_cross_moms_disabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo
)
{
  // Empty method.
}

static void
vmlbo_cross_moms_enabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo
)
{
  // Compute primitive moments for cross-species collisions.
  struct timespec wst = gkyl_wall_clock();

  for (int i = 0; i < lbo->num_cross_collisions; ++i) {
    // Compute the cross-species collision frequency.
    lbo->cross_nu_func(app, vms, lbo, i);

    // Compute alpha_E.
    lbo->alpha_E_func(app, vms, lbo, i);

    // Multiply moments and boundary corrections by cross nu.
    for (int d = 0; d < app->vdim + 2; d++) {
      gkyl_dg_mul_op(&app->basis, d, lbo->nu_moms, d, lbo->moms.marr, 0, lbo->cross_nu[i]);
    }
    for (int d = 0; d < app->vdim + 1; d++) {
      gkyl_dg_mul_op(
        &app->basis, d, lbo->nu_boundary_corrections, d, lbo->boundary_corrections, 0,
        lbo->cross_nu[i]
      );
    }

    // Compute cross primitive moments.
    // Recycle the boundary_corrections array because we don't need those anymore.
    struct gkyl_array *cross_prim_moms = lbo->nu_boundary_corrections;
    gkyl_prim_lbo_cross_calc_advance(
      lbo->cross_calc, &app->local, lbo->alpha_E, vms->mass, lbo->nu_moms, lbo->prim_moms,
      lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
      lbo->nu_boundary_corrections, lbo->cross_nu[i], cross_prim_moms
    );

    // Scale u_{sr} and vtSq_{sr} by nu_{sr}.
    for (int d = 0; d < app->vdim + 1; d++) {
      gkyl_dg_mul_op(&app->basis, d, cross_prim_moms, d, cross_prim_moms, 0, lbo->cross_nu[i]);
    }

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, cross_prim_moms);
  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
vmlbo_rhs_disabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  // Empty method.
}

static void
vmlbo_rhs_enabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();

  // Accumulate update due to collisions onto rhs.
  gkyl_dg_updater_lbo_vlasov_advance(lbo->coll_slvr, &vms->local, fin, vms->cflrate, rhs);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
vmlbo_write_mom_disabled(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame)
{
  // Empty method.
}

static void
vmlbo_write_mom_enabled(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wtm = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });

  // Write out nu_sum and nu_prim_moms.
  const char *fmt = "%s-%s_lbo_nu_sum_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);

  const char *fmt_nu_prim = "%s-%s_lbo_nu_prim_moms_%d.gkyl";
  int sz_nu_prim = gkyl_calc_strlen(fmt_nu_prim, app->name, vms->name, frame);
  char fileNm_nu_prim[sz_nu_prim + 1]; // ensures no buffer overflow
  snprintf(fileNm_nu_prim, sizeof fileNm_nu_prim, fmt_nu_prim, app->name, vms->name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(vms->lbo.nu_sum_host, vms->lbo.nu_sum);
    gkyl_array_copy(vms->lbo.nu_prim_moms_host, vms->lbo.nu_prim_moms);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, vms->lbo.nu_sum_host, fileNm);
  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, vms->lbo.nu_prim_moms_host, fileNm_nu_prim
  );
  app->stat.n_diag_io += 2;

  vlasov_array_meta_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void
vm_species_lbo_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_lbo_collisions *lbo
)
{
  lbo->collision_id = vms->info.collisions.collision_id;
  lbo->write_coll_diagnostics = vms->info.collisions.write_coll_diagnostics;

  // Empty methods.
  lbo->moms_func = vmlbo_moms_disabled;
  lbo->rhs_func = vmlbo_rhs_disabled;
  lbo->write_mom_func = vmlbo_write_mom_disabled;

  if (lbo->collision_id == GKYL_LBO_COLLISIONS) {
    lbo->num_cross_collisions = vms->info.collisions.num_cross_collisions;

    int cdim = app->cdim, vdim = app->vdim;
    double v_bounds[2 * GKYL_MAX_DIM] = {0.0};
    for (int d = 0; d < vdim; ++d) {
      v_bounds[d] = vms->info.lower[d];
      v_bounds[d + vdim] = vms->info.upper[d];
    }
    // Allocate self-species collision frequency and sum of collision frequencies.
    lbo->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    lbo->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    double nu_frac = vms->info.collisions.nu_frac ? vms->info.collisions.nu_frac : 1.0;

    // The Spitzer collision-frequency updater (and the vtsq_min floor it needs)
    // is required whenever any collision frequency is computed in time: a computed
    // self_nu, or a computed cross_nu (cross collisions for which no explicit
    // cross_nu was provided). It is independent of the self_nu mode, so a species
    // may use an explicit self_nu together with a computed cross_nu.
    bool computed_cross = vms->info.collisions.num_cross_collisions > 0 &&
                          !vms->info.collisions.cross_nu[0];
    if (!vms->info.collisions.self_nu || computed_cross) {
      // Compute a minimum representable temperature based on the smallest dv in the grid.
      double vtsq_min = 0.0;
      for (int d = 0; d < vdim; ++d) {
        vtsq_min += (1.0 / 6.0) * pow(vms->grid.dx[cdim + d], 2);
      }
      lbo->vtsq_min = vtsq_min / vdim;

      lbo->spitzer_calc =
        gkyl_spitzer_coll_freq_new(&app->basis, app->poly_order + 1, 1.0, 1.0, 1.0, app->use_gpu);
    }

    if (vms->info.collisions.self_nu) {
      // Project user's self-species collision frequency.
      lbo->norm_nu_self = false;

      struct gkyl_array *self_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(
        &app->grid, &app->basis, app->poly_order + 1, 1, vms->info.collisions.self_nu,
        vms->info.collisions.self_nu_ctx
      );
      gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu_ho);
      gkyl_proj_on_basis_release(proj);
      gkyl_array_copy(lbo->self_nu, self_nu_ho);
      gkyl_array_release(self_nu_ho);

      gkyl_array_scale(lbo->self_nu, nu_frac);
      gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);

      // Set pointers to functions chosen at runtime.
      lbo->self_nu_func = vmlbo_self_nu_calc_constNu;
    } else {
      // Self-collision frequency computed in time.
      lbo->norm_nu_self = true;

      double eps0 = vms->info.collisions.eps0 ? vms->info.collisions.eps0 : GKYL_EPSILON0;
      double hbar = vms->info.collisions.hbar ? vms->info.collisions.hbar :
                                                GKYL_PLANCKS_CONSTANT_H / 2 / M_PI;
      double eV = vms->info.collisions.eV ? vms->info.collisions.eV : GKYL_ELEMENTARY_CHARGE;
      // Vlasov does not use reference magnetic field for cyclotron frequency contribution to log(Lambda)
      double bmag_ref = 0.0;
      // vtsq_min and spitzer_calc are set up above (needed for any computed frequency).

      // We define nu_ss = nu_sr(r=s) = alpha_E/((delta_ss * (1+beta))*n_s), with delta_ss = 2,
      // beta = 0. This gives a nu_ss that is arguably 2X smaller than it should be, but it's
      // cheaper and yields an electron isotropization rate that agrees better with the FPO's.
      lbo->norm_nu_fac_self = nu_frac * gkyl_calc_Morse_alpha_E_const(
                                          vms->info.collisions.den_ref,
                                          vms->info.collisions.den_ref, vms->mass, vms->mass,
                                          vms->charge, vms->charge, vms->info.collisions.temp_ref,
                                          vms->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV
                                        );

      // Set pointers to functions chosen at runtime.
      lbo->self_nu_func = vmlbo_self_nu_calc_normNu;
      lbo->needs_lte_moms = true; // Spitzer self nu reads this species' LTE moments.
    }

    // Create moment calculator to get M0, M1, M2 for primitive moments.
    vm_species_moment_init(app, vms, &lbo->moms, GKYL_F_MOMENT_M0M1M2, false);
    lbo->nu_moms = mkarr(app->use_gpu, lbo->moms.marr->ncomp, lbo->moms.marr->size);

    // Edge of velocity space corrections to momentum and energy.
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_new(
      &vms->grid, &app->basis, &vms->basis, v_bounds, app->use_gpu
    );

    // Primitive moment calculator.
    lbo->coll_pcalc = gkyl_prim_lbo_vlasov_calc_new(
      &vms->grid, &app->basis, &vms->basis, &app->local, app->use_gpu
    );

    // Allocate boundary corrections for primitive mom calculation.
    lbo->boundary_corrections =
      mkarr(app->use_gpu, (vdim + 1) * app->basis.num_basis, app->local_ext.volume);
    lbo->nu_boundary_corrections =
      mkarr(app->use_gpu, lbo->boundary_corrections->ncomp, lbo->boundary_corrections->size);

    // Primitive moments.
    lbo->prim_moms = mkarr(app->use_gpu, (vdim + 1) * app->basis.num_basis, app->local_ext.volume);
    lbo->nu_prim_moms = mkarr(app->use_gpu, lbo->prim_moms->ncomp, lbo->prim_moms->size);

    // Host-side copy for I/O.
    if (lbo->write_coll_diagnostics) {
      if (app->use_gpu) {
        lbo->nu_sum_host = mkarr(false, lbo->nu_sum->ncomp, lbo->nu_sum->size);
        lbo->nu_prim_moms_host = mkarr(false, lbo->nu_prim_moms->ncomp, lbo->nu_prim_moms->size);
      } else {
        lbo->nu_sum_host = lbo->nu_sum;
        lbo->nu_prim_moms_host = lbo->nu_prim_moms;
      }
    }

    // LBO updater.
    struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = {
      .nuSum = lbo->nu_sum,
      .nuPrimMomsSum = lbo->nu_prim_moms,
    };
    struct gkyl_dg_lbo_vlasov_diff_auxfields diff_inp = {
      .nuSum = lbo->nu_sum,
      .nuPrimMomsSum = lbo->nu_prim_moms,
    };
    lbo->coll_slvr = gkyl_dg_updater_lbo_vlasov_new(
      &vms->grid, &app->basis, &vms->basis, &app->local, &drag_inp, &diff_inp, app->use_gpu
    );

    // Methods chosen at runtime.
    lbo->moms_func = vmlbo_moms_enabled;
    lbo->rhs_func = vmlbo_rhs_enabled;
    if (lbo->write_coll_diagnostics) {
      lbo->write_mom_func = vmlbo_write_mom_enabled;
    }
  }
}

void
vm_species_lbo_cross_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_lbo_collisions *lbo
)
{
  // Empty methods.
  lbo->cross_nu_func = vmlbo_cross_nu_calc_constNu;
  lbo->cross_moms_func = vmlbo_cross_moms_disabled;

  if (lbo->collision_id == GKYL_LBO_COLLISIONS) {
    if (vms->lbo.num_cross_collisions) {
      lbo->betaGreenep1 = 1.0; // Greene's beta factor + 1.
      lbo->delta_sr = 2.0; // delta_sr free parameter.

      // Set pointers to species we cross-collide with.
      int my_idx_in_other[GKYL_MAX_SPECIES];
      for (int i = 0; i < lbo->num_cross_collisions; ++i) {
        lbo->collide_with[i] = vm_find_species(app, vms->info.collisions.collide_with[i]);
        // collide_with must name an existing *kinetic* species (a typo, or a
        // fluid species, returns NULL and would segfault below without a message).
        assert(lbo->collide_with[i]);
        my_idx_in_other[i] = -1;
        for (int j = 0; j < lbo->collide_with[i]->lbo.num_cross_collisions; ++j) {
          if (0 == strcmp(vms->name, lbo->collide_with[i]->info.collisions.collide_with[j])) {
            my_idx_in_other[i] = j;
            break;
          }
        }
        // Cross collisions must be specified symmetrically: the partner must
        // list this species in its own collide_with (my_idx_in_other indexes the
        // partner's cross arrays below).
        assert(my_idx_in_other[i] >= 0);
      }

      // Morse's alpha_E.
      lbo->alpha_E = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      for (int i = 0; i < lbo->num_cross_collisions; ++i) {
        // Cross-species collision frequency, nu_sr.
        lbo->cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
        lbo->other_m[i] = lbo->collide_with[i]->mass;
        lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
      }

      double nu_frac = vms->info.collisions.nu_frac ? vms->info.collisions.nu_frac : 1.0;

      // Compute the time-independent part of alpha_E.
      double alpha_E_norm[GKYL_MAX_SPECIES] = {0.0};
      for (int i = 0; i < lbo->num_cross_collisions; ++i) {
        double eps0 = vms->info.collisions.eps0 ? vms->info.collisions.eps0 : GKYL_EPSILON0;
        double hbar = vms->info.collisions.hbar ? vms->info.collisions.hbar :
                                                  GKYL_PLANCKS_CONSTANT_H / 2 / M_PI;
        double eV = vms->info.collisions.eV ? vms->info.collisions.eV : GKYL_ELEMENTARY_CHARGE;
        // Vlasov does not use reference magnetic field for cyclotron frequency contribution to log(Lambda)
        double bmag_ref = 0.0;
        double mass_self = vms->mass, mass_other = lbo->collide_with[i]->mass;

        alpha_E_norm[i] =
          nu_frac * gkyl_calc_Morse_alpha_E_const(
                      vms->info.collisions.den_ref, lbo->collide_with[i]->info.collisions.den_ref,
                      mass_self, mass_other, vms->charge, lbo->collide_with[i]->charge,
                      vms->info.collisions.temp_ref, lbo->collide_with[i]->info.collisions.temp_ref,
                      bmag_ref, eps0, hbar, eV
                    );
      }

      if (vms->info.collisions.cross_nu[0]) {
        // Project user's cross-species collision frequency.
        lbo->norm_nu_cross = false;

        // Cross-collision frequency must be specified symmetrically: if this
        // species provides an explicit cross_nu, its partner must too. The self
        // collision frequency mode is independent (a species may use explicit
        // self_nu with explicit or computed cross_nu), so it is not checked here.
        for (int i = 0; i < lbo->num_cross_collisions; ++i) {
          assert(lbo->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]);
        }

        struct gkyl_array *cross_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
        for (int i = 0; i < lbo->num_cross_collisions; ++i) {
          gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(
            &app->grid, &app->basis, app->poly_order + 1, 1, vms->info.collisions.cross_nu[i],
            vms->info.collisions.cross_nu_ctx[i]
          );
          gkyl_proj_on_basis_advance(proj, 0.0, &app->local, cross_nu_ho);
          gkyl_proj_on_basis_release(proj);
          gkyl_array_copy(lbo->cross_nu[i], cross_nu_ho);

          gkyl_array_scale(lbo->cross_nu[i], nu_frac);
          gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

          // Compute alpha_E using reference parameters.
          assert(vms->info.collisions.den_ref);
          assert(lbo->collide_with[i]->info.collisions.den_ref);
          assert(vms->info.collisions.temp_ref);
          assert(lbo->collide_with[i]->info.collisions.temp_ref);
          double mass_self = vms->mass, mass_other = lbo->collide_with[i]->mass;
          double den_s = vms->info.collisions.den_ref;
          double den_r = lbo->collide_with[i]->info.collisions.den_ref;
          double vtsq_s = vms->info.collisions.temp_ref / mass_self;
          double vtsq_r = lbo->collide_with[i]->info.collisions.temp_ref / mass_other;

          lbo->alpha_E_fac[i] =
            (alpha_E_norm[i] * den_s * den_r / pow(sqrt(vtsq_s + vtsq_r), 3.0)) *
            pow(sqrt(2.0), app->cdim);
        }
        gkyl_array_release(cross_nu_ho);

        // Set pointers to functions chosen at runtime.
        lbo->cross_nu_func = vmlbo_cross_nu_calc_constNu;
        lbo->alpha_E_func = vmlbo_alpha_E_constNu;
      } else {
        // Cross-collision frequency computed in time.
        lbo->norm_nu_cross = true;

        // Cross-collision frequency must be specified symmetrically: if this
        // species uses a computed (Spitzer) cross_nu, its partner must too (i.e.
        // the partner must not provide an explicit cross_nu). The self collision
        // frequency mode is independent and is not checked here.
        for (int i = 0; i < lbo->num_cross_collisions; ++i) {
          assert(!(lbo->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]));
        }

        for (int i = 0; i < lbo->num_cross_collisions; ++i) {
          double mass_self = vms->mass, mass_other = lbo->collide_with[i]->mass;

          lbo->norm_nu_fac_cross[i] = alpha_E_norm[i] * (mass_self + mass_other) /
                                      (lbo->delta_sr * lbo->betaGreenep1 * mass_self);

          lbo->alpha_E_fac[i] =
            (lbo->delta_sr * lbo->betaGreenep1 * mass_self) / (mass_self + mass_other);
        }

        // Set pointers to functions chosen at runtime.
        lbo->cross_nu_func = vmlbo_cross_nu_calc_normNu;
        lbo->alpha_E_func = vmlbo_alpha_E_normNu;

        // The Spitzer cross nu reads both this species' and each partner's LTE
        // moments, so mark both for the per-stage LTE moment calculation in
        // vmlbo_moms_enabled. (Safe here: all lbo_init calls complete before
        // any cross_init runs, and flags are only ever set, never cleared.)
        lbo->needs_lte_moms = true;
        for (int i = 0; i < lbo->num_cross_collisions; ++i) {
          lbo->collide_with[i]->lbo.needs_lte_moms = true;
        }
      }

      // Cross-primitive moment calculator.
      lbo->cross_calc = gkyl_prim_lbo_vlasov_cross_calc_new(
        &vms->grid, &app->basis, &vms->basis, &app->local, app->use_gpu
      );

      // Methods chosen at runtime.
      lbo->cross_moms_func = vmlbo_cross_moms_enabled;
    }
  }
}

void
vm_species_lbo_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  lbo->moms_func(app, species, lbo, fin);
}

void
vm_species_lbo_cross_moms(
  gkyl_vlasov_app *app, const struct vm_species *species, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin
)
{
  lbo->cross_moms_func(app, species, lbo);
}

void
vm_species_lbo_rhs(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  lbo->rhs_func(app, vms, lbo, fin, rhs);
}

void
vm_species_lbo_write_mom(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame)
{
  vms->lbo.write_mom_func(app, vms, tm, frame);
}

void
vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo)
{
  if (lbo->collision_id == GKYL_LBO_COLLISIONS) {
    if (lbo->num_cross_collisions) {
      gkyl_prim_lbo_cross_calc_release(lbo->cross_calc);

      for (int i = 0; i < lbo->num_cross_collisions; ++i) {
        gkyl_array_release(lbo->cross_nu[i]);
      }

      gkyl_array_release(lbo->alpha_E);
    }

    gkyl_dg_updater_lbo_vlasov_release(lbo->coll_slvr);

    if (lbo->write_coll_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(lbo->nu_sum_host);
        gkyl_array_release(lbo->nu_prim_moms_host);
      }
    }

    gkyl_array_release(lbo->nu_prim_moms);
    gkyl_array_release(lbo->prim_moms);

    gkyl_array_release(lbo->nu_boundary_corrections);
    gkyl_array_release(lbo->boundary_corrections);

    gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
    gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

    gkyl_array_release(lbo->nu_moms);
    vm_species_moment_release(app, &lbo->moms);

    // spitzer_calc was allocated if any frequency is computed (self or cross).
    if (lbo->norm_nu_self || lbo->norm_nu_cross) {
      gkyl_spitzer_coll_freq_release(lbo->spitzer_calc);
    }

    gkyl_array_release(lbo->nu_sum);
    gkyl_array_release(lbo->self_nu);
  }
}
