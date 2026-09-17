#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

static void
gklbo_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Empty method.
}

static void
gklbo_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute J*M0, J*M1, J*M2 moments and separate our M0 and M2.
  gk_species_moment_calc(&lbo->moms, species->local, app->local, fin);
  gkyl_array_set_offset_range(lbo->m2self, 1.0, lbo->moms.marr, 2*app->basis.num_basis, &app->local);
  
  // Construct boundary corrections.
  gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
    &species->local, &app->local, fin, lbo->boundary_corrections);

  // Calculate nu_ss (and multibly moms and corrections by it for norm_nu).
  lbo->self_nu_func(app, species, lbo, fin);

  // Multiply moments and boundary corrections by self nu.
  for (int d=0; d<3; d++)
    gkyl_dg_mul_op(&app->basis, d, lbo->nu_moms, d, lbo->moms.marr, 0, lbo->self_nu);
  for (int d=0; d<2; d++)
    gkyl_dg_mul_op(&app->basis, d, lbo->nu_boundary_corrections, d, lbo->boundary_corrections, 0, lbo->self_nu);

  // Construct primitive moments.
  gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local, 
    lbo->nu_moms, lbo->nu_boundary_corrections, lbo->self_nu, lbo->prim_moms);

  // Scale upar and vtSq by self nu.
  for (int d=0; d<2; d++)
    gkyl_dg_mul_op(&app->basis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
gklbo_self_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Empty method.
}

static void
gklbo_self_nu_calc_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Calculate nu_ss(x,t).
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, species->lte.moms.marr,
    0, species->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);

  gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, species->lte.moms.marr, lbo->vtsq_min,
    species->lte.moms.marr, lbo->vtsq_min, lbo->norm_nu_fac_self, lbo->self_nu);

  gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);
}

static void
gklbo_cross_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo, int coll_idx)
{
  // Empty method.
}

static void
gklbo_cross_nu_calc_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo, int coll_idx)
{
  // Calculate nu_sr(x,t).
  gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, s->lte.moms.marr, lbo->vtsq_min,
    lbo->collide_with[coll_idx]->lte.moms.marr, lbo->collide_with[coll_idx]->lbo.vtsq_min,
    lbo->norm_nu_fac_cross[coll_idx], lbo->cross_nu[coll_idx]);

  gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[coll_idx]);
}

static void
gklbo_alpha_E_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo, int coll_idx)
{
  gkyl_array_clear(lbo->alpha_E, 0.0);
  gkyl_array_shiftc(lbo->alpha_E, lbo->alpha_E_fac[coll_idx], 0);
}

static void
gklbo_alpha_E_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo, int coll_idx)
{
  gkyl_dg_mul_op_range(&app->basis, 0, lbo->alpha_E, 0, lbo->cross_nu[coll_idx], 0, s->lte.moms.marr, &app->local);
  gkyl_array_scale_range(lbo->alpha_E, lbo->alpha_E_fac[coll_idx], &app->local);
}

static void
gklbo_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo)
{
  // Empty method.
}

static void
gklbo_cross_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo)
{
  // Compute primitive moments for cross-species collisions.
  struct timespec wst = gkyl_wall_clock();
  
  for (int i=0; i<lbo->num_cross_collisions; ++i) {

    // Compute the cross-species collision frequency.
    lbo->cross_nu_func(app, gks, lbo, i);

    // Compute alpha_E.
    lbo->alpha_E_func(app, gks, lbo, i);

    // Multiply moments and boundary corrections by cross nu.
    for (int d=0; d<3; d++)
      gkyl_dg_mul_op(&app->basis, d, lbo->nu_moms, d, lbo->moms.marr, 0, lbo->cross_nu[i]);
    for (int d=0; d<2; d++)
      gkyl_dg_mul_op(&app->basis, d, lbo->nu_boundary_corrections, d, lbo->boundary_corrections, 0, lbo->cross_nu[i]);

    // Compute cross primitive moments.
    // Recycle the boundary_corrections array because we don't need those anymore.
    struct gkyl_array *cross_prim_moms = lbo->nu_boundary_corrections;
    gkyl_prim_lbo_cross_calc_advance(lbo->cross_calc, &app->local, lbo->alpha_E, 
      gks->info.mass, lbo->nu_moms, lbo->prim_moms,
      lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
      lbo->nu_boundary_corrections, lbo->cross_nu[i], cross_prim_moms);

    // Scale upar_{sr} and vtSq_{sr} by nu_{sr}.
    for (int d=0; d<2; d++)
      gkyl_dg_mul_op(&app->basis, d, cross_prim_moms, d, cross_prim_moms, 0, lbo->cross_nu[i]);

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, cross_prim_moms);

  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
gklbo_rhs_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Empty method.
}

static void
gklbo_rhs_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // Accumulate update due to collisions onto rhs.
  gkyl_dg_updater_lbo_gyrokinetic_advance(lbo->coll_slvr, &gks->local,
    fin, gks->cflrate, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
gklbo_write_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // Empty method.
}

static void
gklbo_write_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wtm = gkyl_wall_clock();
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_conf_len, gks->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_conf_len, gks->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc_nu_sum[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Sum of collision frequencies." }
  };
  int io_meta_nu_sum_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta_nu_sum[] = {gks->io_meta_conf, app->gk_geom->io_meta_basic, desc_nu_sum};
  struct gkyl_msgpack_data *mt_nu_sum = gkyl_msgpack_create_union(sizeof(io_meta_nu_sum_len)/sizeof(int), io_meta_nu_sum_len, io_meta_nu_sum);

  struct gkyl_msgpack_map_elem desc_nu_prim[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING,
      .cval = "Drift velocity and thermal speed squared, times collision frequency, summed over colliding species." }
  };
  int io_meta_nu_prim_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta_nu_prim[] = {gks->io_meta_conf, app->gk_geom->io_meta_basic, desc_nu_prim};
  struct gkyl_msgpack_data *mt_nu_prim = gkyl_msgpack_create_union(sizeof(io_meta_nu_prim_len)/sizeof(int), io_meta_nu_prim_len, io_meta_nu_prim);

  // Write out nu_sum and nu_prim_moms.
  const char *fmt = "%s-%s_lbo_nu_sum_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  const char *fmt_nu_prim = "%s-%s_lbo_nu_prim_moms_%d.gkyl";
  int sz_nu_prim = gkyl_calc_strlen(fmt_nu_prim, app->name, gks->info.name, frame);
  char fileNm_nu_prim[sz_nu_prim+1]; // ensures no buffer overflow
  snprintf(fileNm_nu_prim, sizeof fileNm_nu_prim, fmt_nu_prim, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gks->lbo.nu_sum_host, gks->lbo.nu_sum);
    gkyl_array_copy(gks->lbo.nu_prim_moms_host, gks->lbo.nu_prim_moms);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt_nu_sum, gks->lbo.nu_sum_host, fileNm);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt_nu_prim, gks->lbo.nu_prim_moms_host, fileNm_nu_prim);
  app->stat.n_diag_io += 2;

  gkyl_msgpack_data_release(mt_nu_sum);
  gkyl_msgpack_data_release(mt_nu_prim);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void 
gk_species_lbo_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_lbo_collisions *lbo)
{
  lbo->collision_id = gks->info.collisions.collision_id;
  lbo->write_diagnostics = gks->info.collisions.write_diagnostics;

  // Empty methods.
  lbo->moms_func = gklbo_moms_disabled;
  lbo->rhs_func = gklbo_rhs_disabled;
  lbo->write_mom_func = gklbo_write_mom_disabled;

  if (lbo->collision_id == GKYL_LBO_COLLISIONS) {
    lbo->num_cross_collisions = gks->info.collisions.num_cross_collisions;
    
    int cdim = app->cdim, vdim = gks->info.vdim;
  
    // Allocate self-species collision frequency and sum of collision frequencies.
    lbo->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    lbo->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
    double nu_frac = gks->info.collisions.nu_frac ? gks->info.collisions.nu_frac : 1.0;
  
    if (gks->info.collisions.self_nu) {
      // Project user's self-species collision frequency.
      lbo->norm_nu_self = false;
  
      struct gkyl_array *self_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, gks->info.collisions.self_nu, gks->info.collisions.self_nu_ctx);
      gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu_ho);
      gkyl_proj_on_basis_release(proj);
      gkyl_array_copy(lbo->self_nu, self_nu_ho);
      gkyl_array_release(self_nu_ho);
  
      gkyl_array_scale(lbo->self_nu, nu_frac);
      gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);
  
      // Set pointers to functions chosen at runtime.
      lbo->self_nu_func = gklbo_self_nu_calc_constNu;
    }
    else {
      // Self-collision frequency computed in time.
      lbo->norm_nu_self = true;
  
      double eps0 = gks->info.collisions.eps0 ? gks->info.collisions.eps0 : GKYL_EPSILON0;
      double hbar = gks->info.collisions.hbar ? gks->info.collisions.hbar : GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = gks->info.collisions.eV ? gks->info.collisions.eV : GKYL_ELEMENTARY_CHARGE;
      double bmag_ref = gks->info.collisions.bmag_ref ? gks->info.collisions.bmag_ref : app->bmag_ref;
  
      // Compute a minimum representable temperature based on the smallest dv in the grid.
      double dv_min[vdim];
      gkyl_velocity_map_reduce_dv_range(gks->vel_map, GKYL_MIN, dv_min, gks->vel_map->local_vel);
  
      double tpar_min = (gks->info.mass/6.0)*pow(dv_min[0],2);
      double tperp_min = vdim>1 ? (bmag_ref/3.0)*dv_min[1] : tpar_min;
      lbo->vtsq_min = (tpar_min + 2.0*tperp_min)/(3.0*gks->info.mass);
  
      lbo->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->basis, app->poly_order+1,
        1.0, 1.0, 1.0, app->use_gpu);

      // We define nu_ss = nu_sr(r=s) = alpha_E/((delta_ss * (1+beta))*n_s), with delta_ss = 2, 
      // beta = 0. This gives a nu_ss that is arguably 2X smaller than it should be, but it's
      // cheaper and yields an electron isotropization rate that agrees better with the FPO's.
      lbo->norm_nu_fac_self = nu_frac * gkyl_calc_Morse_alpha_E_const(
        gks->info.collisions.den_ref, gks->info.collisions.den_ref, 
        gks->info.mass, gks->info.mass, gks->info.charge, gks->info.charge,
        gks->info.collisions.temp_ref, gks->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV);
  
      // Set pointers to functions chosen at runtime.
      lbo->self_nu_func = gklbo_self_nu_calc_normNu;
    }
  
    // Create moment calculator to get M0, M1, M2 for primitive moments.
    gk_species_moment_init(app, gks, &lbo->moms, GKYL_F_MOMENT_M0M1M2, false);
    lbo->nu_moms = mkarr(app->use_gpu, lbo->moms.marr->ncomp, lbo->moms.marr->size);
  
    // Edge of velocity space corrections to momentum and energy.
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(&gks->grid, 
      &app->basis, &gks->basis, gks->info.mass, gks->vel_map, app->use_gpu);
    
    // Primitive moment calculator.
    lbo->coll_pcalc = gkyl_prim_lbo_gyrokinetic_calc_new(&gks->grid, 
      &app->basis, &gks->basis, &app->local, app->use_gpu);
  
    // Allocate boundary corrections for primitive mom calculation.
    lbo->boundary_corrections = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
    lbo->nu_boundary_corrections = mkarr(app->use_gpu, lbo->boundary_corrections->ncomp, lbo->boundary_corrections->size);
  
    // Primitive moments.
    lbo->prim_moms = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
    lbo->nu_prim_moms = mkarr(app->use_gpu, lbo->prim_moms->ncomp, lbo->prim_moms->size);

    // M2 of this species.
    lbo->m2self = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
    // Host-side copy for I/O.
    if (lbo->write_diagnostics) {
      if (app->use_gpu) {
        lbo->nu_sum_host = mkarr(false, lbo->nu_sum->ncomp, lbo->nu_sum->size);
        lbo->nu_prim_moms_host = mkarr(false, lbo->nu_prim_moms->ncomp, lbo->nu_prim_moms->size);
      }
      else {
        lbo->nu_sum_host = lbo->nu_sum;
        lbo->nu_prim_moms_host = lbo->nu_prim_moms;
      }
    }
  
    // LBO updater.
    struct gkyl_dg_lbo_gyrokinetic_drag_auxfields drag_inp = { .nuSum = lbo->nu_sum, 
      .nuPrimMomsSum = lbo->nu_prim_moms, .m2self = lbo->m2self };
    struct gkyl_dg_lbo_gyrokinetic_diff_auxfields diff_inp = { .nuSum = lbo->nu_sum, 
      .nuPrimMomsSum = lbo->nu_prim_moms, .m2self = lbo->m2self };
    lbo->coll_slvr = gkyl_dg_updater_lbo_gyrokinetic_new(&gks->grid, 
      &app->basis, &gks->basis, &app->local, &drag_inp, &diff_inp, gks->info.mass, 
      app->gk_geom, gks->vel_map,  app->use_gpu);

    // Methods chosen at runtime.
    lbo->moms_func = gklbo_moms_enabled;
    if (! gks->info.collisions.not_in_dfdt) {;
      lbo->rhs_func = gklbo_rhs_enabled;
    }
    if (lbo->write_diagnostics)
      lbo->write_mom_func = gklbo_write_mom_enabled;
  }
}

void 
gk_species_lbo_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_lbo_collisions *lbo)
{
  // Empty methods.
  lbo->cross_nu_func = gklbo_cross_nu_calc_constNu;
  lbo->cross_moms_func = gklbo_cross_moms_disabled;

  if (lbo->collision_id == GKYL_LBO_COLLISIONS) {
    if (gks->lbo.num_cross_collisions) {
      lbo->betaGreenep1 = 1.0; // Greene's beta factor + 1.
      lbo->delta_sr = 2.0; // delta_sr free parameter.
        
      // Set pointers to species we cross-collide with.
      int my_idx_in_other[GKYL_MAX_SPECIES];
      for (int i=0; i<lbo->num_cross_collisions; ++i) {
        lbo->collide_with[i] = gk_find_species(app, gks->info.collisions.collide_with[i]);
        my_idx_in_other[i] = -1;
        for (int j=0; j<lbo->collide_with[i]->lbo.num_cross_collisions; ++j) {
          if (0 == strcmp(gks->info.name, lbo->collide_with[i]->info.collisions.collide_with[j])) {
            my_idx_in_other[i] = j;
            break;
          }
        }
      }

      // Morse's alpha_E.
      lbo->alpha_E = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      for (int i=0; i<lbo->num_cross_collisions; ++i) {
        // Cross-species collision frequency, nu_sr.
        lbo->cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
        lbo->other_m[i] = lbo->collide_with[i]->info.mass;
        lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
      }

      double nu_frac = gks->info.collisions.nu_frac ? gks->info.collisions.nu_frac : 1.0;

      // Compute the time-independent part of alpha_E.
      double alpha_E_norm[GKYL_MAX_SPECIES] = {0.0};
      for (int i=0; i<lbo->num_cross_collisions; ++i) {
        double eps0 = gks->info.collisions.eps0 ? gks->info.collisions.eps0: GKYL_EPSILON0;
        double hbar = gks->info.collisions.hbar ? gks->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
        double eV = gks->info.collisions.eV ? gks->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
        double bmag_ref = gks->info.collisions.bmag_ref ? gks->info.collisions.bmag_ref : app->bmag_ref;
        double mass_self = gks->info.mass, mass_other = lbo->collide_with[i]->info.mass;

        alpha_E_norm[i] = nu_frac * gkyl_calc_Morse_alpha_E_const(
          gks->info.collisions.den_ref, lbo->collide_with[i]->info.collisions.den_ref,
          mass_self, mass_other, gks->info.charge, lbo->collide_with[i]->info.charge,
          gks->info.collisions.temp_ref, lbo->collide_with[i]->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV);
      }

      if (gks->info.collisions.cross_nu[0]) {
        // Project user's cross-species collision frequency.
        lbo->norm_nu_cross = false;

        // Ensure the other species this collides with also provided self_nu and cross_nu.
        for (int i=0; i<lbo->num_cross_collisions; ++i) {
          assert(lbo->collide_with[i]->info.collisions.self_nu);
          assert(lbo->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]);
        }

        struct gkyl_array *cross_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
        for (int i=0; i<lbo->num_cross_collisions; ++i) {
          gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
            app->poly_order+1, 1, gks->info.collisions.cross_nu[i], gks->info.collisions.cross_nu_ctx);
          gkyl_proj_on_basis_advance(proj, 0.0, &app->local, cross_nu_ho);
          gkyl_proj_on_basis_release(proj);
          gkyl_array_copy(lbo->cross_nu[i], cross_nu_ho);

          gkyl_array_scale(lbo->cross_nu[i], nu_frac);
          gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

          // Compute alpha_E using reference parameters.
          assert(gks->info.collisions.den_ref);
          assert(lbo->collide_with[i]->info.collisions.den_ref);
          assert(gks->info.collisions.temp_ref);
          assert(lbo->collide_with[i]->info.collisions.temp_ref);
          double mass_self = gks->info.mass, mass_other = lbo->collide_with[i]->info.mass;
          double den_s = gks->info.collisions.den_ref;
          double den_r = lbo->collide_with[i]->info.collisions.den_ref;
          double vtsq_s = gks->info.collisions.temp_ref/mass_self;
          double vtsq_r = lbo->collide_with[i]->info.collisions.temp_ref/mass_other;

          lbo->alpha_E_fac[i] = ( alpha_E_norm[i] * den_s * den_r / pow(sqrt(vtsq_s+vtsq_r),3.0) ) * pow(sqrt(2.0),app->cdim);
        }
        gkyl_array_release(cross_nu_ho);

        // Set pointers to functions chosen at runtime.
        lbo->cross_nu_func = gklbo_cross_nu_calc_constNu;
        lbo->alpha_E_func = gklbo_alpha_E_constNu;
      }
      else {
        // Cross-collision frequency computed in time.
        lbo->norm_nu_cross = true;

        // Ensure the other species this collides with didn't provide self_nu nor cross_nu.
        for (int i=0; i<lbo->num_cross_collisions; ++i) {
          assert(!(lbo->collide_with[i]->info.collisions.self_nu));
          assert(!(lbo->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]));
        }

        for (int i=0; i<lbo->num_cross_collisions; ++i) {
          double mass_self = gks->info.mass, mass_other = lbo->collide_with[i]->info.mass;

          lbo->norm_nu_fac_cross[i] = alpha_E_norm[i]
            * (mass_self+mass_other)/(lbo->delta_sr*lbo->betaGreenep1*mass_self);

          lbo->alpha_E_fac[i] = (lbo->delta_sr*lbo->betaGreenep1*mass_self)/(mass_self+mass_other);
        }
    
        // Set pointers to functions chosen at runtime.
        lbo->cross_nu_func = gklbo_cross_nu_calc_normNu;
        lbo->alpha_E_func = gklbo_alpha_E_normNu;
      }

      // Cross-primitive moment calculator.
      lbo->cross_calc = gkyl_prim_lbo_gyrokinetic_cross_calc_new(&gks->grid, 
        &app->basis, &gks->basis, &app->local, app->use_gpu);

      // Methods chosen at runtime.
      lbo->cross_moms_func = gklbo_cross_moms_enabled;
    }
  }
}

void
gk_species_lbo_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  lbo->moms_func(app, species, lbo, fin);
}

void
gk_species_lbo_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  lbo->cross_moms_func(app, species, lbo);
}

void
gk_species_lbo_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  lbo->rhs_func(app, gks, lbo, fin, rhs);
}

void
gk_species_lbo_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->lbo.write_mom_func(app, gks, tm, frame);
}

void 
gk_species_lbo_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lbo_collisions *lbo)
{
  if (lbo->collision_id == GKYL_LBO_COLLISIONS) {
    if (lbo->num_cross_collisions) {
      gkyl_prim_lbo_cross_calc_release(lbo->cross_calc);

      for (int i=0; i<lbo->num_cross_collisions; ++i)
        gkyl_array_release(lbo->cross_nu[i]);

      gkyl_array_release(lbo->alpha_E);
    }

    gkyl_dg_updater_lbo_gyrokinetic_release(lbo->coll_slvr);

    if (lbo->write_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(lbo->nu_sum_host);
        gkyl_array_release(lbo->nu_prim_moms_host);    
      }
    }

    gkyl_array_release(lbo->m2self);

    gkyl_array_release(lbo->nu_prim_moms);
    gkyl_array_release(lbo->prim_moms);

    gkyl_array_release(lbo->nu_boundary_corrections);
    gkyl_array_release(lbo->boundary_corrections);

    gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
    gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

    gkyl_array_release(lbo->nu_moms);
    gk_species_moment_release(app, &lbo->moms);

    if (lbo->norm_nu_self)
      gkyl_spitzer_coll_freq_release(lbo->spitzer_calc);

    gkyl_array_release(lbo->nu_sum);
    gkyl_array_release(lbo->self_nu);
  }
}
