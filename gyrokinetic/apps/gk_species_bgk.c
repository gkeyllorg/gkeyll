#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

static void
gkbgk_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Empty method.
}

static void
gkbgk_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_par, T/m).
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, species->lte.moms.marr,
    0, species->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);
  
  // Calculate nu_ss.
  bgk->self_nu_func(app, species, bgk);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
gkbgk_self_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk)
{
  // Empty method.
}

static void
gkbgk_self_nu_calc_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk)
{
  // Calculate nu_ss(x,t).
  gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, species->lte.moms.marr, bgk->vtsq_min,
    species->lte.moms.marr, bgk->vtsq_min, bgk->norm_nu_fac_self, bgk->self_nu);

  gkyl_array_set(bgk->nu_sum, 1.0, bgk->self_nu);
}

static void
gkbgk_cross_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_bgk_collisions *bgk, int coll_idx)
{
  // Empty method.
}

static void
gkbgk_cross_nu_calc_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_bgk_collisions *bgk, int coll_idx)
{
  // Calculate nu_sr(x,t).
  gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, s->lte.moms.marr, bgk->vtsq_min,
    bgk->collide_with[coll_idx]->lte.moms.marr, bgk->collide_with[coll_idx]->bgk.vtsq_min,
    bgk->norm_nu_fac_cross[coll_idx], bgk->cross_nu[coll_idx]);

  gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[coll_idx]);
}

static void
gkbgk_alpha_E_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_bgk_collisions *bgk, int coll_idx)
{
  gkyl_array_clear(bgk->alpha_E, 0.0);
  gkyl_array_shiftc(bgk->alpha_E, bgk->alpha_E_fac[coll_idx], 0);
}

static void
gkbgk_alpha_E_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_bgk_collisions *bgk, int coll_idx)
{
  gkyl_dg_mul_op_range(&app->basis, 0, bgk->alpha_E, 0, bgk->cross_nu[coll_idx], 0, s->lte.moms.marr, &app->local);
  gkyl_array_scale_range(bgk->alpha_E, bgk->alpha_E_fac[coll_idx], &app->local);
}

static void
gkbgk_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_bgk_collisions *bgk, int coll_idx)
{
  // Empty method.
}

static void
gkbgk_cross_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_bgk_collisions *bgk, int coll_idx)
{
  // Compute primitive moments for cross-species collisions.
  struct timespec wst = gkyl_wall_clock();
  
  // Compute the cross-species collision frequency.
  bgk->cross_nu_func(app, gks, bgk, coll_idx);

  // Compute alpha_E.
  bgk->alpha_E_func(app, gks, bgk, coll_idx);

  // Compute cross primitive moments.
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance(bgk->cross_calc, &app->local, bgk->delta_sr, bgk->betaGreenep1,
    gks->info.mass, gks->lte.moms.marr, bgk->other_m[coll_idx], bgk->other_prim_moms[coll_idx],
    bgk->cross_prim_moms);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
gkbgk_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  // Empty method.
}

static void
gkbgk_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // Compute the self-collisions Maxwellian.
  struct gk_lte *lte = &gks->lte;
  gk_species_lte_from_moms(app, gks, lte, lte->moms.marr);

  // Multiply the Maxwellian by self-species collision frequency.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, bgk->nu_fmax, 
    bgk->self_nu, lte->f_lte, &app->local, &gks->local);

  // Cross-collisions nu*fmax.
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    // Compute cross-primitive moments
    bgk->cross_moms_func(app, gks, bgk, i);

    // Compute the cross-collisions Maxwellian.
    gk_species_lte_from_moms(app, gks, lte, bgk->cross_prim_moms);

    // Compute and accumulate nu*fmax.
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, lte->f_lte, 
      bgk->cross_nu[i], lte->f_lte, &app->local, &gks->local);
    gkyl_array_accumulate(bgk->nu_fmax, 1.0, lte->f_lte);
  }

  // Multiply the Maxwellian by the configuration-space Jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, bgk->nu_fmax,
    app->gk_geom->geo_int.jacobgeo, bgk->nu_fmax, &app->local, &gks->local);

  // Apply BGK collisions.
  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &gks->local, 
    bgk->nu_sum, bgk->nu_fmax, fin, bgk->implicit_step, dt, rhs, gks->cflrate);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
gkbgk_write_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // Empty method.
}

static void
gkbgk_write_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wtm = gkyl_wall_clock();
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_conf_len, gks->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_conf_len, gks->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Sum of collision frequencies." }
  };
  int io_meta_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_conf, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write out nu_sum.
  const char *fmt = "%s-%s_bgk_nu_sum_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  
  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gks->bgk.nu_sum_host, gks->bgk.nu_sum);
  }
  
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->bgk.nu_sum_host, fileNm);
  app->stat.n_diag_io += 2;

  gkyl_msgpack_data_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void 
gk_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_bgk_collisions *bgk)
{
  bgk->collision_id = gks->info.collisions.collision_id;
  bgk->write_diagnostics = gks->info.collisions.write_diagnostics;

  // Empty methods.
  bgk->moms_func = gkbgk_moms_disabled;
  bgk->rhs_func = gkbgk_rhs_disabled;
  bgk->moms_func_implicit = gkbgk_moms_disabled;
  bgk->rhs_func_implicit = gkbgk_rhs_disabled;
  bgk->write_mom_func = gkbgk_write_mom_disabled;

  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    bgk->implicit_step = gks->info.collisions.is_implicit;
    bgk->num_cross_collisions = gks->info.collisions.num_cross_collisions;
    
    int cdim = app->cdim, vdim = gks->info.vdim;
  
    // Allocate self-species collision frequency and sum of collision frequencies.
    bgk->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    bgk->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
    double nu_frac = gks->info.collisions.nu_frac ? gks->info.collisions.nu_frac : 1.0;
  
    if (gks->info.collisions.self_nu) {
      // Project user's self-species collision frequency.
      bgk->norm_nu_self = false;
  
      struct gkyl_array *self_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, gks->info.collisions.self_nu, gks->info.collisions.self_nu_ctx);
      gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu_ho);
      gkyl_proj_on_basis_release(proj);
      gkyl_array_copy(bgk->self_nu, self_nu_ho);
      gkyl_array_release(self_nu_ho);
  
      gkyl_array_scale(bgk->self_nu, nu_frac);
      gkyl_array_set(bgk->nu_sum, 1.0, bgk->self_nu);
  
      // Set pointers to functions chosen at runtime.
      bgk->self_nu_func = gkbgk_self_nu_calc_constNu;
    }
    else {
      // Self-collision frequency computed in time.
      bgk->norm_nu_self = true;
  
      double eps0 = gks->info.collisions.eps0 ? gks->info.collisions.eps0 : GKYL_EPSILON0;
      double hbar = gks->info.collisions.hbar ? gks->info.collisions.hbar : GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = gks->info.collisions.eV ? gks->info.collisions.eV : GKYL_ELEMENTARY_CHARGE;
      double bmag_ref = gks->info.collisions.bmag_ref ? gks->info.collisions.bmag_ref : app->bmag_ref;
  
      // Compute a minimum representable temperature based on the smallest dv in the grid.
      double dv_min[vdim];
      gkyl_velocity_map_reduce_dv_range(gks->vel_map, GKYL_MIN, dv_min, gks->vel_map->local_vel);
  
      double tpar_min = (gks->info.mass/6.0)*pow(dv_min[0],2);
      double tperp_min = vdim>1 ? (bmag_ref/3.0)*dv_min[1] : tpar_min;
      bgk->vtsq_min = (tpar_min + 2.0*tperp_min)/(3.0*gks->info.mass);
  
      bgk->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->basis, app->poly_order+1,
        1.0, 1.0, 1.0, app->use_gpu);

      // We define nu_ss = nu_sr(r=s) = alpha_E/((delta_ss * (1+beta))*n_s), with delta_ss = 2,
      // beta = 0. This gives a nu_ss that is arguably 2X smaller than it should be, but it's
      // cheaper and yields an electron isotropization rate that agrees better with the FPO's.
      bgk->norm_nu_fac_self = nu_frac * gkyl_calc_Morse_alpha_E_const(
        gks->info.collisions.den_ref, gks->info.collisions.den_ref, 
        gks->info.mass, gks->info.mass, gks->info.charge, gks->info.charge,
        gks->info.collisions.temp_ref, gks->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV);
  
      // Set pointers to functions chosen at runtime.
      bgk->self_nu_func = gkbgk_self_nu_calc_normNu;
    }
  
    // Host-side copy for I/O.
    if (bgk->write_diagnostics) {
      if (app->use_gpu) {
        bgk->nu_sum_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      }
      else {
        bgk->nu_sum_host = bgk->nu_sum;
      }
    }
  
    // Collision frequency times Maxwellian.
    bgk->nu_fmax = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);

    // BGK updater.
    bgk->up_bgk = gkyl_bgk_collisions_new(&app->basis, &gks->basis, app->use_gpu);

    // Methods chosen at runtime.
    if (bgk->implicit_step) {
      bgk->moms_func = gkbgk_moms_disabled;
      bgk->rhs_func = gkbgk_rhs_disabled;
      bgk->moms_func_implicit = gkbgk_moms_enabled;
      if (!gks->info.collisions.not_in_dfdt) {
        bgk->rhs_func_implicit = gkbgk_rhs_enabled;
      }
    }
    else {
      bgk->moms_func = gkbgk_moms_enabled;
      if (!gks->info.collisions.not_in_dfdt) {
        bgk->rhs_func = gkbgk_rhs_enabled;
      }
      bgk->moms_func_implicit = gkbgk_moms_disabled;
      bgk->rhs_func_implicit = gkbgk_rhs_disabled;
    }
    if (bgk->write_diagnostics)
      bgk->write_mom_func = gkbgk_write_mom_enabled;
  }
}

void 
gk_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_bgk_collisions *bgk)
{  
  // Empty methods.
  bgk->cross_nu_func = gkbgk_cross_nu_calc_constNu;
  bgk->cross_moms_func = gkbgk_cross_moms_disabled;

  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    if (gks->bgk.num_cross_collisions) {
      bgk->betaGreenep1 = 1.0; // Greene's beta factor + 1.
      bgk->delta_sr = 2.0; // delta_sr free parameter.
        
      // Set pointers to species we cross-collide with.
      int my_idx_in_other[GKYL_MAX_SPECIES];
      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        bgk->collide_with[i] = gk_find_species(app, gks->info.collisions.collide_with[i]);
        my_idx_in_other[i] = -1;
        for (int j=0; j<bgk->collide_with[i]->bgk.num_cross_collisions; ++j) {
          if (0 == strcmp(gks->info.name, bgk->collide_with[i]->info.collisions.collide_with[j])) {
            my_idx_in_other[i] = j;
            break;
          }
        }
      }

      // Morse's alpha_E.
      bgk->alpha_E = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      // Cross primitive moments (n_sr, u_{parallel sr}, v_{t,sr}^2).
      bgk->cross_prim_moms = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        // Cross-species collision frequency, nu_sr.
        bgk->cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
        bgk->other_m[i] = bgk->collide_with[i]->info.mass;
        bgk->other_prim_moms[i] = bgk->collide_with[i]->lte.moms.marr;
      }

      double nu_frac = gks->info.collisions.nu_frac ? gks->info.collisions.nu_frac : 1.0;

      // Compute the time-independent part of alpha_E.
      double alpha_E_norm[GKYL_MAX_SPECIES] = {0.0};
      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        double eps0      = gks->info.collisions.eps0 ? gks->info.collisions.eps0: GKYL_EPSILON0;
        double hbar      = gks->info.collisions.hbar ? gks->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
        double eV        = gks->info.collisions.eV ? gks->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
        double bmag_ref  = gks->info.collisions.bmag_ref ? gks->info.collisions.bmag_ref : app->bmag_ref;
        double mass_self = gks->info.mass, mass_other = bgk->collide_with[i]->info.mass;

        alpha_E_norm[i] = nu_frac * gkyl_calc_Morse_alpha_E_const(
          gks->info.collisions.den_ref, bgk->collide_with[i]->info.collisions.den_ref,
          mass_self, mass_other, gks->info.charge, bgk->collide_with[i]->info.charge,
          gks->info.collisions.temp_ref, bgk->collide_with[i]->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV);
      }

      if (gks->info.collisions.cross_nu[0]) {
        // Project user's cross-species collision frequency.
        bgk->norm_nu_cross = false;

        // Ensure the other species this collides with also provided self_nu and cross_nu.
        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          assert(bgk->collide_with[i]->info.collisions.self_nu);
          assert(bgk->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]);
        }

        struct gkyl_array *cross_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
            app->poly_order+1, 1, gks->info.collisions.cross_nu[i], gks->info.collisions.cross_nu_ctx);
          gkyl_proj_on_basis_advance(proj, 0.0, &app->local, cross_nu_ho);
          gkyl_proj_on_basis_release(proj);
          gkyl_array_copy(bgk->cross_nu[i], cross_nu_ho);

          gkyl_array_scale(bgk->cross_nu[i], nu_frac);
          gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);

          // Compute alpha_E using reference parameters.
          assert(gks->info.collisions.den_ref);
          assert(bgk->collide_with[i]->info.collisions.den_ref);
          assert(gks->info.collisions.temp_ref);
          assert(bgk->collide_with[i]->info.collisions.temp_ref);
          double mass_self = gks->info.mass, mass_other = bgk->collide_with[i]->info.mass;
          double den_s = gks->info.collisions.den_ref;
          double den_r = bgk->collide_with[i]->info.collisions.den_ref;
          double vtsq_s = gks->info.collisions.temp_ref/mass_self;
          double vtsq_r = bgk->collide_with[i]->info.collisions.temp_ref/mass_other;

          bgk->alpha_E_fac[i] = ( alpha_E_norm[i] * den_s * den_r / pow(sqrt(vtsq_s+vtsq_r),3.0) ) * pow(sqrt(2.0),app->cdim);
        }
        gkyl_array_release(cross_nu_ho);

        // Set pointers to functions chosen at runtime.
        bgk->cross_nu_func = gkbgk_cross_nu_calc_constNu;
        bgk->alpha_E_func = gkbgk_alpha_E_constNu;
      }
      else {
        // Cross-collision frequency computed in time.
        bgk->norm_nu_cross = true;

        // Ensure the other species this collides with didn't provide self_nu nor cross_nu.
        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          assert(!(bgk->collide_with[i]->info.collisions.self_nu));
          assert(!(bgk->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]));
        }

        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          double mass_self = gks->info.mass, mass_other = bgk->collide_with[i]->info.mass;

          bgk->norm_nu_fac_cross[i] = alpha_E_norm[i]
            * (mass_self+mass_other)/(bgk->delta_sr*bgk->betaGreenep1*mass_self);

          bgk->alpha_E_fac[i] = (bgk->delta_sr*bgk->betaGreenep1*mass_self)/(mass_self+mass_other);
        }

        // Set pointers to functions chosen at runtime.
        bgk->cross_nu_func = gkbgk_cross_nu_calc_normNu;
        bgk->alpha_E_func = gkbgk_alpha_E_normNu;
      }

      // Cross-primitive moment calculator.
      bgk->cross_calc = gkyl_gyrokinetic_cross_prim_moms_bgk_new(&gks->basis, &app->basis, app->use_gpu);

      // Methods chosen at runtime.
      bgk->cross_moms_func = gkbgk_cross_moms_enabled;
    }
  }
}

void
gk_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->moms_func(app, species, bgk, fin);
}

void
gk_species_bgk_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->moms_func_implicit(app, species, bgk, fin);
}

void
gk_species_bgk_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Compute this in bgk_rhs
}

void
gk_species_bgk_cross_moms_implicit(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Compute this in bgk_rhs
}

void
gk_species_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  bgk->rhs_func(app, gks, bgk, fin, 0.0, rhs);
}

void
gk_species_bgk_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  bgk->rhs_func_implicit(app, gks, bgk, fin, dt, rhs);
}

void
gk_species_bgk_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->bgk.write_mom_func(app, gks, tm, frame);
}

void 
gk_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk)
{
  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    if (bgk->num_cross_collisions) {
      gkyl_gyrokinetic_cross_prim_moms_bgk_release(bgk->cross_calc);

      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        gkyl_array_release(bgk->cross_nu[i]);
      }
      gkyl_array_release(bgk->alpha_E);
      gkyl_array_release(bgk->cross_prim_moms);
    }

    gkyl_array_release(bgk->nu_fmax);
    gkyl_bgk_collisions_release(bgk->up_bgk);

    if (bgk->write_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(bgk->nu_sum_host);
      }
    }

    if (bgk->norm_nu_self)
      gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);

    gkyl_array_release(bgk->nu_sum);
    gkyl_array_release(bgk->self_nu);
  }
}
