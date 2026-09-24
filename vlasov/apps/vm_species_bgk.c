#include <assert.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_const.h>

static void
vmbgk_moms_disabled(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Empty method.
}

static void
vmbgk_moms_enabled(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_i, T/m).
  vm_species_moment_calc(&vms->lte.moms, vms->local, app->local, fin);
  
  // Calculate nu_ss.
  bgk->self_nu_func(app, vms, bgk);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
vmbgk_self_nu_calc_constNu(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk)
{
  // Zero out collision frequency in cells with n, T < 0.
  gkyl_bgk_collisions_correct_nu(bgk->up_bgk, &app->local, vms->lte.moms.marr, 
    bgk->ref_self_nu, bgk->self_nu); 

  gkyl_array_set(bgk->nu_sum, 1.0, bgk->self_nu);    
}

static void
vmbgk_self_nu_calc_normNu(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk)
{
  // Calculate nu_ss(x,t). Automatically zeros out nu_ss if n, T < 0. 
  gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, vms->lte.moms.marr, bgk->vtsq_min,
    vms->lte.moms.marr, bgk->vtsq_min, bgk->norm_nu_fac_self, bgk->self_nu);

  gkyl_array_set(bgk->nu_sum, 1.0, bgk->self_nu);
}

static void
vmbgk_cross_nu_calc_constNu(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, int coll_idx)
{
  // Zero out cross-collision frequency in cells with n, T < 0.
  gkyl_bgk_collisions_correct_nu(bgk->up_bgk, &app->local, bgk->other_prim_moms[coll_idx], 
    bgk->ref_cross_nu[coll_idx], bgk->cross_nu[coll_idx]);   

  gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[coll_idx]);    
}

static void
vmbgk_cross_nu_calc_normNu(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, int coll_idx)
{
  // Calculate nu_sr(x,t). Automatically zeros out nu_sr if n, T < 0. 
  gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, vms->lte.moms.marr, bgk->vtsq_min,
    bgk->other_prim_moms[coll_idx], bgk->collide_with[coll_idx]->bgk.vtsq_min,
    bgk->norm_nu_fac_cross[coll_idx], bgk->cross_nu[coll_idx]);

  gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[coll_idx]);
}

static void
vmbgk_alpha_E_constNu(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_bgk_collisions *bgk, int coll_idx)
{
  gkyl_array_clear(bgk->alpha_E, 0.0);
  gkyl_array_shiftc(bgk->alpha_E, bgk->alpha_E_fac[coll_idx], 0);
}

static void
vmbgk_alpha_E_normNu(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_bgk_collisions *bgk, int coll_idx)
{
  gkyl_dg_mul_op_range(&app->basis, 0, bgk->alpha_E, 0, bgk->cross_nu[coll_idx], 0, s->lte.moms.marr, &app->local);
  gkyl_array_scale_range(bgk->alpha_E, bgk->alpha_E_fac[coll_idx], &app->local);
}

static void
vmbgk_cross_moms_disabled(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, int coll_idx)
{
  // Empty method.
}

static void
vmbgk_cross_moms_enabled(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, int coll_idx)
{
  // Compute primitive moments for cross-species collisions.
  struct timespec wst = gkyl_wall_clock();
  
  // Compute the cross-species collision frequency.
  bgk->cross_nu_func(app, vms, bgk, coll_idx);

  // Compute alpha_E.
  bgk->alpha_E_func(app, vms, bgk, coll_idx);

  // Compute cross primitive moments.
  gkyl_vlasov_cross_prim_moms_bgk_advance(bgk->cross_calc, &app->local, bgk->delta_sr, bgk->betaGreenep1,
    vms->mass, vms->lte.moms.marr, bgk->other_m[coll_idx], bgk->other_prim_moms[coll_idx],
    bgk->cross_prim_moms);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

static void
vmbgk_rhs_disabled(gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  // Empty method.
}

static void
vmbgk_rhs_enabled(gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // Compute the self-collisions Maxwellian.
  struct vm_lte *lte = &vms->lte;
  if (bgk->fixed_temp_relax) { 
    // Set the temperature in the moment array to the pre-computed fixed value
    gkyl_array_set_offset_range(lte->moms.marr, 1.0, bgk->fixed_temp, 
      (app->vdim+1)*app->basis.num_basis, &app->local);
  }
  vm_species_lte_from_moms(app, vms, lte, lte->moms.marr);

  // Multiply the Maxwellian by self-species collision frequency.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &vms->basis, bgk->nu_f_lte, 
    bgk->self_nu, lte->f_lte, &app->local, &vms->local);

  // Cross-collisions nu*fmax.
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    // Compute cross-primitive moments
    bgk->cross_moms_func(app, vms, bgk, i);

    // Compute the cross-collisions Maxwellian.
    vm_species_lte_from_moms(app, vms, lte, bgk->cross_prim_moms);

    // Compute and accumulate nu*fmax.
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &vms->basis, lte->f_lte, 
      bgk->cross_nu[i], lte->f_lte, &app->local, &vms->local);
    gkyl_array_accumulate(bgk->nu_f_lte, 1.0, lte->f_lte);
  }

  // Apply BGK collisions.
  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &vms->local, 
    bgk->nu_sum, bgk->nu_f_lte, fin, bgk->implicit_step, dt, rhs, vms->cflrate);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
vmbgk_fixed_temp_disabled(gkyl_vlasov_app* app, const struct vm_species *vms, 
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  // Empty method.
}

static void
vmbgk_fixed_temp_enabled(gkyl_vlasov_app* app, const struct vm_species *vms, 
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  vm_species_moment_calc(&vms->lte.moms, vms->local, app->local, fin);

  // Set the temperature to the fixed value 
  gkyl_array_set_offset_range(bgk->fixed_temp, 1.0, vms->lte.moms.marr, 
    (app->vdim+1)*app->basis.num_basis, &app->local);
}

static void
vmbgk_write_mom_disabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // Empty method.
}

static void
vmbgk_write_mom_enabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wtm = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    } 
  );

  // Write out nu_sum.
  const char *fmt = "%s-%s_bgk_nu_sum_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);
  
  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(vms->bgk.nu_sum_host, vms->bgk.nu_sum);
  }
  
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, vms->bgk.nu_sum_host, fileNm);
  app->stat.n_diag_io += 1;

  vlasov_array_meta_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void 
vm_species_bgk_init(struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_bgk_collisions *bgk)
{
  bgk->collision_id = vms->info.collisions.collision_id;
  bgk->write_coll_diagnostics = vms->info.collisions.write_coll_diagnostics;

  // Empty methods.
  bgk->moms_func = vmbgk_moms_disabled;
  bgk->moms_func_implicit = vmbgk_moms_disabled;
  bgk->rhs_func = vmbgk_rhs_disabled;
  bgk->rhs_func_implicit = vmbgk_rhs_disabled;
  bgk->fixed_temp_calc_func = vmbgk_fixed_temp_disabled; 
  bgk->write_mom_func = vmbgk_write_mom_disabled;

  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    bgk->implicit_step = vms->info.collisions.is_implicit;
    bgk->num_cross_collisions = vms->info.collisions.num_cross_collisions;
    
    int cdim = app->cdim, vdim = app->vdim;
  
    // Allocate self-species collision frequency and sum of collision frequencies.
    bgk->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    bgk->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
    double nu_frac = vms->info.collisions.nu_frac ? vms->info.collisions.nu_frac : 1.0;

    // The Spitzer collision-frequency updater (and the vtsq_min floor it needs)
    // is required whenever any collision frequency is computed in time: a computed
    // self_nu, or a computed cross_nu (cross collisions for which no explicit
    // cross_nu was provided). It is independent of the self_nu mode, so a species
    // may use an explicit self_nu together with a computed cross_nu.
    bool computed_cross = vms->info.collisions.num_cross_collisions > 0
      && !vms->info.collisions.cross_nu[0];
    if (!vms->info.collisions.self_nu || computed_cross) {
      // Compute a minimum representable temperature based on the smallest dv in the grid.
      double vtsq_min = 0.0;
      for (int d=0; d<vdim; ++d) {
        vtsq_min += (1.0/6.0)*pow(vms->grid.dx[cdim+d],2);
      }
      bgk->vtsq_min = vtsq_min/vdim;

      bgk->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->basis, app->poly_order+1,
        1.0, 1.0, 1.0, app->use_gpu);
    }

    if (vms->info.collisions.self_nu) {
      // Project user's self-species collision frequency.
      bgk->norm_nu_self = false;
  
      struct gkyl_array *self_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, vms->info.collisions.self_nu, vms->info.collisions.self_nu_ctx);
      gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu_ho);
      gkyl_proj_on_basis_release(proj);
      gkyl_array_copy(bgk->self_nu, self_nu_ho);
      gkyl_array_release(self_nu_ho);
  
      gkyl_array_scale(bgk->self_nu, nu_frac);
      gkyl_array_set(bgk->nu_sum, 1.0, bgk->self_nu);

      // Allocate additional array for checking n, T < 0 before applying collisions. 
      // Stores initial projected collision frequency so we do not need to reproject every time step. 
      bgk->ref_self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      gkyl_array_set(bgk->ref_self_nu, 1.0, bgk->self_nu);
  
      // Set pointers to functions chosen at runtime.
      bgk->self_nu_func = vmbgk_self_nu_calc_constNu;
    }
    else {
      // Self-collision frequency computed in time.
      bgk->norm_nu_self = true;
  
      double eps0 = vms->info.collisions.eps0 ? vms->info.collisions.eps0 : GKYL_EPSILON0;
      double hbar = vms->info.collisions.hbar ? vms->info.collisions.hbar : GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = vms->info.collisions.eV ? vms->info.collisions.eV : GKYL_ELEMENTARY_CHARGE;
      // Vlasov does not use reference magnetic field for cyclotron frequency contribution to log(Lambda)
      double bmag_ref = 0.0;
      // vtsq_min and spitzer_calc are set up above (needed for any computed frequency).

      // We define nu_ss = nu_sr(r=s) = alpha_E/((delta_ss * (1+beta))*n_s), with delta_ss = 2,
      // beta = 0. This gives a nu_ss that is arguably 2X smaller than it should be, but it's
      // cheaper and yields an electron isotropization rate that agrees better with the FPO's.
      bgk->norm_nu_fac_self = nu_frac * gkyl_calc_Morse_alpha_E_const(
        vms->info.collisions.den_ref, vms->info.collisions.den_ref, 
        vms->mass, vms->mass, vms->charge, vms->charge,
        vms->info.collisions.temp_ref, vms->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV);
  
      // Set pointers to functions chosen at runtime.
      bgk->self_nu_func = vmbgk_self_nu_calc_normNu;
    }
  
    // Host-side copy for I/O.
    if (bgk->write_coll_diagnostics) {
      if (app->use_gpu) {
        bgk->nu_sum_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      }
      else {
        bgk->nu_sum_host = bgk->nu_sum;
      }
    }
  
    // Collision frequency times Maxwellian.
    bgk->nu_f_lte = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);

    // BGK updater.
    bgk->up_bgk = gkyl_bgk_collisions_new(&app->basis, &vms->basis, app->use_gpu);

    // Is the temperature being relaxed to fixed in time?
    bgk->fixed_temp_relax = vms->info.collisions.fixed_temp_relax;
    if (bgk->fixed_temp_relax) {
      bgk->fixed_temp = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      bgk->fixed_temp_calc_func = vmbgk_fixed_temp_enabled; 
    }  

    // Methods chosen at runtime.
    if (bgk->implicit_step) {
      bgk->moms_func = vmbgk_moms_disabled;
      bgk->rhs_func = vmbgk_rhs_disabled;
      bgk->moms_func_implicit = vmbgk_moms_enabled;
      bgk->rhs_func_implicit = vmbgk_rhs_enabled;
    }
    else {
      bgk->moms_func = vmbgk_moms_enabled;
      bgk->rhs_func = vmbgk_rhs_enabled;
      bgk->moms_func_implicit = vmbgk_moms_disabled;
      bgk->rhs_func_implicit = vmbgk_rhs_disabled;
    }
    if (bgk->write_coll_diagnostics) {
      bgk->write_mom_func = vmbgk_write_mom_enabled;
    }
  }
}

void 
vm_species_bgk_cross_init(struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_bgk_collisions *bgk)
{  
  // Empty methods.
  bgk->cross_nu_func = vmbgk_cross_nu_calc_constNu;
  bgk->cross_moms_func = vmbgk_cross_moms_disabled;

  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    if (vms->bgk.num_cross_collisions) {
      bgk->betaGreenep1 = 1.0; // Greene's beta factor + 1.
      bgk->delta_sr = 2.0; // delta_sr free parameter.
        
      // Set pointers to species we cross-collide with.
      int my_idx_in_other[GKYL_MAX_SPECIES];
      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        bgk->collide_with[i] = vm_find_species(app, vms->info.collisions.collide_with[i]);
        // collide_with must name an existing *kinetic* species (a typo, or a
        // fluid species, returns NULL and would segfault below without a message).
        assert(bgk->collide_with[i]);
        my_idx_in_other[i] = -1;
        for (int j=0; j<bgk->collide_with[i]->bgk.num_cross_collisions; ++j) {
          if (0 == strcmp(vms->name, bgk->collide_with[i]->info.collisions.collide_with[j])) {
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
      bgk->alpha_E = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      // Cross primitive moments (n_sr, u_{parallel sr}, v_{t,sr}^2).
      bgk->cross_prim_moms = mkarr(app->use_gpu, (app->vdim+2)*app->basis.num_basis, app->local_ext.volume);
      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        // Cross-species collision frequency, nu_sr.
        bgk->cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
        bgk->other_m[i] = bgk->collide_with[i]->mass;
        bgk->other_prim_moms[i] = bgk->collide_with[i]->lte.moms.marr;
      }

      double nu_frac = vms->info.collisions.nu_frac ? vms->info.collisions.nu_frac : 1.0;

      // Compute the time-independent part of alpha_E.
      double alpha_E_norm[GKYL_MAX_SPECIES] = {0.0};
      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        double eps0      = vms->info.collisions.eps0 ? vms->info.collisions.eps0: GKYL_EPSILON0;
        double hbar      = vms->info.collisions.hbar ? vms->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
        double eV        = vms->info.collisions.eV ? vms->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
        // Vlasov does not use reference magnetic field for cyclotron frequency contribution to log(Lambda)
        double bmag_ref = 0.0;
        double mass_self = vms->mass, mass_other = bgk->collide_with[i]->mass;

        alpha_E_norm[i] = nu_frac * gkyl_calc_Morse_alpha_E_const(
          vms->info.collisions.den_ref, bgk->collide_with[i]->info.collisions.den_ref,
          mass_self, mass_other, vms->charge, bgk->collide_with[i]->charge,
          vms->info.collisions.temp_ref, bgk->collide_with[i]->info.collisions.temp_ref, bmag_ref, eps0, hbar, eV);
      }

      if (vms->info.collisions.cross_nu[0]) {
        // Project user's cross-species collision frequency.
        bgk->norm_nu_cross = false;

        // Cross-collision frequency must be specified symmetrically: if this
        // species provides an explicit cross_nu, its partner must too. The self
        // collision frequency mode is independent (a species may use explicit
        // self_nu with explicit or computed cross_nu), so it is not checked here.
        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          assert(bgk->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]);
        }

        struct gkyl_array *cross_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
            app->poly_order+1, 1, vms->info.collisions.cross_nu[i], vms->info.collisions.cross_nu_ctx[i]);
          gkyl_proj_on_basis_advance(proj, 0.0, &app->local, cross_nu_ho);
          gkyl_proj_on_basis_release(proj);
          gkyl_array_copy(bgk->cross_nu[i], cross_nu_ho);

          gkyl_array_scale(bgk->cross_nu[i], nu_frac);
          gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);

          // Allocate additional array for checking n, T < 0 before applying cross-collisions. 
          // Stores initial projected cross-collision frequency so we do not need to reproject every time step. 
          bgk->ref_cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
          gkyl_array_set(bgk->ref_cross_nu[i], 1.0, bgk->cross_nu[i]);

          // Compute alpha_E using reference parameters.
          assert(vms->info.collisions.den_ref);
          assert(bgk->collide_with[i]->info.collisions.den_ref);
          assert(vms->info.collisions.temp_ref);
          assert(bgk->collide_with[i]->info.collisions.temp_ref);
          double mass_self = vms->mass, mass_other = bgk->collide_with[i]->mass;
          double den_s = vms->info.collisions.den_ref;
          double den_r = bgk->collide_with[i]->info.collisions.den_ref;
          double vtsq_s = vms->info.collisions.temp_ref/mass_self;
          double vtsq_r = bgk->collide_with[i]->info.collisions.temp_ref/mass_other;

          bgk->alpha_E_fac[i] = ( alpha_E_norm[i] * den_s * den_r / pow(sqrt(vtsq_s+vtsq_r),3.0) ) * pow(sqrt(2.0),app->cdim);
        }
        gkyl_array_release(cross_nu_ho);

        // Set pointers to functions chosen at runtime.
        bgk->cross_nu_func = vmbgk_cross_nu_calc_constNu;
        bgk->alpha_E_func = vmbgk_alpha_E_constNu;
      }
      else {
        // Cross-collision frequency computed in time.
        bgk->norm_nu_cross = true;

        // Cross-collision frequency must be specified symmetrically: if this
        // species uses a computed (Spitzer) cross_nu, its partner must too (i.e.
        // the partner must not provide an explicit cross_nu). The self collision
        // frequency mode is independent and is not checked here.
        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          assert(!(bgk->collide_with[i]->info.collisions.cross_nu[my_idx_in_other[i]]));
        }

        for (int i=0; i<bgk->num_cross_collisions; ++i) {
          double mass_self = vms->mass, mass_other = bgk->collide_with[i]->mass;

          bgk->norm_nu_fac_cross[i] = alpha_E_norm[i]
            * (mass_self+mass_other)/(bgk->delta_sr*bgk->betaGreenep1*mass_self);

          bgk->alpha_E_fac[i] = (bgk->delta_sr*bgk->betaGreenep1*mass_self)/(mass_self+mass_other);
        }

        // Set pointers to functions chosen at runtime.
        bgk->cross_nu_func = vmbgk_cross_nu_calc_normNu;
        bgk->alpha_E_func = vmbgk_alpha_E_normNu;
      }

      // Cross-primitive moment calculator.
      bgk->cross_calc = gkyl_vlasov_cross_prim_moms_bgk_new(&vms->basis, &app->basis, app->use_gpu);

      // Methods chosen at runtime.
      bgk->cross_moms_func = vmbgk_cross_moms_enabled;
    }
  }
}

void
vm_species_bgk_moms(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->moms_func(app, vms, bgk, fin);
}

void
vm_species_bgk_moms_implicit(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->moms_func_implicit(app, vms, bgk, fin);
}

void
vm_species_bgk_rhs(gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  bgk->rhs_func(app, vms, bgk, fin, 0.0, rhs);
}

void
vm_species_bgk_rhs_implicit(gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, double dt, struct gkyl_array *rhs)
{
  bgk->rhs_func_implicit(app, vms, bgk, fin, dt, rhs);
}

void
vm_species_bgk_moms_fixed_temp(gkyl_vlasov_app *app, const struct vm_species *vms,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  bgk->fixed_temp_calc_func(app, vms, bgk, fin);
}

void
vm_species_bgk_write_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  vms->bgk.write_mom_func(app, vms, tm, frame);
}

void 
vm_species_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_bgk_collisions *bgk)
{
  if (bgk->collision_id == GKYL_BGK_COLLISIONS) {
    if (bgk->num_cross_collisions) {
      gkyl_vlasov_cross_prim_moms_bgk_release(bgk->cross_calc);

      for (int i=0; i<bgk->num_cross_collisions; ++i) {
        gkyl_array_release(bgk->cross_nu[i]);
        if (!bgk->norm_nu_cross) {
          gkyl_array_release(bgk->ref_cross_nu[i]);
        }
      }
      gkyl_array_release(bgk->alpha_E);
      gkyl_array_release(bgk->cross_prim_moms);
    }

    gkyl_array_release(bgk->nu_f_lte);
    gkyl_bgk_collisions_release(bgk->up_bgk);

    if (bgk->write_coll_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(bgk->nu_sum_host);
      }
    }

    // ref_self_nu is allocated only for an explicit self_nu; spitzer_calc is
    // allocated whenever any frequency is computed (self or cross). These are
    // now independent (e.g. explicit self_nu with computed cross_nu allocates both).
    if (!bgk->norm_nu_self) {
      gkyl_array_release(bgk->ref_self_nu);
    }
    if (bgk->norm_nu_self || bgk->norm_nu_cross) {
      gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
    }

    if (bgk->fixed_temp_relax) {
      gkyl_array_release(bgk->fixed_temp);
    }

    gkyl_array_release(bgk->nu_sum);
    gkyl_array_release(bgk->self_nu);
  }
}
