#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_moment_priv.h>
#include <gkyl_gr_einstein_excision_priv.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_sources_explicit_priv.h>

static bool
conformal_state_in_excision_for_species(double excision_threshold, const double *q)
{
  return gkyl_gr_einstein_is_excised(excision_threshold, q, 77);
}


static void
conformal_excision_fill_for_species(int meqn, double *q)
{
  for (int m = 0; m < meqn; ++m) {
    q[m] = 0.0;
  }
}

static bool
moment_species_conformal_excision_cleanup(gkyl_moment_app *app,
  const struct moment_species *sp, struct gkyl_array *f)
{
  if (!sp->has_vacuum_einstein_conformal) {
    return false;
  }

  double excision_threshold = sp->vacuum_einstein_conformal_excision_threshold;
  if (!(excision_threshold > 0.0)) {
    return false;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(f, loc);
    if (conformal_state_in_excision_for_species(excision_threshold, q)) {
      conformal_excision_fill_for_species(sp->num_equations, q);
    }
  }
  return true;
}


static void
conformal_gradient_increment(double threshold,
  const double *old, const double *next, double delta[3])
{
  bool live = !conformal_state_in_excision_for_species(threshold, old)
    && !conformal_state_in_excision_for_species(threshold, next);
  for (int j = 0; j < 3; ++j)
    delta[j] = live ? next[65+j]-old[65+j] : 0.0;
}

void
moment_species_conformal_hessian_advance(const gkyl_moment_app *app,
  const struct moment_species *sp, const struct gkyl_array *fin, struct gkyl_array *fout)
{
  if (!sp->has_vacuum_einstein_conformal) return;
  assert(fin != fout);
  double threshold = sp->vacuum_einstein_conformal_excision_threshold;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *old = gkyl_array_cfetch(fin, loc);
    double *next = gkyl_array_fetch(fout, loc);
    if (conformal_state_in_excision_for_species(threshold, next)
        || conformal_state_in_excision_for_species(threshold, old)) {
      for (int m = 68; m < 77; ++m) next[m] = 0.0;
      continue;
    }

    double derivative[3][3] = {{0.0}};
    for (int d = 0; d < app->ndim; ++d) {
      int lo[GKYL_MAX_DIM], hi[GKYL_MAX_DIM];
      gkyl_copy_int_arr(app->ndim, iter.idx, lo);
      gkyl_copy_int_arr(app->ndim, iter.idx, hi);
      lo[d]--; hi[d]++;
      long lm = gkyl_range_idx(&app->local_ext, lo);
      long lp = gkyl_range_idx(&app->local_ext, hi);
      double dm[3], dp[3];
      conformal_gradient_increment(threshold,
        gkyl_array_cfetch(fin, lm), gkyl_array_cfetch(fout, lm), dm);
      conformal_gradient_increment(threshold,
        gkyl_array_cfetch(fin, lp), gkyl_array_cfetch(fout, lp), dp);
      for (int j = 0; j < 3; ++j)
        derivative[d][j] = (dp[j]-dm[j])/(2.0*app->grid.dx[d]);
    }
    // H_t = partial_(i P_j),t. The centered derivative is a shared-face
    // difference: the integrated i-face flux is -average(delta P_j),
    // symmetrized in i,j. On live stencils this preserves the initial
    // H-sym(D_h P), not P-D_h chi. It is not a projection H <- derivative(P).
    // Replace the provisional H transport; adding it would count that update
    // twice. The difference uses Cartesian grid spacings. Excision is a sink.
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        next[68+3*i+j] = old[68+3*i+j]+0.5*(derivative[i][j]+derivative[j][i]);
  }
}

// Differentiate the completed scalar increment, preserving the initial P-D_h chi.
void
moment_species_conformal_derivatives_advance(gkyl_moment_app *app, const struct moment_species *sp,
  const struct gkyl_array *old, struct gkyl_array *next, double time)
{
  if (!sp->has_vacuum_einstein_conformal) return;
  assert(old != next);
  double threshold = sp->vacuum_einstein_conformal_excision_threshold;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long loc = gkyl_range_idx(&app->local_ext, it.idx);
    const double *a = gkyl_array_cfetch(old, loc);
    double *b = gkyl_array_fetch(next, loc);
    bool live = !conformal_state_in_excision_for_species(threshold, a)
      && !conformal_state_in_excision_for_species(threshold, b);
    for (int d=0; d<3; ++d) {
      double increment = 0;
      if (live && d < app->ndim) {
        int lo[GKYL_MAX_DIM], hi[GKYL_MAX_DIM];
        gkyl_copy_int_arr(app->ndim, it.idx, lo);
        gkyl_copy_int_arr(app->ndim, it.idx, hi);
        --lo[d]; ++hi[d];
        double dc[2] = {0};
        const int *indices[2] = {lo, hi};
        for (int side=0; side<2; ++side) {
          long l = gkyl_range_idx(&app->local_ext, indices[side]);
          const double *qa = gkyl_array_cfetch(old, l), *qb = gkyl_array_cfetch(next, l);
          if (!conformal_state_in_excision_for_species(threshold, qa)
              && !conformal_state_in_excision_for_species(threshold, qb))
            dc[side] = qb[64]-qa[64];
        }
        increment = (dc[1]-dc[0])/(2*app->grid.dx[d]);
      }
      b[65+d] = live ? a[65+d]+increment : 0;
    }
  }
  moment_species_apply_bc(app, time, sp, next);
  moment_species_conformal_hessian_advance(app, sp, old, next);
  moment_species_apply_bc(app, time, sp, next);
}

// initialize species
void
moment_species_init(const struct gkyl_moment *mom, const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app, struct moment_species *sp)
{
  sp->ndim = mom->ndim;
  sp->gamma_driver_reference = NULL;
  sp->spacetime_equilibrium = sp->spacetime_equilibrium_source = NULL;
  strcpy(sp->name, mom_sp->name);
  sp->charge = mom_sp->charge;
  sp->mass = mom_sp->mass;
  sp->is_static = mom_sp->is_static; 
  sp->ctx = mom_sp->ctx;
  sp->init = mom_sp->init;

  sp->eqn_type = mom_sp->equation->type;
  sp->num_equations = mom_sp->equation->num_equations;
  sp->equation = gkyl_wv_eqn_acquire(mom_sp->equation);

  // Do we need to update source terms for this fluid?
  // Sources can be electromagnetic fields, closure-related, applied accelerations
  // volume expansion, friction, reactivity, etc. 
  sp->update_sources = false; 
  if (app->has_field) {
    sp->update_sources = true; 
  }

  sp->k0 = 0.0;
  sp->has_grad_closure = false;
  sp->has_nn_closure = false;
  if (mom_sp->equation->type == GKYL_EQN_TEN_MOMENT) {
    sp->update_sources = true; 
    sp->k0 = gkyl_wv_ten_moment_k0(mom_sp->equation);
    sp->has_grad_closure = gkyl_wv_ten_moment_use_grad_closure(mom_sp->equation);
    sp->has_nn_closure = gkyl_wv_ten_moment_use_nn_closure(mom_sp->equation);
    sp->poly_order = gkyl_wv_ten_moment_poly_order(mom_sp->equation);
    sp->ann = gkyl_wv_ten_moment_ann(mom_sp->equation);
  }

  // check if we are running with Braginskii transport and fetch Braginskii type
  if (app->has_braginskii) {
    sp->update_sources = true; 
    sp->type_brag = mom_sp->type_brag;
  }

  if (mom_sp->has_friction) {
    sp->update_sources = true; 
    sp->has_friction = true;
    sp->use_explicit_friction = mom_sp->use_explicit_friction;

    sp->friction_Z = mom_sp->friction_Z;
    sp->friction_T_elc = mom_sp->friction_T_elc;
    sp->friction_Lambda_ee = mom_sp->friction_Lambda_ee;
  }
  else {
    sp->has_friction = false;
    sp->use_explicit_friction = false;
  }

  sp->has_volume_sources = false;
  if (mom_sp->has_volume_sources) {
    sp->update_sources = true; 
    sp->has_volume_sources = true;

    sp->volume_gas_gamma = mom_sp->volume_gas_gamma;
    sp->volume_U0 = mom_sp->volume_U0;
    sp->volume_R0 = mom_sp->volume_R0;
  }

  sp->has_reactivity = false;
  if (mom_sp->has_reactivity) {
    sp->update_sources = true; 
    sp->has_reactivity = true;

    sp->reactivity_gas_gamma = mom_sp->reactivity_gas_gamma;
    sp->reactivity_specific_heat_capacity = mom_sp->reactivity_specific_heat_capacity;
    sp->reactivity_energy_of_formation = mom_sp->reactivity_energy_of_formation;
    sp->reactivity_ignition_temperature = mom_sp->reactivity_ignition_temperature;
    sp->reactivity_reaction_rate = mom_sp->reactivity_reaction_rate;
  }

  sp->has_einstein_medium = false;
  if (mom_sp->has_einstein_medium) {
    sp->update_sources = true; 
    sp->has_einstein_medium = true;

    sp->medium_gas_gamma = mom_sp->medium_gas_gamma;
    sp->medium_kappa = mom_sp->medium_kappa;
  }

  sp->has_gr_ultra_rel = false;
  if (mom_sp->has_gr_ultra_rel) {
    sp->update_sources = true; 
    sp->has_gr_ultra_rel = true;

    sp->gr_ultra_rel_gas_gamma = mom_sp->gr_ultra_rel_gas_gamma;
  }

  sp->has_gr_euler = false;
  if (mom_sp->has_gr_euler) {
    sp->update_sources = true;
    sp->has_gr_euler = true;

    sp->gr_euler_gas_gamma = mom_sp->gr_euler_gas_gamma;
    // The mixed scheme always applies the atmosphere floor in the primitive recovery (fluxes and sources only).
    if (mom->scheme_type == GKYL_MOMENT_MIXED) gkyl_wv_gr_euler_set_atmosphere_floor(true);
  }

  // No metric coupling partner until forward_euler wires it up each RK stage.
  sp->coupling_partner_fin = NULL;

  sp->has_gr_twofluid = false;
  if (mom_sp->has_gr_twofluid) {
    sp->update_sources = true; 
    sp->has_gr_twofluid = true;

    sp->gr_twofluid_mass_elc = mom_sp->gr_twofluid_mass_elc;
    sp->gr_twofluid_mass_ion = mom_sp->gr_twofluid_mass_ion;
    sp->gr_twofluid_charge_elc = mom_sp->gr_twofluid_charge_elc;
    sp->gr_twofluid_charge_ion = mom_sp->gr_twofluid_charge_ion;
    sp->gr_twofluid_gas_gamma_elc = mom_sp->gr_twofluid_gas_gamma_elc;
    sp->gr_twofluid_gas_gamma_ion = mom_sp->gr_twofluid_gas_gamma_ion;
    sp->gr_twofluid_e_fact = mom_sp->gr_twofluid_e_fact;
  }

  sp->has_vacuum_einstein = false;
  if (mom_sp->has_vacuum_einstein) {
    sp->update_sources = true;
    sp->has_vacuum_einstein = true;

    sp->vacuum_einstein_excision_threshold = mom_sp->vacuum_einstein_excision_threshold;
    sp->vacuum_einstein_spacetime_slicing = mom_sp->vacuum_einstein_spacetime_slicing;
    sp->vacuum_einstein_spacetime_evolution = mom_sp->vacuum_einstein_spacetime_evolution;
  }

  sp->has_vacuum_einstein_conformal = false;
  if (mom_sp->has_vacuum_einstein_conformal) {
    sp->update_sources = true;
    sp->has_vacuum_einstein_conformal = true;

    sp->vacuum_einstein_conformal_excision_threshold = mom_sp->vacuum_einstein_conformal_excision_threshold;
    sp->vacuum_einstein_conformal_spacetime_slicing = mom_sp->vacuum_einstein_conformal_spacetime_slicing;
    sp->vacuum_einstein_conformal_spacetime_evolution = mom_sp->vacuum_einstein_conformal_spacetime_evolution;
  }
  
  sp->has_gr_mhd = false;
  if (mom_sp->has_gr_mhd) {
    sp->update_sources = true; 
    sp->has_gr_mhd = true;

    sp->gr_mhd_gas_gamma = mom_sp->gr_mhd_gas_gamma;
  }

  // Per-species scheme override is only consulted (and required) under GKYL_MOMENT_MIXED;
  // otherwise every species inherits the single app-wide scheme, as before.
  sp->scheme_type = (mom->scheme_type == GKYL_MOMENT_MIXED) ? mom_sp->scheme_type : mom->scheme_type;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  enum gkyl_wave_split_type split_type = mom_sp->split_type;    

  int ndim = mom->ndim;
  int meqn = sp->num_equations;  

  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    // create updaters for each directional update
    for (int d=0; d<ndim; ++d)
      sp->slvr[d] = gkyl_wave_prop_new( &(struct gkyl_wave_prop_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .split_type = split_type,
          .limiter = limiter,
          .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
          .force_low_order_flux = mom_sp->force_low_order_flux,
          .check_inv_domain = true,
          .update_dirs = { d },
          .cfl = app->cfl,
          .geom = app->geom,
          .comm = app->comm
        }
      );
      
    sp->fdup = mkarr(false, meqn, app->local_ext.volume);
    // allocate arrays
    for (int d=0; d<ndim+1; ++d)
      sp->f[d] = mkarr(false, meqn, app->local_ext.volume);
  
    // set current solution so ICs and IO work properly
    sp->fcurr = sp->f[0];
  }
  else if ( sp->scheme_type == GKYL_MOMENT_MP || sp->scheme_type == GKYL_MOMENT_KEP ) {
    // determine directions to update
    int num_up_dirs = 0, update_dirs[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      if (!app->is_dir_skipped[d]) {
        update_dirs[num_up_dirs] = d;
        num_up_dirs += 1;
      }

    if (sp->scheme_type == GKYL_MOMENT_MP)
      // single MP updater updates all directions
      sp->mp_slvr = gkyl_mp_scheme_new( &(struct gkyl_mp_scheme_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .mp_recon = app->mp_recon,
          .skip_mp_limiter = mom->skip_mp_limiter,
          .num_up_dirs = num_up_dirs,
          .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
          .cfl = app->cfl,
          .geom = app->geom,
        }
      );
    else
      // single KEP updater updates all directions
      sp->kep_slvr = gkyl_kep_scheme_new( &(struct gkyl_kep_scheme_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .use_hybrid_flux = app->use_hybrid_flux_kep,
          .num_up_dirs = num_up_dirs,
          .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
          .cfl = app->cfl,
          .geom = app->geom,
        }
      );
    
    // allocate arrays
    sp->f0 = mkarr(false, meqn, app->local_ext.volume);
    sp->f1 = mkarr(false, meqn, app->local_ext.volume);
    sp->fnew = mkarr(false, meqn, app->local_ext.volume);
    sp->cflrate = mkarr(false, 1, app->local_ext.volume);
    sp->alpha = mkarr(false, 1, app->local_ext.volume);
    
    // set current solution so ICs and IO work properly
    sp->fcurr = sp->f0;
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    sp->lower_bc[i] = 0;
    sp->upper_bc[i] = 0;
  }

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = mom_sp->bcx;
      else if (dir == 1)
        bc = mom_sp->bcy;
      else
        bc = mom_sp->bcz;

      wv_bc_func_t bc_lower_func;
      if (dir == 0)
        bc_lower_func = mom_sp->bcx_func[0];
      else if (dir == 1)
        bc_lower_func = mom_sp->bcy_func[0];
      else
        bc_lower_func = mom_sp->bcz_func[0];

      wv_bc_func_t bc_upper_func;
      if (dir == 0)
        bc_upper_func = mom_sp->bcx_func[1];
      else if (dir == 1)
        bc_upper_func = mom_sp->bcy_func[1];
      else
        bc_upper_func = mom_sp->bcz_func[1];

      sp->lower_bct[dir] = bc[0];
      sp->upper_bct[dir] = bc[1];

      // lower BCs
      switch (bc[0]) {
        case GKYL_SPECIES_REFLECT:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_lower_func, mom_sp->ctx);
          break;
        
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE: // wedge also uses bc_copy
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_copy, 0);
          break;

        case GKYL_SPECIES_SKIP:
        case GKYL_SPECIES_RADIAL_FALLOFF:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_skip, 0);
          break;

        default:
          assert(false);
          break;
      }
      
      // upper BCs
      switch (bc[1]) {
        case GKYL_SPECIES_REFLECT:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_upper_func, mom_sp->ctx);
          break;
          
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_copy, 0);
          break;

        case GKYL_SPECIES_SKIP:
        case GKYL_SPECIES_RADIAL_FALLOFF:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_skip, 0);
          break;

        default:
          assert(false);
          break;
      }
    }
  }

  // allocate array for applied acceleration/forces for each species
  sp->app_accel = mkarr(false, 3, app->local_ext.volume);
  gkyl_array_clear(sp->app_accel, 0.0);
  sp->has_app_accel = false;
  sp->app_accel_evolve = false;
  if (mom_sp->app_accel) {
    sp->update_sources = true; 
    sp->has_app_accel = true;
    if (mom_sp->app_accel_evolve) {
      sp->app_accel_evolve = mom_sp->app_accel_evolve;
    }
    sp->app_accel_proj = gkyl_fv_proj_new(&app->grid, 2, GKYL_MOM_APP_NUM_APPLIED_ACCELERATION,
      mom_sp->app_accel, mom_sp->app_accel_ctx);  
  }

  sp->embed_mask = mkarr(false, 1, app->local_ext.volume);
  gkyl_array_clear(sp->embed_mask, 1.0);

  if (mom_sp->equation->embed_geo) {
    gkyl_wv_embed_geo_new_mask(mom_sp->equation->embed_geo, &app->grid,
      &app->local, sp->embed_mask);
  }

  sp->nT_source = mkarr(false, 2, app->local_ext.volume);
  sp->nT_source_is_set = false;
  sp->proj_nT_source = 0;
  if (mom_sp->nT_source_func) {
    sp->update_sources = true; 
    void *ctx = sp->ctx;
    if (mom_sp->nT_source_ctx)
      ctx = mom_sp->nT_source_ctx;
    
    sp->proj_nT_source = gkyl_fv_proj_new(&app->grid, 2, GKYL_MOM_APP_NUM_NT_SOURCE,
        mom_sp->nT_source_func, ctx);
    sp->nT_source_set_only_once = mom_sp->nT_source_set_only_once;
  }

  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  sp->bc_buffer = mkarr(false, meqn, buff_sz);

  if (mom_sp->equation->type == GKYL_EQN_EULER)
    sp->integ_q = gkyl_dynvec_new(GKYL_DOUBLE, 6); // KE and PE are stored independently
  else
    sp->integ_q = gkyl_dynvec_new(GKYL_DOUBLE, meqn);
  
  sp->is_first_q_write_call = true;
}

static double
flat_vacuum_einstein_state(enum gkyl_eqn_type type, int c)
{
  double val = 0.0;
  if (c == 0 || c == 4 || c == 8 || c == 9) {
    val = 1.0;
  }
  else if (type == GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL && c == 64) {
    val = 1.0;
  }

  return val;
}



static double
radial_falloff_ghost(double qinf, double skin, double rs, double rg)
{
  double radial_fac = 1.0;
  if (rg > 1.0e-14) {
    radial_fac = rs / rg;
  }

  return qinf + radial_fac * (skin - qinf);
}

static void
moment_species_apply_radial_falloff_bc_corner(gkyl_moment_app *app,
  const struct moment_species *sp, struct gkyl_array *f)
{
  if (sp->equation->type != GKYL_EQN_VACUUM_EINSTEIN &&
      sp->equation->type != GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL) {
    return;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    // Is this a physical-boundary ghost? (any component outside the interior range)
    bool is_ghost = false;
    int donor[GKYL_MAX_DIM];
    for (int d = 0; d < app->ndim; d++) {
      donor[d] = iter.idx[d];
      if (iter.idx[d] < app->local.lower[d]) {
        donor[d] = app->local.lower[d];
        is_ghost = true;
      }
      else if (iter.idx[d] > app->local.upper[d]) {
        donor[d] = app->local.upper[d];
        is_ghost = true;
      }
    }
    if (!is_ghost) {
      continue;
    }

    double xg[3] = { 0.0, 0.0, 0.0 }, xd[3] = { 0.0, 0.0, 0.0 };
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xg);
    gkyl_rect_grid_cell_center(&app->grid, donor, xd);
    double rg = 0.0, rd = 0.0;
    for (int d = 0; d < app->ndim; d++) {
      rg += xg[d] * xg[d];
      rd += xd[d] * xd[d];
    }
    rg = sqrt(rg);
    rd = sqrt(rd);

    const double *donor_q = gkyl_array_cfetch(f, gkyl_range_idx(&app->local_ext, donor));
    double *ghost = gkyl_array_fetch(f, gkyl_range_idx(&app->local_ext, iter.idx));
    for (int c = 0; c < sp->equation->num_equations; c++) {
      double qinf = flat_vacuum_einstein_state(sp->equation->type, c);
      ghost[c] = radial_falloff_ghost(qinf, donor_q[c], rd, rg);
    }
  }
}

static void
moment_species_apply_radial_falloff_bc(gkyl_moment_app *app,
  const struct moment_species *sp, int dir, enum gkyl_edge_loc edge,
  struct gkyl_array *f)
{
  if (sp->equation->type != GKYL_EQN_VACUUM_EINSTEIN &&
      sp->equation->type != GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL) {
    return;
  }

  const struct gkyl_range *ghost_rng = 0;
  if (edge == GKYL_LOWER_EDGE) {
    ghost_rng = &app->skin_ghost.lower_ghost[dir];
  }
  else {
    ghost_rng = &app->skin_ghost.upper_ghost[dir];
  }

  int edge_idx = app->local.lower[dir];
  int fact = -1;
  if (edge == GKYL_UPPER_EDGE) {
    edge_idx = app->local.upper[dir];
    fact = 1;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost_rng);
  while (gkyl_range_iter_next(&iter)) {
    int sidx[GKYL_MAX_DIM];
    for (int d = 0; d < app->ndim; d++) {
      sidx[d] = iter.idx[d];
    }
    sidx[dir] = 2 * edge_idx - iter.idx[dir] + fact;

    double xg[3] = { 0.0, 0.0, 0.0 };
    double xs[3] = { 0.0, 0.0, 0.0 };
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xg);
    gkyl_rect_grid_cell_center(&app->grid, sidx, xs);

    double rg = 0.0, rs = 0.0;
    for (int d = 0; d < app->ndim; d++) {
      rg += xg[d] * xg[d];
      rs += xs[d] * xs[d];
    }
    rg = sqrt(rg);
    rs = sqrt(rs);

    const double *skin = gkyl_array_cfetch(f, gkyl_range_idx(&app->local_ext, sidx));
    double *ghost = gkyl_array_fetch(f, gkyl_range_idx(&app->local_ext, iter.idx));
    for (int c = 0; c < sp->equation->num_equations; c++) {
      double qinf = flat_vacuum_einstein_state(sp->equation->type, c);
      ghost[c] = radial_falloff_ghost(qinf, skin[c], rs, rg);

    }
  }
}

// apply BCs to species
void
moment_species_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_species *sp, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};

  for (int d=0; d<num_periodic_dir; ++d)
    is_non_periodic[app->periodic_dirs[d]] = 0;

  int comm_sz = 1;
  gkyl_comm_get_size(app->comm, &comm_sz);
  bool use_corner_falloff = (comm_sz == 1)
    && (sp->equation->type == GKYL_EQN_VACUUM_EINSTEIN || sp->equation->type == GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL);

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (sp->lower_bct[d] == GKYL_SPECIES_RADIAL_FALLOFF) {
        if (!use_corner_falloff) {
          moment_species_apply_radial_falloff_bc(app, sp, d, GKYL_LOWER_EDGE, f);
        }
      }
      else if (sp->lower_bct[d] != GKYL_SPECIES_WEDGE) {
        gkyl_wv_apply_bc_advance(sp->lower_bc[d], tcurr, &app->local, f);
      }
      if (sp->upper_bct[d] == GKYL_SPECIES_RADIAL_FALLOFF) {
        if (!use_corner_falloff) {
          moment_species_apply_radial_falloff_bc(app, sp, d, GKYL_UPPER_EDGE, f);
        }
      }
      else if (sp->upper_bct[d] != GKYL_SPECIES_WEDGE) {
        gkyl_wv_apply_bc_advance(sp->upper_bc[d], tcurr, &app->local, f);
      }

      // wedge BCs for upper/lower must be handled in one shot
      if (sp->lower_bct[d] == GKYL_SPECIES_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          sp->bc_buffer, d, sp->lower_bc[d], sp->upper_bc[d], f);
    }

  if (use_corner_falloff) {
    moment_species_apply_radial_falloff_bc_corner(app, sp, f);
  }

  // sync interior ghost cells
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);
  // sync periodic ghost cells
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext, num_periodic_dir,
    app->periodic_dirs, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// maximum stable time-step
double
moment_species_max_dt(const gkyl_moment_app *app, const struct moment_species *sp)
{
  double max_dt = DBL_MAX;
  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    for (int d=0; d<app->ndim; ++d)
      max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  }
  else if (sp->scheme_type == GKYL_MOMENT_MP) {
    max_dt = fmin(max_dt, gkyl_mp_scheme_max_dt(sp->mp_slvr, &app->local, sp->f0));
  }
  else if (sp->scheme_type == GKYL_MOMENT_KEP) {
    max_dt = fmin(max_dt, gkyl_kep_scheme_max_dt(sp->kep_slvr, &app->local, sp->f0));
  }  
  return max_dt;
}

// update solution: initial solution is in sp->f[0] and updated
// solution in sp->f[ndim]
static struct gkyl_update_status
moment_species_update_split(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt, bool symmetric)
{
  int ndim = sp->ndim;
  double dt_suggested = DBL_MAX;
  double max_speed = 0.0;
  struct gkyl_wave_prop_status stat;

  if (moment_species_conformal_excision_cleanup(app, sp, sp->f[0])) {
    moment_species_apply_bc(app, tcurr, sp, sp->f[0]);
  }

  int nsweep = symmetric ? 2*ndim-1 : ndim;
  struct gkyl_array *input = sp->f[0];
  for (int sweep=0; sweep<nsweep; ++sweep) {
    int d = sweep < ndim ? sweep : 2*ndim-2-sweep;
    double fraction = symmetric && d != ndim-1 ? 0.5 : 1.0;
    double h = fraction*dt;
    double ts = tcurr + (symmetric && sweep >= ndim ? 0.5*dt : 0.0);
    // Keep f[0] intact: callers use it for rollback and WB subtraction.
    struct gkyl_array *output = symmetric
      ? (input == sp->f[1] ? sp->f[ndim] : sp->f[1]) : sp->f[d+1];
    if (symmetric)
      moment_species_apply_bc(app, ts, sp, input);
    stat = gkyl_wave_prop_advance_equilibrium(sp->slvr[d], ts, h, &app->local,
      sp->embed_mask, input, output, sp->spacetime_equilibrium);

    double my_max_speed = stat.max_speed;
    max_speed = max_speed > my_max_speed ? max_speed : my_max_speed;

    if (!stat.success) {
      return (struct gkyl_update_status) {
        .success = false,
        .dt_actual = 0.0,
        .dt_suggested = stat.dt_suggested/fraction
      };
    }
    
    dt_suggested = fmin(dt_suggested, stat.dt_suggested/fraction);
    moment_species_conformal_excision_cleanup(app, sp, output);
    moment_species_apply_bc(app, symmetric ? ts+h : tcurr, sp, output);
    if (sp->has_vacuum_einstein_conformal) {
      moment_species_conformal_derivatives_advance(app, sp, input, output, symmetric ? ts+h : tcurr);
      moment_species_apply_bc(app, symmetric ? ts+h : tcurr, sp, output);
    }
    input = output;
  }
  if (symmetric && input != sp->f[ndim])
    gkyl_array_copy_range(sp->f[ndim], input, &app->local_ext);

  for (int d=0; d<ndim; ++d) {
    struct gkyl_wv_eqn *eqn = sp->equation;
    if (eqn->type==GKYL_EQN_MHD) {
      gkyl_wv_mhd_set_glm_ch(eqn, max_speed);
    }
  }
  if (app->update_mhd_source) {
    // app->mhd_source.slvr->glm_ch = max_speed;
    // gkyl_mhd_src_set_glm_ch(app->mhd_source.slvr, max_speed);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested
  };
}

struct gkyl_update_status
moment_species_update(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt)
{
  return moment_species_update_split(app, sp, tcurr, dt, false);
}

struct gkyl_update_status
moment_species_update_symmetric(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt)
{
  return moment_species_update_split(app, sp, tcurr, dt, true);
}

// Add the fluid matter source to the vacuum-Einstein RHS (method-of-lines form of the standard ADM matter coupling)
// (Non-static: also reused by the mixed scheme in mom_update_mixed.c; declared in gkyl_moment_priv.h.)
void
add_einstein_matter_source(const double *qe, const double *qf, double gas_gamma, double excision_threshold,
  bool is_conformal, enum gkyl_spacetime_evolution evolution, double *rhs)
{
  double lapse = qe[9];
  if (is_conformal && conformal_state_in_excision_for_species(excision_threshold, qe)) {
    return; // vacuum-Einstein cell is excised: no matter source.
  }
  if (!is_conformal && lapse < excision_threshold) {
    return; // vacuum-Einstein cell is excised: no matter source.
  }
  if (qf[27] < pow(10.0, -8.0)) {
    return; // fluid cell is excised: stress-energy tensor is not meaningful.
  }

  double shift[3] = { qe[52], qe[53], qe[54] };

  double metric_scale = 1.0;
  if (is_conformal) {
    // conformal Einstein state stores chi = psi^-2 and tilde(gamma)_ij = psi^-4 gamma_ij => gamma_ij = chi^-2 tilde(gamma)_ij
    double chi = qe[64];
    if (!(isfinite(chi) && chi > 0.0)) {
      return;
    }
    metric_scale = 1.0 / (chi * chi);
  }

  double spatial_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (is_conformal) {
        spatial_metric[i][j] = metric_scale * qe[(3 * i) + j];
      }
      else {
        spatial_metric[i][j] = qe[(3 * i) + j];
      }
    }
  }

  double shift_lower[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      shift_lower[i] += spatial_metric[i][j] * shift[j];
    }
  }

  double beta_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    beta_sq += shift[i] * shift_lower[i];
  }

  // Explicit contravariant stress-energy tensor T^{mu nu} of the fluid (uses the fluid metric, which
  // has been synced to the vacuum-Einstein metric for this stage, so it is consistent with qe).
  double **T = gkyl_malloc(sizeof(double*[4]));
  for (int a = 0; a < 4; a++) {
    T[a] = gkyl_malloc(sizeof(double[4]));
  }
  gkyl_gr_euler_stress_energy_tensor(gas_gamma, qf, &T);

  // 4D trace T = g_{mu nu} T^{mu nu}, with g_{00} = -alpha^2 + beta^2, g_{0i} = beta_i, g_{ij} = gamma_ij
  double T_trace = (-(lapse * lapse) + beta_sq) * T[0][0];
  for (int a = 0; a < 3; a++) {
    T_trace += 2.0 * shift_lower[a] * T[0][a + 1];
  }
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      T_trace += spatial_metric[a][b] * T[a + 1][b + 1];
    }
  }

  // Spatial projection T_ij = g_{i mu} g_{j nu} T^{mu nu} (both indices lowered) = S_ij
  double T_lower[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double t = shift_lower[i] * shift_lower[j] * T[0][0];
      for (int b = 0; b < 3; b++) {
        t += shift_lower[i] * spatial_metric[j][b] * T[0][b + 1];
      }
      for (int a = 0; a < 3; a++) {
        t += spatial_metric[i][a] * shift_lower[j] * T[a + 1][0];
      }
      for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
          t += spatial_metric[i][a] * spatial_metric[j][b] * T[a + 1][b + 1];
        }
      }
      T_lower[i][j] = t;
    }
  }

  // Mixed projection T^0_k = g_{k nu} T^{0 nu} = beta_k T^{00} + gamma_{kb} T^{0b} = S_k / alpha
  double T_0k[3];
  for (int k = 0; k < 3; k++) {
    double t = shift_lower[k] * T[0][0];
    for (int b = 0; b < 3; b++) {
      t += spatial_metric[k][b] * T[0][b + 1];
    }
    T_0k[k] = t;
  }

  // For conformal, the factors are already present in gamma_ij and therefore in T_ij and T above => no change here
  // dK_ij/dt += -8*pi*alpha*(T_ij - 1/2 gamma_ij T) (standard ADM K_ij matter source)
  // Bona et al. (1997), Eq. (23): n*alpha*gamma_ij/4 * 2*alpha^2*G^{00}.
  // G^{00}=8*pi*T^{00}, E=alpha^2*T^{00}; use the PHYSICAL spatial metric.
  double constraint_matter = evolution == GKYL_EINSTEIN_EVOLUTION
    ? 4.0 * M_PI * lapse * lapse * lapse * T[0][0] : 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rhs[10 + (3 * i) + j] += -8.0 * M_PI * lapse * (T_lower[i][j] - (0.5 * spatial_metric[i][j] * T_trace));
      rhs[10 + (3 * i) + j] += constraint_matter * spatial_metric[i][j];
    }
  }

  // dV_k/dt += +8*pi*alpha^2*T^0_k (momentum-density source for the Bona-Masso auxiliary vector).
  for (int k = 0; k < 3; k++) {
    rhs[49 + k] += 8.0 * M_PI * (lapse * lapse) * T_0k[k];
  }

  for (int a = 0; a < 4; a++) {
    gkyl_free(T[a]);
  }
  gkyl_free(T);
}

// Compute RHS of moment equations
double
moment_species_rhs(gkyl_moment_app *app, struct moment_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec tm = gkyl_wall_clock();
  
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (species->scheme_type == GKYL_MOMENT_MP)
    gkyl_mp_scheme_advance(species->mp_slvr, &app->local, fin,
      app->ql, app->qr, app->amdq, app->apdq,
      species->cflrate, species->embed_mask, rhs);
  else
    gkyl_kep_scheme_advance(species->kep_slvr, &app->local, fin, species->alpha,
      species->cflrate, rhs);

  // Method-of-lines geometric source for GR Euler. 
  if (species->has_gr_euler) {
    const double gamma = species->gr_euler_gas_gamma;
    const int meqn = species->num_equations;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local_ext, iter.idx);
      const double *q = gkyl_array_cfetch(fin, loc);
      double *rhs_c = gkyl_array_fetch(rhs, loc);

      double q_in[71], q_new[71];
      for (int m = 0; m < 71; m++) q_in[m] = q[m];
      explicit_gr_euler_source_update_euler(0, gamma, app->tcurr, 1.0, q_in, q_new); // q_new = q + S(q)
      for (int m = 1; m < 5; m++) rhs_c[m] += (q_new[m] - q_in[m]); // full geometric source

      if (meqn > 71 && !gkyl_gr_euler_wb_disabled(species->equation)) { // WB: subtract the equilibrium source S(q_eq)
        // Use the SAME equilibrium reference the flux subtracts, so the geometric source and the flux are consistent at the equilibrium (matched reference).
        double q_eq[71], q_eq_new[71];
        gkyl_gr_euler_equilibrium(species->equation, q, q_eq);
        explicit_gr_euler_source_update_euler(0, gamma, app->tcurr, 1.0, q_eq, q_eq_new);
        for (int m = 1; m < 5; m++) rhs_c[m] -= (q_eq_new[m] - q_eq[m]);
      }
    }
  }

  if (species->has_vacuum_einstein || species->has_vacuum_einstein_conformal) {
    const int meqn = species->num_equations;

    // If coupled to a GR-Euler fluid, fetch the fluid adiabatic index for the matter source.
    double matter_gas_gamma = 0.0;
    bool is_conformal = species->has_vacuum_einstein_conformal;
    double matter_excision_threshold = 0.0;
    if (is_conformal) {
      matter_excision_threshold = species->vacuum_einstein_conformal_excision_threshold;
    }
    else {
      matter_excision_threshold = species->vacuum_einstein_excision_threshold;
    }
    if (species->coupling_partner_fin) {
      for (int s = 0; s < app->num_species; ++s) {
        if (app->species[s].has_gr_euler) {
          matter_gas_gamma = app->species[s].gr_euler_gas_gamma;
        }
      }
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local_ext, iter.idx);

      const double *q = gkyl_array_cfetch(fin, loc);
      double *rhs_c = gkyl_array_fetch(rhs, loc);

      double sout[meqn];
      gkyl_wv_eqn_source(species->equation, q, sout);

      for (int m = 0; m < meqn; ++m) {
        rhs_c[m] += sout[m];
      }

      // Fluid -> spacetime matter coupling (standard Bona-Masso vacuum-Einstein and conformal share physical lapse, shift, K_ij, and the covariant momentum projection; the conformal branch reconstructs gamma_ij
      if (species->coupling_partner_fin) {
        const double *q_fluid = gkyl_array_cfetch(species->coupling_partner_fin, loc);
        enum gkyl_spacetime_evolution evolution = is_conformal
          ? species->vacuum_einstein_conformal_spacetime_evolution : species->vacuum_einstein_spacetime_evolution;
        add_einstein_matter_source(q, q_fluid, matter_gas_gamma, matter_excision_threshold, is_conformal, evolution, rhs_c);
      }
    }
  }

  double omegaCfl[1];
  gkyl_array_reduce_range(omegaCfl, species->cflrate, GKYL_MAX, &(app->local));


  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl[0];
}

// free species
void
moment_species_release(const struct moment_species *sp)
{
  if (sp->gamma_driver_reference) gkyl_array_release(sp->gamma_driver_reference);
  if (sp->spacetime_equilibrium) gkyl_array_release(sp->spacetime_equilibrium);
  if (sp->spacetime_equilibrium_source) gkyl_array_release(sp->spacetime_equilibrium_source);
  gkyl_wv_eqn_release(sp->equation);
  
  for (int d=0; d<sp->ndim; ++d) {
    if (sp->lower_bc[d])
      gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d])    
      gkyl_wv_apply_bc_release(sp->upper_bc[d]);
  }

  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    for (int d=0; d<sp->ndim; ++d)
      gkyl_wave_prop_release(sp->slvr[d]);
    
    gkyl_array_release(sp->fdup);
    for (int d=0; d<sp->ndim+1; ++d)
      gkyl_array_release(sp->f[d]);
  }
  else if (sp->scheme_type == GKYL_MOMENT_MP || sp->scheme_type == GKYL_MOMENT_KEP) {

    if (sp->scheme_type == GKYL_MOMENT_MP)
      gkyl_mp_scheme_release(sp->mp_slvr);
    else
      gkyl_kep_scheme_release(sp->kep_slvr);
    
    gkyl_array_release(sp->f0);
    gkyl_array_release(sp->f1);
    gkyl_array_release(sp->fnew);
    gkyl_array_release(sp->cflrate);
    gkyl_array_release(sp->alpha);
  }

  gkyl_array_release(sp->app_accel);
  if (sp->has_app_accel) {
    gkyl_fv_proj_release(sp->app_accel_proj);
  }

  gkyl_array_release(sp->nT_source);
  if (sp->proj_nT_source) {
    gkyl_fv_proj_release(sp->proj_nT_source);
  }

  gkyl_array_release(sp->embed_mask);

  gkyl_array_release(sp->bc_buffer);

  gkyl_dynvec_release(sp->integ_q);
}

/** mhd_src functions */

void
mhd_src_init(const struct gkyl_moment_app *app,
             const struct gkyl_moment_species *sp,
             struct mhd_src *src)
{
  double dxyz_min = DBL_MAX;
  for (int d=0; d<app->grid.ndim; ++d) {
    double dx = app->grid.dx[d];
    dxyz_min = dx < dxyz_min ? dx : dxyz_min;
  }

  struct gkyl_mhd_src_inp src_inp = {
    .grid = &app->grid,
    .divergence_constraint = gkyl_wv_mhd_divergence_constraint(sp->equation),
    .glm_ch = gkyl_wv_mhd_glm_ch(sp->equation),
    .glm_alpha = gkyl_wv_mhd_glm_ch(sp->equation),
    .dxyz_min = dxyz_min,
  };

  src->slvr = gkyl_mhd_src_new(src_inp, &app->local_ext);
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
void
mhd_src_update(gkyl_moment_app *app, struct mhd_src *src, int nstrang,
               double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  int i = 0; // mhd has only one 'species'
  struct gkyl_array *fluid = app->species[i].f[sidx[nstrang]];

  if (app->species[i].app_accel_proj)
    gkyl_fv_proj_advance(app->species[i].app_accel_proj, tcurr, &app->local,
                         app->species[i].app_accel);

  // FIXME presently needed for computing divB etc
  moment_species_apply_bc(app, tcurr, &app->species[i], fluid);

  gkyl_mhd_src_advance(src->slvr, dt, &app->local, fluid,
                       app->species[i].app_accel);

  moment_species_apply_bc(app, tcurr, &app->species[i], fluid);

}

void
mhd_src_release(const struct mhd_src *src)
{
  gkyl_mhd_src_release(src->slvr);
}
