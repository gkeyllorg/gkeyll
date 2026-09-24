#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_pkpm_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

// initialize field object
struct pkpm_field* 
pkpm_field_new(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app)
{
  struct pkpm_field *f = gkyl_malloc(sizeof(struct pkpm_field));

  f->info = pkpm->field;

  // allocate EM arrays
  f->em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->em1 = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->emnew = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->em_energy = mkarr(app->use_gpu, 6, app->local_ext.volume);

  // allocate a total field variable for methods which require ext_em + em such as b_hat calculation
  f->tot_em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

  f->em_host = f->em;  
  if (app->use_gpu) {
    f->em_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[6]));
  }

  // Duplicate copy of EM data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  f->em_dup = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_energy_write_call = true;

  // Initialize external EM fields (always used by implicit fluid sources, so always initialize) 
  f->ext_em = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->ext_em, 0.0);
  f->has_ext_em = false;
  f->ext_em_evolve = false;
  // setup external electromagnetic field
  if (f->info.ext_em) {
    f->has_ext_em = true;
    if (f->info.ext_em_evolve) {
      f->ext_em_evolve = f->info.ext_em_evolve;
    }

    f->ext_em_host = f->ext_em;
    if (app->use_gpu) {
      f->ext_em_host = mkarr(false, 6*app->confBasis.num_basis, app->local_ext.volume);
    }
    f->ext_em_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      6, f->info.ext_em, f->info.ext_em_ctx);
  }

  // Initialize applied currents (always used by implicit fluid sources, so always initialize) 
  f->app_current = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->app_current, 0.0);
  f->has_app_current = false;
  f->app_current_evolve = false;
  // setup external currents
  if (f->info.app_current) {
    f->has_app_current = true;
    if (f->info.app_current_evolve) {
      f->app_current_evolve = f->info.app_current_evolve;
    }

    f->app_current_host = f->app_current;
    if (app->use_gpu) {
      f->app_current_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
    }
    f->app_current_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      3, f->info.app_current, f->info.app_current_ctx);
  }

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // equation object
  double c = 1/sqrt(f->info.epsilon0*f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;

  struct gkyl_dg_eqn *eqn;
  eqn = gkyl_dg_maxwell_new(&app->confBasis, c, ef, mf, app->use_gpu);

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[2*GKYL_MAX_DIM] = {0, 0, 0, 0, 0, 0};

  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, eqn,
    app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu);

  // Allocate arrays for diagonstics/parallel-kinetic-perpendicular-moment arrays:
  // bvar = magnetic field unit vector (first 3 components) and unit tensor (last 6 components)
  // ExB = E x B velocity, E x B/|B|^2
  f->cell_avg_magB2 = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);
  f->bvar = mkarr(app->use_gpu, 9*app->confBasis.num_basis, app->local_ext.volume);
  // Surface magnetic field vector organized as:
  // [bx_xl, bx_xr, bxbx_xl, bxbx_xr, bxby_xl, bxby_xr, bxbz_xl, bxbz_xr,
  //  by_yl, by_yr, bxby_yl, bxby_yr, byby_yl, byby_yr, bybz_yl, bybz_yr,
  //  bz_zl, bz_zr, bxbz_zl, bxbz_zr, bybz_zl, bybz_zr, bzbz_zl, bzbz_zr] 
  int cdim = app->cdim;
  int Ncomp_surf = 2*cdim*4;
  int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
  f->bvar_surf = mkarr(app->use_gpu, Ncomp_surf*Nbasis_surf, app->local_ext.volume);
  // Volume expansion of div(b)
  f->div_b = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  // Surface expansion of max b penalization for streaming in PKPM system max(|b_i_l|, |b_i_r|)
  f->max_b = mkarr(app->use_gpu, 2*cdim*Nbasis_surf, app->local_ext.volume);

  // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_em_vars.h)
  double limiter_fac = f->info.limiter_fac == 0 ? 0.0 : f->info.limiter_fac;
  f->limit_em = f->info.limit_em;
  
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, ef, mf, app->use_gpu);
  // Create updaters for bvar (needed by PKPM model)
  f->calc_bvar = gkyl_dg_calc_em_vars_new(&app->grid, &app->confBasis, &app->local_ext, 
    maxwell, app->geom, limiter_fac, 0, app->use_gpu);
  // Create updaters for limiting EM fields
  f->calc_em_vars = gkyl_dg_calc_em_vars_new(&app->grid, &app->confBasis, &app->local_ext, 
    maxwell, app->geom, limiter_fac, 0, app->use_gpu);
  gkyl_wv_eqn_release(maxwell);

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_FIELD_COPY;
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0)
        bc = f->info.bcx;
      else if (dir == 1)
        bc = f->info.bcy;
      else
        bc = f->info.bcz;

      f->lower_bc[dir] = bc[0];
      f->upper_bc[dir] = bc[1];
    }
  }

  // allocate buffer for applying BCs 
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir=0; dir<app->cdim; ++dir) {
    long vol = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, 8*app->confBasis.num_basis, buff_sz);
  
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_FIELD_COPY)
      bctype = GKYL_BC_COPY;
    else if (f->lower_bc[d] == GKYL_FIELD_PEC_WALL)
      bctype = GKYL_BC_MAXWELL_PEC;
    else if (f->lower_bc[d] == GKYL_FIELD_SYM_WALL)
      bctype = GKYL_BC_MAXWELL_SYM;
    else if (f->lower_bc[d] == GKYL_FIELD_RESERVOIR)
      bctype = GKYL_BC_MAXWELL_RESERVOIR;

    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->lower_skin[d], &app->lower_ghost[d], f->em->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_FIELD_COPY)
      bctype = GKYL_BC_COPY;
    else if (f->upper_bc[d] == GKYL_FIELD_PEC_WALL)
      bctype = GKYL_BC_MAXWELL_PEC;
    else if (f->upper_bc[d] == GKYL_FIELD_SYM_WALL)
      bctype = GKYL_BC_MAXWELL_SYM;
    else if (f->upper_bc[d] == GKYL_FIELD_RESERVOIR)
      bctype = GKYL_BC_MAXWELL_RESERVOIR;

    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->upper_skin[d], &app->upper_ghost[d], f->em->ncomp, app->cdim, app->use_gpu);
  }

  gkyl_dg_eqn_release(eqn);

  return f;
}

void
pkpm_field_apply_ic(gkyl_pkpm_app *app, struct pkpm_field *field, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 8, field->info.init, field->info.ctx);

  // run updater; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as magnetic field unit vector computation
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  gkyl_proj_on_basis_advance(proj, t0, &app->local_ext, field->em_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(field->em, field->em_host);

  // pre-compute external EM field and applied current if present
  // pre-computation necessary in case external EM field or applied current
  // are time-independent and not computed in the time-stepping loop
  pkpm_field_calc_ext_em(app, field, t0);
  pkpm_field_calc_app_current(app, field, t0);
}

void
pkpm_field_calc_ext_em(gkyl_pkpm_app *app, struct pkpm_field *field, double tm)
{
  if (field->has_ext_em) {
    gkyl_proj_on_basis_advance(field->ext_em_proj, tm, &app->local_ext, field->ext_em_host);
    if (app->use_gpu) { 
      // note: ext_em_host is same as ext_em when not on GPUs
      gkyl_array_copy(field->ext_em, field->ext_em_host);
    }
  }
}

void
pkpm_field_calc_app_current(gkyl_pkpm_app *app, struct pkpm_field *field, double tm)
{
  if (field->has_app_current) {
    gkyl_proj_on_basis_advance(field->app_current_proj, tm, &app->local_ext, field->app_current_host);
    if (app->use_gpu) { 
      // note: app_current_host is same as app_current when not on GPUs
      gkyl_array_copy(field->app_current, field->app_current_host); 
    }
  }
}

void
pkpm_field_calc_bvar(gkyl_pkpm_app *app, struct pkpm_field *field,
  const struct gkyl_array *em)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(field->tot_em, 0.0);
  gkyl_array_set(field->tot_em, 1.0, em);
  if (field->has_ext_em) {
    gkyl_array_accumulate_range(field->tot_em, 1.0, field->ext_em, &app->local_ext);
  }
  // Assumes magnetic field boundary conditions applied so magnetic field 
  // unit vector and unit tensor are defined everywhere in the domain
  gkyl_dg_calc_em_vars_advance(field->calc_bvar, field->tot_em, 
    field->cell_avg_magB2, field->bvar, field->bvar_surf);

  // Compute div(b) and max_b = max(|b_i_l|, |b_i_r|)
  gkyl_array_clear(field->div_b, 0.0); // Incremented in each dimension, so clear beforehand
  gkyl_dg_calc_em_vars_div_b(field->calc_bvar, &app->local, 
    field->bvar_surf, field->bvar, 
    field->max_b, field->div_b); 

  app->stat.field_em_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_field_limiter(gkyl_pkpm_app *app, struct pkpm_field *field, struct gkyl_array *em)
{
  if (field->limit_em) {
    struct timespec tm = gkyl_wall_clock();

    // Limit the slopes of the solution
    gkyl_dg_calc_em_vars_limiter(field->calc_em_vars, &app->local, em);

    app->stat.field_em_vars_tm += gkyl_time_diff_now_sec(tm);

    // Apply boundary conditions after limiting solution
    pkpm_field_apply_bc(app, field, em);
  }
}

void
pkpm_field_explicit_accumulate_current(gkyl_pkpm_app *app, struct pkpm_field *field, 
  const struct gkyl_array *fluidin[], struct gkyl_array *emout)
{
  struct timespec wst = gkyl_wall_clock();  
  if (!field->info.is_static) {  
    for (int i=0; i<app->num_species; ++i) {
      struct pkpm_species *s = &app->species[i];
      // Need to divide out the mass in pkpm model since we evolve momentum
      double qbymeps = s->info.charge/(s->info.mass*field->info.epsilon0); 
      gkyl_array_accumulate_range(emout, -qbymeps, fluidin[i], &app->local);   
    } 
    // Accumulate applied current to electric field terms
    if (field->has_app_current) {
      gkyl_array_accumulate_range(emout, -1.0/field->info.epsilon0, field->app_current, &app->local);
    }
  }
  app->stat.current_tm += gkyl_time_diff_now_sec(wst);  
}

// Compute the RHS for field update, returning maximum stable
// time-step.
double
pkpm_field_rhs(gkyl_pkpm_app *app, struct pkpm_field *field,
  const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  double omegaCfl = 1/DBL_MAX;
  
  gkyl_array_clear(field->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (!field->info.is_static) {
    gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);
    
    gkyl_array_reduce_range(field->omegaCfl_ptr, field->cflrate, GKYL_MAX, &app->local);

    app->stat.n_field_omega_cfl += 1;
    struct timespec tm = gkyl_wall_clock();
    
    double omegaCfl_ho[1];
    if (app->use_gpu)
      gkyl_cu_memcpy(omegaCfl_ho, field->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
    else
      omegaCfl_ho[0] = field->omegaCfl_ptr[0];
    omegaCfl = omegaCfl_ho[0];

    app->stat.field_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  }

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
  
  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for EM fields
void
pkpm_field_apply_bc(gkyl_pkpm_app *app, const struct pkpm_field *field, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();  
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, app->periodic_dirs, f);
  
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (field->lower_bc[d]) {
        case GKYL_FIELD_COPY:
        case GKYL_FIELD_PEC_WALL:
        case GKYL_FIELD_SYM_WALL:
        case GKYL_FIELD_RESERVOIR:
          gkyl_bc_basic_advance(field->bc_lo[d], field->bc_buffer, f);
          break;

        default:
          break;
      }

      switch (field->upper_bc[d]) {
        case GKYL_FIELD_COPY:
        case GKYL_FIELD_PEC_WALL:
        case GKYL_FIELD_SYM_WALL:
        case GKYL_FIELD_RESERVOIR:
          gkyl_bc_basic_advance(field->bc_up[d], field->bc_buffer, f);
          break;
          
        default:
          break;
      }   
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);

  app->stat.field_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_field_calc_energy(gkyl_pkpm_app *app, double tm, const struct pkpm_field *field)
{
  for (int i=0; i<6; ++i)
    gkyl_dg_calc_l2_range(&app->confBasis, i, field->em_energy, i, field->em, app->local);
  gkyl_array_scale_range(field->em_energy, app->grid.cellVolume, &app->local);
  
  double energy[6] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(field->em_energy_red, field->em_energy, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(energy, field->em_energy_red, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(energy, field->em_energy, GKYL_SUM, &app->local);
  }

  double energy_global[6] = { 0.0 };
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);
  
  gkyl_dynvec_append(field->integ_energy, tm, energy_global);
}

// release resources for field
void
pkpm_field_release(const gkyl_pkpm_app* app, struct pkpm_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
  gkyl_array_release(f->tot_em);
  gkyl_array_release(f->em_dup);
  
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  gkyl_array_release(f->em_energy);
  gkyl_dynvec_release(f->integ_energy);

  gkyl_array_release(f->cell_avg_magB2);
  gkyl_array_release(f->bvar);
  gkyl_array_release(f->bvar_surf);
  gkyl_array_release(f->div_b);
  gkyl_array_release(f->max_b);
  gkyl_dg_calc_em_vars_release(f->calc_bvar);
  gkyl_dg_calc_em_vars_release(f->calc_em_vars);

  gkyl_array_release(f->ext_em);
  if (f->has_ext_em) {
    if (app->use_gpu) {
      gkyl_array_release(f->ext_em_host);
    }
    gkyl_proj_on_basis_release(f->ext_em_proj);
  }
  gkyl_array_release(f->app_current);
  if (f->has_app_current) {
    if (app->use_gpu) {
      gkyl_array_release(f->app_current_host);
    }
    gkyl_proj_on_basis_release(f->app_current_proj);
  }

  gkyl_hyper_dg_release(f->slvr);

  if (app->use_gpu) {
    gkyl_array_release(f->em_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    gkyl_cu_free(f->em_energy_red);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }

  gkyl_free(f);
}

