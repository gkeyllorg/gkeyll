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
struct pkpm_field *
pkpm_field_new(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app)
{
  struct pkpm_field *f = gkyl_malloc(sizeof(struct pkpm_field));

  f->info = pkpm->field;

  // allocate EM arrays
  f->em = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->em1 = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->emnew = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->em_energy = mkarr(app->use_gpu, 6, app->local_ext.volume);

  // allocate a total field variable for methods which require ext_em + em such as b_hat calculation
  f->tot_em = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);

  f->em_host = f->em;
  if (app->use_gpu) {
    f->em_host = mkarr(false, 8 * app->confBasis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[6]));
  }
  // The evolved field carries the position-map Jacobian (J E, J B: the convention of the
  // mapped Maxwell kernels); em_no_J holds the physical fields of the stage for the
  // magnetic unit vector, the sources and the output. Identity map: a copy.
  f->em_no_J = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->em_no_J_host = f->em_no_J;
  if (app->use_gpu) {
    f->em_no_J_host = mkarr(false, 8 * app->confBasis.num_basis, app->local_ext.volume);
  }
  f->pos_map = gkyl_vlasov_position_map_acquire(app->pos_map);

  // Duplicate copy of EM data in case time step fails.
  // Needed because of implicit source split which modifies solution and
  // is always successful, so if a time step fails due to the SSP RK3
  // we must restore the old solution before restarting the time step
  f->em_dup = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_energy_write_call = true;

  // Initialize external EM fields (always used by implicit fluid sources, so always initialize)
  f->ext_em = mkarr(app->use_gpu, 6 * app->confBasis.num_basis, app->local_ext.volume);
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
      f->ext_em_host = mkarr(false, 6 * app->confBasis.num_basis, app->local_ext.volume);
    }
    f->ext_em_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->confBasis,
      .num_quad = app->confBasis.poly_order + 1,
      .num_ret_vals = 6,
      .eval = f->info.ext_em,
      .ctx = f->info.ext_em_ctx,
      .c2p_func = pkpm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // Initialize applied currents (always used by implicit fluid sources, so always initialize)
  f->app_current = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->app_current, 0.0);
  // app_current drives the J-weighted field equation and carries J; the implicit coupling
  // works with the physical current.
  f->app_current_no_J = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->app_current_no_J, 0.0);
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
      f->app_current_host = mkarr(false, 3 * app->confBasis.num_basis, app->local_ext.volume);
    }
    f->app_current_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->confBasis,
      .num_quad = app->confBasis.poly_order + 1,
      .num_ret_vals = 3,
      .eval = f->info.app_current,
      .ctx = f->info.app_current_ctx,
      .c2p_func = pkpm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu) {
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  } else {
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));
  }

  // equation object
  f->use_ghost_current = false;
  if (f->info.use_ghost_current) {
    // 1D periodic electrostatics only: the mean current is a single number.
    assert(app->cdim == 1 && app->num_periodic_dir == 1);
    f->use_ghost_current = true;
    f->ghost_current = mkarr(app->use_gpu, 1, app->local_ext.volume);
    gkyl_array_clear(f->ghost_current, 0.0);
    f->app_current_tot = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu) {
      f->red_ghost_current = gkyl_cu_malloc(sizeof(double[1]));
    }
  }

  double c = 1 / sqrt(f->info.epsilon0 * f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;

  struct gkyl_dg_eqn *eqn;

  // Input structure for building the dg eqn object: the Maxwell curl kernels take the
  // position-map Jacobian of the app (1 everywhere for the identity map).
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &app->confBasis,
    .crange = &app->local,
    .jacob_pos = f->pos_map->jacob_pos,
    .conf_flux_surf = 0,
    .lightSpeed = c,
    .field_id = GKYL_FIELD_E_B,
    .elcErrorSpeedFactor = ef,
    .mgnErrorSpeedFactor = mf,
    .use_gpu = app->use_gpu,
  };
  eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[2 * GKYL_MAX_DIM] = {0, 0, 0, 0, 0, 0};

  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(
    &app->grid, &app->confBasis, eqn, app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu
  );

  // Magnetic field unit vector b = B/|B| (p=1) by nodal normalization, for the PKPM model.
  f->bvar = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  f->grad_b = mkarr(app->use_gpu, app->cdim * 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->grad_b, 0.0);
  // operator split: the frozen transport right-hand side and the curl B current of the source
  // half steps (plan 17.29)
  f->rhs_split = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->omega_cfl_split = 1.0 / DBL_MAX;
  f->curl_rhs = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->curl_rhs_no_J = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  f->curl_current = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  // physical time derivative of the fields of the stage (zero for a static field)
  f->dt_em_no_J = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->dt_em_no_J, 0.0);
  f->dtB_over_magB = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->dtB_over_magB, 0.0);
  f->calc_bvar = gkyl_dg_calc_pkpm_bvar_inew(&(struct gkyl_dg_calc_pkpm_bvar_inp){
    .conf_grid = &app->grid,
    .conf_basis = &app->confBasis,
    .pos_map = app->pos_map,
    .use_gpu = app->use_gpu,
  });

  // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_em_vars.h)
  double limiter_fac = f->info.limiter_fac == 0 ? 0.0 : f->info.limiter_fac;
  f->limit_em = f->info.limit_em;
  if (f->limit_em) {
    f->em_lim = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);
  }

  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, ef, mf, app->use_gpu);
  // Create updaters for limiting EM fields
  f->calc_em_vars = gkyl_dg_calc_em_vars_new(
    &app->grid, &app->confBasis, &app->local_ext, maxwell, app->geom, limiter_fac, 0, app->use_gpu
  );
  gkyl_wv_eqn_release(maxwell);

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  for (int dir = 0; dir < app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_FIELD_COPY;
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0) {
        bc = f->info.bcx;
      } else if (dir == 1) {
        bc = f->info.bcy;
      } else {
        bc = f->info.bcz;
      }

      f->lower_bc[dir] = bc[0];
      f->upper_bc[dir] = bc[1];
    }
  }

  // allocate buffer for applying BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir = 0; dir < app->cdim; ++dir) {
    long vol = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, buff_sz);

  for (int d = 0; d < app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_FIELD_COPY) {
      bctype = GKYL_BC_COPY;
    } else if (f->lower_bc[d] == GKYL_FIELD_PEC_WALL) {
      bctype = GKYL_BC_MAXWELL_PEC;
    } else if (f->lower_bc[d] == GKYL_FIELD_SYM_WALL) {
      bctype = GKYL_BC_MAXWELL_SYM;
    } else if (f->lower_bc[d] == GKYL_FIELD_RESERVOIR) {
      bctype = GKYL_BC_MAXWELL_RESERVOIR;
    }

    f->bc_lo[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.confBasis, &app->lower_skin[d],
      &app->lower_ghost[d], f->em->ncomp, app->cdim, app->use_gpu
    );

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_FIELD_COPY) {
      bctype = GKYL_BC_COPY;
    } else if (f->upper_bc[d] == GKYL_FIELD_PEC_WALL) {
      bctype = GKYL_BC_MAXWELL_PEC;
    } else if (f->upper_bc[d] == GKYL_FIELD_SYM_WALL) {
      bctype = GKYL_BC_MAXWELL_SYM;
    } else if (f->upper_bc[d] == GKYL_FIELD_RESERVOIR) {
      bctype = GKYL_BC_MAXWELL_RESERVOIR;
    }

    f->bc_up[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.confBasis, &app->upper_skin[d],
      &app->upper_ghost[d], f->em->ncomp, app->cdim, app->use_gpu
    );
  }

  gkyl_dg_eqn_release(eqn);

  return f;
}

void
pkpm_field_apply_ic(gkyl_pkpm_app *app, struct pkpm_field *field, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->confBasis,
    .num_quad = poly_order + 1,
    .num_ret_vals = 8,
    .eval = field->info.init,
    .ctx = field->info.ctx,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });

  // run updater; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as magnetic field unit vector computation
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  // The user gives the physical fields at the physical coordinates of the map; the
  // evolved field carries the position-map Jacobian.
  gkyl_proj_on_basis_advance(proj, t0, &app->local_ext, field->em_no_J_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu) {
    gkyl_array_copy(field->em_no_J, field->em_no_J_host);
  }
  gkyl_vlasov_position_map_rescale_jacobpos_conf(app->pos_map, &app->local_ext, field->em_no_J, field->em);

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
    gkyl_proj_on_basis_advance(
      field->app_current_proj, tm, &app->local_ext, field->app_current_host
    );
    if (app->use_gpu) {
      // note: app_current_host is same as app_current when not on GPUs
      gkyl_array_copy(field->app_current, field->app_current_host);
    }
    // Physical current for the implicit coupling; the field equation (J E) is driven by J J_app.
    gkyl_array_copy(field->app_current_no_J, field->app_current);
    gkyl_vlasov_position_map_rescale_jacobpos_conf(
      app->pos_map, &app->local_ext, field->app_current_no_J, field->app_current
    );
  }
}

void
pkpm_field_calc_bvar(gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *em)
{
  struct timespec tm = gkyl_wall_clock();

  // Physical fields of the stage (the evolved em carries the position-map Jacobian).
  gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &app->local_ext, em->ncomp, em, field->em_no_J);
  gkyl_array_clear(field->tot_em, 0.0);
  gkyl_array_set(field->tot_em, 1.0, field->em_no_J);
  if (field->has_ext_em) {
    gkyl_array_accumulate_range(field->tot_em, 1.0, field->ext_em, &app->local_ext);
  }
  // Assumes magnetic field boundary conditions applied so the unit vector is
  // defined everywhere in the extended domain.
  gkyl_dg_calc_pkpm_bvar_advance(field->calc_bvar, &app->local_ext, field->tot_em, field->bvar);
  // Weak gradient of b for the mirror terms (div b), from b in the ghost cells.
  gkyl_dg_calc_pkpm_bvar_grad(field->calc_bvar, &app->local, field->bvar, field->grad_b);

  app->stat.field_em_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_field_limiter(gkyl_pkpm_app *app, struct pkpm_field *field, struct gkyl_array *em)
{
  if (field->limit_em) {
    struct timespec tm = gkyl_wall_clock();

    // Limit the slopes of the PHYSICAL field (on a position map the J-weighted field jumps with J
    // between cells; on the identity map J = 1 exactly and this is the field itself)
    gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &app->local_ext, em->ncomp, em, field->em_lim);
    gkyl_dg_calc_em_vars_limiter(field->calc_em_vars, &app->local, field->em_lim);
    gkyl_vlasov_position_map_rescale_jacobpos_conf(app->pos_map, &app->local, field->em_lim, em);

    app->stat.field_em_vars_tm += gkyl_time_diff_now_sec(tm);

    // Apply boundary conditions after limiting solution
    pkpm_field_apply_bc(app, field, em);
  }
}

void
pkpm_field_calc_ghost_current(
  gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *fluidin[]
)
{
  int nc = app->confBasis.num_basis;
  double mean_local = 0.0, mean_global = 0.0;
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    double qbym_nx = s->info.charge / (s->info.mass * app->grid.cells[0]);
    // ghost_current has one component: coefficient 0 of M_x, scaled by q/(m N_x).
    gkyl_array_set_offset(s->pkpm_mom3, 1.0, fluidin[i], nc);
    gkyl_array_set_range(field->ghost_current, qbym_nx, s->pkpm_mom3, &app->local);
    double red[1];
    if (app->use_gpu) {
      gkyl_array_reduce_range(field->red_ghost_current, field->ghost_current, GKYL_SUM, &app->local);
      gkyl_cu_memcpy(red, field->red_ghost_current, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      gkyl_array_reduce_range(red, field->ghost_current, GKYL_SUM, &app->local);
    }
    mean_local += red[0];
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &mean_local, &mean_global);
  gkyl_array_clear(field->ghost_current, mean_global);
}

void
pkpm_field_explicit_accumulate_current(
  gkyl_pkpm_app *app, struct pkpm_field *field, const struct gkyl_array *fluidin[],
  struct gkyl_array *emout
)
{
  struct timespec wst = gkyl_wall_clock();
  if (!field->info.is_static) {
    for (int i = 0; i < app->num_species; ++i) {
      struct pkpm_species *s = &app->species[i];
      // Need to divide out the mass in pkpm model since we evolve momentum. The fluid
      // array is [rho, rho ux, rho uy, rho uz, E]: only its momentum block is the current.
      double qbymeps = s->info.charge / (s->info.mass * field->info.epsilon0);
      gkyl_array_set_offset(s->pkpm_mom3, 1.0, fluidin[i], app->confBasis.num_basis);
      gkyl_array_accumulate_range(emout, -qbymeps, s->pkpm_mom3, &app->local);
    }
    // Ghost current: add back the domain-mean current (coefficient 0 of E_x), J-weighted like
    // the evolved field (a uniform current is then no source of E on a position map).
    if (field->use_ghost_current) {
      pkpm_field_calc_ghost_current(app, field, fluidin);
      gkyl_vlasov_position_map_rescale_jacobpos_conf(app->pos_map, &app->local, field->ghost_current, field->ghost_current);
      gkyl_array_accumulate_range(emout, 1.0 / field->info.epsilon0, field->ghost_current, &app->local);
    }
    // Accumulate applied current to electric field terms
    if (field->has_app_current) {
      gkyl_array_accumulate_range(
        emout, -1.0 / field->info.epsilon0, field->app_current, &app->local
      );
    }
  }
  app->stat.current_tm += gkyl_time_diff_now_sec(wst);
}

// Compute the RHS for field update, returning maximum stable
// time-step.
double
pkpm_field_rhs(
  gkyl_pkpm_app *app, struct pkpm_field *field, int stage, const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();

  double omegaCfl = 1 / DBL_MAX;

  gkyl_array_clear(field->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (!field->info.is_static && !app->use_explicit_source && stage > 1) {
    // operator split: B drifts with the E of the source half kick over the whole transport step
    gkyl_array_copy(rhs, field->rhs_split);
    gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &app->local, rhs->ncomp, rhs, field->dt_em_no_J);
    gkyl_dg_calc_pkpm_bvar_dtB_over_magB(field->calc_bvar, &app->local, field->tot_em, field->dt_em_no_J, field->dtB_over_magB);
    app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
    return app->cfl / field->omega_cfl_split;
  }

  if (!field->info.is_static) {
    gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);
    if (!app->use_explicit_source) {
      // operator split: E is advanced only by the implicit source half steps (Stormer-Verlet in the
      // curl terms, plan 17.29); zero the E rows (the first 3 components) and keep the right-hand
      // side for stages 2 and 3
      gkyl_array_set_offset(rhs, 0.0, field->curl_current, 0);
      gkyl_array_copy(field->rhs_split, rhs);
    }
    // physical time derivative of the fields (the evolved em carries the position-map Jacobian)
    // and d_t B/|B| at the nodes of b (tot_em holds the physical fields of the stage)
    gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &app->local, rhs->ncomp, rhs, field->dt_em_no_J);
    gkyl_dg_calc_pkpm_bvar_dtB_over_magB(field->calc_bvar, &app->local, field->tot_em, field->dt_em_no_J, field->dtB_over_magB);

    gkyl_array_reduce_range(field->omegaCfl_ptr, field->cflrate, GKYL_MAX, &app->local);

    app->stat.n_field_omega_cfl += 1;
    struct timespec tm = gkyl_wall_clock();

    double omegaCfl_ho[1];
    if (app->use_gpu) {
      gkyl_cu_memcpy(omegaCfl_ho, field->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
    } else {
      omegaCfl_ho[0] = field->omegaCfl_ptr[0];
    }
    omegaCfl = omegaCfl_ho[0];
    field->omega_cfl_split = omegaCfl;

    app->stat.field_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  }

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);

  return app->cfl / omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for EM fields
void
pkpm_field_apply_bc(gkyl_pkpm_app *app, const struct pkpm_field *field, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    app->comm, &app->local, &app->local_ext, num_periodic_dir, app->periodic_dirs, f
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
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
  for (int i = 0; i < 6; ++i) {
    gkyl_dg_calc_l2_range(&app->confBasis, i, field->em_energy, i, field->em, app->local);
  }
  gkyl_array_scale_range(field->em_energy, app->grid.cellVolume, &app->local);
  // The field stores J E, J B: the L2 carries J^2 and the computational cell volume
  // under-counts the physical one by J; divide the net extra J out (identity map: no-op).
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local, field->em_energy->ncomp, field->em_energy, field->em_energy
  );

  double energy[6] = {0.0};
  if (app->use_gpu) {
    gkyl_array_reduce_range(field->em_energy_red, field->em_energy, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(energy, field->em_energy_red, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(energy, field->em_energy, GKYL_SUM, &app->local);
  }

  double energy_global[6] = {0.0};
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);
}

// release resources for field
void
pkpm_field_release(const gkyl_pkpm_app *app, struct pkpm_field *f)
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

  gkyl_array_release(f->bvar);
  gkyl_array_release(f->grad_b);
  gkyl_array_release(f->rhs_split);
  gkyl_array_release(f->curl_rhs);
  gkyl_array_release(f->curl_rhs_no_J);
  gkyl_array_release(f->curl_current);
  gkyl_array_release(f->dt_em_no_J);
  gkyl_array_release(f->dtB_over_magB);
  gkyl_dg_calc_pkpm_bvar_release(f->calc_bvar);
  gkyl_dg_calc_em_vars_release(f->calc_em_vars);

  gkyl_array_release(f->ext_em);
  if (f->has_ext_em) {
    if (app->use_gpu) {
      gkyl_array_release(f->ext_em_host);
    }
    gkyl_proj_on_basis_release(f->ext_em_proj);
  }
  gkyl_array_release(f->app_current);
  if (f->use_ghost_current) {
    gkyl_array_release(f->ghost_current);
    gkyl_array_release(f->app_current_tot);
    if (app->use_gpu) {
      gkyl_cu_free(f->red_ghost_current);
    }
  }
  gkyl_array_release(f->app_current_no_J);
  gkyl_array_release(f->em_no_J);
  if (f->limit_em) {
    gkyl_array_release(f->em_lim);
  }
  if (app->use_gpu) {
    gkyl_array_release(f->em_no_J_host);
  }
  if (f->has_app_current) {
    if (app->use_gpu) {
      gkyl_array_release(f->app_current_host);
    }
    gkyl_proj_on_basis_release(f->app_current_proj);
  }

  gkyl_hyper_dg_release(f->slvr);
  gkyl_vlasov_position_map_release(f->pos_map);

  if (app->use_gpu) {
    gkyl_array_release(f->em_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    gkyl_cu_free(f->em_energy_red);
  } else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }

  gkyl_free(f);
}
