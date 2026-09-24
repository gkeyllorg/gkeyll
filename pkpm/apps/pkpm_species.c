#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_pkpm_priv.h>

#include <assert.h>
#include <time.h>

// initialize species object
void
pkpm_species_init(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app, struct pkpm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = pkpm->cells[d];
    lower[d] = pkpm->lower[d];
    upper[d] = pkpm->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->info.cells[d];
    lower_vel[d] = s->info.lower[d];
    upper_vel[d] = s->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->global_ext, &s->global);
  
  // velocity space grid
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  // phase-space communicator
  s->comm = gkyl_comm_extend_comm(app->comm, &s->local_vel);

  // create local and local_ext from app local range
  struct gkyl_range local;
  // local = conf-local X local_vel
  gkyl_range_ten_prod(&local, &app->local, &s->local_vel);
  gkyl_create_ranges(&local, ghost, &s->local_ext, &s->local);

  // allocate distribution function arrays
  // PKPM has two distribution functions F_0 and T_perp*G, first two Laguerre moments
  s->f = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
  // allocate momentum arrays (rho ux, rho uy, rho uz)
  s->fluid = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  s->fluid1 = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  s->fluidnew = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);

  // create host arrays if on GPUs so initialization occurs host-side
  s->f_host = s->f;
  s->fluid_host = s->fluid;
  if (app->use_gpu) {
    s->f_host = mkarr(false, 2*app->basis.num_basis, s->local_ext.volume);
    s->fluid_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
  }

  // Duplicate copy of momentum data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  s->fluid_dup = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);

  // Wave equation object for upwinding fluid equations
  // We use the 10 moment system since PKPM model generates a full pressure tensor
  // k0 = 0.0 and use_grad_closure = false because we do not need a closure 
  // (the kinetic equations are the closure)
  s->equation = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, app->use_gpu); 

  // Distribution function arrays for coupling different Laguerre moments
  // g_dist_source has two components: 
  // [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
  // (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0]
  // First component is mirror force source *distribution*, second component is *total* vperp characteristics source.
  s->g_dist_source = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume); 
  s->F_k_p_1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume); 
  s->F_k_m_1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume); 

  // allocate cflrate for distribution function and momentum updates (scalar arrays)
  s->cflrate_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  s->cflrate_fluid = mkarr(app->use_gpu, 1, app->local_ext.volume);

  if (app->use_gpu) {
    s->omegaCfl_ptr_dist = gkyl_cu_malloc(sizeof(double));
    s->omegaCfl_ptr_fluid = gkyl_cu_malloc(sizeof(double));
  }
  else {
    s->omegaCfl_ptr_dist = gkyl_malloc(sizeof(double));
    s->omegaCfl_ptr_fluid = gkyl_malloc(sizeof(double));
  }

  // allocate array to store q/m*(E,B) if using a fully explicit update
  s->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

  // pkpm moments for update (rho, p_par, p_perp, M1)
  pkpm_species_moment_init(app, s, &s->pkpm_moms, false);
  // pkpm moments for diagnostics (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
  // For simplicity, is_integrated flag also used by PKPM to turn on diagnostics
  pkpm_species_moment_init(app, s, &s->pkpm_moms_diag, true);

  // div(p_par b_hat), for self-consistent total pressure force
  s->pkpm_div_ppar = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  // allocate array to store primitive moments : 
  // [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp, 3*T_xx/m, 3*T_yy/m, 3*T_zz/m]
  // pressure p_ij : (p_par - p_perp) b_i b_j + p_perp g_ij
  s->pkpm_prim = mkarr(app->use_gpu, 9*app->confBasis.num_basis, app->local_ext.volume);
  s->pkpm_u = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  s->pkpm_p_ij = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
  // boolean array for primitive variables [rho, p_par, p_perp] is negative at control points
  // *only used for diagnostic purposes*
  s->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

  int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
  // Surface primitive variables (2*cdim*4 components). Ordered as:
  // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 3.0*Txx_xl/m, 3.0*Txx_xr/m, 
  //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 
  //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
  s->pkpm_prim_surf = mkarr(app->use_gpu, 2*cdim*4*Nbasis_surf, app->local_ext.volume);
  // Surface expansions of:
  // 1. Lax penalization maximum speed lambda_i = |u_i| + sqrt(3*P_ii/rho)
  // 2. Momentum solve penalization (either 10 moment Roe solve or Lax penalization)
  //    Momentum penalization returns the full penalization, e.g., the modal surface expansion of
  //    jump in the fluctuations from the 10 moment Roe solve 1/2 (A^+ Delta Q - A^- Delta Q)
  s->pkpm_lax = mkarr(app->use_gpu, cdim*Nbasis_surf, app->local_ext.volume);
  s->pkpm_penalization = mkarr(app->use_gpu, 3*cdim*Nbasis_surf, app->local_ext.volume);

  // allocate array for pkpm acceleration variables, stored in pkpm_accel: 
  // 0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
  // 1: bb_grad_u (bb : grad(u))
  // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
  // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2*nu)
  s->pkpm_accel = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume); 

  // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_pkpm_vars.h)
  double limiter_fac = s->info.limiter_fac == 0 ? 0.0 : s->info.limiter_fac;
  s->limit_fluid = s->info.limit_fluid;

  // updater for computing pkpm variables 
  // pressure, primitive variables, and acceleration variables
  // also stores kernels for computing source terms, integrated variables
  // Two instances, one over extended range and one over local range for ease of handling boundary conditions
  s->calc_pkpm_vars_ext = gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, &app->local_ext, 
    s->equation, app->geom, limiter_fac, app->use_gpu);
  s->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, &app->local, 
    s->equation, app->geom, limiter_fac, app->use_gpu); 
  // updater for computing pkpm distribution function variables
  // div(p_par b_hat) for self-consistent total pressure force and distribution function sources for
  // Laguerre couplings and vperp characteristics 
  s->calc_pkpm_dist_vars = gkyl_dg_calc_pkpm_dist_vars_new(&s->grid, &app->confBasis, app->use_gpu);

  // array for storing integrated fluid variables in each cell
  // integ_pkpm_mom : integral[rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_par, p_perp]
  s->integ_pkpm_mom = mkarr(app->use_gpu, 9, app->local_ext.volume);
  // arrays for I/O, fluid_io and pkpm_vars_io
  // fluid_io : [rho, rhoux, rhouy, rhouz, P_xx+rhoux^2, P_xy+rhouxuy, P_xz+rhouxuz, P_yy+rhouy^2, P_yz+rhouyuz, P_zz+rhouz^2]
  // pkpm_vars_io : [ux, uy, uz, T_perp/m, m/T_perp, 1/rho div(p_par b), p_perp/rho div(b), bb : grad(u)]
  s->fluid_io = mkarr(app->use_gpu, 10*app->confBasis.num_basis, app->local_ext.volume);
  s->pkpm_vars_io = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  s->fluid_io_host = s->fluid_io;
  s->pkpm_vars_io_host = s->pkpm_vars_io;
  if (app->use_gpu) {
    s->fluid_io_host = mkarr(false, 10*app->confBasis.num_basis, app->local_ext.volume);
    s->pkpm_vars_io_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);
  }

  // by default, we do not have zero-flux boundary conditions in any direction
  bool is_zero_flux[GKYL_MAX_DIM] = {false};

  struct gkyl_dg_vlasov_pkpm_auxfields vlasov_pkpm_inp = {.bvar = app->field->bvar, .bvar_surf = app->field->bvar_surf, 
    .pkpm_prim = s->pkpm_prim, .pkpm_prim_surf = s->pkpm_prim_surf, 
    .max_b = app->field->max_b, .pkpm_lax = s->pkpm_lax, 
    .div_b = app->field->div_b, .pkpm_accel_vars = s->pkpm_accel, 
    .g_dist_source = s->g_dist_source};
  struct gkyl_dg_euler_pkpm_auxfields euler_pkpm_inp = {.vlasov_pkpm_moms = s->pkpm_moms.marr, 
    .pkpm_prim = s->pkpm_prim, .pkpm_prim_surf = s->pkpm_prim_surf, .pkpm_p_ij = s->pkpm_p_ij, 
    .pkpm_lax = s->pkpm_lax, .pkpm_penalization = s->pkpm_penalization};
  // create solver
  s->slvr = gkyl_dg_updater_pkpm_new(&app->grid, &s->grid, 
    &app->confBasis, &app->basis, 
    &app->local, &s->local, 
    is_zero_flux, &vlasov_pkpm_inp, &euler_pkpm_inp, app->use_gpu);

  // array for storing F_0^2 (0th Laguerre coefficient) in each cell
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  if (app->use_gpu) {
    s->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[9]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and F_0^2
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 9);
  s->is_first_integ_L2_write_call = true;
  s->is_first_integ_write_call = true;

  // Initialize applied acceleration for use in force update. 
  s->app_accel = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(s->app_accel, 0.0);
  s->has_app_accel = false;
  s->app_accel_evolve = false;
  // setup applied acceleration
  if (s->info.app_accel) {
    s->has_app_accel = true;
    if (s->info.app_accel_evolve) {
      s->app_accel_evolve = s->info.app_accel_evolve;
    }

    s->app_accel_host = s->app_accel;
    if (app->use_gpu) {
      s->app_accel_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
    }
    s->app_accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      3, s->info.app_accel, s->info.app_accel_ctx);
  }
  
  // determine collision type to use in PKPM update
  s->collision_id = s->info.collisions.collision_id;
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_init(app, s, &s->lbo);
  }

  // initialize diffusion if present
  s->has_diffusion = false;  
  if (s->info.diffusion.D) {
    s->has_diffusion = true;
    s->info.diffusion.order = s->info.diffusion.order<2? 2 : s->info.diffusion.order;
    int num_eqn = 3;

    int szD = cdim;
    s->diffD = mkarr(app->use_gpu, szD, 1);
    struct gkyl_array *diffD_host = s->diffD;
    if (app->use_gpu)
      diffD_host = mkarr(false, szD, 1);
    // Set diffusion coefficient in each direction to input value.
    gkyl_array_clear(diffD_host, 0.);
    for (int d=0; d<cdim; d++) gkyl_array_shiftc(diffD_host, s->info.diffusion.D, d);

    if (app->use_gpu) {// note: diffD_host is same as diffD when not on GPUs
      gkyl_array_copy(s->diffD, diffD_host);
      gkyl_array_release(diffD_host);
    }

    bool is_zero_flux[GKYL_MAX_CDIM] = {false};
    s->diff_slvr = gkyl_dg_updater_diffusion_fluid_new(&app->grid, &app->confBasis,
      true, num_eqn, NULL, s->info.diffusion.order, &app->local, is_zero_flux, app->use_gpu);
  }

  // determine which directions are not periodic and get non-periodic BC info
  // also create local skin-ghost ranges for ease of handling non-periodic BCs
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<cdim; ++dir) {
    s->lower_bc[dir] = s->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = s->info.bcx;
      else if (dir == 1)
        bc = s->info.bcy;
      else
        bc = s->info.bcz;

      s->lower_bc[dir] = bc[0];
      s->upper_bc[dir] = bc[1];
    }
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin_dist[dir], &s->lower_ghost_dist[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin_dist[dir], &s->upper_ghost_dist[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
  }

  // allocate buffer for applying BCs
  long buff_sz_dist = 0;
  long buff_sz_fluid = 0;
  // compute buffer size needed
  for (int dir=0; dir<cdim; ++dir) {
    long vol_dist = GKYL_MAX2(s->lower_skin_dist[dir].volume, s->upper_skin_dist[dir].volume);
    buff_sz_dist = buff_sz_dist > vol_dist ? buff_sz_dist : vol_dist;
    long vol_fluid = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz_fluid = buff_sz_fluid > vol_fluid ? buff_sz_fluid : vol_fluid;
  }

  s->bc_buffer_dist = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz_dist);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed_dist = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz_dist);
  s->bc_buffer_up_fixed_dist = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz_dist);

  s->bc_buffer_fluid = mkarr(app->use_gpu, 3*app->confBasis.num_basis, buff_sz_fluid);
  // buffer arrays for fixed function boundary conditions on momentum
  s->bc_buffer_lo_fixed_fluid = mkarr(app->use_gpu, 3*app->confBasis.num_basis, buff_sz_fluid);
  s->bc_buffer_up_fixed_fluid = mkarr(app->use_gpu, 3*app->confBasis.num_basis, buff_sz_fluid);

  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  s->bc_is_absorb = false;
  for (int d=0; d<cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype_dist = GKYL_BC_COPY;
    enum gkyl_bc_basic_type bctype_fluid = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_COPY;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype_dist = GKYL_BC_ABSORB;
      bctype_fluid = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype_dist = GKYL_BC_PKPM_SPECIES_REFLECT;
      bctype_fluid = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype_dist = GKYL_BC_FIXED_FUNC;
      bctype_fluid = GKYL_BC_FIXED_FUNC;
    }

    // Distribution function non-periodic lower boundary conditions
    s->bc_lo_dist[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype_dist, app->basis_on_dev.basis,
      &s->lower_skin_dist[d], &s->lower_ghost_dist[d], s->f->ncomp, app->cdim, app->use_gpu);
    // Momentum non-periodic lower boundary conditions
    s->bc_lo_fluid[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype_fluid, app->basis_on_dev.confBasis,
      &app->lower_skin[d], &app->lower_ghost[d], s->fluid->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_COPY;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype_dist = GKYL_BC_ABSORB;
      bctype_fluid = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype_dist = GKYL_BC_PKPM_SPECIES_REFLECT;
      bctype_fluid = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype_dist = GKYL_BC_FIXED_FUNC;
      bctype_fluid = GKYL_BC_FIXED_FUNC;
    }

    // Distribution function non-periodic upper boundary conditions
    s->bc_up_dist[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype_dist, app->basis_on_dev.basis,
      &s->upper_skin_dist[d], &s->upper_ghost_dist[d], s->f->ncomp, app->cdim, app->use_gpu);
    // Momentum non-periodic upper boundary conditions
    s->bc_up_fluid[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype_fluid, app->basis_on_dev.confBasis,
      &app->upper_skin[d], &app->upper_ghost[d], s->fluid->ncomp, app->cdim, app->use_gpu);
  }
}

void
pkpm_species_apply_ic(gkyl_pkpm_app *app, struct pkpm_species *species, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj_dist;
  proj_dist = gkyl_proj_on_basis_new(&species->grid, &app->basis,
    8, 2, species->info.init_dist, species->info.ctx_dist);
  gkyl_proj_on_basis *proj_fluid;
  proj_fluid = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 3, species->info.init_fluid, species->info.ctx_fluid);

  // run updaters; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as primitive variable computations
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  gkyl_proj_on_basis_advance(proj_dist, t0, &species->local_ext, species->f_host);
  gkyl_proj_on_basis_release(proj_dist);    
  gkyl_proj_on_basis_advance(proj_fluid, t0, &app->local_ext, species->fluid_host);
  gkyl_proj_on_basis_release(proj_fluid);    

  // note: f_host and fluid_host are the same as f and fluid respectively when not on GPUs
  if (app->use_gpu) {
    gkyl_array_copy(species->f, species->f_host);
    gkyl_array_copy(species->fluid, species->fluid_host);
  } 

  // we are pre-computing acceleration for now in case is time-independent
  pkpm_species_calc_app_accel(app, species, t0);

  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  gkyl_bc_basic_buffer_fixed_func(species->bc_lo_dist[0], species->bc_buffer_lo_fixed_dist, species->f);
  gkyl_bc_basic_buffer_fixed_func(species->bc_up_dist[0], species->bc_buffer_up_fixed_dist, species->f);
  gkyl_bc_basic_buffer_fixed_func(species->bc_lo_fluid[0], species->bc_buffer_lo_fixed_fluid, species->fluid);
  gkyl_bc_basic_buffer_fixed_func(species->bc_up_fluid[0], species->bc_buffer_up_fixed_fluid, species->fluid);
}

void
pkpm_species_calc_app_accel(gkyl_pkpm_app *app, struct pkpm_species *species, double tm)
{
  if (species->has_app_accel) {
    gkyl_proj_on_basis_advance(species->app_accel_proj, tm, &app->local_ext, species->app_accel_host);
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(species->app_accel, species->app_accel_host); 
    }
  }
}

void
pkpm_species_calc_pkpm_vars(gkyl_pkpm_app *app, struct pkpm_species *species, 
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(species->pkpm_div_ppar, 0.0);
  gkyl_array_clear(species->pkpm_p_ij, 0.0);
  gkyl_array_clear(species->pkpm_prim, 0.0);
  gkyl_array_clear(species->pkpm_prim_surf, 0.0);

  // Compute rho, p_par, & p_perp
  pkpm_species_moment_calc(&species->pkpm_moms, species->local_ext,
    app->local_ext, fin);

  // Compute div(p_par b_hat) for consistent pressure force in vpar acceleration
  gkyl_dg_calc_pkpm_dist_vars_div_ppar(species->calc_pkpm_dist_vars, 
    &app->local, &species->local, 
    app->field->bvar_surf, app->field->bvar, fin, 
    app->field->max_b, species->pkpm_div_ppar);
  gkyl_array_scale(species->pkpm_div_ppar, species->info.mass);

  // Compute p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in the volume
  gkyl_dg_calc_pkpm_vars_pressure(species->calc_pkpm_vars, &app->local_ext, 
    app->field->bvar, species->pkpm_moms.marr, species->pkpm_p_ij);

  // Compute primitive variables in both the volume and on surfaces
  if (species->bc_is_absorb) {
    gkyl_dg_calc_pkpm_vars_advance(species->calc_pkpm_vars,
      species->pkpm_moms.marr, fluidin, 
      species->pkpm_p_ij, species->pkpm_div_ppar, 
      species->cell_avg_prim, species->pkpm_prim, species->pkpm_prim_surf); 
  }
  else {
    gkyl_dg_calc_pkpm_vars_advance(species->calc_pkpm_vars_ext,
      species->pkpm_moms.marr, fluidin, 
      species->pkpm_p_ij, species->pkpm_div_ppar, 
      species->cell_avg_prim, species->pkpm_prim, species->pkpm_prim_surf); 
  }

  // Compute the penalization terms, both the maximum speed used in Lax fluxes
  // and the total momentum penalization (either 10 moment Roe solve or Lax penalization)
  gkyl_array_clear(species->pkpm_lax, 0.0);
  gkyl_array_clear(species->pkpm_penalization, 0.0);
  gkyl_dg_calc_pkpm_vars_penalization(species->calc_pkpm_vars, 
    &app->local, &app->local_ext, 
    species->pkpm_moms.marr, species->pkpm_p_ij, species->pkpm_prim, fluidin, 
    species->pkpm_lax, species->pkpm_penalization);

  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_species_calc_pkpm_update_vars(gkyl_pkpm_app *app, struct pkpm_species *species, 
  const struct gkyl_array *fin)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(species->pkpm_accel, 0.0); // Incremented in each dimension, so clear beforehand
  gkyl_dg_calc_pkpm_vars_accel(species->calc_pkpm_vars, &app->local, 
    species->pkpm_prim_surf, species->pkpm_prim, 
    app->field->bvar, app->field->div_b, species->lbo.nu_sum, 
    species->pkpm_accel); 

  // Calculate distrbution functions for coupling different Laguerre moments
  gkyl_dg_calc_pkpm_dist_vars_mirror_force(species->calc_pkpm_dist_vars, 
    &app->local, &species->local, 
    species->pkpm_prim, species->lbo.nu_prim_moms, 
    app->field->div_b, species->pkpm_accel, 
    fin, species->F_k_p_1, 
    species->g_dist_source, species->F_k_m_1);

  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_fluid_species_limiter(gkyl_pkpm_app *app, struct pkpm_species *species,
  struct gkyl_array *fin, struct gkyl_array *fluid)
{
  if (species->limit_fluid) {  
    struct timespec tm = gkyl_wall_clock();

    // Compute the PKPM moments from the kinetic equation 
    pkpm_species_moment_calc(&species->pkpm_moms, species->local, app->local, fin);
    // Compute the flow velocity 
    gkyl_dg_calc_pkpm_vars_u(species->calc_pkpm_vars,
      species->pkpm_moms.marr, fluid, 
      species->cell_avg_prim, species->pkpm_u);

    // Limit the slopes of the solution of the fluid system
    gkyl_dg_calc_pkpm_vars_limiter(species->calc_pkpm_vars, &app->local, species->pkpm_u, 
      species->pkpm_moms.marr, species->pkpm_p_ij, fluid);

    app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);

    // Apply boundary conditions after limiting solution
    pkpm_fluid_species_apply_bc(app, species, fluid);
  }
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
pkpm_species_rhs(gkyl_pkpm_app *app, struct pkpm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *em, 
  struct gkyl_array *rhs_f, struct gkyl_array *rhs_fluid)
{
  gkyl_array_clear(species->cflrate_f, 0.0);
  gkyl_array_clear(species->cflrate_fluid, 0.0);
  gkyl_array_clear(rhs_f, 0.0);
  gkyl_array_clear(rhs_fluid, 0.0);

  gkyl_dg_updater_pkpm_advance(species->slvr, 
    &species->local, &app->local, 
    fin, fluidin, 
    species->cflrate_f, species->cflrate_fluid, 
    rhs_f, rhs_fluid);

  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_rhs(app, species, &species->lbo, fin, rhs_f);
  }

  if (species->has_diffusion) {
    gkyl_dg_updater_diffusion_fluid_advance(species->diff_slvr,
      &app->local, species->diffD, fluidin, species->cflrate_fluid, rhs_fluid);
  }

  // If PKPM update is fully explicit, include the fluid-EM coupling in fluid update
  if (app->use_explicit_source) { 
    double qbym = species->info.charge/species->info.mass;
    gkyl_array_set(species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (species->has_app_accel) {
      gkyl_array_accumulate_range(species->qmem, 1.0, species->app_accel, &app->local);
    }
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(species->qmem, qbym, app->field->ext_em, &app->local);  
    }

    gkyl_dg_calc_pkpm_vars_source(species->calc_pkpm_vars, &app->local, 
      species->qmem, species->pkpm_moms.marr, fluidin, rhs_fluid);        
  }

  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr_dist, species->cflrate_f, GKYL_MAX, &species->local);
  gkyl_array_reduce_range(species->omegaCfl_ptr_fluid, species->cflrate_fluid, GKYL_MAX, &app->local);

  double omegaCfl_ho_dist[1];
  double omegaCfl_ho_fluid[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omegaCfl_ho_dist, species->omegaCfl_ptr_dist, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(omegaCfl_ho_fluid, species->omegaCfl_ptr_fluid, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omegaCfl_ho_dist[0] = species->omegaCfl_ptr_dist[0];
    omegaCfl_ho_fluid[0] = species->omegaCfl_ptr_fluid[0];
  }
  double omegaCfl = fmax(omegaCfl_ho_dist[0], omegaCfl_ho_fluid[0]);

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);

  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_species_apply_bc(gkyl_pkpm_app *app, const struct pkpm_species *species, 
  struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo_dist[d], species->bc_buffer_dist, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo_dist[d], species->bc_buffer_lo_fixed_dist, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up_dist[d], species->bc_buffer_dist, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up_dist[d], species->bc_buffer_up_fixed_dist, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_fluid_species_apply_bc(gkyl_pkpm_app *app, const struct pkpm_species *species, 
  struct gkyl_array *fluid)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, app->periodic_dirs, fluid); 

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo_fluid[d], species->bc_buffer_fluid, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo_fluid[d], species->bc_buffer_lo_fixed_fluid, fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up_fluid[d], species->bc_buffer_fluid, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up_fluid[d], species->bc_buffer_up_fixed_fluid, fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, fluid);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_calc_L2(gkyl_pkpm_app *app, double tm, const struct pkpm_species *species)
{
  gkyl_dg_calc_l2_range(&app->basis, 0, species->L2_f, 0, species->f, species->local);
  gkyl_array_scale_range(species->L2_f, species->grid.cellVolume, &species->local);
  
  double L2[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(species->red_L2_f, species->L2_f, GKYL_SUM, &species->local);
    gkyl_cu_memcpy(L2, species->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(L2, species->L2_f, GKYL_SUM, &species->local);
  }
  double L2_global[1] = { 0.0 };
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);
  
  gkyl_dynvec_append(species->integ_L2_f, tm, L2_global);  
}

void
pkpm_species_coll_tm(gkyl_pkpm_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_pkpm_tm tm =
        gkyl_dg_updater_lbo_pkpm_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
pkpm_species_tm(gkyl_pkpm_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  app->stat.fluid_species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_pkpm_tm tm =
      gkyl_dg_updater_pkpm_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.vlasov_tm;
    app->stat.fluid_species_rhs_tm += tm.fluid_tm;
  }
}

// release resources for PKPM species
void
pkpm_species_release(const gkyl_pkpm_app* app, const struct pkpm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->fluid1);
  gkyl_array_release(s->fluidnew);
  gkyl_array_release(s->fluid_dup);

  gkyl_array_release(s->cflrate_f);
  gkyl_array_release(s->cflrate_fluid);

  gkyl_wv_eqn_release(s->equation);

  gkyl_array_release(s->bc_buffer_dist);
  gkyl_array_release(s->bc_buffer_lo_fixed_dist);
  gkyl_array_release(s->bc_buffer_up_fixed_dist);
  gkyl_array_release(s->bc_buffer_fluid);
  gkyl_array_release(s->bc_buffer_lo_fixed_fluid);
  gkyl_array_release(s->bc_buffer_up_fixed_fluid);

  gkyl_comm_release(s->comm);

  if (app->use_gpu) {
    gkyl_array_release(s->f_host);
    gkyl_array_release(s->fluid_host);
  }

  gkyl_array_release(s->qmem); 

  // release moment data
  pkpm_species_moment_release(app, &s->pkpm_moms);
  pkpm_species_moment_release(app, &s->pkpm_moms_diag);

  // release pkpm arrays data and updaters
  gkyl_array_release(s->g_dist_source);
  gkyl_array_release(s->F_k_m_1);
  gkyl_array_release(s->F_k_p_1);
  gkyl_array_release(s->pkpm_div_ppar);
  gkyl_array_release(s->pkpm_prim);
  gkyl_array_release(s->pkpm_u);
  gkyl_array_release(s->pkpm_p_ij);
  gkyl_array_release(s->cell_avg_prim);
  gkyl_array_release(s->pkpm_prim_surf);
  gkyl_array_release(s->pkpm_lax);
  gkyl_array_release(s->pkpm_penalization);
  gkyl_array_release(s->pkpm_accel);
  gkyl_array_release(s->integ_pkpm_mom);
  gkyl_array_release(s->fluid_io);
  gkyl_array_release(s->pkpm_vars_io);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars_ext);
  gkyl_dg_calc_pkpm_dist_vars_release(s->calc_pkpm_dist_vars);

  gkyl_dg_updater_pkpm_release(s->slvr);

  gkyl_array_release(s->L2_f);
  gkyl_dynvec_release(s->integ_L2_f);
  gkyl_dynvec_release(s->integ_diag);
  
  gkyl_array_release(s->app_accel);
  if (s->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(s->app_accel_host);
    }
    gkyl_proj_on_basis_release(s->app_accel_proj);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_release(app, &s->lbo);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(s->bc_lo_dist[d]);
    gkyl_bc_basic_release(s->bc_up_dist[d]);
    gkyl_bc_basic_release(s->bc_lo_fluid[d]);
    gkyl_bc_basic_release(s->bc_up_fluid[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr_dist);
    gkyl_cu_free(s->omegaCfl_ptr_fluid);
    gkyl_cu_free(s->red_L2_f);
    gkyl_cu_free(s->red_integ_diag);
  }
  else {
    gkyl_free(s->omegaCfl_ptr_dist);
    gkyl_free(s->omegaCfl_ptr_fluid);
  }
}
