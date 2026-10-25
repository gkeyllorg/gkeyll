#include <gkyl_app_priv.h>
#include <gkyl_flash_priv.h>
#include <gkyl_fv_proj.h>
#include <gkyl_wv_euler.h>

void flash_fluid_species_init(struct gkyl_flash_inp *inp, 
  struct gkyl_flash_app *app, struct flash_fluid_species *s)
{
  int cdim = app->cdim;

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) ghost[d] = app->nghost_fluid[d];

  // Set up finite volume grid.
  gkyl_rect_grid_init(&s->grid, cdim, app->grid_fluid.lower, app->grid_fluid.upper, app->grid_fluid.cells);

  // Local and global ranges over finite volume grid.
  gkyl_create_grid_ranges(&s->grid, ghost, &s->global_ext, &s->global);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  
  // Fluid variables, [rho, rho*ux, rho*uy, rho*uz, E]
  s->fluid = mkarr(app->use_gpu, 5, app->local_fluid_ext.volume);
  s->fluid1 = mkarr(app->use_gpu, 5, app->local_fluid_ext.volume);
  s->fluidnew = mkarr(app->use_gpu, 5, app->local_fluid_ext.volume);

  // Primitive variables, [rho, ux, uy, uz, p]
  s->fluid_prim = mkarr(app->use_gpu, 5, app->local_fluid_ext.volume);

  s->fluid_phi = mkarr(app->use_gpu, 1, app->local_fluid_ext.volume);

  // Flux contributions from the kinetic equation, [div(P)-grad(p), (P - pI) : grad(u)]
  s->flash_coupling_fluxes = mkarr(app->use_gpu, 3, app->local_fluid_ext.volume);

  // Euler wave equation object
  s->gas_gamma = s->info.gas_gamma == 0 ? 5.0/3.0 : s->info.gas_gamma;
  s->wv_eqn = gkyl_wv_euler_new(s->gas_gamma, app->use_gpu);

  // Geometry object for Euler wave_prop solver (on fluid grid)
  s->wv_geom = gkyl_wave_geom_new(&app->grid_fluid, &app->local_fluid_ext, 0, 0, false);

  // Flux limiter to use. Defaults to Monotized Centered
  enum gkyl_wave_limiter limiter = s->info.limiter == 0 ? GKYL_MONOTONIZED_CENTERED : s->info.limiter;

  // Euler wave propagation updater input object
  struct gkyl_wave_prop_inp wave_prop_inp = {
    .grid = &app->grid_fluid,
    .equation = s->wv_eqn,
    .limiter = limiter,
    .num_up_dirs = cdim,
    .update_dirs = { 0 }, // HARDCODE: 1x only
    .cfl = app->cfl,
    .geom = s->wv_geom,
    .comm = app->comm_fluid,
    .force_low_order_flux = true
  };

  // Euler equations wave propagation updater
  s->up_wave_prop = gkyl_wave_prop_new(&wave_prop_inp);

  // Boundary conditions
  // Determine periodic boundaries
  int num_periodic_dir = app->num_periodic_dir;
  int periodic_dirs[3];
  for (int d=0; d<cdim; ++d) 
    periodic_dirs[d] = app->periodic_dirs[d];
  
  // Store non-periodic boundary directions
  int is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    is_np[periodic_dirs[d]] = 0;
  }

  // Get non-periodic BC types from the input file
  for (int dir = 0; dir<cdim; ++dir) {
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
  }

  // Allocate buffer for applying BCs
  long buff_sz_dist = 0;
  long buff_sz_fluid = 0;
  // compute buffer size needed
  for (int dir=0; dir<cdim; ++dir) {
    long vol_fluid = GKYL_MAX2(app->lower_skin_fluid[dir].volume, app->upper_skin_fluid[dir].volume);
    buff_sz_fluid = buff_sz_fluid > vol_fluid ? buff_sz_fluid : vol_fluid;
  }

  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  s->bc_is_absorb = false;
  for (int d=0; d<cdim; ++d) {
    s->bc_lo[d] = s->bc_up[d] = 0;

    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype_dist = GKYL_BC_COPY;
    enum gkyl_bc_basic_type bctype_fluid = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      s->bc_lo[d] = gkyl_wv_apply_bc_new(&s->grid, s->wv_eqn, s->wv_geom, d, 
        GKYL_LOWER_EDGE, ghost, bc_copy, 0);
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      // bctype_dist = GKYL_BC_ABSORB;
      // bctype_fluid = GKYL_BC_ABSORB;
      // s->bc_is_absorb = true;
      gkyl_exit("Absorb BCs not supported.\n");
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype_dist = GKYL_BC_PKPM_SPECIES_REFLECT;
      s->bc_lo[d] = gkyl_wv_apply_bc_new(&s->grid, s->wv_eqn, s->wv_geom, d, 
        GKYL_LOWER_EDGE, ghost, s->wv_eqn->wall_bc_func, 0);
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      // bctype_dist = GKYL_BC_FIXED_FUNC;
      // bctype_fluid = GKYL_BC_FIXED_FUNC;
      gkyl_exit("Fixed function BCs not supported.\n");
    }

    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d] == GKYL_SPECIES_COPY) {
      s->bc_up[d] = gkyl_wv_apply_bc_new(&s->grid, s->wv_eqn, s->wv_geom, d, 
        GKYL_UPPER_EDGE, ghost, bc_copy, 0);
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      // bctype_fluid = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
      gkyl_exit("Absorb BCs not supported.\n");
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      s->bc_up[d] = gkyl_wv_apply_bc_new(&s->grid, s->wv_eqn, s->wv_geom, d, 
        GKYL_UPPER_EDGE, ghost, s->wv_eqn->wall_bc_func, 0);
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      // bctype_fluid = GKYL_BC_FIXED_FUNC;
      gkyl_exit("Fixed function BCs not supported.\n");
    }
  }

  // Skin and ghost ranges used in applying PKPM pressure fluxes at the boundaries
  for (int d=0; d<cdim; ++d) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin[d], &s->lower_ghost[d], d, GKYL_LOWER_EDGE, &app->local_fluid_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin[d], &s->upper_ghost[d], d, GKYL_UPPER_EDGE, &app->local_fluid_ext, ghost);
  }
}

// Apply initial conditions for finite volume Euler equation moments.
void
flash_fluid_species_apply_ic(gkyl_flash_app *app, struct flash_fluid_species *fs, double t0) {
  int num_quad = 2;
  int num_ret_vals = 5; // Number of moments for Euler equations

  gkyl_fv_proj *proj_fluid;
  proj_fluid = gkyl_fv_proj_new(&app->grid_fluid, num_quad, num_ret_vals, 
    fs->info.init_fluid, fs->info.ctx);

  gkyl_fv_proj_advance(proj_fluid, t0, &app->local_fluid_ext, fs->fluid);
  gkyl_fv_proj_release(proj_fluid);
}

// Determine periodic directions and apply boundary conditions to fluid variables
void
flash_fluid_species_apply_bc(gkyl_flash_app *app, const struct flash_fluid_species *species, 
  struct gkyl_array *fluid)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  gkyl_comm_array_per_sync(app->comm_fluid, &app->local_fluid, &app->local_fluid_ext, num_periodic_dir,
   app->periodic_dirs, fluid);

  double tm = 0.0; // HARDCODE: No time-dependent boundary conditions.

  // Apply boundary conditions in non-periodic directions with definied BCs.
  for (int d=0; d<app->cdim; ++d) {
    if (is_np_bc[d]) {
      if (species->bc_lo[d])
        gkyl_wv_apply_bc_advance(species->bc_lo[d], tm, &species->local_ext, fluid);
      if (species->bc_up[d])
        gkyl_wv_apply_bc_advance(species->bc_up[d], tm, &species->local_ext, fluid);
    }
  }

  gkyl_comm_array_sync(app->comm_fluid, &app->local_fluid, &app->local_fluid_ext, fluid);

  app->stat.fluid_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Compute (rho, ux, uy, uz, p) for diagnostic writing.
void
flash_fluid_species_prim_vars(gkyl_flash_app *app, 
  struct flash_fluid_species *s, const struct gkyl_array *fluid, struct gkyl_array *fluid_prim)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_fluid);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&app->local_fluid, iter.idx);

    const double *fluid_d = (const double *)gkyl_array_cfetch(fluid, linc);
    double *prim_d = (double *)gkyl_array_fetch(fluid_prim, linc);

    double rho = fluid_d[0];
    double u = fluid_d[1]/rho;
    double v = fluid_d[2]/rho;
    double w = fluid_d[3]/rho;
    double q_sq = u*u + v*v + w*w;

    prim_d[0] = rho;
    prim_d[1] = u;
    prim_d[2] = v;
    prim_d[3] = w;
    prim_d[4] = (fluid_d[4] - 0.5*rho*q_sq)*(s->gas_gamma - 1.0);
  }
}

// Advance Euler equations in time by step dt0 with wave propagation scheme.
// Return suggested timestep.
double
flash_fluid_species_update(gkyl_flash_app *app, double tcurr, double dt0,
  struct flash_species *s, struct flash_fluid_species *fs, 
  const struct gkyl_array *fluidin, struct gkyl_array *fluidout)
{
  struct gkyl_wave_prop_status wv_status;
  struct timespec wst = gkyl_wall_clock();

  // Compute fluid update to next time step with wave propagation updater
  wv_status = gkyl_wave_prop_advance(fs->up_wave_prop, tcurr, dt0,
      &app->local_fluid, fs->fluid_phi, fluidin, fluidout);

  app->stat.fluid_rhs_tm += gkyl_time_diff_now_sec(wst);
  wst = gkyl_wall_clock();

  // Interpolate kinetic flux contributions onto fluid grid
  gkyl_flash_interp_vars_coupling_fluxes_to_fv(app->interp_vars,
    &app->local_kinetic, &app->local_fluid,
    s->flash_coupling_fluxes, fs->flash_coupling_fluxes);

  app->stat.coupling_tm += gkyl_time_diff_now_sec(wst);
  wst = gkyl_wall_clock();

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_fluid);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&app->local_fluid, iter.idx);

    const double *flux_d = (const double *)gkyl_array_cfetch(fs->flash_coupling_fluxes, linc);
    const double *fluidin_d = (const double *)gkyl_array_cfetch(fluidin, linc);
    double *fluidout_d = (double *)gkyl_array_fetch(fluidout, linc);

    // Momentum flux contribution, div(Pi)
    double div_pi = flux_d[0];

    // Energy flux contribution, div(Pi.u) + div(q) = Pi:grad(u) + u.div(Pi) + div(q)
    double pi_gradu = flux_d[1];
    double u_divpi = fluidin_d[1]/fluidin_d[0]*div_pi;
    double div_q = flux_d[2];

    fluidout_d[1] -= dt0*div_pi;
    fluidout_d[4] -= dt0*(pi_gradu + u_divpi + div_q);
  }

  app->stat.fluid_rhs_tm += gkyl_time_diff_now_sec(wst);

  return wv_status.dt_suggested;
}

void
flash_fluid_species_release(struct flash_fluid_species *fs)
{
  gkyl_array_release(fs->fluid);
  gkyl_array_release(fs->fluid1);
  gkyl_array_release(fs->fluidnew);
  gkyl_array_release(fs->fluid_prim);
  gkyl_array_release(fs->fluid_phi);
  gkyl_array_release(fs->flash_coupling_fluxes);

  gkyl_wave_geom_release(fs->wv_geom);
  gkyl_wv_eqn_release(fs->wv_eqn);
  gkyl_wave_prop_release(fs->up_wave_prop);

  // BC updaters only exist in the cdim configuration-space directions.
  for (int d=0; d<fs->grid.ndim; ++d) {
    if (fs->bc_lo[d]) 
      gkyl_wv_apply_bc_release(fs->bc_lo[d]);
    if (fs->bc_up[d]) 
      gkyl_wv_apply_bc_release(fs->bc_up[d]);
  }
}