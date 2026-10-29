#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_flash_priv.h>
#include <gkyl_flash_proj_pkpm_lte_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_util.h>

void flash_species_init(struct gkyl_flash_inp *inp, 
  struct gkyl_flash_app *app, struct flash_species *s) 
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  int pdim = cdim + vdim;

  s->mass = s->info.mass;

  int cells_kinetic[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_kinetic[GKYL_MAX_DIM], upper_kinetic[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    ghost[d] = app->nghost_kinetic[d];
    cells_kinetic[d] = (int)(inp->cells_fluid[d]/2);
    lower_kinetic[d] = inp->lower_conf[d];
    upper_kinetic[d] = inp->upper_conf[d];
  }
  for (int d=cdim; d<pdim; ++d) {
    ghost[d] = 0;
    cells_kinetic[d] = s->info.cells[d-cdim];
    lower_kinetic[d] = s->info.lower[d-cdim];
    upper_kinetic[d] = s->info.upper[d-cdim];

    ghost_vel[d-cdim] = 0;
    cells_vel[d-cdim] = s->info.cells[d-cdim];
    lower_vel[d-cdim] = s->info.lower[d-cdim];
    upper_vel[d-cdim] = s->info.upper[d-cdim];
  }

  // Phase space grid and ranges
  gkyl_rect_grid_init(&s->grid, pdim, lower_kinetic, upper_kinetic, cells_kinetic);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->global_ext, &s->global);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);

  // Velocity space grid and ranges
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  s->comm = gkyl_comm_extend_comm(app->comm_kinetic, &s->local_vel);

  // Distribution function arrays
  // PKPM has two distribution functions F_0 and T_perp*G, first two Laguerre moments
  s->f = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);

  // Storage for computing L2 norm of F_0.
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);

  // Array to store momentum density from fluid solve (rho u_x, rho u_y, rho u_z)
  s->fluid_rhou = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_kinetic_ext.volume);

  // Distribution function arrays for coupling different Laguerre moments
  // g_dist_source has two components: 
  // [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
  // (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0]
  // First component is mirror force source *distribution*, second component is *total* vperp characteristics source.
  s->g_dist_source = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
  s->F_k_p_1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume); 
  s->F_k_m_1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume); 

  // CFL rate per cell for distribution function
  s->cflrate_f = mkarr(app->use_gpu, 1, s->local_ext.volume);

  // div(p_par b_hat), for self-consistent total pressure force
  s->pkpm_div_ppar = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_kinetic_ext.volume);
  // allocate array to store primitive moments : 
  // [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp, 3*T_xx/m, 3*T_yy/m, 3*T_zz/m]
  // pressure p_ij : (p_par - p_perp) b_i b_j + p_perp g_ij
  s->pkpm_prim = mkarr(app->use_gpu, 9*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  s->pkpm_p_ij = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  // boolean array for primitive variables [rho, p_par, p_perp] is negative at control points
  // *only used for diagnostic purposes*
  s->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_kinetic_ext.volume);

  int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
  // Surface primitive variables (2*cdim*4 components). Ordered as:
  // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 3.0*Txx_xl/m, 3.0*Txx_xr/m, 
  //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 
  //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
  s->pkpm_prim_surf = mkarr(app->use_gpu, 2*cdim*4*Nbasis_surf, app->local_kinetic_ext.volume);
  // Surface expansions of:
  // 1. Lax penalization maximum speed lambda_i = |u_i| + sqrt(3*P_ii/rho)
  // 2. Momentum solve penalization (either 10 moment Roe solve or Lax penalization)
  //    Momentum penalization returns the full penalization, e.g., the modal surface expansion of
  //    jump in the fluctuations from the 10 moment Roe solve 1/2 (A^+ Delta Q - A^- Delta Q)
  s->pkpm_lax = mkarr(app->use_gpu, cdim*Nbasis_surf, app->local_kinetic_ext.volume);
  s->pkpm_penalization = mkarr(app->use_gpu, 3*cdim*Nbasis_surf, app->local_kinetic_ext.volume);

  // allocate array for pkpm acceleration variables, stored in pkpm_accel: 
  // 0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
  // 1: bb_grad_u (bb : grad(u))
  // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
  // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2*nu)
  s->pkpm_accel = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_kinetic_ext.volume);

  // Updater to compute LBO collisions.
  flash_species_lbo_init(app, s, &s->lbo);

  // Wave equation object for upwinding fluid equations
  // We use the 10 moment system since PKPM model generates a full pressure tensor
  // k0 = 0.0 and use_grad_closure = false because we do not need a closure 
  // (the kinetic equations are the closure)
  s->eqn_ten_mom = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, app->use_gpu);

  // Geometry object for ten moment equations.
  s->ten_mom_geom = gkyl_wave_geom_new(&app->grid_kinetic, &app->local_kinetic_ext, 0, 0, app->use_gpu);

  // Updater to compute PKPM vars:
  // pressure, primitive variables, and acceleration variables
  // also stores kernels for computing source terms, integrated variables
  double limiter_fac = 1.0/sqrt(3); // default 1/sqrt(3)
  s->calc_pkpm_vars_ext = gkyl_dg_calc_pkpm_vars_new(&app->grid_kinetic,
    &app->confBasis, &app->local_kinetic_ext, s->eqn_ten_mom, s->ten_mom_geom, limiter_fac, app->use_gpu);
  s->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_new(&app->grid_kinetic,
    &app->confBasis, &app->local_kinetic, s->eqn_ten_mom, s->ten_mom_geom, limiter_fac, app->use_gpu);

  // Updater to compute PKPM dist vars
  // div(p_par b_hat) for self-consistent total pressure force and distribution function sources for
  // Laguerre couplings and vperp characteristics 
  s->calc_pkpm_dist_vars = gkyl_dg_calc_pkpm_dist_vars_new(&s->grid, &app->confBasis, app->use_gpu);

  // Updater to compute PKPM moments (rho, p_par, p_perp, M1) and diagnostics
  s->pkpm_moms = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  s->pkpm_mom_calc = gkyl_dg_updater_moment_pkpm_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local_kinetic, &s->local_vel, s->info.mass, false, app->use_gpu);

  // Diagnostic moments (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
  s->pkpm_moms_diag = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  s->pkpm_mom_calc_diag = gkyl_dg_updater_moment_pkpm_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local_kinetic, &s->local_vel, s->info.mass, true, app->use_gpu);

  // Updater to compute flux contributions to Euler equations, [div(P)-grad(p), (P - pI) : grad(u), div(q)]
  s->calc_flash_coupling = gkyl_dg_calc_flash_coupling_new(&s->grid, &app->grid_kinetic, &app->confBasis);
  s->flash_coupling_fluxes = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_kinetic_ext.volume);

  // Flag for initializing using input arrays rather than function
  s->init_dist_from_arrays = s->info.init_dist_from_arrays;
  if (s->init_dist_from_arrays) {
    s->init_dist_moms = mkarr(app->use_gpu, 5*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  }

  // Flag to skip initializing the distribution function.
  s->skip_init_dist = s->info.skip_init_dist;

  // Updater for kinetic RHS.
  bool is_zero_flux[GKYL_MAX_DIM] = {false};
  struct gkyl_dg_vlasov_pkpm_auxfields aux = {
    .bvar = app->field->bvar, .bvar_surf = app->field->bvar_surf, 
    .pkpm_prim = s->pkpm_prim, .pkpm_prim_surf = s->pkpm_prim_surf, 
    .max_b = app->field->max_b, .pkpm_lax = s->pkpm_lax, 
    .div_b = app->field->div_b, .pkpm_accel_vars = s->pkpm_accel, 
    .g_dist_source = s->g_dist_source};
  s->updater = gkyl_dg_updater_flash_new(&app->grid_kinetic, &s->grid,
    &app->confBasis, &app->basis, &app->local_kinetic, &s->local, is_zero_flux, &aux);

  // Pointer to double used in computing CFL contribution
  s->omegaCFL_ptr = gkyl_malloc(sizeof(double));

  // Boundary conditions
  // Determine periodic boundaries
  int num_periodic_dir = app->num_periodic_dir;
  int periodic_dirs[3];
  for (int d=0; d<cdim; ++d) 
    periodic_dirs[d] = app->periodic_dirs[d];

  int is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    is_np[periodic_dirs[d]] = 0;
  }

  // Get non-periodic BC info from input
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
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin[dir], &s->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin[dir], &s->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);

    // Create local lower skin and ghost ranges for configuration space
    gkyl_skin_ghost_ranges(&s->lower_skin_conf[dir], &s->lower_ghost_conf[dir], dir, GKYL_LOWER_EDGE, &app->local_kinetic_ext, ghost);
    // Create local upper skin and ghost ranges for configuration space
    gkyl_skin_ghost_ranges(&s->upper_skin_conf[dir], &s->upper_ghost_conf[dir], dir, GKYL_UPPER_EDGE, &app->local_kinetic_ext, ghost);
  }

  // Allocate buffer for applying BCs
  long buff_sz_dist = 0;
  long buff_sz_conf = 0;
  // compute buffer size needed
  for (int dir=0; dir<cdim; ++dir) {
    long vol_dist = GKYL_MAX2(s->lower_skin[dir].volume, s->upper_skin[dir].volume);
    buff_sz_dist = buff_sz_dist > vol_dist ? buff_sz_dist : vol_dist;
    long vol_conf = GKYL_MAX2(s->lower_skin_conf[dir].volume, s->upper_skin_conf[dir].volume);
    buff_sz_conf = buff_sz_conf > vol_conf ? buff_sz_conf : vol_conf;
  }

  s->bc_buffer = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz_dist);
  // buffer arrays for fixed function boundary conditions on distribution function and interpolated rhou
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz_dist);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz_dist);

  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  s->bc_is_absorb = false;
  for (int d=0; d<cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY, bctype_conf = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      // bctype_dist = GKYL_BC_ABSORB;
      // bctype_fluid = GKYL_BC_ABSORB;
      // s->bc_is_absorb = true;
      gkyl_exit("Absorb BCs not supported.\n");
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_PKPM_SPECIES_REFLECT;
      bctype_conf = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      // bctype_dist = GKYL_BC_FIXED_FUNC;
      // bctype_fluid = GKYL_BC_FIXED_FUNC;
      gkyl_exit("Fixed function BCs not supported.\n");
    }
    // Distribution function non-periodic lower boundary conditions
    s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, &app->basis,
      &s->lower_skin[d], &s->lower_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    s->bc_conf_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype_conf, &app->confBasis,
      &s->lower_skin_conf[d], &s->lower_ghost_conf[d], s->fluid_rhou->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    bctype = bctype_conf = GKYL_BC_COPY;
    if (s->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
      gkyl_exit("Absorb BCs not supported.\n");
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_PKPM_SPECIES_REFLECT;
      bctype_conf = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype = GKYL_BC_FIXED_FUNC;
      gkyl_exit("Fixed function BCs not supported.\n");
    }
    // Distribution function non-periodic upper boundary conditions
    s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, &app->basis,
      &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    s->bc_conf_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype_conf, &app->confBasis,
      &s->upper_skin_conf[d], &s->upper_ghost_conf[d], s->fluid_rhou->ncomp, app->cdim, app->use_gpu);
  }
}

// Compute initial distribution function as Maxwellian using input fluid moment data.
// Moments are [n, ux, uy, uz, T/m].
void
flash_species_calc_distf_from_moms(gkyl_flash_app *app, struct flash_species *s, struct gkyl_array *init_dist_moms)
{
  struct gkyl_flash_proj_pkpm_lte_on_basis_inp inp = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .use_gpu = app->use_gpu,
    .quad_type = GKYL_GAUSS_QUAD
  };

  struct gkyl_flash_proj_pkpm_lte_on_basis *proj = gkyl_flash_proj_pkpm_lte_on_basis_new(&inp);

  gkyl_flash_proj_pkpm_lte_on_basis_advance(proj, &s->local, &app->local_kinetic, init_dist_moms, s->f);
  gkyl_flash_proj_pkpm_lte_on_basis_release(proj);
}


// Apply initial conditions for distribution function, projecting ICs on p=1 modal basis.
void
flash_species_apply_ic(gkyl_flash_app *app, struct flash_species *s, double t0) 
{
  // Return early if we aren't initializing the distribution function.
  if (s->skip_init_dist) return; 

  // If user specified input moments to initialize distribution functions, use them.
  // Otherwise, project from provided function.
  if (s->init_dist_from_arrays) {
    flash_species_calc_distf_from_moms(app, s, s->init_dist_moms);
    return;
  }

  int poly_order = app->basis.poly_order;
  int num_quad = poly_order+1;

  gkyl_proj_on_basis *proj_dist;
  proj_dist = gkyl_proj_on_basis_new(&s->grid, &app->basis, 
    num_quad, 2, s->info.init_dist, s->info.ctx);

  gkyl_proj_on_basis_advance(proj_dist, t0, &s->local_ext, s->f);
  gkyl_proj_on_basis_release(proj_dist);
}

// Compute pkpm moments, [rho, p_par, p_perp].
void
flash_species_calc_pkpm_moms(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(species->pkpm_moms, 0.0);

  // Compute rho, p_par, and p_perp
  gkyl_dg_updater_moment_pkpm_advance(species->pkpm_mom_calc,
    &species->local_ext, &app->local_kinetic_ext, fin, species->pkpm_moms);

  app->stat.species_mom_tm += gkyl_time_diff_now_sec(wst);
}

// Compute quantities needed by the kinetic update:
// pkpm_div_ppar, pkpm_p_ij, pkpm_prim, pkpm_prim_surf,
// pkpm_lax, pkpm_penalization.
void
flash_species_calc_pkpm_vars(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(species->pkpm_div_ppar, 0.0);
  gkyl_array_clear(species->pkpm_p_ij, 0.0);
  gkyl_array_clear(species->pkpm_prim, 0.0);
  gkyl_array_clear(species->pkpm_prim_surf, 0.0);

  // Compute div(p_par b_hat)
  gkyl_dg_calc_pkpm_dist_vars_div_ppar(species->calc_pkpm_dist_vars,
    &app->local_kinetic, &species->local, 
    app->field->bvar_surf, app->field->bvar, fin, 
    app->field->max_b, species->pkpm_div_ppar);
  gkyl_array_scale(species->pkpm_div_ppar, species->info.mass);

  // Compute p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in the volume
  gkyl_dg_calc_pkpm_vars_pressure(species->calc_pkpm_vars, &app->local_kinetic_ext, 
    app->field->bvar, species->pkpm_moms, species->pkpm_p_ij);

  // Compute primitive variables in both the volume and on surfaces
  gkyl_dg_calc_pkpm_vars_advance(species->calc_pkpm_vars_ext,
    species->pkpm_moms, species->fluid_rhou, 
    species->pkpm_p_ij, species->pkpm_div_ppar,
    species->cell_avg_prim, species->pkpm_prim, species->pkpm_prim_surf);

  // Compute the penalization terms, both the maximum speed used in Lax fluxes
  // and the total momentum penalization (either 10 moment Roe solve or Lax penalization)
  gkyl_array_clear(species->pkpm_lax, 0.0);
  gkyl_array_clear(species->pkpm_penalization, 0.0);
  gkyl_dg_calc_pkpm_vars_penalization(species->calc_pkpm_vars, 
    &app->local_kinetic, &app->local_kinetic_ext, 
    species->pkpm_moms, species->pkpm_p_ij, species->pkpm_prim, species->fluid_rhou,
    species->pkpm_lax, species->pkpm_penalization);

  app->stat.species_mom_tm += gkyl_time_diff_now_sec(wst);
}

// Compute variables used in the kinetic update:
// LBO moments, pkpm_accel, F_k_p_1, g_dist_source, F_k_m_1.
void
flash_species_calc_pkpm_update_vars(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin)
{
  // Compute LBO moments first: nu_prim_moms is needed by the mirror force below.
  // pkpm_moms should be computed prior to this call.
  flash_species_lbo_moms(app, species, &species->lbo, fin);

  struct timespec wst = gkyl_wall_clock();

  // Calculate PKPM acceleration vars.
  gkyl_array_clear(species->pkpm_accel, 0.0);
  gkyl_dg_calc_pkpm_vars_accel(species->calc_pkpm_vars, &app->local_kinetic,
    species->pkpm_prim_surf, species->pkpm_prim, 
    app->field->bvar, app->field->div_b, species->lbo.nu_sum, 
    species->pkpm_accel);

  // Calculate distribution function vars for coupling Laguerre moments through mirror force.
  gkyl_dg_calc_pkpm_dist_vars_mirror_force(species->calc_pkpm_dist_vars,
    &app->local_kinetic, &species->local,
    species->pkpm_prim, species->lbo.nu_prim_moms,
    app->field->div_b, species->pkpm_accel, 
    fin, species->F_k_p_1, 
    species->g_dist_source, species->F_k_m_1);

  app->stat.species_mom_tm += gkyl_time_diff_now_sec(wst);
}

// Compute the RHS for the PKPM Vlasov update. Return maximum
// stable timestep.
double 
flash_species_rhs(gkyl_flash_app *app, 
  struct flash_species *s, const struct gkyl_array *fin,
  struct gkyl_array *rhs_f)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(s->cflrate_f, 0.0);
  gkyl_array_clear(rhs_f, 0.0);

  // Call updater to compute RHS of PKPM kinetic equations.
  gkyl_dg_updater_flash_advance(s->updater, &s->local,
    &app->local_kinetic, fin, s->cflrate_f, rhs_f);

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);

  // Call updater to compute RHS of PKPM LBO collision operator.
  flash_species_lbo_rhs(app, s, &s->lbo, fin, rhs_f);

  // Compute timestep
  gkyl_array_reduce_range(s->omegaCFL_ptr, s->cflrate_f, GKYL_MAX, &s->local);

  return app->cfl/s->omegaCFL_ptr[0];
}

double
flash_species_compute_dt(gkyl_flash_app *app, struct flash_species *s)
{ 
  gkyl_array_clear(s->cflrate_f, 0.0);

  // Compute PKPM moments, [rho, p_par, p_perp].
  flash_species_calc_pkpm_moms(app, s, s->f);

  // Compute coupling moments, self-consistent pressure
  // force (div(p_par b_hat)), and primitive variables
  flash_species_calc_pkpm_vars(app, s, s->f);

  // Compute remaining PKPM variables: acceleration variables in kinetic equation and 
  // source distribution functions for Laguerre couplings.
  flash_species_calc_pkpm_update_vars(app, s, s->f);

  // Vlasov minimum timestep
  gkyl_dg_updater_flash_compute_cflrate(s->updater, &s->grid, &s->local, &app->local_kinetic,
    s->pkpm_accel, s->pkpm_lax, app->field->div_b, app->field->max_b, s->cflrate_f);

  // LBO minimum timestep
  flash_species_lbo_compute_cflrate(app, s, &s->lbo);

  // Compute timestep
  double omegaCFL;
  gkyl_array_reduce_range(&omegaCFL, s->cflrate_f, GKYL_MAX, &s->local);
  // printf("\nGKEYLL: Computed dt = %.6e\n", app->cfl/omegaCFL);

  return app->cfl/omegaCFL;
}

// Apply boundary conditions on distribution function 
void
flash_species_apply_bc(gkyl_flash_app *app, const struct flash_species *species, 
  struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  // Fill ghost cells in periodic directions.
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, app->periodic_dirs, f);

  // For non-periodic directions, advance specified BCs
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {
      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed,f);
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
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
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

// Apply boundary conditions on interpolated momentum density vector 
static void
flash_species_apply_bc_rhou(gkyl_flash_app *app, const struct flash_species *species, 
  struct gkyl_array *fluid_rhou)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  // Fill ghost cells in periodic directions.
  gkyl_comm_array_per_sync(app->comm_kinetic, &app->local_kinetic, &app->local_kinetic_ext,
    num_periodic_dir, app->periodic_dirs, fluid_rhou);

  // For non-periodic directions, advance specified BCs
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {
      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_conf_lo[d], species->bc_buffer, fluid_rhou);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_conf_lo[d], species->bc_buffer_lo_fixed, fluid_rhou);
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
          gkyl_bc_basic_advance(species->bc_conf_up[d], species->bc_buffer, fluid_rhou);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_conf_up[d], species->bc_buffer_up_fixed, fluid_rhou);
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

  gkyl_comm_array_sync(app->comm_kinetic, &app->local_kinetic, &app->local_kinetic_ext, fluid_rhou);
  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Interpolate cell-average finite volume Euler fluid variables to p=1 modal
// DG on coarsened mesh for use in PKPM kinetic update.
void
flash_species_interp_rhou_to_dg(struct gkyl_flash_app *app,
  struct flash_species *species, const struct gkyl_array *fv_fluidin)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_flash_interp_vars_rhou_to_dg(app->interp_vars,
    &app->local_kinetic, &app->local_fluid, fv_fluidin, species->fluid_rhou);

  for (int d=0; d<3; ++d)
    gkyl_dg_mul_op(&app->confBasis, d, species->fluid_rhou, d, species->fluid_rhou, 0, species->pkpm_moms);

  app->stat.coupling_tm += gkyl_time_diff_now_sec(wst);

  // Apply boundary conditions to rho*u
  flash_species_apply_bc_rhou(app, species, species->fluid_rhou);
}

// Compute flux contributions to Euler equations from kinetic closure.
void
flash_species_calc_coupling_fluxes(struct gkyl_flash_app *app,
  struct flash_species *species, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(species->flash_coupling_fluxes, 0.0);

  // Compute div(P) - grad(p) and (P-pI):grad(u)
  gkyl_dg_calc_flash_coupling_advance(species->calc_flash_coupling, 
    &species->local, &app->local_kinetic, species->mass, fin,
    species->pkpm_moms, species->pkpm_prim, species->pkpm_prim_surf, 
    species->pkpm_div_ppar, app->field->bvar, 
    app->field->bvar_surf, app->field->max_b,
    species->flash_coupling_fluxes);

  app->stat.coupling_tm += gkyl_time_diff_now_sec(wst);
}

void
flash_species_set_vel_space_extents(gkyl_flash_app *app, struct flash_species *s, double lower_vel[GKYL_MAX_VDIM], double upper_vel[GKYL_MAX_VDIM])
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  int pdim = cdim+vdim;
  double lower[GKYL_MAX_DIM];
  double upper[GKYL_MAX_DIM];

  int cells[GKYL_MAX_DIM], cells_vel[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    cells[d] = s->grid.cells[d];
    lower[d] = s->grid.lower[d];
    upper[d] = s->grid.upper[d];
  }
  for (int d=cdim; d<pdim; ++d) {
    cells[d] = s->grid.cells[d];
    lower[d] = lower_vel[d-cdim];
    upper[d] = upper_vel[d-cdim];

    cells_vel[d-cdim] = s->grid_vel.cells[d-cdim];
  }

  for (int d=0; d<s->grid.ndim; ++d) cells[d] = s->grid.cells[d];
  for (int d=0; d<s->grid_vel.ndim; ++d) cells_vel[d] = s->grid_vel.cells[d];

  // Update extents of phase grid
  gkyl_rect_grid_init(&s->grid, s->grid.ndim, lower, upper, cells);

  // Update extents of velocity grid
  gkyl_rect_grid_init(&s->grid_vel, s->grid_vel.ndim, lower_vel, upper_vel, cells_vel);

  // printf("Set extents: %.8e\n", *(double *)(s->info.collisions.ctx+0x90));
  // Re-initialize updaters?
  // Updater to compute LBO collisions.
  flash_species_lbo_release(&s->lbo);
  flash_species_lbo_init(app, s, &s->lbo);

  // Wave equation object for upwinding fluid equations
  // We use the 10 moment system since PKPM model generates a full pressure tensor
  // k0 = 0.0 and use_grad_closure = false because we do not need a closure 
  // (the kinetic equations are the closure)
  gkyl_wv_eqn_release(s->eqn_ten_mom);
  s->eqn_ten_mom = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, app->use_gpu);

  // Geometry object for ten moment equations.
  gkyl_wave_geom_release(s->ten_mom_geom);
  s->ten_mom_geom = gkyl_wave_geom_new(&app->grid_kinetic, &app->local_kinetic_ext, 0, 0, app->use_gpu);

  // Updater to compute PKPM vars:
  // pressure, primitive variables, and acceleration variables
  // also stores kernels for computing source terms, integrated variables
  double limiter_fac = 1.0/sqrt(3); // default 1/sqrt(3)
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars_ext);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars);
  s->calc_pkpm_vars_ext = gkyl_dg_calc_pkpm_vars_new(&app->grid_kinetic,
    &app->confBasis, &app->local_kinetic_ext, s->eqn_ten_mom, s->ten_mom_geom, limiter_fac, app->use_gpu);
  s->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_new(&app->grid_kinetic,
    &app->confBasis, &app->local_kinetic, s->eqn_ten_mom, s->ten_mom_geom, limiter_fac, app->use_gpu);

  // Updater to compute PKPM dist vars
  // div(p_par b_hat) for self-consistent total pressure force and distribution function sources for
  // Laguerre couplings and vperp characteristics 
  gkyl_dg_calc_pkpm_dist_vars_release(s->calc_pkpm_dist_vars);
  s->calc_pkpm_dist_vars = gkyl_dg_calc_pkpm_dist_vars_new(&s->grid, &app->confBasis, app->use_gpu);

  // Updater to compute PKPM moments (rho, p_par, p_perp, M1) and diagnostics
  gkyl_array_release(s->pkpm_moms);
  gkyl_dg_updater_moment_pkpm_release(s->pkpm_mom_calc);
  s->pkpm_moms = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  s->pkpm_mom_calc = gkyl_dg_updater_moment_pkpm_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local_kinetic, &s->local_vel, s->info.mass, false, app->use_gpu);

  // Diagnostic moments (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
  gkyl_array_release(s->pkpm_moms_diag);
  gkyl_dg_updater_moment_pkpm_release(s->pkpm_mom_calc_diag);
  s->pkpm_moms_diag = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  s->pkpm_mom_calc_diag = gkyl_dg_updater_moment_pkpm_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local_kinetic, &s->local_vel, s->info.mass, true, app->use_gpu);

  // Updater to compute flux contributions to Euler equations, [div(P)-grad(p), (P - pI) : grad(u), div(q)]
  gkyl_dg_calc_flash_coupling_release(s->calc_flash_coupling);
  s->calc_flash_coupling = gkyl_dg_calc_flash_coupling_new(&s->grid, &app->grid_kinetic, &app->confBasis);

  // Updater for kinetic RHS.
  gkyl_dg_updater_flash_release(s->updater);
  bool is_zero_flux[GKYL_MAX_DIM] = {false};
  struct gkyl_dg_vlasov_pkpm_auxfields aux = {
    .bvar = app->field->bvar, .bvar_surf = app->field->bvar_surf, 
    .pkpm_prim = s->pkpm_prim, .pkpm_prim_surf = s->pkpm_prim_surf, 
    .max_b = app->field->max_b, .pkpm_lax = s->pkpm_lax, 
    .div_b = app->field->div_b, .pkpm_accel_vars = s->pkpm_accel, 
    .g_dist_source = s->g_dist_source};
  s->updater = gkyl_dg_updater_flash_new(&app->grid_kinetic, &s->grid,
    &app->confBasis, &app->basis, &app->local_kinetic, &s->local, is_zero_flux, &aux);
}

// Compute PKPM diagnostic moments.
void
flash_species_calc_pkpm_diags(gkyl_flash_app *app, struct flash_species *species,
  const struct gkyl_array *fin)
{
  gkyl_array_clear(species->pkpm_moms, 0.0);
  gkyl_array_clear(species->pkpm_moms_diag, 0.0);

  // Compute rho, p_par, and p_perp
  gkyl_dg_updater_moment_pkpm_advance(species->pkpm_mom_calc, 
    &species->local_ext, &app->local_kinetic_ext, fin, species->pkpm_moms);

  // Compute diagnostic moments (rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
  gkyl_dg_updater_moment_pkpm_advance(species->pkpm_mom_calc_diag, 
    &species->local, &app->local_kinetic_ext, fin, species->pkpm_moms_diag);
    
  // Compute all vars for pkpm_prim (ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp, 3*T_xx/m, 3*T_yy/m, 3*T_zz/m)
  flash_species_calc_pkpm_vars(app, species, species->f);
}

void
flash_species_calc_L2(gkyl_flash_app *app, double tm, const struct flash_species *species)
{
  gkyl_dg_calc_l2_range(&app->basis, 0, species->L2_f, 0, species->f, species->local);
  gkyl_array_scale_range(species->L2_f, species->grid.cellVolume, &species->local);
  
  double L2[1] = { 0.0 };
  gkyl_array_reduce_range(L2, species->L2_f, GKYL_SUM, &species->local);

  double L2_global[1] = { 0.0 };
  gkyl_comm_allreduce_host(app->comm_kinetic, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);
  
  gkyl_dynvec_append(species->integ_L2_f, tm, L2_global);  
}

// Release memory
void
flash_species_release(gkyl_flash_app *app, struct flash_species *s)
{
  gkyl_array_release(s->f); 
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->pkpm_div_ppar);
  gkyl_array_release(s->pkpm_prim);
  gkyl_array_release(s->pkpm_p_ij);
  gkyl_array_release(s->cell_avg_prim);
  gkyl_array_release(s->pkpm_prim_surf);
  gkyl_array_release(s->pkpm_lax);
  gkyl_array_release(s->pkpm_penalization);
  gkyl_array_release(s->pkpm_accel);
  gkyl_array_release(s->pkpm_moms);
  gkyl_array_release(s->pkpm_moms_diag);
  gkyl_array_release(s->flash_coupling_fluxes);
  gkyl_array_release(s->fluid_rhou);
  gkyl_array_release(s->g_dist_source);
  gkyl_array_release(s->F_k_m_1);
  gkyl_array_release(s->F_k_p_1);
  gkyl_array_release(s->cflrate_f);
  gkyl_array_release(s->L2_f);
  gkyl_dynvec_release(s->integ_L2_f);

  flash_species_lbo_release(&s->lbo);
  gkyl_wave_geom_release(s->ten_mom_geom);
  gkyl_wv_eqn_release(s->eqn_ten_mom);
  gkyl_dg_calc_pkpm_dist_vars_release(s->calc_pkpm_dist_vars);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars_ext);
  gkyl_dg_calc_flash_coupling_release(s->calc_flash_coupling);
  gkyl_dg_updater_moment_pkpm_release(s->pkpm_mom_calc);
  gkyl_dg_updater_moment_pkpm_release(s->pkpm_mom_calc_diag);
  gkyl_dg_updater_flash_release(s->updater);
  gkyl_free(s->omegaCFL_ptr);
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(s->bc_lo[d]);
    gkyl_bc_basic_release(s->bc_up[d]);
    gkyl_bc_basic_release(s->bc_conf_lo[d]);
    gkyl_bc_basic_release(s->bc_conf_up[d]);
  }
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  gkyl_comm_release(s->comm);

  if (s->init_dist_from_arrays) {
    gkyl_array_release(s->init_dist_moms);
  }
}