#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

// Euler and isothermal Euler function pointers for primitive/auxiliary variables, 
// including computing the flow velocity, u, from rhou and rho with weak division, and 
// p, the pressure, from the energy in Euler (and just vth*rho in isothermal Euler),
// and also write method, release method, and method for calculating integrated quantities. 
static void
vm_fluid_species_euler_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *f,
  const struct gkyl_array *fluid)
{
  struct timespec tm = gkyl_wall_clock();

  // Compute flow velocity in both the volume and on surfaces
  if (f->bc_is_absorb) {
    gkyl_dg_calc_fluid_vars_advance(f->calc_fluid_vars,
      fluid, f->cell_avg_prim, f->u, f->u_surf); 
  }
  else {
    gkyl_dg_calc_fluid_vars_advance(f->calc_fluid_vars_ext,
      fluid, f->cell_avg_prim, f->u, f->u_surf); 
  }

  // Compute scalar pressure in the volume and at needed surfaces
  gkyl_dg_calc_fluid_vars_pressure(f->calc_fluid_vars, 
    &app->local_ext, fluid, f->u, f->p, f->p_surf);

  app->stat.fluid_species_vars_tm += gkyl_time_diff_now_sec(tm);  
}

static void
vm_fluid_species_euler_calc_integrated(struct gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm)
{
  double avals_fluid[6], avals_fluid_global[6];  
  gkyl_array_clear(f->integ_mom, 0.0);
  vm_fluid_species_prim_vars(app, f, f->fluid);
  // Euler and isothermal Euler integrated quantities: rho, rhoux, rhouy, rhouz, ke, ie
  // where ke is the kinetic energy 1/2 rhou^2 and ie is the internal energy (for ideal Euler = p/(gas_gamma - 1))
  gkyl_dg_calc_fluid_integrated_vars(f->calc_fluid_vars, &app->local, 
    f->fluid, f->u, f->p, f->integ_mom);
  gkyl_array_scale_range(f->integ_mom, app->grid.cellVolume, &app->local);
  if (app->use_gpu) {
    gkyl_array_reduce_range(f->red_integ_diag, f->integ_mom, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals_fluid, f->red_integ_diag, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(avals_fluid, f->integ_mom, GKYL_SUM, &app->local);
  }

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, avals_fluid, avals_fluid_global);
  gkyl_dynvec_append(f->integ_diag, tm, avals_fluid_global);
}

static void
vm_fluid_species_euler_write(gkyl_vlasov_app *app, struct vm_fluid_species *f, 
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, f->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, f->info.name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(f->fluid_host, f->fluid);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, f->fluid_host, fileNm); 

  // Also write out the primitive variables (u, p)
  vm_fluid_species_prim_vars(app, f, f->fluid);
  const char *fmt_prim = "%s-%s_prim_vars_%d.gkyl";
  int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, f->info.name, frame);
  char fileNm_prim[sz_prim+1]; // ensures no buffer overflow
  snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, f->info.name, frame);

  // copy data to single array and then from device to host (if on GPUs) before writing it out
  gkyl_array_set(f->prim_vars, 1.0, f->u); 
  gkyl_array_set_offset(f->prim_vars, 1.0, f->p, 3*app->confBasis.num_basis); 
  if (app->use_gpu) {
    gkyl_array_copy(f->prim_vars_host, f->prim_vars);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, f->prim_vars_host, fileNm_prim); 

  vlasov_array_meta_release(mt);       
}

static void
vm_fluid_species_euler_release(const gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  gkyl_array_release(f->u);
  gkyl_array_release(f->p);
  gkyl_array_release(f->prim_vars); 
  gkyl_array_release(f->prim_vars_host);
  gkyl_array_release(f->cell_avg_prim);
  gkyl_array_release(f->u_surf);
  gkyl_array_release(f->p_surf);
  gkyl_array_release(f->ke_old);
  gkyl_array_release(f->ke_new);
  gkyl_dg_calc_fluid_vars_release(f->calc_fluid_vars);
  gkyl_dg_calc_fluid_vars_release(f->calc_fluid_vars_ext);
}

// Initialize the necessary pieces to solve Euler's or isothermal Euler's equations with DG. 
static void
vm_fluid_species_euler_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  // allocate array to store fluid velocity (ux, uy, uz) and pressure
  // For isothermal Euler, p : (vth*rho)
  // For Euler, p : (gamma - 1)*(E - 1/2 rho u^2)
  f->u = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  f->p = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  // Combined (u, p) primitive variables array for I/O. 
  f->prim_vars = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);
  f->prim_vars_host = app->use_gpu ? mkarr(false, 4*app->confBasis.num_basis, app->local_ext.volume)
                              : gkyl_array_acquire(f->prim_vars);

  // boolean array for if we are only using the cell average for primitive variables
  f->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

  // Allocate arrays for kinetic energy at old and new time steps.
  // These are used because implicit source solve updates momentum but does not affect 
  // the pressure, so we can construct the updated energy from the updated momentum.
  f->ke_old = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    
  f->ke_new = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    

  int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
  // Surface primitive variables (2*cdim*3 components). Ordered as:
  // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
  //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
  f->u_surf = mkarr(app->use_gpu, 2*app->cdim*3*Nbasis_surf, app->local_ext.volume);
  // Surface pressure tensor (2*cdim components). Ordered as:
  // [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr]
  f->p_surf = mkarr(app->use_gpu, 2*app->cdim*Nbasis_surf, app->local_ext.volume);

  // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_fluid_vars.h)
  double limiter_fac = f->info.limiter_fac == 0 ? 0.0 : f->info.limiter_fac;

  // updater for computing fluid variables: flow velocity and pressure
  // also stores kernels for computing source terms, integrated variables
  // Two instances, one over extended range and one over local range for ease of handling boundary conditions
  f->calc_fluid_vars_ext = gkyl_dg_calc_fluid_vars_new(f->equation, app->geom, &app->confBasis, &app->local_ext, limiter_fac, app->use_gpu);
  f->calc_fluid_vars = gkyl_dg_calc_fluid_vars_new(f->equation, app->geom, &app->confBasis, &app->local, limiter_fac, app->use_gpu); 

  struct gkyl_dg_euler_auxfields aux_inp = {.u = f->u, .p = f->p, 
    .u_surf = f->u_surf, .p_surf = f->p_surf};
  f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
    &app->local, f->equation, app->geom, &aux_inp, app->use_gpu);

  // Euler and isothermal Euler integrated quantities: rho, rhoux, rhouy, rhouz, ke, ie
  // where ke is the kinetic energy 1/2 rhou^2 and ie is the internal energy (for ideal Euler = p/(gas_gamma - 1))
  f->integ_mom = mkarr(app->use_gpu, 6, app->local_ext.volume);
  if (app->use_gpu) {
    f->red_integ_diag = gkyl_cu_malloc(sizeof(double[6]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  f->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_integ_write_call = true;

  f->prim_vars_func = vm_fluid_species_euler_prim_vars; 
  f->calc_integrated_mom_func = vm_fluid_species_euler_calc_integrated;  
  f->write_func = vm_fluid_species_euler_write; 
  f->release_func = vm_fluid_species_euler_release; 
}

// Advection function pointers for primitive/auxiliary variables, 
// and also write method, release method, and method for calculating integrated quantities. 
static void
vm_fluid_species_advect_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *f,
  const struct gkyl_array *fluid)
{
  // No primitive variables in advection equation. 
}

static void
vm_fluid_species_advect_calc_integrated(struct gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm)
{
  double avals_fluid[1], avals_fluid_global[1];  
  // Advection equation integrated quantity: f, f^2
  double advect_int[2];

  // First calculate f
  gkyl_array_clear(f->integ_mom, 0.0); 
  gkyl_dg_calc_average_range(&app->confBasis, 0, f->integ_mom, 0, f->fluid, app->local);
  gkyl_array_scale_range(f->integ_mom, app->grid.cellVolume, &app->local);
  if (app->use_gpu) {
    gkyl_array_reduce_range(f->red_integ_diag, f->integ_mom, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals_fluid, f->red_integ_diag, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(avals_fluid, f->integ_mom, GKYL_SUM, &app->local);
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_fluid, avals_fluid_global);
  advect_int[0] = avals_fluid_global[0];

  // Now calculate f^2
  gkyl_array_clear(f->integ_mom, 0.0); 
  gkyl_dg_calc_l2_range(&app->confBasis, 0, f->integ_mom, 0, f->fluid, app->local);
  gkyl_array_scale_range(f->integ_mom, app->grid.cellVolume, &app->local);
  if (app->use_gpu) {
    gkyl_array_reduce_range(f->red_integ_diag, f->integ_mom, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals_fluid, f->red_integ_diag, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(avals_fluid, f->integ_mom, GKYL_SUM, &app->local);
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_fluid, avals_fluid_global);
  advect_int[1] = avals_fluid_global[0];

  gkyl_dynvec_append(f->integ_diag, tm, advect_int);
}

static void
vm_fluid_species_advect_write(gkyl_vlasov_app *app, struct vm_fluid_species *f, 
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, f->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, f->info.name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(f->fluid_host, f->fluid);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, f->fluid_host, fileNm);  

  // If frame = 0, as part of initial conditions also write out the applied advection. 
  if (frame == 0) {
    const char *fmt_advect = "%s-%s_advect_%d.gkyl";
    int sz_advect = gkyl_calc_strlen(fmt_advect, app->name, f->info.name, frame);
    char fileNm_advect[sz_advect+1]; // ensures no buffer overflow
    snprintf(fileNm_advect, sizeof fileNm_advect, fmt_advect, app->name, f->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(f->app_advect_host, f->app_advect);
    }
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
      mt, f->app_advect_host, fileNm_advect);      
  }

  vlasov_array_meta_release(mt);   
}

static void
vm_fluid_species_advect_release(const gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  gkyl_array_release(f->app_advect);
  if (app->use_gpu) {
    gkyl_array_release(f->app_advect_host);
  }
}

// Initialize the necessary pieces to solve an advection equation with DG. 
static void
vm_fluid_species_advect_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  // setup FEM representation of applied advection 
  // 
  f->app_advect = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  f->app_advect_host = f->app_advect;
  if (app->use_gpu) {
    f->app_advect_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
  }

  // Evaluate specified advection function at nodes to insure continuity of advection velocity
  struct gkyl_eval_on_nodes* app_advect_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis, 3, 
    f->info.advection.velocity, f->info.advection.velocity_ctx);
  gkyl_eval_on_nodes_advance(app_advect_proj, 0.0, &app->local_ext, f->app_advect_host);
  if (app->use_gpu) {
    gkyl_array_copy(f->app_advect, f->app_advect_host);
  }
  gkyl_eval_on_nodes_release(app_advect_proj);

  struct gkyl_dg_advection_auxfields aux_inp = {.u_i = f->app_advect};
  f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
    &app->local, f->equation, app->geom, &aux_inp, app->use_gpu);

  // array for storing integrated quantities in each cell = f, f^2
  // We compute each separate and then write them to the dynvector f->integ_diag
  f->integ_mom = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu) {
    f->red_integ_diag = gkyl_cu_malloc(sizeof(double[1]));
  }
  // allocate dynamic-vector to store all-reduced integrated quantities
  f->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 2);
  f->is_first_integ_write_call = true;

  f->prim_vars_func = vm_fluid_species_advect_prim_vars; 
  f->calc_integrated_mom_func = vm_fluid_species_advect_calc_integrated;  
  f->write_func = vm_fluid_species_advect_write; 
  f->release_func = vm_fluid_species_advect_release; 
}

// Canonical Poisson Bracket fluid function pointers for primitive/auxiliary variables, 
// including computing grad^2 phi = f, where f is (one of) the evolved quantities by 
// the canonical Poisson bracket system such as vorticity in incompressible Euler,  
// and also write method, release method, and method for calculating integrated quantities. 
static void
vm_fluid_species_can_pb_fluid_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *f,
  const struct gkyl_array *fluid)
{
  struct timespec wst = gkyl_wall_clock();

  // Gather RHS of Poisson solve into global array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, 
    fluid, f->poisson_rhs_global);

  // Solve the Poisson problem.
  gkyl_fem_poisson_set_rhs(f->fem_poisson, f->poisson_rhs_global, 0);
  gkyl_fem_poisson_solve(f->fem_poisson, f->phi_global);

  // Copy the portion of global potential corresponding to this MPI process to the local potential.
  gkyl_array_copy_range_to_range(f->phi, f->phi_global, 
    &app->local, &f->global_sub_range);

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst); 
}

static void
vm_fluid_species_can_pb_fluid_calc_integrated(struct gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm)
{ 
  double avals_fluid[1], avals_fluid_global[1];  
  // Canonical Poisson bracket for fluid equations integrated quantities: f, f^2, E
  // where E is the integrated energy = integral(phi*f) since f = grad^2 phi -> E = integral(|grad phi|^2)
  double can_pb_int[3];

  // First calculate f
  gkyl_array_clear(f->integ_mom, 0.0); 
  gkyl_dg_calc_average_range(&app->confBasis, 0, f->integ_mom, 0, f->fluid, app->local);
  gkyl_array_scale_range(f->integ_mom, app->grid.cellVolume, &app->local);
  if (app->use_gpu) {
    gkyl_array_reduce_range(f->red_integ_diag, f->integ_mom, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals_fluid, f->red_integ_diag, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(avals_fluid, f->integ_mom, GKYL_SUM, &app->local);
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_fluid, avals_fluid_global);
  can_pb_int[0] = avals_fluid_global[0];

  // Now calculate f^2
  gkyl_array_clear(f->integ_mom, 0.0); 
  gkyl_dg_calc_l2_range(&app->confBasis, 0, f->integ_mom, 0, f->fluid, app->local);
  gkyl_array_scale_range(f->integ_mom, app->grid.cellVolume, &app->local);
  if (app->use_gpu) {
    gkyl_array_reduce_range(f->red_integ_diag, f->integ_mom, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals_fluid, f->red_integ_diag, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(avals_fluid, f->integ_mom, GKYL_SUM, &app->local);
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_fluid, avals_fluid_global);
  can_pb_int[1] = avals_fluid_global[0];

  // Now calculate E 
  // First solve the Poisson equation for phi 
  vm_fluid_species_prim_vars(app, f, f->fluid);
  // integrate |grad phi|^2
  gkyl_array_integrate_advance(f->calc_can_pb_energy, f->phi,
    app->grid.cellVolume, f->can_pb_energy_fac, &app->local, &app->local, f->red_can_pb_energy);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, f->red_can_pb_energy, f->red_can_pb_energy_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_fluid_global, f->red_can_pb_energy_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    avals_fluid_global[0] = f->red_can_pb_energy_global[0];
  }
  can_pb_int[2] = avals_fluid_global[0];

  gkyl_dynvec_append(f->integ_diag, tm, can_pb_int);

}

static void
vm_fluid_species_can_pb_fluid_write(gkyl_vlasov_app *app, struct vm_fluid_species *f, 
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, f->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, f->info.name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(f->fluid_host, f->fluid);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, f->fluid_host, fileNm);

  // If fluid species is a canonical PB fluid and we are also solving a
  // Poisson equation, also write out phi, the potential. 
  vm_fluid_species_prim_vars(app, f, f->fluid);

  const char *fmt_phi = "%s-phi_%d.gkyl";
  int sz_phi = gkyl_calc_strlen(fmt_phi, app->name, frame);
  char fileNm_phi[sz_phi+1]; // ensures no buffer overflow
  snprintf(fileNm_phi, sizeof fileNm_phi, fmt_phi, app->name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(f->phi_host, f->phi);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, f->phi_host, fileNm_phi);    

  vlasov_array_meta_release(mt);  
}

static void
vm_fluid_species_can_pb_fluid_release(const gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  gkyl_array_release(f->phi);
  gkyl_array_release(f->phi_global);
  gkyl_array_release(f->phi_host);
  gkyl_array_release(f->poisson_rhs_global);
  gkyl_array_release(f->can_pb_n0); 
  gkyl_array_release(f->epsilon);
  if (f->eqn_type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA) {
    gkyl_array_release(f->kSq);
  }
  gkyl_fem_poisson_release(f->fem_poisson);

  gkyl_array_release(f->alpha_surf); 
  gkyl_array_release(f->sgn_alpha_surf); 
  gkyl_array_release(f->const_sgn_alpha); 
  gkyl_dg_calc_canonical_pb_fluid_vars_release(f->calc_can_pb_fluid_vars); 
  gkyl_array_integrate_release(f->calc_can_pb_energy);
  gkyl_array_release(f->can_pb_energy_fac);
  if (app->use_gpu) {
    gkyl_cu_free(f->red_can_pb_energy);
    gkyl_cu_free(f->red_can_pb_energy_global);
  } 
  else {
    gkyl_free(f->red_can_pb_energy);
    gkyl_free(f->red_can_pb_energy_global);
  }
}

// Initialize the necessary pieces to solve an canonical Poisson bracket fluid system with DG 
// such as incompressible Euler, Hasegawa-Mima, and (modified) Hasegawa-Wakatani.
static void
vm_fluid_species_can_pb_fluid_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  assert(app->cdim == 2); // Canonical PB fluid only works with cdim = 2. 

  // canonical PB fluids couple to a Poisson equation; initialize that
  // here for now, though it might be better to utilize the vp_field 
  // infrastructure once that becomes more flexible (JJ: 02/03/2025)
  // Allocate arrays for potential and global rhs of Poisson solve.
  f->phi = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->phi_global = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
  f->poisson_rhs_global = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->global_ext.volume);

  // Host potential for  I/O.
  f->phi_host = app->use_gpu ? mkarr(false, app->confBasis.num_basis, app->local_ext.volume)
                              : gkyl_array_acquire(f->phi);

  // Initialize background gradient which drives turbulence in some fluid systems such as
  // Hasegawa-Mima and Hasegawa-Wakatani. This background is included as a source via 
  // the Poisson bracket: source = {phi, n0} where {., .} is the canonical bracket.
  f->can_pb_n0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->can_pb_n0, 0.0); 
  if (f->eqn_type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA || f->eqn_type == GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI) {
    struct gkyl_array* can_pb_n0_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    // Evaluate specified background gradient function at nodes to insure continuity of gradient
    struct gkyl_eval_on_nodes* can_pb_n0_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis, 1, 
      f->info.can_pb_n0, f->info.can_pb_n0_ctx);
    gkyl_eval_on_nodes_advance(can_pb_n0_proj, 0.0, &app->local_ext, can_pb_n0_host);
    gkyl_array_copy(f->can_pb_n0, can_pb_n0_host);
    gkyl_eval_on_nodes_release(can_pb_n0_proj); 
    gkyl_array_release(can_pb_n0_host); 
  }

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // Set the permittivity in the Poisson equation. Just a constant value of -1.0
  // so that overall Poisson equation becomes grad^2 phi = RHS since FEM Poisson
  // solver by default assumes Poisson equation of the form -nabla . (epsilon . grad phi) = RHS. 
  f->epsilon = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
  gkyl_array_clear(f->epsilon, 0.0);
  gkyl_array_shiftc(f->epsilon, -1.0*pow(sqrt(2.0),app->cdim), 0);

  // Create Poisson solver. Only supports periodic boundary conditions for now. 
  struct gkyl_poisson_bc poisson_bcs = {
    .lo_type = { GKYL_POISSON_PERIODIC, GKYL_POISSON_PERIODIC },
    .up_type = { GKYL_POISSON_PERIODIC, GKYL_POISSON_PERIODIC },
  };
  if (f->eqn_type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA) {
    // If Hasegawa-Mima, we solve a Helmholtz equations (grad^2 - 1)phi = zeta
    // where zeta is vorticity (the quantity we are evolving). 
    f->kSq = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
    gkyl_array_clear(f->kSq, 0.0);
    gkyl_array_shiftc(f->kSq, pow(sqrt(2.0),app->cdim), 0);
    f->fem_poisson = gkyl_fem_poisson_new(&app->global, &app->grid, app->confBasis,
      &poisson_bcs, 0, f->epsilon, f->kSq, true, app->use_gpu);
  }
  else {
    f->fem_poisson = gkyl_fem_poisson_new(&app->global, &app->grid, app->confBasis,
      &poisson_bcs, 0, f->epsilon, NULL, true, app->use_gpu);      
  }
  f->has_poisson = true; 

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis; 
  if (app->confBasis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&surf_basis, app->cdim-1, app->poly_order); 
  }
  else if (app->confBasis.b_type == GKYL_BASIS_MODAL_TENSOR) {
    gkyl_cart_modal_tensor(&surf_basis, app->cdim-1, app->poly_order); 
  }
  else {
    assert(false);
  }
  gkyl_cart_modal_tensor(&surf_quad_basis, app->cdim-1, app->poly_order);    
  int alpha_surf_sz = app->cdim*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = app->cdim*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface configuration space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  f->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, app->local_ext.volume);
  f->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, app->local_ext.volume);
  f->const_sgn_alpha = mk_int_arr(app->use_gpu, app->cdim, app->local_ext.volume);

  // Initialize alpha_surf calculator and source functions for different canonical PB fluids. 
  // NOTE: If solving the modified Hasegawa-Wakatani system which subtracts the zonal components
  // of n and phi, the simulation *must not* have any decomposition in y. 
  f->calc_can_pb_fluid_vars = gkyl_dg_calc_canonical_pb_fluid_vars_new(&app->grid, 
    &app->confBasis, &app->local, &app->local_ext, f->equation, app->use_gpu); 

  struct gkyl_dg_canonical_pb_fluid_auxfields aux_inp = {.phi = f->phi, .alpha_surf = f->alpha_surf, 
    .sgn_alpha_surf = f->sgn_alpha_surf, .const_sgn_alpha = f->const_sgn_alpha};
  f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
    &app->local, f->equation, app->geom, &aux_inp, app->use_gpu);

  // Canonical Poisson bracket for fluid equations integrated quantities: f, f^2, E
  // where E is the integrated energy = integral(phi*f) since f = grad^2 phi -> E = integral(|grad phi|^2)
  // We compute each separate and then write them to the dynvector f->integ_diag
  f->integ_mom = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu) {
    f->red_integ_diag = gkyl_cu_malloc(sizeof(double[1]));
    f->red_can_pb_energy = gkyl_cu_malloc(sizeof(double[1]));
    f->red_can_pb_energy_global = gkyl_cu_malloc(sizeof(double[1]));
  }
  else {
    f->red_can_pb_energy = gkyl_malloc(sizeof(double[1]));
    f->red_can_pb_energy_global = gkyl_malloc(sizeof(double[1]));
  }
  // allocate dynamic-vector to store all-reduced integrated quantities
  f->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 3);
  f->is_first_integ_write_call = true; 

  // Specialized updater for integrating |grad phi|^2
  f->can_pb_energy_fac = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_shiftc(f->can_pb_energy_fac, pow(sqrt(2.0),app->cdim), 0); // Sets can_pb_energy_fac = 1.
  f->calc_can_pb_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis,
    1, GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ, app->use_gpu);

  f->prim_vars_func = vm_fluid_species_can_pb_fluid_prim_vars; 
  f->calc_integrated_mom_func = vm_fluid_species_can_pb_fluid_calc_integrated;  
  f->write_func = vm_fluid_species_can_pb_fluid_write; 
  f->release_func = vm_fluid_species_can_pb_fluid_release; 
}

// initialize fluid species object
void
vm_fluid_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  int cdim = app->cdim;
  // Setup equation-specific memory and equation type/number of equations based on input table
  f->eqn_type = f->info.equation->type;
  f->num_equations = f->info.equation->num_equations;
  f->equation = gkyl_wv_eqn_acquire(f->info.equation);
  if (f->eqn_type == GKYL_EQN_ADVECTION) {
    vm_fluid_species_advect_init(vm, app, f);
  }
  else if (f->eqn_type == GKYL_EQN_EULER || f->eqn_type == GKYL_EQN_ISO_EULER) {
    vm_fluid_species_euler_init(vm, app, f);
  }
  else {
    // Equation type is a Canonical Poisson Bracket (PB) fluid such as
    // incompressible Euler, Hasegawa-Mima, or Hasegawa-Wakatani. 
    vm_fluid_species_can_pb_fluid_init(vm, app, f);
  }

  // allocate fluid arrays
  f->fluid = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);
  f->fluid1 = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);
  f->fluidnew = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);

  f->fluid_host = f->fluid;
  if (app->use_gpu) {
    f->fluid_host = mkarr(false, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);
  }

  // Duplicate copy of fluid data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  f->fluid_dup = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu) {
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  }
  else {
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));
  }

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  f->has_diffusion = false;
  f->diffD = NULL;
  if (f->info.diffusion.Dij) {
    f->has_diffusion = true;
    // allocate space for full diffusion tensor 
    int szD = cdim*(cdim+1)/2;

    f->diffD = mkarr(app->use_gpu, szD*app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *diffD_host = f->diffD;
    if (app->use_gpu)
      diffD_host = mkarr(false, szD*app->confBasis.num_basis, app->local_ext.volume);

    gkyl_proj_on_basis *diff_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = szD,
        .eval = f->info.diffusion.Dij,
        .ctx = f->info.diffusion.Dij_ctx
      }
    );
    gkyl_proj_on_basis_advance(diff_proj, 0.0, &app->local_ext, diffD_host);
    if (app->use_gpu) {// note: diffD_host is same as diffD when not on GPUs
      gkyl_array_copy(f->diffD, diffD_host);
      gkyl_array_release(diffD_host);
    }
    // Free projection object
    gkyl_proj_on_basis_release(diff_proj);

    f->diff_slvr_gen = gkyl_dg_updater_diffusion_gen_new(&app->grid, &app->confBasis,
      &app->local, app->use_gpu);
  }
  else if (f->info.diffusion.D) {
    f->has_diffusion = true;
    f->info.diffusion.order = f->info.diffusion.order<2? 2 : f->info.diffusion.order;

    int szD = cdim;
    f->diffD = mkarr(app->use_gpu, szD, 1);
    struct gkyl_array *diffD_host = f->diffD;
    if (app->use_gpu)
      diffD_host = mkarr(false, szD, 1);
    // Set diffusion coefficient in each direction to input value.
    gkyl_array_clear(diffD_host, 0.);
    for (int d=0; d<cdim; d++) gkyl_array_shiftc(diffD_host, f->info.diffusion.D, d);

    if (app->use_gpu) {// note: diffD_host is same as diffD when not on GPUs
      gkyl_array_copy(f->diffD, diffD_host);
      gkyl_array_release(diffD_host);
    }

    const bool is_zero_flux[GKYL_MAX_CDIM] = {false};

    f->diff_slvr = gkyl_dg_updater_diffusion_fluid_new(&app->grid, &app->confBasis,
      true, f->num_equations, NULL, f->info.diffusion.order, &app->local, is_zero_flux, app->use_gpu);
  }

  // Initialize applied acceleration for use in force update. 
  // Always used by fluid implicit sources, so always initialize.
  f->app_accel = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->app_accel, 0.0);
  f->has_app_accel = false;
  f->app_accel_evolve = false;
  // setup applied acceleration
  if (f->info.app_accel) {
    f->has_app_accel = true;
    if (f->info.app_accel_evolve) {
      f->app_accel_evolve = f->info.app_accel_evolve;
    }

    f->app_accel_host = f->app_accel;
    if (app->use_gpu) {
      f->app_accel_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
    }
    f->app_accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      3, f->info.app_accel, f->info.app_accel_ctx);
  }

  // set species source id
  f->source_id = f->info.source.source_id;

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
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
  f->bc_buffer = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, buff_sz);

  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  f->bc_is_absorb = false;
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_EULER_REFLECT;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_NO_SLIP) {
      bctype = GKYL_BC_EULER_NO_SLIP;
    }

    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->lower_skin[d], &app->lower_ghost[d], f->fluid->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_EULER_REFLECT;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_NO_SLIP) {
      bctype = GKYL_BC_EULER_NO_SLIP;
    }

    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->upper_skin[d], &app->upper_ghost[d], f->fluid->ncomp, app->cdim, app->use_gpu);
  }
}

void
vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0)
{
  int poly_order = app->poly_order;

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, fluid_species->num_equations, fluid_species->info.init, fluid_species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &app->local_ext, fluid_species->fluid_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu) {
    gkyl_array_copy(fluid_species->fluid, fluid_species->fluid_host);
  }
  // Apply limiter at t=0 to insure slopes are well-behaved at beginning of simulation
  vm_fluid_species_limiter(app, fluid_species, fluid_species->fluid);

  // Pre-compute applied acceleration in case it's time-independent
  vm_fluid_species_calc_app_accel(app, fluid_species, t0);

  // we are pre-computing source for now as it is time-independent
  vm_fluid_species_source_calc(app, fluid_species, t0);
}

void
vm_fluid_species_calc_app_accel(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm)
{
  if (fluid_species->has_app_accel) {
    gkyl_proj_on_basis_advance(fluid_species->app_accel_proj, tm, &app->local_ext, fluid_species->app_accel_host);
    // note: app_accel_host is same as app_accel when not on GPUs
    if (app->use_gpu) {
      gkyl_array_copy(fluid_species->app_accel, fluid_species->app_accel_host);
    }
  }
}

void
vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid)
{
  fluid_species->prim_vars_func(app, fluid_species, fluid);
}

void
vm_fluid_species_limiter(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  struct gkyl_array *fluid)
{
  if (fluid_species->eqn_type == GKYL_EQN_EULER || fluid_species->eqn_type == GKYL_EQN_ISO_EULER) {  
    struct timespec tm = gkyl_wall_clock();

    // Limit the slopes of the solution
    gkyl_dg_calc_fluid_vars_limiter(fluid_species->calc_fluid_vars, &app->local, fluid);

    app->stat.fluid_species_vars_tm += gkyl_time_diff_now_sec(tm);

    // Apply boundary conditions after limiting solution
    vm_fluid_species_apply_bc(app, fluid_species, fluid);
  }
}

// Compute the RHS for fluid species update, returning maximum stable
// time-step.
double
vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  double omegaCfl = 1/DBL_MAX;

  gkyl_array_clear(fluid_species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // If we are solving a Poisson equation, need to compute 
  // surface characteristics from potential and source update.
  if (fluid_species->has_poisson) {
    struct timespec tm = gkyl_wall_clock();

    // Compute the surface characteristics from the potential. 
    gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf(fluid_species->calc_can_pb_fluid_vars, 
      &app->local, &app->local_ext, fluid_species->phi, 
      fluid_species->alpha_surf, fluid_species->sgn_alpha_surf, fluid_species->const_sgn_alpha); 

    // Increment the source contribution for certain canonical PB fluids onto the RHS. 
    gkyl_canonical_pb_fluid_vars_source(fluid_species->calc_can_pb_fluid_vars, 
      &app->local, fluid_species->phi, fluid_species->can_pb_n0, fluid, rhs); 

    app->stat.fluid_species_vars_tm += gkyl_time_diff_now_sec(tm); 
  }

  gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, 
    &app->local, fluid, fluid_species->cflrate, rhs);

  // Accumulate explicit source contribution, e.g., external forces
  // Only done if there are external forces and no EM fields, as fluid-EM coupling
  // is handled by implicit source solve, see vm_fluid_em_coupling.c. 
  if (fluid_species->has_app_accel && !app->has_field) {
    gkyl_dg_calc_fluid_vars_source(fluid_species->calc_fluid_vars, &app->local, 
      fluid_species->app_accel, fluid, rhs); 
  }

  if (fluid_species->has_diffusion) {
    if (fluid_species->info.diffusion.Dij) {
      gkyl_dg_updater_diffusion_gen_advance(fluid_species->diff_slvr_gen,
        &app->local, fluid_species->diffD, fluid, fluid_species->cflrate, rhs);
    }
    else if (fluid_species->info.diffusion.D) {
      gkyl_dg_updater_diffusion_fluid_advance(fluid_species->diff_slvr,
        &app->local, fluid_species->diffD, fluid, fluid_species->cflrate, rhs);
    }
  }

  gkyl_array_reduce_range(fluid_species->omegaCfl_ptr, fluid_species->cflrate, GKYL_MAX, &app->local);

  double omegaCfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omegaCfl_ho, fluid_species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omegaCfl_ho[0] = fluid_species->omegaCfl_ptr[0];
  }
  omegaCfl = omegaCfl_ho[0];

  app->stat.fluid_species_rhs_tm += gkyl_time_diff_now_sec(wst);

  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for fluid species
void
vm_fluid_species_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f)
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

      switch (fluid_species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
          gkyl_bc_basic_advance(fluid_species->bc_lo[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }

      switch (fluid_species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
          gkyl_bc_basic_advance(fluid_species->bc_up[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);

  app->stat.fluid_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Integrated quantities calculator for different fluid equations.
void 
vm_fluid_species_calc_integrated_mom(gkyl_vlasov_app* app, struct vm_fluid_species *f, double tm)
{
  f->calc_integrated_mom_func(app, f, tm); 
}

// Write method for different fluid equations. 
void
vm_fluid_species_write(gkyl_vlasov_app *app, struct vm_fluid_species *f, double tm, int frame)
{
  f->write_func(app, f, tm, frame);
}

// Release resources for fluid species.
void
vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f)
{
  // Release acquired equation object pointer 
  // along with other equation specific data. 
  gkyl_wv_eqn_release(f->equation); 
  f->release_func(app, f); 

  // Release solver and fluid arrays for update. 
  gkyl_dg_updater_fluid_release(f->advect_slvr);
  gkyl_array_release(f->fluid);
  gkyl_array_release(f->fluid1);
  gkyl_array_release(f->fluidnew);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  gkyl_array_release(f->fluid_dup);

  if (f->has_diffusion) {
    gkyl_array_release(f->diffD);
    if (f->info.diffusion.Dij) {
      gkyl_dg_updater_diffusion_gen_release(f->diff_slvr_gen);
    }
    else if (f->info.diffusion.D) {
      gkyl_dg_updater_diffusion_fluid_release(f->diff_slvr);
    }
  }

  gkyl_array_release(f->integ_mom);
  gkyl_dynvec_release(f->integ_diag);

  gkyl_array_release(f->app_accel);
  if (f->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(f->app_accel_host);
    }
    gkyl_proj_on_basis_release(f->app_accel_proj);
  }

  if (f->source_id) {
    vm_fluid_species_source_release(app, &f->src);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->fluid_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    gkyl_cu_free(f->red_integ_diag);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }
}
