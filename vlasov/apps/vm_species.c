#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <time.h>

// initialize species object
void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
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

  s->model_id = s->info.model_id; // Save model type.
  s->field_id = app->has_field ? app->field->field_id : GKYL_FIELD_NULL; // Save field type.

  // allocate distribution function arrays
  s->f = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu)
    s->f_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

  // allocate cflrate (scalar array)
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);

  if (app->use_gpu)
    s->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else 
    s->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // allocate array to store q/m*(E,B) or potentials (phi, A) depending on equation system
  s->qbym = s->info.charge/s->info.mass;
  if (s->field_id == GKYL_FIELD_E_B)
    s->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  else if (s->field_id != GKYL_FIELD_NULL) {
    s->qmem = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);
    s->qmem_ext = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
  }

  // Determine which directions are not periodic.
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  bool is_zero_flux[2*GKYL_MAX_DIM] = {false}; // Default: no zero-flux BCs in any direction.
  for (int dir=0; dir<app->cdim; ++dir) {
    s->lower_bc[dir].type = s->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const struct gkyl_vlasov_bcs *bc;
      if (dir == 0)
        bc = &s->info.bcx;
      else if (dir == 1)
        bc = &s->info.bcy;
      else
        bc = &s->info.bcz;

      s->lower_bc[dir] = bc->lower;
      s->upper_bc[dir] = bc->upper;
      if (s->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir] = true;
      }
      if (s->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir+pdim] = true;
      }
    }
  }

  if (s->model_id  == GKYL_MODEL_SR) {
    // Allocate special relativistic variables gamma and its inverse
    s->gamma = mkarr(app->use_gpu, app->velBasis.num_basis, s->local_vel.volume);
    s->gamma_host = s->gamma;
    s->gamma_inv = mkarr(app->use_gpu, app->velBasis.num_basis, s->local_vel.volume);
    s->gamma_inv_host = s->gamma_inv;

    if (app->use_gpu) {
      s->gamma_host = mkarr(false, app->velBasis.num_basis, s->local_vel.volume);
      s->gamma_inv_host = mkarr(false, app->velBasis.num_basis, s->local_vel.volume);
    }
    s->sr_vars = gkyl_dg_calc_sr_vars_new(&s->grid, &s->grid_vel,
      &app->confBasis,  &app->velBasis, &app->local, &s->local_vel, app->use_gpu);
    // Project gamma and its inverse
    gkyl_calc_sr_vars_init_p_vars(s->sr_vars, s->gamma, s->gamma_inv);

    struct gkyl_dg_vlasov_sr_auxfields aux_inp = {.qmem = s->qmem, .gamma = s->gamma};

    // create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }
  else if (s->model_id == GKYL_MODEL_CANONICAL_PB || s->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    // Allocate arrays for specified hamiltonian
    s->hamil = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
    s->hamil_host = s->hamil;
    if (app->use_gpu){
      s->hamil_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    }

    // Allocate arrays for specified metric inverse
    s->h_ij = mkarr(app->use_gpu, app->confBasis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    s->h_ij_host = s->h_ij;
    if (app->use_gpu){
      s->h_ij_host = mkarr(false, app->confBasis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    }

    // Allocate arrays for specified metric inverse
    s->h_ij_inv = mkarr(app->use_gpu, app->confBasis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    s->h_ij_inv_host = s->h_ij_inv;
    if (app->use_gpu){
      s->h_ij_inv_host = mkarr(false, app->confBasis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    }

    // Allocate arrays for specified metric determinant
    s->det_h = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    s->det_h_host = s->det_h;
    if (app->use_gpu){
      s->det_h_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    }

    // Evaluate specified hamiltonian function at nodes to insure continuity of hamiltoniam
    struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&s->grid, &app->basis, 1, s->info.hamil, s->info.hamil_ctx);
    gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &s->local_ext, s->hamil_host);
    if (app->use_gpu){
      gkyl_array_copy(s->hamil, s->hamil_host);
    }
    gkyl_eval_on_nodes_release(hamil_proj);

    // Evaluate specified metric function at nodes to insure continuity
    struct gkyl_eval_on_nodes* h_ij_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis, vdim*(vdim+1)/2, s->info.h_ij, s->info.h_ij_ctx);
    gkyl_eval_on_nodes_advance(h_ij_proj, 0.0, &app->local, s->h_ij_host);
    if (app->use_gpu){
      gkyl_array_copy(s->h_ij, s->h_ij_host);
    }
    gkyl_eval_on_nodes_release(h_ij_proj);

    // Evaluate specified inverse metric function at nodes to insure continuity of the inverse 
    struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis, vdim*(vdim+1)/2, s->info.h_ij_inv, s->info.h_ij_inv_ctx);
    gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &app->local, s->h_ij_inv_host);
    if (app->use_gpu){
      gkyl_array_copy(s->h_ij_inv, s->h_ij_inv_host);
    }
    gkyl_eval_on_nodes_release(h_ij_inv_proj);

    // Evaluate specified determinant metric function at nodes to insure continuity of the determinant
    struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis, 1, s->info.det_h, s->info.det_h_ctx);
    gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &app->local, s->det_h_host);
    if (app->use_gpu){
      gkyl_array_copy(s->det_h, s->det_h_host);
    }
    gkyl_eval_on_nodes_release(det_h_proj);

    // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
    struct gkyl_basis surf_basis, surf_quad_basis;
    if (app->basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
      // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
      // sizes in each direction.
      gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
      gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
    } 
    else {
      gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
      gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);
    }

    int alpha_surf_sz = (cdim + vdim)*surf_basis.num_basis; 
    int sgn_alpha_surf_sz = (cdim + vdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

    // allocate arrays to store fields: 
    // 1. alpha_surf (surface phase space velocity)
    // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
    // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
    s->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, s->local_ext.volume);
    s->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, s->local_ext.volume);
    s->const_sgn_alpha = mk_int_arr(app->use_gpu, (cdim + vdim), s->local_ext.volume);

    // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
    struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&s->grid, 
      &app->confBasis, &app->basis, app->use_gpu);
    gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &app->local, &s->local, &s->local_ext, s->hamil,
      s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
    gkyl_dg_calc_canonical_pb_vars_release(calc_vars);

    struct gkyl_dg_canonical_pb_auxfields aux_inp = {.hamil = s->hamil, .alpha_surf = s->alpha_surf, 
      .sgn_alpha_surf = s->sgn_alpha_surf, .const_sgn_alpha = s->const_sgn_alpha};

    //create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }
  else {
    if (s->field_id == GKYL_FIELD_NULL || s->field_id == GKYL_FIELD_E_B) {
      struct gkyl_dg_vlasov_auxfields aux_inp = {.field = s->qmem, .cot_vec = 0, 
        .alpha_surf = 0, .sgn_alpha_surf = 0, .const_sgn_alpha = 0 };
      s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
        &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);
    }
    else {
      struct gkyl_dg_vlasov_poisson_auxfields aux_inp = {.potentials = s->qmem, .fields_ext = s->qmem_ext};
      s->slvr = gkyl_dg_updater_vlasov_poisson_new(&s->grid, &app->confBasis, &app->basis,
        &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);
    }
  }

  // acquire equation object
  if (s->field_id == GKYL_FIELD_NULL || s->field_id == GKYL_FIELD_E_B)
    s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);
  else
    s->eqn_vlasov = gkyl_dg_updater_vlasov_poisson_acquire_eqn(s->slvr);

  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, GKYL_F_MOMENT_M1, false);
  // allocate date for density (for use in charge density accumulation and weak division for V_drift)
  vm_species_moment_init(app, s, &s->m0, GKYL_F_MOMENT_M0, false);
  // allocate data for integrated moments
  vm_species_moment_init(app, s, &s->integ_moms,
    s->model_id == GKYL_MODEL_SR? GKYL_F_MOMENT_M0ENERGYM3 : GKYL_F_MOMENT_M0M1M2, true);

  // allocate data for diagnostic moments
  int ndm = s->info.num_diag_moments;
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m], false);

  // array for storing f^2 in each cell
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  if (app->use_gpu) {
    s->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and f^2
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
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

  // initialize projection routine for initial conditions
  s->num_init = s->info.num_init;
  for (int k=0; k<s->num_init; k++) {
    vm_species_projection_init(app, s, s->info.projection[k], &s->proj_init[k]);
  }
  // If the number of initial condition functions > 1, make a temporary array for accumulation
  if (s->num_init > 1) {
    s->f_tmp = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  }

  // set species source id
  s->source_id = s->info.source.source_id;
  if (s->source_id == GKYL_BFLUX_SOURCE) {
    s->calc_bflux = true;
  }
  
  // determine collision type to use in vlasov update
  s->collision_id = s->info.collisions.collision_id;
  s->lte = (struct vm_lte) { };
  // initialize empty collision structs so inputs of structs are set to 0
  s->lbo = (struct vm_lbo_collisions) { };
  s->bgk = (struct vm_bgk_collisions) { };
  if (s->info.output_f_lte){
    // Always have correct moments on for the f_lte output
    struct correct_all_moms_inp corr_inp = { .correct_all_moms = true, 
      .max_iter = s->info.max_iter, .iter_eps = s->info.iter_eps, 
      .use_last_converged = s->info.use_last_converged };
    vm_species_lte_init(app, s, &s->lte, corr_inp);
  }
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_init(app, s, &s->lbo);
  }
  else if (s->collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_init(app, s, &s->bgk);
  }

  // determine radiation type to use in vlasov update
  s->radiation_id = s->info.radiation.radiation_id;
  s->rad = (struct vm_rad_drag) { };
  if (s->radiation_id == GKYL_VM_COMPTON_RADIATION) {
    vm_species_radiation_init(app, s, &s->rad);
  }

  // Local lower/upper skin and ghost ranges (for applying BCs).
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&s->lower_skin[dir], &s->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    gkyl_skin_ghost_ranges(&s->upper_skin[dir], &s->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
  }

  // Allocate buffer for applying BCs.
  long buff_sz = 0;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(s->lower_skin[dir].volume, s->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  // Buffers for fixed function BCs on distribution function.
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {

    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      s->emit_lo = true;
      s->calc_bflux = true;
      vm_species_emission_init(app, &s->bc_emission_lo, d, GKYL_LOWER_EDGE, s->lower_bc[d].aux_ctx);
    }
    else {
      if (s->lower_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (s->lower_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (s->lower_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_DISTF_REFLECT;
      else if (s->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.basis,
        &s->lower_skin[d], &s->lower_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }

    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      s->emit_up = true;
      s->calc_bflux = true;
      vm_species_emission_init(app, &s->bc_emission_up, d, GKYL_UPPER_EDGE, s->upper_bc[d].aux_ctx);
    }
    else {
      if (s->upper_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (s->upper_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (s->upper_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_DISTF_REFLECT;
      else if (s->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.basis,
        &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }
  }
  if (s->calc_bflux) {
    // intitalize boundary flux updater if we need boundary fluxes
    vm_species_bflux_init(app, s, &s->bflux);
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0)
{
  if (species->num_init > 1) {
    gkyl_array_clear(species->f, 0.0); 
    for (int k=0; k<species->num_init; k++) {
      vm_species_projection_calc(app, species, &species->proj_init[k], species->f_tmp, t0);
      gkyl_array_accumulate(species->f, 1.0, species->f_tmp);
    }
    // Free the temporary array now that initial conditions are complete
    gkyl_array_release(species->f_tmp);
  }
  else {
    vm_species_projection_calc(app, species, &species->proj_init[0], species->f, t0);
  }

  // Pre-compute applied acceleration in case it's time-independent
  vm_species_calc_app_accel(app, species, t0);

  // we are pre-computing source for now as it is time-independent
  vm_species_source_calc(app, species, &species->src, t0);

  if (species->calc_bflux) {
    vm_species_bflux_rhs(app, species, &species->bflux, species->f, species->f1);
  }

  // Optional runtime configuration to use BGK collisions but with fixed input 
  // temperature relaxation based on the initial temperature value. 
  if (species->bgk.fixed_temp_relax) {
    vm_species_bgk_moms_fixed_temp(app, species, &species->bgk, species->f);
  }
  
  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  if (species->lower_bc[0].type == GKYL_SPECIES_FIXED_FUNC)
    gkyl_bc_basic_buffer_fixed_func(species->bc_lo[0], species->bc_buffer_lo_fixed, species->f);
  if (species->upper_bc[0].type == GKYL_SPECIES_FIXED_FUNC)
    gkyl_bc_basic_buffer_fixed_func(species->bc_up[0], species->bc_buffer_up_fixed, species->f);
}

void
vm_species_calc_app_accel(gkyl_vlasov_app *app, struct vm_species *species, double tm)
{
  if (species->has_app_accel) {
    gkyl_proj_on_basis_advance(species->app_accel_proj, tm, &app->local_ext, species->app_accel_host);
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(species->app_accel, species->app_accel_host);
    }
  }
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (species->field_id  == GKYL_FIELD_NULL || species->field_id  == GKYL_FIELD_E_B) {
    if (species->field_id  == GKYL_FIELD_E_B) {
      gkyl_array_set(species->qmem, species->qbym, em);

      // Accumulate applied acceleration and/or q/m*(external electromagnetic)
      // fields onto qmem to get the total acceleration
      if (species->has_app_accel) {
        gkyl_array_accumulate_range(species->qmem, 1.0, species->app_accel, &app->local);
      }
      if (app->field->has_ext_em) {
        gkyl_array_accumulate_range(species->qmem, species->qbym, app->field->ext_em, &app->local);
      }
    }

    gkyl_dg_updater_vlasov_advance(species->slvr, &species->local, 
      fin, species->cflrate, rhs);
  }
  else {
    if (app->field->has_ext_pot) {
      gkyl_array_set_range(species->qmem, species->qbym, app->field->ext_pot, &app->local);
      gkyl_array_accumulate_offset(species->qmem, species->qbym, app->field->phi, 0);
    }
    else {
      gkyl_array_set_offset(species->qmem, species->qbym, app->field->phi, 0);
      if (app->field->ext_em_evolve) {
        gkyl_array_set_range(species->qmem_ext, species->qbym, app->field->ext_em, &app->local);
      }
    }

    gkyl_dg_updater_vlasov_poisson_advance(species->slvr, &species->local,
      fin, species->cflrate, rhs);
  }

  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  }
  else if (species->collision_id == GKYL_BGK_COLLISIONS && !app->has_implicit_coll_scheme) {
    species->bgk.implicit_step = false;
    vm_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  }

  if (species->calc_bflux) {
    vm_species_bflux_rhs(app, species, &species->bflux, fin, rhs);
  }

  if (species->radiation_id == GKYL_VM_COMPTON_RADIATION) {
    vm_species_radiation_rhs(app, species, &species->rad, fin, rhs);
  }
  
  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, &species->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = species->omegaCfl_ptr[0];
  double omegaCfl = omegaCfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl;
}


// Compute the implicit RHS for species update, returning maximum stable
// time-step.
double
vm_species_rhs_implicit(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (species->collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  }

  if (species->calc_bflux) {
    vm_species_bflux_rhs(app, species, &species->bflux, fin, rhs);
  }
  
  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, &species->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = species->omegaCfl_ptr[0];
  double omegaCfl = omegaCfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl;
}


// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f,
  double tcurr)
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

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(app, &species->bc_emission_lo, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(app, &species->bc_emission_up, f, tcurr);
          break;
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


void
vm_species_calc_L2(gkyl_vlasov_app *app, double tm, const struct vm_species *species)
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
vm_species_coll_tm(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_vlasov_tm tm =
        gkyl_dg_updater_lbo_vlasov_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
vm_species_bgk_niter(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      app->stat.niter_self_bgk_corr[i] = app->species[i].bgk.lte.niter;
    }
  }
}


void
vm_species_tm(gkyl_vlasov_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    struct gkyl_dg_updater_vlasov_tm tm;
    if (s->field_id == GKYL_FIELD_NULL || s->field_id == GKYL_FIELD_E_B)
      tm = gkyl_dg_updater_vlasov_get_tm(s->slvr);
    else
      tm = gkyl_dg_updater_vlasov_poisson_get_tm(s->slvr);
    app->stat.species_rhs_tm += tm.vlasov_tm;
  }
}

void
vm_species_rad_tm(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].radiation_id == GKYL_VM_COMPTON_RADIATION) {
      struct gkyl_dg_updater_rad_vlasov_tm tm =
        gkyl_dg_updater_rad_vlasov_get_tm(app->species[i].rad.rad_slvr);
      app->stat.species_rad_tm += tm.drag_tm;
    }
  }
}

// release resources for species
void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  for (int k=0; k<s->num_init; k++) {
    vm_species_projection_release(app, &s->proj_init[k]);
  }

  if (s->calc_bflux) {
    vm_species_bflux_release(app, &s->bflux);
  }

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  // Release arrays for different types of Vlasov equations
  if (s->model_id  == GKYL_MODEL_SR) {
    gkyl_dg_calc_sr_vars_release(s->sr_vars);
    // release relativistic arrays data
    gkyl_array_release(s->gamma);
    gkyl_array_release(s->gamma_inv);
    if (app->use_gpu) {
      gkyl_array_release(s->gamma_host);
      gkyl_array_release(s->gamma_inv_host);
    }
  }
  else if (s->model_id == GKYL_MODEL_CANONICAL_PB || s->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    gkyl_array_release(s->hamil);
    gkyl_array_release(s->h_ij);
    gkyl_array_release(s->h_ij_inv);
    gkyl_array_release(s->det_h);
    gkyl_array_release(s->alpha_surf);
    gkyl_array_release(s->sgn_alpha_surf);
    gkyl_array_release(s->const_sgn_alpha);
    if (app->use_gpu){
      gkyl_array_release(s->hamil_host);
      gkyl_array_release(s->h_ij_host);
      gkyl_array_release(s->h_ij_inv_host);
      gkyl_array_release(s->det_h_host);
    }
  }

  gkyl_array_release(s->qmem);

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_vlasov);
  if (s->field_id == GKYL_FIELD_NULL || s->field_id == GKYL_FIELD_E_B)
    gkyl_dg_updater_vlasov_release(s->slvr);
  else {
    gkyl_array_release(s->qmem_ext);
    gkyl_dg_updater_vlasov_poisson_release(s->slvr);
  }

  // release moment data
  vm_species_moment_release(app, &s->m1i);
  vm_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  vm_species_moment_release(app, &s->integ_moms); 

  gkyl_array_release(s->L2_f);
  gkyl_dynvec_release(s->integ_L2_f);
  gkyl_dynvec_release(s->integ_diag);
  
  gkyl_array_release(s->app_accel);
  if (s->has_app_accel) {
    if (app->use_gpu)
      gkyl_array_release(s->app_accel_host);

    gkyl_proj_on_basis_release(s->app_accel_proj);
  }

  if (s->source_id) {
    vm_species_source_release(app, &s->src);
  }
  if (s->info.output_f_lte){
    vm_species_lte_release(app, &s->lte);
  }
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_release(app, &s->lbo);
  }
  else if (s->collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_release(app, &s->bgk);
  }

  if (s->radiation_id == GKYL_VM_COMPTON_RADIATION) {
    vm_species_radiation_release(app, &s->rad);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d].type == GKYL_SPECIES_EMISSION)
      vm_species_emission_release(&s->bc_emission_lo);
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if (s->upper_bc[d].type == GKYL_SPECIES_EMISSION)
      vm_species_emission_release(&s->bc_emission_up);
    else 
      gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr);
    gkyl_cu_free(s->red_L2_f);
    gkyl_cu_free(s->red_integ_diag);
  }
  else {
    gkyl_free(s->omegaCfl_ptr);
  }
}
