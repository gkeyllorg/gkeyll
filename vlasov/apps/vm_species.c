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

// Configuration-space c2p: sample geometry/metric functions at physical conf
// coordinates on non-uniform meshes via the position map (identity map =>
// identity c2p, so uniform grids are unaffected).
static void
vm_species_c2p_conf(const double *xcomp, double *xphys, void *ctx)
{
  struct vm_proj_c2p_ctx *c = (struct vm_proj_c2p_ctx *) ctx;
  gkyl_vlasov_position_map_eval_mc2p(c->pos_map, xcomp, xphys);
}

// Phase-space c2p: conf coords via the position map, velocity coords via the
// velocity map. Used for the phase-space Hamiltonian nodal projection so H is
// sampled at the physical coordinates the mapped mesh represents (nodes shared
// between cells map to the same physical point, preserving C^0 continuity).
static void
vm_species_c2p_phase(const double *xcomp, double *xphys, void *ctx)
{
  struct vm_proj_c2p_ctx *c = (struct vm_proj_c2p_ctx *) ctx;
  gkyl_vlasov_position_map_eval_mc2p(c->pos_map, xcomp, xphys);
  gkyl_vlasov_velocity_map_eval_c2p(c->vel_map, &xcomp[c->cdim], &xphys[c->cdim]);
}

static void
vm_species_new_hamil(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int vdim = app->vdim;  

  // Allocate arrays for configuration space Poisson tensor (unused)
  int num_pt_indices[3] = { 1 , 6, 18 };
  vms->conf_poisson_tensor = mkarr(app->use_gpu, app->basis.num_basis*num_pt_indices[vdim-1], app->local_ext.volume);
  vms->conf_poisson_tensor_host = vms->conf_poisson_tensor;
  if (app->use_gpu){
    vms->conf_poisson_tensor_host = mkarr(false, app->basis.num_basis*num_pt_indices[vdim-1], app->local_ext.volume);
  }

  // Coordinate-map context for sampling geometry/Hamiltonians at physical
  // coordinates on non-uniform meshes. All uses complete within this function,
  // so a stack-local context is safe.
  struct vm_proj_c2p_ctx c2p_ctx = {
    .cdim = app->cdim,
    .pos_map = vms->pos_map,
    .vel_map = vms->vel_map,
  };

  if (vms->model_id == GKYL_MODEL_DEFAULT || vms->model_id == GKYL_MODEL_SR) {
    // H = v^2/2 is separable (sparse velocity-space expansion); the relativistic
    // H = sqrt(1 + p^2) is a dense velocity-space expansion.
    vms->hamil_id = gkyl_hamil_id_from_model_id(vms->model_id);
    vms->mom_hamil_id = vms->hamil_id;

    // Hamiltonain for computing the moments is only a function of velocity space.
    vms->mom_hamil_range = vms->local_vel;
    vms->mom_hamil = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    vms->gamma_inv = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    gkyl_dg_vlasov_calc_hamil(&vms->grid_vel, &vms->basis_vel, &vms->local_vel, 
      vms->model_id, vms->vel_map, vms->mom_hamil, vms->gamma_inv, app->use_gpu); 

    // As the mom_hamil and hamil are identical, aquire hamil and hamil_range here
    vms->hamil_range = vms->mom_hamil_range;
    vms->hamil = gkyl_array_acquire(vms->mom_hamil);

    // If relativistic, allocate additional updater for computing derived relativistic moments,
    // such as the spatial component of the four-velocity and pressure.
    if (vms->model_id == GKYL_MODEL_SR) {
      vms->sr_vars = gkyl_dg_calc_sr_vars_new(&vms->grid, &vms->grid_vel,
        &app->basis,  &vms->basis_vel, &app->local, &vms->local_vel,
        vms->vel_map, app->use_gpu);
    }
  } 
  else if (vms->model_id == GKYL_MODEL_TRIAD) {

    int cdim = app->cdim;
    int vdim = app->vdim;

    // Triad bracket with the separable H = v^2/2 (sparse velocity-space expansion).
    vms->hamil_id = GKYL_HAMIL_VEL_SPARSE;
    vms->mom_hamil_id = GKYL_HAMIL_VEL_SPARSE;

    // Hamiltonain is only a function of velocity space, non-relativistic only using the
    // same infrastructure as GKYL_MODEL_DEFAULT
    vms->mom_hamil_range = vms->local_vel; 
    vms->mom_hamil = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    vms->gamma_inv = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    gkyl_dg_vlasov_calc_hamil(&vms->grid_vel, &vms->basis_vel, &vms->local_vel, 
      GKYL_MODEL_DEFAULT, vms->vel_map, vms->mom_hamil, vms->gamma_inv, app->use_gpu);

    // As the mom_hamil and hamil are identical, aquire hamil and hamil_range here
    vms->hamil_range = vms->mom_hamil_range;
    vms->hamil = gkyl_array_acquire(vms->mom_hamil);

    // Triad geometry context if needed for preset geometries.
    struct gkyl_triad_geom_ctx preset_geom_ctx = {
      .mass_bh = vms->geom->mass_bh,
      .spin_bh = vms->geom->spin_bh,
    };

    struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
    inp_triad_geom.use_vierbein = vms->info.use_vierbein;
    inp_triad_geom.use_preset_geom = vms->geom->use_preset_geom;
    if ( (vms->info.use_vierbein) && (vms->info.vierbein) && (vms->info.vierbein_gradient) ) {
      inp_triad_geom.eval_vierbein = vms->info.vierbein;
      inp_triad_geom.eval_vierbein_gradient = vms->info.vierbein_gradient;
      inp_triad_geom.eval_vierbein_ctx = vms->info.vierbein_ctx;
      inp_triad_geom.eval_vierbein_gradient_ctx = vms->info.vierbein_gradient_ctx;
    }
    else if (vms->geom->use_preset_geom) {
      inp_triad_geom.triad_preset_geom_type = vms->geom->triad_preset_geom_type;
      inp_triad_geom.eval_vierbein_ctx = &preset_geom_ctx;
      inp_triad_geom.eval_vierbein_gradient_ctx = &preset_geom_ctx;
    }
    else if ((vms->info.triad_basis) && (vms->info.triad_basis_gradient) && (vms->info.cov_tangent_basis))  {
      inp_triad_geom.eval_cov_tangent_basis = vms->info.cov_tangent_basis; 
      inp_triad_geom.eval_triad_basis = vms->info.triad_basis; 
      inp_triad_geom.eval_triad_basis_gradient = vms->info.triad_basis_gradient; 
      inp_triad_geom.eval_cov_tangent_basis_ctx = vms->info.cov_tangent_basis_ctx; 
      inp_triad_geom.eval_triad_basis_ctx = vms->info.triad_basis_ctx; 
      inp_triad_geom.eval_triad_basis_gradient_ctx = vms->info.triad_basis_gradient_ctx; 
    } 
    else {
      // Missing valid geometry inputs for triads
      assert(false);
    }

    // Sample the triad geometry at physical conf coordinates on non-uniform meshes.
    inp_triad_geom.c2p_func = vm_species_c2p_conf;
    inp_triad_geom.c2p_func_ctx = &c2p_ctx;

    // The geometry comes from the tangents and triads
    gkyl_vlasov_triad_geom_new(&app->grid, &app->local, app->basis,
      &vms->grid, &vms->local, vms->basis, inp_triad_geom, vms->conf_poisson_tensor_host);

    // Copy Pi_conf onto the device.
    if (app->use_gpu) {
      gkyl_array_copy(vms->conf_poisson_tensor, vms->conf_poisson_tensor_host); 
    }

  }
  else if (vms->model_id == GKYL_MODEL_TRIAD_GR) {

    int cdim = app->cdim;
    int vdim = app->vdim;  
  

    // NOTE: Two hamiltonians are built here. mom_hamil (1.) is used for local moments and hamil (2.) is 
    // a function of the full phase space used for chateristics.
    // 1. (mom_hamil)  H(\hat{p}_i) = \gamma
    // 2. (hamil)      H(x^i,\hat{p}_i) = \alpha \gamma - \hat{p}_a e_i^a \beta^i 
    
    // The characteristics Hamiltonian is a phase-space expansion; the moment
    // Hamiltonian is the dense velocity-space (SR) expansion.
    vms->hamil_id = GKYL_HAMIL_PHASE;
    vms->mom_hamil_id = GKYL_HAMIL_VEL_DENSE;

    // 1. Build a velocity space only hamiltonian for moments.
    // Hamiltonain is only a function of velocity space, non-relativistic only using the
    // same infrastructure as GKYL_MODEL_DEFAULT
    vms->mom_hamil_range = vms->local_vel;
    vms->mom_hamil = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    vms->gamma_inv = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    gkyl_dg_vlasov_calc_hamil(&vms->grid_vel, &vms->basis_vel, &vms->local_vel, 
      GKYL_MODEL_SR, vms->vel_map, vms->mom_hamil, vms->gamma_inv, app->use_gpu);

    // Hamiltonian (GR) is a full phase-space array. 
    vms->hamil_range = vms->local; 
    vms->hamil = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
    vms->hamil_host = vms->hamil;
    if (app->use_gpu){
      vms->hamil_host = mkarr(false, vms->basis.num_basis, vms->local_ext.volume);
    }

    // Triad geometry context if needed for preset geometries.
    struct gkyl_triad_geom_ctx preset_geom_ctx = {
      .mass_bh = vms->geom->mass_bh,
      .spin_bh = vms->geom->spin_bh,
    };

    // Evaluate specified hamiltonian function at nodes to ensure continuity of
    // hamiltonian. On non-uniform meshes H(x,p) is sampled at the physical
    // phase-space coordinates via the position/velocity maps (nodes shared by
    // adjacent cells map to the same physical point, so continuity holds).
    struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
        .grid = &vms->grid, .basis = &vms->basis, .num_ret_vals = 1,
        .eval = gkyl_vlasov_triad_preset_hamil(cdim, vdim, vms->geom->triad_preset_geom_type),
        .ctx = &preset_geom_ctx,
        .c2p_func = vm_species_c2p_phase, .c2p_func_ctx = &c2p_ctx,
      }
    );
    gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &vms->local_ext, vms->hamil_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->hamil, vms->hamil_host);
    }
    gkyl_eval_on_nodes_release(hamil_proj);

    struct gkyl_vlasov_triad_geom_inp inp_triad_geom = { 0 };
    inp_triad_geom.use_vierbein = vms->info.use_vierbein;
    inp_triad_geom.use_preset_geom = vms->geom->use_preset_geom;
    if (vms->geom->use_preset_geom) {
      inp_triad_geom.triad_preset_geom_type = vms->geom->triad_preset_geom_type;
      inp_triad_geom.eval_vierbein_ctx = &preset_geom_ctx;
      inp_triad_geom.eval_vierbein_gradient_ctx = &preset_geom_ctx;
    }
    else {
      // Only Preset Geometries are allowed for GKYL_MODEL_TRIAD_GR
      assert(false);
    }

    // Sample the triad geometry at physical conf coordinates on non-uniform meshes.
    inp_triad_geom.c2p_func = vm_species_c2p_conf;
    inp_triad_geom.c2p_func_ctx = &c2p_ctx;

    // The geometry comes from the tangents and triads
    gkyl_vlasov_triad_geom_new(&app->grid, &app->local, app->basis,
      &vms->grid, &vms->local, vms->basis, inp_triad_geom, vms->conf_poisson_tensor_host);

    // Copy Pi_conf onto the device.
    if (app->use_gpu) {
      gkyl_array_copy(vms->conf_poisson_tensor, vms->conf_poisson_tensor_host); 
    }

  }
  else {
    // Canonical-PB models carry a general phase-space Hamiltonian.
    vms->hamil_id = GKYL_HAMIL_PHASE;
    vms->mom_hamil_id = GKYL_HAMIL_PHASE;

    // Canonical bracket = identity Poisson tensor: streaming through the
    // configuration-space flux machinery uses alpha_dir = P . grad_v H, which
    // reduces to dH/dv_dir with P[dir][j] = delta_dir,j. Set the cell-constant
    // coefficient of each diagonal block to the modal representation of 1.
    gkyl_array_clear(vms->conf_poisson_tensor, 0.0);
    for (int d = 0; d < app->cdim; ++d) {
      gkyl_array_shiftc(vms->conf_poisson_tensor, pow(sqrt(2.0), app->cdim),
        app->basis.num_basis*(d*vdim + d));
    }

    // Hamiltonian is a full phase-space array.
    vms->hamil_range = vms->local;
    vms->hamil = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
    vms->hamil_host = vms->hamil;
    if (app->use_gpu){
      vms->hamil_host = mkarr(false, vms->basis.num_basis, vms->local_ext.volume);
    }

    // Allocate arrays for specified metric inverse
    vms->h_ij = mkarr(app->use_gpu, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    vms->h_ij_host = vms->h_ij;
    if (app->use_gpu){
      vms->h_ij_host = mkarr(false, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    }

    // Allocate arrays for specified metric inverse
    vms->h_ij_inv = mkarr(app->use_gpu, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    vms->h_ij_inv_host = vms->h_ij_inv;
    if (app->use_gpu){
      vms->h_ij_inv_host = mkarr(false, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    }

    // Allocate arrays for specified metric determinant
    vms->det_h = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    vms->det_h_host = vms->det_h;
    if (app->use_gpu){
      vms->det_h_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    // Determine if we are using the extended Hamiltonian defintion
    if (vms->info.use_extended_hamil_def) {
      // Allocate arrays for background flows
      vms->background_flows = mkarr(app->use_gpu, app->basis.num_basis*vdim, app->local_ext.volume);
      vms->background_flows_host = vms->background_flows;
      if (app->use_gpu){
        vms->background_flows_host = mkarr(false, app->basis.num_basis*vdim, app->local_ext.volume);
      }

      // Allocate arrays for the effective potential
      vms->effective_potential = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      vms->effective_potential_host = vms->effective_potential;
      if (app->use_gpu){
        vms->effective_potential_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      }
    }

    // Evaluate specified hamiltonian function at nodes to ensure continuity of
    // hamiltonian, sampling at physical phase-space coordinates on mapped meshes.
    struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
        .grid = &vms->grid, .basis = &vms->basis, .num_ret_vals = 1,
        .eval = vms->info.hamil, .ctx = vms->info.hamil_ctx,
        .c2p_func = vm_species_c2p_phase, .c2p_func_ctx = &c2p_ctx,
      }
    );
    gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &vms->local_ext, vms->hamil_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->hamil, vms->hamil_host);
    }
    gkyl_eval_on_nodes_release(hamil_proj);

    // Evaluate specified metric function at nodes to ensure continuity
    struct gkyl_eval_on_nodes* h_ij_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
        .grid = &app->grid, .basis = &app->basis, .num_ret_vals = vdim*(vdim+1)/2,
        .eval = vms->info.h_ij, .ctx = vms->info.h_ij_ctx,
        .c2p_func = vm_species_c2p_conf, .c2p_func_ctx = &c2p_ctx,
      }
    );
    gkyl_eval_on_nodes_advance(h_ij_proj, 0.0, &app->local, vms->h_ij_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->h_ij, vms->h_ij_host);
    }
    gkyl_eval_on_nodes_release(h_ij_proj);

    // Evaluate specified inverse metric function at nodes to ensure continuity of the inverse
    struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
        .grid = &app->grid, .basis = &app->basis, .num_ret_vals = vdim*(vdim+1)/2,
        .eval = vms->info.h_ij_inv, .ctx = vms->info.h_ij_inv_ctx,
        .c2p_func = vm_species_c2p_conf, .c2p_func_ctx = &c2p_ctx,
      }
    );
    gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &app->local, vms->h_ij_inv_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->h_ij_inv, vms->h_ij_inv_host);
    }
    gkyl_eval_on_nodes_release(h_ij_inv_proj);

    // Evaluate specified determinant metric function at nodes to ensure continuity of the determinant
    struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
        .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 1,
        .eval = vms->info.det_h, .ctx = vms->info.det_h_ctx,
        .c2p_func = vm_species_c2p_conf, .c2p_func_ctx = &c2p_ctx,
      }
    );
    gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &app->local, vms->det_h_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->det_h, vms->det_h_host);
    }
    gkyl_eval_on_nodes_release(det_h_proj);    

    // As the mom_hamil and hamil are identical, aquire hamil and hamil_range here
    vms->mom_hamil_range = vms->hamil_range;
    vms->mom_hamil = gkyl_array_acquire(vms->hamil);

    // Evaluate specified determinant metric function at nodes to ensure continuity of the background flows
    if (vms->info.use_extended_hamil_def) {
      struct gkyl_eval_on_nodes* background_flows_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
          .grid = &app->grid, .basis = &app->basis, .num_ret_vals = vdim,
          .eval = vms->info.background_flows, .ctx = vms->info.background_flows_ctx,
          .c2p_func = vm_species_c2p_conf, .c2p_func_ctx = &c2p_ctx,
        }
      );
      gkyl_eval_on_nodes_advance(background_flows_proj, 0.0, &app->local, vms->background_flows_host);
      if (app->use_gpu){
        gkyl_array_copy(vms->background_flows, vms->background_flows_host);
      }
      gkyl_eval_on_nodes_release(background_flows_proj);    

      // Evaluate specified determinant metric function at nodes to ensure continuity of the effective potential
      struct gkyl_eval_on_nodes* effective_potential_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp) {
          .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 1,
          .eval = vms->info.effective_potential, .ctx = vms->info.effective_potential_ctx,
          .c2p_func = vm_species_c2p_conf, .c2p_func_ctx = &c2p_ctx,
        }
      );
      gkyl_eval_on_nodes_advance(effective_potential_proj, 0.0, &app->local, vms->effective_potential_host);
      if (app->use_gpu){
        gkyl_array_copy(vms->effective_potential, vms->effective_potential_host);
      }
      gkyl_eval_on_nodes_release(effective_potential_proj); 
    }   
  }
}

static void 
vm_species_new_radiation(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int vdim = app->vdim;    
  vms->has_rad = false;
  enum gkyl_vlasov_radiation_id radiation_id = vms->info.radiation.radiation_id;
  vms->rad = mkarr(app->use_gpu, vdim*vms->basis_vel.num_basis, vms->local_vel.volume);
  if (radiation_id ==  GKYL_VM_COMPTON_RADIATION || radiation_id == GKYL_VM_CURVATURE_RADIATION) {
    vms->has_rad = true;
    gkyl_dg_vlasov_calc_radiation(&vms->grid_vel, &vms->basis_vel, &vms->local_vel,
      radiation_id, vms->vel_map, vms->info.radiation.t_cool, vms->info.radiation.p0,
      vms->rad, app->use_gpu);
  }
}

static void
vm_species_apply_ic_dynamic(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  if (vms->num_init > 1) {
    gkyl_array_clear(vms->f, 0.0);  
    for (int k=0; k<vms->num_init; k++) {
      vm_species_projection_calc(app, vms, &vms->proj_init[k], vms->f1, tm);
      gkyl_array_accumulate(vms->f, 1.0, vms->f1);
    }
  }
  else {
    vm_species_projection_calc(app, vms, &vms->proj_init[0], vms->f, tm);
  }

  // Pre-compute applied acceleration in case it's time-independent
  vm_species_collisionless_app_accel(app, &vms->collisionless, tm);

  // we are pre-computing source for now as it is time-independent
  vm_species_source_calc(app, vms, &vms->src, tm);

  if (vms->calc_bflux) {
    vm_species_bflux_rhs(app, vms, &vms->bflux, vms->f, vms->f1);
  }

  // Optional runtime configuration to use BGK collisions but with fixed input 
  // temperature relaxation based on the initial temperature value. 
  if (vms->bgk.fixed_temp_relax) {
    vm_species_bgk_moms_fixed_temp(app, vms, &vms->bgk, vms->f);
  }
  
  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  if (vms->lower_bc[0].type == GKYL_SPECIES_FIXED_FUNC) {
    gkyl_bc_basic_buffer_fixed_func(vms->bc_lo[0], vms->bc_buffer_lo_fixed, vms->f);
  }
  if (vms->upper_bc[0].type == GKYL_SPECIES_FIXED_FUNC) {
    gkyl_bc_basic_buffer_fixed_func(vms->bc_up[0], vms->bc_buffer_up_fixed, vms->f);  
  }
}

static void
vm_species_apply_ic_static(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  if (vms->num_init > 1) {
    gkyl_array_clear(vms->f, 0.0);  
    struct gkyl_array *f_tmp = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
    for (int k=0; k<vms->num_init; k++) {
      vm_species_projection_calc(app, vms, &vms->proj_init[k], f_tmp, tm);
      gkyl_array_accumulate(vms->f, 1.0, f_tmp);
    }
    gkyl_array_release(f_tmp);
  }
  else {
    vm_species_projection_calc(app, vms, &vms->proj_init[0], vms->f, tm);
  }  
}

static double
vm_species_rhs_dynamic(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  gkyl_array_clear(vms->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Update RHS due to collisionless terms. 
  vm_species_collisionless_rhs(app, vms, &vms->collisionless, fin, em, rhs);

  // Update RHS due to collisions if collisions are present. 
  vm_species_lbo_rhs(app, vms, &vms->lbo, fin, rhs);
  vm_species_bgk_rhs(app, vms, &vms->bgk, fin, rhs);

  if (vms->calc_bflux) {
    vm_species_bflux_rhs(app, vms, &vms->bflux, fin, rhs);
  }

  // Reduce the CFL frequency anc compute stable dt needed by this species.
  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(vms->omega_cfl, vms->cflrate, GKYL_MAX, &vms->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, vms->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = vms->omega_cfl[0];
  }
  double dt_out = app->cfl/omega_cfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return dt_out;
}

double
vm_species_rhs_implicit_dynamic(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  gkyl_array_clear(vms->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute implicit update and update rhs to new time with time step size dt. 
  vm_species_bgk_rhs_implicit(app, vms, &vms->bgk, fin, dt, rhs);
  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);

  if (vms->calc_bflux) {
    vm_species_bflux_rhs(app, vms, &vms->bflux, fin, rhs);
  }

  // Implicit step does not affect the stable time step. 
  double dt_out = DBL_MAX;  
  return dt_out;
}

static double
vm_species_rhs_static(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  double dt_out = DBL_MAX;  
  return dt_out;
}

static double
vm_species_rhs_implicit_static(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  double dt_out = DBL_MAX;  
  return dt_out;
}

static void
vm_species_apply_bc_dynamic(gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f,
  double tcurr)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(vms->comm, &vms->local, &vms->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 
  
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (vms->lower_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(app, vms, &vms->bc_emission_lo, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(vms->bc_lo[d], vms->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(vms->bc_lo[d], vms->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (vms->upper_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(app, vms, &vms->bc_emission_up, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(vms->bc_up[d], vms->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(vms->bc_up[d], vms->bc_buffer_up_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(vms->comm, &vms->local, &vms->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_apply_bc_static(gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f, 
  double tcurr)
{
  // do nothing
}

static void
vm_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

static void
vm_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // do nothing
}

static void
vm_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

static void
vm_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // do nothing
}

static void
vm_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

static void
vm_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // do nothing
}

static void
vm_species_write_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  vms->write_cfl_func(app, vms, tm, frame);
  vms->write_cell_avg_func(app, vms, tm, frame);
  vms->write_lte_func(app, vms, tm, frame);

  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = vms->basis.id
    }
  );
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);
  
  // Divide out the velocity space Jacobian if present
  // We do the division before I/O to increase the accuracy since we know
  // the velocity-space Jacobian at specific quadrature points. 
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, vms->f, vms->f_no_J);
  // Also divide out the (per-conf-cell constant) configuration-space Jacobian so
  // the written distribution is the physical f (J_x*J_v*f -> f).
  gkyl_vlasov_position_map_divide_jacobpos(vms->pos_map, &vms->basis,
    &vms->local, vms->f_no_J, vms->f_no_J);

  // If we are on device, copy the distribution function without the velocity-space
  // Jacobian to the host, otherwise just write out the f_no_J array.
  if (app->use_gpu) {
    gkyl_array_copy(vms->f_host, vms->f_no_J);
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, vms->f_host, fileNm);
  }
  else {
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, vms->f_no_J, fileNm);
  }
    
  vlasov_array_meta_release(mt);  

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_cfl_enabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = vms->basis.id,
    }
  );

  const char *fmt = "%s-%s-cflrate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);
  gkyl_array_copy(vms->cflrate_host, vms->cflrate);
  gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt,
    vms->cflrate_host, fileNm);

  vlasov_array_meta_release(mt);  

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_cfl_disabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_cell_avg_enabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = vms->basis.id,
    }
  );

  const char *fmt = "%s-%s_avg_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);

  // Divide out the velocity space Jacobian if present
  // We do the division before I/O to increase the accuracy since we know
  // the velocity-space Jacobian at specific quadrature points. 
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, vms->f, vms->f_no_J);
  gkyl_vlasov_position_map_divide_jacobpos(vms->pos_map, &vms->basis,
    &vms->local, vms->f_no_J, vms->f_no_J);

  // Copy the cell average into a temporary array and re-scale
  gkyl_array_set_offset(vms->cflrate, 1.0/pow(2.0, (app->cdim+app->vdim)/2.0), vms->f_no_J, 0); 
  gkyl_array_copy(vms->cflrate_host, vms->cflrate);
  gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt,
    vms->cflrate_host, fileNm);

  vlasov_array_meta_release(mt);  

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_cell_avg_disabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_lte_enabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = vms->basis.id
    }
  );

  const char *fmt = "%s-%s_%d_lte.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, frame);

  vm_species_lte(app, vms, &vms->lte, vms->f);
  
  // Divide out the velocity space Jacobian from LTE distribution if present
  // We do the division before I/O to increase the accuracy since we know
  // the velocity-space Jacobian at specific quadrature points. 
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, vms->lte.f_lte, vms->f_no_J);
  gkyl_vlasov_position_map_divide_jacobpos(vms->pos_map, &vms->basis,
    &vms->local, vms->f_no_J, vms->f_no_J);

  // If we are on device, copy the LTE distribution function without the velocity-space
  // Jacobian to the host, otherwise just write out the f_no_J array. 
  if (app->use_gpu) {
    gkyl_array_copy(vms->f_host, vms->f_no_J);
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, vms->f_host, fileNm);
  }
  else {
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, vms->f_no_J, fileNm);
  }

  vlasov_array_meta_release(mt);  

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_lte_disabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_mom_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  for (int m=0; m<vms->info.num_diag_moments; ++m) {
    struct timespec wtm = gkyl_wall_clock();
    vm_species_moment_calc(&vms->moms[m], vms->local, app->local, vms->f);
    app->stat.n_mom += 1;
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wtm);
      
    struct timespec wst = gkyl_wall_clock();
    if (app->use_gpu) {
      gkyl_array_copy(vms->moms[m].marr_host, vms->moms[m].marr);
    }

    // Moments of the stored J_x J_v f carry the configuration-space Jacobian J
    // (it factors out of the velocity integral). Divide it out so the written
    // moment field is physical. The map is C^0 linear, so this is an exact
    // per-cell scalar division. For the LTE moment diagnostic (n, V_drift, T/m)
    // only the density n (the first num_basis coefficients) carries J; V_drift
    // and T/m are velocity ratios in which J cancels, so they are left alone.
    // Pure moments (M0, M1i, M2, ...) carry J in every component. (Identity
    // map: skip, leaving output unchanged.)
    if (!vms->pos_map->is_identity) {
      int num_coeff_divide = vms->moms[m].is_vlasov_lte_moms ?
        app->basis.num_basis : vms->moms[m].marr_host->ncomp;
      gkyl_vlasov_position_map_divide_jacobpos_conf(vms->pos_map, &app->local,
        num_coeff_divide, vms->moms[m].marr_host, vms->moms[m].marr_host);
    }

    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->name,
      gkyl_distribution_moments_strs[vms->info.diag_moments[m]], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name,
      gkyl_distribution_moments_strs[vms->info.diag_moments[m]], frame);
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      vms->moms[m].marr_host, fileNm);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;
  }
  
  vlasov_array_meta_release(mt);

  // Collision diagnostics (nu_sum, nu_prim_moms), no-ops unless the species
  // has collisions with write_coll_diagnostics set.
  vm_species_lbo_write_mom(app, vms, tm, frame);
  vm_species_bgk_write_mom(app, vms, tm, frame);

  app->stat.n_diag += 1;
}

static void
vm_species_write_mom_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_calc_integrated_mom_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  int vdim = app->vdim;
  double avals[2+vdim], avals_global[2+vdim];
  
  vm_species_moment_calc(&vms->integ_moms, vms->local, app->local, vms->f);
  app->stat.n_mom += 1;

  // reduce to compute sum over whole domain, append to diagnostics
  if (app->use_gpu) {
    gkyl_array_reduce_range(vms->red_integ_diag, vms->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals, vms->red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    gkyl_array_reduce_range(avals, vms->integ_moms.marr_host, GKYL_SUM, &app->local);
  }

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
  gkyl_dynvec_append(vms->integ_diag, tm, avals_global);  

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_calc_integrated_mom_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  // do nothing
}

static void
vm_species_write_integrated_mom_dynamic(gkyl_vlasov_app *app, struct vm_species *vms)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->name, "imom");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, "imom");
    
    if (vms->is_first_integ_write_call) {
      gkyl_dynvec_write(vms->integ_diag, fileNm);
      vms->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(vms->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(vms->integ_diag);
  app->stat.n_diag_io += 1;
  
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_write_integrated_mom_static(gkyl_vlasov_app* app, struct vm_species *vms)
{
  // do nothing
}

static void
vm_species_calc_L2_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  // L^2 energy with nonuniform meshes is (J_x J_v f)*f. Dividing both the
  // velocity- and configuration-space Jacobians out of one factor leaves the
  // physical f; multiplied by the stored J_x J_v f this gives the Jacobian-
  // weighted f^2 integrated over the computational grid.
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, vms->f, vms->f_no_J);
  gkyl_vlasov_position_map_divide_jacobpos(vms->pos_map, &vms->basis,
    &vms->local, vms->f_no_J, vms->f_no_J);
  gkyl_dg_calc_prod_op_range(vms->basis, 0, vms->L2_f, 0, vms->f_no_J, 0, vms->f, vms->local);
  gkyl_array_scale_range(vms->L2_f, vms->grid.cellVolume, &vms->local); 
  
  double L2[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(vms->red_L2_f, vms->L2_f, GKYL_SUM, &vms->local);
    gkyl_cu_memcpy(L2, vms->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(L2, vms->L2_f, GKYL_SUM, &vms->local);
  }
  double L2_global[1] = { 0.0 };
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);
  
  gkyl_dynvec_append(vms->integ_L2_f, tm, L2_global);  

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

static void
vm_species_calc_L2_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  // do nothing
}

static void
vm_species_write_L2_dynamic(gkyl_vlasov_app* app, struct vm_species *vms)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write the L2 norm.
    const char *fmt = "%s-%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->name, "L2");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->name, "L2");

    if (vms->is_first_integ_L2_write_call) {
      gkyl_dynvec_write(vms->integ_L2_f, fileNm);
      vms->is_first_integ_L2_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(vms->integ_L2_f, fileNm);
    }
  }
  gkyl_dynvec_clear(vms->integ_L2_f);

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

static void
vm_species_write_L2_static(gkyl_vlasov_app* app, struct vm_species *vms)
{
  // do nothing
}

static void
vm_species_release_dynamic(const gkyl_vlasov_app* app, const struct vm_species *vms)
{
  // Release various arrays used in time stepping and boundary conditions. 
  gkyl_array_release(vms->f1);
  gkyl_array_release(vms->fnew);
  gkyl_array_release(vms->cflrate);
  gkyl_array_release(vms->cflrate_host);
  gkyl_array_release(vms->bc_buffer);
  gkyl_array_release(vms->bc_buffer_lo_fixed);
  gkyl_array_release(vms->bc_buffer_up_fixed);

  vm_species_moment_release(app, &vms->integ_moms); 

  gkyl_array_release(vms->L2_f);
  gkyl_dynvec_release(vms->integ_L2_f);
  gkyl_dynvec_release(vms->integ_diag);

  if (vms->source_id) {
    vm_species_source_release(app, &vms->src);
  }

  if (vms->lbo.collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_release(app, &vms->lbo);
  }
  else if (vms->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_release(app, &vms->bgk);
  }

  if (vms->calc_bflux) {
    vm_species_bflux_release(app, &vms->bflux);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (vms->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      vm_species_emission_release(&vms->bc_emission_lo);
    }
    else {
      gkyl_bc_basic_release(vms->bc_lo[d]);
    }
    
    if (vms->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      vm_species_emission_release(&vms->bc_emission_up);
    }
    else {
      gkyl_bc_basic_release(vms->bc_up[d]);
    }
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(vms->omega_cfl);
    gkyl_cu_free(vms->red_L2_f);
    gkyl_cu_free(vms->red_integ_diag);
  }
  else {
    gkyl_free(vms->omega_cfl);
  }  
}

static void
vm_species_release_static(const gkyl_vlasov_app* app, const struct vm_species *vms)
{
}

// Initialize dynamic species object.
static void
vm_species_new_dynamic(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  // Allocate remaining distribution function arrays for RK time stepping. 
  vms->f1 = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
  vms->fnew = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);

  // Allocate cflrate (scalar array) and maximum frequency for stable time step.
  vms->cflrate = mkarr(app->use_gpu, 1, vms->local_ext.volume);
  // Host-side cflrate for I/O on GPUs.
  vms->cflrate_host = app->use_gpu ? mkarr(false, vms->cflrate->ncomp, vms->cflrate->size)
                                   : gkyl_array_acquire(vms->cflrate);    
  if (app->use_gpu) {
    vms->omega_cfl = gkyl_cu_malloc(sizeof(double));
  }
  else {
    vms->omega_cfl = gkyl_malloc(sizeof(double));
  }

  // Allocate data for integrated moments. 
  vm_species_moment_init(app, vms, &vms->integ_moms, GKYL_F_MOMENT_M0M1M2, true);  
  // Array for storing f^2 in each cell.
  vms->L2_f = mkarr(app->use_gpu, 1, vms->local_ext.volume);
  if (app->use_gpu) {
    vms->red_L2_f = gkyl_cu_malloc(sizeof(double));
    vms->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments and f^2.
  vms->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  vms->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  vms->is_first_integ_L2_write_call = true;
  vms->is_first_integ_write_call = true; 
  
  // Allocate buffer for applying BCs.
  long buff_sz = 0;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(vms->lower_skin[dir].volume, vms->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  vms->bc_buffer = mkarr(app->use_gpu, vms->basis.num_basis, buff_sz);
  // Buffers for fixed function BCs on distribution function.
  vms->bc_buffer_lo_fixed = mkarr(app->use_gpu, vms->basis.num_basis, buff_sz);
  vms->bc_buffer_up_fixed = mkarr(app->use_gpu, vms->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (vms->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      vms->emit_lo = true;
      vms->calc_bflux = true;
      vm_species_emission_init(app, &vms->bc_emission_lo, d, GKYL_LOWER_EDGE, vms->lower_bc[d].aux_ctx);
    }
    else {
      if (vms->lower_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (vms->lower_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (vms->lower_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_DISTF_REFLECT;
      else if (vms->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      vms->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, vms->basis_on_dev,
        &vms->lower_skin[d], &vms->lower_ghost[d], vms->f->ncomp, cdim, app->use_gpu);
    }

    // Upper BC updater. Copy BCs by default.
    if (vms->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      vms->emit_up = true;
      vms->calc_bflux = true;
      vm_species_emission_init(app, &vms->bc_emission_up, d, GKYL_UPPER_EDGE, vms->upper_bc[d].aux_ctx);
    }
    else {
      if (vms->upper_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (vms->upper_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (vms->upper_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_DISTF_REFLECT;
      else if (vms->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      vms->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, vms->basis_on_dev,
        &vms->upper_skin[d], &vms->upper_ghost[d], vms->f->ncomp, cdim, app->use_gpu);
    }
  }

  // Intitalize boundary flux updater if we need boundary fluxes.
  if (vms->calc_bflux) {
    vm_species_bflux_init(app, vms, &vms->bflux);
  }

  // Set function pointers.
  vms->apply_ic_func = vm_species_apply_ic_dynamic; 
  vms->rhs_func = vm_species_rhs_dynamic;
  vms->rhs_implicit_func = vm_species_rhs_implicit_dynamic;
  vms->bc_func = vm_species_apply_bc_dynamic;
  vms->release_func = vm_species_release_dynamic;
  vms->step_f_func = vm_species_step_f_dynamic;
  vms->combine_func = vm_species_combine_dynamic;
  vms->copy_func = vm_species_copy_range_dynamic;
  vms->write_func = vm_species_write_dynamic;
  if (vms->info.correct.output_f_lte) {
    vms->write_lte_func = vm_species_write_lte_enabled;
  }
  else {
    vms->write_lte_func = vm_species_write_lte_disabled;
  }
  if (vms->info.write_omega_cfl) {
    vms->write_cfl_func = vm_species_write_cfl_enabled;
  }
  else {
    vms->write_cfl_func = vm_species_write_cfl_disabled;
  }
  if (vms->info.write_cell_avg) {
    vms->write_cell_avg_func = vm_species_write_cell_avg_enabled;
  }
  else {
    vms->write_cell_avg_func = vm_species_write_cell_avg_disabled;
  }
  vms->write_mom_func = vm_species_write_mom_dynamic;
  vms->calc_integrated_mom_func = vm_species_calc_integrated_mom_dynamic;
  vms->write_integrated_mom_func = vm_species_write_integrated_mom_dynamic;
  vms->calc_L2_func = vm_species_calc_L2_dynamic;
  vms->write_L2_func = vm_species_write_L2_dynamic;
}

// Initialize static species object.
static void
vm_species_new_static(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  // When using static species, only need one allocation.
  vms->f1 = vms->f;
  vms->fnew = vms->f;
  
  // Set function pointers.
  vms->apply_ic_func = vm_species_apply_ic_static; 
  vms->rhs_func = vm_species_rhs_static;
  vms->rhs_implicit_func = vm_species_rhs_implicit_static;
  vms->bc_func = vm_species_apply_bc_static;
  vms->release_func = vm_species_release_static;
  vms->step_f_func = vm_species_step_f_static;
  vms->combine_func = vm_species_combine_static;
  vms->copy_func = vm_species_copy_range_static;
  vms->write_func = vm_species_write_static;
  vms->write_lte_func = vm_species_write_lte_disabled;
  vms->write_cfl_func = vm_species_write_cfl_disabled;
  vms->write_cell_avg_func = vm_species_write_cell_avg_disabled;
  vms->write_mom_func = vm_species_write_mom_static;
  vms->calc_integrated_mom_func = vm_species_calc_integrated_mom_static;
  vms->write_integrated_mom_func = vm_species_write_integrated_mom_static;
  vms->calc_L2_func = vm_species_calc_L2_static;
  vms->write_L2_func = vm_species_write_L2_static;
}

// End static function definitions.

// Initialize species object.
// The 2x3v and 3x3v tensor p=1 hybrid kernels are optional build sets (see
// ./configure --help and the top-level Makefile). Refuse a configuration whose
// kernels were not built with an actionable message instead of a NULL kernel
// pointer assert deep in a constructor.
static void
vm_species_check_hyb_build(struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  if (vms->basis.b_type != GKYL_BASIS_MODAL_HYBRID) return;
  int cdim = app->cdim, vdim = vms->basis.ndim - cdim;
#ifndef GKYL_BUILD_VLASOV_HYB_2X3V
  if (cdim == 2 && vdim == 3)
    gkyl_exit("vm_species: the 2x3v tensor p=1 hybrid kernels were not built. Reconfigure with --build-vlasov-hyb-2x3v=yes.");
#endif
#ifndef GKYL_BUILD_VLASOV_HYB_3X3V
  if (cdim == 3 && vdim == 3)
    gkyl_exit("vm_species: the 3x3v tensor p=1 hybrid kernels were not built. Reconfigure with --build-vlasov-hyb-3x3v=yes.");
#endif
#ifndef GKYL_BUILD_VLASOV_HYB_3X3V_PHASE
  if (cdim == 3 && vdim == 3 && vms->hamil_id == GKYL_HAMIL_PHASE)
    gkyl_exit("vm_species: the 3x3v tensor p=1 hybrid phase-space Hamiltonian kernels were not built. Reconfigure with --build-vlasov-hyb-3x3v-phase=yes (and --build-vlasov-hyb-3x3v=yes).");
#endif
}

void
vm_species_init(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm_app_inp->cells[d];
    lower[d] = vm_app_inp->lower[d];
    upper[d] = vm_app_inp->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // Full phase-space grid.
    cells[cdim+d] = vms->info.cells[d];
    lower[cdim+d] = vms->info.lower[d];
    upper[cdim+d] = vms->info.upper[d];
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.

    // Only velocity space.
    cells_vel[d] = vms->info.cells[d];
    lower_vel[d] = vms->info.lower[d];
    upper_vel[d] = vms->info.upper[d];
    ghost_vel[d] = 0; // No ghost-cells in velocity space.
  }

  // Allocate device basis if we are using GPUs.
  if (app->use_gpu) {
    vms->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    vms->basis_on_dev = &vms->basis;
  }  

  // Determine basis type from configuration-space basis and create phase-space basis. 
  enum gkyl_basis_type b_type = app->basis.b_type;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&vms->basis, pdim, app->poly_order);
      gkyl_cart_modal_serendip(&vms->basis_surf, pdim-1, app->poly_order);
      gkyl_cart_modal_serendip(&vms->basis_vel, vdim, app->poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(vms->basis_on_dev, pdim, app->poly_order);
      }
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (app->poly_order == 1) {
        // Tensor p=1 is the tensor hybrid basis: p=1 in configuration space
        // tensored with p=2 in velocity space (the pure p=1 tensor phase
        // basis has no kernels). The velocity basis is the p=2 tensor basis,
        // which also selects the C^1 cubic velocity map and the p=2 velocity
        // representation of the Hamiltonian.
        gkyl_cart_modal_hybrid(&vms->basis, cdim, vdim);
        gkyl_cart_modal_tensor(&vms->basis_surf, pdim-1, app->poly_order);
        gkyl_cart_modal_tensor(&vms->basis_vel, vdim, 2);
        if (app->use_gpu) {
          gkyl_cart_modal_hybrid_cu_dev(vms->basis_on_dev, cdim, vdim);
        }
      }
      else {
        gkyl_cart_modal_tensor(&vms->basis, pdim, app->poly_order);
        gkyl_cart_modal_tensor(&vms->basis_surf, pdim-1, app->poly_order);
        gkyl_cart_modal_tensor(&vms->basis_vel, vdim, app->poly_order);
        if (app->use_gpu) {
          gkyl_cart_modal_tensor_cu_dev(vms->basis_on_dev, pdim, app->poly_order);
        }
      }
      break;
    default:
      assert(false);
      break;
  }  

  // Full phase-space grid.
  gkyl_rect_grid_init(&vms->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vms->grid, ghost, &vms->global_ext, &vms->global);
  
  // Velocity-space grid.
  gkyl_rect_grid_init(&vms->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&vms->grid_vel, ghost_vel, &vms->local_ext_vel, &vms->local_vel);

  // Phase-space communicator.
  vms->comm = gkyl_comm_extend_comm(app->comm, &vms->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // local = conf-local X local_vel
  gkyl_range_ten_prod(&local, &app->local, &vms->local_vel);
  gkyl_create_ranges(&local, ghost, &vms->local_ext, &vms->local);

  // Determine which directions are not periodic.
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }
  // Set non-periodic boundary conditions. 
  for (int dir=0; dir<cdim; ++dir) {
    vms->lower_bc[dir].type = vms->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const struct gkyl_vlasov_bcs *bc;
      if (dir == 0) {
        bc = &vms->info.bcx;
      }
      else if (dir == 1) {
        bc = &vms->info.bcy;
      }
      else {
        bc = &vms->info.bcz;
      }

      vms->lower_bc[dir] = bc->lower;
      vms->upper_bc[dir] = bc->upper;
    }
  }

  // Store model type.
  vms->model_id = vms->info.model_id; 
  // Store field type (a field object always exists; GKYL_FIELD_NULL if none).
  vms->field_id = app->field->field_id;
  
  // Allocate distribution function array.
  vms->f = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
  // Host-side distribution function for I/O on GPUs.
  vms->f_host = app->use_gpu ? mkarr(false, vms->f->ncomp, vms->f->size)
                             : gkyl_array_acquire(vms->f);  

  vms->write_cell_avg = vms->info.write_cell_avg; // Write out only the cell averages?

  // Velocity-space mapping object. Always created: directions without a
  // user-specified mapping get the identity map, so uniform and non-uniform
  // velocity grids are handled transparently downstream.
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  for (int v=0; v<vdim; ++v) {
    if (vms->info.mapc2p_vel[v].mapc2p_vel_func) {
      inp_vmap[v].eval_vmap = vms->info.mapc2p_vel[v].mapc2p_vel_func;
      inp_vmap[v].ctx = vms->info.mapc2p_vel[v].mapc2p_vel_ctx;
    }
  }
  vms->vel_map = gkyl_vlasov_velocity_map_new(&vms->grid_vel, &vms->local_vel,
    &vms->basis_vel, inp_vmap, vms->info.use_lo, app->use_gpu);

  // The velocity map is static in time, so write it (uniform grids included)
  // only with the first frame. Both the p=3 map and its p=0 cell average are
  // written, with metadata built from the map's I/O basis.
  gkyl_vlasov_velocity_map_write(vms->vel_map, vms->comm,
    app->name, vms->name);

  // Configuration-space mapping object. Created once on the app and shared by
  // all species (configuration space is common to every species); acquire a
  // reference here.
  vms->pos_map = gkyl_vlasov_position_map_acquire(app->pos_map);

  // Allocate array for dividing out velocity-space Jacobian. 
  // If the mesh is uniform, we simply copy the distribution function at that RK stage
  // into this array for use in the velocity-space surface flux computation and 
  // the magnetic field volume update if magnetic fields are present. 
  vms->f_no_J = mkarr(app->use_gpu, vms->f->ncomp, vms->f->size); 

  // Copy the geometry pointer
  vms->geom = app->vm_geom;

  // Construct Hamiltonian. 
  vm_species_new_hamil(vm_app_inp, app, vms); 
  vm_species_check_hyb_build(app, vms);

  // Determine whether we have radiation. 
  vm_species_new_radiation(vm_app_inp, app, vms); 

  // Initialize the collisionless solver.
  vms->collisionless = (struct vm_collisionless) { };
  vm_species_collisionless_init(app, vms, &vms->collisionless);

  // Allocate data for momentum (for use in current accumulation).
  vm_species_moment_init(app, vms, &vms->m1i, GKYL_F_MOMENT_M1_FROM_H, false);
  // Allocate date for density (for use in charge density accumulation and weak division for V_drift).
  vm_species_moment_init(app, vms, &vms->m0, GKYL_F_MOMENT_M0, false);

  // Allocate data for diagnostic moments.
  int ndm = vms->info.num_diag_moments;
  vms->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m) {
    vm_species_moment_init(app, vms, &vms->moms[m], vms->info.diag_moments[m], false); 
  }

  // Create skin/ghost ranges for applying BCs.
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&vms->lower_skin[dir], &vms->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &vms->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&vms->upper_skin[dir], &vms->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &vms->local_ext, ghost);
  }

  // Global skin and ghost ranges, only valid (i.e. volume>0) in ranges abutting boundaries.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&vms->global_lower_skin[dir], &vms->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE, &vms->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&vms->global_upper_skin[dir], &vms->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE, &vms->global_ext, ghost);

    gkyl_sub_range_intersect(&vms->global_lower_skin[dir], &vms->local_ext, &vms->global_lower_skin[dir]);
    gkyl_sub_range_intersect(&vms->global_upper_skin[dir], &vms->local_ext, &vms->global_upper_skin[dir]);

    gkyl_sub_range_intersect(&vms->global_lower_ghost[dir], &vms->local_ext, &vms->global_lower_ghost[dir]);
    gkyl_sub_range_intersect(&vms->global_upper_ghost[dir], &vms->local_ext, &vms->global_upper_ghost[dir]);
  }

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine.
  // Projection routine optionally corrects all the Maxwellian/LTE moments.
  // This routine is utilized by BGK collisions.
  vms->lte = (struct vm_lte) { };
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = vms->info.correct.correct_all_moms,
    .max_iter = vms->info.correct.max_iter > 0 ? vms->info.correct.max_iter : 100,
    .iter_eps = vms->info.correct.iter_eps > 0 ? vms->info.correct.iter_eps : 1e-12,
    .use_last_converged = vms->info.correct.use_last_converged,
  };
  vm_species_lte_init(app, vms, &vms->lte, corr_inp);

  // Initialize projection routine for initial conditions. 
  vms->num_init = vms->info.num_init;
  for (int k=0; k<vms->num_init; k++) {
    vm_species_projection_init(app, vms, vms->info.projection[k], &vms->proj_init[k]);
  }

  // Initialize empty structs. New methods will fill them if specified.
  vms->src = (struct vm_source) { };
  vms->bflux = (struct vm_boundary_fluxes) { };

  // Initialize collision objects. Init method checks if collision
  // table has been specified and allocates appropriate arrays and
  // sets relevant function pointers. 
  vms->lbo = (struct vm_lbo_collisions) { };
  vm_species_lbo_init(app, vms, &vms->lbo);
  vms->bgk = (struct vm_bgk_collisions) { };
  vm_species_bgk_init(app, vms, &vms->bgk);

  // Set species source id. 
  vms->source_id = vms->info.source.source_id;
  if (vms->source_id == GKYL_BFLUX_SOURCE) {
    vms->calc_bflux = true;
  }

  if (!vms->info.is_static) {
    vm_species_new_dynamic(vm_app_inp, app, vms);
  }
  else {
    vm_species_new_static(vm_app_inp, app, vms);
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  return vms->apply_ic_func(app, vms, tm); 
}

// Compute the RHS for species update, returning maximum stable time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  return vms->rhs_func(app, vms, fin, em, rhs);
}

// Compute the implicit RHS for species update; returns DBL_MAX because implicit
// update does not affect the stable time step.
double
vm_species_rhs_implicit(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  return vms->rhs_implicit_func(app, vms, fin, rhs, dt);
}

// Accummulate function for forward Euler method.
void
vm_species_step_f(struct vm_species *vms, struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  vms->step_f_func(out, a, inp);
}

// Combine function for SSP-RK3 updates.
void
vm_species_combine(struct vm_species *vms, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  vms->combine_func(out, c1, arr1, c2, arr2, rng);
}

// Copy function for SSP-RK3 updates.
void
vm_species_copy_range(struct vm_species *vms, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  vms->copy_func(out, inp, range);
}

// Apply boundary conditions to the distribution function.
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f, 
  double tcurr)
{
  vms->bc_func(app, vms, f, tcurr);
}

// Write distribution function. 
void
vm_species_write(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  if (frame == 0) {
    vm_species_write_dynamic(app, vms, tm, frame);
  }
  else {
    vms->write_func(app, vms, tm, frame);
  }
}

// Write distribution function moments. 
void
vm_species_write_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  if (frame == 0) {
    vm_species_write_mom_dynamic(app, vms, tm, frame);
  }
  else {
    vms->write_mom_func(app, vms, tm, frame);
  }
}

// Calculate integrated moments. 
void
vm_species_calc_integrated_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  vms->calc_integrated_mom_func(app, vms, tm);
}

// Write integrated moments. 
void
vm_species_write_integrated_mom(gkyl_vlasov_app* app, struct vm_species *vms)
{
  vms->write_integrated_mom_func(app, vms);
}

// Calculate integrated f^2. 
void
vm_species_calc_L2(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  vms->calc_L2_func(app, vms, tm);
}

// Write integrated f^2.
void
vm_species_write_L2(gkyl_vlasov_app* app, struct vm_species *vms)
{
  vms->write_L2_func(app, vms);
}

// Write the number of correction iterations to output statistics struct. 
void
vm_species_n_iter_corr(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (!app->species[i].kinetic) continue;
    app->stat.num_corr[i] = app->species[i].kinetic->lte.num_corr;
    app->stat.n_iter_corr[i] = app->species[i].kinetic->lte.n_iter;
  }
}

void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *vms)
{
  // Release resources for Vlasov species.
  gkyl_array_release(vms->f);
  gkyl_comm_release(vms->comm);
  gkyl_array_release(vms->f_host);
  if (app->use_gpu) {
    gkyl_cu_free(vms->basis_on_dev);
  }

  for (int k=0; k<vms->num_init; k++) {
    vm_species_projection_release(app, &vms->proj_init[k]);
  }

  gkyl_array_release(vms->rad); 
  gkyl_array_release(vms->f_no_J); 

  gkyl_vlasov_velocity_map_release(vms->vel_map);
  gkyl_vlasov_position_map_release(vms->pos_map);

  // Release arrays for different types of Vlasov equations.
  if (vms->model_id  == GKYL_MODEL_DEFAULT || vms->model_id  == GKYL_MODEL_SR 
      || vms->model_id  == GKYL_MODEL_TRIAD || vms->model_id  == GKYL_MODEL_TRIAD_GR) {
    if (vms->model_id == GKYL_MODEL_SR) {
      gkyl_dg_calc_sr_vars_release(vms->sr_vars);
    }
    gkyl_array_release(vms->gamma_inv);   
  }
  else  { 
    gkyl_array_release(vms->h_ij);
    gkyl_array_release(vms->h_ij_inv);
    gkyl_array_release(vms->det_h);
    if (vms->info.use_extended_hamil_def) {
      gkyl_array_release(vms->background_flows);
      gkyl_array_release(vms->effective_potential);
    }
    if (app->use_gpu){
      gkyl_array_release(vms->hamil_host);
      gkyl_array_release(vms->h_ij_host);
      gkyl_array_release(vms->h_ij_inv_host);
      gkyl_array_release(vms->det_h_host);
      if (vms->info.use_extended_hamil_def) {
        gkyl_array_release(vms->background_flows_host);
        gkyl_array_release(vms->effective_potential_host);
      }
    }
  }
  gkyl_array_release(vms->hamil);
  if (app->use_gpu && vms->model_id == GKYL_MODEL_TRIAD_GR) {
    gkyl_array_release(vms->hamil_host);
  }
  gkyl_array_release(vms->mom_hamil);

  gkyl_array_release(vms->conf_poisson_tensor);
  if (app->use_gpu) {
    gkyl_array_release(vms->conf_poisson_tensor_host);
  }

  // Release collisionless solver.
  vm_species_collisionless_release(app, vms, &vms->collisionless);

  // Release moment data.
  vm_species_moment_release(app, &vms->m1i);
  vm_species_moment_release(app, &vms->m0);
  for (int i=0; i<vms->info.num_diag_moments; ++i) {
    vm_species_moment_release(app, &vms->moms[i]);
  }
  gkyl_free(vms->moms);
  
  vm_species_lte_release(app, &vms->lte);

  vms->release_func(app, vms);
}
