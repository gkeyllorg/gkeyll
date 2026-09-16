#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Configuration-space c2p: map computational conf coords to physical via the
// position map. Used by conf-space projections (e.g. LTE moments).
static void
vm_proj_c2p_conf(const double *xcomp, double *xphys, void *ctx)
{
  struct vm_proj_c2p_ctx *c = ctx;
  gkyl_vlasov_position_map_eval_mc2p(c->pos_map, xcomp, xphys);
}

// Phase-space c2p: map conf coords via the position map and velocity coords via
// the velocity map (the DG maps the solver kernels assume). Used by the FUNC
// initial-condition projection.
static void
vm_proj_c2p_phase(const double *xcomp, double *xphys, void *ctx)
{
  struct vm_proj_c2p_ctx *c = ctx;
  gkyl_vlasov_position_map_eval_mc2p(c->pos_map, xcomp, xphys);
  gkyl_vlasov_velocity_map_eval_c2p(c->vel_map, &xcomp[c->cdim], &xphys[c->cdim]);
}

void
vm_species_projection_init(struct gkyl_vlasov_app *app, struct vm_species *vms,
  struct gkyl_vlasov_projection inp, struct vm_proj *proj)
{
  proj->proj_id = inp.proj_id;
  proj->model_id = vms->model_id;
  // Coordinate-map context for projecting on non-uniform meshes (identity maps
  // make the c2p the identity, so this is transparent for uniform grids).
  proj->c2p_ctx = (struct vm_proj_c2p_ctx) {
    .cdim = app->cdim,
    .pos_map = vms->pos_map,
    .vel_map = vms->vel_map,
  };
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &vms->grid,
        .basis = &vms->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = vms->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.func,
        .ctx = inp.ctx_func,
        // Sample the IC at physical phase-space coordinates on non-uniform meshes.
        .c2p_func = vm_proj_c2p_phase,
        .c2p_func_ctx = &proj->c2p_ctx,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, vms->basis.num_basis, vms->local_ext.volume);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_VLASOV_LTE) {
    int vdim = app->vdim; 
    proj->dens = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->V_drift = mkarr(false, vdim*app->basis.num_basis, app->local_ext.volume);
    proj->T_over_m = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->vlasov_lte_moms_host = mkarr(false, (vdim+2)*app->basis.num_basis, app->local_ext.volume);

    // Sample the LTE moments at physical configuration coordinates on
    // non-uniform conf meshes (the conf c2p uses the position map).
    proj->proj_dens = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid, .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD, .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1, .eval = inp.density, .ctx = inp.ctx_density,
        .c2p_func = vm_proj_c2p_conf, .c2p_func_ctx = &proj->c2p_ctx,
      }
    );
    proj->proj_V_drift = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid, .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD, .num_quad = app->basis.poly_order+1,
        .num_ret_vals = vdim, .eval = inp.V_drift, .ctx = inp.ctx_V_drift,
        .c2p_func = vm_proj_c2p_conf, .c2p_func_ctx = &proj->c2p_ctx,
      }
    );
    proj->proj_temp = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid, .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD, .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1, .eval = inp.temp, .ctx = inp.ctx_temp,
        .c2p_func = vm_proj_c2p_conf, .c2p_func_ctx = &proj->c2p_ctx,
      }
    );

    proj->vlasov_lte_moms = mkarr(app->use_gpu, (vdim+2)*app->basis.num_basis, app->local_ext.volume);

    struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
      .phase_grid = &vms->grid,
      .vel_grid = &vms->grid_vel, 
      .conf_basis = &app->basis,
      .vel_basis = &vms->basis_vel, 
      .phase_basis = &vms->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &vms->local_vel,
      .phase_range = &vms->local,
      .vel_map = vms->vel_map,
      .hamil_range = &vms->mom_hamil_range,
      .hamil = vms->mom_hamil,
      .model_id = vms->model_id,
      .hamil_id = vms->mom_hamil_id,
      .gamma_inv = vms->gamma_inv,
      .h_ij = vms->h_ij,
      .h_ij_inv = vms->h_ij_inv,
      .det_h = vms->det_h,
      .use_extended_hamil_def = vms->info.use_extended_hamil_def,
      .background_flows = vms->background_flows,
      .effective_potential = vms->effective_potential,
      .use_gpu = app->use_gpu,
      .quad_type = inp.quad_type
    };
    proj->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

    proj->correct_all_moms = false; 
    if (inp.correct_all_moms) {
      proj->correct_all_moms = true;

      int max_iter = inp.max_iter > 0 ? inp.max_iter : 100;
      double iter_eps = inp.iter_eps > 0 ? inp.iter_eps  : 1e-12;
      bool use_last_converged = inp.use_last_converged; 

      struct gkyl_vlasov_lte_correct_inp inp_corr = {
        .phase_grid = &vms->grid,
        .vel_grid = &vms->grid_vel, 
        .conf_basis = &app->basis,
        .vel_basis = &vms->basis_vel, 
        .phase_basis = &vms->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &vms->local_vel,
        .phase_range = &vms->local,
        .vel_map = vms->vel_map,
        .hamil_range = &vms->mom_hamil_range,
        .hamil = vms->mom_hamil,
        .model_id = vms->model_id,
        .hamil_id = vms->mom_hamil_id,
        .gamma_inv = vms->gamma_inv,
        .h_ij = vms->h_ij,
        .h_ij_inv = vms->h_ij_inv,
        .det_h = vms->det_h,
        .use_extended_hamil_def = vms->info.use_extended_hamil_def,
        .background_flows = vms->background_flows,
        .effective_potential = vms->effective_potential,
        .use_gpu = app->use_gpu,
        .quad_type = inp.quad_type,
        .max_iter = max_iter,
        .eps = iter_eps,
        .use_last_converged = use_last_converged, 
      };
      proj->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );
    }
  }
}

void
vm_species_projection_calc(gkyl_vlasov_app *app, const struct vm_species *vms, 
  struct vm_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    if (app->use_gpu) {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &vms->local_ext, proj->proj_host);
      gkyl_array_copy(f, proj->proj_host);
    }
    else {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &vms->local_ext, f);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_VLASOV_LTE) {
    int vdim = app->vdim;
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_V_drift, tm, &app->local_ext, proj->V_drift);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local_ext, proj->T_over_m);
    gkyl_array_scale(proj->T_over_m, 1.0/vms->info.mass);

    // Projection routines expect the LTE moments as a single array.
    gkyl_array_set_offset(proj->vlasov_lte_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
    gkyl_array_set_offset(proj->vlasov_lte_moms_host, 1.0, proj->V_drift, 1*app->basis.num_basis);
    gkyl_array_set_offset(proj->vlasov_lte_moms_host, 1.0, proj->T_over_m, (vdim+1)*app->basis.num_basis);

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->vlasov_lte_moms, proj->vlasov_lte_moms_host);

    // Project the LTE distribution function.
    // Projection routine also corrects the density of the projected distribution function.
    gkyl_vlasov_lte_proj_on_basis_advance(proj->proj_lte, &vms->local, &app->local, 
      proj->vlasov_lte_moms, f);

    // Correct all the moments of the projected LTE distribution function.
    if (proj->correct_all_moms) {
      struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(proj->corr_lte,
        f, proj->vlasov_lte_moms, &vms->local, &app->local);
    }
  }

  // Weight by the configuration-space (position-map) Jacobian so the stored
  // distribution carries J_x*J_v*f. This must come AFTER the projection is
  // complete: after gkyl_proj_on_basis for GKYL_PROJ_FUNC, and after the LTE
  // projection *and* its moment correction for GKYL_PROJ_VLASOV_LTE (the
  // correction matches moments of the J_v-weighted f to the target moments, so
  // the J_x weight must not be present during the correction). For the identity
  // map J_x=1 and this is a no-op.
  gkyl_vlasov_position_map_rescale_jacobpos(vms->pos_map, &vms->basis, &vms->local, f, f);
}

void
vm_species_projection_release(const struct gkyl_vlasov_app *app, const struct vm_proj *proj)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (app->use_gpu) {
      gkyl_array_release(proj->proj_host);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_VLASOV_LTE) { 
    gkyl_array_release(proj->dens);
    gkyl_array_release(proj->V_drift); 
    gkyl_array_release(proj->T_over_m);
    gkyl_array_release(proj->vlasov_lte_moms_host);
    gkyl_array_release(proj->vlasov_lte_moms);

    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_V_drift);
    gkyl_proj_on_basis_release(proj->proj_temp);

    gkyl_vlasov_lte_proj_on_basis_release(proj->proj_lte);
    if (proj->correct_all_moms) {
      gkyl_vlasov_lte_correct_release(proj->corr_lte);
    }
  } 
}
