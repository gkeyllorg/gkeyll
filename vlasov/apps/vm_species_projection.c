#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_projection_init(struct gkyl_vlasov_app *app, struct vm_species *s, 
  struct gkyl_vlasov_projection inp, struct vm_proj *proj)
{
  proj->proj_id = inp.proj_id;
  proj->model_id = s->model_id;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.func,
        .ctx = inp.ctx_func,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_VLASOV_LTE) {
    int vdim = app->vdim; 
    proj->dens = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->V_drift = mkarr(false, vdim*app->confBasis.num_basis, app->local_ext.volume);
    proj->T_over_m = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vlasov_lte_moms_host = mkarr(false, (vdim+2)*app->confBasis.num_basis, app->local_ext.volume);

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_V_drift = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, vdim, inp.V_drift, inp.ctx_V_drift);
    proj->proj_temp = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.temp, inp.ctx_temp);

    proj->vlasov_lte_moms = mkarr(app->use_gpu, (vdim+2)*app->confBasis.num_basis, app->local_ext.volume);

    struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
      .phase_grid = &s->grid,
      .vel_grid = &s->grid_vel, 
      .conf_basis = &app->confBasis,
      .vel_basis = &app->velBasis, 
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel,
      .phase_range = &s->local,
      .gamma = s->gamma,
      .gamma_inv = s->gamma_inv,
      .h_ij = s->h_ij, 
      .h_ij_inv = s->h_ij_inv, 
      .det_h = s->det_h,
      .hamil = s->hamil,
      .model_id = s->model_id,
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
        .phase_grid = &s->grid,
        .vel_grid = &s->grid_vel, 
        .conf_basis = &app->confBasis,
        .vel_basis = &app->velBasis, 
        .phase_basis = &app->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &s->local_vel,
        .phase_range = &s->local,
        .gamma = s->gamma,
        .gamma_inv = s->gamma_inv,
        .h_ij = s->h_ij,
        .h_ij_inv = s->h_ij_inv,
        .det_h = s->det_h,
        .hamil = s->hamil,
        .model_id = s->model_id,
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
vm_species_projection_calc(gkyl_vlasov_app *app, const struct vm_species *s, 
  struct vm_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    if (app->use_gpu) {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local_ext, proj->proj_host);
      gkyl_array_copy(f, proj->proj_host);
    }
    else {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local_ext, f);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_VLASOV_LTE) {
    int vdim = app->vdim;
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_V_drift, tm, &app->local_ext, proj->V_drift);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local_ext, proj->T_over_m);
    gkyl_array_scale(proj->T_over_m, 1.0/s->info.mass);

    // Projection routines expect the LTE moments as a single array.
    gkyl_array_set_offset(proj->vlasov_lte_moms_host, 1.0, proj->dens, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->vlasov_lte_moms_host, 1.0, proj->V_drift, 1*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->vlasov_lte_moms_host, 1.0, proj->T_over_m, (vdim+1)*app->confBasis.num_basis);

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->vlasov_lte_moms, proj->vlasov_lte_moms_host);

    // Project the LTE distribution function.
    // Projection routine also corrects the density of the projected distribution function.
    gkyl_vlasov_lte_proj_on_basis_advance(proj->proj_lte, &s->local, &app->local, 
      proj->vlasov_lte_moms, f);

    // Correct all the moments of the projected LTE distribution function.
    if (proj->correct_all_moms) {
      struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(proj->corr_lte, 
        f, proj->vlasov_lte_moms, &s->local, &app->local);
    } 
  } 
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
