#include <assert.h>
#include <gkyl_vlasov_priv.h>

static void
vm_species_collisionless_rhs_enabled(gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Set values of q/m*EM and the total potentials based on field and external forces configuration. 
  gkyl_array_clear(cls->qmem, 0.0); 
  if (cls->has_app_accel) {
    gkyl_array_accumulate_range(cls->qmem, 1.0, cls->app_accel, &app->local);
  }

  // A field object always exists; for the null field (GKYL_FIELD_NULL) none of
  // these force terms fire (has_ext_em is false and the field type matches none).
  if (app->field->has_ext_em) {
    gkyl_array_accumulate_range(cls->qmem, cls->qbym, app->field->ext_em, &app->local);
  }

  if (vms->field_id == GKYL_FIELD_E_B) {
    // The Lorentz force uses the physical E, B. On a mapped grid the field stores
    // J*E, J*B, so divide out the conf Jacobian first (identity map => no buffer).
    const struct gkyl_array *em_force = em;
    if (cls->em_no_J) {
      gkyl_vlasov_position_map_divide_jacobpos_conf(vms->pos_map, &app->local,
        cls->em_no_J->ncomp, em, cls->em_no_J);
      em_force = cls->em_no_J;
    }
    gkyl_array_accumulate_range(cls->qmem, cls->qbym, em_force, &app->local);
  }
  else if (cls->has_gr_em_triad_coupling) {
    gkyl_dg_gr_maxwell_lorentz_conf_advance(cls->calc_lorentz, &app->local,
      app->vm_geom->lapse, app->vm_geom->shift, app->vm_geom->h_ij, app->vm_geom->h_ij_inv,
      app->vm_geom->det_h, app->vm_geom->vierb_cov, app->vm_geom->vierb_con, em, cls->qmem);
  }
  else if (vms->field_id == GKYL_FIELD_PHI) {
    gkyl_array_set_offset(cls->pot_tot, cls->qbym, app->field->phi, 0);
    if (app->field->has_ext_pot) {
      gkyl_array_accumulate_offset(cls->pot_tot, cls->qbym, app->field->ext_pot, 0);
    }
  }

  // Divide out velocity-space Jacobian.
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, fin, vms->f_no_J);

  // Compute the surface expansion of the phase space flux in configuration space. 
  if (vms->model_id == GKYL_MODEL_TRIAD || vms->hamil_id == GKYL_HAMIL_PHASE) {
    gkyl_dg_vlasov_conf_flux_surf_advance(cls->calc_conf_flux, &app->local, &vms->local, &vms->local_ext, 
      vms->conf_poisson_tensor, vms->hamil, fin, vms->cflrate, cls->conf_flux_surf);
  }

  // Compute the surface expansion of the phase space flux in velocity space.
  gkyl_dg_vlasov_vel_flux_surf_advance(cls->calc_vel_flux, &app->local, &vms->local,
    vms->conf_poisson_tensor, vms->hamil, cls->qmem, cls->pot_tot, vms->rad,
    vms->f_no_J, vms->cflrate, cls->vel_flux_surf);

  gkyl_hyper_dg_advance(cls->slvr, &vms->local, fin, vms->cflrate, rhs);

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_collisionless_rhs_disabled(gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs)
{
  // Do nothing.
}

void 
vm_species_collisionless_init(struct gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;
  enum gkyl_basis_type b_type = app->basis.b_type;

  // Allocate array to store q/m*(E,B) or potentials (q/m*phi + m*phi_g, q/m*A) depending on equation system. 
  // Note: the potentials are the total potentials and thus can include both (or either) gravitational
  // or electrostatic interactions. 
  cls->qbym = vms->info.charge/vms->info.mass;
  cls->has_gr_em_triad_coupling = app->vm_geom->has_gr_em_triad_coupling &&
    (vms->field_id == GKYL_FIELD_GR_D_B) &&
    (vms->model_id == GKYL_MODEL_TRIAD || vms->model_id == GKYL_MODEL_TRIAD_GR);
  cls->qmem = mkarr(app->use_gpu, 8*app->basis.num_basis, app->local_ext.volume);
  cls->calc_lorentz = 0;
  cls->calc_current_dep = 0;
  cls->pot_tot = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);
  // Buffer for the physical E, B (field stores J*E, J*B on a mapped grid) used in
  // the Lorentz force; only allocated when the standard Maxwell field is on a
  // non-identity position map.
  cls->em_no_J = ((vms->field_id == GKYL_FIELD_E_B) && (!vms->pos_map->is_identity)) ?
    mkarr(app->use_gpu, 8*app->basis.num_basis, app->local_ext.volume) : 0;

  // Initialize applied acceleration for use in force update. 
  cls->app_accel = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(cls->app_accel, 0.0);
  cls->has_app_accel = false;
  cls->app_accel_evolve = false;
  // setup applied acceleration
  if (vms->info.app_accel) {
    cls->has_app_accel = true;
    if (vms->info.app_accel_evolve) {
      cls->app_accel_evolve = vms->info.app_accel_evolve;
    }

    // Host-side distribution function for projection on GPUs.
    cls->app_accel_host = app->use_gpu ? mkarr(false, cls->app_accel->ncomp, cls->app_accel->size)
                                       : gkyl_array_acquire(cls->app_accel);      
    cls->app_accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->basis.poly_order+1,
      3, vms->info.app_accel, vms->info.app_accel_ctx);
  }

  // Determine which forces we need based on combination of field ID and presence 
  // of applied accelerations and external fields/potentials. 
  // A field object always exists (a null field with field_id == GKYL_FIELD_NULL
  // when none is present), so dispatch on the field type. With no field force
  // (GKYL_FIELD_NULL), an applied acceleration still acts as an electric force.
  cls->has_E = false;
  cls->has_B = false;
  cls->has_phi = false;
  if (cls->has_gr_em_triad_coupling) {
    cls->has_E = true;
    cls->has_B = true;
  }
  else if (vms->field_id == GKYL_FIELD_E_B || app->field->has_ext_em) {
    cls->has_E = true;
    cls->has_B = true;
  }
  if (vms->field_id == GKYL_FIELD_PHI) {
    cls->has_phi = true;
  }
  if (vms->field_id == GKYL_FIELD_NULL && cls->has_app_accel) {
    cls->has_E = true;
  }

  if (cls->has_gr_em_triad_coupling) {
    struct gkyl_dg_gr_maxwell_lorentz_conf_inp inp_lorentz = {
      .pos_map = vms->pos_map,
      .conf_grid = &app->grid,
      .conf_basis = &app->basis,
      .vdim = app->vdim,
      .qbym = cls->qbym,
      .chi = app->field->info.elcErrorSpeedFactor,
      .gamma = app->field->info.mgnErrorSpeedFactor,
      .K_phi = app->field->info.K_phi,
      .K_psi = app->field->info.K_psi,
      .use_gpu = app->use_gpu,
    };
    cls->calc_lorentz = gkyl_dg_gr_maxwell_lorentz_conf_inew(&inp_lorentz);
    struct gkyl_dg_gr_maxwell_current_deposition_inp inp_current_dep = {
      .conf_basis = &app->basis,
      .vdim = vdim,
      .use_gpu = app->use_gpu,
    };
    cls->calc_current_dep = gkyl_dg_gr_maxwell_current_deposition_inew(&inp_current_dep);
  }

  cls->use_lo = false; 
  if (vms->info.use_lo == true) {
    cls->use_lo = true; 
  }
  cls->use_vierbein = false; 
  if (vms->info.use_vierbein == true) {
    cls->use_vierbein = true; 
  }
  cls->use_extended_hamil_def = false; 
  if (vms->info.use_extended_hamil_def == true) {
    cls->use_extended_hamil_def = true; 
  }
  cls->use_preset_geom = false;
  if (vms->geom && vms->geom->use_preset_geom) {
    cls->use_preset_geom = true;
  }

  // Select the number of nodes, with case for hybrid-tensor.
  int highorder = cls->use_lo ? 0 : 1;

  // p + 1 is equivalent to p + 2 for ser p1
  if ((app->poly_order == 1) && (b_type == GKYL_BASIS_MODAL_SERENDIPITY)) {
    highorder = 0;
  }
  cls->num_surf_vel_nodes = pow(app->poly_order+1+highorder, pdim - 1);
  if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (app->poly_order == 1)) {
    // Tensor p=1 hybrid: a velocity-direction surface has 2 nodes per
    // configuration direction and 3 (lo) or 4 (ho) nodes per remaining
    // velocity direction. Must match the vel_flux updater's
    // num_nodes_conf*num_nodes_vel.
    int nq_vel = cls->use_lo ? 3 : 4;
    cls->num_surf_vel_nodes = (int) ( pow(2, cdim) * pow(nq_vel, vdim - 1) );
  }

  // Allocate nodal surface expansion of velocity space flux array (conf). 
  if (vms->model_id == GKYL_MODEL_TRIAD || vms->hamil_id == GKYL_HAMIL_PHASE) {

    // Compute the number of configuration space nodes, with case for hybrid-tensor.
    cls->num_surf_conf_nodes = pow(app->poly_order+1+highorder,pdim - 1);
    if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (app->poly_order == 1)) {
      // Tensor p=1 hybrid: a configuration-direction surface has 2 nodes per
      // remaining configuration direction and 3 (lo) or 4 (ho) nodes per
      // velocity direction.
      int nq_vel = cls->use_lo ? 3 : 4;
      cls->num_surf_conf_nodes = (int) ( pow(2, cdim - 1) * pow(nq_vel, vdim) );
    }

    cls->conf_flux_surf = mkarr(app->use_gpu, cdim*cls->num_surf_conf_nodes, vms->local_ext.volume);
    struct gkyl_dg_vlasov_conf_flux_surf_inp inp_conf_flux = {
      .phase_grid = &vms->grid, 
      .conf_basis = &app->basis,
      .phase_basis = &vms->basis,
      .vel_range = &vms->local_vel,
      .vel_map = vms->vel_map,
      .pos_map = vms->pos_map,
      .hamil_range = &vms->hamil_range,
      .skip_cell_thresh = vms->info.skip_cell_thresh > 0.0 ? vms->info.skip_cell_thresh : 0.0, 
      .model_id = vms->model_id,
      .hamil_id = vms->hamil_id,
      .use_lo = cls->use_lo,
      .use_gpu = app->use_gpu,
    }; 
    cls->calc_conf_flux = gkyl_dg_vlasov_conf_flux_surf_inew(&inp_conf_flux); 
  }

  // Allocate nodal surface expansion of velocity space flux array (vel).
  cls->vel_flux_surf = mkarr(app->use_gpu, vdim*cls->num_surf_vel_nodes, vms->local_ext.volume);
  struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
    .phase_grid = &vms->grid,
    .conf_basis = &app->basis,
    .phase_basis = &vms->basis,
    .vel_map = vms->vel_map,
    .pos_map = vms->pos_map,
    .hamil_range = &vms->hamil_range,
    .skip_cell_thresh = vms->info.skip_cell_thresh > 0.0 ? vms->info.skip_cell_thresh : 0.0, 
    .model_id = vms->model_id,
    .hamil_id = vms->hamil_id,
    .has_E = cls->has_E, 
    .has_phi = cls->has_phi, 
    .has_B = cls->has_B, 
    .has_rad = vms->has_rad, 
    .use_lo = cls->use_lo,
    .use_gpu = app->use_gpu,
  }; 
  cls->calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux);

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &app->basis,
    .phase_basis = &vms->basis,
    .conf_range =  &app->local,
    .hamil_range = &vms->hamil_range,
    .phase_range = &vms->local,
    .vel_map = vms->vel_map,
    .pos_map = vms->pos_map,
    .skip_cell_thresh = vms->info.skip_cell_thresh > 0.0 ? vms->info.skip_cell_thresh : 0.0,
    .model_id = vms->model_id,
    .hamil_id = vms->hamil_id,
    .has_E = cls->has_E, 
    .has_phi = cls->has_phi, 
    .has_B = cls->has_B, 
    .has_rad = vms->has_rad, 
    .poisson_tensor_conf = vms->conf_poisson_tensor,
    .hamil = vms->hamil,
    .qmem = cls->qmem, 
    .pot_tot = cls->pot_tot, 
    .conf_flux_surf = cls->conf_flux_surf,
    .vel_flux_surf = cls->vel_flux_surf, 
    .f_no_J = vms->f_no_J, 
    .rad = vms->rad, 
    .use_lo = cls->use_lo,
    .use_gpu = app->use_gpu,
  };  
  // Construct Vlasov equation and Hyper DG object for updating equation. 
  cls->eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  int up_dirs[GKYL_MAX_DIM];
  for (int d=0; d<pdim; ++d) {
    up_dirs[d] = d;
  }
  int num_up_dirs = pdim;

  // Default: no zero-flux BCs in any configuration-space direction.
  int zero_flux_flags[2*GKYL_MAX_DIM] = {false}; 
  for (int dir=0; dir<cdim; ++dir) {
    if (vms->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
      zero_flux_flags[dir] = true;
    }
    if (vms->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
      zero_flux_flags[dir+pdim] = true;
    }
  }
  // Default: zero-flux BCs in velocity space
  for (int dir=cdim; dir<pdim; ++dir) {
    zero_flux_flags[dir] = zero_flux_flags[dir+pdim] = 1;
  }

  cls->slvr = gkyl_hyper_dg_new(&vms->grid, &vms->basis, cls->eqn, 
    num_up_dirs, up_dirs, zero_flux_flags, 1, app->use_gpu);

  cls->rhs_func = vm_species_collisionless_rhs_enabled;
  if (vms->info.no_collisionless_terms) {
    cls->rhs_func = vm_species_collisionless_rhs_disabled;
  }
}

void
vm_species_collisionless_app_accel(gkyl_vlasov_app *app, struct vm_collisionless *cls, double tm)
{
  if (cls->has_app_accel) {
    gkyl_proj_on_basis_advance(cls->app_accel_proj, tm, &app->local_ext, cls->app_accel_host);
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(cls->app_accel, cls->app_accel_host);
    }
  }
}

void
vm_species_collisionless_rhs(gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs)
{
  cls->rhs_func(app, vms, cls, fin, em, rhs);
}

void
vm_species_collisionless_release(const struct gkyl_vlasov_app *app, 
  const struct vm_species *vms, const struct vm_collisionless *cls)
{
  if (vms->model_id == GKYL_MODEL_TRIAD || vms->hamil_id == GKYL_HAMIL_PHASE) {
    gkyl_dg_vlasov_conf_flux_surf_release(cls->calc_conf_flux);
    gkyl_array_release(cls->conf_flux_surf);
  }
  gkyl_dg_vlasov_vel_flux_surf_release(cls->calc_vel_flux);
  gkyl_array_release(cls->qmem);
  if (cls->has_gr_em_triad_coupling) {
    gkyl_dg_gr_maxwell_lorentz_conf_release(cls->calc_lorentz);
    gkyl_dg_gr_maxwell_current_deposition_release(cls->calc_current_dep);
  }
  if (cls->em_no_J) gkyl_array_release(cls->em_no_J);
  gkyl_array_release(cls->pot_tot); 
  gkyl_array_release(cls->app_accel);
  if (cls->has_app_accel) {
    gkyl_array_release(cls->app_accel_host);
    gkyl_proj_on_basis_release(cls->app_accel_proj);
  } 
  gkyl_array_release(cls->vel_flux_surf);  
  gkyl_hyper_dg_release(cls->slvr);
  gkyl_dg_eqn_release(cls->eqn);
}
