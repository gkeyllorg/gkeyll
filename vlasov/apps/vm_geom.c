#include <gkyl_app.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_vlasov_triad_geom.h>


// Configuration-space c2p for sampling the GR geometry at physical coordinates
// on a non-uniform conf mesh via the position map (identity map => identity
// c2p, so uniform grids are unaffected).
static void
vm_geom_pos_c2p(const double *xcomp, double *xphys, void *ctx)
{
  gkyl_vlasov_position_map_eval_mc2p((struct gkyl_vlasov_position_map *) ctx, xcomp, xphys);
}

// Create the dg-gr-maxwell geometry
static void 
vm_dg_maxwell_geom_new(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_geom *vmg)
{

  // populate the geometry context struct with mass and spin
  struct gkyl_dg_gr_maxwell_geom_ctx ctx = {
    .mass_bh = vmg->mass_bh,
    .spin_bh = vmg->spin_bh,
  };

  int is_np[3] = { 1, 1, 1 };
  for (int d=0; d<app->num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  // Record which configured field boundaries use the theta-pole BC.
  for (int i=0; i<app->cdim; ++i) {
    const enum gkyl_field_bc_type *bc = i == 0 ? vm_app_inp->field.bcx :
      i == 1 ? vm_app_inp->field.bcy : vm_app_inp->field.bcz;

    vmg->theta_pole_lo[i] = is_np[i] && bc[0] == GKYL_FIELD_THETA_POLE;
    vmg->theta_pole_up[i] = is_np[i] && bc[1] == GKYL_FIELD_THETA_POLE;
  }

  // Evaluation of geometry at surface and volume nodal points.
  // Lapse - \alpha in the ADM split
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* lapse_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 1, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_dg_gr_maxwell_preset_lapse(vmg->triad_preset_geom_type), .ctx = &ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
  vmg->lapse_init = gkyl_surf_and_vol_node_arrays_new(lapse_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(lapse_proj, 0.0, &app->local_ext, vmg->lapse_init);
  if (app->use_gpu) {
    vmg->lapse = gkyl_surf_and_vol_node_copy_to_device(vmg->lapse_init, app->cdim);
  }
  else {
    vmg->lapse = gkyl_surf_and_vol_node_arrays_acquire(vmg->lapse_init);
  }
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(lapse_proj);

  // shift - \beta^i components in the ADM split (contravariant)
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* shift_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 3, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_dg_gr_maxwell_preset_shift(vmg->triad_preset_geom_type), .ctx = &ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
  vmg->shift_init = gkyl_surf_and_vol_node_arrays_new(shift_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(shift_proj, 0.0, &app->local_ext, vmg->shift_init);
  if (app->use_gpu) {
    vmg->shift = gkyl_surf_and_vol_node_copy_to_device(vmg->shift_init, app->cdim);
  }
  else {
    vmg->shift = gkyl_surf_and_vol_node_arrays_acquire(vmg->shift_init);
  }
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(shift_proj);

  // geom_factor_con - contravariant factors for geometric source terms
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* geom_factor_con_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 3, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_dg_gr_maxwell_preset_geom_factor_con(vmg->triad_preset_geom_type), .ctx = &ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
  vmg->geom_factor_con_init = gkyl_surf_and_vol_node_arrays_new(geom_factor_con_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(geom_factor_con_proj, 0.0, &app->local_ext, vmg->geom_factor_con_init);
  if (app->use_gpu) {
    vmg->geom_factor_con = gkyl_surf_and_vol_node_copy_to_device(vmg->geom_factor_con_init, app->cdim);
  }
  else {
    vmg->geom_factor_con = gkyl_surf_and_vol_node_arrays_acquire(vmg->geom_factor_con_init);
  }
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(geom_factor_con_proj);

  // h_ij - Covariant components of the spatial metric (assumed to allways be a upper 
  // triangular matrix of 6 unique elements)
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* h_ij_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 6, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_dg_gr_maxwell_preset_h_ij(vmg->triad_preset_geom_type), .ctx = &ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
  vmg->h_ij_init = gkyl_surf_and_vol_node_arrays_new(h_ij_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(h_ij_proj, 0.0, &app->local_ext, vmg->h_ij_init);
  if (app->use_gpu) {
    vmg->h_ij = gkyl_surf_and_vol_node_copy_to_device(vmg->h_ij_init, app->cdim);
  }
  else {
    vmg->h_ij = gkyl_surf_and_vol_node_arrays_acquire(vmg->h_ij_init);
  }
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(h_ij_proj);

  // h_ij_inv - Contravariant components of the spatial metric (upper triangular
  // matrix of 6 unique elements).
  // Allocate arrays for specified metric inverse
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* h_ij_inv_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 6, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_dg_gr_maxwell_preset_h_ij_inv(vmg->triad_preset_geom_type), .ctx = &ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
  vmg->h_ij_inv_init = gkyl_surf_and_vol_node_arrays_new(h_ij_inv_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(h_ij_inv_proj, 0.0, &app->local_ext, vmg->h_ij_inv_init);
  if (app->use_gpu) {
    vmg->h_ij_inv = gkyl_surf_and_vol_node_copy_to_device(vmg->h_ij_inv_init, app->cdim);
  }
  else {
    vmg->h_ij_inv = gkyl_surf_and_vol_node_arrays_acquire(vmg->h_ij_inv_init);
  }
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(h_ij_inv_proj);

  // Allocate arrays for the metric determinant (computed from J = sqrt(det(h_ij)))
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* det_h_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = 1, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_dg_gr_maxwell_preset_det_h(vmg->triad_preset_geom_type), .ctx = &ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
  vmg->det_h_init = gkyl_surf_and_vol_node_arrays_new(det_h_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(det_h_proj, 0.0, &app->local_ext, vmg->det_h_init);
  if (app->use_gpu) {
    vmg->det_h = gkyl_surf_and_vol_node_copy_to_device(vmg->det_h_init, app->cdim);
  }
  else {
    vmg->det_h = gkyl_surf_and_vol_node_arrays_acquire(vmg->det_h_init);
  }
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(det_h_proj);

  // Vierbeins used to transform GR-Maxwell fields into local Lorentz-force fields.
  if (vmg->has_gr_em_triad_coupling) {
    // vdim needed for vierbein size
    int vdim = app->vdim;
    struct gkyl_triad_geom_ctx triad_ctx = {
      .mass_bh = vmg->mass_bh,
      .spin_bh = vmg->spin_bh,
    };

    struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* vierb_cov_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = vdim*vdim, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_vlasov_triad_preset_vierbein(vdim, vmg->triad_preset_geom_type), .ctx = &triad_ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
    vmg->vierb_cov_init = gkyl_surf_and_vol_node_arrays_new(vierb_cov_proj, app->local_ext.volume, app->use_gpu);
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(vierb_cov_proj, 0.0, &app->local_ext, vmg->vierb_cov_init);
    if (app->use_gpu) {
      vmg->vierb_cov = gkyl_surf_and_vol_node_copy_to_device(vmg->vierb_cov_init, app->cdim);
    }
    else {
      vmg->vierb_cov = gkyl_surf_and_vol_node_arrays_acquire(vmg->vierb_cov_init);
    }
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(vierb_cov_proj);

    struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* vierb_con_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_inew( &(struct gkyl_dg_gr_maxwell_surf_and_vol_nodes_inp) {
      .grid = &app->grid, .basis = &app->basis, .num_ret_vals = vdim*vdim, .polyorder = vm_app_inp->poly_order,
      .eval = gkyl_vlasov_triad_preset_vierbein_inv(vdim, vmg->triad_preset_geom_type), .ctx = &triad_ctx,
      // Sample the GR geometry at physical conf coordinates on mapped meshes.
      .c2p_func = vm_geom_pos_c2p, .c2p_func_ctx = app->pos_map,
    });
    vmg->vierb_con_init = gkyl_surf_and_vol_node_arrays_new(vierb_con_proj, app->local_ext.volume, app->use_gpu);
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(vierb_con_proj, 0.0, &app->local_ext, vmg->vierb_con_init);
    if (app->use_gpu) {
      vmg->vierb_con = gkyl_surf_and_vol_node_copy_to_device(vmg->vierb_con_init, app->cdim);
    }
    else {
      vmg->vierb_con = gkyl_surf_and_vol_node_arrays_acquire(vmg->vierb_con_init);
    }
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(vierb_con_proj);
  }

}


// Initialize geom object.
void
vm_geom_init(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_geom *vmg)
{

  *vmg = (struct vm_geom) { 0 };
  vmg->info = vm_app_inp->geom;

  // Set the black hole parameters
  vmg->mass_bh = vmg->info.mass_bh;
  vmg->spin_bh = vmg->info.spin_bh;

  // Set triad preset-geometry parameters (used by TRIAD/TRIAD_GR models).
  vmg->use_preset_geom = vmg->info.use_preset_geom;
  vmg->triad_preset_geom_type = vmg->info.triad_preset_geom_type;

  // Determine if gr-maxwell are used
  vmg->has_gr_fields = false;
  if ( (vm_app_inp->skip_field == false) && ( vm_app_inp->field.field_id == GKYL_FIELD_GR_D_B ) ) {
    vmg->has_gr_fields = true;
  }
  // (GR only) All coupled species are required to be either triad / triad_gr
  vmg->has_gr_em_triad_coupling = vmg->has_gr_fields && vm_app_inp->num_species > 0;
  for (int i=0; i<vm_app_inp->num_species; ++i) {
    enum gkyl_model_id model_id = vm_app_inp->species[i].model_id;
    vmg->has_gr_em_triad_coupling = vmg->has_gr_em_triad_coupling &&
      (model_id == GKYL_MODEL_TRIAD || model_id == GKYL_MODEL_TRIAD_GR);
  }

  // If the fields are included, and are GR, then build the gr maxwell geometry
  if ( vmg->has_gr_fields ) {
    vm_dg_maxwell_geom_new(vm_app_inp, app, vmg);
  }

}

void
vm_geom_release(const gkyl_vlasov_app* app, const struct vm_geom *vmg)
{

  // Release all variables associated with vm geom for dg-gr-maxwells
  if ( vmg->has_gr_fields ) {
    gkyl_surf_and_vol_node_arrays_release(vmg->lapse);
    gkyl_surf_and_vol_node_arrays_release(vmg->shift);
    gkyl_surf_and_vol_node_arrays_release(vmg->geom_factor_con);
    gkyl_surf_and_vol_node_arrays_release(vmg->h_ij);
    gkyl_surf_and_vol_node_arrays_release(vmg->h_ij_inv);
    gkyl_surf_and_vol_node_arrays_release(vmg->det_h);
    if (vmg->has_gr_em_triad_coupling) {
      gkyl_surf_and_vol_node_arrays_release(vmg->vierb_cov);
      gkyl_surf_and_vol_node_arrays_release(vmg->vierb_con);
    }
    gkyl_surf_and_vol_node_arrays_release(vmg->lapse_init);
    gkyl_surf_and_vol_node_arrays_release(vmg->shift_init);
    gkyl_surf_and_vol_node_arrays_release(vmg->geom_factor_con_init);
    gkyl_surf_and_vol_node_arrays_release(vmg->h_ij_init);
    gkyl_surf_and_vol_node_arrays_release(vmg->h_ij_inv_init);
    gkyl_surf_and_vol_node_arrays_release(vmg->det_h_init);
    if (vmg->has_gr_em_triad_coupling) {
      gkyl_surf_and_vol_node_arrays_release(vmg->vierb_cov_init);
      gkyl_surf_and_vol_node_arrays_release(vmg->vierb_con_init);
    }
  }
}
