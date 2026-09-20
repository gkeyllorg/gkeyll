#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_neut_species_collisionless_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_neut_species_collisionless_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Divide out the velocity-space Jacobian (identity velocity map => copy).
  gkyl_vlasov_velocity_map_divide_jacobvel(species->vlasov_vel_map, &app->basis, &species->basis,
    &species->local, fin, species->f_no_J);

  // Compute the surface expansion of the phase space flux in configuration space.
  gkyl_dg_vlasov_conf_flux_surf_advance(gkcls->calc_conf_flux, &app->local, &species->local, &species->local_ext,
    species->conf_poisson_tensor, species->hamil, fin, species->cflrate, gkcls->conf_flux_surf);

  // Compute the surface expansion of the phase space flux in velocity space.
  gkyl_dg_vlasov_vel_flux_surf_advance(gkcls->calc_vel_flux, &app->local, &species->local,
    species->conf_poisson_tensor, species->hamil, 0, 0, 0,
    species->f_no_J, species->cflrate, gkcls->vel_flux_surf);

  gkyl_hyper_dg_advance(gkcls->slvr_vlasov, &species->local, fin, species->cflrate, rhs);

  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_collisionless_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame)
{
}

static void
gk_neut_species_collisionless_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_neut_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls)
{
  gkcls->collisionless_id = gkns->info.collisionless.type;
  gkcls->write_diagnostics = gkns->info.collisionless.write_diagnostics;

  gkcls->write_diags_func_neut = gk_neut_species_collisionless_write_diags_disabled;
  gkcls->rhs_func_neut = gk_neut_species_collisionless_rhs_disabled;

  if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL) {

    int cdim = app->cdim, vdim = gkns->info.vdim;
    int pdim = cdim+vdim;

    // Surface node counts for the nodal flux expansions. Tensor p=1 hybrid:
    // a velocity-direction surface has 2 nodes per configuration direction and
    // 4 nodes per remaining velocity direction; a configuration-direction
    // surface has 2 nodes per remaining configuration direction and 4 per
    // velocity direction. Must match the flux updaters' node counts.
    gkcls->num_surf_vel_nodes = (int) (pow(2, cdim) * pow(4, vdim - 1));
    gkcls->num_surf_conf_nodes = (int) (pow(2, cdim - 1) * pow(4, vdim));

    // Surface expansion of the phase-space flux in configuration space.
    gkcls->conf_flux_surf = mkarr(app->use_gpu, cdim*gkcls->num_surf_conf_nodes, gkns->local_ext.volume);
    struct gkyl_dg_vlasov_conf_flux_surf_inp inp_conf_flux = {
      .phase_grid = &gkns->grid,
      .conf_basis = &app->basis,
      .phase_basis = &gkns->basis,
      .vel_range = &gkns->local_vel,
      .vel_map = gkns->vlasov_vel_map,
      .pos_map = gkns->vlasov_pos_map,
      .hamil_range = &gkns->hamil_range,
      .model_id = gkns->model_id,
      .hamil_id = gkns->hamil_id,
      .use_gpu = app->use_gpu,
    };
    gkcls->calc_conf_flux = gkyl_dg_vlasov_conf_flux_surf_inew(&inp_conf_flux);

    // Surface expansion of the phase-space flux in velocity space.
    gkcls->vel_flux_surf = mkarr(app->use_gpu, vdim*gkcls->num_surf_vel_nodes, gkns->local_ext.volume);
    struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
      .phase_grid = &gkns->grid,
      .conf_basis = &app->basis,
      .phase_basis = &gkns->basis,
      .vel_map = gkns->vlasov_vel_map,
      .pos_map = gkns->vlasov_pos_map,
      .hamil_range = &gkns->hamil_range,
      .model_id = gkns->model_id,
      .hamil_id = gkns->hamil_id,
      .use_gpu = app->use_gpu,
    };
    gkcls->calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux);

    struct gkyl_dg_vlasov_inp inp_eqn = {
      .conf_basis = &app->basis,
      .phase_basis = &gkns->basis,
      .conf_range = &app->local,
      .hamil_range = &gkns->hamil_range,
      .phase_range = &gkns->local,
      .vel_map = gkns->vlasov_vel_map,
      .pos_map = gkns->vlasov_pos_map,
      .model_id = gkns->model_id,
      .hamil_id = gkns->hamil_id,
      .poisson_tensor_conf = gkns->conf_poisson_tensor,
      .hamil = gkns->hamil,
      .conf_flux_surf = gkcls->conf_flux_surf,
      .vel_flux_surf = gkcls->vel_flux_surf,
      .f_no_J = gkns->f_no_J,
      .use_gpu = app->use_gpu,
    };
    gkcls->eqn_vlasov = gkyl_dg_vlasov_inew(&inp_eqn);

    int up_dirs[GKYL_MAX_DIM];
    for (int d=0; d<pdim; ++d) {
      up_dirs[d] = d;
    }

    // Zero-flux BCs from the input in configuration space, always in velocity space.
    int zero_flux_flags[2*GKYL_MAX_DIM] = {false};
    for (int dir=0; dir<cdim; ++dir) {
      if (gkns->lower_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        zero_flux_flags[dir] = true;
      if (gkns->upper_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        zero_flux_flags[dir+pdim] = true;
    }
    for (int dir=cdim; dir<pdim; ++dir) {
      zero_flux_flags[dir] = zero_flux_flags[dir+pdim] = 1;
    }

    gkcls->slvr_vlasov = gkyl_hyper_dg_new(&gkns->grid, &gkns->basis, gkcls->eqn_vlasov,
      pdim, up_dirs, zero_flux_flags, 1, app->use_gpu);

    // Methods chosen at runtime.
    gkcls->rhs_func_neut = gk_neut_species_collisionless_rhs_enabled;
    if (gkcls->write_diagnostics) {
      gkcls->write_diags_func_neut = gk_neut_species_collisionless_write_diags_enabled;
    }
  }
}

void
gk_neut_species_collisionless_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkcls->rhs_func_neut(app, species, gkcls, fin, rhs);
}

void
gk_neut_species_collisionless_write_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  gkcls->write_diags_func_neut(app, gkns, gkcls, tm, frame);
}

void
gk_neut_species_collisionless_release(const struct gkyl_gyrokinetic_app *app, const struct gk_collisionless *gkcls)
{
  if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL) {
    gkyl_dg_vlasov_conf_flux_surf_release(gkcls->calc_conf_flux);
    gkyl_array_release(gkcls->conf_flux_surf);
    gkyl_dg_vlasov_vel_flux_surf_release(gkcls->calc_vel_flux);
    gkyl_array_release(gkcls->vel_flux_surf);

    gkyl_hyper_dg_release(gkcls->slvr_vlasov);
    gkyl_dg_eqn_release(gkcls->eqn_vlasov);

    if (gkcls->write_diagnostics) {
    }
  }
}
