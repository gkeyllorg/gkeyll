#include <assert.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_dg_updater_moment.h>

void 
vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_boundary_fluxes *bflux)
{ 
  // allocate solver
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov, app->cdim, app->use_gpu);
  int ndim = app->cdim + app->vdim;
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    cells[d] = s->grid.cells[d];
    lower[d] = s->grid.lower[d];
    upper[d] = s->grid.upper[d];
    ghost[d] = 0;
  }

  // initialize moment solver
  for (int i=0; i<app->cdim; ++i) {
    cells[i] = 1;

    bflux->flux_arr[2*i] = mkarr(app->use_gpu, app->basis.num_basis, s->lower_ghost[i].volume);
    bflux->flux_arr[2*i+1] = mkarr(app->use_gpu, app->basis.num_basis, s->upper_ghost[i].volume);

    gkyl_range_init(&bflux->flux_r[2*i], ndim, s->lower_ghost[i].lower, s->lower_ghost[i].upper);
    gkyl_range_init(&bflux->flux_r[2*i+1], ndim, s->upper_ghost[i].lower, s->upper_ghost[i].upper);

    gkyl_range_init(&bflux->conf_r[2*i], app->cdim, s->lower_ghost[i].lower,
      s->lower_ghost[i].upper);
    gkyl_range_init(&bflux->conf_r[2*i+1], app->cdim, s->upper_ghost[i].lower,
      s->upper_ghost[i].upper);

    upper[i] = s->grid.lower[i] + s->grid.dx[i];

    gkyl_rect_grid_init(&bflux->boundary_grid[2*i], ndim, lower, upper, cells);

    upper[i] = s->grid.upper[i];
    lower[i] = s->grid.upper[i] - s->grid.dx[i];

    gkyl_rect_grid_init(&bflux->boundary_grid[2*i+1], ndim, lower, upper, cells);

    bflux->integ_moms[2*i] = gkyl_dg_updater_moment_new(&bflux->boundary_grid[2*i],
      &app->confBasis, &app->basis, &bflux->conf_r[2*i], &s->local_vel, &s->local, 
      s->model_id, 0, GKYL_F_MOMENT_M0M1M2, true, app->use_gpu);
    bflux->integ_moms[2*i+1] = gkyl_dg_updater_moment_new(&bflux->boundary_grid[2*i+1],
      &app->confBasis, &app->basis, &bflux->conf_r[2*i+1], &s->local_vel, &s->local, 
      s->model_id, 0, GKYL_F_MOMENT_M0M1M2, true, app->use_gpu);

    cells[i] = s->grid.cells[i];

    bflux->mom_arr[2*i] = mkarr(app->use_gpu, app->confBasis.num_basis, bflux->conf_r[2*i].volume);
    bflux->mom_arr[2*i+1] = mkarr(app->use_gpu, app->confBasis.num_basis, bflux->conf_r[2*i+1].volume);
  }
}

// computes rhs of the boundary flux
void
vm_species_bflux_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // zero ghost cells before calculation to ensure there's no residual data
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_clear_range(rhs, 0.0, &(species->lower_ghost[j]));
    gkyl_array_clear_range(rhs, 0.0, &(species->upper_ghost[j]));
  }
  // ghost cells of the rhs array are filled with the bflux
  // This is overwritten by the boundary conditions and is not being stored,
  // it is only currently used to calculate moments for other applications
  if (app->use_gpu) {
    gkyl_ghost_surf_calc_advance_cu(bflux->flux_slvr, &species->local_ext, fin, rhs);
  } else {
    gkyl_ghost_surf_calc_advance(bflux->flux_slvr, &species->local_ext, fin, rhs);
  }

  // only calculating integrated moments for use in the bflux source for now,
  // others can be added if applications require
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*j], rhs, &bflux->flux_r[2*j],
      &species->lower_ghost[j]);
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*j+1], rhs, &bflux->flux_r[2*j+1],
      &species->upper_ghost[j]);
    
    gkyl_dg_updater_moment_advance(bflux->integ_moms[2*j], &bflux->flux_r[2*j],
      &bflux->conf_r[2*j], bflux->flux_arr[2*j], bflux->mom_arr[2*j]);
    gkyl_dg_updater_moment_advance(bflux->integ_moms[2*j+1], &bflux->flux_r[2*j+1],
      &bflux->conf_r[2*j+1], bflux->flux_arr[2*j+1], bflux->mom_arr[2*j+1]);
  }
}

void
vm_species_bflux_release(const struct gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux)
{
  gkyl_ghost_surf_calc_release(bflux->flux_slvr);
  for (int i=0; i<2*app->cdim; ++i) {
    gkyl_array_release(bflux->mom_arr[i]);
    gkyl_array_release(bflux->flux_arr[i]);
    gkyl_dg_updater_moment_release(bflux->integ_moms[i]);
  }
}
