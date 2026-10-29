#include <gkyl_alloc.h>
#include <gkyl_flash_interp_vars_priv.h>

struct gkyl_flash_interp_vars*
gkyl_flash_interp_vars_new(struct gkyl_basis *basis, 
  struct gkyl_rect_grid *grid_kinetic, struct gkyl_rect_grid *grid_fluid,
  struct gkyl_range *range_fluid_ext) 
{
  struct gkyl_flash_interp_vars *up = gkyl_malloc(sizeof(struct gkyl_flash_interp_vars));

  up->cdim = grid_fluid->ndim;
  up->basis = basis;

  up->grid_fluid = grid_fluid;
  up->grid_kinetic = grid_kinetic;

  // Interpolation updater
  up->interp = gkyl_dg_fv_interpolate_new(basis, grid_kinetic, grid_fluid);

  // Fluid velocity buffer
  up->fluid_vel = gkyl_array_new(GKYL_DOUBLE, 3, range_fluid_ext->volume);

  return up;
}

// Interpolate momentum density from finite volume grid to kinetic grid.
void
gkyl_flash_interp_vars_rhou_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_fluid_in, struct gkyl_array *dg_rhou_out)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range_fluid);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(range_fluid, iter.idx);
    const double *fluid_d = (const double *)gkyl_array_cfetch(fv_fluid_in, linc);
    double *vel_d = (double *)gkyl_array_fetch(up->fluid_vel, linc);

    for (int d=0; d<3; ++d) {
      vel_d[d] = fluid_d[d+1]/fluid_d[0];
    }
  }

  // Interpolate fluid velocity vector onto DG grid.
  gkyl_dg_fv_interpolate_fv_to_modal_advance(up->interp,
    range_kinetic, range_fluid, up->fluid_vel, 0, dg_rhou_out, 0, 3);
}

// Interpolate fluid velocity vector from finite volume grid to kinetic grid.
void
gkyl_flash_interp_vars_vel_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_fluid_vel_in, struct gkyl_array *dg_vel_out)
{
  // Interpolate fluid velocity vector onto DG grid.
  gkyl_dg_fv_interpolate_fv_to_modal_advance(up->interp,
    range_kinetic, range_fluid, fv_fluid_vel_in, 0, dg_vel_out, 0, 3);
}

// Interpolate initial condition arrays [n, ux, uy, uz, T/m] onto kinetic grid.
void
gkyl_flash_interp_vars_init_dist_moms_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_moms_in, struct gkyl_array *dg_moms_out)
{
  // Interpolate moments onto DG grid.
  gkyl_dg_fv_interpolate_fv_to_modal_advance(up->interp,
    range_kinetic, range_fluid, fv_moms_in, 0, dg_moms_out, 0, 5);
}

void
gkyl_flash_interp_vars_coupling_fluxes_to_fv(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *dg_coupling_fluxes, struct gkyl_array *fv_coupling_fluxes) 
{
  gkyl_dg_fv_interpolate_modal_to_fv_advance(up->interp,
    range_kinetic, range_fluid, dg_coupling_fluxes, 0, fv_coupling_fluxes, 0, 3);
}

void
gkyl_flash_interp_vars_release(struct gkyl_flash_interp_vars *up)
{
  gkyl_dg_fv_interpolate_release(up->interp);
  gkyl_array_release(up->fluid_vel);
  gkyl_free(up);
}