#pragma once
#include <gkyl_basis.h>
#include <gkyl_dg_fv_interpolate.h>
#include <gkyl_rect_grid.h>

struct gkyl_flash_interp_vars {
    int cdim;
    struct gkyl_basis *basis;
    struct gkyl_rect_grid *grid_fluid, *grid_kinetic;
    gkyl_dg_fv_interpolate *interp; // Updater object for interpolation.
    struct gkyl_array *fluid_vel; // Fluid velocity on the extended fluid range.
};
