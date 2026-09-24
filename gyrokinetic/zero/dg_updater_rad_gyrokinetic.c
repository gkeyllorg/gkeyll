#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_collisions_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_collisions*
gkyl_dg_updater_rad_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_velocity_map *vel_map, void *aux_inp, bool use_gpu)
{

  struct gkyl_dg_updater_collisions *up = gkyl_malloc(sizeof(gkyl_dg_updater_collisions));
  up->use_gpu = use_gpu;
  up->coll_drag = gkyl_dg_rad_gyrokinetic_drag_new(conf_basis, phase_basis, phase_range, conf_range, vel_map, use_gpu);
  struct gkyl_dg_rad_gyrokinetic_auxfields *rad_inp = aux_inp;
  gkyl_rad_gyrokinetic_drag_set_auxfields(up->coll_drag, *rad_inp);

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + phase_basis->ndim - vdim;

  int zero_flux_flags[2*GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1;

  up->drag = gkyl_hyper_dg_new(grid, phase_basis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);
  
  up->drag_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_rad_gyrokinetic_advance(struct gkyl_dg_updater_collisions *rad,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(rad->drag, update_rng, fIn, cflrate, rhs);
  rad->drag_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_rad_gyrokinetic_tm
gkyl_dg_updater_rad_gyrokinetic_get_tm(const struct gkyl_dg_updater_collisions *coll)
{
  return (struct gkyl_dg_updater_rad_gyrokinetic_tm) {
    .drag_tm = coll->drag_tm
  };
}

void
gkyl_dg_updater_rad_gyrokinetic_release(struct gkyl_dg_updater_collisions* coll)
{
  gkyl_dg_eqn_release(coll->coll_drag);
  gkyl_hyper_dg_release(coll->drag);
  gkyl_free(coll);
}
