#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>
#include <gkyl_dg_lbo_gyrokinetic_diff.h>
#include <gkyl_dg_updater_lbo_gyrokinetic.h>
#include <gkyl_dg_updater_collisions_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_collisions* 
gkyl_dg_updater_lbo_gyrokinetic_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range, 
  struct gkyl_dg_lbo_gyrokinetic_drag_auxfields *drag_inp, struct gkyl_dg_lbo_gyrokinetic_diff_auxfields *diff_inp, 
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map,
  bool use_gpu)
{
  struct gkyl_dg_updater_collisions *up = gkyl_malloc(sizeof(gkyl_dg_updater_collisions));

  up->use_gpu = use_gpu;

  up->coll_drag = gkyl_dg_lbo_gyrokinetic_drag_new(conf_basis, phase_basis, conf_range,
    phase_grid, mass, gk_geom, vel_map, up->use_gpu);
  gkyl_lbo_gyrokinetic_drag_set_auxfields(up->coll_drag, *drag_inp);

  up->coll_diff = gkyl_dg_lbo_gyrokinetic_diff_new(conf_basis, phase_basis, conf_range,
    phase_grid, mass, gk_geom, vel_map, up->use_gpu); 
  gkyl_lbo_gyrokinetic_diff_set_auxfields(up->coll_diff, *diff_inp);

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + phase_basis->ndim - vdim;

  int zero_flux_flags[2*GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1;
  
  up->drag = gkyl_hyper_dg_new(phase_grid, phase_basis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);
  up->diff = gkyl_hyper_dg_new(phase_grid, phase_basis, up->coll_diff, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->diff_tm = 0.0; 
  up->drag_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_lbo_gyrokinetic_advance(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(lbo->drag, update_rng, fIn, cflrate, rhs);
  lbo->drag_tm += gkyl_time_diff_now_sec(wst);

  wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(lbo->diff, update_rng, fIn, cflrate, rhs);
  lbo->diff_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_lbo_gyrokinetic_tm
gkyl_dg_updater_lbo_gyrokinetic_get_tm(const gkyl_dg_updater_collisions *coll)
{
  return (struct gkyl_dg_updater_lbo_gyrokinetic_tm) {
    .drag_tm = coll->drag_tm,
    .diff_tm = coll->diff_tm
  };
}

void
gkyl_dg_updater_lbo_gyrokinetic_release(gkyl_dg_updater_collisions* coll)
{
  gkyl_dg_eqn_release(coll->coll_drag);
  gkyl_dg_eqn_release(coll->coll_diff);
  gkyl_hyper_dg_release(coll->drag);
  gkyl_hyper_dg_release(coll->diff);
  gkyl_free(coll);
}
