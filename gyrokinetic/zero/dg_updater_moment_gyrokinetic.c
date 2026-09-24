#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_dg_updater_moment_priv.h>
#include <gkyl_mom_calc.h>
#include <gkyl_util.h>

struct gkyl_mom_type*
gkyl_dg_updater_moment_gyrokinetic_acquire_type(const gkyl_dg_updater_moment* moment)
{
  return gkyl_mom_type_acquire(moment->type);
}

int
gkyl_dg_updater_moment_gyrokinetic_num_mom(const gkyl_dg_updater_moment* moment)
{
  return gkyl_mom_type_num_mom(moment->type);
}

struct gkyl_dg_updater_moment*
gkyl_dg_updater_moment_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, double mass, double charge,
  const struct gkyl_velocity_map *vel_map, const struct gk_geometry *gk_geom,
  struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool is_integrated, bool use_gpu)
{
  gkyl_dg_updater_moment *up = gkyl_malloc(sizeof(gkyl_dg_updater_moment));
  up->use_gpu = use_gpu;

  if (is_integrated)
    up->type = gkyl_int_mom_gyrokinetic_new(conf_basis, phase_basis, conf_range, mass, charge, vel_map, gk_geom, phi, mom_type, use_gpu);
  else
    up->type = gkyl_mom_gyrokinetic_new(conf_basis, phase_basis, conf_range, mass, charge, vel_map, gk_geom, phi, mom_type, use_gpu);

  up->up_moment = gkyl_mom_calc_new(grid, up->type, use_gpu);

  up->moment_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_moment_gyrokinetic_advance(struct gkyl_dg_updater_moment *moment,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT mout)
{  
  struct timespec wst = gkyl_wall_clock();
  if (moment->use_gpu)
    gkyl_mom_calc_advance_cu(moment->up_moment, update_phase_rng, update_conf_rng, fIn, mout);
  else 
    gkyl_mom_calc_advance(moment->up_moment, update_phase_rng, update_conf_rng, fIn, mout);
  moment->moment_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_moment_tm
gkyl_dg_updater_moment_gyrokinetic_get_tm(const gkyl_dg_updater_moment *moment)
{
  return (struct gkyl_dg_updater_moment_tm) {
    .moment_tm = moment->moment_tm,
  };
}

void
gkyl_dg_updater_moment_gyrokinetic_release(gkyl_dg_updater_moment* moment)
{
  gkyl_mom_type_release(moment->type);
  gkyl_mom_calc_release(moment->up_moment);
  gkyl_free(moment);
}
