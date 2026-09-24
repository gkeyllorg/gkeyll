#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_dg_updater_gyrokinetic_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_gyrokinetic_acquire_eqn(const gkyl_dg_updater_gyrokinetic* gyrokinetic)
{
  return gkyl_dg_eqn_acquire(gyrokinetic->eqn_gyrokinetic);
}

struct gkyl_dg_updater_gyrokinetic*
gkyl_dg_updater_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc, double charge, double mass, 
  enum gkyl_gk_collisionless_type collless_type, const struct gk_geometry *gk_geom, 
  const struct gkyl_velocity_map *vel_map, void *aux_inp, bool use_gpu)
{
  struct gkyl_dg_updater_gyrokinetic *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_gyrokinetic));

  up->use_gpu = use_gpu;

  up->eqn_gyrokinetic = gkyl_dg_gyrokinetic_new(cbasis, pbasis, conf_range, phase_range, 
    charge, mass, collless_type, gk_geom, vel_map, up->use_gpu);

  struct gkyl_dg_gyrokinetic_auxfields *gk_inp = aux_inp;
  gkyl_gyrokinetic_set_auxfields(up->eqn_gyrokinetic, *gk_inp);

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int up_dirs[GKYL_MAX_DIM] = {0};
  int num_up_dirs = cdim+1;
  for (int d=0; d<num_up_dirs; ++d) up_dirs[d] = d;

  int zero_flux_flags[2*GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; ++d) {
    zero_flux_flags[d] = is_zero_flux_bc[d]? 1 : 0;
    zero_flux_flags[d+pdim] = is_zero_flux_bc[d+pdim]? 1 : 0;
  }
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1; // zero-flux BCs in vel-space

  up->up_gyrokinetic = gkyl_hyper_dg_new(grid, pbasis, up->eqn_gyrokinetic,
    num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->gyrokinetic_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_gyrokinetic_advance(struct gkyl_dg_updater_gyrokinetic *gyrokinetic,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(gyrokinetic->up_gyrokinetic, update_rng, fIn, cflrate, rhs);
  gyrokinetic->gyrokinetic_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_gyrokinetic_tm
gkyl_dg_updater_gyrokinetic_get_tm(const gkyl_dg_updater_gyrokinetic *gyrokinetic)
{
  return (struct gkyl_dg_updater_gyrokinetic_tm) {
    .gyrokinetic_tm = gyrokinetic->gyrokinetic_tm,
  };
}

void
gkyl_dg_updater_gyrokinetic_release(struct gkyl_dg_updater_gyrokinetic* gyrokinetic)
{
  gkyl_dg_eqn_release(gyrokinetic->eqn_gyrokinetic);
  gkyl_hyper_dg_release(gyrokinetic->up_gyrokinetic);
  gkyl_free(gyrokinetic);
}
