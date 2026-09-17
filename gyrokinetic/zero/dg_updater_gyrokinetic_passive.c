#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_gyrokinetic_passive.h>
#include <gkyl_dg_updater_gyrokinetic_passive.h>
#include <gkyl_dg_updater_gyrokinetic_passive_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_gyrokinetic_passive_acquire_eqn(const gkyl_dg_updater_gyrokinetic_passive* up)
{
  return gkyl_dg_eqn_acquire(up->eqn_gyrokinetic_passive);
}

gkyl_dg_updater_gyrokinetic_passive*
gkyl_dg_updater_gyrokinetic_passive_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom,
  const struct gkyl_velocity_map *vel_map,
  void *aux_inp, bool use_gpu)
{
  struct gkyl_dg_updater_gyrokinetic_passive *up =
    gkyl_malloc(sizeof(struct gkyl_dg_updater_gyrokinetic_passive));

  up->use_gpu = use_gpu;

  // Create the passive GK equation object.
  up->eqn_gyrokinetic_passive = gkyl_dg_gyrokinetic_passive_new(cbasis, pbasis, conf_range, phase_range,
    charge, mass, gk_geom, vel_map, use_gpu);

  // Set auxiliary fields.
  struct gkyl_dg_gyrokinetic_passive_auxfields *gkp_inp = aux_inp;
  gkyl_gyrokinetic_passive_set_auxfields(up->eqn_gyrokinetic_passive, *gkp_inp);

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;

  // Update only conf-space directions (no vpar advection for passive type).
  int up_dirs[GKYL_MAX_DIM] = {0};
  int num_up_dirs = cdim;
  for (int d = 0; d < num_up_dirs; d++) up_dirs[d] = d;

  // Zero-flux flags: conf-space BCs from input.
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0};
  for (int d = 0; d < cdim; d++) {
    zero_flux_flags[d]      = is_zero_flux_bc[d]      ? 1 : 0;
    zero_flux_flags[d+pdim] = is_zero_flux_bc[d+pdim] ? 1 : 0;
  }
  for (int d = cdim; d < pdim; d++)
    zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1; // zero-flux BCs in vel-space

  up->up_gyrokinetic_passive = gkyl_hyper_dg_new(grid, pbasis, up->eqn_gyrokinetic_passive,
    num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->gyrokinetic_passive_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_gyrokinetic_passive_advance(gkyl_dg_updater_gyrokinetic_passive *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *fIn,
  struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->up_gyrokinetic_passive, update_rng, fIn, cflrate, rhs);
  up->gyrokinetic_passive_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_gyrokinetic_passive_tm
gkyl_dg_updater_gyrokinetic_passive_get_tm(const gkyl_dg_updater_gyrokinetic_passive *up)
{
  return (struct gkyl_dg_updater_gyrokinetic_passive_tm) {
    .gyrokinetic_passive_tm = up->gyrokinetic_passive_tm,
  };
}

void
gkyl_dg_updater_gyrokinetic_passive_release(gkyl_dg_updater_gyrokinetic_passive *up)
{
  gkyl_dg_eqn_release(up->eqn_gyrokinetic_passive);
  gkyl_hyper_dg_release(up->up_gyrokinetic_passive);
  gkyl_free(up);
}
