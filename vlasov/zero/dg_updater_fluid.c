#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_canonical_pb_fluid.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_updater_fluid.h>
#include <gkyl_dg_updater_fluid_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn *
gkyl_dg_updater_fluid_acquire_eqn(const gkyl_dg_updater_fluid *fluid)
{
  return gkyl_dg_eqn_acquire(fluid->eqn_fluid);
}

gkyl_dg_updater_fluid *
gkyl_dg_updater_fluid_inew(const struct gkyl_dg_updater_fluid_inp *inp)
{
  gkyl_dg_updater_fluid *up = gkyl_malloc(sizeof(gkyl_dg_updater_fluid));
  up->eqn_id = inp->wv_eqn->type;
  up->use_gpu = inp->use_gpu;
  bool has_wall_flux = false;
  if (up->eqn_id == GKYL_EQN_ADVECTION) {
    up->eqn_fluid = gkyl_dg_advection_new(inp->cbasis, inp->conf_range, up->use_gpu);
    struct gkyl_dg_advection_auxfields *adv_in = inp->aux_inp;
    gkyl_advection_set_auxfields(up->eqn_fluid, *adv_in);
  } else if (up->eqn_id == GKYL_EQN_EULER || up->eqn_id == GKYL_EQN_ISO_EULER) {
    struct gkyl_dg_euler_auxfields *euler_in = inp->aux_inp;
    up->eqn_fluid = gkyl_dg_euler_inew(&(struct gkyl_dg_euler_inp){
      .cbasis = inp->cbasis,
      .conf_range = inp->conf_range,
      .wv_eqn = inp->wv_eqn,
      .geom = inp->geom,
      .pos_map = inp->pos_map,
      .u = euler_in->u,
      .p = euler_in->p,
      .use_gpu = up->use_gpu,
    });
    has_wall_flux = true;
  } else {
    up->eqn_fluid =
      gkyl_dg_canonical_pb_fluid_new(inp->cbasis, inp->conf_range, inp->wv_eqn, up->use_gpu);
    struct gkyl_dg_canonical_pb_fluid_auxfields *dg_can_pb_fluid_inp = inp->aux_inp;
    gkyl_canonical_pb_fluid_set_auxfields(up->eqn_fluid, *dg_can_pb_fluid_inp);
  }

  // The absorbing-wall edges take the boundary surface term of the equation (hyper_dg's
  // zero-flux flags); only the Euler equations have one.
  int cdim = inp->cbasis->ndim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2 * GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = has_wall_flux && inp->wall_flux[d];
    zero_flux_flags[d + cdim] = has_wall_flux && inp->wall_flux[d + cdim];
    assert(has_wall_flux || !(inp->wall_flux[d] || inp->wall_flux[d + cdim]));
  }
  int num_up_dirs = cdim;

  up->up_fluid = gkyl_hyper_dg_new(
    inp->grid, inp->cbasis, up->eqn_fluid, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu
  );

  up->fluid_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_fluid_advance(
  gkyl_dg_updater_fluid *fluid, const struct gkyl_range *update_rng,
  const struct gkyl_array *GKYL_RESTRICT fluidIn, struct gkyl_array *GKYL_RESTRICT cflrate,
  struct gkyl_array *GKYL_RESTRICT rhs
)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(fluid->up_fluid, update_rng, fluidIn, cflrate, rhs);
  fluid->fluid_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_fluid_tm
gkyl_dg_updater_fluid_get_tm(const gkyl_dg_updater_fluid *fluid)
{
  return (struct gkyl_dg_updater_fluid_tm){.fluid_tm = fluid->fluid_tm};
}

void
gkyl_dg_updater_fluid_release(gkyl_dg_updater_fluid *fluid)
{
  gkyl_dg_eqn_release(fluid->eqn_fluid);
  gkyl_hyper_dg_release(fluid->up_fluid);
  gkyl_free(fluid);
}
