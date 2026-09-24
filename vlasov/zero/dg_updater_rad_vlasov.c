#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_rad_vlasov.h>
#include <gkyl_dg_updater_rad_vlasov_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_rad_vlasov*
gkyl_dg_updater_rad_vlasov_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, struct gkyl_dg_lbo_vlasov_drag_auxfields *drag_inp, 
  bool use_gpu)
{
  struct gkyl_dg_updater_rad_vlasov *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;

  up->rad_drag = gkyl_dg_lbo_vlasov_drag_new(conf_basis, phase_basis, conf_range, phase_grid, use_gpu);
  gkyl_lbo_vlasov_drag_set_auxfields(up->rad_drag, *drag_inp);

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + phase_basis->ndim - vdim;

  int zero_flux_flags[GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = 1;

  up->drag = gkyl_hyper_dg_new(phase_grid, phase_basis, up->rad_drag, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);  

  up->drag_tm = 0.0; 
  
  return up;
}

void
gkyl_dg_updater_rad_vlasov_advance(struct gkyl_dg_updater_rad_vlasov *rad,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(rad->drag, update_rng, fIn, cflrate, rhs);
  rad->drag_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_rad_vlasov_tm
gkyl_dg_updater_rad_vlasov_get_tm(const gkyl_dg_updater_rad_vlasov *rad)
{
  return (struct gkyl_dg_updater_rad_vlasov_tm) {
    .drag_tm = rad->drag_tm
  };
}

void
gkyl_dg_updater_rad_vlasov_release(gkyl_dg_updater_rad_vlasov *rad)
{
  gkyl_dg_eqn_release(rad->rad_drag);
  gkyl_hyper_dg_release(rad->drag);
  gkyl_free(rad);
}
