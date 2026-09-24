#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_vlasov_poisson.h>
#include <gkyl_dg_updater_vlasov_poisson.h>
#include <gkyl_dg_updater_vlasov_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_vlasov_poisson_acquire_eqn(const gkyl_dg_updater_vlasov* vlasov)
{
  return gkyl_dg_eqn_acquire(vlasov->eqn_vlasov);
}

struct gkyl_dg_updater_vlasov_tm
gkyl_dg_updater_vlasov_poisson_get_tm(const gkyl_dg_updater_vlasov *vlasov)
{
  return (struct gkyl_dg_updater_vlasov_tm) {
    .vlasov_tm = vlasov->vlasov_tm,
  };
}

gkyl_dg_updater_vlasov*
gkyl_dg_updater_vlasov_poisson_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc, enum gkyl_model_id model_id, enum gkyl_field_id field_id, void *aux_inp, bool use_gpu)
{
  gkyl_dg_updater_vlasov *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;

  up->eqn_vlasov = gkyl_dg_vlasov_poisson_new(cbasis, pbasis, conf_range, phase_range, model_id, field_id, up->use_gpu);
  struct gkyl_dg_vlasov_poisson_auxfields *vlasov_inp = aux_inp;
  gkyl_vlasov_poisson_set_auxfields(up->eqn_vlasov, *vlasov_inp); 

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;

  int up_dirs[GKYL_MAX_DIM] = {0};
  int num_up_dirs = pdim;
  for (int d=0; d<num_up_dirs; ++d) up_dirs[d] = d;

  int zero_flux_flags[2*GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    zero_flux_flags[d] = is_zero_flux_bc[d]? 1 : 0;
    zero_flux_flags[d+pdim] = is_zero_flux_bc[d+pdim]? 1 : 0;
  }
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1; // Zero-flux BCs in v-space.

  up->hdg_vlasov = gkyl_hyper_dg_new(grid, pbasis, up->eqn_vlasov,
    num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->vlasov_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_vlasov_poisson_advance(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(vlasov->hdg_vlasov, update_rng, fIn, cflrate, rhs);
  vlasov->vlasov_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_dg_updater_vlasov_poisson_release(gkyl_dg_updater_vlasov* vlasov)
{
  gkyl_dg_eqn_release(vlasov->eqn_vlasov);
  gkyl_hyper_dg_release(vlasov->hdg_vlasov);
  gkyl_free(vlasov);
}
