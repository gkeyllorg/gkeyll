#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_sr.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_updater_vlasov_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_vlasov_acquire_eqn(const gkyl_dg_updater_vlasov* vlasov)
{
  return gkyl_dg_eqn_acquire(vlasov->eqn_vlasov);
}

struct gkyl_dg_updater_vlasov_tm
gkyl_dg_updater_vlasov_get_tm(const gkyl_dg_updater_vlasov *vlasov)
{
  return (struct gkyl_dg_updater_vlasov_tm) {
    .vlasov_tm = vlasov->vlasov_tm,
  };
}

gkyl_dg_updater_vlasov*
gkyl_dg_updater_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_bc, enum gkyl_model_id model_id, enum gkyl_field_id field_id, void *aux_inp, bool use_gpu)
{
  gkyl_dg_updater_vlasov *up = gkyl_malloc(sizeof(gkyl_dg_updater_vlasov));
  up->model_id = model_id;
  up->field_id = field_id;
  up->use_gpu = use_gpu;
  if (up->model_id == GKYL_MODEL_SR) {
    up->eqn_vlasov = gkyl_dg_vlasov_sr_new(cbasis, pbasis, conf_range, vel_range, up->field_id, up->use_gpu);
    struct gkyl_dg_vlasov_sr_auxfields *sr_inp = aux_inp;
    gkyl_vlasov_sr_set_auxfields(up->eqn_vlasov, *sr_inp);
  } 
  else if (up->model_id == GKYL_MODEL_CANONICAL_PB || up->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    up->eqn_vlasov = gkyl_dg_canonical_pb_new(cbasis, pbasis, phase_range, up->use_gpu);
    struct gkyl_dg_canonical_pb_auxfields *canonical_pb_inp = aux_inp;
    gkyl_canonical_pb_set_auxfields(up->eqn_vlasov, *canonical_pb_inp); 
  }
  else {
    up->eqn_vlasov = gkyl_dg_vlasov_new(cbasis, pbasis, conf_range, phase_range, up->model_id, up->field_id, up->use_gpu);
    struct gkyl_dg_vlasov_auxfields *vlasov_inp = aux_inp;
    gkyl_vlasov_set_auxfields(up->eqn_vlasov, *vlasov_inp); 
  }

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2*GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = is_zero_flux_bc[d]? 1 : 0;
    zero_flux_flags[d+pdim] = is_zero_flux_bc[d+pdim]? 1 : 0;
  }
  int num_up_dirs = cdim;
  // update velocity space only when field is present 
  // Need to include Canonical_pb to update velocity space directions
  if (field_id != GKYL_FIELD_NULL || up->model_id == GKYL_MODEL_CANONICAL_PB || up->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    for (int d=cdim; d<pdim; ++d) {
      up_dirs[d] = d;
      zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1; // zero-flux BCs in vel-space
    }
    num_up_dirs = pdim;
  }
  up->hdg_vlasov = gkyl_hyper_dg_new(grid, pbasis, up->eqn_vlasov, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->vlasov_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_vlasov_advance(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(vlasov->hdg_vlasov, update_rng, fIn, cflrate, rhs);
  vlasov->vlasov_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_dg_updater_vlasov_release(gkyl_dg_updater_vlasov* vlasov)
{
  gkyl_dg_eqn_release(vlasov->eqn_vlasov);
  gkyl_hyper_dg_release(vlasov->hdg_vlasov);
  gkyl_free(vlasov);
}
