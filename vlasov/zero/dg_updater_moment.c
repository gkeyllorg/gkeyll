#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_dg_updater_moment_priv.h>
#include <gkyl_mom_calc.h>
#include <gkyl_util.h>

struct gkyl_mom_type*
gkyl_dg_updater_moment_acquire_type(const gkyl_dg_updater_moment* moment)
{
  return gkyl_mom_type_acquire(moment->type);
}

int
gkyl_dg_updater_moment_num_mom(const gkyl_dg_updater_moment* moment)
{
  return gkyl_mom_type_num_mom(moment->type);
}

struct gkyl_dg_updater_moment*
gkyl_dg_updater_moment_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  const struct gkyl_range *phase_range,
  enum gkyl_model_id model_id, void *aux_inp, 
  enum gkyl_distribution_moments mom_type, bool is_integrated, bool use_gpu)
{
  gkyl_dg_updater_moment *up = gkyl_malloc(sizeof(gkyl_dg_updater_moment));
  up->model_id = model_id;
  up->use_gpu = use_gpu;
  if (up->model_id == GKYL_MODEL_SR) {
    if (is_integrated)
      up->type = gkyl_int_mom_vlasov_sr_new(cbasis, pbasis, conf_range, vel_range, mom_type, use_gpu);
    else
      up->type = gkyl_mom_vlasov_sr_new(cbasis, pbasis, conf_range, vel_range, mom_type, use_gpu);

    struct gkyl_mom_vlasov_sr_auxfields *sr_inp = aux_inp;
    gkyl_mom_vlasov_sr_set_auxfields(up->type, *sr_inp);

  } 
  else if ((up->model_id == GKYL_MODEL_CANONICAL_PB || up->model_id == GKYL_MODEL_CANONICAL_PB_GR) 
    && (mom_type == GKYL_F_MOMENT_M1_FROM_H || mom_type == GKYL_F_MOMENT_ENERGY
    || (is_integrated && mom_type == GKYL_F_MOMENT_M0M1M2))) {
    if (is_integrated)
      up->type = gkyl_int_mom_canonical_pb_new(cbasis, pbasis, phase_range, mom_type, use_gpu);
    else
      up->type = gkyl_mom_canonical_pb_new(cbasis, pbasis, phase_range, mom_type, use_gpu);
    
    struct gkyl_mom_canonical_pb_auxfields *can_pb_inp = aux_inp;
    gkyl_mom_canonical_pb_set_auxfields(up->type, *can_pb_inp);

  } 
  else {
    if (is_integrated)
      up->type = gkyl_int_mom_vlasov_new(cbasis, pbasis, mom_type, use_gpu);
    else
      up->type = gkyl_mom_vlasov_new(cbasis, pbasis, mom_type, use_gpu);
  }

  up->up_moment = gkyl_mom_calc_new(grid, up->type, use_gpu);

  up->moment_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_moment_advance(struct gkyl_dg_updater_moment *moment,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT mout)
{  
  struct timespec wst = gkyl_wall_clock();
  if (moment->use_gpu) {
    gkyl_mom_calc_advance_cu(moment->up_moment, update_phase_rng, update_conf_rng, fIn, mout);
  }
  else {
    gkyl_mom_calc_advance(moment->up_moment, update_phase_rng, update_conf_rng, fIn, mout);
  }
  moment->moment_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_moment_tm
gkyl_dg_updater_moment_get_tm(const gkyl_dg_updater_moment *moment)
{
  return (struct gkyl_dg_updater_moment_tm) {
    .moment_tm = moment->moment_tm,
  };
}

void
gkyl_dg_updater_moment_release(gkyl_dg_updater_moment* moment)
{
  gkyl_mom_type_release(moment->type);
  gkyl_mom_calc_release(moment->up_moment);
  gkyl_free(moment);
}
