#include <gkyl_gyrokinetic_priv.h>

int
gk_stepped_moms_add(
  struct gk_stepped_moms *sm, bool use_gpu, long ncomp, long size, const struct gkyl_range *range,
  struct gkyl_array *rate
)
{
  int idx = sm->num;
  sm->num += 1;
  for (int s = 0; s < GK_STEP_STAGE_COUNT; s++) {
    sm->stage[s] = gkyl_realloc(sm->stage[s], sizeof(struct gkyl_array * [sm->num]));
    sm->stage[s][idx] = gkyl_array_clear(mkarr(use_gpu, ncomp, size), 0.0);
  }
  sm->rate = gkyl_realloc(sm->rate, sizeof(struct gkyl_array * [sm->num]));
  sm->range = gkyl_realloc(sm->range, sizeof(const struct gkyl_range * [sm->num]));
  sm->rate[idx] = rate ? gkyl_array_acquire(rate) : 0;
  sm->range[idx] = range;
  return idx;
}

void
gk_stepped_moms_forward_euler(
  const struct gk_stepped_moms *sm, enum gk_step_stage out, enum gk_step_stage in
)
{
  for (int i = 0; i < sm->num; i++) {
    struct gkyl_array *arr_out = sm->stage[out][i];
    if (sm->rate[i]) {
      // The rhs put the rate in the rate array, not in the output stage.
      gkyl_array_copy_range(arr_out, sm->rate[i], sm->range[i]);
    }
    gkyl_array_accumulate_range(arr_out, 1.0, sm->stage[in][i], sm->range[i]);
  }
}

void
gk_stepped_moms_set(
  const struct gk_stepped_moms *sm, enum gk_step_stage out, double fac, enum gk_step_stage in
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_set_range(sm->stage[out][i], fac, sm->stage[in][i], sm->range[i]);
  }
}

void
gk_stepped_moms_accumulate(
  const struct gk_stepped_moms *sm, enum gk_step_stage out, double fac, enum gk_step_stage in
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_accumulate_range(sm->stage[out][i], fac, sm->stage[in][i], sm->range[i]);
  }
}

void
gk_stepped_moms_release(const struct gk_stepped_moms *sm)
{
  for (int i = 0; i < sm->num; i++) {
    for (int s = 0; s < GK_STEP_STAGE_COUNT; s++) {
      gkyl_array_release(sm->stage[s][i]);
    }
    if (sm->rate[i]) {
      gkyl_array_release(sm->rate[i]);
    }
  }
  if (sm->num > 0) {
    for (int s = 0; s < GK_STEP_STAGE_COUNT; s++) {
      gkyl_free(sm->stage[s]);
    }
    gkyl_free(sm->rate);
    gkyl_free(sm->range);
  }
}
