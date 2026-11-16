#include <gkyl_gyrokinetic_priv.h>

int
gk_stepped_moms_add(
  struct gk_stepped_moms *sm, bool use_gpu, long ncomp, long size, const struct gkyl_range *range
)
{
  int idx = sm->num;
  sm->num += 1;
  sm->f = gkyl_realloc(sm->f, sizeof(struct gkyl_array * [sm->num]));
  sm->f1 = gkyl_realloc(sm->f1, sizeof(struct gkyl_array * [sm->num]));
  sm->fnew = gkyl_realloc(sm->fnew, sizeof(struct gkyl_array * [sm->num]));
  sm->range = gkyl_realloc(sm->range, sizeof(const struct gkyl_range * [sm->num]));

  sm->f[idx] = gkyl_array_clear(mkarr(use_gpu, ncomp, size), 0.0);
  sm->f1[idx] = gkyl_array_clear(mkarr(use_gpu, ncomp, size), 0.0);
  sm->fnew[idx] = gkyl_array_clear(mkarr(use_gpu, ncomp, size), 0.0);
  sm->range[idx] = range;
  return idx;
}

void
gk_stepped_moms_step_f(
  const struct gk_stepped_moms *sm, struct gkyl_array **fout, double dt, struct gkyl_array **fin
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_accumulate_range(
      gkyl_array_scale_range(fout[i], dt, sm->range[i]), 1.0, fin[i], sm->range[i]
    );
  }
}

void
gk_stepped_moms_set(
  const struct gk_stepped_moms *sm, struct gkyl_array **fout, double fac, struct gkyl_array **fin
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_set_range(fout[i], fac, fin[i], sm->range[i]);
  }
}

void
gk_stepped_moms_accumulate(
  const struct gk_stepped_moms *sm, struct gkyl_array **fout, double fac, struct gkyl_array **fin
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_accumulate_range(fout[i], fac, fin[i], sm->range[i]);
  }
}

void
gk_stepped_moms_combine(
  const struct gk_stepped_moms *sm, struct gkyl_array **fout, double fac1, struct gkyl_array **fin1,
  double fac2, struct gkyl_array **fin2
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_accumulate_range(
      gkyl_array_set_range(fout[i], fac1, fin1[i], sm->range[i]), fac2, fin2[i], sm->range[i]
    );
  }
}

void
gk_stepped_moms_copy(
  const struct gk_stepped_moms *sm, struct gkyl_array **fout, struct gkyl_array **fin
)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_copy_range(fout[i], fin[i], sm->range[i]);
  }
}

void
gk_stepped_moms_release(const struct gk_stepped_moms *sm)
{
  for (int i = 0; i < sm->num; i++) {
    gkyl_array_release(sm->f[i]);
    gkyl_array_release(sm->f1[i]);
    gkyl_array_release(sm->fnew[i]);
  }
  if (sm->num > 0) {
    gkyl_free(sm->f);
    gkyl_free(sm->f1);
    gkyl_free(sm->fnew);
    gkyl_free(sm->range);
  }
}
