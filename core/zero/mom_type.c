#include <gkyl_mom_type.h>
#include <gkyl_alloc_flags_priv.h>

bool
gkyl_mom_type_is_cu_dev(const struct gkyl_mom_type* momt)
{
  return GKYL_IS_CU_ALLOC(momt->flags);
}

struct gkyl_mom_type*
gkyl_mom_type_acquire(const struct gkyl_mom_type* momt)
{
  gkyl_ref_count_inc(&momt->ref_count);
  return (struct gkyl_mom_type*) momt;
}

void
gkyl_mom_type_release(const struct gkyl_mom_type* momt)
{
  gkyl_ref_count_dec(&momt->ref_count);
}

void
gkyl_mom_type_calc(const struct gkyl_mom_type* momt,
  const double *xc, const double *dx, const int *idx,
  const double *f, double* GKYL_RESTRICT out, void *param)
{
  momt->kernel(momt, xc, dx, idx, f, out, param);
}

int
gkyl_mom_type_num_mom(const struct gkyl_mom_type* momt)
{
  return momt->num_mom;
}
