#include <gkyl_dist_type.h>
#include <gkyl_alloc_flags_priv.h>

bool
gkyl_dist_proj_type_is_cu_dev(const struct gkyl_dist_proj_type *dist_proj)
{
  return GKYL_IS_CU_ALLOC(dist_proj->flags);
}

struct gkyl_dist_proj_type*
gkyl_dist_proj_type_acquire(const struct gkyl_dist_proj_type *dist_proj)
{
  gkyl_ref_count_inc(&dist_proj->ref_count);
  return (struct gkyl_dist_proj_type*) dist_proj;
}

void
gkyl_dist_proj_type_release(const struct gkyl_dist_proj_type *dist_proj)
{
  gkyl_ref_count_dec(&dist_proj->ref_count);
}


void
gkyl_dist_proj_type_calc(const struct gkyl_dist_proj_type *dist_proj,
  const double *xc, const double *dx,
  const int *idx, const double *moms,
  double *GKYL_RESTRICT out)
{
  dist_proj->proj_kernel(dist_proj, xc, dx, idx, moms, out);
}

int
gkyl_dist_proj_type_num_mom(const struct gkyl_dist_proj_type *dist_proj)
{
  return dist_proj->num_mom;
}