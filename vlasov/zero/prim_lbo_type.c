#include <gkyl_prim_lbo_type.h>

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_type_acquire(const struct gkyl_prim_lbo_type* prim)
{
  gkyl_ref_count_inc(&prim->ref_count);
  return (struct gkyl_prim_lbo_type*) prim;
}

void
gkyl_prim_lbo_type_release(const struct gkyl_prim_lbo_type* prim)
{
  gkyl_ref_count_dec(&prim->ref_count);
}
