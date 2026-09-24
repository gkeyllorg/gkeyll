#include <gkyl_dg_eqn.h>
#include <gkyl_alloc_flags_priv.h>

bool
gkyl_dg_eqn_is_cu_dev(const struct gkyl_dg_eqn *eqn)
{
  return GKYL_IS_CU_ALLOC(eqn->flags);
}

struct gkyl_dg_eqn*
gkyl_dg_eqn_acquire(const struct gkyl_dg_eqn* eqn)
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_dg_eqn*) eqn;
}

void
gkyl_dg_eqn_release(const struct gkyl_dg_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
