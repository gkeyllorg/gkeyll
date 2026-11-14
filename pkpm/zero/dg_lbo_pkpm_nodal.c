#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_lbo_pkpm_nodal.h>
#include <gkyl_dg_lbo_pkpm_nodal_priv.h>
#include <gkyl_util.h>

void
gkyl_lbo_pkpm_nodal_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_lbo_pkpm_nodal *lbo = container_of(base, struct dg_lbo_pkpm_nodal, eqn);
  gkyl_free(lbo);
}

void
gkyl_lbo_pkpm_nodal_set_auxfields(
  const struct gkyl_dg_eqn *eqn, const struct gkyl_dg_lbo_pkpm_nodal_auxfields auxin
)
{
  struct dg_lbo_pkpm_nodal *lbo = container_of(eqn, struct dg_lbo_pkpm_nodal, eqn);
  lbo->auxfields.nu = auxin.nu;
  lbo->auxfields.prim_moms = auxin.prim_moms;
}

struct gkyl_dg_eqn *
gkyl_dg_lbo_pkpm_nodal_new(
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_rect_grid *pgrid, bool use_gpu
)
{
  assert(!use_gpu);
  assert(cbasis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && cbasis->poly_order == 1);
  assert(pbasis->b_type == GKYL_BASIS_MODAL_HYBRID);

  struct dg_lbo_pkpm_nodal *lbo = gkyl_malloc(sizeof(struct dg_lbo_pkpm_nodal));

  int cdim = cbasis->ndim, vdim = pbasis->ndim - cbasis->ndim;
  lbo->cdim = cdim;
  lbo->num_nodes = cbasis->num_basis;
  lbo->num_vel_basis = pbasis->num_basis / cbasis->num_basis;
  lbo->stride = pbasis->num_basis;
  lbo->vMaxSq = pgrid->upper[cdim] * pgrid->upper[cdim];

  lbo->eqn.num_equations = PKPM_DIST_NUM_BLOCKS;
  lbo->eqn.vol_term = vol;
  lbo->eqn.surf_term = surf;
  lbo->eqn.boundary_surf_term = boundary_surf;

  lbo->vol = lbo_pkpm_nodal_vol_kernels.kernels[vdim - 1];
  lbo->surf = lbo_pkpm_nodal_surf_vpar_kernels.kernels[vdim - 1];
  lbo->boundary_surf = lbo_pkpm_nodal_boundary_surf_vpar_kernels.kernels[vdim - 1];

  // ensure non-NULL pointers
  assert(lbo->vol);
  assert(lbo->surf);
  assert(lbo->boundary_surf);

  lbo->auxfields.nu = 0;
  lbo->auxfields.prim_moms = 0;
  lbo->conf_range = *conf_range;

  lbo->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(lbo->eqn.flags);
  lbo->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_pkpm_nodal_free);
  lbo->eqn.on_dev = &lbo->eqn;

  return &lbo->eqn;
}
