#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_pkpm_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_pkpm_fluid *pkpm = container_of(base->on_dev, struct dg_pkpm_fluid, eqn);
    gkyl_cu_free(pkpm);
  }

  struct dg_pkpm_fluid *pkpm = container_of(base, struct dg_pkpm_fluid, eqn);
  gkyl_vlasov_position_map_release(pkpm->pos_map);
  gkyl_array_release(pkpm->vlasov_pkpm_moms);
  gkyl_array_release(pkpm->bvar);
  gkyl_array_release(pkpm->u);
  gkyl_array_release(pkpm->fluid_flux_surf);
  gkyl_free(pkpm);
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_fluid_inew(const struct gkyl_dg_pkpm_fluid_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_fluid_cu_dev_inew(inp);
  }
#endif
  struct dg_pkpm_fluid *pkpm = gkyl_malloc(sizeof(struct dg_pkpm_fluid));

  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;

  pkpm->cdim = cdim;
  pkpm->num_surf_basis = inp->conf_basis->num_basis / 2;

  pkpm->eqn.num_equations = PKPM_FLUID_NUM_COMP;
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;

  choose_pkpm_fluid_kerns(pkpm, inp->conf_basis->b_type, cdim, poly_order);

  pkpm->vlasov_pkpm_moms = gkyl_array_acquire(inp->vlasov_pkpm_moms);
  pkpm->bvar = gkyl_array_acquire(inp->bvar);
  pkpm->u = gkyl_array_acquire(inp->u);
  pkpm->fluid_flux_surf = gkyl_array_acquire(inp->fluid_flux_surf);
  pkpm->conf_range = *inp->conf_range;
  assert(inp->pos_map);
  pkpm->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  pkpm->jacob_pos = inp->pos_map->jacob_pos;

  pkpm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(pkpm->eqn.flags);

  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_fluid_free);
  pkpm->eqn.on_dev = &pkpm->eqn; // CPU eqn obj points to itself

  return &pkpm->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn *
gkyl_dg_pkpm_fluid_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_inp *inp)
{
  assert(false);
  return 0;
}

#endif
