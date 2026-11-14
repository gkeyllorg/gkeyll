#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_kinetic_vel_priv.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_util.h>

void
gkyl_dg_pkpm_kinetic_vel_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_pkpm_kinetic_vel *pkpm = container_of(base->on_dev, struct dg_pkpm_kinetic_vel, eqn);
    gkyl_cu_free(pkpm);
  }

  struct dg_pkpm_kinetic_vel *pkpm = container_of(base, struct dg_pkpm_kinetic_vel, eqn);
  gkyl_free(pkpm);
}

void
gkyl_dg_pkpm_kinetic_vel_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_vel_auxfields auxin
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_dg_pkpm_kinetic_vel_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  pkpm->auxfields = auxin;
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_kinetic_vel_inew(const struct gkyl_dg_pkpm_kinetic_vel_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_kinetic_vel_cu_dev_inew(inp);
  }
#endif
  struct dg_pkpm_kinetic_vel *pkpm = gkyl_malloc(sizeof(struct dg_pkpm_kinetic_vel));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;
  assert(pdim == cdim + 1);

  pkpm->cdim = cdim;
  pkpm->pdim = pdim;

  pkpm->eqn.num_equations = PKPM_DIST_NUM_BLOCKS;
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;

  choose_pkpm_kinetic_vel_kerns(pkpm, inp->conf_basis->b_type, cdim, poly_order);

  pkpm->auxfields = inp->auxfields;
  pkpm->conf_range = *inp->conf_range;
  pkpm->phase_range = *inp->phase_range;

  pkpm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(pkpm->eqn.flags);

  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_kinetic_vel_free);
  pkpm->eqn.on_dev = &pkpm->eqn; // CPU eqn obj points to itself

  return &pkpm->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn *
gkyl_dg_pkpm_kinetic_vel_cu_dev_inew(const struct gkyl_dg_pkpm_kinetic_vel_inp *inp)
{
  assert(false);
  return 0;
}

#endif
