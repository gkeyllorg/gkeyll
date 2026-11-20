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
  gkyl_array_release(pkpm->c_0);
  gkyl_array_release(pkpm->alpha_0);
  gkyl_array_release(pkpm->c_h);
  gkyl_array_release(pkpm->a_perp);
  gkyl_array_release(pkpm->vel_flux_surf);
  gkyl_free(pkpm);
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

  pkpm->c_0 = gkyl_array_acquire(inp->c_0);
  pkpm->alpha_0 = gkyl_array_acquire(inp->alpha_0);
  pkpm->c_h = gkyl_array_acquire(inp->c_h);
  pkpm->a_perp = gkyl_array_acquire(inp->a_perp);
  pkpm->vel_flux_surf = gkyl_array_acquire(inp->vel_flux_surf);
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
