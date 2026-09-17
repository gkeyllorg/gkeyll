#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gyrokinetic_passive.h>
#include <gkyl_dg_gyrokinetic_passive_priv.h>
#include <gkyl_util.h>

void
gkyl_gyrokinetic_passive_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_gyrokinetic_passive *gkp = container_of(base, struct dg_gyrokinetic_passive, eqn);
  gkyl_gk_geometry_release(gkp->gk_geom);
  gkyl_velocity_map_release(gkp->vel_map);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_gyrokinetic_passive *gkp_cu = container_of(base->on_dev, struct dg_gyrokinetic_passive, eqn);
    gkyl_cu_free(gkp_cu);
  }

  gkyl_free(gkp);
}

void
gkyl_gyrokinetic_passive_set_auxfields(const struct gkyl_dg_eqn *eqn,
  struct gkyl_dg_gyrokinetic_passive_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_gyrokinetic_passive_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  struct dg_gyrokinetic_passive *gkp = container_of(eqn, struct dg_gyrokinetic_passive, eqn);
  gkp->auxfields.flux_surf = auxin.flux_surf;
  gkp->auxfields.speeds    = auxin.speeds;
}

struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_passive_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map,
  bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_dg_gyrokinetic_passive_cu_dev_new(cbasis, pbasis, conf_range, phase_range,
      charge, mass, gk_geom, vel_map);
#endif

  struct dg_gyrokinetic_passive *gkp = gkyl_malloc(sizeof(struct dg_gyrokinetic_passive));

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = cbasis->poly_order;

  gkp->cdim = cdim;
  gkp->pdim = pdim;

  gkp->charge = charge;
  gkp->mass   = mass;

  gkp->eqn.num_equations      = 1;
  gkp->eqn.surf_term          = surf_passive;
  gkp->eqn.boundary_surf_term = boundary_surf_passive;
  gkp->eqn.boundary_diag_term = boundary_diag_passive;

  // Volume kernel.
  gkp->eqn.vol_term = CK_PASSIVE(ser_passive_vol_kernels, cdim, vdim, poly_order);

  // Conf-space surface and boundary-surface kernels.
  if (cdim >= 1) {
    gkp->surf[0]          = CK_PASSIVE(ser_passive_surf_x_kernels,          cdim, vdim, poly_order);
    gkp->boundary_surf[0] = CK_PASSIVE(ser_passive_boundary_surf_x_kernels, cdim, vdim, poly_order);
  }
  if (cdim >= 2) {
    gkp->surf[1]          = CK_PASSIVE(ser_passive_surf_y_kernels,          cdim, vdim, poly_order);
    gkp->boundary_surf[1] = CK_PASSIVE(ser_passive_boundary_surf_y_kernels, cdim, vdim, poly_order);
  }
  if (cdim >= 3) {
    gkp->surf[2]          = CK_PASSIVE(ser_passive_surf_z_kernels,          cdim, vdim, poly_order);
    gkp->boundary_surf[2] = CK_PASSIVE(ser_passive_boundary_surf_z_kernels, cdim, vdim, poly_order);
  }

  // Ensure non-NULL pointers.
  for (int i=0; i<cdim; ++i) {
    assert(gkp->surf[i]);
    assert(gkp->boundary_surf[i]);
  }

  gkp->conf_range  = *conf_range;
  gkp->phase_range = *phase_range;

  gkp->gk_geom  = gkyl_gk_geometry_acquire(gk_geom);
  gkp->vel_map  = gkyl_velocity_map_acquire(vel_map);

  gkp->auxfields.flux_surf = 0;
  gkp->auxfields.speeds    = 0;

  gkp->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gkp->eqn.flags);
  gkp->eqn.ref_count = gkyl_ref_count_init(gkyl_gyrokinetic_passive_free);
  gkp->eqn.on_dev    = &gkp->eqn;

  return &gkp->eqn;
}
