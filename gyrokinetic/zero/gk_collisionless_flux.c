#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_gk_collisionless_flux.h>
#include <gkyl_gk_collisionless_flux_priv.h>
#include <gkyl_util.h>

gkyl_gk_collisionless_flux*
gkyl_gk_collisionless_flux_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass,
  enum gkyl_gk_collisionless_type type,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom, 
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_gk_collisionless_flux_cu_dev_new(phase_grid, conf_basis, phase_basis,
      charge, mass, type, gk_geom, dg_geom, gk_dg_geom, vel_map, bctype_conf);
#endif     

  gkyl_gk_collisionless_flux *up = gkyl_malloc(sizeof(gkyl_gk_collisionless_flux));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;

  up->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  up->dg_geom = gkyl_dg_geom_acquire(dg_geom);
  up->gk_dg_geom = gkyl_gk_dg_geom_acquire(gk_dg_geom);
  up->vel_map = gkyl_velocity_map_acquire(vel_map);

  if (type == GKYL_GK_COLLISIONLESS_ES) {
    for (int d=0; d<cdim; ++d) {
      // BC option in ->flux_surf kernel doesn't matter as long as it's not SKIP.
      up->flux_surf[d] = choose_gk_collisionless_flux_surf_conf_kern(d, cdim, vdim, poly_order, GKYL_BC_GK_SPECIES_ABSORB);
      up->flux_surf_edge_lo[d] = choose_gk_collisionless_flux_surf_conf_kern(d, cdim, vdim,
        poly_order, bctype_conf[d]);
      up->flux_surf_edge_up[d] = choose_gk_collisionless_flux_edge_surf_conf_kern(d, cdim, vdim,
        poly_order, bctype_conf[GKYL_MAX_CDIM+d]);
    }
    up->flux_surfvpar[0] = choose_gk_collisionless_flux_surf_vpar_kern(cdim, vdim, poly_order);
  }
  else if (type == GKYL_GK_COLLISIONLESS_ES_NO_BY) {
    for (int d=0; d<cdim; ++d) {
      // BC option in ->flux_surf kernel doesn't matter as long as it's not SKIP.
      up->flux_surf[d] = choose_gk_collisionless_flux_no_by_surf_conf_kern(d, cdim, vdim, poly_order, GKYL_BC_GK_SPECIES_ABSORB);
      up->flux_surf_edge_lo[d] = choose_gk_collisionless_flux_no_by_surf_conf_kern(d, cdim, vdim,
        poly_order, bctype_conf[d]);
      up->flux_surf_edge_up[d] = choose_gk_collisionless_flux_no_by_edge_surf_conf_kern(d, cdim, vdim,
        poly_order, bctype_conf[GKYL_MAX_CDIM+d]);
    }
    up->flux_surfvpar[0] = choose_gk_collisionless_flux_no_by_surf_vpar_kern(cdim, vdim, poly_order);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_gk_collisionless_flux_surf(struct gkyl_gk_collisionless_flux *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array *fin, struct gkyl_array* flux_surf, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(flux_surf)) {
    return gkyl_gk_collisionless_flux_surf_cu(up, conf_range, phase_range,
      conf_ext_range, phase_ext_range, phi, fin, flux_surf, cflrate);
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM], idx_vel[2];
  int idxL[GKYL_MAX_DIM];
  int idx_velL[2];
  double xc[GKYL_MAX_DIM];

  struct gkyl_range vpar_range;
  int extend_lo[GKYL_MAX_DIM] = {0};
  int extend_up[GKYL_MAX_DIM] = {0};
  extend_lo[cdim] = -1;
  gkyl_range_extend(&vpar_range, phase_range, extend_lo, extend_up);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = iter.idx[d];

    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(phase_range, idx);

    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, loc_conf);
    const double *phi_d = gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double *flux_surf_d = gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d = gkyl_array_fetch(cflrate, loc_phase);

    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, idx, idxL);
      idxL[dir] = idx[dir] - 1;
      long loc_confL = gkyl_range_idx(conf_range, idxL);
      long locL = gkyl_range_idx(phase_range, idxL);
      const double *fL = gkyl_array_cfetch(fin, locL);
      const double *fR = gkyl_array_cfetch(fin, loc_phase);

      const struct gkyl_dg_surf_geom *dgs = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx);
      const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx);

      const double *jacgeo_rat_surfL_d = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_confL);
      const double *jacgeo_rat_surfR_d = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);

      if (idx[dir] == phase_range->lower[dir]) {
        // Lower domain/block boundary.
        cflrate_d[0] += up->flux_surf_edge_lo[dir](xc, up->phase_grid.dx, vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d, phi_d, fL, fR, flux_surf_d);
      } else {
        // Interior, lower cell surface.
        cflrate_d[0] += up->flux_surf[dir](xc, up->phase_grid.dx, vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d, phi_d, fL, fR, flux_surf_d);
      }

      // Upper domain/block boundary.
      // If the phase space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == phase_range->upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_conf_ext = gkyl_range_idx(conf_ext_range, idx_edge);
        long loc_phase_ext = gkyl_range_idx(phase_ext_range, idx_edge);

        double *cflrate_ext_d = gkyl_array_fetch(cflrate, loc_phase_ext);
        const double *fL = gkyl_array_cfetch(fin, loc_phase);
        const double *fR = gkyl_array_cfetch(fin, loc_phase_ext);
        const struct gkyl_dg_surf_geom *dgs = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx_edge);
        const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx_edge);

        const double *jacgeo_rat_surfL_d = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);
        const double *jacgeo_rat_surfR_d = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf_ext);

        double* flux_surf_ext_d = gkyl_array_fetch(flux_surf, loc_phase_ext);

        cflrate_ext_d[0] += up->flux_surf_edge_up[dir](xc, up->phase_grid.dx, vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d, phi_d, fL, fR, flux_surf_ext_d);
      }  
    }
  }

  gkyl_range_iter_init(&iter, &vpar_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = iter.idx[d];

    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(phase_range, idx);

    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, loc_conf);
    const double *phi_d = gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double *flux_surf_d = gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d = gkyl_array_fetch(cflrate, loc_phase);

    int dir=cdim;
    gkyl_copy_int_arr(pdim, idx, idxL);
    idxL[dir] = idx[dir] - 1;
    gkyl_copy_int_arr(pdim-cdim, idx_vel, idx_velL);
    idx_velL[0] = idx_velL[0]-1;
    long locL = gkyl_range_idx(phase_range, idxL);
    long loc_velL = gkyl_range_idx(&up->vel_map->local_vel, idx_velL);
    const double *fL = gkyl_array_cfetch(fin, locL);
    const double *fR = gkyl_array_cfetch(fin, loc_phase);

    const double *vpL = gkyl_array_cfetch(up->vel_map->vmap_prime, loc_velL);
    const double *vpR = gkyl_array_cfetch(up->vel_map->vmap_prime, loc_vel);

    const struct gkyl_dg_vol_geom *dgv = gkyl_dg_geom_get_vol(up->dg_geom, idx);
    const struct gkyl_gk_dg_vol_geom *gkdgv = gkyl_gk_dg_geom_get_vol(up->gk_dg_geom, idx);
    cflrate_d[0] += up->flux_surfvpar[0](xc, up->phase_grid.dx, vpL, vpR, vmap_d, vmapSq_d, up->charge, up->mass,
      dgv, gkdgv, bmag_d, phi_d,  fL, fR, flux_surf_d);
  }
}

void gkyl_gk_collisionless_flux_release(gkyl_gk_collisionless_flux *up)
{
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_dg_geom_release(up->dg_geom);
  gkyl_gk_dg_geom_release(up->gk_dg_geom);
  gkyl_velocity_map_release(up->vel_map);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
