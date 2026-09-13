#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_collisionless_passive_flux.h>
#include <gkyl_gk_collisionless_passive_flux_priv.h>
#include <gkyl_util.h>

gkyl_gk_collisionless_passive_flux*
gkyl_gk_collisionless_passive_flux_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_array *passive_speeds,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom,
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_gk_collisionless_passive_flux_cu_dev_new(phase_grid, conf_basis, phase_basis,
      passive_speeds, charge, mass, gk_geom, dg_geom, gk_dg_geom, vel_map, bctype_conf);
#endif

  gkyl_gk_collisionless_passive_flux *up = gkyl_malloc(sizeof(gkyl_gk_collisionless_passive_flux));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;

  // Store pointers to the projected passive-speed arrays (not owned here).
  up->passive_speeds = passive_speeds;

  up->gk_geom  = gkyl_gk_geometry_acquire(gk_geom);
  up->dg_geom  = gkyl_dg_geom_acquire(dg_geom);
  up->gk_dg_geom = gkyl_gk_dg_geom_acquire(gk_dg_geom);
  up->vel_map  = gkyl_velocity_map_acquire(vel_map);

  // Select kernels for each conf-space direction.
  for (int d = 0; d < cdim; d++) {
    up->flux_surf[d]         = choose_gk_collisionless_passive_flux_surf_conf_kern(d, cdim, vdim, poly_order);
    up->flux_surf_edge_lo[d] = choose_gk_collisionless_passive_flux_edge_lo_surf_conf_kern(d, cdim, vdim, poly_order, bctype_conf[d]);
    up->flux_surf_edge_up[d] = choose_gk_collisionless_passive_flux_edge_up_surf_conf_kern(d, cdim, vdim, poly_order, bctype_conf[GKYL_MAX_CDIM+d]);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_gk_collisionless_passive_flux_surf(gkyl_gk_collisionless_passive_flux *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range,
  const struct gkyl_array *fin,
  struct gkyl_array *flux_surf, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_gk_collisionless_passive_flux_surf_cu(up, conf_range, phase_range,
      conf_ext_range, phase_ext_range, fin, flux_surf, cflrate);
    return;
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_ghost[GKYL_MAX_DIM];
  int idxL[GKYL_MAX_DIM];
  int idx_vel[2];
  double xc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    for (int d = cdim; d < pdim; d++) idx_vel[d-cdim] = iter.idx[d];

    long loc_conf  = gkyl_range_idx(conf_range, idx);
    long loc_vel   = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(phase_range, idx);

    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *bmag_d    = gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, loc_conf);
    const double *vmap_d    = gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d  = gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double *flux_surf_d = gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d   = gkyl_array_fetch(cflrate, loc_phase);

    for (int dir = 0; dir < cdim; dir++) {
      gkyl_copy_int_arr(pdim, idx, idxL);
      idxL[dir] = idx[dir] - 1;
      long loc_confL  = gkyl_range_idx(conf_range, idxL);
      long locL       = gkyl_range_idx(phase_range, idxL);
      const double *fL = gkyl_array_cfetch(fin, locL);
      const double *fR = gkyl_array_cfetch(fin, loc_phase);

      const struct gkyl_dg_surf_geom    *dgs   = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx);
      const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx);

      const double *jacgeo_rat_surfL_d = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_confL);
      const double *jacgeo_rat_surfR_d = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);

      // Fetch passive speed DG coefficients at the left and right conf cells.
      const double *speeds_L = gkyl_array_cfetch(up->passive_speeds, loc_confL);
      const double *speeds_R = gkyl_array_cfetch(up->passive_speeds, loc_conf);

      if (idx[dir] == phase_range->lower[dir]) {
        // Lower domain boundary.
        cflrate_d[0] += up->flux_surf_edge_lo[dir](xc, up->phase_grid.dx,
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d,
          speeds_L, speeds_R, fL, fR, flux_surf_d);
      } else {
        // Interior lower surface.
        cflrate_d[0] += up->flux_surf[dir](xc, up->phase_grid.dx,
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d,
          speeds_L, speeds_R, fL, fR, flux_surf_d);
      }

      // Upper domain boundary: also compute the upper-edge surface expansion
      // (stored in the ghost cell) to avoid evaluating geometry in ghost cells.
      if (idx[dir] == phase_range->upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_ghost);
        idx_ghost[dir] += 1;
        long loc_conf_ghost  = gkyl_range_idx(conf_ext_range, idx_ghost);
        long loc_phase_ghost = gkyl_range_idx(phase_ext_range, idx_ghost);

        double *cflrate_ghost_d  = gkyl_array_fetch(cflrate, loc_phase_ghost);
        const double *f_skin  = gkyl_array_cfetch(fin, loc_phase);
        const double *f_ghost = gkyl_array_cfetch(fin, loc_phase_ghost);

        const struct gkyl_dg_surf_geom    *dgs_ghost   = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx_ghost);
        const struct gkyl_gk_dg_surf_geom *gkdgs_ghost = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx_ghost);

        const double *jacgeo_rat_surf_skin  = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);
        const double *jacgeo_rat_surf_ghost = gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf_ghost);

        // Passive speeds at the upper edge: right cell is the ghost cell,
        // so use the skin-cell speed for both sides (speeds are not defined in ghosts).
        const double *speeds_skin  = gkyl_array_cfetch(up->passive_speeds, loc_conf);
        const double *speeds_ghost = gkyl_array_cfetch(up->passive_speeds, loc_conf_ghost);

        double *flux_surf_ghost_d = gkyl_array_fetch(flux_surf, loc_phase_ghost);

        cflrate_ghost_d[0] += up->flux_surf_edge_up[dir](xc, up->phase_grid.dx,
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs_ghost, gkdgs_ghost, bmag_d,
          jacgeo_rat_surf_skin, jacgeo_rat_surf_ghost,
          speeds_skin, speeds_ghost,
          f_skin, f_ghost, flux_surf_ghost_d);

      }
    }
    // No vpar loop: passive advection has no velocity-space component.
  }
}

void
gkyl_gk_collisionless_passive_flux_release(gkyl_gk_collisionless_passive_flux *up)
{
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_dg_geom_release(up->dg_geom);
  gkyl_gk_dg_geom_release(up->gk_dg_geom);
  gkyl_velocity_map_release(up->vel_map);

  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
