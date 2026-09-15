/* -*- c++ -*- */

#include <math.h>
#include <time.h>
#include <float.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_collisionless_passive_flux.h>
#include <gkyl_gk_collisionless_passive_flux_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_gk_collisionless_passive_flux_surf_conf_cu_kernel(struct gkyl_gk_collisionless_passive_flux *up,
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_range conf_ext_range, struct gkyl_range phase_ext_range,
  const struct gkyl_array *fin,
  struct gkyl_array *flux_surf, struct gkyl_array *cflrate)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_ghost[GKYL_MAX_DIM], idx_vel[2];
  int idxL[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    for (int d = cdim; d < pdim; d++) idx_vel[d-cdim] = idx[d];

    long loc_conf  = gkyl_range_idx(&conf_range, idx);
    long loc_vel   = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d   = (const double*) gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, loc_conf);
    const double *vmap_d   = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double *flux_surf_d = (double*) gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d   = (double*) gkyl_array_fetch(cflrate, loc_phase);

    for (int dir = 0; dir < cdim; dir++) {
      gkyl_copy_int_arr(pdim, idx, idxL);
      idxL[dir] = idx[dir] - 1;
      long loc_confL = gkyl_range_idx(&conf_range, idxL);
      long locL      = gkyl_range_idx(&phase_range, idxL);
      const double *fL = (const double*) gkyl_array_cfetch(fin, locL);
      const double *fR = (const double*) gkyl_array_cfetch(fin, loc_phase);

      const struct gkyl_dg_surf_geom    *dgs   = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx);
      const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx);

      const double *jacgeo_rat_surfL_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_confL);
      const double *jacgeo_rat_surfR_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);

      const double *speeds_L = (const double*) gkyl_array_cfetch(up->passive_speeds, loc_confL);
      const double *speeds_R = (const double*) gkyl_array_cfetch(up->passive_speeds, loc_conf);

      if (idx[dir] == phase_range.lower[dir]) {
        // Lower domain/block boundary.
        cflrate_d[0] += up->flux_surf_edge_lo[dir](xc, up->phase_grid.dx,
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d,
          speeds_L, speeds_R, fL, fR, flux_surf_d);
      } else {
        // Interior, lower cell surface.
        cflrate_d[0] += up->flux_surf[dir](xc, up->phase_grid.dx,
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d,
          speeds_L, speeds_R, fL, fR, flux_surf_d);
      }

      // Upper domain boundary: compute the upper-edge surface expansion stored
      // in the ghost cell so we avoid reading geometry from ghost cells.
      if (idx[dir] == phase_range.upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_ghost);
        idx_ghost[dir] = idx_ghost[dir] + 1;
        long loc_conf_ghost  = gkyl_range_idx(&conf_ext_range, idx_ghost);
        long loc_phase_ghost = gkyl_range_idx(&phase_ext_range, idx_ghost);

        // Write into the skin cell's own cflrate (not the ghost cell's, which is excluded
        // from the CFL reduction range). Use a max instead of accumulating since cflrate_d
        // already holds this cell's lower-surface contribution from earlier in this dir loop.
        double *cflrate_ghost_d = (double*) gkyl_array_fetch(cflrate, loc_phase);
        const double *f_skin  = (const double*) gkyl_array_cfetch(fin, loc_phase);
        const double *f_ghost = (const double*) gkyl_array_cfetch(fin, loc_phase_ghost);

        const struct gkyl_dg_surf_geom    *dgs_ghost   = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx_ghost);
        const struct gkyl_gk_dg_surf_geom *gkdgs_ghost = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx_ghost);

        const double *jacgeo_rat_surf_skin  = (const double*) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);
        const double *jacgeo_rat_surf_ghost = (const double*) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf_ghost);

        const double *speeds_skin  = (const double*) gkyl_array_cfetch(up->passive_speeds, loc_conf);
        const double *speeds_ghost = (const double*) gkyl_array_cfetch(up->passive_speeds, loc_conf_ghost);

        double *flux_surf_ghost_d = (double*) gkyl_array_fetch(flux_surf, loc_phase_ghost);

        cflrate_ghost_d[0] = GKYL_MAX2(cflrate_ghost_d[0], up->flux_surf_edge_up[dir](xc, up->phase_grid.dx,
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs_ghost, gkdgs_ghost, bmag_d,
          jacgeo_rat_surf_skin, jacgeo_rat_surf_ghost,
          speeds_skin, speeds_ghost,
          f_skin, f_ghost, flux_surf_ghost_d));
      }
    }
    // No vpar loop: passive advection is conf-space only.
  }
}

// Host-side wrapper for the passive surface flux kernel.
void
gkyl_gk_collisionless_passive_flux_surf_cu(gkyl_gk_collisionless_passive_flux *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range,
  const struct gkyl_array *fin,
  struct gkyl_array *flux_surf, struct gkyl_array *cflrate)
{
  gkyl_gk_collisionless_passive_flux_surf_conf_cu_kernel<<<phase_range->volume, GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, *conf_range, *phase_range, *conf_ext_range, *phase_ext_range,
    fin->on_dev, flux_surf->on_dev, cflrate->on_dev);
}

// CUDA kernel to set device function pointers.
// Doing function-pointer work here avoids troublesome cudaMemcpyFromSymbol.
__global__ static void
gk_collisionless_passive_flux_set_cu_dev_ptrs(struct gkyl_gk_collisionless_passive_flux *up,
  int cdim, int vdim, int poly_order, const enum gkyl_gyrokinetic_bc_type *bctype_conf)
{
  for (int d = 0; d < cdim; d++) {
    up->flux_surf[d]         = choose_gk_collisionless_passive_flux_surf_conf_kern(d, cdim, vdim, poly_order);
    up->flux_surf_edge_lo[d] = choose_gk_collisionless_passive_flux_edge_lo_surf_conf_kern(d, cdim, vdim, poly_order, bctype_conf[d]);
    up->flux_surf_edge_up[d] = choose_gk_collisionless_passive_flux_edge_up_surf_conf_kern(d, cdim, vdim, poly_order, bctype_conf[GKYL_MAX_CDIM+d]);
  }
}

gkyl_gk_collisionless_passive_flux*
gkyl_gk_collisionless_passive_flux_cu_dev_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_array *passive_speeds,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom,
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf)
{
  gkyl_gk_collisionless_passive_flux *up =
    (gkyl_gk_collisionless_passive_flux*) gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass   = mass;

  // Store on_dev pointers so the memcpy below carries device-side array references.
  up->passive_speeds = passive_speeds->on_dev;

  // Acquire geometry/vel_map and use their on_dev pointers for the device copy.
  struct gk_geometry       *geom_ho        = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_dg_geom      *dg_geom_ho     = gkyl_dg_geom_acquire(dg_geom);
  struct gkyl_gk_dg_geom   *gk_dg_geom_ho  = gkyl_gk_dg_geom_acquire(gk_dg_geom);
  struct gkyl_velocity_map *vel_map_ho     = gkyl_velocity_map_acquire(vel_map);

  up->gk_geom    = geom_ho->on_dev;
  up->dg_geom    = dg_geom_ho->on_dev;
  up->gk_dg_geom = gk_dg_geom_ho->on_dev;
  up->vel_map    = vel_map_ho->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  // Temporarily copy array of BCs to device for kernel selection.
  enum gkyl_gyrokinetic_bc_type *bctype_conf_dev = (enum gkyl_gyrokinetic_bc_type *) gkyl_cu_malloc(2*GKYL_MAX_CDIM*sizeof(enum gkyl_gyrokinetic_bc_type));
  gkyl_cu_memcpy(bctype_conf_dev, bctype_conf, 2*GKYL_MAX_CDIM*sizeof(enum gkyl_gyrokinetic_bc_type), GKYL_CU_MEMCPY_H2D);

  // Copy host struct to device and select kernels.
  struct gkyl_gk_collisionless_passive_flux *up_cu =
    (struct gkyl_gk_collisionless_passive_flux*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);

  gk_collisionless_passive_flux_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, vdim, poly_order, bctype_conf_dev);

  gkyl_cu_free(bctype_conf_dev);

  up->on_dev = up_cu;

  // Restore host-side pointers for the host updater.
  up->passive_speeds = passive_speeds;

  up->gk_geom    = geom_ho;
  up->dg_geom    = dg_geom_ho;
  up->gk_dg_geom = gk_dg_geom_ho;
  up->vel_map    = vel_map_ho;

  return up;
}
