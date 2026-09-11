/* -*- c++ -*- */

#include <math.h>
#include <time.h>
#include <float.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_gk_collisionless_flux.h>
#include <gkyl_gk_collisionless_flux_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_gk_collisionless_flux_surf_conf_cu_kernel(struct gkyl_gk_collisionless_flux *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_range conf_ext_range, struct gkyl_range phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array *fin, struct gkyl_array* flux_surf, struct gkyl_array *cflrate)
{ 
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM], idx_vel[2];
  int idxL[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  // 2D thread grid
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = idx[d];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, loc_conf);
    const double *phi_d = (const double*) gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double* flux_surf_d = (double*) gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, loc_phase);

    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, idx, idxL);
      idxL[dir] = idx[dir] - 1;
      long loc_confL = gkyl_range_idx(&conf_range, idxL);
      long locL = gkyl_range_idx(&phase_range, idxL);
      const double *fL = (const double*) gkyl_array_cfetch(fin, locL);
      const double *fR = (const double*) gkyl_array_cfetch(fin, loc_phase);

      const struct gkyl_dg_surf_geom *dgs = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx);
      const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx);

      const double *jacgeo_rat_surfL_d = (const double *) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_confL);
      const double *jacgeo_rat_surfR_d = (const double *) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);

      if (idx[dir] == phase_range.lower[dir]) {
        // Lower domain/block boundary.
        cflrate_d[0] += up->flux_surf_edge_lo[dir](xc, up->phase_grid.dx, vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d, phi_d, fL, fR, flux_surf_d);
      } else {
        // Interior, lower cell surface.
        cflrate_d[0] += up->flux_surf[dir](xc, up->phase_grid.dx, vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d, phi_d, fL, fR, flux_surf_d);
      }

      // If the phase space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == phase_range.upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_conf_ext = gkyl_range_idx(&conf_ext_range, idx_edge);
        long loc_phase_ext = gkyl_range_idx(&phase_ext_range, idx_edge);

        double *cflrate_ext_d = (double*) gkyl_array_fetch(cflrate, loc_phase);
        const double *fL = (const double*)  gkyl_array_cfetch(fin, loc_phase);
        const double *fR = (const double*)  gkyl_array_cfetch(fin, loc_phase_ext);
        const struct gkyl_dg_surf_geom *dgs = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx_edge);
        const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx_edge);

        const double *jacgeo_rat_surfL_d = (const double *) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf);
        const double *jacgeo_rat_surfR_d = (const double *) gkyl_array_cfetch(up->gk_geom->geo_surf[dir].jacobgeo_ratio, loc_conf_ext);

        double* flux_surf_ext_d = (double*) gkyl_array_fetch(flux_surf, loc_phase_ext);

        // Write into the skin cell's own cflrate (not the ghost cell's, which is excluded
        // from the CFL reduction range). Use a max instead of accumulating since cflrate_d
        // already holds this cell's lower-surface contribution from earlier in this dir loop.
        cflrate_ext_d[0] = GKYL_MAX2(cflrate_ext_d[0], up->flux_surf_edge_up[dir](xc, up->phase_grid.dx, vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs, bmag_d, jacgeo_rat_surfL_d, jacgeo_rat_surfR_d, phi_d, fL, fR, flux_surf_ext_d));
      }  
    }
  }
}

__global__ void
gkyl_gk_collisionless_flux_surf_surfvpar_cu_kernel(struct gkyl_gk_collisionless_flux *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_range conf_ext_range, struct gkyl_range phase_ext_range, struct gkyl_range vpar_range,
  const struct gkyl_array *phi, const struct gkyl_array *fin, 
  struct gkyl_array* flux_surf, struct gkyl_array *cflrate)
{ 
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_vel[2];
  int idxL[GKYL_MAX_DIM];
  int idx_velL[2];
  double xc[GKYL_MAX_DIM];

  // 2D thread grid
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < vpar_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&vpar_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = idx[d];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_corn.bmag, loc_conf);
    const double *phi_d = (const double*) gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double* flux_surf_d = (double*) gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, loc_phase);

    int dir = cdim;
    gkyl_copy_int_arr(pdim, idx, idxL);
    idxL[dir] = idx[dir] - 1;
    gkyl_copy_int_arr(pdim-cdim, idx_vel, idx_velL);
    idx_velL[0] = idx_velL[0]-1;
    long locL = gkyl_range_idx(&phase_range, idxL);
    long loc_velL = gkyl_range_idx(&up->vel_map->local_vel, idx_velL);
    const double *fL = (const double*) gkyl_array_cfetch(fin, locL);
    const double *fR = (const double*) gkyl_array_cfetch(fin, loc_phase);

    const double *vpL = (const double*) gkyl_array_cfetch(up->vel_map->vmap_prime, loc_velL);
    const double *vpR = (const double*) gkyl_array_cfetch(up->vel_map->vmap_prime, loc_vel);

    const struct gkyl_dg_vol_geom *dgv = gkyl_dg_geom_get_vol(up->dg_geom, idx);
    const struct gkyl_gk_dg_vol_geom *gkdgv = gkyl_gk_dg_geom_get_vol(up->gk_dg_geom, idx);

    cflrate_d[0] += up->flux_surfvpar[0](xc, up->phase_grid.dx, 
      vpL, vpR,
      vmap_d, vmapSq_d, up->charge, up->mass,
      dgv, gkdgv, bmag_d, phi_d,  fL, fR, flux_surf_d);
  }
}

// Host-side wrapper for gyrokinetic surface alpha calculation
void gkyl_gk_collisionless_flux_surf_cu(struct gkyl_gk_collisionless_flux *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array *fin, struct gkyl_array* flux_surf, struct gkyl_array* cflrate)
{
  gkyl_gk_collisionless_flux_surf_conf_cu_kernel<<<phase_range->volume, GKYL_DEFAULT_NUM_THREADS>>>(up->on_dev, 
    *conf_range, *phase_range, *conf_ext_range, *phase_ext_range, phi->on_dev, fin->on_dev,
    flux_surf->on_dev, cflrate->on_dev);

  struct gkyl_range vpar_range;
  int sublower[GKYL_MAX_DIM];
  int subupper[GKYL_MAX_DIM];
  for(int i = 0; i < up->pdim; i++) {
     sublower[i] = phase_range->lower[i];
     subupper[i] = phase_range->upper[i];
  }
  sublower[up->cdim] += 1;
  gkyl_sub_range_init(&vpar_range, phase_ext_range, sublower, subupper);
  gkyl_gk_collisionless_flux_surf_surfvpar_cu_kernel<<<vpar_range.volume, GKYL_DEFAULT_NUM_THREADS>>>(up->on_dev, 
    *conf_range, *phase_range, *conf_ext_range, *phase_ext_range, vpar_range, phi->on_dev, fin->on_dev,
    flux_surf->on_dev, cflrate->on_dev);

}

// CUDA kernel to set device pointers to gyrokinetic vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
gk_collisionless_flux_set_cu_dev_ptrs(struct gkyl_gk_collisionless_flux *up, 
  int cdim, int vdim, int poly_order, enum gkyl_gk_collisionless_type type,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf)
{
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
}

gkyl_gk_collisionless_flux*
gkyl_gk_collisionless_flux_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass,
  enum gkyl_gk_collisionless_type type,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom, 
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf)
{
  struct gkyl_gk_collisionless_flux *up = (struct gkyl_gk_collisionless_flux*) gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_dg_geom *dg_geom_ho = gkyl_dg_geom_acquire(dg_geom);
  struct gkyl_gk_dg_geom *gk_dg_geom_ho = gkyl_gk_dg_geom_acquire(gk_dg_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  
  up->gk_geom = geom_ho->on_dev;
  up->dg_geom = dg_geom_ho->on_dev;
  up->gk_dg_geom = gk_dg_geom_ho->on_dev;
  up->vel_map = vel_map_ho->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  // Temporarily copy array of BCs to device for kernel selection.
  enum gkyl_gyrokinetic_bc_type *bctype_conf_dev = (enum gkyl_gyrokinetic_bc_type *) gkyl_cu_malloc(2*GKYL_MAX_CDIM*sizeof(enum gkyl_gyrokinetic_bc_type));
  gkyl_cu_memcpy(bctype_conf_dev, bctype_conf, 2*GKYL_MAX_CDIM*sizeof(enum gkyl_gyrokinetic_bc_type), GKYL_CU_MEMCPY_H2D);

  struct gkyl_gk_collisionless_flux *up_cu = (struct gkyl_gk_collisionless_flux*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_gk_collisionless_flux), GKYL_CU_MEMCPY_H2D);

  gk_collisionless_flux_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, vdim, poly_order, type, bctype_conf_dev);

  gkyl_cu_free(bctype_conf_dev);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  // Updater should store host pointers.
  up->gk_geom = geom_ho; 
  up->dg_geom = dg_geom_ho; 
  up->gk_dg_geom = gk_dg_geom_ho; 
  up->vel_map = vel_map_ho; 
  
  return up;
}
