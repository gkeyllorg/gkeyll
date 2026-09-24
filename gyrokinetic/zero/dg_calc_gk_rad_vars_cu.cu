/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gk_rad_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_calc_gk_rad_vars_nu_advance_cu_kernel(struct gkyl_dg_calc_gk_rad_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  double a, double alpha, double beta, double gamma, double v0, 
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_vel[2];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = idx[d];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_int.bmag, loc_conf);

    double* vnu_surf_d = (double*) gkyl_array_fetch(vnu_surf, loc_phase);
    double* vnu_d = (double*) gkyl_array_fetch(vnu, loc_phase);
    double* vsqnu_surf_d = (double*) gkyl_array_fetch(vsqnu_surf, loc_phase);  
    double* vsqnu_d = (double*) gkyl_array_fetch(vsqnu, loc_phase);   
    const double *vmap_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    up->rad_nu_vpar(vmap_d, vmapSq_d, up->charge, up->mass, 
      a, alpha, beta, gamma, v0, bmag_d, vnu_surf_d, vnu_d);
    up->rad_nu_mu(vmap_d, vmapSq_d, up->charge, up->mass, 
      a, alpha, beta, gamma, v0, bmag_d, vsqnu_surf_d, vsqnu_d);
  }  
}

// Host-side wrapper for radiation drag coefficient calculation
void 
gkyl_dg_calc_gk_rad_vars_nu_advance_cu(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  double a, double alpha, double beta, double gamma, double v0, 
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_gk_rad_vars_nu_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, a, alpha, beta, gamma, v0,
    vnu_surf->on_dev, vnu->on_dev, vsqnu_surf->on_dev, vsqnu->on_dev);
}

__global__ void
gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu_kernel(struct gkyl_dg_calc_gk_rad_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_gk_rad_drag* vnu_surf, const struct gkyl_gk_rad_drag* vnu, 
  const struct gkyl_gk_rad_drag* vsqnu_surf, const struct gkyl_gk_rad_drag* vsqnu,
  const struct gkyl_array* n_elc_rad, const struct gkyl_array* n_elc,
  const struct gkyl_array *nI, 
  struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu,
  struct gkyl_array* vtsq_min_normalized, struct gkyl_array* vtsq)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM];
  
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_phase = gkyl_range_idx(&phase_range, idx);
    const double* vtsq_d = (const double*) gkyl_array_cfetch(vtsq, loc_conf);

    const double* ne = (const double*)gkyl_array_cfetch(n_elc, loc_conf);
    double ne_cell_avg = ne[0]/pow(2.0, cdim/2.0);
    
    // Find nearest index
    int ne_idx = gkyl_dg_rad_gyrokinetic_find_nearest_idx(n_elc_rad, ne_cell_avg);

    const double* vtsq_min_d = (const double*) gkyl_array_cfetch(vtsq_min_normalized, ne_idx);
    if ( vtsq_d[0] > vtsq_min_d[0] ) {            
      const double* vnu_surf_d = (const double*) gkyl_array_cfetch(vnu_surf[ne_idx].arr, loc_phase);
      const double* vnu_d = (const double*) gkyl_array_cfetch(vnu[ne_idx].arr, loc_phase);
      const double* vsqnu_surf_d = (const double*) gkyl_array_cfetch(vsqnu_surf[ne_idx].arr, loc_phase);  
      const double* vsqnu_d = (const double*) gkyl_array_cfetch(vsqnu[ne_idx].arr, loc_phase);   
      
      const double *nI_d = (const double*) gkyl_array_cfetch(nI, loc_conf);
      
      double* nvnu_surf_d = (double*) gkyl_array_fetch(nvnu_surf, loc_phase);
      double* nvnu_d = (double*) gkyl_array_fetch(nvnu, loc_phase);
      double* nvsqnu_surf_d = (double*) gkyl_array_fetch(nvsqnu_surf, loc_phase);  
      double* nvsqnu_d = (double*) gkyl_array_fetch(nvsqnu, loc_phase);   
      
      up->rad_nI_nu(vnu_surf_d, vnu_d, vsqnu_surf_d, vsqnu_d, nI_d, 
		    nvnu_surf_d, nvnu_d, nvsqnu_surf_d, nvsqnu_d);
    }
  }  
}

// Host-side wrapper for sum_s n_{i_s} nu_s(v) calculation for a given input n_{i_s} and nu_s(v)
void 
gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_gk_rad_drag* vnu_surf, const struct gkyl_gk_rad_drag* vnu, 
  const struct gkyl_gk_rad_drag* vsqnu_surf, const struct gkyl_gk_rad_drag* vsqnu,
  const struct gkyl_array* n_elc_rad, const struct gkyl_array* n_elc,
  const struct gkyl_array *nI, 
  struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu,
  struct gkyl_array* vtsq_min_normalized, struct gkyl_array* vtsq)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range,
    vnu_surf->on_dev, vnu->on_dev, vsqnu_surf->on_dev, vsqnu->on_dev, 
    n_elc_rad->on_dev, n_elc->on_dev, nI->on_dev, 
    nvnu_surf->on_dev, nvnu->on_dev, nvsqnu_surf->on_dev, nvsqnu->on_dev,
    vtsq_min_normalized->on_dev, vtsq->on_dev);
}

// CUDA kernel to set device pointers to gyrokinetic radiation vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_gk_rad_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_gk_rad_vars *up, 
  int cdim, int vdim, int poly_order)
{
  up->rad_nu_vpar = choose_rad_gyrokinetic_nu_vpar_kern(cdim, vdim, poly_order);
  up->rad_nu_mu = choose_rad_gyrokinetic_nu_mu_kern(cdim, vdim, poly_order);
  up->rad_nI_nu = choose_rad_gyrokinetic_nI_nu_kern(cdim, vdim, poly_order);
}

gkyl_dg_calc_gk_rad_vars*
gkyl_dg_calc_gk_rad_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, double charge,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map)
{
  struct gkyl_dg_calc_gk_rad_vars *up = (struct gkyl_dg_calc_gk_rad_vars*) gkyl_malloc(sizeof(*up));

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
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  up->gk_geom = geom_ho->on_dev;
  up->vel_map = vel_map_ho->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_gk_rad_vars *up_cu = (struct gkyl_dg_calc_gk_rad_vars*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_gk_rad_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_gk_rad_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, vdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  // Updater should store host pointers.
  up->gk_geom = geom_ho; 
  up->vel_map = vel_map_ho; 
  
  return up;
}
