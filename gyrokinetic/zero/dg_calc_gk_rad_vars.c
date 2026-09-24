#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gk_rad_vars_priv.h>
#include <gkyl_util.h>

struct gkyl_gk_rad_drag*
gkyl_dg_calc_gk_rad_vars_drag_new(int num_collisions,
  const int *num_densities, int ncomp, long sz, bool use_gpu)
{
  // Drag coefficient for each species.
  struct gkyl_gk_rad_drag *drag_s = gkyl_malloc(num_collisions*sizeof(struct gkyl_gk_rad_drag));
  for (int i=0; i<num_collisions; i++) {
    drag_s[i].num_dens = num_densities[i];
    drag_s[i].data = gkyl_malloc(num_densities[i]*sizeof(struct gkyl_gk_rad_drag));

    // Drag coefficient for each density.
    struct gkyl_gk_rad_drag *drag_ne = &drag_s[i];
    for (int n=0; n<num_densities[i]; n++) {
      struct gkyl_gk_rad_drag *drag = &drag_ne->data[n];
      drag->arr = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, sz)
                         : gkyl_array_new(GKYL_DOUBLE, ncomp, sz);
    }
  }

  // Create a temporary drag struct to store the on_dev pointers.
  struct gkyl_gk_rad_drag *drag_s_dev = gkyl_malloc(num_collisions*sizeof(struct gkyl_gk_rad_drag));
  for (int i=0; i<num_collisions; i++) {
    drag_s_dev[i].on_dev = gkyl_malloc(num_densities[i]*sizeof(struct gkyl_gk_rad_drag));

    struct gkyl_gk_rad_drag *drag_ne_dev = &drag_s_dev[i];
    struct gkyl_gk_rad_drag *drag_ne_ho = &drag_s[i];
    for (int n=0; n<num_densities[i]; n++) {
      struct gkyl_gk_rad_drag *drag_dev = &drag_ne_dev->on_dev[n];
      struct gkyl_gk_rad_drag *drag_ho = &drag_ne_ho->data[n];
      drag_dev->arr = drag_ho->arr->on_dev;
    }
  }

  // Now allocate the drag coeff as a function of density and assign its array to the on_dev pointers.
  for (int i=0; i<num_collisions; i++) {
    if (use_gpu) {
      drag_s[i].on_dev = gkyl_cu_malloc(num_densities[i]*sizeof(struct gkyl_gk_rad_drag));
      gkyl_cu_memcpy(drag_s[i].on_dev, drag_s_dev[i].on_dev,
        num_densities[i] * sizeof(struct gkyl_gk_rad_drag), GKYL_CU_MEMCPY_H2D);
    }
    else {
      drag_s[i].on_dev = gkyl_malloc(num_densities[i]*sizeof(struct gkyl_gk_rad_drag));
      memcpy(drag_s[i].on_dev, drag_s_dev[i].on_dev,
        num_densities[i] * sizeof(struct gkyl_gk_rad_drag));
    }
  }

  // Release the temporary struct with on_dev pointers.
  for (int i=0; i<num_collisions; i++) {
    gkyl_free(drag_s_dev[i].on_dev);
  }
  gkyl_free(drag_s_dev);
  
  return drag_s;
}

void
gkyl_dg_calc_gk_rad_vars_drag_release(struct gkyl_gk_rad_drag *drag_s, int num_collisions, bool use_gpu)
{
  // Free memory allocated to store a drag coefficient for each collision and each density.
  for (int i=0; i<num_collisions; i++) {
    if (use_gpu)
      gkyl_cu_free(drag_s[i].on_dev);
    else
      gkyl_free(drag_s[i].on_dev);
  }
  for (int i=0; i<num_collisions; i++) {
    struct gkyl_gk_rad_drag *drag_ne = &drag_s[i];
    for (int n=0; n<drag_s[i].num_dens; n++) {
      struct gkyl_gk_rad_drag *drag = &drag_ne->data[n];
      gkyl_array_release(drag->arr);
    }
    gkyl_free(drag_s[i].data);
  }
  gkyl_free(drag_s);
}

gkyl_dg_calc_gk_rad_vars* 
gkyl_dg_calc_gk_rad_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, double charge,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, bool use_gpu) 
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_dg_calc_gk_rad_vars_cu_dev_new(phase_grid, conf_basis, phase_basis, 
      charge, mass, gk_geom, vel_map);
  } 
#endif 

  gkyl_dg_calc_gk_rad_vars *up = gkyl_malloc(sizeof(*up));

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
  up->vel_map = gkyl_velocity_map_acquire(vel_map);

  up->rad_nu_vpar = choose_rad_gyrokinetic_nu_vpar_kern(cdim, vdim, poly_order);
  up->rad_nu_mu = choose_rad_gyrokinetic_nu_mu_kern(cdim, vdim, poly_order);
  up->rad_nI_nu = choose_rad_gyrokinetic_nI_nu_kern(cdim, vdim, poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_gk_rad_vars_nu_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  double a, double alpha, double beta, double gamma, double v0, 
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vnu_surf)) {
    return gkyl_dg_calc_gk_rad_vars_nu_advance_cu(up, conf_range, phase_range, 
      a, alpha, beta, gamma, v0, vnu_surf, vnu, vsqnu_surf, vsqnu);
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_vel[2];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = iter.idx[d];

    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(phase_range, idx);

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->geo_int.bmag, loc_conf);

    double* vnu_surf_d = gkyl_array_fetch(vnu_surf, loc_phase);
    double* vnu_d = gkyl_array_fetch(vnu, loc_phase);
    double* vsqnu_surf_d = gkyl_array_fetch(vsqnu_surf, loc_phase);  
    double* vsqnu_d = gkyl_array_fetch(vsqnu, loc_phase);   
    const double *vmap_d = gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    up->rad_nu_vpar(vmap_d, vmapSq_d, up->charge, up->mass, 
      a, alpha, beta, gamma, v0, bmag_d, vnu_surf_d, vnu_d);
    up->rad_nu_mu(vmap_d, vmapSq_d, up->charge, up->mass, 
      a, alpha, beta, gamma, v0, bmag_d, vsqnu_surf_d, vsqnu_d);
  }
}

void gkyl_dg_calc_gk_rad_vars_nI_nu_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_gk_rad_drag* vnu_surf, const struct gkyl_gk_rad_drag* vnu,
  const struct gkyl_gk_rad_drag* vsqnu_surf, const struct gkyl_gk_rad_drag* vsqnu,
  const struct gkyl_array* n_elc_rad, const struct gkyl_array* n_elc,
  const struct gkyl_array *nI, 
  struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu,
  struct gkyl_array* vtsq_min_normalized, struct gkyl_array* vtsq)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(nI)) {
    return gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu(up, conf_range, phase_range, 
      vnu_surf, vnu, vsqnu_surf, vsqnu, n_elc_rad, n_elc, nI, 
      nvnu_surf, nvnu, nvsqnu_surf, nvsqnu, vtsq_min_normalized, vtsq);
  }
#endif
  int pdim = up->pdim; // pdim and cdim are constant across densities
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_phase = gkyl_range_idx(phase_range, idx);
    const double* vtsq_d = gkyl_array_cfetch(vtsq, loc_conf);

    const double* ne = gkyl_array_cfetch(n_elc, loc_conf);
    double ne_cell_avg = ne[0]/pow(2.0, cdim/2.0);
    int ne_idx = gkyl_dg_rad_gyrokinetic_find_nearest_idx(n_elc_rad, ne_cell_avg);
    const double* vtsq_min_d = gkyl_array_cfetch(vtsq_min_normalized, ne_idx);
    if ( vtsq_d[0] > vtsq_min_d[0] ) {      
      const double* vnu_surf_d = gkyl_array_cfetch(vnu_surf->data[ne_idx].arr, loc_phase);
      const double* vnu_d = gkyl_array_cfetch(vnu->data[ne_idx].arr, loc_phase);
      const double* vsqnu_surf_d = gkyl_array_cfetch(vsqnu_surf->data[ne_idx].arr, loc_phase);  
      const double* vsqnu_d = gkyl_array_cfetch(vsqnu->data[ne_idx].arr, loc_phase);   
      const double *nI_d = gkyl_array_cfetch(nI, loc_conf);
      
      double* nvnu_surf_d = gkyl_array_fetch(nvnu_surf, loc_phase);
      double* nvnu_d = gkyl_array_fetch(nvnu, loc_phase);
      double* nvsqnu_surf_d = gkyl_array_fetch(nvsqnu_surf, loc_phase);  
      double* nvsqnu_d = gkyl_array_fetch(nvsqnu, loc_phase);   
      
      up->rad_nI_nu(vnu_surf_d, vnu_d, vsqnu_surf_d, vsqnu_d, nI_d, 
		    nvnu_surf_d, nvnu_d, nvsqnu_surf_d, nvsqnu_d);
    }
  }
}

void gkyl_dg_calc_gk_rad_vars_release(gkyl_dg_calc_gk_rad_vars *up)
{
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_velocity_map_release(up->vel_map);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
