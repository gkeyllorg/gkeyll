#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_dist_vars.h>
#include <gkyl_dg_calc_pkpm_dist_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_pkpm_dist_vars*
gkyl_dg_calc_pkpm_dist_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis* cbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_pkpm_dist_vars_cu_dev_new(phase_grid, cbasis);
  } 
#endif     
  gkyl_dg_calc_pkpm_dist_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_dist_vars));

  up->phase_grid = *phase_grid;
  enum gkyl_basis_type b_type = cbasis->b_type;
  int cdim = cbasis->ndim;
  up->cdim = cdim;
  int poly_order = cbasis->poly_order;

  up->pkpm_dist_mirror_force = choose_pkpm_dist_mirror_force_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    up->pkpm_dist_div_ppar[d] = choose_pkpm_dist_div_ppar_kern(d, b_type, cdim, poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_pkpm_dist_vars_mirror_force(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* pkpm_prim, const struct gkyl_array* nu_prim_moms_sum, 
  const struct gkyl_array* div_b, const struct gkyl_array* pkpm_accel, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(g_dist_source)) {
    return gkyl_dg_calc_pkpm_dist_vars_mirror_force_cu(up, 
      conf_range, phase_range, 
      pkpm_prim, nu_prim_moms_sum, div_b, pkpm_accel, 
      fIn, F_k_p_1, g_dist_source, F_k_m_1);
  }
#endif  
  // Cell center array
  double xc[GKYL_MAX_DIM];
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->phase_grid, iter.idx, xc);
    long loc_conf = gkyl_range_idx(conf_range, iter.idx);
    long loc_phase = gkyl_range_idx(phase_range, iter.idx);

    const double *pkpm_prim_d = gkyl_array_cfetch(pkpm_prim, loc_conf);
    const double *nu_prim_moms_sum_d = gkyl_array_cfetch(nu_prim_moms_sum, loc_conf);
    const double *div_b_d = gkyl_array_cfetch(div_b, loc_conf);
    const double *pkpm_accel_d = gkyl_array_cfetch(pkpm_accel, loc_conf);
    const double *fIn_d = gkyl_array_cfetch(fIn, loc_phase);
    const double *F_k_p_1_d = gkyl_array_cfetch(F_k_p_1, loc_phase);

    double *g_dist_source_d = gkyl_array_fetch(g_dist_source, loc_phase);
    double *F_k_m_1_d = gkyl_array_fetch(F_k_m_1, loc_phase);

    up->pkpm_dist_mirror_force(xc, up->phase_grid.dx, 
      pkpm_prim_d, nu_prim_moms_sum_d, div_b_d, pkpm_accel_d, 
      fIn_d, F_k_p_1_d, g_dist_source_d, F_k_m_1_d);
  }  
}

void gkyl_dg_calc_pkpm_dist_vars_div_ppar(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, const struct gkyl_array* fIn, 
  const struct gkyl_array* max_b, struct gkyl_array* pkpm_div_ppar)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_div_ppar)) {
    return gkyl_dg_calc_pkpm_dist_vars_div_ppar_cu(up, 
      conf_range, phase_range, 
      bvar_surf, bvar, fIn, 
      max_b, pkpm_div_ppar);
  }
#endif
  int cdim = up->cdim;
  // Cell center array
  double xc[GKYL_MAX_DIM];
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim+1, iter.idx, idxc);
    gkyl_rect_grid_cell_center(&up->phase_grid, idxc, xc);
    long loc_conf_c = gkyl_range_idx(conf_range, idxc);
    long loc_phase_c = gkyl_range_idx(phase_range, idxc);

    const double *bvar_surf_c = gkyl_array_cfetch(bvar_surf, loc_conf_c);
    const double *bvar_c = gkyl_array_cfetch(bvar, loc_conf_c);
    const double *f_c = gkyl_array_cfetch(fIn, loc_phase_c);
    const double *max_b_c = gkyl_array_cfetch(max_b, loc_conf_c);
    double *pkpm_div_ppar_d = gkyl_array_fetch(pkpm_div_ppar, loc_conf_c);
    
    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim+1, iter.idx, idxl);
      gkyl_copy_int_arr(cdim+1, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long loc_conf_l = gkyl_range_idx(conf_range, idxl);
      long loc_phase_l = gkyl_range_idx(phase_range, idxl);
      long loc_conf_r = gkyl_range_idx(conf_range, idxr);
      long loc_phase_r = gkyl_range_idx(phase_range, idxr);

      const double *bvar_surf_l = gkyl_array_cfetch(bvar_surf, loc_conf_l);
      const double *f_l = gkyl_array_cfetch(fIn, loc_phase_l);
      const double *bvar_surf_r = gkyl_array_cfetch(bvar_surf, loc_conf_r);
      const double *f_r = gkyl_array_cfetch(fIn, loc_phase_r);

      up->pkpm_dist_div_ppar[dir](xc, up->phase_grid.dx, 
        bvar_surf_l, bvar_surf_c, bvar_surf_r, 
        f_l, f_c, f_r, 
        bvar_c, max_b_c, pkpm_div_ppar_d);
    }    
  }  
}

void gkyl_dg_calc_pkpm_dist_vars_release(gkyl_dg_calc_pkpm_dist_vars *up)
{ 
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
