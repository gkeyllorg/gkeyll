#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_em_vars*
gkyl_dg_calc_em_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  double limiter_fac, bool is_ExB, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_em_vars_cu_dev_new(conf_grid, cbasis, 
      mem_range, wv_eqn, geom, limiter_fac, is_ExB);
  } 
#endif     
  gkyl_dg_calc_em_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_em_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->mem_range = *mem_range;

  up->wv_eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->geom = gkyl_wave_geom_acquire(geom);  

  if (is_ExB) {
    up->Ncomp = 3;
    up->em_calc_temp = choose_em_calc_num_ExB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_ExB_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_ExB_kern(b_type, cdim, poly_order);
  }
  else {
    up->Ncomp = 6;
    up->em_calc_temp = choose_em_calc_BB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_bvar_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_bvar_kern(b_type, cdim, poly_order);      
    // Fetch the kernels in each direction
    for (int d=0; d<cdim; ++d) {
      up->em_div_b[d] = choose_em_div_b_kern(d, b_type, cdim, poly_order); 
      up->em_limiter[d] = choose_em_limiter_kern(d, b_type, cdim, poly_order);
    }
  }

  // Limiter factor for relationship between slopes and cell average differences
  // By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
  // and a cell slope estimate from two adjacent cells is (for the x variation): 
  // integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
  // where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
  // This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
  if (limiter_fac == 0.0) {
    up->limiter_fac = 0.5773502691896258;
  }
  else {
    up->limiter_fac = limiter_fac;
  }

  // There are Ncomp more linear systems to be solved 
  // 6 components of bb and 3 components of E x B
  up->As = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);
  // 6 component temporary variable for either storing B_i B_j (for computing bb) 
  // or (E x B)_i and B_i^2 (for computing E x B/|B|^2)
  up->temp_var = gkyl_array_new(GKYL_DOUBLE, 6*nc, mem_range->volume);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_em_vars_advance(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, 
  struct gkyl_array* out, struct gkyl_array* out_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_calc_em_vars_advance_cu(up, em, cell_avg_magB2, out, out_surf);
  }
#endif
  gkyl_array_clear(up->temp_var, 0.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);
    int* cell_avg_magB2_d = gkyl_array_fetch(cell_avg_magB2, loc);

    up->em_calc_temp(em_d, gkyl_array_fetch(up->temp_var, loc));
    cell_avg_magB2_d[0] = up->em_set(count, up->As, up->xs, gkyl_array_cfetch(up->temp_var, loc));

    count += up->Ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = gkyl_array_fetch(cell_avg_magB2, loc);
    double *out_d = gkyl_array_fetch(out, loc);
    double *out_surf_d = gkyl_array_fetch(out_surf, loc);

    up->em_copy(count, up->xs, em_d, cell_avg_magB2_d, out_d, out_surf_d);

    count += up->Ncomp;
  }  
}

void gkyl_dg_calc_em_vars_div_b(struct gkyl_dg_calc_em_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(div_b)) {
    return gkyl_dg_calc_em_vars_div_b_cu(up, conf_range, bvar_surf, bvar, max_b, div_b);
  }
#endif
  // Loop over configuration space range to compute div(b) and max(|b_i|) penalization
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);

    const double *bvar_surf_c = gkyl_array_cfetch(bvar_surf, linc);
    const double *bvar_d = gkyl_array_cfetch(bvar, linc);

    double *max_b_d = gkyl_array_fetch(max_b, linc);
    double *div_b_d = gkyl_array_fetch(div_b, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      const double *bvar_surf_l = gkyl_array_cfetch(bvar_surf, linl);
      const double *bvar_surf_r = gkyl_array_cfetch(bvar_surf, linr);

      up->em_div_b[dir](up->conf_grid.dx, 
        bvar_surf_l, bvar_surf_c, bvar_surf_r, 
        bvar_d, max_b_d, div_b_d);
    }
  }
}

void gkyl_dg_calc_em_vars_limiter(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_range *conf_range, struct gkyl_array* em)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(em)) {
    return gkyl_dg_calc_em_vars_limiter_cu(up, conf_range, em);
  }
#endif
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);
    const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(up->geom, idxc);

    double *em_c = gkyl_array_fetch(em, linc);
    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      double *em_l = gkyl_array_fetch(em, linl);
      double *em_r = gkyl_array_fetch(em, linr);

      up->em_limiter[dir](up->limiter_fac, up->wv_eqn, geom, em_l, em_c, em_r);    
    }
  }
}

void gkyl_dg_calc_em_vars_release(gkyl_dg_calc_em_vars *up)
{
  gkyl_wv_eqn_release(up->wv_eqn);
  gkyl_wave_geom_release(up->geom);

  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);
  gkyl_array_release(up->temp_var);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
