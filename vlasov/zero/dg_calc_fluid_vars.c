#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_fluid_vars.h>
#include <gkyl_dg_calc_fluid_vars_priv.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <gkyl_util.h>

gkyl_dg_calc_fluid_vars*
gkyl_dg_calc_fluid_vars_new(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  double limiter_fac, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_fluid_vars_cu_dev_new(wv_eqn, geom, cbasis, mem_range, limiter_fac);
  } 
#endif     
  gkyl_dg_calc_fluid_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_fluid_vars));

  up->eqn_type = wv_eqn->type;
  up->wv_eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->geom = gkyl_wave_geom_acquire(geom);  
  if (up->eqn_type == GKYL_EQN_EULER)
    up->param = gkyl_wv_euler_gas_gamma(up->wv_eqn);

  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->Ncomp = 3;
  up->mem_range = *mem_range;

  up->fluid_set = choose_fluid_set_kern(b_type, cdim, poly_order);
  up->fluid_copy = choose_fluid_copy_kern(b_type, cdim, poly_order);
  up->fluid_pressure = choose_fluid_pressure_kern(b_type, cdim, poly_order);
  up->fluid_ke = choose_fluid_ke_kern(b_type, cdim, poly_order);
  up->fluid_int = choose_fluid_int_kern(b_type, cdim, poly_order);
  up->fluid_source = choose_fluid_source_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    up->fluid_limiter[d] = choose_fluid_limiter_kern(d, b_type, cdim, poly_order);

  // Limiter factor for relationship between slopes and cell average differences
  // By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
  // and a cell slope estimate from two adjacent cells is (for the x variation): 
  // integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
  // where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
  // This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
  if (limiter_fac == 0.0)
    up->limiter_fac = 0.5773502691896258;
  else
    up->limiter_fac = limiter_fac;

  // There are Ncomp*range->volume linear systems to be solved 
  // 3 components: ux, uy, uz, 
  up->As = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_fluid_vars_advance(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_array* fluid, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* u, struct gkyl_array* u_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(u)) {
    return gkyl_dg_calc_fluid_vars_advance_cu(up, 
      fluid, cell_avg_prim, u, u_surf);
  }
#endif

  // First loop over mem_range for solving linear systems to compute primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *fluid_d = gkyl_array_cfetch(fluid, loc);

    int* cell_avg_prim_d = gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->fluid_set(count, up->As, up->xs, fluid_d);

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

    double* u_d = gkyl_array_fetch(u, loc);
    double* u_surf_d = gkyl_array_fetch(u_surf, loc);

    up->fluid_copy(count, up->xs, u_d, u_surf_d);

    count += up->Ncomp;
  }
}

void gkyl_dg_calc_fluid_vars_pressure(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(p)) {
    return gkyl_dg_calc_fluid_vars_pressure_cu(up, conf_range, 
      fluid, u, p, p_surf);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    const double *u_d = gkyl_array_cfetch(u, loc);

    double* p_d = gkyl_array_fetch(p, loc);
    double* p_surf_d = gkyl_array_fetch(p_surf, loc);

    up->fluid_pressure(up->param, fluid_d, u_d, p_d, p_surf_d);
  }
}

void gkyl_dg_calc_fluid_vars_ke(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* ke)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(ke)) {
    return gkyl_dg_calc_fluid_vars_ke_cu(up, conf_range, 
      fluid, u, ke);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    const double *u_d = gkyl_array_cfetch(u, loc);

    double* ke_d = gkyl_array_fetch(ke, loc);

    up->fluid_ke(fluid_d, u_d, ke_d);
  }
}

void gkyl_dg_calc_fluid_vars_limiter(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, struct gkyl_array* fluid)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid)) {
    return gkyl_dg_calc_fluid_vars_limiter_cu(up, conf_range, fluid);
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

    double *fluid_c = gkyl_array_fetch(fluid, linc);
    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      double *fluid_l = gkyl_array_fetch(fluid, linl);
      double *fluid_r = gkyl_array_fetch(fluid, linr);

      up->fluid_limiter[dir](up->limiter_fac, up->wv_eqn, geom, fluid_l, fluid_c, fluid_r);    
    }
  }
}

void gkyl_dg_calc_fluid_integrated_vars(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* fluid, 
  const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* fluid_int_vars)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid_int_vars)) {
    return gkyl_dg_calc_fluid_integrated_vars_cu(up, conf_range, 
      fluid, u_i, p_ij, fluid_int_vars);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    const double *u_i_d = gkyl_array_cfetch(u_i, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    
    double *fluid_int_vars_d = gkyl_array_fetch(fluid_int_vars, loc);
    up->fluid_int(fluid_d, u_i_d, p_ij_d, fluid_int_vars_d);
  }
}

void gkyl_dg_calc_fluid_vars_source(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* app_accel, const struct gkyl_array* fluid, 
  struct gkyl_array* rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_calc_fluid_vars_source_cu(up, conf_range, 
      app_accel, fluid, rhs);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *app_accel_d = gkyl_array_cfetch(app_accel, loc);
    const double *fluid_d = gkyl_array_cfetch(fluid, loc);

    double *rhs_d = gkyl_array_fetch(rhs, loc);
    up->fluid_source(app_accel_d, fluid_d, rhs_d);
  }
}

void gkyl_dg_calc_fluid_vars_release(gkyl_dg_calc_fluid_vars *up)
{
  gkyl_wv_eqn_release(up->wv_eqn);
  gkyl_wave_geom_release(up->geom);  

  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);
  
  if (GKYL_IS_CU_ALLOC(up->flags)) 
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
