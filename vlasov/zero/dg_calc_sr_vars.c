#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_calc_sr_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_sr_vars*
gkyl_dg_calc_sr_vars_new(const struct gkyl_rect_grid *phase_grid, const struct gkyl_rect_grid *vel_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *vel_basis, 
  const struct gkyl_range *mem_range, const struct gkyl_range *vel_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_sr_vars_cu_dev_new(phase_grid, vel_grid, 
      conf_basis, vel_basis, mem_range, vel_range);
  } 
#endif     
  gkyl_dg_calc_sr_vars *up = gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  up->vel_grid = *vel_grid;
  up->vel_range = *vel_range;

  int nc = conf_basis->num_basis;
  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  enum gkyl_basis_type b_type = conf_basis->b_type;
  // store polynomial order and mem_range for linear solve
  up->poly_order = poly_order;
  up->mem_range = *mem_range;

  int vdim = vel_basis->ndim;
  int poly_order_v = vel_basis->poly_order;
  enum gkyl_basis_type b_type_v = vel_basis->b_type;

  up->sr_p_vars = choose_sr_p_vars_kern(b_type_v, vdim, poly_order_v);
  up->sr_n_set = choose_sr_vars_n_set_kern(b_type, cdim, vdim, poly_order);
  up->sr_n_copy = choose_sr_vars_n_copy_kern(b_type, cdim, vdim, poly_order);
  up->sr_GammaV = choose_sr_vars_GammaV_kern(b_type, cdim, vdim, poly_order);
  up->sr_pressure = choose_sr_vars_pressure_kern(b_type, cdim, vdim, poly_order);

  // Linear system for solving for the drift velocity V_drift = M1i/M0 
  // and then computing the rest-frame density n = GammaV_inv*M0 
  // where GammaV_inv = sqrt(1 - |V_drift|^2)
  up->Ncomp = vdim; 
  up->As = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_calc_sr_vars_init_p_vars(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(gamma)) {
    return gkyl_calc_sr_vars_init_p_vars_cu(up, gamma, gamma_inv);
  }
#endif

  // Cell center array
  double xc[GKYL_MAX_DIM];  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->vel_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->vel_grid, iter.idx, xc);
    long loc = gkyl_range_idx(&up->vel_range, iter.idx);

    double *gamma_d = gkyl_array_fetch(gamma, loc);
    double *gamma_inv_d = gkyl_array_fetch(gamma_inv, loc);
    up->sr_p_vars(xc, up->vel_grid.dx, gamma_d, gamma_inv_d);
  }
}

void gkyl_dg_calc_sr_vars_n(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_array* M0, const struct gkyl_array* M1i, struct gkyl_array* n)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(n)) {
    return gkyl_dg_calc_sr_vars_n_cu(up, M0, M1i, n);
  }
#endif

  // First loop over mem_range for setting matrices to solve linear systems for V_drift
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *M0_d = gkyl_array_cfetch(M0, loc);
    const double *M1i_d = gkyl_array_cfetch(M1i, loc);

    up->sr_n_set(count, up->As, up->xs, M0_d, M1i_d);

    count += up->Ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  // Then loop over mem_range to construct 1/Gamma = sqrt(1 - V_drift^2)
  // to solve for the rest-frame density n = M0/Gamma
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *M0_d = gkyl_array_cfetch(M0, loc);
    double* n_d = gkyl_array_fetch(n, loc);

    up->sr_n_copy(count, up->xs, M0_d, n_d);

    count += up->Ncomp;
  }
}

void gkyl_dg_calc_sr_vars_GammaV(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range,
  const struct gkyl_array* u_i, struct gkyl_array* u_i_sq, 
  struct gkyl_array* GammaV, struct gkyl_array* GammaV_sq)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(GammaV)) {
    return gkyl_dg_calc_sr_vars_GammaV_cu(up, conf_range, 
      u_i, u_i_sq, GammaV, GammaV_sq);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc_conf = gkyl_range_idx(conf_range, iter.idx);

    const double *u_i_d = gkyl_array_cfetch(u_i, loc_conf);
    double *u_i_sq_d = gkyl_array_fetch(u_i_sq, loc_conf);
    double *GammaV_d = gkyl_array_fetch(GammaV, loc_conf);
    double *GammaV_sq_d = gkyl_array_fetch(GammaV_sq, loc_conf);

    up->sr_GammaV(u_i_d, u_i_sq_d, GammaV_d, GammaV_sq_d); 
  }
}

void gkyl_dg_calc_sr_vars_pressure(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* gamma, const struct gkyl_array* gamma_inv, 
  const struct gkyl_array* u_i, const struct gkyl_array* u_i_sq, 
  const struct gkyl_array* GammaV, const struct gkyl_array* GammaV_sq, 
  const struct gkyl_array* f, struct gkyl_array* sr_pressure)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(sr_pressure)) {
    return gkyl_dg_calc_sr_vars_pressure_cu(up, 
      conf_range, phase_range, 
      gamma, gamma_inv, u_i, u_i_sq, GammaV, GammaV_sq, 
      f, sr_pressure);
  }
#endif
  gkyl_array_clear(sr_pressure, 0.0); 

  int cdim = conf_range->ndim;
  int pdim = phase_range->ndim;
  int idx_vel[GKYL_MAX_DIM];
  // Cell center array
  double xc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->phase_grid, iter.idx, xc);
    long loc_conf = gkyl_range_idx(conf_range, iter.idx);
    long loc_phase = gkyl_range_idx(phase_range, iter.idx);

    for (int i=0; i<pdim-cdim; ++i) {
      idx_vel[i] = iter.idx[cdim+i];
    }
    long loc_vel = gkyl_range_idx(&up->vel_range, idx_vel);

    const double *gamma_d = gkyl_array_cfetch(gamma, loc_vel);
    const double *gamma_inv_d = gkyl_array_cfetch(gamma_inv, loc_vel);
    const double *u_i_d = gkyl_array_cfetch(u_i, loc_conf);
    const double *u_i_sq_d = gkyl_array_cfetch(u_i_sq, loc_conf);
    const double *GammaV_d = gkyl_array_cfetch(GammaV, loc_conf);
    const double *GammaV_sq_d = gkyl_array_cfetch(GammaV_sq, loc_conf);
    const double *f_d = gkyl_array_cfetch(f, loc_phase);

    double *sr_pressure_d = gkyl_array_fetch(sr_pressure, loc_conf);

    up->sr_pressure(xc, up->phase_grid.dx, 
      gamma_d, gamma_inv_d, u_i_d, u_i_sq_d, GammaV_d, GammaV_sq_d, 
      f_d, sr_pressure_d);   
  }  
}

void gkyl_dg_calc_sr_vars_release(gkyl_dg_calc_sr_vars *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);

  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
