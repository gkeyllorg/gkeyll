/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_cross_calc_priv.h>
#include <gkyl_prim_lbo_kernels.h> 
#include <gkyl_util.h>
#include <stdio.h>
}

__global__ static void
gkyl_prim_lbo_cross_calc_set_cu_ker(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_nmat *As, struct gkyl_nmat *xs,
  const struct gkyl_range conf_rng, const struct gkyl_array *alpha_E,
  double self_m, const struct gkyl_array *self_moms, const struct gkyl_array *self_prim_moms,
  double other_m, const struct gkyl_array *other_moms, const struct gkyl_array *other_prim_moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_rng.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_rng, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&conf_rng, idx);

    struct gkyl_mat lhs = gkyl_nmat_get(As, linc1);
    struct gkyl_mat rhs = gkyl_nmat_get(xs, linc1);

    const double *alpha_E_d = (const double*) gkyl_array_cfetch(alpha_E, start);
    const double *self_moms_d = (const double*) gkyl_array_cfetch(self_moms, start);
    const double *self_prim_moms_d = (const double*) gkyl_array_cfetch(self_prim_moms, start);
    const double *other_moms_d = (const double*) gkyl_array_cfetch(other_moms, start);
    const double *other_prim_moms_d = (const double*) gkyl_array_cfetch(other_prim_moms, start);
    const double *boundary_corrections_d = (const double*) gkyl_array_cfetch(boundary_corrections, start);
    const double *nu_d = (const double*) gkyl_array_cfetch(nu, start);

    gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

    calc->prim->cross_prim(calc->prim, &lhs, &rhs, idx, alpha_E_d, 
      self_m, self_moms_d, self_prim_moms_d,
      other_m, other_moms_d, other_prim_moms_d,
      boundary_corrections_d, nu_d
    );
  }
}

__global__ static void
gkyl_prim_lbo_copy_sol_cu_ker(struct gkyl_nmat *xs,
  const struct gkyl_range conf_rng,
  int nc, int udim, 
  struct gkyl_array* prim_moms_out)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_rng.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_rng, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&conf_rng, idx);

    struct gkyl_mat out_d = gkyl_nmat_get(xs, linc1);
    double *prim_moms_d = (double*) gkyl_array_fetch(prim_moms_out, start);
    
    prim_lbo_copy_sol(&out_d, nc, udim, prim_moms_d);
  }
}

void
gkyl_prim_lbo_cross_calc_advance_cu(struct gkyl_prim_lbo_cross_calc* calc,
  const struct gkyl_range *conf_rng, const struct gkyl_array *alpha_E,
  double self_m, const struct gkyl_array *self_moms, const struct gkyl_array *self_prim_moms,
  double other_m, const struct gkyl_array *other_moms, const struct gkyl_array *other_prim_moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu,  
  struct gkyl_array *prim_moms_out)
{
  int nc = calc->prim->num_config;
  int udim = calc->prim->udim;
  int N = nc*(udim + 1);
  
  if (calc->is_first) {
    calc->As = gkyl_nmat_cu_dev_new(conf_rng->volume, N, N);
    calc->xs = gkyl_nmat_cu_dev_new(conf_rng->volume, N, 1);
    calc->mem = gkyl_nmat_linsolve_lu_cu_dev_new(calc->As->num, calc->As->nr);
    calc->is_first = false;
  }
  
  gkyl_prim_lbo_cross_calc_set_cu_ker<<<conf_rng->nblocks, conf_rng->nthreads>>>(calc->on_dev,
    calc->As->on_dev, calc->xs->on_dev, *conf_rng, alpha_E->on_dev, 
    self_m, self_moms->on_dev, self_prim_moms->on_dev,
    other_m, other_moms->on_dev, other_prim_moms->on_dev,
    boundary_corrections->on_dev, nu->on_dev);
  
  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);
  
  gkyl_prim_lbo_copy_sol_cu_ker<<<conf_rng->nblocks, conf_rng->nthreads>>>(calc->xs->on_dev,
    *conf_rng, nc, udim, 
    prim_moms_out->on_dev);
}

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim)
{
  gkyl_prim_lbo_cross_calc *up = (gkyl_prim_lbo_cross_calc*) gkyl_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  up->grid = *grid;
  up->prim = prim;

  up->is_first = true;
  up->As = up->xs = 0;
  up->mem = 0;

  struct gkyl_prim_lbo_type *pt = gkyl_prim_lbo_type_acquire(prim);
  up->prim = pt->on_dev; // so memcpy below gets dev copy

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_prim_lbo_cross_calc *up_cu = (gkyl_prim_lbo_cross_calc*) gkyl_cu_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_prim_lbo_cross_calc), GKYL_CU_MEMCPY_H2D);

  up->prim = pt; // host portion of struct should have host copy
  up->on_dev = up_cu; // host pointer
  
  return up;
}
