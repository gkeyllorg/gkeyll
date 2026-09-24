#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_calc_priv.h>
#include <gkyl_prim_lbo_kernels.h>
#include <gkyl_mat.h>
#include <assert.h>

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_prim_lbo_calc_cu_dev_new(grid, prim);
  } 
#endif     
  gkyl_prim_lbo_calc *up = gkyl_malloc(sizeof(gkyl_prim_lbo_calc));
  up->grid = *grid;
  up->prim = gkyl_prim_lbo_type_acquire(prim);

  up->is_first = true;
  up->As = up->xs = 0;
  up->mem = 0;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void
gkyl_prim_lbo_calc_advance(struct gkyl_prim_lbo_calc* calc, 
  const struct gkyl_range *conf_rng, const struct gkyl_array *moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu,
  struct gkyl_array *prim_moms_out)
{
#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(calc->flags)) {
    gkyl_prim_lbo_calc_advance_cu(calc, conf_rng, moms, boundary_corrections, nu, prim_moms_out);
    return;
  }
#endif     

  struct gkyl_range_iter conf_iter;

  // allocate memory for use in kernels
  int nc = calc->prim->num_config;
  int udim = calc->prim->udim;
  int N = nc*(udim + 1);

  if (calc->is_first) {
    calc->As = gkyl_nmat_new(conf_rng->volume, N, N);
    calc->xs = gkyl_nmat_new(conf_rng->volume, N, 1);
    calc->mem = gkyl_nmat_linsolve_lu_new(calc->As->num, calc->As->nr);
    calc->is_first = false;
  }

  gkyl_array_clear_range(prim_moms_out, 0.0, conf_rng);

  // loop over configuration space cells.
  gkyl_range_iter_init(&conf_iter, conf_rng);
  long count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);

    struct gkyl_mat lhs = gkyl_nmat_get(calc->As, count);
    struct gkyl_mat rhs = gkyl_nmat_get(calc->xs, count);
    gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

    calc->prim->self_prim(calc->prim, &lhs, &rhs, conf_iter.idx,
      gkyl_array_cfetch(moms, midx), gkyl_array_cfetch(boundary_corrections, midx),
      gkyl_array_cfetch(nu, midx)
    );

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);
  assert(status);

  gkyl_range_iter_init(&conf_iter, conf_rng);
  count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);
    
    struct gkyl_mat out = gkyl_nmat_get(calc->xs, count);
    double *prim_moms = gkyl_array_fetch(prim_moms_out, midx);
    prim_lbo_copy_sol(&out, nc, udim, prim_moms);
    count += 1;
  }
}

const struct gkyl_prim_lbo_type* gkyl_prim_lbo_calc_get_prim(gkyl_prim_lbo_calc* calc)
{
  return calc->prim;
}

void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* up)
{
  gkyl_prim_lbo_type_release(up->prim);

  if (up->As)
    gkyl_nmat_release(up->As);
  if (up->xs)
    gkyl_nmat_release(up->xs);
  if (up->mem)
    gkyl_nmat_linsolve_lu_release(up->mem);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim)
{
  assert(false);
  return 0;
}

void
gkyl_prim_lbo_calc_advance_cu(struct gkyl_prim_lbo_calc* calc, 
  const struct gkyl_range *conf_rng,
  const struct gkyl_array* moms, const struct gkyl_array* boundary_corrections,
  struct gkyl_array* prim_moms_out)
{
  assert(false);
}

#endif
