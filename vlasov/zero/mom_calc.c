#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_priv.h>

#include <assert.h>

struct gkyl_mom_calc*
gkyl_mom_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_calc_cu_dev_new(grid, momt);
  } 
#endif
  gkyl_mom_calc *up = gkyl_malloc(sizeof(gkyl_mom_calc));
  
  up->grid = *grid;
  up->momt = gkyl_mom_type_acquire(momt);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

void
gkyl_mom_calc_advance(const struct gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout)
{
  double xc[GKYL_MAX_DIM];
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;
  
  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_array_clear_range(mout, 0.0, conf_rng);

  // the outer loop is over configuration space cells; for each
  // config-space cell the inner loop walks over the velocity space
  // computing the contribution to the moment
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);

    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);

    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_rng->ndim, phase_rng->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&calc->grid, pidx, xc);
      
      long fidx = gkyl_range_idx(&vel_rng, vel_iter.idx);

      gkyl_mom_type_calc(calc->momt, xc, calc->grid.dx, pidx,
        gkyl_array_cfetch(fin, fidx), gkyl_array_fetch(mout, midx), 0
      );
    }
  }
}

void gkyl_mom_calc_release(gkyl_mom_calc* up)
{
  gkyl_mom_type_release(up->momt);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_mom_calc_advance_cu(const struct gkyl_mom_calc* mcalc,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout)
{
  assert(false);
}

gkyl_mom_calc*
gkyl_mom_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt)
{
  assert(false);
}

#endif
