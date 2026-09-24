#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_calc_bcorr_priv.h>
#include <gkyl_util.h>

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_mom_type *momt, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_mom_calc_bcorr_cu_dev_new(grid, momt);
  } 
#endif
  gkyl_mom_calc_bcorr *up = gkyl_malloc(sizeof(gkyl_mom_calc_bcorr));
  up->grid = *grid;
  up->momt = gkyl_mom_type_acquire(momt);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up;
  
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
gkyl_mom_calc_bcorr_advance(const struct gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out)
{
#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(bcorr->flags)) {
    gkyl_mom_calc_bcorr_advance_cu(bcorr, phase_rng, conf_rng, fIn, out);
    return;
  }
#endif

  double xc[GKYL_MAX_DIM];
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;
  
  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  enum gkyl_vel_edge edge;

  int pdim = phase_rng->ndim;
  int cdim = conf_rng->ndim;
  int vdim = pdim - cdim;
  
  for (int d=0; d<cdim; ++d) rem_dir[d] = 1;

  int defl_vdims[GKYL_MAX_VDIM*(GKYL_MAX_VDIM-1)] = { 0 };
  for (int d=0; d<vdim; ++d) {
    int c = 0;
    for (int i=0; i<vdim; ++i) {
      if (i != d)
        defl_vdims[d*(GKYL_MAX_VDIM-1)+c++] = i;
    }
  }

  gkyl_array_clear_range(out, 0.0, conf_rng);

  // outer loop is over configuration space cells; for each
  // config-space cell inner loop walks over the edges of velocity
  // space
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);
    
    for (int d=0; d<cdim; ++d)
      pidx[d] = conf_iter.idx[d];
    
    for (int d=0; d<vdim; ++d) {
      rem_dir[cdim + d] = 1;
      
      // loop over upper edge of velocity space
      edge = d + GKYL_MAX_CDIM;
      pidx[cdim + d] = phase_rng->upper[cdim + d];
      gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, pidx);
      gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
      
      while (gkyl_range_iter_next(&vel_iter)) {
        for (int i=0; i<vel_rng.ndim; ++i)
          pidx[cdim+defl_vdims[d*(GKYL_MAX_VDIM-1)+i]] = vel_iter.idx[i];

        gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);
      
        long fidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
        gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
          gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx), &edge
        );
      }
      
      // loop over lower edge of velocity space
      edge = d;
      pidx[cdim + d] = phase_rng->lower[cdim + d];
      gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, pidx);
      gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
      
      while (gkyl_range_iter_next(&vel_iter)) {
        for (int i=0; i<vel_rng.ndim; ++i)
          pidx[cdim+defl_vdims[d*(GKYL_MAX_VDIM-1)+i]] = vel_iter.idx[i];
  
        gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);
      
        long fidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
        gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
          gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx), &edge
        );
      }
      rem_dir[cdim + d] = 0;
      pidx[cdim + d] = 0;
    }
  }
}

void
gkyl_mom_calc_bcorr_release(gkyl_mom_calc_bcorr* up)
{
  gkyl_mom_type_release(up->momt);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_mom_type *momt)
{
  assert(false);
}

void
gkyl_mom_calc_bcorr_advance_cu(const struct gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out)
{
  assert(false);
}

#endif
