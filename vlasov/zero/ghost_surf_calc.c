#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_ghost_surf_calc_priv.h>
#include <gkyl_util.h>

void
gkyl_ghost_surf_calc_advance(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs)
{
  // Ghost and skin index and cell center coordinates.
  int idxg[GKYL_MAX_DIM], idxs[GKYL_MAX_DIM];
  double xcg[GKYL_MAX_DIM], xcs[GKYL_MAX_DIM];

  struct gkyl_range edge_rng;
  struct gkyl_range_iter edge_iter;
  
  int clower_idx[GKYL_MAX_DIM], cupper_idx[GKYL_MAX_DIM] = { 0 };

  for (int d=0; d<phase_rng->ndim; ++d) {
    clower_idx[d] = phase_rng->lower[d];
    cupper_idx[d] = phase_rng->upper[d];
  }

  for (int dir=0; dir<gcalc->cdim; ++dir) {
    // Ghost surf at lower boundary.
    clower_idx[dir] = phase_rng->lower[dir];
    cupper_idx[dir] = phase_rng->lower[dir];
    gkyl_sub_range_init(&edge_rng, phase_rng, clower_idx, cupper_idx);
    gkyl_range_iter_no_split_init(&edge_iter, &edge_rng);

    while (gkyl_range_iter_next(&edge_iter)) {
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxg);
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxs);
      idxs[dir] += 1;

      gkyl_rect_grid_cell_center(&gcalc->grid, idxg, xcg);
      gkyl_rect_grid_cell_center(&gcalc->grid, idxs, xcs);

      long ling = gkyl_range_idx(&edge_rng, idxg); 
      long lins = gkyl_range_idx(&edge_rng, idxs);

      gcalc->equation->boundary_surf_term(gcalc->equation, dir, xcs, xcg,
        gcalc->grid.dx, gcalc->grid.dx, idxs, idxg, -1,
	gkyl_array_cfetch(fIn, lins), gkyl_array_cfetch(fIn, ling), gkyl_array_fetch(rhs, ling)
      );
    }

    // Ghost surf at upper boundary.
    clower_idx[dir] = phase_rng->upper[dir];
    cupper_idx[dir] = phase_rng->upper[dir];
    gkyl_sub_range_init(&edge_rng, phase_rng, clower_idx, cupper_idx);
    gkyl_range_iter_no_split_init(&edge_iter, &edge_rng);

    while (gkyl_range_iter_next(&edge_iter)) {
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxs);
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxg);
      idxs[dir] -= 1;

      gkyl_rect_grid_cell_center(&gcalc->grid, idxs, xcs);
      gkyl_rect_grid_cell_center(&gcalc->grid, idxg, xcg);

      long lins = gkyl_range_idx(&edge_rng, idxs);
      long ling = gkyl_range_idx(&edge_rng, idxg);

      gcalc->equation->boundary_surf_term(gcalc->equation, dir, xcs, xcg,
        gcalc->grid.dx, gcalc->grid.dx, idxs, idxg, 1, 
	gkyl_array_cfetch(fIn, lins), gkyl_array_cfetch(fIn, ling), gkyl_array_fetch(rhs, ling)
      );
    }

    // Reset clower and cupper for the next iteration.
    clower_idx[dir] = phase_rng->lower[dir];
    cupper_idx[dir] = phase_rng->upper[dir];
  }
}

gkyl_ghost_surf_calc*
gkyl_ghost_surf_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation, int cdim, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_ghost_surf_calc_cu_dev_new(grid, equation, cdim);
  } 
#endif
  gkyl_ghost_surf_calc *up = gkyl_malloc(sizeof(gkyl_ghost_surf_calc));
  up->grid = *grid;
  up->equation = gkyl_dg_eqn_acquire(equation);
  up->cdim = cdim;
  
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  
  up->on_dev = up;
  
  return up;
}

void
gkyl_ghost_surf_calc_release(gkyl_ghost_surf_calc* up)
{
  gkyl_dg_eqn_release(up->equation);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_ghost_surf_calc_advance_cu(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs)
{
  assert(false);
}

gkyl_ghost_surf_calc*
gkyl_ghost_surf_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation, int cdim)
{
  assert(false);
}

#endif
