#include <assert.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_boundary_flux.h>
#include <gkyl_boundary_flux_priv.h>
#include <gkyl_util.h>

gkyl_boundary_flux*
gkyl_boundary_flux_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  int num_eqns, const struct gkyl_dg_eqn **eqns, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_boundary_flux_cu_dev_new(dir, edge, grid, skin_r, ghost_r, num_eqns, eqns);
  }
#endif

  gkyl_boundary_flux *up = gkyl_malloc(sizeof(gkyl_boundary_flux));

  up->dir = dir;
  up->edge = edge;
  up->grid = *grid;
  up->skin_r = *skin_r;
  up->ghost_r = *ghost_r;
  up->use_gpu = use_gpu;

  up->num_eqns = num_eqns;
  up->eqns = gkyl_malloc(up->num_eqns*sizeof(struct gkyl_dg_eqn *));
  for (int i=0; i<up->num_eqns; i++)
    up->eqns[i] = gkyl_dg_eqn_acquire(eqns[i]);

  up->eqns_ho = up->eqns;
  
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  
  up->on_dev = up;
  
  return up;
}

void
gkyl_boundary_flux_advance(gkyl_boundary_flux *up,
  const struct gkyl_array *fIn, struct gkyl_array *fluxOut)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_boundary_flux_advance_cu(up, fIn, fluxOut);
    return;
  } 
#endif

  int idx_s[GKYL_MAX_DIM];
  double xc_g[GKYL_MAX_DIM], xc_s[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    int *idx_g = iter.idx;
    gkyl_copy_int_arr(up->ghost_r.ndim, idx_g, idx_s);
    idx_s[up->dir] = up->edge == GKYL_LOWER_EDGE? idx_g[up->dir]+1 : idx_g[up->dir]-1;

    gkyl_rect_grid_cell_center(&up->grid, idx_g, xc_g);
    gkyl_rect_grid_cell_center(&up->grid, idx_s, xc_s);

    long linidx_g = gkyl_range_idx(&up->ghost_r, idx_g); 
    long linidx_s = gkyl_range_idx(&up->skin_r, idx_s);

    const double *fIn_s = gkyl_array_cfetch(fIn, linidx_s);
    const double *fIn_g = gkyl_array_cfetch(fIn, linidx_g);
    double *fluxOut_g = gkyl_array_fetch(fluxOut, linidx_g);

    for (int i=0; i<up->num_eqns; i++) {
      up->eqns[i]->boundary_diag_term(up->eqns[i], up->dir, xc_s, xc_g,
        up->grid.dx, up->grid.dx, idx_s, idx_g, up->edge == GKYL_LOWER_EDGE? -1 : 1,
        fIn_s, fIn_g, fluxOut_g);
    }
  }
}

void
gkyl_boundary_flux_release(gkyl_boundary_flux* up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->eqns);
  } 
#endif
  for (int i=0; i<up->num_eqns; i++)
    gkyl_dg_eqn_release(up->eqns_ho[i]);
  gkyl_free(up->eqns_ho);

  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
