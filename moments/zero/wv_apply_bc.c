#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_eqn.h>

struct gkyl_wv_apply_bc {
  struct gkyl_rect_grid grid;
  int dir; // direction to apply BC
  enum gkyl_edge_loc edge; // edge to apply BC
  int nghost[GKYL_MAX_DIM]; // number of ghost cells

  const struct gkyl_wv_eqn *eqn; // equation 
  const struct gkyl_wave_geom *geom; // geometry needed for BCs
  
  wv_bc_func_t bcfunc; // function pointer
  void *ctx; // context to pass to function

  struct gkyl_range ext_range, range; // ranges on grid
  struct gkyl_range skin, ghost; // skin and ghost ranges
};

gkyl_wv_apply_bc*
gkyl_wv_apply_bc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_wv_eqn *eqn, const struct gkyl_wave_geom *geom,
  int dir, enum gkyl_edge_loc edge, const int *nghost,
  wv_bc_func_t bcfunc, void *ctx)
{
  gkyl_wv_apply_bc *up = gkyl_malloc(sizeof(gkyl_wv_apply_bc));

  up->grid = *grid;
  up->dir = dir;
  up->edge = edge;
  gkyl_copy_int_arr(grid->ndim, nghost, up->nghost);

  up->eqn = gkyl_wv_eqn_acquire(eqn);
  up->geom = gkyl_wave_geom_acquire(geom);
  
  up->bcfunc = bcfunc;
  up->ctx = ctx;

  // compute range and extended range over grid ....
  gkyl_create_grid_ranges(grid, nghost, &up->ext_range, &up->range);
  // ... from these compute skin and ghost ranges
  gkyl_skin_ghost_ranges(&up->skin, &up->ghost, dir, edge, &up->ext_range, nghost);

  return up;
}

void
gkyl_wv_apply_bc_advance(const gkyl_wv_apply_bc *bc, double tm,
  const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  enum gkyl_edge_loc edge = bc->edge;
  int dir = bc->dir, ndim = bc->grid.ndim, ncomp = out->ncomp;
  int meqn = bc->eqn->num_equations;

  double skin_local[meqn], ghost_local[meqn];

  // return immediately if update region does not touch boundary
  if ( (edge == GKYL_LOWER_EDGE) && (update_rng->lower[dir] > bc->range.lower[dir]) )
    return;
  if ( (edge == GKYL_UPPER_EDGE) && (update_rng->upper[dir] < bc->range.upper[dir]) )
    return;

  // compute intersection for region to update
  struct gkyl_range up_range;
  gkyl_range_intersect(&up_range, update_rng, &bc->skin);

  int edge_idx = (edge == GKYL_LOWER_EDGE) ? bc->range.lower[dir] : bc->range.upper[dir];
  int fact = (edge == GKYL_LOWER_EDGE) ? -1 : 1;

  int eidx[GKYL_MAX_DIM]; // index into geometry
  int gidx[GKYL_MAX_DIM]; // index into ghost cell
  
  // create iterator to walk over skin cells
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up_range);
  
  while (gkyl_range_iter_next(&iter)) {

    long sloc = gkyl_range_idx(update_rng, iter.idx);

    gkyl_copy_int_arr(ndim, iter.idx, gidx);
    
    // compute linear index into appropriate ghost-cell:
    // this strange indexing ensures that the ghost cell index is
    // essentially "reflection" of the skin cell index; might not be
    // correct for all BCs but I am not sure how else to handle
    // multiple ghost-cell situations
    gidx[dir] = 2*edge_idx-gidx[dir]+fact;
    long gloc = gkyl_range_idx(update_rng, gidx);

    // compute index into geometry: we always use the edge on the
    // physical boundary as the geometry may not be defined in the
    // ghost region. This is likely not right, but there does not seem
    // any clean way to do this
    if (edge == GKYL_LOWER_EDGE) {
      gkyl_copy_int_arr(ndim, iter.idx, eidx);
      eidx[dir] = edge_idx;
    }
    else {
      gkyl_copy_int_arr(ndim, gidx, eidx);
      eidx[dir] = edge_idx+1;
    }

    const struct gkyl_wave_cell_geom *wg = gkyl_wave_geom_get(bc->geom, eidx);

    // rotate skin data to local coordinates
    gkyl_wv_eqn_rotate_to_local(bc->eqn, wg->tau1[dir], wg->tau2[dir], wg->norm[dir],
      gkyl_array_fetch(out, sloc), skin_local);
      
    // apply boundary condition in local coordinates
    bc->bcfunc(bc->eqn, tm, ncomp, skin_local, ghost_local, bc->ctx);

    // rotate back to global
    gkyl_wv_eqn_rotate_to_global(bc->eqn, wg->tau1[dir], wg->tau2[dir], wg->norm[dir],
      ghost_local, gkyl_array_fetch(out, gloc));
  }
}

void
gkyl_wv_apply_bc_to_buff(const gkyl_wv_apply_bc *bc, double tm,
  const struct gkyl_range *update_rng, const struct gkyl_array *inp, double *buffer)
{
  enum gkyl_edge_loc edge = bc->edge;
  int dir = bc->dir, ndim = bc->grid.ndim, ncomp = inp->ncomp;
  int meqn = bc->eqn->num_equations;

  double skin_local[meqn], ghost_local[meqn];

  // return immediately if update region does not touch boundary
  if ( (edge == GKYL_LOWER_EDGE) && (update_rng->lower[dir] > bc->range.lower[dir]) )
    return;
  if ( (edge == GKYL_UPPER_EDGE) && (update_rng->upper[dir] < bc->range.upper[dir]) )
    return;

  // compute intersection for region to update
  struct gkyl_range up_range;
  gkyl_range_intersect(&up_range, update_rng, &bc->skin);

  int edge_idx = (edge == GKYL_LOWER_EDGE) ? bc->range.lower[dir] : bc->range.upper[dir];
  int fact = (edge == GKYL_LOWER_EDGE) ? -1 : 1;

  int ncells = bc->range.upper[dir]-bc->range.lower[dir]+1; // cells in 'dir'

  int sidx[GKYL_MAX_DIM]; // index into skin-cell geometry
  int gidx[GKYL_MAX_DIM]; // index into ghost-cell geometry
  
  // create iterator to walk over skin cells
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up_range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {

    long sloc = gkyl_range_idx(update_rng, iter.idx);

    // compute index into geometry
    if (edge == GKYL_LOWER_EDGE) {
      gkyl_copy_int_arr(ndim, iter.idx, sidx);
      sidx[dir] = iter.idx[dir];
      
      gkyl_copy_int_arr(ndim, iter.idx, gidx);
      gidx[dir] = sidx[dir]+ncells;
    }
    else {
      gkyl_copy_int_arr(ndim, iter.idx, sidx);
      sidx[dir] = iter.idx[dir]+1;
      
      gkyl_copy_int_arr(ndim, iter.idx, gidx);
      gidx[dir] = sidx[dir]-ncells;
    }

    const struct gkyl_wave_cell_geom *wgs = gkyl_wave_geom_get(bc->geom, sidx);

    // rotate skin data to local coordinates of skin-cell edge
    gkyl_wv_eqn_rotate_to_local(bc->eqn, wgs->tau1[dir], wgs->tau2[dir], wgs->norm[dir],
      gkyl_array_cfetch(inp, sloc), skin_local);
      
    // apply boundary condition in local coordinates
    bc->bcfunc(bc->eqn, tm, ncomp, skin_local, ghost_local, bc->ctx);

    // rotate back to global coordinates as defined on ghost cell edge
    const struct gkyl_wave_cell_geom *wgg = gkyl_wave_geom_get(bc->geom, gidx);
    gkyl_wv_eqn_rotate_to_global(bc->eqn, wgg->tau1[dir], wgg->tau2[dir], wgg->norm[dir],
      ghost_local, buffer+meqn*count);

    count += 1;
  }  
}

void
gkyl_wv_apply_bc_release(gkyl_wv_apply_bc* bc)
{
  gkyl_wv_eqn_release(bc->eqn);
  gkyl_wave_geom_release(bc->geom);
  gkyl_free(bc);
}
