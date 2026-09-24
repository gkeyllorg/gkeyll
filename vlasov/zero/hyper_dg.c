#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_hyper_dg_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static void
create_offsets(struct gkyl_hyper_dg *up, const struct gkyl_range *range, long offsets[])
{
  // Construct the offsets *only* in the directions being updated.
  // No need to load the neighbors that are not needed for the update.
  int lower_offset[GKYL_MAX_DIM] = {0};
  int upper_offset[GKYL_MAX_DIM] = {0};
  for (int d=0; d<up->num_up_dirs; ++d) {
    int dir = up->update_dirs[d];
    lower_offset[dir] = -1;
    upper_offset[dir] = 1;
  }  

  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, lower_offset, upper_offset);
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

void
gkyl_hyper_dg_set_update_vol(gkyl_hyper_dg *up, int update_vol_term)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_hyper_dg_set_update_vol_cu(up, update_vol_term);
    return;
  }
#endif
  up->update_vol_term = update_vol_term;
}

void
gkyl_hyper_dg_advance(struct gkyl_hyper_dg *up, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_hyper_dg_advance_cu(up, update_range, fIn, cflrate, rhs);
    return;
  }
#endif

  int ndim = up->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM], xc_edge[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(&up->grid, idxc, xcc);

    long linc = gkyl_range_idx(update_range, idxc);
    if (up->update_vol_term) {
      double cflr = up->equation->vol_term(
        up->equation, xcc, up->grid.dx, idxc,
        gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
      );
      double *cflrate_d = gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cflr; // frequencies are additive
    }
    
    for (int d=0; d<up->num_up_dirs; ++d) {
      int dir = up->update_dirs[d];
      double cfls = 0.0;
      // Assumes update_range owns lower and upper edges of the domain
      if ((up->zero_flux_flags[dir]      && idxc[dir] == update_range->lower[dir]) ||
          (up->zero_flux_flags[dir+ndim] && idxc[dir] == update_range->upper[dir]) ) {
        gkyl_copy_int_arr(ndim, iter.idx, idx_edge);
        edge = (idxc[dir] == update_range->lower[dir]) ? -1 : 1;
        // idx_edge stores interior edge index (first index away from skin cell)
        idx_edge[dir] = idx_edge[dir]-edge;

        gkyl_rect_grid_cell_center(&up->grid, idx_edge, xc_edge);
        long lin_edge = gkyl_range_idx(update_range, idx_edge);

        cfls = up->equation->boundary_surf_term(up->equation,
          dir, xc_edge, xcc, up->grid.dx, up->grid.dx,
          idx_edge, idxc, edge,
          gkyl_array_cfetch(fIn, lin_edge), gkyl_array_cfetch(fIn, linc),
          gkyl_array_fetch(rhs, linc)
        );
      }
      else {
        gkyl_copy_int_arr(ndim, iter.idx, idxl);
        gkyl_copy_int_arr(ndim, iter.idx, idxr);
        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;
        
        gkyl_rect_grid_cell_center(&up->grid, idxl, xcl);
        gkyl_rect_grid_cell_center(&up->grid, idxr, xcr);
        long linl = gkyl_range_idx(update_range, idxl); 
        long linr = gkyl_range_idx(update_range, idxr);

        cfls = up->equation->surf_term(up->equation,
          dir, xcl, xcc, xcr, up->grid.dx, up->grid.dx, up->grid.dx,
          idxl, idxc, idxr,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc), gkyl_array_cfetch(fIn, linr),
          gkyl_array_fetch(rhs, linc)
        );
      }
      double *cflrate_d = gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cfls; // frequencies are additive      
    }
  }
}

void
gkyl_hyper_dg_gen_stencil_advance(gkyl_hyper_dg *up, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = up->ndim;
  long sz[] = { 3, 9, 27 };
  long sz_dim = sz[up->num_up_dirs-1];
  long offsets[sz_dim];
  create_offsets(up, update_range, offsets);

  // idx, xc, and dx for volume update
  int idxc[GKYL_MAX_DIM];
  double xcc[GKYL_MAX_DIM];

  // idx, xc, and dx for generic surface update
  int idx[sz_dim][GKYL_MAX_DIM];
  double xc[sz_dim][GKYL_MAX_DIM];
  double dx[sz_dim][GKYL_MAX_DIM];
  const double* fIn_d[sz_dim];

  // bool for checking if index is in the domain
  int in_grid = 1;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(update_range, iter.idx);

    // Call volume kernel and get CFL rate
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(&up->grid, idxc, xcc);
    double cflr = up->equation->vol_term(
      up->equation, xcc, up->grid.dx, idxc,
      gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
    );
    double *cflrate_d = gkyl_array_fetch(cflrate, linc);
    cflrate_d[0] += cflr; // frequencies are additive

    // Get pointers to all neighbor values (i.e., 9 cells in 2D, 27 cells in 3D)
    for (int i=0; i<sz_dim; ++i) {
      // Get index based on offset (not necessarily a valid index) 
      gkyl_sub_range_inv_idx(update_range, linc+offsets[i], idx[i]);
      
      // Check if the index is in the domain
      // Assumes update_range owns lower and upper edges of the domain
      for (int d=0; d<up->num_up_dirs; ++d) {
        int dir = up->update_dirs[d];
        if (idx[i][dir] < update_range->lower[dir] || idx[i][dir] > update_range->upper[dir]) {
          in_grid = 0;
        }
      }
          
      // Only if the index is in the domain, fetch the pointer (otherwise pointer stays NULL)
      if (in_grid) {
        gkyl_rect_grid_cell_center(&up->grid, idx[i], xc[i]);
        for (int j=0; j<ndim; ++j)
          dx[i][j] = up->grid.dx[j];
        fIn_d[i] = gkyl_array_cfetch(fIn, linc + offsets[i]);
      }
      // reset in_grid for next neighbor value check
      in_grid = 1;
    }

    // Loop over surfaces and update using any/all neighbors needed
    // NOTE: ASSUMES UNIFORM GRIDS FOR NOW
    for (int d1=0; d1<up->num_up_dirs; ++d1) {
      for (int d2=0; d2<up->num_up_dirs; ++d2) {
        double cfls = 0.0;
        int dir1 = up->update_dirs[d1];
        int dir2 = up->update_dirs[d2];
        // Assumes update_range owns lower and upper edges of the domain
        if (idxc[dir1] == update_range->lower[dir1] || idxc[dir1] == update_range->upper[dir1]
             || idxc[dir2] == update_range->lower[dir2] || idxc[dir2] == update_range->upper[dir2]) {
          cfls = up->equation->gen_boundary_surf_term(up->equation,
            dir1, dir2, xcc, up->grid.dx, idxc,
            sz_dim, idx, fIn_d,
            gkyl_array_fetch(rhs, linc)
          );
        }
        else {
          cfls = up->equation->gen_surf_term(up->equation,
            dir1, dir2, xcc, up->grid.dx, idxc,
            sz_dim, idx, fIn_d,
            gkyl_array_fetch(rhs, linc)
          );
        }
        double *cflrate_d = gkyl_array_fetch(cflrate, linc);
        cflrate_d[0] += cfls; // frequencies are additive
      }
    }
  }
}

gkyl_hyper_dg*
gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[2*GKYL_MAX_DIM],
  int update_vol_term, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_hyper_dg_cu_dev_new(grid, basis, equation, num_up_dirs, update_dirs, zero_flux_flags, update_vol_term);
  } 
#endif
  gkyl_hyper_dg *up = gkyl_malloc(sizeof(gkyl_hyper_dg));

  up->grid = *grid;
  up->ndim = basis->ndim;
  up->num_basis = basis->num_basis;
  up->num_up_dirs = num_up_dirs;

  for (int i=0; i<num_up_dirs; ++i)
    up->update_dirs[i] = update_dirs[i];

  for (int i=0; i<2*GKYL_MAX_DIM; ++i)
    up->zero_flux_flags[i] = zero_flux_flags[i];
    
  up->update_vol_term = update_vol_term;
  up->equation = gkyl_dg_eqn_acquire(equation);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  
  up->on_dev = up; // on host, on_dev points to itself
  up->use_gpu = use_gpu;

  return up;
}

void gkyl_hyper_dg_release(struct gkyl_hyper_dg* up)
{
  gkyl_dg_eqn_release(up->equation);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

