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
gkyl_hyper_dg_gen_stencil_advance(gkyl_hyper_dg *up, long offsets[36], const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_hyper_dg_gen_stencil_advance_cu(up, offsets, update_range, fIn, cflrate, rhs);
    return;
  }
#endif  
  int ndim = up->ndim;
  // idxc, xc, and dx for volume update
  int idxc[GKYL_MAX_DIM] = {0};
  double xcc[GKYL_MAX_DIM] = {0.0};

  // idx for generic surface update
  int idx[9][GKYL_MAX_DIM] = {0};
  const double* fIn_d[9] = {0};

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
    cflrate_d[0] += cflr;

    for (int d1=0; d1<up->num_up_dirs; ++d1) {
      for (int d2=0; d2<up->num_up_dirs; ++d2) {
        int dir1 = up->update_dirs[d1];
        int dir2 = up->update_dirs[d2];

        // Sort update directions from least to greatest for domain indexing 
        int update_dirs[2] = {0};
        update_dirs[0] = dir1 < dir2 ? dir1 : dir2;
        update_dirs[1] = dir1 < dir2 ? dir2 : dir1;

        // If dir1=dir2, only need in/lo/up domain stencil,
        // otherwise, need a 9-region domain stencil for cross terms
        int num_up_dirs = (d1 == d2) ? 1 : 2;
        int num_stencil = (d1 == d2) ? 3 : 9;
        int offsets_idx = (d1 == d2) ? d1*3 : (d1+d2)*9;

        // Index into relative offsets and kernel list
        long offsets_arr[9] = {0};
        gkyl_copy_long_arr(num_stencil, &offsets[offsets_idx], offsets_arr);
        int keri = idx_to_inloup_ker(num_up_dirs, idxc, update_dirs, update_range->upper);

        // Get pointers to all neighbor values
        for (int i=0; i<9; ++i) {
          gkyl_range_inv_idx(update_range, linc+offsets_arr[i], idx[i]);
          fIn_d[i] = gkyl_array_cfetch(fIn, linc + offsets_arr[i]);
        }

        // Domain stencil location is handled by the kernel selectors
        // gen_surf_term contains both surf and boundary surf kernels
        double cfls = up->equation->gen_surf_term(up->equation,
          dir1, dir2, xcc, up->grid.dx, idxc,
          keri, idx, fIn_d,
          gkyl_array_fetch(rhs, linc));

        cflrate_d[0] += cfls;
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
