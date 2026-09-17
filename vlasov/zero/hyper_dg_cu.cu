/* -*- c++ -*- */

#include "gkyl_alloc_flags_priv.h"
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_hyper_dg_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_hyper_dg_set_update_vol_cu_kernel(gkyl_hyper_dg *up, int update_vol_term)
{
  up->update_vol_term = update_vol_term;
}

__global__ static void
gkyl_hyper_dg_advance_cu_kernel(gkyl_hyper_dg* up, struct gkyl_range update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  int ndim = up->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < update_range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idxc);
    gkyl_rect_grid_cell_center(&up->grid, idxc, xcc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idxc);

    if (up->update_vol_term) {
      double cflr = up->equation->vol_term(
        up->equation, xcc, up->grid.dx, idxc,
        (const double*) gkyl_array_cfetch(fIn, linc), (double*) gkyl_array_fetch(rhs, linc)
      );
      double *cflrate_d = (double*) gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cflr; // frequencies are additive
    }
    
    for (int d=0; d<up->num_up_dirs; ++d) {
      int dir = up->update_dirs[d];
      double cfls = 0.0;
      gkyl_copy_int_arr(ndim, idxc, idxl);
      gkyl_copy_int_arr(ndim, idxc, idxr);
      // TODO: fix for arbitrary subrange
      if ((up->zero_flux_flags[dir]      && idxc[dir] == update_range.lower[dir]) ||
          (up->zero_flux_flags[dir+ndim] && idxc[dir] == update_range.upper[dir])) {
        edge = (idxc[dir] == update_range.lower[dir]) ? -1 : 1;
        // use idxl to store interior edge index (first index away from skin cell)
        idxl[dir] = idxl[dir]-edge;

        gkyl_rect_grid_cell_center(&up->grid, idxl, xcl);
        long linl = gkyl_range_idx(&update_range, idxl);

        cfls = up->equation->boundary_surf_term(up->equation,
          dir, xcl, xcc, up->grid.dx, up->grid.dx,
          idxl, idxc, edge,
          (const double*) gkyl_array_cfetch(fIn, linl), (const double*) gkyl_array_cfetch(fIn, linc),
          (double*) gkyl_array_fetch(rhs, linc)
        );
      }
      else {
        idxl[dir] = idxl[dir]-1;
        idxr[dir] = idxr[dir]+1;
        gkyl_rect_grid_cell_center(&up->grid, idxl, xcl);
        gkyl_rect_grid_cell_center(&up->grid, idxr, xcr);
        long linl = gkyl_range_idx(&update_range, idxl); 
        long linr = gkyl_range_idx(&update_range, idxr);

        cfls = up->equation->surf_term(up->equation,
          dir, xcl, xcc, xcr, up->grid.dx, up->grid.dx, up->grid.dx,
          idxl, idxc, idxr,
          (const double*) gkyl_array_cfetch(fIn, linl), (const double*) gkyl_array_cfetch(fIn, linc),
          (const double*) gkyl_array_cfetch(fIn, linr), (double*) gkyl_array_fetch(rhs, linc)
        );
      }
      double *cflrate_d = (double*) gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cfls; // frequencies are additive     
    }
  }
}

// wrapper to call advance kernel on device
void
gkyl_hyper_dg_advance_cu(gkyl_hyper_dg* up, const struct gkyl_range *update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  int nblocks = update_range->nblocks;
  int nthreads = update_range->nthreads;

  gkyl_hyper_dg_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *update_range,
    fIn->on_dev, cflrate->on_dev, rhs->on_dev);
}

__global__ static void
gkyl_hyper_dg_gen_stencil_advance_cu_kernel(gkyl_hyper_dg* hdg, long offsets[36], struct gkyl_range update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  // idxc, xc, and dx for volume update
  int idxc[GKYL_MAX_DIM] = {0};
  double xcc[GKYL_MAX_DIM] = {0.0};

  // idx for generic surface update
  int idx[9][GKYL_MAX_DIM] = {0};
  const double* fIn_d[9];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < update_range.volume; tid += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, tid, idxc);
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idxc);

    // Call volume kernel and get CFL rate
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);
    double cflr = hdg->equation->vol_term(
      hdg->equation, xcc, hdg->grid.dx, idxc,
      (const double*)gkyl_array_cfetch(fIn, linc), (double*) gkyl_array_fetch(rhs, linc));
    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, linc);
    cflrate_d[0] += cflr;

    for (int d1=0; d1<hdg->num_up_dirs; ++d1) {
      for (int d2=0; d2<hdg->num_up_dirs; ++d2) {
        int dir1 = hdg->update_dirs[d1];
        int dir2 = hdg->update_dirs[d2];

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
        int keri = idx_to_inloup_ker(num_up_dirs, idxc, update_dirs, update_range.upper);

        // Get pointers to all neighbor values
        for (int i=0; i<9; ++i) {
          gkyl_range_inv_idx(&update_range, linc+offsets_arr[i], idx[i]);
          fIn_d[i] = (const double*)gkyl_array_cfetch(fIn, linc + offsets_arr[i]);
        }

        // Domain stencil location is handled by the kernel selectors
        // gen_surf_term contains both surf and boundary surf kernels
        double cfls = hdg->equation->gen_surf_term(hdg->equation,
          dir1, dir2, xcc, hdg->grid.dx, idxc,
          keri, idx, fIn_d,
          (double *)gkyl_array_fetch(rhs, linc));

        cflrate_d[0] += cfls;
      }
    }
  }
}

// wrapper to call advance kernel on device
void
gkyl_hyper_dg_gen_stencil_advance_cu(gkyl_hyper_dg* hdg, long offsets[36], const struct gkyl_range *update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  int nblocks = update_range->nblocks;
  int nthreads = update_range->nthreads;

  long *offsets_cu = (long *)gkyl_cu_malloc(36*sizeof(long));
  gkyl_cu_memcpy(offsets_cu, offsets, 36*sizeof(long), GKYL_CU_MEMCPY_H2D);

  gkyl_hyper_dg_gen_stencil_advance_cu_kernel<<<nblocks, nthreads>>>(hdg->on_dev,
    offsets_cu, *update_range, fIn->on_dev, cflrate->on_dev, rhs->on_dev);

  gkyl_cu_free(offsets_cu);
}

void
gkyl_hyper_dg_set_update_vol_cu(gkyl_hyper_dg *up, int update_vol_term)
{
  gkyl_hyper_dg_set_update_vol_cu_kernel<<<1,1>>>(up, update_vol_term);
}

gkyl_hyper_dg*
gkyl_hyper_dg_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[2*GKYL_MAX_DIM],
  int update_vol_term)
{
  gkyl_hyper_dg *up = (gkyl_hyper_dg*) gkyl_malloc(sizeof(gkyl_hyper_dg));

  up->ndim = basis->ndim;
  up->num_basis = basis->num_basis;
  up->num_up_dirs = num_up_dirs;
  up->grid = *grid;

  for (int i=0; i<num_up_dirs; ++i)
    up->update_dirs[i] = update_dirs[i];

  for (int i=0; i<2*GKYL_MAX_DIM; ++i)
    up->zero_flux_flags[i] = zero_flux_flags[i];
    
  up->update_vol_term = update_vol_term;

  // aquire pointer to equation object
  struct gkyl_dg_eqn *eqn = gkyl_dg_eqn_acquire(equation);
  up->equation = eqn->on_dev; // this is so the memcpy below has eqn on_dev

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  
  // copy host struct to device struct
  gkyl_hyper_dg *up_cu = (gkyl_hyper_dg*) gkyl_cu_malloc(sizeof(gkyl_hyper_dg));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_hyper_dg), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu; // set parent pointer

  up->equation = eqn; // updater should store host pointer

  up->use_gpu = true;

  return up;
}
