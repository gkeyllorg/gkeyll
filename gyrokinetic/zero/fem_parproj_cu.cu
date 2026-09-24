/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>
}

// CUDA kernel to set device pointers to l2g, RHS src and solution
// kernels. Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_parproj_set_cu_ker_ptrs(struct gkyl_fem_parproj_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order, bool has_weight_lhs, bool has_weight_rhs,
  enum gkyl_fem_parproj_bc_type bctype)
{
  // Set l2g kernels.
  int bckey_periodic = bctype == GKYL_FEM_PARPROJ_PERIODIC? 0 : 1;
  const local2global_kern_list *local2global_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      local2global_kernels = ser_loc2glob_list;
      break;
    default:
      assert(false);
      break;
  }
  for (int k=0; k<2; k++)
    kers->l2g[k] = CK(local2global_kernels, dim, bckey_periodic, poly_order, k);

  // Set RHS stencil kernels.
  int bckey_dirichlet;
  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    bckey_dirichlet = 1;
  else if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN)
    bckey_dirichlet = 2;
  else
    bckey_dirichlet = 0;

  const srcstencil_kern_list *srcstencil_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      srcstencil_kernels = has_weight_rhs? ser_srcstencil_list_weighted : ser_srcstencil_list_noweight;
      break;
    default:
      assert(false);
  }
  for (int k=0; k<3; k++)
    kers->srcker[k] = CK(srcstencil_kernels, dim, bckey_dirichlet, poly_order, k);

  // Set the get solution stencil kernel.
  const solstencil_kern_list *solstencil_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      solstencil_kernels = ser_solstencil_list;
      break;
    default:
      assert(false);
  }
  kers->solker = solstencil_kernels[dim-1].kernels[poly_order-1];

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    kers->get_dirichlet_value = get_dirichlet_value_enabled_ghost;
  else if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN)
    kers->get_dirichlet_value = get_dirichlet_value_enabled_skin;
  else
    kers->get_dirichlet_value = get_dirichlet_value_disabled;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<2; k++)
        kers->bias_src_ker[k] = CK(ser_bias_src_list, dim, bckey_periodic, poly_order, k); 

      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

void
fem_parproj_choose_kernels_cu(const struct gkyl_basis *basis, bool has_weight_lhs, bool has_weight_rhs,
  enum gkyl_fem_parproj_bc_type bctype, struct gkyl_fem_parproj_kernels *kers)
{
  fem_parproj_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, basis->ndim,
    basis->poly_order, has_weight_lhs, has_weight_rhs, bctype);
}

__global__ void
gkyl_fem_parproj_set_rhs_kernel(double *rhs_global, const struct gkyl_array *rhsin,
  const struct gkyl_array *weight, const struct gkyl_array *phibc,
  struct gkyl_range range, struct gkyl_range perp_range2d, struct gkyl_range par_range1d,
  struct gkyl_fem_parproj_kernels *kers, long numnodes_global)
{
  int idx[GKYL_MAX_CDIM];
  long globalidx[32];
  int parnum_cells = range.upper[range.ndim-1]-range.lower[range.ndim-1]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linidx = gkyl_range_idx(&range, idx);

    const double *wgt_p = weight? (const double *) gkyl_array_cfetch(weight, linidx) : NULL;
    const double *rhsin_p = (const double*) gkyl_array_cfetch(rhsin, linidx);
    const double *phibc_p = kers->get_dirichlet_value(range.ndim-1, parnum_cells, idx, &range, phibc);

    int idx1d[] = {idx[range.ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    int keri = idx1d[0] == parnum_cells? 1 : 0;
    kers->l2g[keri](parnum_cells, paridx, globalidx);

    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<range.ndim-1; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(parnum_cells, idx1d[0]);
    kers->srcker[keri](wgt_p, rhsin_p, phibc_p, perpProbOff, globalidx, rhs_global);
  }
}

void
gkyl_fem_parproj_set_rhs_cu(gkyl_fem_parproj *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{
  gkyl_culinsolver_clear_rhs(up->prob_cu, 0);
  double *rhs_cu = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);

  const struct gkyl_array *phibc_cu = phibc? phibc->on_dev : NULL;
  const struct gkyl_array *wgt_cu = up->has_weight_rhs? up->weight_rhs->on_dev : NULL;

  gkyl_fem_parproj_set_rhs_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(rhs_cu, rhsin->on_dev, wgt_cu, phibc_cu,
    *up->solve_range, up->perp_range2d, up->par_range1d, up->kernels, up->numnodes_global);

  // Set the corresponding entries to the biasing potential.
  up->bias_line_src(up, rhsin);
}

__global__ void
gkyl_fem_parproj_bias_src_kernel(double *rhs_global, struct gkyl_rect_grid grid,
  struct gkyl_range range, struct gkyl_range perp_range2d, struct gkyl_range par_range1d,
  struct gkyl_fem_parproj_kernels *kers, long numnodes_global,
  int num_bias_line, struct gkyl_poisson_bias_line *bias_lines)
{
  const int bl_ndim_perp = 2;

  int ndim_perp = range.ndim-1;
  int parnum_cells = range.upper[ndim_perp]-range.lower[ndim_perp]+1;

  int idx[GKYL_MAX_CDIM];
  long globalidx[32];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume; linc1 += gridDim.x*blockDim.x)
  {
    // Inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    int idx1d[] = {idx[ndim_perp]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    int keri = idx1d[0] == parnum_cells? 1 : 0;
    kers->l2g[keri](parnum_cells, paridx, globalidx);

    // Modify the RHS source to enforce biasing of the solution.
    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<ndim_perp; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    for (int i=0; i<num_bias_line; i++) {
      // Index of the cell that abuts the line from below.
      struct gkyl_poisson_bias_line *bl = &bias_lines[i];
      int bl_idx_m[bl_ndim_perp];
      for (int d=0; d<bl_ndim_perp; d++) {
        int perp_dir = bl->perp_dirs[d];
        double dx = grid.dx[perp_dir];
        bl_idx_m[d] = (bl->perp_coords[d]-1e-3*dx - grid.lower[perp_dir])/dx+1;
      }

      if (
          ( idx[bl->perp_dirs[0]] == bl_idx_m[0]   && idx[bl->perp_dirs[1]] == bl_idx_m[1]   ) ||
          ( idx[bl->perp_dirs[0]] == bl_idx_m[0]+1 && idx[bl->perp_dirs[1]] == bl_idx_m[1]   ) ||
          ( idx[bl->perp_dirs[0]] == bl_idx_m[0]   && idx[bl->perp_dirs[1]] == bl_idx_m[1]+1 ) ||
          ( idx[bl->perp_dirs[0]] == bl_idx_m[0]+1 && idx[bl->perp_dirs[1]] == bl_idx_m[1]+1 )
         ) {
        int edge[2] = {
          -1+2*((bl_idx_m[0]+1)-idx[bl->perp_dirs[0]]),
          -1+2*((bl_idx_m[1]+1)-idx[bl->perp_dirs[1]]),
        };
        kers->bias_src_ker[keri](edge, bl->perp_dirs, bl->val, perpProbOff, globalidx, rhs_global);
      }
    }
  }
}

void
gkyl_fem_parproj_bias_src_enabled_cu(gkyl_fem_parproj *up, const struct gkyl_array *rhsin)
{
  double *rhs_cu = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);
  gkyl_fem_parproj_bias_src_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(rhs_cu, up->grid,
    *up->solve_range, up->perp_range2d, up->par_range1d, up->kernels, up->numnodes_global,
    up->num_bias_line, up->bias_lines);
}

__global__ void
gkyl_fem_parproj_get_sol_kernel(struct gkyl_array *phiout, const double *x_global, struct gkyl_range range,
  struct gkyl_range perp_range2d, struct gkyl_range par_range1d, struct gkyl_fem_parproj_kernels *kers,
  long numnodes_global)
{
  int idx[GKYL_MAX_DIM];
  long globalidx[32];
  int parnum_cells = range.upper[range.ndim-1]-range.lower[range.ndim-1]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linidx = gkyl_range_idx(&range, idx);
    double *phiout_p = (double*) gkyl_array_cfetch(phiout, linidx);

    int idx1d[] = {idx[range.ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    int keri = idx1d[0] == parnum_cells? 1 : 0;
    kers->l2g[keri](parnum_cells, paridx, globalidx);

    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<range.ndim-1; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    kers->solker(x_global, perpProbOff, globalidx, phiout_p);
  }
}

void
gkyl_fem_parproj_solve_cu(gkyl_fem_parproj *up, struct gkyl_array *phiout)
{
  gkyl_culinsolver_solve(up->prob_cu);
  double *x_cu = gkyl_culinsolver_get_sol_ptr(up->prob_cu, 0);

  gkyl_fem_parproj_get_sol_kernel<<<phiout->nblocks, phiout->nthreads>>>(phiout->on_dev, x_cu,
    *up->solve_range, up->perp_range2d, up->par_range1d, up->kernels, up->numnodes_global);
}
