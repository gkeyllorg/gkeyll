/* -*- c++ -*- */

extern "C" {
#include <gkyl_dg_eval_at_coord_proj.h>
#include <gkyl_dg_eval_at_coord_proj_priv.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}

__global__ void
dg_eval_at_coord_choose_ker_cu_ker(int cdim, int ndim, struct gkyl_basis basis, int num_eval_dirs,
  dg_evproj_struct_int_t eval_dirs, struct dg_ev_proj_kernels *kers)
{
  int dir_mask = eval_dirs_to_mask(num_eval_dirs, eval_dirs.c);

  int poly_order = basis.poly_order;

  assert(ndim >= 1 && ndim <= 3);
  assert(dir_mask >= 1 && dir_mask < (1 << ndim));
  assert(poly_order >= 1 && poly_order <= 3);

  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->ev_ker = ser_eval_at_coord_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      kers->basis_ker = ser_basis_tar_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      kers->ev_ker = ten_eval_at_coord_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      kers->basis_ker = ten_basis_tar_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_GKHYBRID:
      kers->ev_ker = gkhyb_eval_at_coord_list[ndim-2][dir_mask-1].kernels[poly_order-1];
      kers->basis_ker = gkhyb_basis_tar_list[ndim-2][dir_mask-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
  assert(kers->ev_ker);
}

struct dg_ev_proj_kernels*
dg_eval_at_coord_choose_ker_cu(int cdim, int ndim, const struct gkyl_basis *basis,
  int num_eval_dirs, const int *eval_dirs)
{
  struct dg_ev_proj_kernels *kers = (struct dg_ev_proj_kernels *) gkyl_cu_malloc(sizeof(struct dg_ev_proj_kernels));

  dg_evproj_struct_int_t eval_dirs_st = {0};
  for (int i=0; i<num_eval_dirs; i++)
    eval_dirs_st.c[i] = eval_dirs[i];

  dg_eval_at_coord_choose_ker_cu_ker<<<1,1>>>(cdim, ndim, *basis, num_eval_dirs, eval_dirs_st, kers);

  return kers;
}

__global__ void
dg_eval_at_coord_proj_range_cu_kernel(int num_basis_do, int num_basis_tar, int ncomp,
  dg_evproj_struct_bool_t is_eval, dg_evproj_struct_double_t eval_coords_log, dg_evproj_struct_int_t cell_idx,
  struct dg_ev_proj_kernels *kers,
  struct gkyl_range rng_do, struct gkyl_range rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar)
{
  int idx_tar[GKYL_MAX_DIM];
  int idx_do[GKYL_MAX_DIM];
  const int ndim_do = rng_do.ndim;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < rng_tar.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&rng_tar, linc1, idx_tar);

    eval_at_coord_get_idx_do(is_eval.c, ndim_do, idx_tar, cell_idx.c, idx_do);

    long start_do  = gkyl_range_idx(&rng_do,  idx_do);
    long start_tar = gkyl_range_idx(&rng_tar, idx_tar);

    const double *fdo_c  = (const double *) gkyl_array_cfetch(fdo,  start_do);
    double       *ftar_c = (double *)       gkyl_array_fetch(ftar, start_tar);

    for (int n=0; n<ncomp; n++)
      kers->ev_ker(eval_coords_log.c, fdo_c+n*num_basis_do, ftar_c+n*num_basis_tar);
  }
}

void
gkyl_dg_eval_at_coord_proj_advance_cu(struct gkyl_dg_eval_at_coord_proj *up, const double *eval_coords,
  const struct gkyl_rect_grid *grid, const bool *pick_lower, const int *known_index,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar)
{

  // We assume that if fdo has multiple DG fields (vector components), ftar has the
  // same number of vector components.
  int ncomp = fdo->ncomp / up->num_basis_do;
  int num_basis_tar = fdo->ncomp / ncomp;

  // Build full ndim_do-dimensional point for gkyl_rect_grid_find_cell.
  double point[GKYL_MAX_DIM];
  int eval_ctr = 0;
  for (int d=0; d<up->ndim_do; d++) {
    if (up->is_eval[d])
      point[d] = eval_coords[eval_ctr++];
    else
      point[d] = grid->lower[d] + (rng_do->lower[d] - 0.5) * grid->dx[d];
  }

  dg_evproj_struct_int_t cell_idx = {0};
  gkyl_rect_grid_find_cell(grid, point, pick_lower, known_index, cell_idx.c);

  // Convert comp eval_coords to logical coords.
  dg_evproj_struct_double_t eval_coords_log = {0};
  for (int i=0; i<up->num_eval_dirs; i++) {
    int d = up->eval_dirs[i];
    double xc_d = grid->lower[d] + (cell_idx.c[d] - 0.5) * grid->dx[d];
    eval_coords_log.c[i] = 2.0 * (eval_coords[i] - xc_d) / grid->dx[d];
  }

  dg_evproj_struct_bool_t is_eval = {0};
  for (int d=0; d<GKYL_MAX_DIM; d++)
    is_eval.c[d] = up->is_eval[d];

  int nblocks = rng_tar->nblocks;
  int nthreads = rng_tar->nthreads;

  dg_eval_at_coord_proj_range_cu_kernel<<<nblocks, nthreads>>>(up->num_basis_do, num_basis_tar,
   ncomp, is_eval, eval_coords_log, cell_idx, up->kers,
   *rng_do, *rng_tar, fdo->on_dev, ftar->on_dev);
}
