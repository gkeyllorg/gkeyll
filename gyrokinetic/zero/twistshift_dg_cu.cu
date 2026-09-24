#include <cstdio>
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
#include <gkyl_twistshift_dg.h>
#include <gkyl_twistshift_dg_priv.h>
}

#include <cstdio>
#include <cassert>

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

__global__ void
gkyl_twistshift_dg_set_distf_mats_cu_ker(const struct gkyl_array *fdo, const long *num_numcol_fidx_do, struct gkyl_nmat *fmat)
{
  // Assign the distribution matrices.
  // This assumes that fdo->ncomp = fmat->nr.
  // Recall:
  //   fmat->num = sum_i^Nx num_do(i)
  //   fmat->nr = num_basis = ncomp
  //   fmat->nc = Ny*Nvpar*Nmu

  long nc_idx = START_ID / fmat->nr; // num-num_col index: current num-num_col plane.
  int row_idx = START_ID % fmat->nr; // row index: current DG coeff.

  if ((nc_idx < fmat->num * fmat->nc) && (row_idx < fdo->ncomp)) {
    const double *fdo_c = (const double*) gkyl_array_cfetch(fdo, num_numcol_fidx_do[nc_idx]);
    struct gkyl_mat mcurr = gkyl_nmat_get(fmat, nc_idx % fmat->num);

    gkyl_mat_set(&mcurr, row_idx, nc_idx/fmat->num, fdo_c[row_idx]);
  }
}

__global__ void
gkyl_twistshift_dg_add_contr_cu_ker(struct gkyl_array *ftar, long *num_numcol_fidx_tar, int num_cells_skin,
  struct gkyl_nmat *mm_contr, int *num_do_cum, struct gkyl_range permutted_ghost_r, struct gkyl_rect_grid grid)
{
  long linidx_tar = START_ID / ftar->ncomp;
  int row_idx = START_ID % ftar->ncomp;

  // This if-statement may only be needed in GPU kernel, not for CPUs.
  if ((linidx_tar < num_cells_skin) && (row_idx < ftar->ncomp)) {
    double *ftar_c = (double*) gkyl_array_fetch(ftar, num_numcol_fidx_tar[linidx_tar]);

    int idx[GKYL_MAX_DIM] = {1};
    gkyl_sub_range_inv_idx(&permutted_ghost_r, linidx_tar, idx);

    int ac[GKYL_MAX_DIM] = {1};
    for (int d=2; d<grid.ndim-1; d++)
      ac[d-2] = grid.cells[d+1];
    ac[permutted_ghost_r.ndim-2] = mm_contr->num;

    int start = 0;
    for (int d=0; d<permutted_ghost_r.ndim-1; d++)
      start = (start + (idx[d]-1)) * ac[d];

    int shear_idx = idx[permutted_ghost_r.ndim-1];

    int do_start = num_do_cum[shear_idx-1];
    int do_end   = num_do_cum[shear_idx-1+1];
    for (int j=do_start; j<do_end; j++) { // Only loop over num_do[i] elements.
      int linidx_mm_contr = start + j;
      struct gkyl_mat mat = gkyl_nmat_get(mm_contr, linidx_mm_contr % mm_contr->num);
      ftar_c[row_idx] += gkyl_mat_get(&mat, row_idx, linidx_mm_contr / mm_contr->num);
    }
  }
}

void
gkyl_twistshift_dg_advance_cu(struct gkyl_twistshift_dg *up, struct gkyl_array *fdo, struct gkyl_array *ftar)
{
  // Set the columns of the donor matrix with the donor distributions.
  int num_blocks_set = (up->fmat->num * up->fmat->nr * up->fmat->nc+GKYL_DEFAULT_NUM_THREADS-1)/GKYL_DEFAULT_NUM_THREADS;
  gkyl_twistshift_dg_set_distf_mats_cu_ker<<<num_blocks_set, GKYL_DEFAULT_NUM_THREADS>>>
    (fdo->on_dev, up->num_numcol_fidx_do, up->fmat->on_dev);

  // Perform the mat-mat multiplications.
  gkyl_nmat_mm(1.0, 0.0, GKYL_NO_TRANS, up->scimat, GKYL_NO_TRANS, up->fmat, up->mm_contr);

  // Clear the ghost range.
  gkyl_array_clear_range(ftar, 0.0, &up->ghost_r);

  // Add the contributions of mat-vec multiplications.
  int num_cells_skin = (up->shear_r.upper[0]-up->shear_r.lower[0]+1) * up->fmat->nc;
  int num_blocks_add = (ftar->ncomp * num_cells_skin+GKYL_DEFAULT_NUM_THREADS-1)/GKYL_DEFAULT_NUM_THREADS;
  gkyl_twistshift_dg_add_contr_cu_ker<<<num_blocks_add,GKYL_DEFAULT_NUM_THREADS>>>
    (ftar->on_dev, up->num_numcol_fidx_tar, num_cells_skin, up->mm_contr->on_dev,
     up->num_do_cum, up->permutted_ghost_r, up->grid);
}
