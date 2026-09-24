/* -*- c++ -*- */

extern "C" {
#include <gkyl_translate_dim.h>
#include <gkyl_translate_dim_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_trans_dim_set_cu_ker_ptrs(struct gkyl_translate_dim_kernels *kernels, int cdim_do,
  struct gkyl_basis basis_do, int cdim_tar, struct gkyl_basis basis_tar, int dir, enum gkyl_edge_loc edge)
{
  enum gkyl_basis_type basis_type = basis_tar.b_type;
  int poly_order = basis_tar.poly_order;
  int dir_idx = cdim_tar-1;
  int edge_idx = 0;
  if (cdim_tar < cdim_do) {
    dir_idx = dir;
    switch (edge) {
      case GKYL_LOWER_EDGE:
        edge_idx = 0;
        break;
      case GKYL_NO_EDGE:
        edge_idx = 1;
        break;
      case GKYL_UPPER_EDGE:
        edge_idx = 2;
        break;
    }
  }

  // Choose kernel that translates DG coefficients.
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
      kernels->translate = trans_dim_kern_list_gkhyb[cdim_tar+cdim_do-3].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->translate = trans_dim_kern_list_ser[cdim_do-1].list[dir_idx*3+edge_idx].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }

  // Choose the function that populates the donor index.
  int vdim = basis_do.ndim - cdim_do;
  if (vdim > 0) {
    kernels->get_idx_do = translate_dim_get_idx_do_gk;
  }
  else {
    if (cdim_tar < cdim_do)
      kernels->get_idx_do = translate_dim_get_idx_do_conf_down;
    else
      kernels->get_idx_do = translate_dim_get_idx_do_conf_up;
  }
};

void
trans_dim_choose_kernel_cu(struct gkyl_translate_dim_kernels *kernels, int cdim_do,
  struct gkyl_basis basis_do, int cdim_tar, struct gkyl_basis basis_tar, int dir, enum gkyl_edge_loc edge)
{
  gkyl_trans_dim_set_cu_ker_ptrs<<<1,1>>>(kernels, cdim_do, basis_do, cdim_tar, basis_tar, dir, edge);
}

__global__ static void
gkyl_translate_dim_advance_cu_ker(int cdim_do, int cdim_tar, int vdim_do, int vdim_tar,
  int num_basis_do, int num_basis_tar, int dir, struct gkyl_translate_dim_kernels *kernels,
  const struct gkyl_range rng_do, const struct gkyl_range rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, int ncomp,
  struct gkyl_array *GKYL_RESTRICT ftar)
{
  int idx_do[GKYL_MAX_DIM] = {0};
  int idx_tar[GKYL_MAX_DIM] = {0};

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < rng_tar.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&rng_tar, tid, idx_tar);

    // Translate the target idx to the donor idx:
    kernels->get_idx_do(cdim_tar, vdim_tar, idx_tar, &rng_do, cdim_do, idx_do, dir);

    long linidx_do = gkyl_range_idx(&rng_do, idx_do);
    long linidx_tar = gkyl_range_idx(&rng_tar, idx_tar);

    const double *fdo_c = (const double *) gkyl_array_cfetch(fdo, linidx_do);
    double *ftar_c = (double *) gkyl_array_fetch(ftar, linidx_tar);

    for (int n=0; n<ncomp; n++) {
      kernels->translate(&fdo_c[n*num_basis_do], &ftar_c[n*num_basis_tar]);
    }

  }
}

void
gkyl_translate_dim_advance_cu(gkyl_translate_dim* up,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, int ncomp,
  struct gkyl_array *GKYL_RESTRICT ftar)
{
  int nblocks = rng_tar->nblocks, nthreads = rng_tar->nthreads;

  gkyl_translate_dim_advance_cu_ker<<<nblocks, nthreads>>>
    (up->cdim_do, up->cdim_tar, up->vdim_do, up->vdim_tar, up->num_basis_do, up->num_basis_tar,
     up->dir, up->kernels, *rng_do, *rng_tar, fdo->on_dev, ncomp, ftar->on_dev);
}
