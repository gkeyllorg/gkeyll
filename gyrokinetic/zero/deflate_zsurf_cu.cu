/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_deflate_zsurf_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_deflate_zsurf_advance_cu_kernel(const struct gkyl_deflate_zsurf *up, int zidx, 
  struct gkyl_range range, struct gkyl_range deflated_range, 
  const struct gkyl_array* field, struct gkyl_array* deflated_field, int ncomp) 
{
  int idx[GKYL_MAX_DIM];
  int do_idx[3];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < deflated_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since deflated_range is a subrange
    gkyl_sub_range_inv_idx(&deflated_range, linc1, idx);

    for(int i = 0; i < up->cdim-1; i++)
      do_idx[i] = idx[i];
    do_idx[up->cdim-1] = zidx;
    long linc = gkyl_range_idx(&range, do_idx);
    const double *fld = (const double*) gkyl_array_cfetch(field, linc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc_deflated = gkyl_range_idx(&deflated_range, idx); 
    double *fld_deflated = (double*) gkyl_array_fetch(deflated_field, linc_deflated);
    for(int c = 0; c<ncomp; c++)
      up->kernel(&fld[c*up->num_basis], &fld_deflated[c*up->num_deflated_basis]);
  }
}

// Host-side wrapper for deflating 2d (x,z) modal expansion to a 1d (x) modal expansion
void
gkyl_deflate_zsurf_advance_cu(const struct gkyl_deflate_zsurf *up, int zidx, 
  const struct gkyl_range *range, const struct gkyl_range *deflated_range, 
  const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp)
{
  int nblocks = deflated_range->nblocks;
  int nthreads = deflated_range->nthreads;
  gkyl_deflate_zsurf_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, zidx, 
    *range, *deflated_range, field->on_dev, deflated_field->on_dev, ncomp);
}

// CUDA kernel to set device pointers to em vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
deflate_zsurf_set_cu_dev_ptrs(struct gkyl_deflate_zsurf *up, enum gkyl_basis_type b_type,
  int edge, int poly_order)
{
  up->kernel = deflate_zsurf_choose_kernel(b_type, up->cdim, edge, poly_order); // edge = 0,1 = lo, up
}

struct gkyl_deflate_zsurf* 
gkyl_deflate_zsurf_cu_dev_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis, int edge) 
{
  struct gkyl_deflate_zsurf *up = (struct gkyl_deflate_zsurf*) gkyl_malloc(sizeof(*up));

  up->num_basis = cbasis->num_basis;
  up->num_deflated_basis = deflated_cbasis->num_basis;
  up->cdim = cbasis->ndim;
  up->use_gpu = true;

  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_deflate_zsurf *up_cu = (struct gkyl_deflate_zsurf*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_deflate_zsurf), GKYL_CU_MEMCPY_H2D);

  deflate_zsurf_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, edge, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
