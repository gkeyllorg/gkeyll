/* -*- c++ -*- */
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_gk_neut_hamil.h>
#include <gkyl_dg_calc_gk_neut_hamil_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_dg_calc_gk_neut_hamil_set_cu_dev_ptrs(struct gkyl_dg_calc_gk_neut_hamil *up, 
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  up->calc_hamil = choose_kern(b_type, cdim, vdim, poly_order);
};

__global__ static void
gkyl_dg_calc_gk_neut_hamil_calc_cu_kernel(struct gkyl_dg_calc_gk_neut_hamil *up,
  const struct gkyl_range conf_range, const struct gkyl_range phase_range,
  const struct gkyl_array* gij, struct gkyl_array* hamil)
{
  int idx[GKYL_MAX_DIM];
  // Cell center array
  double xc[GKYL_MAX_DIM];  
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *gij_d = (const double*) gkyl_array_cfetch(gij, loc_conf);
    double *hamil_d = (double*) gkyl_array_fetch(hamil, loc_phase);

    up->calc_hamil(xc, up->phase_grid.dx, gij_d, hamil_d);
  }
}

void gkyl_dg_calc_gk_neut_hamil_calc_cu(struct gkyl_dg_calc_gk_neut_hamil *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* gij, struct gkyl_array* hamil)
{

  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_gk_neut_hamil_calc_cu_kernel<<<nblocks, nthreads>>>(up->on_dev,
    *conf_range, *phase_range, gij->on_dev, hamil->on_dev);
 
}

gkyl_dg_calc_gk_neut_hamil*
gkyl_dg_calc_gk_neut_hamil_cu_dev_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *basis, int cdim)
{    
  gkyl_dg_calc_gk_neut_hamil *up = (struct gkyl_dg_calc_gk_neut_hamil*) gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  int vdim = 3; 
  int poly_order = basis->poly_order;
  enum gkyl_basis_type b_type = basis->b_type;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_gk_neut_hamil *up_cu = (struct gkyl_dg_calc_gk_neut_hamil*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_gk_neut_hamil), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_calc_gk_neut_hamil_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, vdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
