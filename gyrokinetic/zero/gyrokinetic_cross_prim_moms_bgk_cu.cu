/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk_priv.h>
}

__global__ void
gkyl_gyrokinetic_cross_prim_moms_bgk_advance_cu_kernel(gkyl_gyrokinetic_cross_prim_moms_bgk *up,
  struct gkyl_range conf_range, double delta_sr, double betap1,
  double m_self, const struct gkyl_array *prim_moms_self, double m_other, const struct gkyl_array *prim_moms_other,
  struct gkyl_array *prim_moms_cross)
{ 
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);

    const double *prim_moms_self_d = (const double*) gkyl_array_cfetch(prim_moms_self, loc_conf);
    const double *prim_moms_other_d = (const double*) gkyl_array_cfetch(prim_moms_other, loc_conf);
    double *out_d = (double*) gkyl_array_fetch(prim_moms_cross, loc_conf);

    up->cross_prim_moms_calc(delta_sr, betap1, m_self, prim_moms_self_d, m_other, prim_moms_other_d, out_d);
  }
}

// Host-side wrapper for cross BGK moments
void
gkyl_gyrokinetic_cross_prim_moms_bgk_advance_cu(gkyl_gyrokinetic_cross_prim_moms_bgk *up,
  const struct gkyl_range *conf_range, double delta_sr, double betap1,
  double m_self, const struct gkyl_array *prim_moms_self, double m_other, const struct gkyl_array *prim_moms_other,
  struct gkyl_array *prim_moms_cross)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, delta_sr, betap1, m_self, prim_moms_self->on_dev, m_other, prim_moms_other->on_dev, 
    prim_moms_cross->on_dev);
}

__global__
static void
set_gyrokinetic_cross_prim_moms_bgk_cu_ptrs(struct gkyl_gyrokinetic_cross_prim_moms_bgk *up, int cdim, int vdim, int poly_order)
{
  up->cross_prim_moms_calc = choose_gyrokinetic_cross_prim_moms_bgk_kern(cdim, vdim, poly_order); 
}

gkyl_gyrokinetic_cross_prim_moms_bgk* 
gkyl_gyrokinetic_cross_prim_moms_bgk_cu_dev_new(const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis)
{
  struct gkyl_gyrokinetic_cross_prim_moms_bgk *up = (struct gkyl_gyrokinetic_cross_prim_moms_bgk*) gkyl_malloc(sizeof(*up));
  up->use_gpu = true;
  
  // copy struct to device
  struct gkyl_gyrokinetic_cross_prim_moms_bgk *up_cu = (struct gkyl_gyrokinetic_cross_prim_moms_bgk*) gkyl_cu_malloc(sizeof(*up_cu));
  
  set_gyrokinetic_cross_prim_moms_bgk_cu_ptrs<<<1,1>>>(up_cu, conf_basis->ndim,
    phase_basis->ndim-conf_basis->ndim, phase_basis->poly_order);

  up->on_dev = up_cu;  
  return up;   
}
