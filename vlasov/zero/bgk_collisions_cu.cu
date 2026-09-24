/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_util.h>
#include <gkyl_bgk_collisions_priv.h>
}

__global__ void
gkyl_bgk_collisions_advance_cu_kernel(unsigned cdim, unsigned vdim, unsigned poly_order,
  unsigned pnum_basis, enum gkyl_basis_type b_type, double cellav_fac,
  struct gkyl_range crange, struct gkyl_range prange,
  const struct gkyl_array* nu, const struct gkyl_array* nufM, const struct gkyl_array* fin,
  bool implicit_step, double dt, struct gkyl_array* out, struct gkyl_array* cflfreq)
{
  mul_op_t mul_op = choose_mul_conf_phase_kern(b_type, cdim, vdim, poly_order);

  int pidx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < prange.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&prange, linc1, pidx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long pstart = gkyl_range_idx(&prange, pidx);
    long cstart = gkyl_range_idx(&crange, pidx);

    const double *nufM_d = (const double*) gkyl_array_cfetch(nufM, pstart);
    const double *fin_d = (const double*) gkyl_array_cfetch(fin, pstart);
    double *out_d = (double*) gkyl_array_fetch(out, pstart);

    const double *nu_d = (const double*) gkyl_array_cfetch(nu, cstart);

    // Add contribution to CFL frequency.
    if (implicit_step) {

      // Add nu*f_M.
      array_acc1(pnum_basis, out_d, 1.0/(1.0 + nu_d[0]*cellav_fac*dt), nufM_d);

      // Calculate and add -nu*f.
      double incr[160]; // mul_op assigns, but need increment, so use a buffer.
      mul_op(nu_d, fin_d, incr);
      array_acc1(pnum_basis, out_d, -1.0/(1.0 + nu_d[0]*cellav_fac*dt), incr);

      // No CFL contribution in the implicit case
    } 
    else {
      // Add nu*f_M.
      array_acc1(pnum_basis, out_d, 1., nufM_d);

      // Calculate -nu*f.
      double incr[160]; // mul_op assigns, but need increment, so use a buffer.
      mul_op(nu_d, fin_d, incr);
      array_acc1(pnum_basis, out_d, -1., incr);

      // Add contribution to CFL frequency.
      double *cflfreq_d = (double *) gkyl_array_fetch(cflfreq, pstart);
      cflfreq_d[0] += nu_d[0]*cellav_fac;
    }
  }
}

void
gkyl_bgk_collisions_advance_cu(const gkyl_bgk_collisions *up,
  const struct gkyl_range *crange, const struct gkyl_range *prange,
  const struct gkyl_array *nu, const struct gkyl_array *nufM, const struct gkyl_array *fin,
  bool implicit_step, double dt, struct gkyl_array *out, struct gkyl_array *cflfreq)
{
  int nblocks = prange->nblocks;
  int nthreads = prange->nthreads;
  gkyl_bgk_collisions_advance_cu_kernel<<<nblocks, nthreads>>>(up->cdim, up->vdim,
    up->poly_order, up->pnum_basis, up->pb_type, up->cellav_fac, *crange, *prange,
    nu->on_dev, nufM->on_dev, fin->on_dev, implicit_step, dt, out->on_dev, cflfreq->on_dev);
}
