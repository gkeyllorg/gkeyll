/* -*- c++ -*- */

extern "C" {
#include <gkyl_proj_exp_on_basis.h>
#include <gkyl_proj_exp_on_basis_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_proj_exp_on_basis_advance_cu_ker(int num_quad, 
  const struct gkyl_range range, const struct gkyl_array* GKYL_RESTRICT basis_at_ords,
  const struct gkyl_array* GKYL_RESTRICT weights, double alpha, double beta,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT fOut)
{
  int num_basis = basis_at_ords->ncomp;
  int tot_quad = basis_at_ords->size;

  int idx[GKYL_MAX_DIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, tid, idx);
    long linidx = gkyl_range_idx(&range, idx);

    const double *fIn_d = (const double *) gkyl_array_cfetch(fIn, linidx);

    double *fOut_d = (double *) gkyl_array_fetch(fOut, linidx);
    for (int k=0; k<num_basis; ++k) fOut_d[k] = 0.0;

    // Compute expansion coefficients of fOut using quadrature.
    const double *w_d = (const double *) weights->data;
    const double *bo_d = (const double *) basis_at_ords->data;

    for (int n=0; n<tot_quad; ++n) {

      const double *b_ord = (const double *) gkyl_array_cfetch(basis_at_ords, n);

      // Evaluate input function f at quad point.
      double fIn_q=0.;
      for (int k=0; k<num_basis; ++k)
        fIn_q += fIn_d[k]*b_ord[k];

      // Evaluate exp at quad point.
      double fOut_o = alpha*exp(beta*fIn_q);

      // Compute expansion coefficients.
      double tmp = w_d[n]*fOut_o;
      for (int k=0; k<num_basis; ++k)
        fOut_d[k] += tmp*bo_d[k+num_basis*n];
    }

  }
}

void
gkyl_proj_exp_on_basis_advance_cu(const gkyl_proj_exp_on_basis *up,
  const struct gkyl_range *range, double alpha, double beta,
  const struct gkyl_array *fIn, struct gkyl_array *fOut)
{
  int nblocks = range->nblocks, nthreads = range->nthreads;
  gkyl_proj_exp_on_basis_advance_cu_ker<<<nblocks, nthreads>>>
    (up->num_quad, *range, up->basis_at_ords->on_dev,
     up->weights->on_dev, alpha, beta, fIn->on_dev, fOut->on_dev);
}
