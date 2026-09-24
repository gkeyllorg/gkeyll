/* -*- c++ -*- */

extern "C" {
#include <gkyl_positivity_shift_vlasov.h>
#include <gkyl_positivity_shift_vlasov_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_pos_shift_vlasov_set_cu_ker_ptrs(struct gkyl_positivity_shift_vlasov_kernels *kernels,
  struct gkyl_basis cbasis, struct gkyl_basis pbasis, enum gkyl_positivity_shift_type stype)
{
  int cdim = cbasis.ndim, pdim = pbasis.ndim;
  int vdim = pdim-cdim;
  enum gkyl_basis_type cbasis_type = cbasis.b_type, pbasis_type = pbasis.b_type;
  int poly_order = pbasis.poly_order;

  int plin = pos_shift_vlasov_cv_index[cdim].vdim[vdim];

  switch (pbasis_type) {
    case GKYL_BASIS_MODAL_TENSOR:
      kernels->is_m0_positive = pos_shift_vlasov_kern_list_m0_pos_check_tensor[cdim-1].kernels[poly_order-1];
      kernels->shift = stype == GKYL_POSITIVITY_SHIFT_TYPE_SHIFT_ONLY?
        pos_shift_vlasov_kern_list_shift_tensor[plin].kernels[poly_order-1] :
        pos_shift_vlasov_kern_list_MRSlimiter_tensor[plin].kernels[poly_order-1];
      kernels->m0 = pos_shift_vlasov_kern_list_m0_tensor[plin].kernels[poly_order-1];
      kernels->conf_phase_mul_op = choose_mul_conf_phase_kern(pbasis_type, cdim, vdim, poly_order);
      break;
    default:
      assert(false);
      break;
  }

  switch (cbasis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->conf_inv_op = choose_ser_inv_kern(cdim, poly_order);
      kernels->conf_mul_op = choose_ser_mul_kern(cdim, poly_order);
      break;
    default:
      assert(false);
      break;
  }
};

void
pos_shift_vlasov_choose_shift_kernel_cu(struct gkyl_positivity_shift_vlasov_kernels *kernels,
  struct gkyl_basis cbasis, struct gkyl_basis pbasis, enum gkyl_positivity_shift_type stype)
{
  gkyl_pos_shift_vlasov_set_cu_ker_ptrs<<<1,1>>>(kernels, cbasis, pbasis, stype);
}

// Function borrowed from array_reduce_cu.cu.
__device__ static __forceinline__ double
pos_shift_atomicMax_double(double *address, double val)
{
  unsigned long long int ret = __double_as_longlong(*address);
  while(val > __longlong_as_double(ret))
  {
    unsigned long long int old = ret;
    if((ret = atomicCAS((unsigned long long int*)address, old, __double_as_longlong(val))) == old)
      break;
  }
  return __longlong_as_double(ret);
}

__global__ void
gkyl_positivity_shift_vlasov_advance_int_array_clear_cu_ker(struct gkyl_array* out, int val)
{
  int *out_d = (int*) out->data;
  unsigned long start_id = threadIdx.x + blockIdx.x*blockDim.x;
  unsigned long nelm = out->size*out->ncomp;
  for (unsigned long linc = start_id; linc < nelm; linc += blockDim.x*gridDim.x)
    out_d[linc] = val;
}

__global__ static void
gkyl_positivity_shift_vlasov_advance_shift_cu_ker(
  struct gkyl_positivity_shift_vlasov_kernels *kers, const struct gkyl_rect_grid grid,
  const struct gkyl_range conf_range, const struct gkyl_range phase_range,
  double *ffloor, double ffloor_fac, double cellav_fac, struct gkyl_array* GKYL_RESTRICT shiftedf,
  struct gkyl_array* GKYL_RESTRICT distf, struct gkyl_array* GKYL_RESTRICT m0, struct gkyl_array* GKYL_RESTRICT delta_m0)
{
  int pidx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  double distf_max = -DBL_MAX;
  const int num_cbasis = 20; // MF 2024/09/03: Hardcoded to p=2 3x ser for now.

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    long clinidx = gkyl_range_idx(&conf_range, pidx);
    long plinidx = gkyl_range_idx(&phase_range, pidx);

    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    int *shiftedf_c = (int*) gkyl_array_fetch(shiftedf, clinidx);
    double *m0_c = (double*) gkyl_array_fetch(m0, clinidx);
    double *delta_m0_c = (double*) gkyl_array_fetch(delta_m0, clinidx);
    double *distf_c = (double*) gkyl_array_fetch(distf, plinidx);

    // Contribution to the old number density from this v-space cell.
    double m0Local_in[num_cbasis];
    for (unsigned int k=0; k<delta_m0->ncomp; ++k)
      m0Local_in[k] = 0.0;
    kers->m0(xc, grid.dx, pidx, distf_c, m0Local_in);

    // Add to the old number density.
    for (unsigned int k = 0; k < delta_m0->ncomp; ++k)
      atomicAdd(&delta_m0_c[k], m0Local_in[k]);

    // Shift f if needed.
    bool shifted_node = kers->shift(ffloor[0], distf_c);

    if (shifted_node) {
      // Compute the new number density local to this phase-space cell.
      double m0Local_out[num_cbasis];
      for (unsigned int k=0; k<m0->ncomp; ++k)
        m0Local_out[k] = 0.0;
      kers->m0(xc, grid.dx, pidx, distf_c, m0Local_out);

      if (kers->is_m0_positive(m0Local_in)) {
        // Rescale f in this cell so it keeps the same density.
        double m0ratio_c[num_cbasis];
        kers->conf_inv_op(m0Local_out, m0ratio_c);
        kers->conf_mul_op(m0Local_in, m0ratio_c, m0ratio_c);

        kers->conf_phase_mul_op(m0ratio_c, distf_c, distf_c);

        // Add contribution from this phase-space cell to the new number density.
        for (unsigned int k = 0; k < m0->ncomp; ++k)
          atomicAdd(&m0_c[k], m0Local_in[k]);
      }
      else {
        // Add contribution from this phase-space cell to the new number density.
        for (unsigned int k = 0; k < m0->ncomp; ++k)
          atomicAdd(&m0_c[k], m0Local_out[k]);

        atomicOr(shiftedf_c, shifted_node);
      }
    }
    else {
      // Add contribution from this phase-space cell to the new number density.
      for (unsigned int k = 0; k < m0->ncomp; ++k)
        atomicAdd(&m0_c[k], m0Local_in[k]);
    }

    distf_max = fmax(distf_max, distf_c[0]);

  }

  pos_shift_atomicMax_double(ffloor, ffloor_fac * distf_max * cellav_fac);
}

__global__ static void
gkyl_positivity_shift_vlasov_advance_scalef_cu_ker(
  struct gkyl_positivity_shift_vlasov_kernels *kers,
  const struct gkyl_range conf_range, const struct gkyl_range phase_range,
  const struct gkyl_array* GKYL_RESTRICT shiftedf, const struct gkyl_array* GKYL_RESTRICT m0,
  const struct gkyl_array* GKYL_RESTRICT delta_m0, struct gkyl_array* GKYL_RESTRICT distf)
{
  int pidx[GKYL_MAX_DIM];

  const int num_cbasis = 20; // MF 2024/09/03: Hardcoded to p=2 3x ser for now.

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    long clinidx = gkyl_range_idx(&conf_range, pidx);
    const int *shiftedf_c = (const int*) gkyl_array_cfetch(shiftedf, clinidx);

    if (shiftedf_c[0]) {
      const double *delta_m0_c = (const double*) gkyl_array_cfetch(delta_m0, clinidx);
      if (kers->is_m0_positive(delta_m0_c)) {
        // Rescale f so it has the same m0 at this conf-space cell.
        const double *m0_c = (const double*) gkyl_array_cfetch(m0, clinidx);
        double m0ratio_c[num_cbasis];
        kers->conf_inv_op(m0_c, m0ratio_c);
        kers->conf_mul_op(delta_m0_c, m0ratio_c, m0ratio_c);

        long plinidx = gkyl_range_idx(&phase_range, pidx);
        double *distf_c = (double*) gkyl_array_fetch(distf, plinidx);
        kers->conf_phase_mul_op(m0ratio_c, distf_c, distf_c);
      }
    }
  }
}

__global__ static void
gkyl_positivity_shift_vlasov_advance_m0fix_cu_ker(
  struct gkyl_positivity_shift_vlasov_kernels *kers,
  const struct gkyl_range conf_range, const struct gkyl_array* GKYL_RESTRICT shiftedf,
  struct gkyl_array* GKYL_RESTRICT m0, struct gkyl_array* GKYL_RESTRICT delta_m0)
{
  int cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long clinidx = gkyl_range_idx(&conf_range, cidx);

    const int *shiftedf_c = (const int*) gkyl_array_cfetch(shiftedf, clinidx);
    double *delta_m0_c = (double*) gkyl_array_fetch(delta_m0, clinidx);

    if (shiftedf_c[0]) {
      double *m0_c = (double*) gkyl_array_fetch(m0, clinidx);
      if (kers->is_m0_positive(delta_m0_c)) {
	for (int k=0; k<m0->ncomp; k++) {
          m0_c[k] = delta_m0_c[k];
          delta_m0_c[k] = 0.0;
	}
      }
      else {
	for (int k=0; k<m0->ncomp; k++)
          delta_m0_c[k] = m0_c[k] - delta_m0_c[k];
      }
    }
    else {
      for (int k=0; k<m0->ncomp; k++)
        delta_m0_c[k] = 0.0;
    }
  }
}

void
gkyl_positivity_shift_vlasov_advance_cu(gkyl_positivity_shift_vlasov* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT m0,
  struct gkyl_array *GKYL_RESTRICT delta_m0)
{
  int nblocks_phase = phase_rng->nblocks, nthreads_phase = phase_rng->nthreads;
  int nblocks_conf = conf_rng->nblocks, nthreads_conf = conf_rng->nthreads;

  gkyl_array_clear_range(m0, 0.0, conf_rng);
  gkyl_array_clear_range(delta_m0, 0.0, conf_rng);

  // Set shiftedf boolean (int) to 0s.
  gkyl_positivity_shift_vlasov_advance_int_array_clear_cu_ker<<<nblocks_conf, nthreads_conf>>>
    (up->shiftedf->on_dev, 0);

  // Shift f is needed & scale f locally if initial local contribution to M0 was >0.
  gkyl_positivity_shift_vlasov_advance_shift_cu_ker<<<nblocks_phase, nthreads_phase>>>
    (up->kernels, up->grid, *conf_rng, *phase_rng, up->ffloor, up->ffloor_fac,
     up->cellav_fac, up->shiftedf->on_dev, distf->on_dev, m0->on_dev, delta_m0->on_dev);

  // If a shift took place, rescale f so it keeps the same M0.
  gkyl_positivity_shift_vlasov_advance_scalef_cu_ker<<<nblocks_phase, nthreads_phase>>>
    (up->kernels, *conf_rng, *phase_rng, up->shiftedf->on_dev, m0->on_dev, delta_m0->on_dev, distf->on_dev);

  // Ensure m0 and delta_m0 are correct based on whether a shift took place.
  gkyl_positivity_shift_vlasov_advance_m0fix_cu_ker<<<nblocks_conf, nthreads_conf>>>
    (up->kernels, *conf_rng, up->shiftedf->on_dev, m0->on_dev, delta_m0->on_dev);
}
