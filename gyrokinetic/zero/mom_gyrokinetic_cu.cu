/* -*- c++ -*- */

extern "C" {
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
#include <cub/block/block_reduce.cuh>
#endif

// Constrain the indirect call to this interface's registry. The runtime pointer
// still selects the kernel; unrelated functions with the same PTX signature are
// not candidates. Every selected kernel must occur in the registry.
GKYL_CU_D static __forceinline__ void
mom_gk_call(
  wmomf_t kernel, const double *dxv, const double *vmap, double mass, double charge,
  const double *bmag, const double *phi, const double *f, double *out
)
{
  asm volatile("{\n"
               " .param .b64 a0, a1, a2, a3, a4, a5, a6, a7;\n"
               " st.param.b64 [a0], %1;\n"
               " st.param.b64 [a1], %2;\n"
               " st.param.f64 [a2], %3;\n"
               " st.param.f64 [a3], %4;\n"
               " st.param.b64 [a4], %5;\n"
               " st.param.b64 [a5], %6;\n"
               " st.param.b64 [a6], %7;\n"
               " st.param.b64 [a7], %8;\n"
               " targets: .calltargets " GKYL_MOM_GK_CALL_TARGETS ";\n"
               " call %0, (a0, a1, a2, a3, a4, a5, a6, a7), targets;\n"
               "}\n"
               :
               : "l"(kernel), "l"(dxv), "l"(vmap), "d"(mass), "d"(charge), "l"(bmag), "l"(phi),
                 "l"(f), "l"(out)
               : "memory");
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_mom_gk_set_cu_ker_ptrs(
  struct gkyl_mom_gyrokinetic_kernels *kernels, int *num_mom, const struct gkyl_basis cbasis,
  const struct gkyl_basis pbasis, enum gkyl_distribution_moments mom_type, bool is_integrated
)
{
  mom_gk_choose_kernel(kernels, num_mom, &cbasis, &pbasis, mom_type, is_integrated);
};

void
mom_gk_choose_kernel_cu(
  struct gkyl_mom_gyrokinetic_kernels *kernels, int *num_mom, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, enum gkyl_distribution_moments mom_type, bool is_integrated
)
{
  // Allocate num_mom on device so we can copy it back to host.
  int *num_mom_dev = (int *)gkyl_cu_malloc(sizeof(int));

  gkyl_mom_gk_set_cu_ker_ptrs<<<1, 1>>>(
    kernels, num_mom_dev, *cbasis, *pbasis, mom_type, is_integrated
  );

  gkyl_cu_memcpy(num_mom, num_mom_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(num_mom_dev);
}

#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
__global__ static void
gkyl_mom_gk_advance_cu_ker(
  double mass, double charge, struct gkyl_mom_gyrokinetic_kernels *kernels,
  struct gkyl_rect_grid phase_grid, struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  struct gkyl_range phase_rng, const struct gkyl_array *GKYL_RESTRICT bmag,
  const struct gkyl_array *GKYL_RESTRICT vmap, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
)
{
  int pidx[GKYL_MAX_DIM], vidx[GKYL_MAX_VDIM];
  int cdim = conf_rng.ndim;

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < phase_rng.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_rng, tid, pidx);
    if (!gkyl_range_contains_idx(&conf_rng, pidx)) {
      continue;
    }

    for (int d = cdim; d < phase_rng.ndim; d++) {
      vidx[d - cdim] = pidx[d];
    }

    long linidx_conf = gkyl_range_idx(&conf_rng, pidx);
    long linidx_vel = gkyl_range_idx(&vel_rng, vidx);
    long linidx_phase = gkyl_range_idx(&phase_rng, pidx);

    const double *phi_c = 0;
    if (phi) {
      phi_c = (const double *)gkyl_array_cfetch(phi, linidx_conf);
    }
    const double *bmag_c = (const double *)gkyl_array_cfetch(bmag, linidx_conf);
    const double *vmap_c = (const double *)gkyl_array_cfetch(vmap, linidx_vel);
    const double *f_c = (const double *)gkyl_array_cfetch(fin, linidx_phase);

    double mom_local[48]; // hard-coded to 6 moments * 8 configuration coefficients for now.
    for (unsigned int k = 0; k < 48; ++k) {
      mom_local[k] = 0.0;
    }

    // Reduce local f to local mom.
    mom_gk_call(kernels->wmom, phase_grid.dx, vmap_c, mass, charge, bmag_c, phi_c, f_c, mom_local);

    double *mout_c = (double *)gkyl_array_fetch(mout, linidx_conf);
    for (unsigned int k = 0; k < mout->ncomp; ++k) {
      atomicAdd(&mout_c[k], mom_local[k]);
    }
  }
}
#else
enum { MOM_GK_BLOCK_SIZE = 128 };

__global__ static void
mom_gk_advance_deterministic_cu_ker(
  double mass, double charge, struct gkyl_mom_gyrokinetic_kernels *kernels,
  struct gkyl_rect_grid phase_grid, struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  struct gkyl_range phase_rng, const struct gkyl_array *bmag, const struct gkyl_array *vmap,
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *mout, long num_vel
)
{
  typedef cub::BlockReduce<double, MOM_GK_BLOCK_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  int cidx[GKYL_MAX_CDIM], pidx[GKYL_MAX_DIM], vidx[GKYL_MAX_VDIM];
  gkyl_sub_range_inv_idx(&conf_rng, blockIdx.x, cidx);
  for (int d = 0; d < conf_rng.ndim; ++d) {
    pidx[d] = cidx[d];
  }
  long clin = gkyl_range_idx(&conf_rng, cidx);
  const double *bmag_c = (const double *)gkyl_array_cfetch(bmag, clin);
  const double *phi_c = phi ? (const double *)gkyl_array_cfetch(phi, clin) : NULL;
  double local[48] = {0.0};
  for (long v = threadIdx.x; v < num_vel; v += MOM_GK_BLOCK_SIZE) {
    long remaining = v;
    for (int d = conf_rng.ndim; d < phase_rng.ndim; ++d) {
      int extent = phase_rng.upper[d] - phase_rng.lower[d] + 1;
      pidx[d] = phase_rng.lower[d] + remaining % extent;
      vidx[d - conf_rng.ndim] = pidx[d];
      remaining /= extent;
    }
    double moment[48] = {0.0};
    mom_gk_call(
      kernels->wmom, phase_grid.dx,
      (const double *)gkyl_array_cfetch(vmap, gkyl_range_idx(&vel_rng, vidx)), mass, charge, bmag_c,
      phi_c, (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&phase_rng, pidx)), moment
    );
    for (unsigned k = 0; k < mout->ncomp; ++k) {
      local[k] += moment[k];
    }
  }
  double *out = (double *)gkyl_array_fetch(mout, clin);
  for (unsigned k = 0; k < mout->ncomp; ++k) {
    double sum = BlockReduceT(temp).Sum(local[k]);
    if (threadIdx.x == 0) {
      out[k] = sum;
    }
    __syncthreads();
  }
}
#endif

void
gkyl_mom_gyrokinetic_advance_cu(
  struct gkyl_mom_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
)
{
  gkyl_array_clear_range(mout, 0.0, conf_rng);

  if (phase_rng->volume == 0 || conf_rng->volume == 0) {
    return;
  }

  struct gkyl_array *phi_on_dev = NULL;
  if (phi) {
    phi_on_dev = phi->on_dev;
  }

#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  long num_vel = 1;
  for (int d = conf_rng->ndim; d < phase_rng->ndim; ++d) {
    num_vel *= phase_rng->upper[d] - phase_rng->lower[d] + 1;
  }
  mom_gk_advance_deterministic_cu_ker<<<conf_rng->volume, MOM_GK_BLOCK_SIZE>>>(
    up->mass, up->charge, up->kernels, up->phase_grid, *conf_rng, up->vel_map->local_vel,
    *phase_rng, up->gk_geom->geo_corn.bmag->on_dev, up->vel_map->vmap->on_dev, phi_on_dev,
    fin->on_dev, mout->on_dev, num_vel
  );
#else
  gkyl_mom_gk_advance_cu_ker<<<phase_rng->nblocks, phase_rng->nthreads>>>(
    up->mass, up->charge, up->kernels, up->phase_grid, *conf_rng, up->vel_map->local_vel,
    *phase_rng, up->gk_geom->geo_corn.bmag->on_dev, up->vel_map->vmap->on_dev, phi_on_dev,
    fin->on_dev, mout->on_dev
  );
#endif
}
