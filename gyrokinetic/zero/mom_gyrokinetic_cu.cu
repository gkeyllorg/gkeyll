/* -*- c++ -*- */

extern "C" {
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
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
    kernels->wmom(phase_grid.dx, vmap_c, mass, charge, bmag_c, phi_c, f_c, mom_local);

    double *mout_c = (double *)gkyl_array_fetch(mout, linidx_conf);
    for (unsigned int k = 0; k < mout->ncomp; ++k) {
      atomicAdd(&mout_c[k], mom_local[k]);
    }
  }
}

void
gkyl_mom_gyrokinetic_advance_cu(
  struct gkyl_mom_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
)
{
  int nblocks = phase_rng->nblocks, nthreads = phase_rng->nthreads;

  gkyl_array_clear_range(mout, 0.0, conf_rng);

  if (phase_rng->volume == 0 || conf_rng->volume == 0) {
    return;
  }

  struct gkyl_array *phi_on_dev = NULL;
  if (phi) {
    phi_on_dev = phi->on_dev;
  }

  gkyl_mom_gk_advance_cu_ker<<<nblocks, nthreads>>>(
    up->mass, up->charge, up->kernels, up->phase_grid, *conf_rng, up->vel_map->local_vel,
    *phase_rng, up->gk_geom->geo_corn.bmag->on_dev, up->vel_map->vmap->on_dev, phi_on_dev,
    fin->on_dev, mout->on_dev
  );
}
