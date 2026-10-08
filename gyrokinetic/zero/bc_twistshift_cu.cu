/* -*- c++ -*- */
// All of this will be useless when we'll have the DG average updater for phase space.

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_bc_twistshift_shift_dir_avg_cu_ker(
  int sdir, const int *GKYL_RESTRICT shift_indep, struct gkyl_range ghost_r,
  struct gkyl_range ghost_avg_r, const struct gkyl_array *GKYL_RESTRICT fprolong,
  struct gkyl_array *GKYL_RESTRICT favg
)
{
  // One thread per cell of the ghost plane collapsed along shift_dir.
  int idx[GKYL_MAX_DIM];
  int num_cells_shift = ghost_r.upper[sdir] - ghost_r.lower[sdir] + 1;

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < ghost_avg_r.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&ghost_avg_r, tid, idx);

    for (int k = 0; k < fprolong->ncomp; k++) {
      if (!shift_indep[k]) {
        continue;
      }
      double avg = 0.0;
      for (int i = ghost_r.lower[sdir]; i <= ghost_r.upper[sdir]; i++) {
        idx[sdir] = i;
        const double *f_c =
          (const double *)gkyl_array_cfetch(fprolong, gkyl_range_idx(&ghost_r, idx));
        avg += f_c[k];
      }
      avg /= num_cells_shift;
      for (int i = ghost_r.lower[sdir]; i <= ghost_r.upper[sdir]; i++) {
        idx[sdir] = i;
        double *favg_c = (double *)gkyl_array_fetch(favg, gkyl_range_idx(&ghost_r, idx));
        favg_c[k] = avg;
      }
    }
  }
}

void
gkyl_bc_twistshift_shift_dir_avg_cu(struct gkyl_bc_twistshift *up)
{
  int nblocks = up->ghost_avg_r.nblocks, nthreads = up->ghost_avg_r.nthreads;
  gkyl_bc_twistshift_shift_dir_avg_cu_ker<<<nblocks, nthreads>>>(
    up->shift_dir, up->shift_indep_cu, up->ghost_r, up->ghost_avg_r, up->fprolong->on_dev,
    up->favg->on_dev
  );
}
