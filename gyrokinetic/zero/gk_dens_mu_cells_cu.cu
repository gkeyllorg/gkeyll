/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_gk_dens_mu_cells.h>
#include <gkyl_gk_dens_mu_cells_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
}

__global__ static void
gk_dens_mu_cells_set_cu_dev_ptrs(struct gkyl_gk_dens_mu_cells *up)
{
  // Set the pointer to the kernel on the device.
  up->kernel = gk_dens_mu_cells_choose_kernel(up->cdim);
}

void
gk_dens_mu_cells_dev_new(struct gkyl_gk_dens_mu_cells *up)
{
  struct gkyl_gk_dens_mu_cells *up_dev =
    (struct gkyl_gk_dens_mu_cells *)gkyl_cu_malloc(sizeof(*up_dev));
  gkyl_cu_memcpy(up_dev, up, sizeof(*up_dev), GKYL_CU_MEMCPY_H2D);
  gk_dens_mu_cells_set_cu_dev_ptrs<<<1, 1>>>(up_dev);
  up->on_dev = up_dev;
}

__global__ static void
gk_dens_mu_cells_advance_cu_kernel(
  const struct gkyl_gk_dens_mu_cells *up, struct gkyl_range phase_range,
  struct gkyl_range conf_range, const struct gkyl_array *fin, struct gkyl_array *out
)
{
  // Each thread computes the density in one mu cell of one configuration space cell, from the
  // phase space cells along vpar. No two threads write to the same memory.
  int cdim = up->cdim, pdim = up->pdim, num_basis = up->num_basis;
  int num_vpar_cells = phase_range.upper[cdim] - phase_range.lower[cdim] + 1;
  int num_mu_cells = phase_range.upper[pdim - 1] - phase_range.lower[pdim - 1] + 1;
  // Factor from the integral over mu of the vpar-integrated distribution.
  double mu_fac = 2.0 * M_PI / up->mass * 0.5 * up->phase_grid.dx[pdim - 1];

  int idx[GKYL_MAX_DIM];
  long num_tot = conf_range.volume * num_mu_cells;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    long conf_linc = tid / num_mu_cells;
    int mu_cell = tid % num_mu_cells;
    // The configuration space range is a sub-range, so use the inverse index of sub-ranges.
    gkyl_sub_range_inv_idx(&conf_range, conf_linc, idx);
    double *out_d = (double *)gkyl_array_fetch(out, gkyl_range_idx(&conf_range, idx));

    double cell_dens[GK_DENS_MU_MAX_BASIS] = {0.0};
    idx[pdim - 1] = phase_range.lower[pdim - 1] + mu_cell;
    for (int vpar_cell = 0; vpar_cell < num_vpar_cells; vpar_cell++) {
      idx[cdim] = phase_range.lower[cdim] + vpar_cell;
      const double *fin_d =
        (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&phase_range, idx));

      double vpar_int[GK_DENS_MU_MAX_CMU_BASIS] = {0.0};
      up->kernel(up->phase_grid.dx, 0, up->mass, 0, fin_d, vpar_int);

      for (int k = 0; k < num_basis; k++) {
        for (int j = 0; j < up->num_cmu_basis; j++) {
          cell_dens[k] += up->mu_weight[k][j] * vpar_int[j];
        }
      }
    }

    for (int k = 0; k < num_basis; k++) {
      out_d[mu_cell * num_basis + k] = mu_fac * cell_dens[k];
    }
  }
}

void
gkyl_gk_dens_mu_cells_advance_cu(
  const struct gkyl_gk_dens_mu_cells *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *fin, struct gkyl_array *out
)
{
  int num_mu_cells = phase_range->upper[up->pdim - 1] - phase_range->lower[up->pdim - 1] + 1;
  long num_tot = conf_range->volume * num_mu_cells;
  int num_blocks = (int)((num_tot + GKYL_DEFAULT_NUM_THREADS - 1) / GKYL_DEFAULT_NUM_THREADS);
  gk_dens_mu_cells_advance_cu_kernel<<<num_blocks, GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, *phase_range, *conf_range, fin->on_dev, out->on_dev
  );
}
