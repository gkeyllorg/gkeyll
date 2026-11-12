/* -*- c++ -*- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_dg_gaussian_filter.h>
#include <gkyl_dg_gaussian_filter_priv.h>
#include <gkyl_range.h>

#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(
  dim3 *dimGrid, dim3 *dimBlock, gkyl_range range, int ncomp
)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS / ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ static void
gkyl_apply_dg_gaussian_filter_ker(
  struct gkyl_rect_grid conf_grid, struct gkyl_range conf_range, int tot_quad,
  const struct gkyl_array *ordinates, const struct gkyl_array *weights,
  const struct gkyl_array *conf_arr_quad, struct gkyl_array *conf_arr_filter
)
{
  int cdim = conf_range.ndim;
  int cidx_l[GKYL_MAX_CDIM], cidx_c[GKYL_MAX_CDIM], cidx_r[GKYL_MAX_CDIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  double xmu[GKYL_MAX_DIM];
  double xmul[GKYL_MAX_DIM], xmuc[GKYL_MAX_DIM], xmur[GKYL_MAX_DIM];
  double dx[GKYL_MAX_CDIM];
  for (int d = 0; d < cdim; ++d) {
    dx[d] = conf_grid.dx[d];
  }

  // 2D thread grid
  // linc2 goes from 0 to tot_conf_quad
  long linc2 = threadIdx.y + blockIdx.y * blockDim.y;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < conf_range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx_c);

    gkyl_copy_int_arr(cdim, cidx_c, cidx_l);
    gkyl_copy_int_arr(cdim, cidx_c, cidx_r);
    // Assumes filter is only applied in first dimension for now JJ 09/28/25
    cidx_l[0] = cidx_l[0] - 1;
    cidx_r[0] = cidx_r[0] + 1;

    long loc_l = gkyl_range_idx(&conf_range, cidx_l);
    long loc_c = gkyl_range_idx(&conf_range, cidx_c);
    long loc_r = gkyl_range_idx(&conf_range, cidx_r);

    const double *conf_arr_quad_l = (const double *)gkyl_array_cfetch(conf_arr_quad, loc_l);
    const double *conf_arr_quad_c = (const double *)gkyl_array_cfetch(conf_arr_quad, loc_c);
    const double *conf_arr_quad_r = (const double *)gkyl_array_cfetch(conf_arr_quad, loc_r);
    double *conf_arr_filter_d = (double *)gkyl_array_fetch(conf_arr_filter, loc_c);

    // Fetch the cell center in the lower, center, and upper cells.
    gkyl_rect_grid_cell_center(&conf_grid, cidx_l, xcl);
    gkyl_rect_grid_cell_center(&conf_grid, cidx_c, xcc);
    gkyl_rect_grid_cell_center(&conf_grid, cidx_r, xcr);

    // Fetch the quadrature point location for the quadrature point we are updating.
    comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, linc2), dx, xcc, xmu);
    // First construct the filter normalization
    double denominator =
      filter_norm(cdim, tot_quad, ordinates, weights, dx, xmu, xcl, xmul, xcc, xmuc, xcr, xmur);
    for (int j = 0; j < tot_quad; ++j) {
      // Fetch the quadrature point locations for all other quadrature points.
      comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, j), dx, xcl, xmul);
      comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, j), dx, xcc, xmuc);
      comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, j), dx, xcr, xmur);
      const double *wq = (const double *)gkyl_array_cfetch(weights, j);
      // Sum the contributions from each cell in the three-cell filter.
      double expamp_l = -pow((xmu[0] - xmul[0]) / (sqrt(2.0) * dx[0]), 2.0);
      double expamp_c = -pow((xmu[0] - xmuc[0]) / (sqrt(2.0) * dx[0]), 2.0);
      double expamp_r = -pow((xmu[0] - xmur[0]) / (sqrt(2.0) * dx[0]), 2.0);
      conf_arr_filter_d[linc2] +=
        wq[0] / denominator *
        (conf_arr_quad_l[j] * exp(expamp_l) + conf_arr_quad_c[j] * exp(expamp_c) +
         conf_arr_quad_r[j] * exp(expamp_r));
    }
  }
}

void
gkyl_dg_gaussian_filter_advance_cu(
  gkyl_dg_gaussian_filter *up, const struct gkyl_range *conf_range, struct gkyl_array *conf_arr
)
{
  gkyl_array_clear(up->conf_arr_quad, 0.0);
  gkyl_array_clear(up->conf_arr_filter, 0.0);

  // Convert from modal basis to nodal quadrature basis.
  gkyl_mat_mm_array(up->m2n_mem, conf_arr, up->conf_arr_quad);

  dim3 dimGrid_conf, dimBlock_conf;
  gkyl_parallelize_components_kernel_launch_dims(
    &dimGrid_conf, &dimBlock_conf, *conf_range, up->tot_quad
  );
  gkyl_apply_dg_gaussian_filter_ker<<<dimGrid_conf, dimBlock_conf>>>(
    up->conf_grid, *conf_range, up->tot_quad, up->ordinates->on_dev, up->weights->on_dev,
    up->conf_arr_quad->on_dev, up->conf_arr_filter->on_dev
  );

  // Convert back to modal basis
  gkyl_mat_mm_array(up->n2m_mem, up->conf_arr_filter, conf_arr);
}