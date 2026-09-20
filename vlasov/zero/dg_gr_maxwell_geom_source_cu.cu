/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_gr_maxwell_geom_source.h>
#include <gkyl_dg_gr_maxwell_geom_source_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>
}

#include <cassert>

__global__ void
gkyl_dg_gr_maxwell_geom_source_advance_cu_kernel(struct gkyl_dg_gr_maxwell_geom_source *up,
  const struct gkyl_range conf_range,
  const struct gkyl_surf_and_vol_node_arrays *geom_factor_con,
  const struct gkyl_array *field_con, struct gkyl_array *rhs)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < conf_range.volume;
      linc += gridDim.x*blockDim.x) {

    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xc);

    const double *geom_factor_con_d = (const double*) gkyl_array_cfetch(geom_factor_con->nodal_arr_vol, cidx);
    const double *field_con_d = (const double*) gkyl_array_cfetch(field_con, cidx);
    double *rhs_d = (double*) gkyl_array_fetch(rhs, cidx);

    up->geom_source(&up->gr_maxwell_data, xc, up->conf_grid.dx,
      geom_factor_con_d, field_con_d, rhs_d);
  }
}

void
gkyl_dg_gr_maxwell_geom_source_advance_cu(struct gkyl_dg_gr_maxwell_geom_source *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *geom_factor_con,
  const struct gkyl_array *field_con, struct gkyl_array *rhs)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_gr_maxwell_geom_source_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev,
    *conf_range, geom_factor_con->on_dev, field_con->on_dev, rhs->on_dev);
}

// CUDA kernel to set device pointers to kernel functions.
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol.
__global__ static void
gkyl_dg_gr_maxwell_geom_source_set_cu_dev_ptrs(struct gkyl_dg_gr_maxwell_geom_source *up,
  enum gkyl_basis_type b_type, int cdim, int poly_order, enum gkyl_field_id field_id)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->geom_source = ser_geom_source_kernels[cdim-1].kernels[poly_order];
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      up->geom_source = ten_geom_source_kernels[cdim-1].kernels[poly_order];
      break;

    default:
      assert(false);
      break;
  }
}

gkyl_dg_gr_maxwell_geom_source*
gkyl_dg_gr_maxwell_geom_source_cu_dev_inew(const struct gkyl_dg_gr_maxwell_geom_source_inp *inp)
{
  struct gkyl_dg_gr_maxwell_geom_source *up =
    (struct gkyl_dg_gr_maxwell_geom_source*) gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;

  up->conf_grid = *inp->conf_grid;
  up->cdim = cdim;
  up->gr_maxwell_data.chi = inp->chi;
  up->gr_maxwell_data.gamma = inp->gamma;
  up->gr_maxwell_data.K_phi = inp->K_phi;
  up->gr_maxwell_data.K_psi = inp->K_psi;
  up->use_gpu = true;
  up->flags = 0;

  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_gr_maxwell_geom_source *up_cu =
    (struct gkyl_dg_gr_maxwell_geom_source*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_gr_maxwell_geom_source), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_gr_maxwell_geom_source_set_cu_dev_ptrs<<<1,1>>>(up_cu,
    inp->conf_basis->b_type, cdim, poly_order, inp->field_id);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;
}
