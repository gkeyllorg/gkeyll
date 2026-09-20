/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gr_maxwell_current_deposition.h>
#include <gkyl_dg_gr_maxwell_current_deposition_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>
}

#include <cassert>

__global__ void
gkyl_dg_gr_maxwell_current_deposition_advance_cu_kernel(
  struct gkyl_dg_gr_maxwell_current_deposition *up,
  const struct gkyl_range conf_range, double q_over_eps0,
  const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *m0, const struct gkyl_array *m1i,
  struct gkyl_array *rhs)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < conf_range.volume;
      linc += gridDim.x*blockDim.x) {

    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);

    up->current_deposition(q_over_eps0,
      (const double*) gkyl_array_cfetch(lapse->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(shift->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(vierb_con->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(m0, cidx),
      (const double*) gkyl_array_cfetch(m1i, cidx),
      (double*) gkyl_array_fetch(rhs, cidx));
  }
}

void
gkyl_dg_gr_maxwell_current_deposition_advance_cu(
  struct gkyl_dg_gr_maxwell_current_deposition *up,
  const struct gkyl_range *conf_range, double q_over_eps0,
  const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *m0, const struct gkyl_array *m1i,
  struct gkyl_array *rhs)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;

  gkyl_dg_gr_maxwell_current_deposition_advance_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, *conf_range, q_over_eps0, lapse->on_dev, shift->on_dev,
    vierb_con->on_dev, m0->on_dev, m1i->on_dev, rhs->on_dev);
}

__global__ static void
gkyl_dg_gr_maxwell_current_deposition_set_cu_dev_ptrs(
  struct gkyl_dg_gr_maxwell_current_deposition *up,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  int kernel_index = current_dep_cv_index[cdim].vdim[vdim];
  assert(kernel_index != -1);

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->current_deposition = ser_current_deposition_kernels[kernel_index].kernels[poly_order];
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      up->current_deposition = tensor_current_deposition_kernels[kernel_index].kernels[poly_order];
      break;

    default:
      assert(false);
      break;
  }
}

gkyl_dg_gr_maxwell_current_deposition*
gkyl_dg_gr_maxwell_current_deposition_cu_dev_inew(
  const struct gkyl_dg_gr_maxwell_current_deposition_inp *inp)
{
  struct gkyl_dg_gr_maxwell_current_deposition *up =
    (struct gkyl_dg_gr_maxwell_current_deposition*) gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  int vdim = inp->vdim;
  int poly_order = inp->conf_basis->poly_order;
  int kernel_index = current_dep_cv_index[cdim].vdim[vdim];

  assert(kernel_index != -1);

  up->cdim = cdim;
  up->vdim = vdim;
  up->num_basis = inp->conf_basis->num_basis;
  up->use_gpu = true;
  up->flags = 0;

  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_gr_maxwell_current_deposition *up_cu =
    (struct gkyl_dg_gr_maxwell_current_deposition*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_gr_maxwell_current_deposition), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_gr_maxwell_current_deposition_set_cu_dev_ptrs<<<1,1>>>(up_cu,
    inp->conf_basis->b_type, cdim, vdim, poly_order);

  up->on_dev = up_cu;

  return up;
}
