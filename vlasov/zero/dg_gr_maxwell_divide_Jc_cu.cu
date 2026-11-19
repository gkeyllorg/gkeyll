/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
#include <gkyl_dg_gr_maxwell_divide_Jc_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_gr_maxwell_divide_Jc_cu_kernel(
  const struct gkyl_basis conf_basis, const struct gkyl_range conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_with_J_con,
  struct gkyl_array *field_no_J_con
)
{
  int cdim = conf_basis.ndim;
  int poly_order = conf_basis.poly_order;
  divide_Jc_t divide_Jc;
  switch (conf_basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      divide_Jc = choose_ser_divide_Jc_kern(cdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      divide_Jc = choose_tensor_divide_Jc_kern(cdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  int idx[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    long cidx = gkyl_range_idx(&conf_range, idx);

    const double *field_with_J_con_d = (const double *)gkyl_array_cfetch(field_with_J_con, cidx);
    double *field_no_J_con_d = (double *)gkyl_array_fetch(field_no_J_con, cidx);
    divide_Jc(
      det_h ? (const double *)gkyl_array_cfetch(det_h->nodal_arr_vol, cidx) : 0, field_with_J_con_d,
      field_no_J_con_d
    );
  }
}

// Host-side wrapper for dividing out velocity-space Jacobian from Jf.
void
gkyl_dg_gr_maxwell_divide_Jc_cu(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_with_J_con,
  struct gkyl_array *field_no_J_con
)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_gr_maxwell_divide_Jc_cu_kernel<<<nblocks, nthreads>>>(
    *conf_basis, *conf_range, det_h ? det_h->on_dev : 0, field_with_J_con->on_dev,
    field_no_J_con->on_dev
  );
}

__global__ void
gkyl_dg_gr_maxwell_rescale_Jc_cu_kernel(
  const struct gkyl_basis conf_basis, const struct gkyl_range conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_no_J_con,
  struct gkyl_array *field_with_J_con
)
{
  int cdim = conf_basis.ndim;
  int poly_order = conf_basis.poly_order;
  rescale_Jc_t rescale_Jc;
  switch (conf_basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      rescale_Jc = choose_ser_rescale_Jc_kern(cdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      rescale_Jc = choose_tensor_rescale_Jc_kern(cdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  int idx[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    long cidx = gkyl_range_idx(&conf_range, idx);

    const double *field_no_J_con_d = (const double *)gkyl_array_cfetch(field_no_J_con, cidx);
    double *field_with_J_con_d = (double *)gkyl_array_fetch(field_with_J_con, cidx);
    rescale_Jc(
      det_h ? (const double *)gkyl_array_cfetch(det_h->nodal_arr_vol, cidx) : 0, field_no_J_con_d,
      field_with_J_con_d
    );
  }
}

// Host-side wrapper for multiplying conf-space Jacobian by field to obtain Jfield.
void
gkyl_dg_gr_maxwell_rescale_Jc_cu(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_no_J_con,
  struct gkyl_array *field_with_J_con
)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_gr_maxwell_rescale_Jc_cu_kernel<<<nblocks, nthreads>>>(
    *conf_basis, *conf_range, det_h ? det_h->on_dev : 0, field_no_J_con->on_dev,
    field_with_J_con->on_dev
  );
}
