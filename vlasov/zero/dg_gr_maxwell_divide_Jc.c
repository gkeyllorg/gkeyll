#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
#include <gkyl_dg_gr_maxwell_divide_Jc_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>

void
gkyl_dg_gr_maxwell_divide_Jc(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_with_J_con,
  struct gkyl_array *field_no_J_con, bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    // change pointer type for det_h
    gkyl_dg_gr_maxwell_divide_Jc_cu(conf_basis, conf_range, det_h, field_with_J_con, field_no_J_con);
    return;
  }
#endif

  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  divide_Jc_t divide_Jc;
  switch (conf_basis->b_type) {
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

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  int idx_vel[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);

    const double *field_with_J_con_d = gkyl_array_cfetch(field_with_J_con, cidx);
    double *field_no_J_con_d = gkyl_array_fetch(field_no_J_con, cidx);
    divide_Jc(
      det_h ? gkyl_array_cfetch(det_h->nodal_arr_vol, cidx) : 0, field_with_J_con_d,
      field_no_J_con_d
    );
  }
}

void
gkyl_dg_gr_maxwell_rescale_Jc(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_no_J_con,
  struct gkyl_array *field_with_J_con, bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_dg_gr_maxwell_rescale_Jc_cu(
      conf_basis, conf_range, det_h, field_no_J_con, field_with_J_con
    );
    return;
  }
#endif

  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  rescale_Jc_t rescale_Jc;
  switch (conf_basis->b_type) {
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

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  int idx_vel[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);

    const double *field_no_J_con_d = gkyl_array_cfetch(field_no_J_con, cidx);
    double *field_with_J_con_d = gkyl_array_fetch(field_with_J_con, cidx);
    rescale_Jc(
      det_h ? gkyl_array_cfetch(det_h->nodal_arr_vol, cidx) : 0, field_no_J_con_d,
      field_with_J_con_d
    );
  }
}
