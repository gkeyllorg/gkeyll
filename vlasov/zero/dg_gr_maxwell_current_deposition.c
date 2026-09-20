#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gr_maxwell_current_deposition.h>
#include <gkyl_dg_gr_maxwell_current_deposition_priv.h>
#include <gkyl_util.h>

gkyl_dg_gr_maxwell_current_deposition*
gkyl_dg_gr_maxwell_current_deposition_inew(
  const struct gkyl_dg_gr_maxwell_current_deposition_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_gr_maxwell_current_deposition_cu_dev_inew(inp);
  }
#endif

  struct gkyl_dg_gr_maxwell_current_deposition *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  int vdim = inp->vdim;
  int poly_order = inp->conf_basis->poly_order;
  int kernel_index = current_dep_cv_index[cdim].vdim[vdim];

  assert(kernel_index != -1);

  up->cdim = cdim;
  up->vdim = vdim;
  up->num_basis = inp->conf_basis->num_basis;
  up->use_gpu = inp->use_gpu;

  switch (inp->conf_basis->b_type) {
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

  assert(up->current_deposition);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // Self-reference on host.

  return up;
}

void
gkyl_dg_gr_maxwell_current_deposition_advance(
  struct gkyl_dg_gr_maxwell_current_deposition *up,
  const struct gkyl_range *conf_range, double q_over_eps0,
  const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *m0, const struct gkyl_array *m1i,
  struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_gr_maxwell_current_deposition_advance_cu(up, conf_range,
      q_over_eps0, lapse, shift, vierb_con, m0, m1i, rhs);
  }
#endif

  assert(up);
  assert(up->current_deposition);
  assert(m0->ncomp >= up->num_basis);
  assert(m1i->ncomp >= up->vdim*up->num_basis);
  assert(rhs->ncomp >= 3*up->num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);

    up->current_deposition(q_over_eps0,
      gkyl_array_cfetch(lapse->nodal_arr_vol, cidx),
      gkyl_array_cfetch(shift->nodal_arr_vol, cidx),
      gkyl_array_cfetch(vierb_con->nodal_arr_vol, cidx),
      gkyl_array_cfetch(m0, cidx),
      gkyl_array_cfetch(m1i, cidx),
      gkyl_array_fetch(rhs, cidx));
  }
}

void
gkyl_dg_gr_maxwell_current_deposition_release(
  struct gkyl_dg_gr_maxwell_current_deposition *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
