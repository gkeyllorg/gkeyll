#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gr_maxwell_lorentz_conf.h>
#include <gkyl_dg_gr_maxwell_lorentz_conf_priv.h>
#include <gkyl_util.h>

gkyl_dg_gr_maxwell_lorentz_conf*
gkyl_dg_gr_maxwell_lorentz_conf_inew(const struct gkyl_dg_gr_maxwell_lorentz_conf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_gr_maxwell_lorentz_conf_cu_dev_inew(inp);
  }
#endif

  struct gkyl_dg_gr_maxwell_lorentz_conf *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  int vdim = inp->vdim;
  int poly_order = inp->conf_basis->poly_order;
  int kernel_index = lorentz_cv_index[cdim].vdim[vdim];

  assert(kernel_index != -1);

  up->conf_grid = *inp->conf_grid;
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  up->gr_maxwell_data.chi = inp->chi;
  up->gr_maxwell_data.gamma = inp->gamma;
  up->gr_maxwell_data.K_phi = inp->K_phi;
  up->gr_maxwell_data.K_psi = inp->K_psi;
  up->cdim = cdim;
  up->vdim = vdim;
  up->num_basis = inp->conf_basis->num_basis;
  up->qbym = inp->qbym;
  up->use_gpu = inp->use_gpu;

  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->lorentz_conf = ser_lorentz_conf_kernels[kernel_index].kernels[poly_order];
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      up->lorentz_conf = tensor_lorentz_conf_kernels[kernel_index].kernels[poly_order];
      break;

    default:
      assert(false);
      break;
  }

  assert(up->lorentz_conf);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // Self-reference on host.

  return up;
}

void
gkyl_dg_gr_maxwell_lorentz_conf_advance(struct gkyl_dg_gr_maxwell_lorentz_conf *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h,
  const struct gkyl_surf_and_vol_node_arrays *vierb_cov, const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *field_con, struct gkyl_array *qmem)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(qmem)) {
    return gkyl_dg_gr_maxwell_lorentz_conf_advance_cu(up, conf_range,
      lapse, shift, h_ij, h_ij_inv, det_h, vierb_cov, vierb_con, field_con, qmem);
  }
#endif

  int num_basis = up->num_basis;
  double E_conf_for_force[3*num_basis], B_conf_for_force[3*num_basis];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);

    for (int k=0; k<3*num_basis; ++k) E_conf_for_force[k] = 0.0;
    for (int k=0; k<3*num_basis; ++k) B_conf_for_force[k] = 0.0;

    up->lorentz_conf(&up->gr_maxwell_data, up->conf_grid.dx,
      gkyl_array_cfetch(up->jacob_pos, cidx),
      gkyl_array_cfetch(lapse->nodal_arr_vol, cidx),
      gkyl_array_cfetch(shift->nodal_arr_vol, cidx),
      gkyl_array_cfetch(h_ij->nodal_arr_vol, cidx),
      gkyl_array_cfetch(h_ij_inv->nodal_arr_vol, cidx),
      gkyl_array_cfetch(det_h->nodal_arr_vol, cidx),
      gkyl_array_cfetch(field_con, cidx),
      gkyl_array_cfetch(vierb_cov->nodal_arr_vol, cidx),
      gkyl_array_cfetch(vierb_con->nodal_arr_vol, cidx),
      E_conf_for_force, B_conf_for_force);

    // qmem is cleared at the start of the RHS, has already added applied
    // acceleration or external EM contributions at this point, so these are additive
    double *qmem_d = gkyl_array_fetch(qmem, cidx);
    for (int k=0; k<num_basis; ++k) {
      qmem_d[0*num_basis+k] += up->qbym*E_conf_for_force[0*num_basis+k];
      qmem_d[1*num_basis+k] += up->qbym*E_conf_for_force[1*num_basis+k];
      qmem_d[2*num_basis+k] += up->qbym*E_conf_for_force[2*num_basis+k];
      qmem_d[3*num_basis+k] += up->qbym*B_conf_for_force[0*num_basis+k];
      qmem_d[4*num_basis+k] += up->qbym*B_conf_for_force[1*num_basis+k];
      qmem_d[5*num_basis+k] += up->qbym*B_conf_for_force[2*num_basis+k];
    }
  }
}

void
gkyl_dg_gr_maxwell_lorentz_conf_release(struct gkyl_dg_gr_maxwell_lorentz_conf* up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
