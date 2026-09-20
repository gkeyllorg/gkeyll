/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gr_maxwell_lorentz_conf.h>
#include <gkyl_dg_gr_maxwell_lorentz_conf_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>
}

#include <cassert>

__global__ void
gkyl_dg_gr_maxwell_lorentz_conf_advance_cu_kernel(
  struct gkyl_dg_gr_maxwell_lorentz_conf *up,
  const struct gkyl_range conf_range,
  const struct gkyl_array *jacob_pos,
  const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij,
  const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h,
  const struct gkyl_surf_and_vol_node_arrays *vierb_cov,
  const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *field_con, struct gkyl_array *qmem)
{
  int idx[GKYL_MAX_DIM];
  int num_basis = up->num_basis;
  assert(3*num_basis <= 96);

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < conf_range.volume;
      linc += gridDim.x*blockDim.x) {

    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);

    double E_conf_for_force[96], B_conf_for_force[96];
    for (int k=0; k<3*num_basis; ++k) E_conf_for_force[k] = 0.0;
    for (int k=0; k<3*num_basis; ++k) B_conf_for_force[k] = 0.0;

    up->lorentz_conf(&up->gr_maxwell_data, up->conf_grid.dx,
      (const double*) gkyl_array_cfetch(jacob_pos, cidx),
      (const double*) gkyl_array_cfetch(lapse->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(shift->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(h_ij->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(h_ij_inv->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(det_h->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(field_con, cidx),
      (const double*) gkyl_array_cfetch(vierb_cov->nodal_arr_vol, cidx),
      (const double*) gkyl_array_cfetch(vierb_con->nodal_arr_vol, cidx),
      E_conf_for_force, B_conf_for_force);

    double *qmem_d = (double*) gkyl_array_fetch(qmem, cidx);
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
gkyl_dg_gr_maxwell_lorentz_conf_advance_cu(struct gkyl_dg_gr_maxwell_lorentz_conf *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h,
  const struct gkyl_surf_and_vol_node_arrays *vierb_cov, const struct gkyl_surf_and_vol_node_arrays *vierb_con,
  const struct gkyl_array *field_con, struct gkyl_array *qmem)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;

  gkyl_dg_gr_maxwell_lorentz_conf_advance_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, *conf_range, up->jacob_pos->on_dev, lapse->on_dev, shift->on_dev, h_ij->on_dev,
    h_ij_inv->on_dev, det_h->on_dev, vierb_cov->on_dev, vierb_con->on_dev,
    field_con->on_dev, qmem->on_dev);
}

__global__ static void
gkyl_dg_gr_maxwell_lorentz_conf_set_cu_dev_ptrs(
  struct gkyl_dg_gr_maxwell_lorentz_conf *up,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  int kernel_index = lorentz_cv_index[cdim].vdim[vdim];
  assert(kernel_index != -1);

  switch (b_type) {
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
}

gkyl_dg_gr_maxwell_lorentz_conf*
gkyl_dg_gr_maxwell_lorentz_conf_cu_dev_inew(
  const struct gkyl_dg_gr_maxwell_lorentz_conf_inp *inp)
{
  struct gkyl_dg_gr_maxwell_lorentz_conf *up =
    (struct gkyl_dg_gr_maxwell_lorentz_conf*) gkyl_malloc(sizeof(*up));

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
  up->use_gpu = true;
  up->flags = 0;

  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_gr_maxwell_lorentz_conf *up_cu =
    (struct gkyl_dg_gr_maxwell_lorentz_conf*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_gr_maxwell_lorentz_conf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_gr_maxwell_lorentz_conf_set_cu_dev_ptrs<<<1,1>>>(up_cu,
    inp->conf_basis->b_type, cdim, vdim, poly_order);

  up->on_dev = up_cu;

  return up;
}
