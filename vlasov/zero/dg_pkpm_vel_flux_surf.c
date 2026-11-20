#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_pkpm_vel_flux_surf_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_pkpm_vel_flux_surf_setup(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
)
{
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  assert(pdim == cdim + 1);
  assert(inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = inp->use_gpu;
  // The surface basis of the v_par face is the p=1 configuration basis.
  up->num_surf_basis = inp->conf_basis->num_basis;
}

gkyl_dg_pkpm_vel_flux_surf *
gkyl_dg_pkpm_vel_flux_surf_inew(const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_vel_flux_surf *up = gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_vel_flux_surf_setup(up, inp);
  up->vel_nodes = ser_pkpm_vel_flux_nodes_kernels[up->cdim - 1];
  up->vel_flux =
    choose_pkpm_vel_flux_kern(inp->conf_basis->b_type, up->cdim, inp->conf_basis->poly_order);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

void
gkyl_dg_pkpm_vel_flux_surf_advance(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *bvar, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *vel_flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vel_flux_surf)) {
    return gkyl_dg_pkpm_vel_flux_surf_advance_cu(
      up, conf_range, phase_range, wnode, theta, bvar, alpha_0, c_h, fin, vel_flux_surf
    );
  }
#endif
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  double nv[PKPM_VF_MAXNODE * PKPM_HS_NODE_NV];
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(conf_range, citer.idx);
    // the node values of the configuration cell, once for all its velocity faces
    up->vel_nodes(gkyl_array_cfetch(wnode, cidx), gkyl_array_cfetch(theta, cidx),
      gkyl_array_cfetch(bvar, cidx), gkyl_array_cfetch(alpha_0, cidx), gkyl_array_cfetch(c_h, cidx), nv);
    gkyl_copy_int_arr(cdim, citer.idx, idx);
    gkyl_copy_int_arr(cdim, citer.idx, idx_l);
    // the lower face of every cell but the lowest (the velocity domain edge, set by the closure)
    for (int v = phase_range->lower[cdim] + 1; v <= phase_range->upper[cdim]; ++v) {
      idx[cdim] = v;
      idx_l[cdim] = v - 1;
      long pidx = gkyl_range_idx(phase_range, idx), pidx_l = gkyl_range_idx(phase_range, idx_l);
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      up->vel_flux(xc, up->phase_grid.dx, nv, gkyl_array_cfetch(fin, pidx_l),
        gkyl_array_cfetch(fin, pidx), gkyl_array_fetch(vel_flux_surf, pidx));
    }
  }
}

int
gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_vel_flux_surf *up)
{
  return up->num_surf_basis;
}

void
gkyl_dg_pkpm_vel_flux_surf_release(struct gkyl_dg_pkpm_vel_flux_surf *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_vel_flux_surf *
gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_pkpm_vel_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *bvar, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *vel_flux_surf
)
{
  assert(false);
}

#endif
