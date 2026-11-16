/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_pkpm_vel_flux_surf_priv.h>
#include <gkyl_util.h>
}

// One thread per configuration cell: its node values, then the lower v_par faces of its velocity
// cells but the lowest (the lower velocity domain edge, set by the closure).
__global__ void
gkyl_dg_pkpm_vel_flux_surf_advance_cu_kernel(
  struct gkyl_dg_pkpm_vel_flux_surf *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *wnode, const struct gkyl_array *theta, const struct gkyl_array *bvar,
  const struct gkyl_array *alpha_0, const struct gkyl_array *c_h, const struct gkyl_array *fin,
  struct gkyl_array *vel_flux_surf
)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  double nv[PKPM_VF_MAXNODE * PKPM_HS_NODE_NV];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    up->vel_nodes((const double *)gkyl_array_cfetch(wnode, cidx),
      (const double *)gkyl_array_cfetch(theta, cidx), (const double *)gkyl_array_cfetch(bvar, cidx),
      (const double *)gkyl_array_cfetch(alpha_0, cidx), (const double *)gkyl_array_cfetch(c_h, cidx), nv);
    gkyl_copy_int_arr(cdim, idx, idx_l);
    for (int v = phase_range.lower[cdim] + 1; v <= phase_range.upper[cdim]; ++v) {
      idx[cdim] = v;
      idx_l[cdim] = v - 1;
      long pidx = gkyl_range_idx(&phase_range, idx), pidx_l = gkyl_range_idx(&phase_range, idx_l);
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      up->vel_flux(xc, up->phase_grid.dx, nv, (const double *)gkyl_array_cfetch(fin, pidx_l),
        (const double *)gkyl_array_cfetch(fin, pidx), (double *)gkyl_array_fetch(vel_flux_surf, pidx));
    }
  }
}

void
gkyl_dg_pkpm_vel_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *bvar, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *vel_flux_surf
)
{
  gkyl_dg_pkpm_vel_flux_surf_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, wnode->on_dev, theta->on_dev, bvar->on_dev, alpha_0->on_dev,
    c_h->on_dev, fin->on_dev, vel_flux_surf->on_dev
  );
}

__global__ static void
dg_pkpm_vel_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_pkpm_vel_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  up->vel_nodes = ser_pkpm_vel_flux_nodes_kernels[cdim - 1];
  up->vel_flux = choose_pkpm_vel_flux_kern(b_type, cdim, poly_order);
}

struct gkyl_dg_pkpm_vel_flux_surf *
gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp)
{
  struct gkyl_dg_pkpm_vel_flux_surf *up =
    (struct gkyl_dg_pkpm_vel_flux_surf *)gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_vel_flux_surf_setup(up, inp);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_pkpm_vel_flux_surf *up_cu =
    (struct gkyl_dg_pkpm_vel_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_vel_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, inp->conf_basis->b_type, inp->conf_basis->ndim, inp->conf_basis->poly_order
  );

  up->on_dev = up_cu;
  return up;
}
