/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_harm_vol_priv.h>
}

// One thread per configuration cell: its stored operator and raw w-flux node values.
__global__ static void
gkyl_dg_pkpm_harm_vol_coef_cu_kernel(
  const struct gkyl_dg_pkpm_harm_vol *up, struct gkyl_range conf_range, const struct gkyl_array *bvar,
  const struct gkyl_array *grad_b, const struct gkyl_array *u, const struct gkyl_array *grad_u,
  const struct gkyl_array *xi, const struct gkyl_array *grad_xi, const struct gkyl_array *grad_ln_n,
  const struct gkyl_array *theta, const struct gkyl_array *dtb, const struct gkyl_array *em,
  struct gkyl_array *coef, struct gkyl_array *wnode
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    up->build_kern(up->qbym, (const double *)gkyl_array_cfetch(bvar, cidx),
      (const double *)gkyl_array_cfetch(grad_b, cidx), (const double *)gkyl_array_cfetch(u, cidx),
      (const double *)gkyl_array_cfetch(grad_u, cidx), (const double *)gkyl_array_cfetch(xi, cidx),
      (const double *)gkyl_array_cfetch(grad_xi, cidx), (const double *)gkyl_array_cfetch(grad_ln_n, cidx),
      (const double *)gkyl_array_cfetch(theta, cidx), (const double *)gkyl_array_cfetch(dtb, cidx),
      (const double *)gkyl_array_cfetch(em, cidx), (double *)gkyl_array_fetch(coef, cidx),
      (double *)gkyl_array_fetch(wnode, cidx));
  }
}

void
gkyl_dg_pkpm_harm_vol_coef_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef,
  struct gkyl_array *wnode
)
{
  gkyl_dg_pkpm_harm_vol_coef_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, bvar->on_dev, grad_b->on_dev, u->on_dev, grad_u->on_dev, xi->on_dev,
    grad_xi->on_dev, grad_ln_n->on_dev, theta->on_dev, dtB_over_magB->on_dev, em->on_dev, coef->on_dev,
    wnode->on_dev
  );
}

// One thread per configuration cell: its velocity cells from the stored coefficients.
__global__ static void
gkyl_dg_pkpm_harm_vol_advance_cu_kernel(
  const struct gkyl_dg_pkpm_harm_vol *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *coef, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *rhs
)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM];
  double mat[PKPM_HM_MAT];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    // the matrices of the stored operator of the configuration cell and the source rate of its column
    const double *ho = (const double *)gkyl_array_cfetch(coef, cidx);
    up->mat_kern(up->phase_grid.dx, (const double *)gkyl_array_cfetch(up->jacob_pos, cidx), ho, mat);
    double rate = pkpm_hm_src_bound(up, ho);
    for (int v = phase_range.lower[cdim]; v <= phase_range.upper[cdim]; ++v) {
      idx[cdim] = v;
      long pidx = gkyl_range_idx(&phase_range, idx);
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      up->vol_kern(xc, up->phase_grid.dx, mat, (const double *)gkyl_array_cfetch(fin, pidx),
        (double *)gkyl_array_fetch(rhs, pidx));
      double *cfl = (double *)gkyl_array_fetch(cflrate, pidx);
      cfl[0] += rate;
    }
  }
}

void
gkyl_dg_pkpm_harm_vol_advance_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  gkyl_dg_pkpm_harm_vol_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, coef->on_dev, fin->on_dev, cflrate->on_dev, rhs->on_dev
  );
}

// The device kernel pointers of the generated kernels.
__global__ static void
dg_pkpm_harm_vol_set_cu_dev_ptrs(struct gkyl_dg_pkpm_harm_vol *up, int cdim)
{
  up->build_kern = ser_pkpm_harm_build_kernels[cdim - 1];
  up->mat_kern = ser_pkpm_harm_mat_kernels[cdim - 1];
  up->vol_kern = ser_pkpm_harm_vol_kernels[cdim - 1];
}

struct gkyl_dg_pkpm_harm_vol *
gkyl_dg_pkpm_harm_vol_cu_dev_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  struct gkyl_dg_pkpm_harm_vol *up = (struct gkyl_dg_pkpm_harm_vol *)gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_harm_vol_setup(up, inp);
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  // The device copy needs the device pointer of the Jacobian array; the host object keeps the
  // acquired map for its own lifetime management.
  up->jacob_pos = inp->pos_map->jacob_pos->on_dev;
  struct gkyl_dg_pkpm_harm_vol *up_cu =
    (struct gkyl_dg_pkpm_harm_vol *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);
  dg_pkpm_harm_vol_set_cu_dev_ptrs<<<1, 1>>>(up_cu, up->cdim);
  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  return up;
}
