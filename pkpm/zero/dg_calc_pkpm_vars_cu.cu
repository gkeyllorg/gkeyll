/* -*- c++ -*- */

#include <math.h>
#include <string.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_calc_pkpm_vars_closure_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  double w_lo, double w_up, double dv, const struct gkyl_array *norm_moms,
  const struct gkyl_array *rhs_moms, const struct gkyl_array *bvar, const struct gkyl_array *wnode,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *c_0,
  struct gkyl_array *alpha_0, struct gkyl_array *c_h, struct gkyl_array *a_perp,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate, int *num_null,
  unsigned long long *cond_min
)
{
  int cdim = conf_range.ndim;
  int idx[GKYL_MAX_DIM];
  double cfl[2], cond[1];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    idx[cdim] = phase_range.lower[cdim];
    long plo = gkyl_range_idx(&phase_range, idx);
    idx[cdim] = phase_range.upper[cdim];
    long pup = gkyl_range_idx(&phase_range, idx);
    idx[cdim] = phase_range.upper[cdim] + 1; // velocity ghost cell: its slot is the upper edge face
    long pgh = gkyl_range_idx(&phase_range, idx);
    int null = up->pkpm_closure(
      w_lo, w_up, dv, (const double *)gkyl_array_cfetch(norm_moms, loc),
      (const double *)gkyl_array_cfetch(rhs_moms, loc),
      (const double *)gkyl_array_cfetch(bvar, loc), (const double *)gkyl_array_cfetch(wnode, loc),
      (const double *)gkyl_array_cfetch(theta, loc),
      (const double *)gkyl_array_cfetch(alpha_0, loc), (const double *)gkyl_array_cfetch(c_h, loc),
      (const double *)gkyl_array_cfetch(fin, plo), (const double *)gkyl_array_cfetch(fin, pup),
      (double *)gkyl_array_fetch(c_0, loc), (double *)gkyl_array_fetch(alpha_0, loc),
      (double *)gkyl_array_fetch(c_h, loc), (double *)gkyl_array_fetch(a_perp, loc),
      (double *)gkyl_array_fetch(vel_flux_surf, plo),
      (double *)gkyl_array_fetch(vel_flux_surf, pgh), cfl, cond
    );
    // the column rate of the velocity-space terms, the same in every cell of the column
    for (int v = phase_range.lower[cdim]; v <= phase_range.upper[cdim]; ++v) {
      idx[cdim] = v;
      double *c = (double *)gkyl_array_fetch(cflrate, gkyl_range_idx(&phase_range, idx));
      c[0] += cfl[0];
    }
    if (null > 0) {
      atomicAdd(num_null, null);
    }
    // minimum of the conditioning: for non-negative doubles the bit pattern orders like the
    // value; a negative ratio (a nulled node) is recorded as 0
    double cpos = cond[0] > 0.0 ? cond[0] : 0.0;
    atomicMin(cond_min, (unsigned long long)__double_as_longlong(cpos));
  }
}

// Host-side wrapper for the closures
int
gkyl_dg_calc_pkpm_vars_closure_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0,
  struct gkyl_array *c_h, struct gkyl_array *a_perp, struct gkyl_array *vel_flux_surf,
  struct gkyl_array *cflrate
)
{
  int cdim = conf_range->ndim;
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  unsigned long long one = 0x3FF0000000000000ULL; // bit pattern of 1.0
  gkyl_cu_memcpy(up->closure_cond_min_dev, &one, sizeof(unsigned long long), GKYL_CU_MEMCPY_H2D);
  gkyl_dg_calc_pkpm_vars_closure_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, phase_grid->lower[cdim], phase_grid->upper[cdim],
    phase_grid->dx[cdim], norm_moms->on_dev, rhs_moms->on_dev, bvar->on_dev, wnode->on_dev,
    theta->on_dev, fin->on_dev, c_0->on_dev, alpha_0->on_dev, c_h->on_dev, a_perp->on_dev,
    vel_flux_surf->on_dev, cflrate->on_dev, up->num_null_rho_dev, up->closure_cond_min_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  unsigned long long cmin = 0;
  gkyl_cu_memcpy(&cmin, up->closure_cond_min_dev, sizeof(unsigned long long), GKYL_CU_MEMCPY_D2H);
  double cmin_d;
  memcpy(&cmin_d, &cmin, sizeof(double));
  up->closure_cond_min = cmin_d;
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_prim_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  struct gkyl_array *u, struct gkyl_array *xi, int *num_null
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    int null = up->pkpm_prim(
      (const double *)gkyl_array_cfetch(fluid, loc), (double *)gkyl_array_fetch(u, loc),
      (double *)gkyl_array_fetch(xi, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

// Host-side wrapper for the nodal primitive variables
int
gkyl_dg_calc_pkpm_vars_prim_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
)
{
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_prim_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, fluid->on_dev, u->on_dev, xi->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_int_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *u,
  struct gkyl_array *int_pkpm_vars
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    up->pkpm_int(
      (const double *)gkyl_array_cfetch(fluid, loc),
      (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, loc),
      (const double *)gkyl_array_cfetch(u, loc), (double *)gkyl_array_fetch(int_pkpm_vars, loc)
    );
  }
}

// Host-side wrapper for the integrated diagnostics
void
gkyl_dg_calc_pkpm_vars_int_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
)
{
  gkyl_dg_calc_pkpm_vars_int_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, fluid->on_dev, vlasov_pkpm_moms->on_dev, u->on_dev,
    int_pkpm_vars->on_dev
  );
}

__global__ void
gkyl_dg_calc_pkpm_vars_grad_ln_n_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *rho,
  const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n, int *num_null
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);

    int null = up->pkpm_grad_ln_n(
      (const double *)gkyl_array_cfetch(rho, loc), (const double *)gkyl_array_cfetch(grad_rho, loc),
      (double *)gkyl_array_fetch(grad_ln_n, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

int null = up->pkpm_grad_ln_n(
  (const double *)gkyl_array_cfetch(rho, loc), (const double *)gkyl_array_cfetch(grad_rho, loc),
  (double *)gkyl_array_fetch(grad_ln_n, loc)
);
if (null > 0) {
  atomicAdd(num_null, null);
}
}
}
// Host-side wrapper for grad(ln n)
int
gkyl_dg_calc_pkpm_vars_grad_ln_n_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
)
{
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_grad_ln_n_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, rho->on_dev, grad_rho->on_dev, grad_ln_n->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_positivity_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, double eps,
  struct gkyl_array *fluid, struct gkyl_array *count, int *num_pos
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);

    int flag = up->pkpm_fluid_positivity(eps, (double *)gkyl_array_fetch(fluid, loc));
    double *cnt = (double *)gkyl_array_fetch(count, loc);
    if (flag < 0) {
      cnt[3] += 1.0;
      atomicAdd(&num_pos[3], 1);
    }
    if (flag > 0 && (flag & 1)) {
      cnt[0] += 1.0;
      atomicAdd(&num_pos[0], 1);
    }
    if (flag > 0 && (flag & 2)) {
      cnt[1] += 1.0;
      atomicAdd(&num_pos[1], 1);
    }
    if (flag > 0 && (flag & 4)) {
      cnt[2] += 1.0;
      atomicAdd(&num_pos[2], 1);
    }
  }
}
// Host-side wrapper for the fluid positivity limiter
void
gkyl_dg_calc_pkpm_vars_positivity_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, double eps,
  struct gkyl_array *fluid, struct gkyl_array *count, int *num_rho, int *num_e, int *num_xi,
  int *num_avg
)
{
  gkyl_cu_memset(up->num_pos_dev, 0, sizeof(int[4]));
  gkyl_dg_calc_pkpm_vars_positivity_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, eps, fluid->on_dev, count->on_dev, up->num_pos_dev
  );
  int num_pos[4] = {0, 0, 0, 0};
  gkyl_cu_memcpy(num_pos, up->num_pos_dev, sizeof(int[4]), GKYL_CU_MEMCPY_D2H);
  *num_rho = num_pos[0];
  *num_e = num_pos[1];
  *num_xi = num_pos[2];
  *num_avg = num_pos[3];
}

__global__ void
gkyl_dg_calc_pkpm_vars_normalize_dist_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *n, const struct gkyl_array *f, struct gkyl_array *fbar, int *num_null
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);

    int null = up->pkpm_normalize(
      (const double *)gkyl_array_cfetch(n, cidx), (const double *)gkyl_array_cfetch(f, pidx),
      (double *)gkyl_array_fetch(fbar, pidx)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

int null = up->pkpm_normalize(
  (const double *)gkyl_array_cfetch(n, cidx), (const double *)gkyl_array_cfetch(f, pidx),
  (double *)gkyl_array_fetch(fbar, pidx)
);
if (null > 0) {
  atomicAdd(num_null, null);
}
}
}
__global__ void
gkyl_dg_calc_pkpm_vars_transverse_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *bvar, struct gkyl_array *f
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);
    up->pkpm_transverse(
      (const double *)gkyl_array_cfetch(bvar, cidx), (double *)gkyl_array_fetch(f, pidx)
    );
  }
}

// Host-side wrapper for the transversality projection of the first-harmonic blocks
void
gkyl_dg_calc_pkpm_vars_transverse_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, struct gkyl_array *f
)
{
  gkyl_dg_calc_pkpm_vars_transverse_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, bvar->on_dev, f->on_dev
  );
}

// Host-side wrapper for the weak division of the distributions by n
int
gkyl_dg_calc_pkpm_vars_normalize_dist_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *n, const struct gkyl_array *f,
  struct gkyl_array *fbar
)
{
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_normalize_dist_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, n->on_dev, f->on_dev, fbar->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_diag_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    up->pkpm_diag(
      (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, loc),
      (const double *)gkyl_array_cfetch(euler_pkpm, loc),
      (const double *)gkyl_array_cfetch(bvar, loc), (const double *)gkyl_array_cfetch(grad_b, loc),
      (const double *)gkyl_array_cfetch(u, loc), (const double *)gkyl_array_cfetch(grad_u, loc),
      (const double *)gkyl_array_cfetch(theta, loc),
      (const double *)gkyl_array_cfetch(alpha_0, loc), (const double *)gkyl_array_cfetch(xi, loc),
      (const double *)gkyl_array_cfetch(c_h, loc), (double *)gkyl_array_fetch(fluid_io, loc),
      (double *)gkyl_array_fetch(pkpm_vars_io, loc)
    );
  }
}

void
gkyl_dg_calc_pkpm_vars_diag_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
)
{
  gkyl_dg_calc_pkpm_vars_diag_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, bvar->on_dev,
    grad_b->on_dev, u->on_dev, grad_u->on_dev, theta->on_dev, alpha_0->on_dev, xi->on_dev,
    c_h->on_dev, fluid_io->on_dev, pkpm_vars_io->on_dev
  );
}

// One thread per interior cell: the gradient of the field in every direction from the
// cell and its two neighbors.
__global__ void
gkyl_dg_calc_pkpm_vars_in_cell_grad_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *blend,
  const struct gkyl_array *q, struct gkyl_array *grad
)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    const double *q_c = (const double *)gkyl_array_cfetch(q, loc);
    const double *blend_c = (const double *)gkyl_array_cfetch(blend, loc);
    double *gd = (double *)gkyl_array_fetch(grad, loc);
    for (int d = 0; d < cdim; ++d) {
      gkyl_copy_int_arr(cdim, idx, idx_l);
      gkyl_copy_int_arr(cdim, idx, idx_r);
      idx_l[d] = idx_l[d] - 1;
      idx_r[d] = idx_r[d] + 1;
      long loc_r = gkyl_range_idx(&conf_range, idx_r);
      const double *blend_r = (const double *)gkyl_array_cfetch(blend, loc_r);
      up->pkpm_in_cell_grad[d](
        up->conf_grid.dx, (const double *)gkyl_array_cfetch(up->jacob_pos, loc), blend_c[d],
        blend_r[d], (const double *)gkyl_array_cfetch(q, gkyl_range_idx(&conf_range, idx_l)), q_c,
        (const double *)gkyl_array_cfetch(q, loc_r), gd
      );
    }
  }
}

void
gkyl_dg_calc_pkpm_vars_in_cell_grad_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *blend, const struct gkyl_array *q, struct gkyl_array *grad
)
{
  gkyl_dg_calc_pkpm_vars_in_cell_grad_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, blend->on_dev, q->on_dev, grad->on_dev
  );
}

__global__ void
gkyl_dg_calc_pkpm_vars_source_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *qmem,
  const struct gkyl_array *fluid, struct gkyl_array *rhs
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    up->pkpm_source(
      (const double *)gkyl_array_cfetch(qmem, loc), (const double *)gkyl_array_cfetch(fluid, loc),
      (double *)gkyl_array_fetch(rhs, loc)
    );
  }
}

// Host-side wrapper for pkpm source term calculations
void
gkyl_dg_calc_pkpm_vars_source_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *fluid, struct gkyl_array *rhs
)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_source_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, *conf_range, qmem->on_dev, fluid->on_dev, rhs->on_dev
  );
}

// CUDA kernel to set device pointers to pkpm vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_calc_pkpm_vars_set_cu_dev_ptrs(
  struct gkyl_dg_calc_pkpm_vars *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  up->pkpm_u_set = choose_pkpm_u_set_kern(b_type, cdim, poly_order);
  up->pkpm_u_copy = choose_pkpm_u_copy_kern(b_type, cdim, poly_order);
  for (int d = 0; d < cdim; ++d) {
    up->pkpm_in_cell_grad[d] = choose_pkpm_in_cell_grad_kern(b_type, cdim, poly_order, d);
  }
  up->pkpm_closure = choose_pkpm_closure_kern(b_type, cdim, poly_order);
  up->pkpm_prim = choose_pkpm_prim_kern(b_type, cdim, poly_order);
  up->pkpm_int = choose_pkpm_int_kern(b_type, cdim, poly_order);
  up->pkpm_grad_ln_n = choose_pkpm_grad_ln_n_kern(b_type, cdim, poly_order);
  up->pkpm_fluid_positivity = choose_pkpm_fluid_positivity_kern(b_type, cdim, poly_order);
  up->pkpm_normalize = choose_pkpm_normalize_kern(b_type, cdim, poly_order);
  up->pkpm_transverse = choose_pkpm_transverse_kern(b_type, cdim, poly_order);
  up->pkpm_source = choose_pkpm_source_kern(b_type, cdim, poly_order);
  up->pkpm_diag = choose_pkpm_diag_kern(b_type, cdim, poly_order);
}

gkyl_dg_calc_pkpm_vars *
gkyl_dg_calc_pkpm_vars_cu_dev_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis,
  const struct gkyl_vlasov_position_map *pos_map
)
{
  struct gkyl_dg_calc_pkpm_vars *up =
    (struct gkyl_dg_calc_pkpm_vars *)gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  assert(pos_map);
  up->pos_map = 0; // host pointer set after the device memcpy below
  up->jacob_pos = pos_map->jacob_pos->on_dev;

  up->num_null_rho = 0;
  up->num_null_rho_dev = (int *)gkyl_cu_malloc(sizeof(int));
  up->closure_cond_min = 1.0;
  up->closure_cond_min_dev = (unsigned long long *)gkyl_cu_malloc(sizeof(unsigned long long));
  up->num_pos_dev = (int *)gkyl_cu_malloc(sizeof(int[4]));

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_vars *up_cu =
    (struct gkyl_dg_calc_pkpm_vars *)gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_vars_set_cu_dev_ptrs<<<1, 1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(pos_map);
  up->jacob_pos = pos_map->jacob_pos;

  return up;
}
