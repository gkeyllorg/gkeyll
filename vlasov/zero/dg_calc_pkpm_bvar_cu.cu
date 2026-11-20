/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_bvar_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_calc_pkpm_bvar_advance_cu_kernel(
  struct gkyl_dg_calc_pkpm_bvar *up, struct gkyl_range conf_range, const struct gkyl_array *em,
  struct gkyl_array *bvar, int *num_null_b
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

    int null_b = up->pkpm_bvar(
      (const double *)gkyl_array_cfetch(em, loc), (double *)gkyl_array_fetch(bvar, loc)
    );
    if (null_b > 0) {
      atomicAdd(num_null_b, null_b);
    }
  }
}

// Host-side wrapper for the magnetic field unit vector calculation
void
gkyl_dg_calc_pkpm_bvar_advance_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *bvar
)
{
  gkyl_cu_memset(up->num_null_b_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_bvar_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, em->on_dev, bvar->on_dev, up->num_null_b_dev
  );
  gkyl_cu_memcpy(&up->num_null_b, up->num_null_b_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
}

__global__ void
gkyl_dg_calc_pkpm_bvar_dtB_over_magB_cu_kernel(
  struct gkyl_dg_calc_pkpm_bvar *up, struct gkyl_range conf_range, const struct gkyl_array *em,
  const struct gkyl_array *dt_em, struct gkyl_array *out, int *num_null
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    int null = up->pkpm_dtB_over_magB(
      (const double *)gkyl_array_cfetch(em, loc), (const double *)gkyl_array_cfetch(dt_em, loc),
      (double *)gkyl_array_fetch(out, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

// Host-side wrapper for d_t B/|B| (reuses the null-node counter of b)
int
gkyl_dg_calc_pkpm_bvar_dtB_over_magB_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, const struct gkyl_array *dt_em, struct gkyl_array *dtB_over_magB
)
{
  gkyl_cu_memset(up->num_null_b_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_bvar_dtB_over_magB_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, em->on_dev, dt_em->on_dev, dtB_over_magB->on_dev, up->num_null_b_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_b_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

// One thread per interior cell: the weak gradient of b in every direction.
__global__ void
gkyl_dg_calc_pkpm_bvar_grad_cu_kernel(
  struct gkyl_dg_calc_pkpm_bvar *up, struct gkyl_range conf_range, const struct gkyl_array *bvar,
  struct gkyl_array *grad_b
)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    const double *b_c = (const double *)gkyl_array_cfetch(bvar, loc);
    double *gb = (double *)gkyl_array_fetch(grad_b, loc);
    for (int d = 0; d < cdim; ++d) {
      gkyl_copy_int_arr(cdim, idx, idx_l);
      gkyl_copy_int_arr(cdim, idx, idx_r);
      idx_l[d] = idx_l[d] - 1;
      idx_r[d] = idx_r[d] + 1;
      up->pkpm_grad[d](
        up->dx, (const double *)gkyl_array_cfetch(up->jacob_pos, loc),
        (const double *)gkyl_array_cfetch(bvar, gkyl_range_idx(&conf_range, idx_l)), b_c,
        (const double *)gkyl_array_cfetch(bvar, gkyl_range_idx(&conf_range, idx_r)), gb
      );
    }
  }
}

void
gkyl_dg_calc_pkpm_bvar_grad_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, struct gkyl_array *grad_b
)
{
  gkyl_dg_calc_pkpm_bvar_grad_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, bvar->on_dev, grad_b->on_dev
  );
}

// CUDA kernel to set device function pointers
__global__ static void
dg_calc_pkpm_bvar_set_cu_dev_ptrs(
  struct gkyl_dg_calc_pkpm_bvar *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  up->pkpm_bvar = choose_pkpm_bvar_kern(b_type, cdim, poly_order);
  up->pkpm_dtB_over_magB = choose_pkpm_dtB_over_magB_kern(b_type, cdim, poly_order);
  for (int d = 0; d < cdim; ++d) {
    up->pkpm_grad[d] = choose_pkpm_bvar_grad_kern(b_type, cdim, poly_order, d);
  }
}

gkyl_dg_calc_pkpm_bvar *
gkyl_dg_calc_pkpm_bvar_cu_dev_inew(const struct gkyl_dg_calc_pkpm_bvar_inp *inp)
{
  struct gkyl_dg_calc_pkpm_bvar *up =
    (struct gkyl_dg_calc_pkpm_bvar *)gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_bvar));

  up->cdim = inp->conf_basis->ndim;
  for (int d = 0; d < inp->conf_basis->ndim; ++d) {
    up->dx[d] = inp->conf_grid->dx[d];
  }
  up->num_null_b = 0;
  up->num_null_b_dev = (int *)gkyl_cu_malloc(sizeof(int));
  assert(inp->pos_map);
  up->pos_map = 0; // host pointer set after the device memcpy below
  up->jacob_pos = inp->pos_map->jacob_pos->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_bvar *up_cu =
    (struct gkyl_dg_calc_pkpm_bvar *)gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_bvar));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_bvar), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_bvar_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, inp->conf_basis->b_type, inp->conf_basis->ndim, inp->conf_basis->poly_order
  );

  // set parent on_dev pointer
  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  return up;
}
