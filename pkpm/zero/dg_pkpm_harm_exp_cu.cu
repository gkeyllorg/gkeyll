/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_exp_priv.h>
}

// One thread per configuration cell.
__global__ static void
gkyl_dg_pkpm_harm_exp_coef_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, struct gkyl_range conf_range, const struct gkyl_array *em,
  struct gkyl_array *coef
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    pkpm_hx_coef_cell(up, (const double *)gkyl_array_cfetch(em, loc), (double *)gkyl_array_fetch(coef, loc));
  }
}

void
gkyl_dg_pkpm_harm_exp_coef_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *coef
)
{
  gkyl_dg_pkpm_harm_exp_coef_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, em->on_dev, coef->on_dev
  );
}

// One thread per phase cell.
__global__ static void
gkyl_dg_pkpm_harm_exp_advance_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, struct gkyl_range conf_range,
  struct gkyl_range phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  struct gkyl_array *nu, struct gkyl_array *act, struct gkyl_array *out
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < phase_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc, idx);
    long cloc = gkyl_range_idx(&conf_range, idx);
    long ploc = gkyl_range_idx(&phase_range, idx);
    // act is shared by the cells of a column; every writer writes 1
    pkpm_hx_advance_cell(up, stage, h, (const double *)gkyl_array_cfetch(coef, cloc),
      (const double *)gkyl_array_cfetch(f0, ploc), 0, (double *)gkyl_array_fetch(nu, ploc),
      (double *)gkyl_array_fetch(act, cloc), (double *)gkyl_array_fetch(out, ploc));
  }
}

void
gkyl_dg_pkpm_harm_exp_advance_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  struct gkyl_array *nu, struct gkyl_array *act, struct gkyl_array *out
)
{
  gkyl_dg_pkpm_harm_exp_advance_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, stage, h, *conf_range, *phase_range, coef->on_dev, f0->on_dev, nu->on_dev,
    act->on_dev, out->on_dev
  );
}

// One thread per configuration cell: the w integrals of a flagged column.
__global__ static void
gkyl_dg_pkpm_harm_exp_col_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *act, struct gkyl_array *col, const struct gkyl_array *out
)
{
  int cidx[GKYL_MAX_DIM], pidx[GKYL_MAX_DIM];
  int cdim = up->cdim, nc = 4 * up->num_conf_nodes;
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, cidx);
    long cloc = gkyl_range_idx(&conf_range, cidx);
    if (*(const double *)gkyl_array_cfetch(act, cloc) == 0.0) {
      continue;
    }
    double *c = (double *)gkyl_array_fetch(col, cloc);
    for (int i = 0; i < nc; ++i) {
      c[i] = 0.0;
    }
    for (int d = 0; d < cdim; ++d) {
      pidx[d] = cidx[d];
    }
    for (int iv = phase_range.lower[cdim]; iv <= phase_range.upper[cdim]; ++iv) {
      pidx[cdim] = iv;
      pkpm_hx_col_add(up, (const double *)gkyl_array_cfetch(out, gkyl_range_idx(&phase_range, pidx)), c);
    }
  }
}

// One thread per phase cell: the frame shift of the flagged columns.
__global__ static void
gkyl_dg_pkpm_harm_exp_shift_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *act, const struct gkyl_array *col, struct gkyl_array *out
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < phase_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc, idx);
    long cloc = gkyl_range_idx(&conf_range, idx);
    if (*(const double *)gkyl_array_cfetch(act, cloc) == 0.0) {
      continue;
    }
    pkpm_hx_restore_cell(up, (const double *)gkyl_array_cfetch(col, cloc),
      (double *)gkyl_array_fetch(out, gkyl_range_idx(&phase_range, idx)));
  }
}

void
gkyl_dg_pkpm_harm_exp_restore_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *act, struct gkyl_array *col,
  struct gkyl_array *out
)
{
  gkyl_dg_pkpm_harm_exp_col_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, act->on_dev, col->on_dev, out->on_dev
  );
  gkyl_dg_pkpm_harm_exp_shift_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, act->on_dev, col->on_dev, out->on_dev
  );
}

struct gkyl_dg_pkpm_harm_exp *
gkyl_dg_pkpm_harm_exp_cu_dev_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
  struct gkyl_dg_pkpm_harm_exp *up =
    (struct gkyl_dg_pkpm_harm_exp *)gkyl_malloc(sizeof(struct gkyl_dg_pkpm_harm_exp));
  up->cdim = inp->conf_basis->ndim;
  up->num_conf_nodes = 1 << up->cdim;
  up->num_phase_nodes = 3 * up->num_conf_nodes;
  up->qbym = inp->qbym;
  up->mix = inp->mix > 0.0 ? inp->mix : 0.0;
  up->half_dw = 0.5 * inp->phase_grid->dx[up->cdim];
  up->use_gpu = true;
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  struct gkyl_dg_pkpm_harm_exp *up_cu =
    (struct gkyl_dg_pkpm_harm_exp *)gkyl_cu_malloc(sizeof(struct gkyl_dg_pkpm_harm_exp));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_dg_pkpm_harm_exp), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;
  return up;
}
