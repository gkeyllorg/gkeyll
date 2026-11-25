/* -*- c++ -*- */

#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
#include <cub/cub.cuh>
#endif

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_calc_sr_vars_priv.h>
#include <gkyl_util.h>
}

// Restrict each indirect call to its generated family; the runtime pointer
// still selects the kernel. Integer matrix offsets retain their 32-bit ABI.
GKYL_CU_D static __forceinline__ void
sr_vars_n_set_call(
  sr_n_set_t kernel, int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0,
  const double *M1i
)
{
  asm volatile("{\n"
               " .param .b32 a0;\n"
               " .param .b64 a1, a2, a3, a4;\n"
               " st.param.b32 [a0], %1;\n"
               " st.param.b64 [a1], %2;\n"
               " st.param.b64 [a2], %3;\n"
               " st.param.b64 [a3], %4;\n"
               " st.param.b64 [a4], %5;\n"
               " targets: .calltargets " GKYL_SR_VARS_N_SET_CALL_TARGETS ";\n"
               " call %0, (a0, a1, a2, a3, a4), targets;\n"
               "}\n"
               :
               : "l"(kernel), "r"(count), "l"(A), "l"(rhs), "l"(M0), "l"(M1i)
               : "memory");
}

GKYL_CU_D static __forceinline__ void
sr_vars_n_copy_call(sr_n_copy_t kernel, int count, struct gkyl_nmat *x, const double *M0, double *n)
{
  asm volatile("{\n"
               " .param .b32 a0;\n"
               " .param .b64 a1, a2, a3;\n"
               " st.param.b32 [a0], %1;\n"
               " st.param.b64 [a1], %2;\n"
               " st.param.b64 [a2], %3;\n"
               " st.param.b64 [a3], %4;\n"
               " targets: .calltargets " GKYL_SR_VARS_N_COPY_CALL_TARGETS ";\n"
               " call %0, (a0, a1, a2, a3), targets;\n"
               "}\n"
               :
               : "l"(kernel), "r"(count), "l"(x), "l"(M0), "l"(n)
               : "memory");
}

GKYL_CU_D static __forceinline__ void
sr_vars_GammaV_call(
  sr_GammaV_t kernel, const double *u_i, double *u_i_sq, double *GammaV, double *GammaV_sq
)
{
  asm volatile("{\n"
               " .param .b64 a0, a1, a2, a3;\n"
               " st.param.b64 [a0], %1;\n"
               " st.param.b64 [a1], %2;\n"
               " st.param.b64 [a2], %3;\n"
               " st.param.b64 [a3], %4;\n"
               " targets: .calltargets " GKYL_SR_VARS_GAMMAV_CALL_TARGETS ";\n"
               " call %0, (a0, a1, a2, a3), targets;\n"
               "}\n"
               :
               : "l"(kernel), "l"(u_i), "l"(u_i_sq), "l"(GammaV), "l"(GammaV_sq)
               : "memory");
}

GKYL_CU_D static __forceinline__ void
sr_vars_pressure_call(
  sr_pressure_t kernel, const double *w, const double *dxv, const double *vmap, const double *gamma,
  const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV,
  const double *GammaV_sq, const double *f, double *pressure
)
{
  asm volatile("{\n"
               " .param .b64 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;\n"
               " st.param.b64 [a0], %1;\n"
               " st.param.b64 [a1], %2;\n"
               " st.param.b64 [a2], %3;\n"
               " st.param.b64 [a3], %4;\n"
               " st.param.b64 [a4], %5;\n"
               " st.param.b64 [a5], %6;\n"
               " st.param.b64 [a6], %7;\n"
               " st.param.b64 [a7], %8;\n"
               " st.param.b64 [a8], %9;\n"
               " st.param.b64 [a9], %10;\n"
               " st.param.b64 [a10], %11;\n"
               " targets: .calltargets " GKYL_SR_VARS_PRESSURE_CALL_TARGETS ";\n"
               " call %0, (a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10), targets;\n"
               "}\n"
               :
               : "l"(kernel), "l"(w), "l"(dxv), "l"(vmap), "l"(gamma), "l"(gamma_inv), "l"(u_i),
                 "l"(u_i_sq), "l"(GammaV), "l"(GammaV_sq), "l"(f), "l"(pressure)
               : "memory");
}

__global__ static void
gkyl_dg_calc_sr_vars_n_set_cu_kernel(
  gkyl_dg_calc_sr_vars *up, struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_range conf_range, const struct gkyl_array *M0, const struct gkyl_array *M1i
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
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1 * up->Ncomp;

    const double *M0_d = (const double *)gkyl_array_cfetch(M0, loc);
    const double *M1i_d = (const double *)gkyl_array_cfetch(M1i, loc);

    sr_vars_n_set_call(up->sr_n_set, count, As, xs, M0_d, M1i_d);
  }
}

__global__ static void
gkyl_dg_calc_sr_vars_n_copy_cu_kernel(
  gkyl_dg_calc_sr_vars *up, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array *M0, struct gkyl_array *n
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
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1 * up->Ncomp;

    const double *M0_d = (const double *)gkyl_array_cfetch(M0, loc);
    double *n_d = (double *)gkyl_array_fetch(n, loc);

    sr_vars_n_copy_call(up->sr_n_copy, count, xs, M0_d, n_d);
  }
}

// Host-side wrapper for SR rest-frame density calculation
void
gkyl_dg_calc_sr_vars_n_cu(
  struct gkyl_dg_calc_sr_vars *up, const struct gkyl_array *M0, const struct gkyl_array *M1i,
  struct gkyl_array *n
)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_dg_calc_sr_vars_n_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(
    up->on_dev, up->As->on_dev, up->xs->on_dev, conf_range, M0->on_dev, M1i->on_dev
  );

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_dg_calc_sr_vars_n_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(
    up->on_dev, up->xs->on_dev, conf_range, M0->on_dev, n->on_dev
  );
}

__global__ void
gkyl_dg_calc_sr_vars_GammaV_cu_kernel(
  struct gkyl_dg_calc_sr_vars *up, struct gkyl_range conf_range, const struct gkyl_array *u_i,
  struct gkyl_array *u_i_sq, struct gkyl_array *GammaV, struct gkyl_array *GammaV_sq
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
    long loc_conf = gkyl_range_idx(&conf_range, idx);

    const double *u_i_d = (const double *)gkyl_array_cfetch(u_i, loc_conf);
    double *u_i_sq_d = (double *)gkyl_array_fetch(u_i_sq, loc_conf);
    double *GammaV_d = (double *)gkyl_array_fetch(GammaV, loc_conf);
    double *GammaV_sq_d = (double *)gkyl_array_fetch(GammaV_sq, loc_conf);

    sr_vars_GammaV_call(up->sr_GammaV, u_i_d, u_i_sq_d, GammaV_d, GammaV_sq_d);
  }
}

// Host-side wrapper for derived quantities from bulk four-velocity
void
gkyl_dg_calc_sr_vars_GammaV_cu(
  struct gkyl_dg_calc_sr_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *u_i, struct gkyl_array *u_i_sq, struct gkyl_array *GammaV,
  struct gkyl_array *GammaV_sq
)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_sr_vars_GammaV_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, *conf_range, u_i->on_dev, u_i_sq->on_dev, GammaV->on_dev, GammaV_sq->on_dev
  );
}

__device__ static void
dg_calc_sr_vars_pressure_cell(
  struct gkyl_dg_calc_sr_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  bool use_vmap, const struct gkyl_array *vmap, const struct gkyl_array *gamma,
  const struct gkyl_array *gamma_inv, const struct gkyl_array *u_i, const struct gkyl_array *u_i_sq,
  const struct gkyl_array *GammaV, const struct gkyl_array *GammaV_sq, const struct gkyl_array *f,
  const int *idx, double *mom_local
)
{
  int cdim = conf_range.ndim;
  int pdim = phase_range.ndim;

  double xc[GKYL_MAX_DIM];
  int idx_vel[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
  long loc_conf = gkyl_range_idx(&conf_range, idx);
  long loc_phase = gkyl_range_idx(&phase_range, idx);
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long loc_vel = gkyl_range_idx(&up->vel_range, idx_vel);

  sr_vars_pressure_call(
    up->sr_pressure, xc, up->phase_grid.dx,
    use_vmap ? (const double *)gkyl_array_cfetch(vmap, loc_vel) : 0,
    (const double *)gkyl_array_cfetch(gamma, loc_vel),
    (const double *)gkyl_array_cfetch(gamma_inv, loc_vel),
    (const double *)gkyl_array_cfetch(u_i, loc_conf),
    (const double *)gkyl_array_cfetch(u_i_sq, loc_conf),
    (const double *)gkyl_array_cfetch(GammaV, loc_conf),
    (const double *)gkyl_array_cfetch(GammaV_sq, loc_conf),
    (const double *)gkyl_array_cfetch(f, loc_phase), mom_local
  );
}

#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
__global__ static void
dg_calc_sr_vars_pressure_cu_kernel(
  struct gkyl_dg_calc_sr_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  bool use_vmap, const struct gkyl_array *vmap, const struct gkyl_array *gamma,
  const struct gkyl_array *gamma_inv, const struct gkyl_array *u_i, const struct gkyl_array *u_i_sq,
  const struct gkyl_array *GammaV, const struct gkyl_array *GammaV_sq, const struct gkyl_array *f,
  struct gkyl_array *sr_pressure
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    double mom_local[32] = {0.0};
    dg_calc_sr_vars_pressure_cell(
      up, conf_range, phase_range, use_vmap, vmap, gamma, gamma_inv, u_i, u_i_sq, GammaV, GammaV_sq,
      f, idx, mom_local
    );
    double *pressure = (double *)gkyl_array_fetch(sr_pressure, gkyl_range_idx(&conf_range, idx));
    for (unsigned int k = 0; k < sr_pressure->ncomp; ++k) {
      atomicAdd(&pressure[k], mom_local[k]);
    }
  }
}
#else
enum { SR_PRESSURE_TILE_SIZE = 128 };

__global__ static void
dg_calc_sr_vars_pressure_partial_cu_kernel(
  struct gkyl_dg_calc_sr_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  bool use_vmap, const struct gkyl_array *vmap, const struct gkyl_array *gamma,
  const struct gkyl_array *gamma_inv, const struct gkyl_array *u_i, const struct gkyl_array *u_i_sq,
  const struct gkyl_array *GammaV, const struct gkyl_array *GammaV_sq, const struct gkyl_array *f,
  double *partial_dev, long conf_start, long num_tiles, long num_vel, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, SR_PRESSURE_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long conf_local = conf_start + blockIdx.x / num_tiles;
  long tile = blockIdx.x % num_tiles;
  long vel_local = tile * SR_PRESSURE_TILE_SIZE + threadIdx.x;
  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&conf_range, conf_local, idx);

  double mom_local[32] = {0.0};
  if (vel_local < num_vel) {
    long remaining = vel_local;
    for (int d = conf_range.ndim; d < phase_range.ndim; ++d) {
      int extent = phase_range.upper[d] - phase_range.lower[d] + 1;
      idx[d] = phase_range.lower[d] + remaining % extent;
      remaining /= extent;
    }
    dg_calc_sr_vars_pressure_cell(
      up, conf_range, phase_range, use_vmap, vmap, gamma, gamma_inv, u_i, u_i_sq, GammaV, GammaV_sq,
      f, idx, mom_local
    );
  }

  for (unsigned int k = 0; k < ncomp; ++k) {
    double sum = BlockReduceT(temp).Sum(mom_local[k]);
    if (threadIdx.x == 0) {
      long segment = (blockIdx.x / num_tiles) * ncomp + k;
      partial_dev[segment * num_tiles + tile] = sum;
    }
    __syncthreads();
  }
}

__global__ static void
dg_calc_sr_vars_pressure_finish_cu_kernel(
  struct gkyl_range conf_range, struct gkyl_array *sr_pressure, const double *partial_dev,
  long conf_start, long num_tiles, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, SR_PRESSURE_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long segment = blockIdx.x;
  double local_sum = 0.0;
  for (long tile = threadIdx.x; tile < num_tiles; tile += blockDim.x) {
    local_sum += partial_dev[segment * num_tiles + tile];
  }
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    long conf_local = conf_start + segment / ncomp;
    int idx[GKYL_MAX_DIM];
    gkyl_sub_range_inv_idx(&conf_range, conf_local, idx);
    double *pressure = (double *)gkyl_array_fetch(sr_pressure, gkyl_range_idx(&conf_range, idx));
    pressure[segment % ncomp] = sum;
  }
}
#endif
// Host-side wrapper for rest-frame pressure calculation
void
gkyl_dg_calc_sr_vars_pressure_cu(
  struct gkyl_dg_calc_sr_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *gamma,
  const struct gkyl_array *gamma_inv, const struct gkyl_array *u_i, const struct gkyl_array *u_i_sq,
  const struct gkyl_array *GammaV, const struct gkyl_array *GammaV_sq, const struct gkyl_array *f,
  struct gkyl_array *sr_pressure
)
{
  gkyl_array_clear(sr_pressure, 0.0);
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  long num_vel = 1;
  for (int d = conf_range->ndim; d < phase_range->ndim; ++d) {
    num_vel *= phase_range->upper[d] - phase_range->lower[d] + 1;
  }
  if (conf_range->volume == 0 || num_vel == 0) {
    return;
  }
  long num_tiles = (num_vel + SR_PRESSURE_TILE_SIZE - 1) / SR_PRESSURE_TILE_SIZE;
  unsigned int ncomp = sr_pressure->ncomp;
  const long max_partial_values = (64L * 1024 * 1024) / sizeof(double);
  long batch_conf = max_partial_values / (num_tiles * ncomp);
  if (batch_conf < 1) {
    batch_conf = 1;
  }
  if (batch_conf > conf_range->volume) {
    batch_conf = conf_range->volume;
  }
  size_t required_bytes = batch_conf * num_tiles * ncomp * sizeof(double);
  if (required_bytes > up->pressure_buffer_bytes) {
    if (up->pressure_buffer_dev) {
      gkyl_cu_free(up->pressure_buffer_dev);
    }
    up->pressure_buffer_dev = (double *)gkyl_cu_malloc(required_bytes);
    up->pressure_buffer_bytes = required_bytes;
  }
  for (long start = 0; start < conf_range->volume; start += batch_conf) {
    long count = conf_range->volume - start < batch_conf ? conf_range->volume - start : batch_conf;
    long segments = count * ncomp;
    dg_calc_sr_vars_pressure_partial_cu_kernel<<<count * num_tiles, SR_PRESSURE_TILE_SIZE>>>(
      up->on_dev, *conf_range, *phase_range, up->use_vmap, up->use_vmap ? up->vmap->on_dev : 0,
      gamma->on_dev, gamma_inv->on_dev, u_i->on_dev, u_i_sq->on_dev, GammaV->on_dev,
      GammaV_sq->on_dev, f->on_dev, up->pressure_buffer_dev, start, num_tiles, num_vel, ncomp
    );
    dg_calc_sr_vars_pressure_finish_cu_kernel<<<segments, SR_PRESSURE_TILE_SIZE>>>(
      *conf_range, sr_pressure->on_dev, up->pressure_buffer_dev, start, num_tiles, ncomp
    );
  }
#else
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  dg_calc_sr_vars_pressure_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, *conf_range, *phase_range, up->use_vmap, up->use_vmap ? up->vmap->on_dev : 0,
    gamma->on_dev, gamma_inv->on_dev, u_i->on_dev, u_i_sq->on_dev, GammaV->on_dev,
    GammaV_sq->on_dev, f->on_dev, sr_pressure->on_dev
  );
#endif
}

// CUDA kernel to set device pointers to sr vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_calc_sr_vars_set_cu_dev_ptrs(
  struct gkyl_dg_calc_sr_vars *up, enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order
)
{
  up->sr_pressure = choose_sr_vars_pressure_kern(b_type, cdim, vdim, poly_order);
  up->sr_n_set = choose_sr_vars_n_set_kern(b_type, cdim, vdim, poly_order);
  up->sr_n_copy = choose_sr_vars_n_copy_kern(b_type, cdim, vdim, poly_order);
  up->sr_GammaV = choose_sr_vars_GammaV_kern(b_type, cdim, vdim, poly_order);
  up->sr_pressure = choose_sr_vars_pressure_kern(b_type, cdim, vdim, poly_order);
}

gkyl_dg_calc_sr_vars *
gkyl_dg_calc_sr_vars_cu_dev_new(
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_rect_grid *vel_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *vel_basis,
  const struct gkyl_range *mem_range, const struct gkyl_range *vel_range,
  const struct gkyl_vlasov_velocity_map *vel_map
)
{
  struct gkyl_dg_calc_sr_vars *up = (struct gkyl_dg_calc_sr_vars *)gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  up->vel_grid = *vel_grid;
  up->vel_range = *vel_range;

  int nc = conf_basis->num_basis;
  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  enum gkyl_basis_type b_type = conf_basis->b_type;
  // store polynomial order and mem_range for linear solve
  up->poly_order = poly_order;
  up->mem_range = *mem_range;

  int vdim = vel_basis->ndim;
  up->use_vmap = vel_map && vel_map->is_mapped;
  up->vmap = 0;
  if (up->use_vmap) {
    up->vmap = gkyl_array_acquire(vel_map->vmap);
  }

  // Linear system for solving for the drift velocity V_drift = M1i/M0
  // and then computing the rest-frame density n = GammaV_inv*M0
  // where GammaV_inv = sqrt(1 - |V_drift|^2)
  up->Ncomp = vdim;
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp * mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp * mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->pressure_buffer_dev = 0;
  up->pressure_buffer_bytes = 0;

  struct gkyl_dg_calc_sr_vars *up_cu =
    (struct gkyl_dg_calc_sr_vars *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_sr_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_sr_vars_set_cu_dev_ptrs<<<1, 1>>>(up_cu, b_type, cdim, vdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;
}
