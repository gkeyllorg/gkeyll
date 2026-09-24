/* -*- c++ -*- */

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

__global__ void
gkyl_calc_sr_vars_init_p_vars_cu_kernel(gkyl_dg_calc_sr_vars* up, 
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv)
{
  int idx[GKYL_MAX_DIM];
  // Cell center array
  double xc[GKYL_MAX_DIM];  

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < up->vel_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&up->vel_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->vel_grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&up->vel_range, idx);

    double *gamma_d = (double*) gkyl_array_fetch(gamma, loc);
    double *gamma_inv_d = (double*) gkyl_array_fetch(gamma_inv, loc);
    up->sr_p_vars(xc, up->vel_grid.dx, gamma_d, gamma_inv_d);
  }
}

// Host-side wrapper for initialization of momentum variables (gamma, gamma_inv) 
void
gkyl_calc_sr_vars_init_p_vars_cu(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv)
{
  int nblocks = up->vel_range.nblocks;
  int nthreads = up->vel_range.nthreads;
  gkyl_calc_sr_vars_init_p_vars_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    gamma->on_dev, gamma_inv->on_dev);
}

__global__ static void
gkyl_dg_calc_sr_vars_n_set_cu_kernel(gkyl_dg_calc_sr_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* M0, const struct gkyl_array* M1i)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp;

    const double *M0_d = (const double*) gkyl_array_cfetch(M0, loc);
    const double *M1i_d = (const double*) gkyl_array_cfetch(M1i, loc);

    up->sr_n_set(count, As, xs, M0_d, M1i_d);
  }
}

__global__ static void
gkyl_dg_calc_sr_vars_n_copy_cu_kernel(gkyl_dg_calc_sr_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* M0, struct gkyl_array* n)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp;

    const double *M0_d = (const double*) gkyl_array_cfetch(M0, loc);
    double* n_d = (double*) gkyl_array_fetch(n, loc);

    up->sr_n_copy(count, xs, M0_d, n_d);
  }
}

// Host-side wrapper for SR rest-frame density calculation
void gkyl_dg_calc_sr_vars_n_cu(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_array* M0, const struct gkyl_array* M1i, struct gkyl_array* n)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_sr_vars_n_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    M0->on_dev, M1i->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_dg_calc_sr_vars_n_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, M0->on_dev, n->on_dev);
}

__global__ void
gkyl_dg_calc_sr_vars_GammaV_cu_kernel(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_range conf_range,
  const struct gkyl_array* u_i, struct gkyl_array* u_i_sq, 
  struct gkyl_array* GammaV, struct gkyl_array* GammaV_sq)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);

    const double *u_i_d = (const double*) gkyl_array_cfetch(u_i, loc_conf);
    double *u_i_sq_d = (double*) gkyl_array_fetch(u_i_sq, loc_conf);
    double *GammaV_d = (double*) gkyl_array_fetch(GammaV, loc_conf);
    double *GammaV_sq_d = (double*) gkyl_array_fetch(GammaV_sq, loc_conf);

    up->sr_GammaV(u_i_d, u_i_sq_d, GammaV_d, GammaV_sq_d); 
  }
}

// Host-side wrapper for derived quantities from bulk four-velocity
void
gkyl_dg_calc_sr_vars_GammaV_cu(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range,
  const struct gkyl_array* u_i, struct gkyl_array* u_i_sq, 
  struct gkyl_array* GammaV, struct gkyl_array* GammaV_sq)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_sr_vars_GammaV_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    u_i->on_dev, u_i_sq->on_dev, GammaV->on_dev, GammaV_sq->on_dev);
}

__global__ void
gkyl_dg_calc_sr_vars_pressure_cu_kernel(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range, 
  const struct gkyl_array* gamma, const struct gkyl_array* gamma_inv, 
  const struct gkyl_array* u_i, const struct gkyl_array* u_i_sq, 
  const struct gkyl_array* GammaV, const struct gkyl_array* GammaV_sq, 
  const struct gkyl_array* f, struct gkyl_array* sr_pressure)
{
  int cdim = conf_range.ndim;
  int pdim = phase_range.ndim;

  double xc[GKYL_MAX_DIM];
  int idx_vel[GKYL_MAX_DIM];
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    for (int i=0; i<pdim-cdim; ++i) {
      idx_vel[i] = idx[cdim+i];
    }
    long loc_vel = gkyl_range_idx(&up->vel_range, idx_vel);

    const double *gamma_d = (const double*) gkyl_array_cfetch(gamma, loc_vel);
    const double *gamma_inv_d = (const double*) gkyl_array_cfetch(gamma_inv, loc_vel);
    const double *u_i_d = (const double*) gkyl_array_cfetch(u_i, loc_conf);
    const double *u_i_sq_d = (const double*) gkyl_array_cfetch(u_i_sq, loc_conf);
    const double *GammaV_d = (const double*) gkyl_array_cfetch(GammaV, loc_conf);
    const double *GammaV_sq_d = (const double*) gkyl_array_cfetch(GammaV_sq, loc_conf);
    const double *f_d = (const double*) gkyl_array_cfetch(f, loc_phase);

    double momLocal[32]; // hard-coded to max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k=0; k<32; ++k) {
      momLocal[k] = 0.0;
    }

    up->sr_pressure(xc, up->phase_grid.dx, 
      gamma_d, gamma_inv_d, u_i_d, u_i_sq_d, GammaV_d, GammaV_sq_d, 
      f_d, &momLocal[0]);  

    // Accumulate output to output array atomically to avoid race conditions
    double *sr_pressure_d = (double*) gkyl_array_fetch(sr_pressure, loc_conf);
    for (unsigned int k = 0; k < sr_pressure->ncomp; ++k) {
      atomicAdd(&sr_pressure_d[k], momLocal[k]);
    }   
  }  
}
// Host-side wrapper for rest-frame pressure calculation
void 
gkyl_dg_calc_sr_vars_pressure_cu(struct gkyl_dg_calc_sr_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* gamma, const struct gkyl_array* gamma_inv, 
  const struct gkyl_array* u_i, const struct gkyl_array* u_i_sq, 
  const struct gkyl_array* GammaV, const struct gkyl_array* GammaV_sq, 
  const struct gkyl_array* f, struct gkyl_array* sr_pressure)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_array_clear(sr_pressure, 0.0); 
  gkyl_dg_calc_sr_vars_pressure_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, 
    gamma->on_dev, gamma_inv->on_dev, 
    u_i->on_dev, u_i_sq->on_dev, GammaV->on_dev, GammaV_sq->on_dev, 
    f->on_dev, sr_pressure->on_dev);
}

// CUDA kernel to set device pointers to sr vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_sr_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_sr_vars *up, 
  enum gkyl_basis_type b_type, enum gkyl_basis_type b_type_v,
  int cdim, int vdim, int poly_order, int poly_order_v)
{
  up->sr_p_vars = choose_sr_p_vars_kern(b_type_v, vdim, poly_order_v);
  up->sr_n_set = choose_sr_vars_n_set_kern(b_type, cdim, vdim, poly_order);
  up->sr_n_copy = choose_sr_vars_n_copy_kern(b_type, cdim, vdim, poly_order);
  up->sr_GammaV = choose_sr_vars_GammaV_kern(b_type, cdim, vdim, poly_order);
  up->sr_pressure = choose_sr_vars_pressure_kern(b_type, cdim, vdim, poly_order);
}

gkyl_dg_calc_sr_vars*
gkyl_dg_calc_sr_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, const struct gkyl_rect_grid *vel_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *vel_basis, 
  const struct gkyl_range *mem_range, const struct gkyl_range *vel_range)
{
  struct gkyl_dg_calc_sr_vars *up = (struct gkyl_dg_calc_sr_vars*) gkyl_malloc(sizeof(*up));

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
  int poly_order_v = vel_basis->poly_order;
  enum gkyl_basis_type b_type_v = vel_basis->b_type;

  // Linear system for solving for the drift velocity V_drift = M1i/M0 
  // and then computing the rest-frame density n = GammaV_inv*M0 
  // where GammaV_inv = sqrt(1 - |V_drift|^2)
  up->Ncomp = vdim; 
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_sr_vars *up_cu = (struct gkyl_dg_calc_sr_vars*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_sr_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_sr_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, b_type_v, cdim, vdim, poly_order, poly_order_v);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
