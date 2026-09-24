/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_fluid_vars.h>
#include <gkyl_dg_calc_fluid_vars_priv.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_fluid_vars_set_cu_kernel(gkyl_dg_calc_fluid_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* fluid, struct gkyl_array* cell_avg_prim)
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

    const double *fluid_d = (const double*) gkyl_array_cfetch(fluid, loc);

    int* cell_avg_prim_d = (int*) gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->fluid_set(count, As, xs, fluid_d);
  }
}

__global__ static void
gkyl_dg_calc_fluid_vars_copy_cu_kernel(gkyl_dg_calc_fluid_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array* prim, struct gkyl_array* prim_surf)
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

    double* prim_d = (double*) gkyl_array_fetch(prim, loc);
    double* prim_surf_d = (double*) gkyl_array_fetch(prim_surf, loc);

    up->fluid_copy(count, xs, prim_d, prim_surf_d);
  }
}

// Host-side wrapper for primitive variable calculation
void gkyl_dg_calc_fluid_vars_advance_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_array* fluid, struct gkyl_array* cell_avg_prim, 
  struct gkyl_array* prim, struct gkyl_array* prim_surf)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_fluid_vars_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    fluid->on_dev, cell_avg_prim->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_dg_calc_fluid_vars_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, prim->on_dev, prim_surf->on_dev);
}

__global__ void
gkyl_calc_fluid_vars_pressure_cu_kernel(struct gkyl_dg_calc_fluid_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf)
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

    const double *fluid_d = (const double*) gkyl_array_cfetch(fluid, loc);
    const double *u_d = (const double*) gkyl_array_cfetch(u, loc);

    double *p_d = (double*) gkyl_array_fetch(p, loc);
    double *p_surf_d = (double*) gkyl_array_fetch(p_surf, loc);
    up->fluid_pressure(up->param, fluid_d, u_d, p_d, p_surf_d);
  }
}

// Host-side wrapper for pressure calculation
void gkyl_dg_calc_fluid_vars_pressure_cu(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_calc_fluid_vars_pressure_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    fluid->on_dev, u->on_dev, 
    p->on_dev, p_surf->on_dev);
}

__global__ void
gkyl_calc_fluid_vars_ke_cu_kernel(struct gkyl_dg_calc_fluid_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* ke)
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

    const double *fluid_d = (const double*) gkyl_array_cfetch(fluid, loc);
    const double *u_d = (const double*) gkyl_array_cfetch(u, loc);

    double *ke_d = (double*) gkyl_array_fetch(ke, loc);
    up->fluid_ke(fluid_d, u_d, ke_d);
  }
}

// Host-side wrapper for kinetic energy calculation
void gkyl_dg_calc_fluid_vars_ke_cu(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* ke)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_calc_fluid_vars_ke_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    fluid->on_dev, u->on_dev, 
    ke->on_dev);
}

__global__ void
gkyl_dg_calc_fluid_vars_limiter_cu_kernel(struct gkyl_dg_calc_fluid_vars *up, struct gkyl_range conf_range, 
  struct gkyl_array* fluid)
{
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idxc);
    const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(up->geom, idxc);    

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&conf_range, idxc);

    double *fluid_c = (double*) gkyl_array_fetch(fluid, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      double *fluid_l = (double*) gkyl_array_fetch(fluid, linl);
      double *fluid_r = (double*) gkyl_array_fetch(fluid, linr);
      
      up->fluid_limiter[dir](up->limiter_fac, up->wv_eqn, geom, fluid_l, fluid_c, fluid_r);
    }
  }
}

// Host-side wrapper for slope limiter of fluid variables
void
gkyl_dg_calc_fluid_vars_limiter_cu(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_range *conf_range, 
  struct gkyl_array* fluid)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_fluid_vars_limiter_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, fluid->on_dev);
}

__global__ void
gkyl_dg_calc_fluid_integrated_vars_cu_kernel(struct gkyl_dg_calc_fluid_vars *up, 
  struct gkyl_range conf_range, const struct gkyl_array* fluid, 
  const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* int_fluid_vars)
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

    const double *fluid_d = (const double*) gkyl_array_cfetch(fluid, loc);
    const double *u_i_d = (const double*) gkyl_array_cfetch(u_i, loc);
    const double *p_ij_d = (const double*) gkyl_array_cfetch(p_ij, loc);

    double *int_fluid_vars_d = (double*) gkyl_array_fetch(int_fluid_vars, loc);
    up->fluid_int(fluid_d, u_i_d, p_ij_d, int_fluid_vars_d);
  }
}

// Host-side wrapper for fluid integrated variables calculation
void
gkyl_dg_calc_fluid_integrated_vars_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* fluid, 
  const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* int_fluid_vars)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_fluid_integrated_vars_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    fluid->on_dev, u_i->on_dev, p_ij->on_dev, int_fluid_vars->on_dev);
}

__global__ void
gkyl_dg_calc_fluid_vars_source_cu_kernel(struct gkyl_dg_calc_fluid_vars *up, 
  struct gkyl_range conf_range, 
  const struct gkyl_array* app_accel, const struct gkyl_array* fluid, 
  struct gkyl_array* rhs)
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

    const double *app_accel_d = (const double*) gkyl_array_cfetch(app_accel, loc);
    const double *fluid_d = (const double*) gkyl_array_cfetch(fluid, loc);

    double *rhs_d = (double*) gkyl_array_fetch(rhs, loc);
    up->fluid_source(app_accel_d, fluid_d, rhs_d);
  }
}

// Host-side wrapper for fluid source term calculations
void
gkyl_dg_calc_fluid_vars_source_cu(struct gkyl_dg_calc_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* app_accel, const struct gkyl_array* fluid, 
  struct gkyl_array* rhs)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_fluid_vars_source_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    app_accel->on_dev, fluid->on_dev, rhs->on_dev);
}

// CUDA kernel to set device pointers to fluid vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_fluid_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_wv_eqn *wv_eqn, 
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  up->fluid_set = choose_fluid_set_kern(b_type, cdim, poly_order);
  up->fluid_copy = choose_fluid_copy_kern(b_type, cdim, poly_order);
  up->fluid_pressure = choose_fluid_pressure_kern(b_type, cdim, poly_order);
  up->fluid_ke = choose_fluid_ke_kern(b_type, cdim, poly_order);
  up->fluid_int = choose_fluid_int_kern(b_type, cdim, poly_order);
  up->fluid_source = choose_fluid_source_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    up->fluid_limiter[d] = choose_fluid_limiter_kern(d, b_type, cdim, poly_order);
}

gkyl_dg_calc_fluid_vars*
gkyl_dg_calc_fluid_vars_cu_dev_new(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *wg, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  double limiter_fac)
{
  struct gkyl_dg_calc_fluid_vars *up = (struct gkyl_dg_calc_fluid_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_fluid_vars));

  up->eqn_type = wv_eqn->type;
  if (up->eqn_type == GKYL_EQN_EULER)
    up->param = gkyl_wv_euler_gas_gamma(wv_eqn);

  // acquire pointer to wave equation object
  struct gkyl_wv_eqn *eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->wv_eqn = eqn->on_dev; // this is so the memcpy below has eqn on_dev

  // acquire pointer to wave equation object
  struct gkyl_wave_geom *geom = gkyl_wave_geom_acquire(wg);
  up->geom = geom->on_dev; // this is so the memcpy below has geom on_dev

  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->Ncomp = 3;
  up->mem_range = *mem_range;

  // Limiter factor for relationship between slopes and cell average differences
  // By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
  // and a cell slope estimate from two adjacent cells is (for the x variation): 
  // integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
  // where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
  // This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
  if (limiter_fac == 0.0)
    up->limiter_fac = 0.5773502691896258;
  else
    up->limiter_fac = limiter_fac;

  // There are Ncomp*range->volume linear systems to be solved 
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_fluid_vars *up_cu = (struct gkyl_dg_calc_fluid_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_fluid_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_fluid_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_fluid_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, wv_eqn->on_dev, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  up->wv_eqn = eqn; // updater should store host pointer 
  up->geom = geom; 

  return up;
}
