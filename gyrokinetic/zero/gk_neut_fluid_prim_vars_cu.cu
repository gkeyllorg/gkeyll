/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_gk_neut_fluid_prim_vars.h>
#include <gkyl_gk_neut_fluid_prim_vars_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_udrift_set_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* moms)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->udrift_ncomp;

    const double *moms_d = (const double*) gkyl_array_cfetch(moms, loc);

    up->udrift_set_prob_ker(count, As, xs, moms_d);
  }
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_udrift_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[3*20]; // udrift_comp=3 for 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->udrift_ncomp;

    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_get_sol_ker(count, xs, prim_vars_d);

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_udrift_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_udrift_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_udrift_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_pressure_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[20]; // 3D p=2.
  double udrift_d[3*20]; // udrift_comp=3 for 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->udrift_ncomp;

    const double *moms_d = (const double*) gkyl_array_cfetch(moms, loc);
    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_get_sol_ker(count, xs, prim_vars_d);
    up->pressure_ker(up->gas_gamma, moms_d, udrift_d, prim_vars_d);

    for (int i=0; i<up->num_basis; i++)
      prim_vars_d[i] *= up->thermalE_fac;

    if (up->is_integrated) {
      out_d[out_coff] = up->integrated_fac * prim_vars_d[0];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_pressure_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_udrift_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_pressure_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_temp_set_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* moms)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve 1 system in each cell)
    long count = linc1;

    const double *moms_d = (const double*) gkyl_array_cfetch(moms, loc);

    up->temp_set_prob_ker(count, As, xs, moms_d, up->gas_gamma, up->mass);
  }
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_temp_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[20]; // 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // Convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // Fetch the correct count in the matrix (since we solve 1 system in each cell)
    long count = linc1;

    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->temp_get_sol_ker(count, xs, prim_vars_d);

    if (up->is_integrated) {
      out_d[out_coff] = up->integrated_fac * prim_vars_d[0];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_temp_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_temp_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_temp_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_udrift_pressure_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[(3+1)*20]; // udrift_comp=3 for 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->udrift_ncomp;

    const double *moms_d = (const double*) gkyl_array_cfetch(moms, loc);
    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    double* udrift_d = prim_vars_d;
    double* pressure_d = &prim_vars_d[up->udrift_ncomp*up->num_basis];

    up->udrift_get_sol_ker(count, xs, udrift_d);
    up->pressure_ker(up->gas_gamma, moms_d, udrift_d, pressure_d);

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+1; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_udrift_pressure_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_udrift_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_udrift_pressure_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_udrift_temp_set_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* moms)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*(up->udrift_ncomp+1);

    const double *moms_d = (const double*) gkyl_array_cfetch(moms, loc);

    up->udrift_temp_set_prob_ker(count, As, xs, moms_d, up->gas_gamma, up->mass);
  }
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_udrift_temp_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[(3+1)*20]; // udrift_comp=3 for 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*(up->udrift_ncomp+1);

    const double* moms_d = (const double*) gkyl_array_cfetch(moms, loc);
    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_temp_get_sol_ker(count, xs, prim_vars_d);

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+1; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_udrift_temp_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_udrift_temp_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_udrift_temp_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_lte_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[(3+2)*20]; // udrift_comp=3 for 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)

  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*(up->udrift_ncomp+1);

    const double* moms_d = (const double*) gkyl_array_cfetch(moms, loc);
    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_temp_get_sol_ker(count, xs, &prim_vars_d[up->num_basis]);

    // Scale rho and temp by 1/m.
    for (int i=0; i<up->num_basis; i++) {
      prim_vars_d[out_coff+i] = moms_d[i]/up->mass;
      prim_vars_d[out_coff+(up->udrift_ncomp+1)*up->num_basis+i] *= 1.0/up->mass;
    }

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+2; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_lte_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_udrift_temp_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_lte_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_flow_energy_set_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* moms)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // Convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // Fetch the correct count in the matrix (since we solve 1 system in each cell).
    long count = linc1;

    const double *moms_d = (const double*) gkyl_array_cfetch(moms, loc);

    up->flowE_set_prob_ker(count, As, xs, moms_d);
  }
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_flow_energy_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[20]; // 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // Convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // Fetch the correct count in the matrix (since we solve 1 system in each cell)
    long count = linc1;

    const double* moms_d = (const double*) gkyl_array_cfetch(moms, loc);
    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->flowE_get_sol_ker(count, xs, prim_vars_d);

    if (up->is_integrated) {
      out_d[out_coff] = up->integrated_fac * prim_vars_d[0];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_flow_energy_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_flow_energy_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_flow_energy_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

__global__ static void
gkyl_gk_neut_fluid_prim_vars_mass_momentum_flow_thermal_energy_copy_cu_kernel(gkyl_gk_neut_fluid_prim_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range, const struct gkyl_array* moms,
  struct gkyl_array* out, int out_coff)
{
  int idx[GKYL_MAX_DIM];

  double prim_vars_buff[(3+3)*20]; // udrift_ncomp=3 for 3D p=2.

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // Invert index from linc1 to idx.
    // Must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange.
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // Convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // Fetch the correct count in the matrix (since we solve 1 system in each cell)
    long count = linc1;

    const double* moms_d = (const double*) gkyl_array_cfetch(moms, loc);
    double* out_d = (double*) gkyl_array_fetch(out, loc);
    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    int fourth_comp_off = 4*up->num_basis;
    for (int i=0; i<fourth_comp_off; i++)
      prim_vars_d[i] = moms_d[i];

    up->flowE_get_sol_ker(count, xs, &prim_vars_d[fourth_comp_off]);

    for (int i=0; i<up->num_basis; i++)
      prim_vars_d[5*up->num_basis+i] = moms_d[fourth_comp_off+i] - prim_vars_d[fourth_comp_off+i];

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+3; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
  }
}

void gkyl_gk_neut_fluid_prim_vars_mass_momentum_flow_thermal_energy_advance_cu(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_gk_neut_fluid_prim_vars_flow_energy_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range, moms->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_gk_neut_fluid_prim_vars_mass_momentum_flow_thermal_energy_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, moms->on_dev, out->on_dev, out_coff);
}

// CUDA kernel to set device pointers to fluid vars kernel functions.
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol.
__global__ static void
gk_neut_fluid_prim_vars_set_cu_dev_ptrs(struct gkyl_gk_neut_fluid_prim_vars *up,
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  up->udrift_set_prob_ker = choose_udrift_set_prob_ker(b_type, cdim, poly_order);
  up->udrift_get_sol_ker = choose_udrift_get_sol_ker(b_type, cdim, poly_order);
  up->pressure_ker = choose_pressure_ker(b_type, cdim, poly_order);
  up->temp_set_prob_ker = choose_temp_set_prob_ker(b_type, cdim, poly_order);
  up->temp_get_sol_ker = choose_temp_get_sol_ker(b_type, cdim, poly_order);
  up->udrift_temp_set_prob_ker = choose_udrift_temp_set_prob_ker(b_type, cdim, poly_order);
  up->udrift_temp_get_sol_ker = choose_udrift_temp_get_sol_ker(b_type, cdim, poly_order);
  up->flowE_set_prob_ker = choose_flowE_set_prob_ker(b_type, cdim, poly_order);
  up->flowE_get_sol_ker = choose_flowE_get_sol_ker(b_type, cdim, poly_order);
}

gkyl_gk_neut_fluid_prim_vars*
gkyl_gk_neut_fluid_prim_vars_cu_dev_new(double gas_gamma, double mass, const struct gkyl_basis* cbasis,
  struct gkyl_rect_grid *grid, const struct gkyl_range *mem_range,
  enum gkyl_gk_neut_fluid_prim_vars_type prim_vars_type, bool is_integrated)
{

  struct gkyl_gk_neut_fluid_prim_vars *up = (struct gkyl_gk_neut_fluid_prim_vars*) gkyl_malloc(sizeof(gkyl_gk_neut_fluid_prim_vars));

  up->gas_gamma = gas_gamma;
  up->mass = mass;

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->num_basis = cbasis->num_basis;
  up->udrift_ncomp = 3;
  up->mem_range = *mem_range;
  up->is_integrated = is_integrated;
  
  up->integrated_fac = 0.0;
  if (up->is_integrated) {
    up->integrated_fac = 1.0;
    for (int d=0; d<up->cdim; d++)
      up->integrated_fac *= (grid->dx[d]/2.0)*sqrt(2.0);
  }

  int nprob;
  up->thermalE_fac = 0.0;
  if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT) {
    nprob = up->udrift_ncomp;
  }
  else if ( (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_PRESSURE) ||
            (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_THERMAL_ENERGY) ) {
    nprob = up->udrift_ncomp;
    up->thermalE_fac = prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_PRESSURE? 1.0 : 1.0/(up->gas_gamma-1.0);
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_TEMP) {
    nprob = 1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_PRESSURE) {
    nprob = up->udrift_ncomp;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_TEMP) {
    nprob = up->udrift_ncomp+1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_LTE) {
    nprob = up->udrift_ncomp+1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_FLOW_ENERGY) {
    nprob = 1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_MASS_MOMENTUM_FLOW_THERMAL_ENERGY) {
    nprob = 1;
  }

  // There are udrift_ncomp*range->volume linear systems to be solved
  // for 3 components of u: ux, uy, uz.
  up->As = gkyl_nmat_cu_dev_new(nprob*mem_range->volume, up->num_basis, up->num_basis);
  up->xs = gkyl_nmat_cu_dev_new(nprob*mem_range->volume, up->num_basis, 1);
  if (up->poly_order > 1) {
    up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_gk_neut_fluid_prim_vars *up_cu = (struct gkyl_gk_neut_fluid_prim_vars*) gkyl_cu_malloc(sizeof(gkyl_gk_neut_fluid_prim_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_gk_neut_fluid_prim_vars), GKYL_CU_MEMCPY_H2D);

  gk_neut_fluid_prim_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;

}
