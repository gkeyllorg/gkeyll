/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_em_coupling.h>
#include <gkyl_dg_calc_pkpm_em_coupling_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_pkpm_em_coupling_set_one_fluid_cu_kernel(gkyl_dg_calc_pkpm_em_coupling* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range, double dt, 
  const struct gkyl_array* app_accel, 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, 
  struct gkyl_array* em)
{
  const double *app_accels[GKYL_MAX_SPECIES]; 
  const double *pkpm_moms[GKYL_MAX_SPECIES];
  const double *pkpm_flows[GKYL_MAX_SPECIES];
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

    app_accels[0] = (const double*) gkyl_array_cfetch(app_accel, loc);
    pkpm_moms[0] = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    pkpm_flows[0] = (const double*) gkyl_array_cfetch(pkpm_u, loc);

    const double *ext_em_d = (const double*) gkyl_array_cfetch(ext_em, loc);
    const double *app_current_d = (const double*) gkyl_array_cfetch(app_current, loc);
    double *em_d = (double*) gkyl_array_fetch(em, loc);

    up->pkpm_em_coupling_set(linc1, 
      up->num_species, up->qbym, up->epsilon0, up->pkpm_field_static, dt, 
      As, xs, 
      app_accels, ext_em_d, app_current_d, pkpm_moms, pkpm_flows, em_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_em_coupling_copy_one_fluid_cu_kernel(gkyl_dg_calc_pkpm_em_coupling* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u,   
  struct gkyl_array* euler_pkpm, struct gkyl_array* em)
{
  const double *pkpm_moms[GKYL_MAX_SPECIES];
  const double *pkpm_flows[GKYL_MAX_SPECIES];  
  double *fluids[GKYL_MAX_SPECIES];  
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

    pkpm_moms[0] = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    pkpm_flows[0] = (const double*) gkyl_array_cfetch(pkpm_u, loc);
    fluids[0] = (double*) gkyl_array_fetch(euler_pkpm, loc);
    double *em_d = (double*) gkyl_array_fetch(em, loc);

    up->pkpm_em_coupling_copy(linc1, up->num_species, up->qbym, up->epsilon0, 
      xs, pkpm_moms, pkpm_flows, fluids, em_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_em_coupling_set_two_fluids_cu_kernel(gkyl_dg_calc_pkpm_em_coupling* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range, double dt, 
  const struct gkyl_array* app_accel_1, const struct gkyl_array* app_accel_2, 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  const struct gkyl_array* vlasov_pkpm_moms_1, const struct gkyl_array* vlasov_pkpm_moms_2, 
  const struct gkyl_array* pkpm_u_1, const struct gkyl_array* pkpm_u_2, 
  struct gkyl_array* em)
{
  const double *app_accels[GKYL_MAX_SPECIES];  
  const double *pkpm_moms[GKYL_MAX_SPECIES];
  const double *pkpm_flows[GKYL_MAX_SPECIES];
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

    app_accels[0] = (const double*) gkyl_array_cfetch(app_accel_1, loc);
    app_accels[1] = (const double*) gkyl_array_cfetch(app_accel_2, loc);
    pkpm_moms[0] = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms_1, loc);
    pkpm_moms[1] = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms_2, loc);
    pkpm_flows[0] = (const double*) gkyl_array_cfetch(pkpm_u_1, loc);
    pkpm_flows[1] = (const double*) gkyl_array_cfetch(pkpm_u_2, loc);

    const double *ext_em_d = (const double*) gkyl_array_cfetch(ext_em, loc);
    const double *app_current_d = (const double*) gkyl_array_cfetch(app_current, loc);
    double *em_d = (double*) gkyl_array_fetch(em, loc);

    up->pkpm_em_coupling_set(linc1, 
      up->num_species, up->qbym, up->epsilon0, up->pkpm_field_static, dt, 
      As, xs, 
      app_accels, ext_em_d, app_current_d, pkpm_moms, pkpm_flows, em_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_em_coupling_copy_two_fluids_cu_kernel(gkyl_dg_calc_pkpm_em_coupling* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* vlasov_pkpm_moms_1, const struct gkyl_array* vlasov_pkpm_moms_2, 
  const struct gkyl_array* pkpm_u_1, const struct gkyl_array* pkpm_u_2,   
  struct gkyl_array* euler_pkpm_1, struct gkyl_array* euler_pkpm_2, struct gkyl_array* em)
{
  const double *pkpm_moms[GKYL_MAX_SPECIES];
  const double *pkpm_flows[GKYL_MAX_SPECIES];  
  double *fluids[GKYL_MAX_SPECIES];  
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

    pkpm_moms[0] = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms_1, loc);
    pkpm_moms[1] = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms_2, loc);
    pkpm_flows[0] = (const double*) gkyl_array_cfetch(pkpm_u_1, loc);
    pkpm_flows[1] = (const double*) gkyl_array_cfetch(pkpm_u_2, loc);
    fluids[0] = (double*) gkyl_array_fetch(euler_pkpm_1, loc);
    fluids[1] = (double*) gkyl_array_fetch(euler_pkpm_2, loc);
    double *em_d = (double*) gkyl_array_fetch(em, loc);

    up->pkpm_em_coupling_copy(linc1, up->num_species, up->qbym, up->epsilon0, 
      xs, pkpm_moms, pkpm_flows, fluids, em_d);
  }
}

// Host-side wrapper for primitive variable calculation
void gkyl_dg_calc_pkpm_em_coupling_advance_cu(struct gkyl_dg_calc_pkpm_em_coupling *up, double dt, 
  const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  const struct gkyl_array* vlasov_pkpm_moms[GKYL_MAX_SPECIES], const struct gkyl_array* pkpm_u[GKYL_MAX_SPECIES], 
  struct gkyl_array* euler_pkpm[GKYL_MAX_SPECIES], struct gkyl_array* em)
{
  struct gkyl_range conf_range = up->mem_range;

  int num_species = up->num_species;

  if (num_species == 1) {
    gkyl_dg_calc_pkpm_em_coupling_set_one_fluid_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev, 
      up->As->on_dev, up->xs->on_dev, conf_range, dt, 
      app_accel[0]->on_dev, ext_em->on_dev, app_current->on_dev, 
      vlasov_pkpm_moms[0]->on_dev, pkpm_u[0]->on_dev, 
      em->on_dev);
  }
  else if (num_species == 2) {
    gkyl_dg_calc_pkpm_em_coupling_set_two_fluids_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev, 
      up->As->on_dev, up->xs->on_dev, conf_range, dt, 
      app_accel[0]->on_dev, app_accel[1]->on_dev, ext_em->on_dev, app_current->on_dev, 
      vlasov_pkpm_moms[0]->on_dev, vlasov_pkpm_moms[1]->on_dev, 
      pkpm_u[0]->on_dev, pkpm_u[1]->on_dev, 
      em->on_dev);
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  if (num_species == 1) {
    gkyl_dg_calc_pkpm_em_coupling_copy_one_fluid_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
      up->xs->on_dev, conf_range, 
      vlasov_pkpm_moms[0]->on_dev, pkpm_u[0]->on_dev, 
      euler_pkpm[0]->on_dev, em->on_dev);
  }
  else if (num_species == 2) {
    gkyl_dg_calc_pkpm_em_coupling_copy_two_fluids_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
      up->xs->on_dev, conf_range, 
      vlasov_pkpm_moms[0]->on_dev, vlasov_pkpm_moms[1]->on_dev, 
      pkpm_u[0]->on_dev, pkpm_u[1]->on_dev,       
      euler_pkpm[0]->on_dev, euler_pkpm[1]->on_dev, em->on_dev);
  }
}

// CUDA kernel to set device pointers to pkpm-em coupling kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_pkpm_em_coupling_set_cu_dev_ptrs(struct gkyl_dg_calc_pkpm_em_coupling *up, 
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  up->pkpm_em_coupling_set = choose_pkpm_em_coupling_set_kern(b_type, cdim, poly_order);
  up->pkpm_em_coupling_copy = choose_pkpm_em_coupling_copy_kern(b_type, cdim, poly_order);
}

gkyl_dg_calc_pkpm_em_coupling*
gkyl_dg_calc_pkpm_em_coupling_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *mem_range, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  bool pkpm_field_static)
{
  struct gkyl_dg_calc_pkpm_em_coupling *up = (struct gkyl_dg_calc_pkpm_em_coupling*) gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_em_coupling));

  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->mem_range = *mem_range;

  // Linear system size is nc*(3*num_species + 3)
  up->num_species = num_species;
  up->As = gkyl_nmat_cu_dev_new(mem_range->volume, nc*(3*up->num_species + 3), nc*(3*up->num_species + 3));
  up->xs = gkyl_nmat_cu_dev_new(mem_range->volume, nc*(3*up->num_species + 3), 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  // Boolean for whether or not self-consistent EM fields are static
  up->pkpm_field_static = pkpm_field_static;

  // Needed constants for the source solve
  up->epsilon0 = epsilon0;
  for (int n = 0; n < num_species; ++n) {
    up->qbym[n] = qbym[n];
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_em_coupling *up_cu = (struct gkyl_dg_calc_pkpm_em_coupling*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_em_coupling));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_em_coupling), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_em_coupling_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;
}
