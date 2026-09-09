/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars_priv.h>
#include <gkyl_wv_canonical_pb_fluid.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf_cu_kernel(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range conf_range, const struct gkyl_range conf_ext_range, 
  struct gkyl_array *phi,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);  
    long loc = gkyl_range_idx(&conf_range, idx);
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xc);

    double* alpha_surf_d = (double*) gkyl_array_fetch(alpha_surf, loc);
    double* sgn_alpha_surf_d = (double*) gkyl_array_fetch(sgn_alpha_surf, loc);
    int* const_sgn_alpha_d = (int*) gkyl_array_fetch(const_sgn_alpha, loc);
    for (int dir = 0; dir<cdim; ++dir) {
      const_sgn_alpha_d[dir] = up->alpha_surf[dir](xc, up->conf_grid.dx, 
        (const double*) gkyl_array_cfetch(phi, loc),
        alpha_surf_d, sgn_alpha_surf_d);

      // If the configuration space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == conf_range.upper[dir]) {
        gkyl_copy_int_arr(cdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_ext = gkyl_range_idx(&conf_ext_range, idx_edge);

        double* alpha_surf_ext_d = (double*) gkyl_array_fetch(alpha_surf, loc_ext);
        double* sgn_alpha_surf_ext_d = (double*) gkyl_array_fetch(sgn_alpha_surf, loc_ext);
        int* const_sgn_alpha_ext_d = (int*) gkyl_array_fetch(const_sgn_alpha, loc_ext);
        const_sgn_alpha_ext_d[dir] = up->alpha_edge_surf[dir](xc, up->conf_grid.dx, 
          (const double*) gkyl_array_fetch(phi, loc),
          alpha_surf_ext_d, sgn_alpha_surf_ext_d);
      }  
    }
  }
}
// Host-side wrapper for configuration-space surface alpha computation. 
void 
gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf_cu(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_array *phi,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *conf_ext_range, 
    phi->on_dev, alpha_surf->on_dev, sgn_alpha_surf->on_dev, const_sgn_alpha->on_dev);
}

__global__ void
gkyl_canonical_pb_fluid_vars_subtract_zonal_cu_kernel(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range x_range, 
  const struct gkyl_array *phi_zonal, const struct gkyl_array *n_zonal, 
  struct gkyl_array *adiabatic_coupling_phi_n)
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
    long loc_1d = gkyl_range_idx(&x_range, idx);  

    const double *phi_zonal_d = (const double*) gkyl_array_cfetch(phi_zonal, loc_1d);
    const double *n_zonal_d = (const double*) gkyl_array_cfetch(n_zonal, loc_1d);

    double *adiabatic_coupling_phi_n_d = (double*) gkyl_array_fetch(adiabatic_coupling_phi_n, loc); 
    up->subtract_zonal(phi_zonal_d, n_zonal_d, adiabatic_coupling_phi_n_d); 
  }
}

__global__ void
gkyl_canonical_pb_fluid_vars_source_cu_kernel(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  struct gkyl_range conf_range, 
  const struct gkyl_array *phi, const struct gkyl_array *n0, 
  const struct gkyl_array *adiabatic_coupling_phi_n, struct gkyl_array *rhs)
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

    const double *phi_d = (const double*) gkyl_array_cfetch(phi, loc);
    const double *n0_d = (const double*) gkyl_array_cfetch(n0, loc);

    double* rhs_d = (double*) gkyl_array_fetch(rhs, loc);
    up->canonical_pb_fluid_source(up->conf_grid.dx, up->alpha, phi_d, n0_d, 
      adiabatic_coupling_phi_n ? (const double*) gkyl_array_cfetch(adiabatic_coupling_phi_n, loc) : 0, 
      rhs_d);
  }
}

// Host-side wrapper for source update of canonical PB fluid systems.
void 
gkyl_canonical_pb_fluid_vars_source_cu(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array *phi, const struct gkyl_array *n0, 
  const struct gkyl_array *fluid, struct gkyl_array *rhs)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;

  // If alpha is specified, we are solving Hasegawa-Wakatani and need to check 
  // whether we are solving the modified version of Hasegawa-Wakatani which 
  // requires computing the zonal components of phi, n: f_zonal = 1/Ly int f dy
  // and subtracting the zonal components off the adiabatic coupling, f_tilde = f - f_zonal. 
  // Otherwise, we just copy n and phi into a temporary array for use in the updater.
  if (up->alpha > 0.0) {
    gkyl_array_set_offset(up->n, 1.0, fluid, up->conf_basis.num_basis);
    gkyl_array_set_offset(up->adiabatic_coupling_phi_n, 1.0, phi, 0);
    gkyl_array_set_offset(up->adiabatic_coupling_phi_n, 1.0, up->n, up->conf_basis.num_basis);
    if (up->is_modified) {
      // Compute the zonal components of phi and n. 
      gkyl_array_average_advance(up->int_y, phi, up->phi_zonal);
      gkyl_array_average_advance(up->int_y, up->n, up->n_zonal);
      gkyl_canonical_pb_fluid_vars_subtract_zonal_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
        *conf_range, up->x_local, up->phi_zonal->on_dev, up->n_zonal->on_dev, 
        up->adiabatic_coupling_phi_n->on_dev);
    }
  }

  gkyl_canonical_pb_fluid_vars_source_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, phi->on_dev, n0->on_dev, 
    up->adiabatic_coupling_phi_n ? up->adiabatic_coupling_phi_n->on_dev : 0, 
    rhs->on_dev);
}


// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_canoncial_pb_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_canonical_pb_fluid_vars *up,
  enum gkyl_basis_type b_type, int cdim, int poly_order, enum gkyl_eqn_type eqn_type, bool is_modified)
{
  if (eqn_type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA) {
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_hasegawa_mima_source_kern(b_type, cdim, poly_order);
  }
  else if (eqn_type == GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI) {
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_hasegawa_wakatani_source_kern(b_type, cdim, poly_order);
    if (is_modified) {
      up->subtract_zonal = choose_canonical_pb_fluid_subtract_zonal_kern(b_type, cdim, poly_order);
    }
  }
  else {
    // Default source kernel; immediately returns and does not do anything. 
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_default_source_kern(b_type, cdim, poly_order);
  }
  for (int d=0; d<cdim; ++d) {
    up->alpha_surf[d] = choose_canonical_pb_fluid_alpha_surf_kern(b_type, d, cdim, poly_order);
    up->alpha_edge_surf[d] = choose_canonical_pb_fluid_alpha_edge_surf_kern(b_type, d, cdim, poly_order);
  }
}


gkyl_dg_calc_canonical_pb_fluid_vars*
gkyl_dg_calc_canonical_pb_fluid_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_wv_eqn *wv_eqn)
{   
  struct gkyl_dg_calc_canonical_pb_fluid_vars *up = (struct gkyl_dg_calc_canonical_pb_fluid_vars *) gkyl_malloc(sizeof(gkyl_dg_calc_canonical_pb_fluid_vars));

  up->conf_grid = *conf_grid;
  up->conf_basis = *conf_basis; 
  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  up->cdim = cdim;
  up->alpha = 0.0;
  up->is_modified = 0; 

  up->adiabatic_coupling_phi_n = 0;
  if (wv_eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI) {
    up->alpha = gkyl_wv_can_pb_hasegawa_wakatani_alpha(wv_eqn); 
    up->is_modified = gkyl_wv_can_pb_hasegawa_wakatani_is_modified(wv_eqn); 

    // Temporary array for holding the density and combined potential and density for computing the adiabatic coupling.
    // These are stored separately from the input phi and fluid arrays in case we are solving the
    // modified Hasegawa-Wakatani system and need to subtract the zonal components. 
    up->n = gkyl_array_cu_dev_new(GKYL_DOUBLE, conf_basis->num_basis, conf_ext_range->volume);
    // Component 0 is n, Component 1 is phi. 
    up->adiabatic_coupling_phi_n = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*conf_basis->num_basis, conf_ext_range->volume);

    // Set up the array averaging to correctly subtact the zonal component of fluctuations
    // Currently assumes that updater has *no decomposition* in y and thus owns the whole y range.
    if (up->is_modified) {
      // Make the one-dimensional x basis, and ranges for constructing the average. 
      struct gkyl_basis basis_x;
      gkyl_cart_modal_serendip(&basis_x, 1, poly_order);
      gkyl_range_init(&up->x_local, 1, &conf_range->lower[0], &conf_range->upper[0]);
      gkyl_range_init(&up->x_local_ext, 1, &conf_ext_range->lower[0], &conf_ext_range->upper[0]);
      // Integration over y only, (x,y) to (x).
      int int_dim_y[] = {0,1,0};
      struct gkyl_array_average_inp inp_int_y = {
        .grid = &up->conf_grid,
        .basis = up->conf_basis,
        .basis_avg = basis_x,
        .local = conf_range,
        .local_avg = &up->x_local,
        .local_avg_ext = &up->x_local_ext,
        .weight = NULL,
        .avg_dim = int_dim_y,
        .use_gpu = true, // We will perform the average on GPUs
      };
      up->int_y = gkyl_array_average_inew(&inp_int_y);
      up->phi_zonal = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis_x.num_basis, up->x_local_ext.volume);
      up->n_zonal = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis_x.num_basis, up->x_local_ext.volume);
    }    
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_canonical_pb_fluid_vars *up_cu = (struct gkyl_dg_calc_canonical_pb_fluid_vars *) gkyl_cu_malloc(sizeof(gkyl_dg_calc_canonical_pb_fluid_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_canonical_pb_fluid_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_canoncial_pb_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, conf_basis->b_type, cdim, poly_order, wv_eqn->type, up->is_modified);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
