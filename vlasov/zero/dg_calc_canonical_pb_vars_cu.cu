/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars_priv.h>
#include <gkyl_util.h>
}


__global__ void
gkyl_dg_calc_canonical_pb_vars_alpha_surf_cu_kernel(struct gkyl_dg_calc_canonical_pb_vars *up, 
  const struct gkyl_range conf_range, const struct gkyl_range phase_range,  const struct gkyl_range phase_ext_range, 
  struct gkyl_array *hamil,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);  
    long loc_phase = gkyl_range_idx(&phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    double* alpha_surf_d = (double*) gkyl_array_fetch(alpha_surf, loc_phase);
    double* sgn_alpha_surf_d = (double*) gkyl_array_fetch(sgn_alpha_surf, loc_phase);
    int* const_sgn_alpha_d = (int*) gkyl_array_fetch(const_sgn_alpha, loc_phase);
    for (int dir = 0; dir<vdim; ++dir) {
      const_sgn_alpha_d[dir+cdim] = up->alpha_surf[dir+cdim](xc, up->phase_grid.dx, 
        (const double*) gkyl_array_cfetch(hamil, loc_phase),
        alpha_surf_d, sgn_alpha_surf_d);
    }
    for (int dir = 0; dir<cdim; ++dir) {
      const_sgn_alpha_d[dir] = up->alpha_surf[dir](xc, up->phase_grid.dx, 
        (const double*) gkyl_array_cfetch(hamil, loc_phase),
        alpha_surf_d, sgn_alpha_surf_d);

      // If the phase space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == conf_range.upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_phase_ext = gkyl_range_idx(&phase_ext_range, idx_edge);

        double* alpha_surf_ext_d = (double*) gkyl_array_fetch(alpha_surf, loc_phase_ext);
        double* sgn_alpha_surf_ext_d = (double*) gkyl_array_fetch(sgn_alpha_surf, loc_phase_ext);
        int* const_sgn_alpha_ext_d = (int*) gkyl_array_fetch(const_sgn_alpha, loc_phase_ext);
        const_sgn_alpha_ext_d[dir] = up->alpha_edge_surf[dir](xc, up->phase_grid.dx, 
          (const double*) gkyl_array_fetch(hamil, loc_phase),
          alpha_surf_ext_d, sgn_alpha_surf_ext_d);
      }  
    }
  }
}
// Host-side wrapper
void 
gkyl_dg_calc_canonical_pb_vars_alpha_surf_cu(struct gkyl_dg_calc_canonical_pb_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,  const struct gkyl_range *phase_ext_range, 
  struct gkyl_array *hamil,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_canonical_pb_vars_alpha_surf_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, *phase_ext_range, 
    hamil->on_dev, alpha_surf->on_dev, sgn_alpha_surf->on_dev, const_sgn_alpha->on_dev);
}

/* Convert the bulk velocity from contravariant to covaraint components for can-pb*/
__global__ void
gkyl_canonical_pb_contra_to_covariant_m1i_cu_kernel(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range conf_range,
 const struct gkyl_array *h_ij,  const struct gkyl_array *V_drift, const struct gkyl_array *M1i, struct gkyl_array *V_drift_cov,
 struct gkyl_array *M1i_cov)
{
  int cdim = up->cdim;
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

    const double *h_ij_d = (const double*) gkyl_array_cfetch(h_ij, loc);
    const double *v_j_d = (const double*) gkyl_array_cfetch(V_drift, loc);
    const double *nv_i_d = (const double*) gkyl_array_cfetch(M1i, loc);
    double *v_i_cov_d = (double*) gkyl_array_fetch(V_drift_cov, loc);
    double *nv_i_cov_d = (double*) gkyl_array_fetch(M1i_cov, loc);

    up->canonical_pb_covariant_u_i(h_ij_d, v_j_d, nv_i_d, v_i_cov_d, nv_i_cov_d);
  }
}
// Host-side wrapper
void 
gkyl_canonical_pb_contra_to_covariant_m1i_cu(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij, const struct gkyl_array *V_drift, const struct gkyl_array *M1i, struct gkyl_array *V_drift_cov, 
 struct gkyl_array *M1i_cov)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_canonical_pb_contra_to_covariant_m1i_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, h_ij->on_dev, V_drift->on_dev, M1i->on_dev, V_drift_cov->on_dev, M1i_cov->on_dev);
}


/* Compute the pressure for can-pb*/
__global__ void
gkyl_canonical_pb_pressure_cu_kernel(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range conf_range,
 const struct gkyl_array *h_ij_inv, 
 const struct gkyl_array *MEnergy, const struct gkyl_array *V_drift, const struct gkyl_array *M1i,
 struct gkyl_array *pressure)
{
  int cdim = up->cdim;
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

    const double *h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv, loc);
    const double *MEnergy_d = (const double*) gkyl_array_cfetch(MEnergy, loc);
    const double *v_j_d = (const double*) gkyl_array_cfetch(V_drift, loc);
    const double *nv_i_d = (const double*) gkyl_array_cfetch(M1i, loc);

    double* d_Jv_P_d = (double*) gkyl_array_fetch(pressure, loc);

    up->canonical_pb_pressure(h_ij_inv_d, MEnergy_d, v_j_d, nv_i_d, d_Jv_P_d);
  }
}
// Host-side wrapper
void 
gkyl_canonical_pb_pressure_cu(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij_inv, 
 const struct gkyl_array *MEnergy, const struct gkyl_array *V_drift, const struct gkyl_array *M1i,
 struct gkyl_array *pressure)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_canonical_pb_pressure_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, h_ij_inv->on_dev, MEnergy->on_dev, V_drift->on_dev, M1i->on_dev, pressure->on_dev);
}


// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
  dg_calc_canoncial_pb_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_canonical_pb_vars *up,
  enum gkyl_basis_type b_type, int cv_index, int cdim, int vdim, int poly_order)
{
  up->canonical_pb_pressure = choose_canonical_pb_pressure_kern(b_type, cv_index, cdim, poly_order);
  up->canonical_pb_covariant_u_i = choose_canonical_pb_m1i_contra_to_cov_kern(b_type, cv_index, cdim, poly_order);
  for (int d=0; d<cdim; ++d) {  
    up->alpha_surf[d] = choose_canonical_pb_alpha_surf_kern(b_type, d, cv_index, cdim, vdim, poly_order);
    up->alpha_edge_surf[d] = choose_canonical_pb_alpha_edge_surf_kern(b_type, d, cv_index, cdim, vdim, poly_order);
  }
  for (int d=0; d<vdim; ++d) {  
    up->alpha_surf[d+cdim] = choose_canonical_pb_alpha_surf_v_kern(b_type, d, cv_index, cdim, vdim, poly_order);
  }
}


gkyl_dg_calc_canonical_pb_vars*
gkyl_dg_calc_canonical_pb_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis)
{   
  struct gkyl_dg_calc_canonical_pb_vars *up = (struct gkyl_dg_calc_canonical_pb_vars *) gkyl_malloc(sizeof(gkyl_dg_calc_canonical_pb_vars));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_canonical_pb_vars *up_cu = (struct gkyl_dg_calc_canonical_pb_vars *) gkyl_cu_malloc(sizeof(gkyl_dg_calc_canonical_pb_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_canonical_pb_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_canoncial_pb_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, phase_basis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
