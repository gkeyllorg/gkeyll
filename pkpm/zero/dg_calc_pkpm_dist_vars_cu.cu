/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_dist_vars.h>
#include <gkyl_dg_calc_pkpm_dist_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_calc_pkpm_dist_vars_mirror_force_cu_kernel(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array* pkpm_prim, const struct gkyl_array* nu_prim_moms_sum, 
  const struct gkyl_array* div_b, const struct gkyl_array* pkpm_accel, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1)
{
  double xc[GKYL_MAX_DIM] = {0.0};
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

    const double *pkpm_prim_d = (const double*) gkyl_array_cfetch(pkpm_prim, loc_conf);
    const double *nu_prim_moms_sum_d = (const double*) gkyl_array_cfetch(nu_prim_moms_sum, loc_conf);
    const double *div_b_d = (const double*) gkyl_array_cfetch(div_b, loc_conf);
    const double *pkpm_accel_d = (const double*) gkyl_array_cfetch(pkpm_accel, loc_conf);
    const double *fIn_d = (const double*) gkyl_array_cfetch(fIn, loc_phase);
    const double *F_k_p_1_d = (const double*) gkyl_array_cfetch(F_k_p_1, loc_phase);

    double *g_dist_source_d = (double*) gkyl_array_fetch(g_dist_source, loc_phase);
    double *F_k_m_1_d = (double*) gkyl_array_fetch(F_k_m_1, loc_phase);

    up->pkpm_dist_mirror_force(xc, up->phase_grid.dx, 
      pkpm_prim_d, nu_prim_moms_sum_d, div_b_d, pkpm_accel_d, 
      fIn_d, F_k_p_1_d, g_dist_source_d, F_k_m_1_d);
  }  
}
// Host-side wrapper for pkpm mirror force source distribution function calculation
void 
gkyl_dg_calc_pkpm_dist_vars_mirror_force_cu(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* pkpm_prim, const struct gkyl_array* nu_prim_moms_sum, 
  const struct gkyl_array* div_b, const struct gkyl_array* pkpm_accel, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_pkpm_dist_vars_mirror_force_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, 
    pkpm_prim->on_dev, nu_prim_moms_sum->on_dev, 
    div_b->on_dev, pkpm_accel->on_dev, 
    fIn->on_dev, F_k_p_1->on_dev, 
    g_dist_source->on_dev, F_k_m_1->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_dist_vars_div_ppar_cu_kernel(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, const struct gkyl_array* fIn, 
  const struct gkyl_array* max_b, struct gkyl_array* pkpm_div_ppar)
{
  double xc[GKYL_MAX_DIM] = {0.0};
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idxc);
    gkyl_rect_grid_cell_center(&up->phase_grid, idxc, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc_conf = gkyl_range_idx(&conf_range, idxc);
    long linc_phase = gkyl_range_idx(&phase_range, idxc);

    const double *bvar_surf_c = (const double*) gkyl_array_cfetch(bvar_surf, linc_conf);
    const double *bvar_c = (const double*) gkyl_array_cfetch(bvar, linc_conf);
    const double *f_c = (const double*) gkyl_array_cfetch(fIn, linc_phase);
    const double *max_b_c = (const double*) gkyl_array_cfetch(max_b, linc_conf);

    double momLocal[96]; // hard-coded to 3 * max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k=0; k<96; ++k)
      momLocal[k] = 0.0;

    for (int dir=0; dir<up->cdim; ++dir) {
      gkyl_copy_int_arr(up->cdim+1, idxc, idxl);
      gkyl_copy_int_arr(up->cdim+1, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl_conf = gkyl_range_idx(&conf_range, idxl); 
      long linl_phase = gkyl_range_idx(&phase_range, idxl); 
      long linr_conf = gkyl_range_idx(&conf_range, idxr); 
      long linr_phase = gkyl_range_idx(&phase_range, idxr); 

      const double *bvar_surf_l = (const double*) gkyl_array_cfetch(bvar_surf, linl_conf);
      const double *f_l = (const double*) gkyl_array_cfetch(fIn, linl_phase);
      const double *bvar_surf_r = (const double*) gkyl_array_cfetch(bvar_surf, linr_conf);
      const double *f_r = (const double*) gkyl_array_cfetch(fIn, linr_phase);

      up->pkpm_dist_div_ppar[dir](xc, up->phase_grid.dx, 
        bvar_surf_l, bvar_surf_c, bvar_surf_r, 
        f_l, f_c, f_r, 
        bvar_c, max_b_c, &momLocal[0]);
    }
    // Accumulate output to output array atomically to avoid race conditions
    double *pkpm_div_ppar_d = (double*) gkyl_array_fetch(pkpm_div_ppar, linc_conf);
    for (unsigned int k = 0; k < pkpm_div_ppar->ncomp; ++k) {
       atomicAdd(&pkpm_div_ppar_d[k], momLocal[k]);
    }    
  }  
}
// Host-side wrapper for pkpm div(p_parallel b_hat) calculation
void 
gkyl_dg_calc_pkpm_dist_vars_div_ppar_cu(struct gkyl_dg_calc_pkpm_dist_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, const struct gkyl_array* fIn, 
  const struct gkyl_array* max_b, struct gkyl_array* pkpm_div_ppar)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_pkpm_dist_vars_div_ppar_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, 
    bvar_surf->on_dev, bvar->on_dev, fIn->on_dev, 
    max_b->on_dev, pkpm_div_ppar->on_dev);
}

// CUDA kernel to set device pointers to pkpm dist vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_pkpm_dist_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_pkpm_dist_vars *up, enum gkyl_basis_type b_type,
  int cdim,int poly_order)
{
  up->pkpm_dist_mirror_force = choose_pkpm_dist_mirror_force_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    up->pkpm_dist_div_ppar[d] = choose_pkpm_dist_div_ppar_kern(d, b_type, cdim, poly_order);
}

gkyl_dg_calc_pkpm_dist_vars*
gkyl_dg_calc_pkpm_dist_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis* cbasis)
{
  struct gkyl_dg_calc_pkpm_dist_vars *up = (struct gkyl_dg_calc_pkpm_dist_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_dist_vars));

  up->phase_grid = *phase_grid;
  enum gkyl_basis_type b_type = cbasis->b_type;
  int cdim = cbasis->ndim;
  up->cdim = cdim;
  int poly_order = cbasis->poly_order;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_dist_vars *up_cu = (struct gkyl_dg_calc_pkpm_dist_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_dist_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_dist_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_dist_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
