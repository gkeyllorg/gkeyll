#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis_priv.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ static void
gkyl_proj_maxwellian_pots_on_basis_conf_quad_ker(struct gkyl_range conf_range,
  const struct gkyl_array *conf_basis_at_ords, const struct gkyl_array *prim_moms,
  struct gkyl_array *prim_moms_conf_quad)
{
  int vdim = 3;
  int num_conf_basis = conf_basis_at_ords->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_conf_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long linc = gkyl_range_idx(&conf_range, cidx);

    const double *prim_moms_d = (const double *)gkyl_array_cfetch(prim_moms, linc);
    double *prim_moms_conf_quad_d = (double *)gkyl_array_fetch(prim_moms_conf_quad, linc);

    // Sum over basis of primitive moments (n, u_drift, T/m) at configuration space points
    const double *conf_basis_at_ord_q = (const double *) gkyl_array_cfetch(conf_basis_at_ords, linq);
    for (int k=0; k<num_conf_basis; ++k) {
      for (int d=0; d<vdim+2; ++d) {
        prim_moms_conf_quad_d[d*tot_conf_quad+linq] += prim_moms_d[num_conf_basis*d+k]*conf_basis_at_ord_q[k];
      }
    }
  }
}

__global__ static void
gkyl_proj_h_on_basis_phase_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, const struct gkyl_array *ordinates,
  const struct gkyl_array *conf_basis_at_ords,
  const struct gkyl_array *prim_moms_conf_quad, int *p2c_qidx, struct gkyl_array *pot_phase_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_conf_quad = conf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_phase_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *pot_q = (double *)gkyl_array_fetch(pot_phase_quad, linp);

    comp_to_phys(pdim, (const double*)gkyl_array_cfetch(ordinates, linq),
      phase_grid.dx, xc, &xmu[0]);
    int cqidx = p2c_qidx[linq];

    double rel_speedsq_d = 0.0;
    for (int d=0; d<vdim; ++d) {
      rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
    }
    double rel_speed_d = sqrt(rel_speedsq_d);

    pot_q[linq] += eval_fpo_h(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
  }
}

__global__ static void
gkyl_proj_g_on_basis_phase_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, const struct gkyl_array *ordinates,
  const struct gkyl_array *conf_basis_at_ords,
  const struct gkyl_array *prim_moms_conf_quad, int *p2c_qidx, struct gkyl_array *pot_phase_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_conf_quad = conf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_phase_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *pot_q = (double *)gkyl_array_fetch(pot_phase_quad, linp);

    comp_to_phys(pdim, (const double*)gkyl_array_cfetch(ordinates, linq),
      phase_grid.dx, xc, &xmu[0]);
    int cqidx = p2c_qidx[linq];

    double rel_speedsq_d = 0.0;
    for (int d=0; d<vdim; ++d) {
      rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
    }
    double rel_speed_d = sqrt(rel_speedsq_d);

    pot_q[linq] += eval_fpo_g(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
  }
}

__global__ static void
gkyl_proj_h_on_basis_surf_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, int tot_conf_quad,
  const struct gkyl_array *surf_basis_at_ords, const struct gkyl_array *surf_ordinates,
  const struct gkyl_array *prim_moms_conf_quad, int *surf2c_qidx,
  struct gkyl_array *pot_surf_quad, struct gkyl_array *pot_deriv_surf_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_surf_quad = surf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_surf_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *pot_q = (double *)gkyl_array_fetch(pot_surf_quad, linp);
    double *pot_deriv_q = (double *)gkyl_array_fetch(pot_deriv_surf_quad, linp);

    // Iterate over velocity directions for H, dH/dv, G, d2G/dv2 at each dir1 domain boundary
    // Derivatives dH/dv and d2G/dv2 are normal to surface
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.upper[dir1]);

      const double* surf_ord = (const double*)gkyl_array_cfetch(surf_ordinates, linq);
      int cqidx = surf2c_qidx[linq];

      if (is_edge_in_dir1) {
        // Velocity value at dir1 boundary
        double vbound = pidx[dir1] == phase_range.lower[dir1] ?
          phase_grid.lower[dir1] : phase_grid.upper[dir1];

        surf_comp_to_phys(dir1, pdim, surf_ord, phase_grid.dx, xc, xmu);
        xmu[dir1] = vbound;

        double rel_speedsq_d = 0.0;
        double rel_vel_in_dir1_d = xmu[dir1] - u_drift_d[tot_conf_quad*d1+cqidx];
        for (int d=0; d<vdim; ++d) {
          rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
        }
        double rel_speed_d = sqrt(rel_speedsq_d);

        pot_q[d1*tot_surf_quad + linq] += eval_fpo_h(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
        pot_deriv_q[d1*tot_surf_quad + linq] += eval_fpo_dhdv(den_d[cqidx], rel_vel_in_dir1_d, vtsq_d[cqidx], rel_speed_d);
      }
    }
  }
}

__global__ static void
gkyl_proj_g_on_basis_surf_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, int tot_conf_quad,
  const struct gkyl_array *surf_basis_at_ords, const struct gkyl_array *surf_ordinates,
  const struct gkyl_array *prim_moms_conf_quad, int *surf2c_qidx,
  struct gkyl_array *pot_surf_quad, struct gkyl_array *pot_deriv_surf_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_surf_quad = surf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_surf_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *pot_q = (double *)gkyl_array_fetch(pot_surf_quad, linp);
    double *pot_deriv_q = (double *)gkyl_array_fetch(pot_deriv_surf_quad, linp);

    // Iterate over velocity directions for H, dH/dv, G, d2G/dv2 at each dir1 domain boundary
    // Derivatives dH/dv and d2G/dv2 are normal to surface
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.upper[dir1]);

      const double* surf_ord = (const double*)gkyl_array_cfetch(surf_ordinates, linq);
      int cqidx = surf2c_qidx[linq];

      if (is_edge_in_dir1) {
        // Velocity value at dir1 boundary
        double vbound = pidx[dir1] == phase_range.lower[dir1] ?
          phase_grid.lower[dir1] : phase_grid.upper[dir1];

        surf_comp_to_phys(dir1, pdim, surf_ord, phase_grid.dx, xc, xmu);
        xmu[dir1] = vbound;

        double rel_speedsq_d = 0.0;
        double rel_vel_in_dir1_d = xmu[dir1] - u_drift_d[tot_conf_quad*d1+cqidx];
        for (int d=0; d<vdim; ++d) {
          rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
        }
        double rel_speed_d = sqrt(rel_speedsq_d);

        pot_q[d1*tot_surf_quad + linq] += eval_fpo_g(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
        pot_deriv_q[d1*tot_surf_quad + linq] += eval_fpo_d2gdv2(den_d[cqidx], rel_vel_in_dir1_d, vtsq_d[cqidx], rel_speed_d);
      }
    }
  }
}

__global__ static void
gkyl_proj_dgdv_on_basis_surf_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, int tot_conf_quad,
  const struct gkyl_array *surf_basis_at_ords, const struct gkyl_array *surf_ordinates,
  const struct gkyl_array *prim_moms_conf_quad, int *surf2c_qidx,
  struct gkyl_array *fpo_dgdv_surf_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_surf_quad = surf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_surf_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *dgdv_surf_q = (double *)gkyl_array_fetch(fpo_dgdv_surf_quad, linp);

    // Compute dG/dv in transverse directions at each boundary surface
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.upper[dir1]);

      const double* surf_ord = (const double*)gkyl_array_cfetch(surf_ordinates, linq);
      int cqidx = surf2c_qidx[linq];

      // Velocity value at dir1 boundary
      double vbound = pidx[dir1] == phase_range.lower[dir1] ?
        phase_grid.lower[dir1] : phase_grid.upper[dir1];

      surf_comp_to_phys(dir1, pdim, surf_ord, phase_grid.dx, xc, xmu);
      xmu[dir1] = vbound;

      double rel_speedsq_d = 0.0;
      for (int d=0; d<vdim; ++d) {
        rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
      }
      double rel_speed_d = sqrt(rel_speedsq_d);

      if (is_edge_in_dir1) {
        for (int d2=0; d2<vdim; ++d2) {
          if (d1 == d2) continue;
          int dir2 = d2+cdim;
          double rel_vel_in_dir2_d = xmu[dir2] - u_drift_d[tot_conf_quad*d2+cqidx];

          int dgdv_index = (d1 < d2) ? d2-1 : d2;
          dgdv_index += (vdim-1)*d1;
          dgdv_surf_q[dgdv_index*tot_surf_quad + linq] += eval_fpo_dgdv(den_d[cqidx],
            rel_vel_in_dir2_d, vtsq_d[cqidx], rel_speedsq_d);
        }
      }
    }
  }
}


void
gkyl_proj_maxwellian_pots_on_basis_advance_cu(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* prim_moms,
  struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
  struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
  struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf,
  struct gkyl_array *fpo_d2gdv2_surf)
{
  // Zero out input arrays
  gkyl_array_clear(fpo_h, 0.0);
  gkyl_array_clear(fpo_g, 0.0);
  gkyl_array_clear(fpo_h_surf, 0.0);
  gkyl_array_clear(fpo_g_surf, 0.0);
  gkyl_array_clear(fpo_dhdv_surf, 0.0);
  gkyl_array_clear(fpo_dgdv_surf, 0.0);
  gkyl_array_clear(fpo_d2gdv2_surf, 0.0);

  // Evaluate primitive moments at configuration space quadrature points and nodes
  dim3 dimGrid_conf, dimBlock_conf;
  int tot_conf_quad = up->tot_conf_quad;
  gkyl_array_clear(up->prim_moms_conf_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, tot_conf_quad);
  gkyl_proj_maxwellian_pots_on_basis_conf_quad_ker<<<dimGrid_conf, dimBlock_conf>>>(*conf_range,
    up->conf_basis_at_ords->on_dev, prim_moms->on_dev, up->prim_moms_conf_quad->on_dev);

  // Phase space quadrature quantities. First compute and project H, then the same for G.
  // Computed individually like this to avoid having to allocate more arrays
  dim3 dimGrid, dimBlock;
  int tot_quad = up->tot_quad;
  gkyl_array_clear(up->pot_phase_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_quad);
  gkyl_proj_h_on_basis_phase_quad_ker<<<dimGrid, dimBlock>>>(up->grid,
    *phase_range, *conf_range, up->ordinates->on_dev, up->conf_basis_at_ords->on_dev,
    up->prim_moms_conf_quad->on_dev, up->p2c_qidx, up->pot_phase_quad->on_dev);

  // Matrix multiplication for nodal to modal conversion
  gkyl_mat_mm_array(up->phase_quad_nodal_to_modal_mem, up->pot_phase_quad, fpo_h);

  // Compute G
  gkyl_array_clear(up->pot_phase_quad, 0.0);
  gkyl_proj_g_on_basis_phase_quad_ker<<<dimGrid, dimBlock>>>(up->grid,
    *phase_range, *conf_range, up->ordinates->on_dev, up->conf_basis_at_ords->on_dev,
    up->prim_moms_conf_quad->on_dev, up->p2c_qidx, up->pot_phase_quad->on_dev);
  gkyl_mat_mm_array(up->phase_quad_nodal_to_modal_mem, up->pot_phase_quad, fpo_g);

  // Phase space surface quadrature quantities.
  // Compute H, dH/dv surface expansions
  dim3 dimGrid_surf, dimBlock_surf;
  int tot_surf_quad = up->tot_surf_quad;
  gkyl_array_clear(up->pot_surf_quad, 0.0);
  gkyl_array_clear(up->pot_deriv_surf_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_surf, &dimBlock_surf, *phase_range, tot_surf_quad);
  gkyl_proj_h_on_basis_surf_quad_ker<<<dimGrid_surf, dimBlock_surf>>>(up->grid,
    *phase_range, *conf_range, up->tot_conf_quad, up->surf_basis_at_ords->on_dev,
    up->surf_ordinates->on_dev, up->prim_moms_conf_quad->on_dev, up->surf2c_qidx,
    up->pot_surf_quad->on_dev, up->pot_deriv_surf_quad->on_dev);

  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_h_surf, 1.0, up->sol_pot_surf_modal);
  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_deriv_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_dhdv_surf, 1.0, up->sol_pot_surf_modal);

  // Compute G, d2G/dv2 surface expansions
  gkyl_array_clear(up->pot_surf_quad, 0.0);
  gkyl_array_clear(up->pot_deriv_surf_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_surf, &dimBlock_surf, *phase_range, tot_surf_quad);
  gkyl_proj_g_on_basis_surf_quad_ker<<<dimGrid_surf, dimBlock_surf>>>(up->grid,
    *phase_range, *conf_range, up->tot_conf_quad, up->surf_basis_at_ords->on_dev,
    up->surf_ordinates->on_dev, up->prim_moms_conf_quad->on_dev, up->surf2c_qidx,
    up->pot_surf_quad->on_dev, up->pot_deriv_surf_quad->on_dev);

  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_g_surf, 1.0, up->sol_pot_surf_modal);
  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_deriv_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_d2gdv2_surf, 1.0, up->sol_pot_surf_modal);

  // Compute dG/dv surface expansions 
  gkyl_array_clear(up->fpo_dgdv_at_surf_ords, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_surf, &dimBlock_surf, *phase_range, tot_surf_quad);
  gkyl_proj_dgdv_on_basis_surf_quad_ker<<<dimGrid_surf, dimBlock_surf>>>(up->grid,
    *phase_range, *conf_range, up->tot_conf_quad, up->surf_basis_at_ords->on_dev,
    up->surf_ordinates->on_dev, up->prim_moms_conf_quad->on_dev, up->surf2c_qidx,
    up->fpo_dgdv_at_surf_ords->on_dev);
  gkyl_mat_mm_array(up->dgdv_surf_quad_nodal_to_modal_mem, up->fpo_dgdv_at_surf_ords, fpo_dgdv_surf);
}

struct gkyl_proj_maxwellian_pots_on_basis*
gkyl_proj_maxwellian_pots_on_basis_cu_ho_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, int num_quad)
{
  gkyl_proj_maxwellian_pots_on_basis *up = (gkyl_proj_maxwellian_pots_on_basis *)
    gkyl_malloc(sizeof(gkyl_proj_maxwellian_pots_on_basis));
  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  int vdim = up->pdim-up->cdim;
  up->num_quad = num_quad;

  up->conf_basis = *conf_basis;

  up->use_gpu = true;

  if (phase_basis->poly_order == 1) {
    gkyl_cart_modal_hybrid(&up->surf_basis, up->cdim, vdim-1);
    up->surf_basis_dev = gkyl_cart_modal_hybrid_cu_dev_new(up->cdim, vdim-1);
  } else {
    gkyl_cart_modal_serendip(&up->surf_basis, up->pdim-1, phase_basis->poly_order);
    up->surf_basis_dev = gkyl_cart_modal_serendip_cu_dev_new(up->pdim-1, phase_basis->poly_order);
  }
  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->num_surf_basis = up->surf_basis.num_basis;
 
  // Quadrature points for phase space expansion
  up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, false);

  // Initialize quadrature for surface expansion
  up->tot_surf_quad = init_quad_values(up->cdim, &up->surf_basis, num_quad,
    &up->surf_ordinates, &up->surf_weights, &up->surf_basis_at_ords, false);

  // Hybrid basis support: uses p=2 in velocity space
  int num_quad_v = num_quad;
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    // Maxwellian potentials are always 3V
    is_vdim_p2[0] = true, is_vdim_p2[1] = true, is_vdim_p2[2] = true;
  }

  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);
  up->surf_qrange = get_qrange(up->cdim, up->pdim-1, num_quad, num_quad_v, is_vdim_p2);

  // Evaluate conf basis at nodes
  struct gkyl_array *conf_nodes_ho = gkyl_array_new(GKYL_DOUBLE,
    up->cdim, up->surf_basis.num_basis);
  up->conf_basis.node_list((double *)gkyl_array_fetch(conf_nodes_ho, 0));

  // Allocate the memory for computing the specific phase and surface nodal to modal calculations
  struct gkyl_mat_mm_array_mem *phase_quad_nodal_to_modal_mem_ho;
  phase_quad_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(up->num_phase_basis,
    up->tot_quad, 1.0, 0.0,
    GKYL_NO_TRANS, GKYL_NO_TRANS, false);

  struct gkyl_mat_mm_array_mem *surf_quad_nodal_to_modal_mem_ho;
  surf_quad_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(vdim*up->num_surf_basis,
    vdim*up->tot_surf_quad, 1.0, 0.0,
    GKYL_NO_TRANS, GKYL_NO_TRANS, false);

  struct gkyl_mat_mm_array_mem *dgdv_surf_quad_nodal_to_modal_mem_ho;
  dgdv_surf_quad_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(2*vdim*up->num_surf_basis, 
    2*vdim*up->tot_surf_quad, 1.0, 0.0,
    GKYL_NO_TRANS, GKYL_NO_TRANS, false);

  // Compute the matrix A for the phase and surface nodal to modal memory
  const double *phase_w = (const double*) up->weights->data;
  const double *phaseb_o = (const double*) up->basis_at_ords->data;
  for (int n=0; n<up->tot_quad; ++n){
    for (int k=0; k<up->num_phase_basis; ++k){
      gkyl_mat_set(phase_quad_nodal_to_modal_mem_ho->A, k, n, phase_w[n]*phaseb_o[k+up->num_phase_basis*n]);
    }
  }

  // Block matrix for vdim component quad nodal to modal matrix multiplication
  const double *surf_w = (const double*) up->surf_weights->data;
  const double *surfb_o = (const double*) up->surf_basis_at_ords->data;
  for (int n=0; n<vdim*up->tot_surf_quad; ++n){
    for (int k=0; k<vdim*up->num_surf_basis; ++k){
      bool block = !((n-n%up->tot_surf_quad)/up->tot_surf_quad - (k-k%up->num_surf_basis)/up->num_surf_basis);
      if (block) {
        gkyl_mat_set(surf_quad_nodal_to_modal_mem_ho->A, k, n,
          surf_w[n%up->tot_surf_quad]*surfb_o[(k%up->num_surf_basis)+up->num_surf_basis*(n%up->tot_surf_quad)]);
      }
    }
  }
  for (int n=0; n<2*vdim*up->tot_surf_quad; ++n){
    for (int k=0; k<2*vdim*up->num_surf_basis; ++k){
      bool block = !((n-n%up->tot_surf_quad)/up->tot_surf_quad - (k-k%up->num_surf_basis)/up->num_surf_basis);
      if (block) {
        gkyl_mat_set(dgdv_surf_quad_nodal_to_modal_mem_ho->A, k, n,
          surf_w[n%up->tot_surf_quad]*surfb_o[(k%up->num_surf_basis)+up->num_surf_basis*(n%up->tot_surf_quad)]);
      }
    }
  }

  // Copy matrix memory to device
  up->phase_quad_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(up->num_phase_basis, up->tot_quad, 1.0, 0.0,
    GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
  gkyl_mat_copy(up->phase_quad_nodal_to_modal_mem->A, phase_quad_nodal_to_modal_mem_ho->A);
  gkyl_mat_mm_array_mem_release(phase_quad_nodal_to_modal_mem_ho);

  up->surf_quad_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(vdim*up->num_surf_basis, vdim*up->tot_surf_quad, 1.0, 0.0,
    GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
  gkyl_mat_copy(up->surf_quad_nodal_to_modal_mem->A, surf_quad_nodal_to_modal_mem_ho->A);
  gkyl_mat_mm_array_mem_release(surf_quad_nodal_to_modal_mem_ho);

  up->dgdv_surf_quad_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(2*vdim*up->num_surf_basis, 2*vdim*up->tot_surf_quad, 1.0, 0.0,
    GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
  gkyl_mat_copy(up->dgdv_surf_quad_nodal_to_modal_mem->A, dgdv_surf_quad_nodal_to_modal_mem_ho->A);
  gkyl_mat_mm_array_mem_release(dgdv_surf_quad_nodal_to_modal_mem_ho);

  // Initialize quadrature data on device
  // Quadrature points for configuration space expansion
  up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, up->use_gpu);

  // Quadrature points for phase space expansion
  up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, up->use_gpu);

  // Initialize quadrature for surface expansion
  up->tot_surf_quad = init_quad_values(up->cdim, &up->surf_basis, num_quad,
    &up->surf_ordinates, &up->surf_weights, &up->surf_basis_at_ords, up->use_gpu);

  // Allocate arrays to store quantities at quadrature points
  up->prim_moms_conf_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    (vdim+2)*up->tot_conf_quad, conf_range->volume);
  up->pot_phase_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    up->tot_quad, phase_range->volume);

  // Primitive moments evaluated at nodes
  up->prim_moms_conf_nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    (vdim+2)*up->num_conf_basis, conf_range->volume);

  // H and G surface quadrature
  up->pot_surf_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    vdim*up->tot_surf_quad, phase_range->volume);
  // dH/dv and d2G/dv2 surface quadrature
  up->pot_deriv_surf_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    vdim*up->tot_surf_quad, phase_range->volume);
  // Nodal evaluation of dG/dv at surfaces in transverse directions
  up->fpo_dgdv_at_surf_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
    2*vdim*up->tot_surf_quad, phase_range->volume);

  // Arrays to store computed potentials for copying
  up->sol_pot_surf_modal = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    vdim*up->num_surf_basis, phase_range->volume);

  // Mappings between phase and configuration space quadrature points
  int p2c_qidx_ho[up->phase_qrange.volume];
  int surf2c_qidx_ho[up->surf_qrange.volume];
  up->p2c_qidx = (int *) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);
  up->surf2c_qidx = (int *) gkyl_cu_malloc(sizeof(int)*up->surf_qrange.volume);

  struct gkyl_range conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);

  int pidx[GKYL_MAX_DIM];
  for (int n=0; n<up->tot_quad; ++n) {
    gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
    int cqidx = gkyl_range_idx(&conf_qrange, pidx);
    p2c_qidx_ho[n] = cqidx;
  }
  gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);

  for (int n=0; n<up->tot_surf_quad; ++n) {
    gkyl_range_inv_idx(&up->surf_qrange, n, pidx);
    int cqidx = gkyl_range_idx(&conf_qrange, pidx);
    surf2c_qidx_ho[n] = cqidx;
  }
  gkyl_cu_memcpy(up->surf2c_qidx, surf2c_qidx_ho, sizeof(int)*up->surf_qrange.volume, GKYL_CU_MEMCPY_H2D);

  return up;
}

