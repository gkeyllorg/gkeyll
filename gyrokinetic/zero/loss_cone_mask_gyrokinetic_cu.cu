/* -*- c++ -*- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>

#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
}

__global__ static void
gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu_ker(struct gkyl_range conf_range,
  const struct gkyl_array* basis_at_ords_conf, const struct gkyl_array* bmag, const double *bmag_max,
  struct gkyl_array* Dbmag_quad_d)
{    
  int num_basis_conf = basis_at_ords_conf->ncomp;
  int tot_quad_conf = basis_at_ords_conf->size;

  int cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);
    long linidx = gkyl_range_idx(&conf_range, cidx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(bmag, linidx);

    double *bmag_quad = (double*) gkyl_array_fetch(Dbmag_quad_d, linidx);

    for (int n=0; n<tot_quad_conf; ++n) {
      const double *b_ord = (const double*) gkyl_array_cfetch(basis_at_ords_conf, n);

      for (int k=0; k<num_basis_conf; ++k)
        bmag_quad[n] += bmag_d[k]*b_ord[k];

      bmag_quad[n] = bmag_max[0] - bmag_quad[n];
    }
  }
}

void 
gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const double *bmag_max)
{
  int nblocks = conf_range->nblocks, nthreads = conf_range->nthreads;
  gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu_ker<<<nblocks, nthreads>>>(*conf_range, 
    up->basis_at_ords_conf->on_dev, bmag->on_dev, bmag_max, up->Dbmag_quad->on_dev);
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
gkyl_loss_cone_mask_gyrokinetic_qDphiDbmag_quad_ker(struct gkyl_range conf_range, 
  const struct gkyl_array* basis_at_ords_conf, double charge, const struct gkyl_array* phi,
  const double *phi_m, const struct gkyl_array* Dbmag_quad, struct gkyl_array* qDphiDbmag_quad)
{
  int num_basis_conf = basis_at_ords_conf->ncomp;

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linc2 goes from 0 to tot_quad_conf= basis_at_ords_conf->size.
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long linidx = gkyl_range_idx(&conf_range, cidx);

    const double *phi_d = (const double*) gkyl_array_cfetch(phi, linidx);
    const double *Dbmag_quad_d = (const double*) gkyl_array_cfetch(Dbmag_quad, linidx);

    // Sum over basis at configuration-space quadrature points. 
    const double *b_ord = (const double*) gkyl_array_cfetch(basis_at_ords_conf, linc2);
    double phi_quad = 0;
    for (int k=0; k<num_basis_conf; ++k)
      phi_quad += phi_d[k]*b_ord[k];

    // Potential energy term at each quadrature point.
    double *qDphiDbmag_quad_d = (double*) gkyl_array_fetch(qDphiDbmag_quad, linidx);
    if (Dbmag_quad_d[linc2] > 0.0)
      qDphiDbmag_quad_d[linc2] = charge*(phi_quad-phi_m[0])/Dbmag_quad_d[linc2];
    else
      qDphiDbmag_quad_d[linc2] = 0.0;
  }
}

__global__ static void
gkyl_loss_cone_mask_gyrokinetic_ker(struct gkyl_rect_grid grid_phase,
  struct gkyl_range phase_range, struct gkyl_range conf_range, struct gkyl_range vel_range,
  double mass, const struct gkyl_array* phase_ordinates, 
  const double *bmag_max_loc, const struct gkyl_array* qDphiDbmag_quad, const struct gkyl_array* Dbmag_quad,
  const int *p2c_qidx, struct gkyl_array* vmap, struct gkyl_basis* vmap_basis, struct gkyl_array* mask_out)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.0};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  int tot_phase_quad = phase_ordinates->size;

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // Get configuration-space linear index.
    for (unsigned int k = 0; k < cdim; k++) cidx[k] = pidx[k];
    long linidx_conf = gkyl_range_idx(&conf_range, cidx);

    const double *Dbmag_quad_d = (const double*) gkyl_array_cfetch(Dbmag_quad, linidx_conf);
    const double *qDphiDbmag_quad_d = (const double*) gkyl_array_cfetch(qDphiDbmag_quad, linidx_conf);

    gkyl_rect_grid_cell_center(&grid_phase, pidx, xc);
    long linidx_phase = gkyl_range_idx(&phase_range, pidx);
    double *mask_d = (double*) gkyl_array_fetch(mask_out, linidx_phase);

    for (int d = cdim; d < pdim; d++) vidx[d-cdim] = pidx[d];
    long linidx_vel = gkyl_range_idx(&vel_range, vidx);
    const double *vmap_d = (const double*) gkyl_array_cfetch(vmap, linidx_vel);

    mask_d[0] = 1.0; // In this case the mask has ncomp=1.

    for (int n=0; n<tot_phase_quad; ++n) {
      int cqidx = p2c_qidx[n];

      const double *xcomp_d = (const double*) gkyl_array_cfetch(phase_ordinates, n);

      // Convert comp position coordinate to phys pos coord.
      log_to_comp(cdim, xcomp_d, grid_phase.dx, xc, xmu);
  
      // Convert comp velocity coordinate to phys velocity coord.
      double xcomp[1];
      for (int vd = 0; vd < vdim; vd++) {
        xcomp[0] = xcomp_d[cdim+vd];
        xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
      }
  
      // KEparDbmag = 0.5*mass*pow(vpar,2)/(bmag_max-bmag[0]).
      double KEparDbmag = 0.0;
      if (Dbmag_quad_d[cqidx] > 0.0)
        KEparDbmag = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_d[cqidx];
      else
        KEparDbmag = 0.0;
  
      double mu_bound = GKYL_MAX2(0.0, KEparDbmag+qDphiDbmag_quad_d[cqidx]);
  
      if ( !(mu_bound < xmu[cdim+1] && fabs(xmu[cdim-1]) < fabs(bmag_max_loc[cdim-1])) ) {
        mask_d[0] = 0.0;
        break;
      }
    }
  }
}

__global__ static void
gkyl_loss_cone_mask_gyrokinetic_quad_ker(struct gkyl_rect_grid grid_phase,
  struct gkyl_range phase_range, struct gkyl_range conf_range, struct gkyl_range vel_range,
  double mass, double norm_fac, const struct gkyl_array* phase_ordinates, 
  const double *bmag_max_loc, const struct gkyl_array* qDphiDbmag_quad, const struct gkyl_array* Dbmag_quad,
  const int *p2c_qidx, struct gkyl_array* vmap, struct gkyl_basis* vmap_basis, struct gkyl_array* mask_out_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.0};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  // 2D thread grid
  // linc2 goes from 0 to tot_quad_phase
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // Get configuration-space linear index.
    for (unsigned int k = 0; k < cdim; k++) cidx[k] = pidx[k];

    long linidx_conf = gkyl_range_idx(&conf_range, cidx);

    const double *Dbmag_quad_d = (const double*) gkyl_array_cfetch(Dbmag_quad, linidx_conf);
    const double *qDphiDbmag_quad_d = (const double*) gkyl_array_cfetch(qDphiDbmag_quad, linidx_conf);

    gkyl_rect_grid_cell_center(&grid_phase, pidx, xc);
    long linidx_phase = gkyl_range_idx(&phase_range, pidx);

    int cqidx = p2c_qidx[linc2];
    for (int d = cdim; d < pdim; d++) vidx[d-cdim] = pidx[d];

    long linidx_vel = gkyl_range_idx(&vel_range, vidx);
    const double *vmap_d = (const double*) gkyl_array_cfetch(vmap, linidx_vel);
    const double *xcomp_d = (const double*) gkyl_array_cfetch(phase_ordinates, linc2);

    // Convert comp position coordinate to phys pos coord.
    gkyl_rect_grid_cell_center(&grid_phase, pidx, xc);
    log_to_comp(cdim, xcomp_d, grid_phase.dx, xc, xmu);

    // Convert comp velocity coordinate to phys velocity coord.
    double xcomp[1];
    for (int vd = 0; vd < vdim; vd++) {
      xcomp[0] = xcomp_d[cdim+vd];
      xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
    }

    // KEparDbmag = 0.5*mass*pow(vpar,2)/(bmag_max-bmag[0]).
    double KEparDbmag = 0.0;
    if (Dbmag_quad_d[cqidx] > 0.0)
      KEparDbmag = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_d[cqidx];
    else
      KEparDbmag = 0.0;

    double mu_bound = GKYL_MAX2(0.0, KEparDbmag+qDphiDbmag_quad_d[cqidx]);

    double *fq = (double*) gkyl_array_fetch(mask_out_quad, linidx_phase);
    if (mu_bound < xmu[cdim+1] && fabs(xmu[cdim-1]) < fabs(bmag_max_loc[cdim-1])) 
      fq[linc2] = norm_fac;
    else
      fq[linc2] = 0.0;
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_advance_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *phi, const double *phi_m, struct gkyl_array *mask_out)
{
  dim3 dimGrid_conf, dimBlock_conf;
  int tot_quad_conf = up->basis_at_ords_conf->size;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, tot_quad_conf);

  gkyl_loss_cone_mask_gyrokinetic_qDphiDbmag_quad_ker<<<dimGrid_conf, dimBlock_conf>>>(*conf_range, 
    up->basis_at_ords_conf->on_dev, up->charge, phi->on_dev, phi_m, up->Dbmag_quad->on_dev,
    up->qDphiDbmag_quad->on_dev);

  const struct gkyl_velocity_map *gvm = up->vel_map;

  if (up->cellwise_trap_loss) {
    // Don't do quadrature.
    int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;
    gkyl_loss_cone_mask_gyrokinetic_ker<<<nblocks, nthreads>>>(*up->grid_phase, *phase_range, *conf_range,
      gvm->local_ext_vel, up->mass, up->ordinates_phase->on_dev,
      up->bmag_max_loc, up->qDphiDbmag_quad->on_dev, up->Dbmag_quad->on_dev, up->p2c_qidx, gvm->vmap->on_dev,
      gvm->vmap_basis, mask_out->on_dev);
  }
  else {
    // Use quadrature.
    dim3 dimGrid, dimBlock;
    int tot_quad_phase = up->basis_at_ords_phase->size;
    gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_quad_phase);

    gkyl_loss_cone_mask_gyrokinetic_quad_ker<<<dimGrid, dimBlock>>>(*up->grid_phase, *phase_range, *conf_range,
      gvm->local_ext_vel, up->mass, up->norm_fac, up->ordinates_phase->on_dev,
      up->bmag_max_loc, up->qDphiDbmag_quad->on_dev, up->Dbmag_quad->on_dev, up->p2c_qidx, gvm->vmap->on_dev,
      gvm->vmap_basis, up->mask_out_quad->on_dev);

    // Call cublas to do the matrix multiplication nodal to modal conversion
    gkyl_mat_mm_array(up->phase_nodal_to_modal_mem, up->mask_out_quad, mask_out);
  }
}
