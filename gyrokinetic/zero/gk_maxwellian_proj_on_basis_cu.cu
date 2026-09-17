/* -*- c++ -*- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_gk_maxwellian_proj_on_basis_priv.h>
#include <gkyl_range.h>

#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
}

__global__ static void
gkyl_gk_maxwellian_proj_on_basis_geom_quad_vars_cu_ker(struct gkyl_range conf_range,
  const struct gkyl_array* conf_basis_at_ords, 
  const struct gkyl_array* bmag, const struct gkyl_array* jacobtot, 
  struct gkyl_array* bmag_quad_d, struct gkyl_array* jacobtot_quad_d)
{    
  int num_conf_basis = conf_basis_at_ords->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;

  int cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);
    long lincC = gkyl_range_idx(&conf_range, cidx);
    const double *bmag_d = (const double*) gkyl_array_cfetch(bmag, lincC);
    const double *jacobtot_d = (const double*) gkyl_array_cfetch(jacobtot, lincC);

    double *bmag_quad = (double*) gkyl_array_fetch(bmag_quad_d, lincC);
    double *jacobtot_quad = (double*) gkyl_array_fetch(jacobtot_quad_d, lincC);

    // Sum over basis for the geometric quantities at configuration-space quadrature points. 
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double*) gkyl_array_cfetch(conf_basis_at_ords, n);

      for (int k=0; k<num_conf_basis; ++k) {
        bmag_quad[n] += bmag_d[k]*b_ord[k];
        jacobtot_quad[n] += jacobtot_d[k]*b_ord[k];
      }
    }
  }
}

void 
gkyl_gk_maxwellian_proj_on_basis_geom_quad_vars_cu(gkyl_gk_maxwellian_proj_on_basis *up,
  const struct gkyl_range *conf_range, 
  const struct gkyl_array *bmag, const struct gkyl_array *jacobtot)
{
  int nblocks = conf_range->nblocks, nthreads = conf_range->nthreads;
  gkyl_gk_maxwellian_proj_on_basis_geom_quad_vars_cu_ker<<<nblocks, nthreads>>>(*conf_range, 
    up->conf_basis_at_ords->on_dev, 
    bmag->on_dev, jacobtot->on_dev, 
    up->bmag_quad->on_dev, up->jacobtot_quad->on_dev);
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
gkyl_gk_maxwellian_proj_on_basis_moms_quad_ker(struct gkyl_range conf_range, 
  int vdim_phys, int num_comp, bool bimaxwellian, bool use_jacobtot, 
  const struct gkyl_array* conf_basis_at_ords, 
  const struct gkyl_array* moms_maxwellian, 
  const struct gkyl_array* bmag_quad, const struct gkyl_array* jacobtot_quad, 
  struct gkyl_array* moms_maxwellian_quad, struct gkyl_array* expamp_quad)
{
  int num_conf_basis = conf_basis_at_ords->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;  

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linc2 goes from 0 to tot_conf_quad
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long lincC = gkyl_range_idx(&conf_range, cidx);

    const double *moms_maxwellian_d = (const double*) gkyl_array_cfetch(moms_maxwellian, lincC);
    const double *bmag_quad_d = (const double*) gkyl_array_cfetch(bmag_quad, lincC);
    const double *jacobtot_quad_d = (const double*) gkyl_array_cfetch(jacobtot_quad, lincC);

    double *moms_maxwellian_quad_d = (double*) gkyl_array_fetch(moms_maxwellian_quad, lincC);

    // Sum over basis for given Maxwellian moments (n, upar, T/m) or (n, upar, Tpar/m, Tperp/m) 
    // at configuration-space quadrature points. 
    const double *b_ord = (const double*) gkyl_array_cfetch(conf_basis_at_ords, linc2);
    for (int k=0; k<num_conf_basis; ++k) {
      for (int d=0; d<num_comp; ++d) {
        moms_maxwellian_quad_d[tot_conf_quad*d+linc2] += moms_maxwellian_d[num_conf_basis*d+k]*b_ord[k];
      }
    }

    // Amplitude of the exponential at each quadrature point.
    const double *n_quad = moms_maxwellian_quad_d;
    const double *T_over_m_quad = &moms_maxwellian_quad_d[tot_conf_quad*2]; 
    double *expamp_quad_d = (double*) gkyl_array_fetch(expamp_quad, lincC);
    if ((n_quad[linc2] > 0.0) && (T_over_m_quad[linc2] > 0.0)) {
      if (bimaxwellian) {
        const double *Tperp_over_m_quad = &moms_maxwellian_quad_d[tot_conf_quad*3];
        expamp_quad_d[linc2] = n_quad[linc2]/(sqrt(pow(2.0*GKYL_PI, 3.0)*T_over_m_quad[linc2])*Tperp_over_m_quad[linc2]);
      }
      else {
        expamp_quad_d[linc2] = n_quad[linc2]/(sqrt(pow(2.0*GKYL_PI*T_over_m_quad[linc2], vdim_phys)));
      }      
    }
    else {
      expamp_quad_d[linc2] = 0.0;
    }

    // Scale amplitude of the exponential by desired Jacobian factor 
    // Either the total Jacobian or just the velocity-space Jacobian bmag
    if (use_jacobtot) {
      expamp_quad_d[linc2] *= jacobtot_quad_d[linc2];
    }  
    else {
      expamp_quad_d[linc2] *= bmag_quad_d[linc2];
    }    
  }
}

__global__ static void
gkyl_gk_maxwellian_proj_on_basis_f_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, struct gkyl_range vel_range,
  bool bimaxwellian, double mass, 
  const struct gkyl_array* conf_basis_at_ords, const struct gkyl_array* phase_ordinates, 
  const struct gkyl_array* moms_maxwellian_quad, const struct gkyl_array* expamp_quad, 
  const struct gkyl_array* bmag_quad, const int *p2c_qidx, 
  struct gkyl_array* vmap, struct gkyl_array* jacobvel, struct gkyl_basis* vmap_basis, 
  struct gkyl_array* f_maxwellian_quad)
{
  double f_floor = 1.0e-40;
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;
  int tot_conf_quad = conf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.0};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  // 2D thread grid
  // linc2 goes from 0 to tot_phase_quad
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // get configuration-space linear index.
    for (unsigned int k = 0; k < conf_range.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_range, cidx);

    // Fetch upar and T/m; density dependence already included in expamp_quad
    const double *moms_maxwellian_quad_d = (const double*) gkyl_array_cfetch(moms_maxwellian_quad, lincC);
    const double *upar_quad = &moms_maxwellian_quad_d[tot_conf_quad];
    const double *T_over_m_quad = &moms_maxwellian_quad_d[tot_conf_quad*2];

    const double *bmag_quad_d = (const double*) gkyl_array_cfetch(bmag_quad, lincC);

    const double *expamp_quad_d = (const double*) gkyl_array_cfetch(expamp_quad, lincC);

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long lidx = gkyl_range_idx(&phase_range, pidx);

    // Select for a phase space index fq
    double *fq = (double*) gkyl_array_fetch(f_maxwellian_quad, lidx);

    int cqidx = p2c_qidx[linc2];
    for (int d = cdim; d < pdim; d++) {
      vidx[d-cdim] = pidx[d];
    }
    long vlinidx = gkyl_range_idx(&vel_range, vidx);
    const double *vmap_d = (const double*) gkyl_array_cfetch(vmap, vlinidx);
    const double *xcomp_d = (const double*) gkyl_array_cfetch(phase_ordinates, linc2);
    // Convert comp velocity coordinate to phys velocity coord.
    double xcomp[1];
    for (int vd = 0; vd < vdim; vd++) {
      xcomp[0] = xcomp_d[cdim+vd];
      xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
    }
    // Fetch velocity space Jacobian for scaling distribution function
    const double *jacobvel_d = (const double*) gkyl_array_cfetch(jacobvel, lidx);

    double efact = 0.0;
    // vpar term.
    efact += pow(xmu[cdim]-upar_quad[cqidx],2)/(2.0*T_over_m_quad[cqidx]);
    // mu term (only for 2v).
    if (bimaxwellian) {
      const double *Tperp_over_m_quad = &moms_maxwellian_quad_d[tot_conf_quad*3];
      efact += xmu[cdim+1]*bmag_quad_d[cqidx]/(mass*Tperp_over_m_quad[cqidx]);
    }
    else {
      efact += xmu[cdim+1]*bmag_quad_d[cqidx]/(mass*T_over_m_quad[cqidx]);
    }
    fq[linc2] = T_over_m_quad[cqidx] > 0.0 ? f_floor + 
      jacobvel_d[0]*expamp_quad_d[cqidx]*exp(-efact) : f_floor; 
  }
}

void
gkyl_gk_maxwellian_proj_on_basis_advance_cu(gkyl_gk_maxwellian_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_maxwellian, bool use_jacobtot, 
  struct gkyl_array *f_maxwellian)
{
  int vdim = up->pdim - up->cdim;
  int vdim_phys = vdim == 1 ? 1 : 3;

  gkyl_array_clear(up->moms_maxwellian_quad, 0.0); 
  dim3 dimGrid_conf, dimBlock_conf;
  int tot_conf_quad = up->conf_basis_at_ords->size;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, tot_conf_quad);
  gkyl_gk_maxwellian_proj_on_basis_moms_quad_ker<<<dimGrid_conf, dimBlock_conf>>>(*conf_range, 
    vdim_phys, up->num_comp, up->bimaxwellian, use_jacobtot, up->conf_basis_at_ords->on_dev, 
    moms_maxwellian->on_dev, 
    up->bmag_quad->on_dev, up->jacobtot_quad->on_dev, 
    up->moms_maxwellian_quad->on_dev, up->expamp_quad->on_dev);


  const struct gkyl_velocity_map *gvm = up->vel_map;
  dim3 dimGrid, dimBlock;
  int tot_phase_quad = up->basis_at_ords->size;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_phase_quad);
  gkyl_gk_maxwellian_proj_on_basis_f_quad_ker<<<dimGrid, dimBlock>>>(up->phase_grid, 
    *phase_range, *conf_range, gvm->local_ext_vel, 
    up->bimaxwellian, up->mass, 
    up->conf_basis_at_ords->on_dev, up->ordinates->on_dev,
    up->moms_maxwellian_quad->on_dev, up->expamp_quad->on_dev, 
    up->bmag_quad->on_dev, 
    up->p2c_qidx, gvm->vmap->on_dev, gvm->jacobvel->on_dev, gvm->vmap_basis, 
    up->f_maxwellian_quad->on_dev);

  // Call cublas to do the matrix multiplication nodal to modal conversion
  gkyl_mat_mm_array(up->phase_nodal_to_modal_mem, up->f_maxwellian_quad, f_maxwellian);

  // Correct the density of the projected Maxwellian (or bi-Maxwellian) 
  // distribution function through rescaling.  
  gkyl_gk_maxwellian_density_moment_advance(up->moments_up, phase_range, conf_range, 
    f_maxwellian, up->num_ratio);

  // compute number density ratio: num_ratio = n/n0
  // 0th component of moms_target is the target density
  gkyl_dg_div_op_range(up->mem, &up->conf_basis, 0, up->num_ratio,
    0, moms_maxwellian, 0, up->num_ratio, conf_range);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis,
    f_maxwellian, up->num_ratio, f_maxwellian, conf_range, phase_range);  
}
