/* -*- c++ -*- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_lte_proj_on_basis_priv.h>
#include <gkyl_range.h>

#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
}

__global__ static void
gkyl_vlasov_lte_proj_on_basis_geom_quad_vars_cu_ker(struct gkyl_range conf_range,
  const struct gkyl_array* conf_basis_at_ords, int vdim, const struct gkyl_array* h_ij,
  const struct gkyl_array* h_ij_inv, const struct gkyl_array* det_h, 
  struct gkyl_array* h_ij_quad_d, struct gkyl_array* h_ij_inv_quad_d, struct gkyl_array* det_h_quad_d)
{    
  int num_conf_basis = conf_basis_at_ords->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;

  int cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);
    long lincC = gkyl_range_idx(&conf_range, cidx);
    const double *h_ij_d = (const double*) gkyl_array_cfetch(h_ij, lincC);
    const double *h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv, lincC);
    const double *det_h_d = (const double*) gkyl_array_cfetch(det_h, lincC);

    double *h_ij_quad = (double*) gkyl_array_fetch(h_ij_quad_d, lincC);
    double *h_ij_inv_quad = (double*) gkyl_array_fetch(h_ij_inv_quad_d, lincC);
    double *det_h_quad = (double*) gkyl_array_fetch(det_h_quad_d, lincC);

    // Sum over basis for the geometric quantities at configuration-space quadrature points. 
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double*) gkyl_array_cfetch(conf_basis_at_ords, n);

      for (int k=0; k<num_conf_basis; ++k) {
        det_h_quad[n] += det_h_d[k]*b_ord[k];
        for (int j=0; j<vdim*(vdim+1)/2; ++j) {
          h_ij_quad[tot_conf_quad*j + n] += h_ij_d[num_conf_basis*j+k]*b_ord[k];
          h_ij_inv_quad[tot_conf_quad*j + n] += h_ij_inv_d[num_conf_basis*j+k]*b_ord[k];
        }
      }
    }
  }
}

void 
gkyl_vlasov_lte_proj_on_basis_geom_quad_vars_cu(gkyl_vlasov_lte_proj_on_basis *up,
  const struct gkyl_range *conf_range, const struct gkyl_array *h_ij,
  const struct gkyl_array *h_ij_inv, const struct gkyl_array *det_h)
{
  int vdim = up->pdim - up->cdim;
  int nblocks = conf_range->nblocks, nthreads = conf_range->nthreads;
  gkyl_vlasov_lte_proj_on_basis_geom_quad_vars_cu_ker<<<nblocks, nthreads>>>(*conf_range, 
    up->conf_basis_at_ords->on_dev, vdim, h_ij->on_dev,
    h_ij_inv->on_dev, det_h->on_dev, up->h_ij_quad->on_dev,
    up->h_ij_inv_quad->on_dev, up->det_h_quad->on_dev);
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
gkyl_vlasov_lte_proj_on_basis_moms_lte_quad_ker(struct gkyl_range conf_range, int vdim, 
  const struct gkyl_array* conf_basis_at_ords, 
  const struct gkyl_array* moms_lte, const struct gkyl_array* det_h_quad, 
  bool is_relativistic, bool is_canonical_pb, 
  struct gkyl_array* moms_lte_quad, struct gkyl_array* expamp_quad)
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

    const double *moms_lte_d = (const double*) gkyl_array_cfetch(moms_lte, lincC);

    double *moms_lte_quad_d = (double*) gkyl_array_fetch(moms_lte_quad, lincC);

    // Sum over basis for given LTE moments (n, V_drift, T/m) in the stationary frame
    // at configuration-space quadrature points. 
    const double *b_ord = (const double*) gkyl_array_cfetch(conf_basis_at_ords, linc2);
    for (int k=0; k<num_conf_basis; ++k) {
      for (int d=0; d<vdim+2; ++d) {
        moms_lte_quad_d[tot_conf_quad*d+linc2] += moms_lte_d[num_conf_basis*d+k]*b_ord[k];
      }
    }
    // Amplitude of the exponential at each quadrature point.
    const double *n_quad = moms_lte_quad_d;
    const double *T_over_m_quad = &moms_lte_quad_d[tot_conf_quad*(vdim+1)]; 
    double *expamp_quad_d = (double*) gkyl_array_fetch(expamp_quad, lincC);
    if ((n_quad[linc2] > 0.0) && (T_over_m_quad[linc2] > 0.0)) {
      if (is_relativistic) {
        expamp_quad_d[linc2] = n_quad[linc2]*(1.0/(4.0*GKYL_PI*T_over_m_quad[linc2]))*(sqrt(2.0*T_over_m_quad[linc2]/GKYL_PI));
      }
      else if (is_canonical_pb) { 
        const double *det_h_quad_d = (const double*) gkyl_array_cfetch(det_h_quad, lincC);
        expamp_quad_d[linc2] = (1.0/det_h_quad_d[linc2])*n_quad[linc2]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[linc2], vdim));
      }
      else {
        expamp_quad_d[linc2] = n_quad[linc2]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[linc2], vdim));
      }
    }
    else {
      expamp_quad_d[linc2] = 0.0;
    }  
  }
}

__global__ static void
gkyl_vlasov_lte_proj_on_basis_f_lte_quad_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, 
  const struct gkyl_array* conf_basis_at_ords, const struct gkyl_array* phase_ordinates, 
  const struct gkyl_array* moms_lte_quad, const struct gkyl_array* expamp_quad, 
  const struct gkyl_array* h_ij_inv_quad, 
  const int *p2c_qidx, bool is_relativistic, bool is_canonical_pb, 
  struct gkyl_array* f_lte_quad)
{
  double f_floor = 1.0e-40;
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;
  int tot_conf_quad = conf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

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

    // Fetch V_drift and T/m; density dependence already included in expamp_quad
    const double *moms_lte_quad_d = (const double*) gkyl_array_cfetch(moms_lte_quad, lincC);
    const double *V_drift_quad = &moms_lte_quad_d[tot_conf_quad];
    const double *T_over_m_quad = &moms_lte_quad_d[tot_conf_quad*(vdim+1)];

    const double *expamp_quad_d = (const double*) gkyl_array_cfetch(expamp_quad, lincC);

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long lidx = gkyl_range_idx(&phase_range, pidx);

    // Select for a phase space index fq
    double *fq = (double*) gkyl_array_fetch(f_lte_quad, lidx);

    int cqidx = p2c_qidx[linc2];
    comp_to_phys(pdim, (const double*) gkyl_array_cfetch(phase_ordinates, linc2),
      phase_grid.dx, xc, &xmu[0]);

    fq[linc2] = f_floor;
    if (T_over_m_quad[cqidx] > 0.0) {
      if (is_relativistic) {
        double vv = 0.0;
        double vu = 0.0;
        double uu = 0.0;
        // V_drift_quad is the spatial component of the four-velocity u_i = GammaV*V_drift
        for (int d=0; d<vdim; ++d) {
          vv += (V_drift_quad[tot_conf_quad*d + cqidx]*V_drift_quad[tot_conf_quad*d + cqidx]);
          vu += (V_drift_quad[tot_conf_quad*d + cqidx]*xmu[cdim+d]);
          uu += (xmu[cdim+d]*xmu[cdim+d]);
        }
        double GammaV_quad = sqrt(1.0 + vv);
        fq[linc2] += expamp_quad_d[cqidx]*exp((1.0/T_over_m_quad[cqidx]) 
          - (1.0/T_over_m_quad[cqidx])*(GammaV_quad*sqrt(1 + uu) - vu));
      }
      else if (is_canonical_pb) {
        // Assumes a (particle) hamiltonian in canocial form: g = 1/2 g^{ij} w_i_w_j
        const double *h_ij_inv_quad_d = (const double*) gkyl_array_cfetch(h_ij_inv_quad, lincC);
        double efact = 0.0;
        for (int d0=0; d0<vdim; ++d0) {
          for (int d1=d0; d1<vdim; ++d1) {
            int sym_tensor_index = (d0*(2*vdim - d0 + 1))/2 + (d1-d0);
            // Grab the spatial metric component, the ctx includes geometry that isn't 
            // part of the canonical set of variables, like R on the surf of a sphere
            // q_can includes the canonical variables list
            double h_ij_inv_loc = h_ij_inv_quad_d[tot_conf_quad*sym_tensor_index + cqidx]; 
            // For off-diagonal components, we need to count these twice, due to symmetry
            int sym_fact = (d0 == d1) ? 1 : 2;
            efact += sym_fact*h_ij_inv_loc*(xmu[cdim+d0]-V_drift_quad[tot_conf_quad*d0 + cqidx])*(xmu[cdim+d1]-V_drift_quad[tot_conf_quad*d1 + cqidx]);
          }
        }
        fq[linc2] += expamp_quad_d[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
      }
      else {
        double efact = 0.0;        
        for (int d=0; d<vdim; ++d) {
          efact += (xmu[cdim+d]-V_drift_quad[tot_conf_quad*d + cqidx])*(xmu[cdim+d]-V_drift_quad[tot_conf_quad*d + cqidx]);
        }
        fq[linc2] += expamp_quad_d[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
      }
    }
  }
}

void
gkyl_vlasov_lte_proj_on_basis_advance_cu(gkyl_vlasov_lte_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_lte, struct gkyl_array *f_lte)
{
  int vdim = up->pdim - up->cdim;

  gkyl_array_clear(up->moms_lte_quad, 0.0); 
  dim3 dimGrid_conf, dimBlock_conf;
  int tot_conf_quad = up->conf_basis_at_ords->size;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, tot_conf_quad);
  gkyl_vlasov_lte_proj_on_basis_moms_lte_quad_ker<<<dimGrid_conf, dimBlock_conf>>>(*conf_range, 
    vdim, up->conf_basis_at_ords->on_dev, 
    moms_lte->on_dev, 
    up->is_canonical_pb ? up->det_h_quad->on_dev : 0, 
    up->is_relativistic, up->is_canonical_pb, 
    up->moms_lte_quad->on_dev, up->expamp_quad->on_dev);

  dim3 dimGrid, dimBlock;
  int tot_phase_quad = up->basis_at_ords->size;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_phase_quad);
  gkyl_vlasov_lte_proj_on_basis_f_lte_quad_ker<<<dimGrid, dimBlock>>>(up->phase_grid, 
    *phase_range, *conf_range, 
    up->conf_basis_at_ords->on_dev, up->ordinates->on_dev,
    up->moms_lte_quad->on_dev, up->expamp_quad->on_dev, 
    up->is_canonical_pb ? up->h_ij_inv_quad->on_dev : 0, 
    up->p2c_qidx, up->is_relativistic, up->is_canonical_pb, 
    up->f_lte_quad->on_dev);

  // Call cublas to do the matrix multiplication nodal to modal conversion
  gkyl_mat_mm_array(up->phase_nodal_to_modal_mem, up->f_lte_quad, f_lte);

  // Correct the density of the projected LTE distribution function through rescaling.
  // This correction is needed especially for the relativistic LTE, whose pre-factor
  // we construct through an expansion of the Bessel functions to avoid finite 
  // precision effects in such a way that we can recover arbitrary temperature  
  // relativistic LTE distributions by rescaling the distribution to the desired density.  
  gkyl_vlasov_lte_density_moment_advance(up->moments_up, phase_range, conf_range, 
    f_lte, up->num_ratio);

  // compute number density ratio: num_ratio = n/n0
  // 0th component of moms_target is the target density
  gkyl_dg_div_op_range(up->mem, &up->conf_basis, 0, up->num_ratio,
    0, moms_lte, 0, up->num_ratio, conf_range);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis,
    f_lte, up->num_ratio, f_lte, conf_range, phase_range);  
}
