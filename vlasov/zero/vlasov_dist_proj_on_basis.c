#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mat.h>
#include <gkyl_vlasov_dist_proj_on_basis.h>
#include <gkyl_vlasov_dist_proj_on_basis_priv.h>
#include <gkyl_range.h>
#include <assert.h>

// create range to loop over quadrature points.
static inline struct gkyl_range
get_qrange(int cdim, int dim, int num_quad, int num_quad_v, bool *is_vdim_p2)
{
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  for (int i=cdim; i<dim; ++i) qshape[i] = is_vdim_p2[i-cdim] ? num_quad_v : num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

// Sets ordinates, weights and basis functions at ords.
// Returns the total number of quadrature nodes
static int
init_quad_values(int cdim, const struct gkyl_basis *basis,
  enum gkyl_quad_type quad_type, int num_quad,
  struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  int vdim = ndim-cdim;
  int num_quad_v = num_quad;
  // hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    for (int d=0; d<vdim; d++) {
      is_vdim_p2[d] = true;
    }
  }

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];

  if (quad_type == GKYL_POSITIVITY_QUAD) {
    // Positivity nodes are -1/3, 1/3 and have weights 1.0. All other
    // weights are set to zero.
    for (int i=0; i<num_quad; ++i) {
      ordinates1[i] = 0.0;
      weights1[i] = 0.0;
    }
    for (int i=0; i<num_quad_v; ++i) {
      ordinates1_v[i] = 0.0;
      weights1_v[i] = 0.0;
    }
    ordinates1[0] = -1.0/3.0;
    ordinates1[1] = 1.0/3.0;
    weights1[0] = 1.0;
    weights1[1] = 1.0;

    ordinates1_v[0] = -1.0/3.0;
    ordinates1_v[1] = 0.0;
    ordinates1_v[2] = 1.0/3.0;
    weights1_v[0] = 3.0;
    weights1_v[1] = -4.0;
    weights1_v[2] = 3.0;
  }
  else if (quad_type == GKYL_GAUSS_QUAD) {
    if (num_quad <= gkyl_gauss_max) {
      // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
      memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
      memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
    } 
    else {
      gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
    }
    if (num_quad_v <= gkyl_gauss_max) {
      memcpy(ordinates1_v, gkyl_gauss_ordinates[num_quad_v], sizeof(double[num_quad_v]));
      memcpy(weights1_v, gkyl_gauss_weights[num_quad_v], sizeof(double[num_quad_v]));
    } 
    else {
      gkyl_gauleg(-1, 1, ordinates1_v, weights1_v, num_quad_v);
    }
  }

  struct gkyl_range qrange = get_qrange(cdim, ndim, num_quad, num_quad_v, is_vdim_p2);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *ordinates = gkyl_array_cu_dev_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  } 
  else {
    *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<cdim; ++i) {
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    for (int i=cdim; i<ndim; ++i) {
      ord[i] = is_vdim_p2[i-cdim] ? 
        ordinates1_v[iter.idx[i]-qrange.lower[i]] : ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    
    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<cdim; ++i) {
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
    }
    for (int i=cdim; i<ndim; ++i) {
      wgt[0] *= is_vdim_p2[i-cdim] ? 
        weights1_v[iter.idx[i]-qrange.lower[i]] : weights1[iter.idx[i]-qrange.lower[i]];
    }
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu) {
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  }
  else {
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  }

  for (int n=0; n<tot_quad; ++n) {
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));
  }

  // copy host array to device array
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

static void
gkyl_vlasov_dist_proj_on_basis_geom_quad_vars(gkyl_vlasov_dist_proj_on_basis *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array *h_ij,
  const struct gkyl_array *h_ij_inv, const struct gkyl_array *det_h)
{
// Setup the intial geometric vars, on GPU 
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_vlasov_dist_proj_on_basis_geom_quad_vars_cu(up, conf_range, h_ij,
      h_ij_inv, det_h);
#endif

  // Otherwise run the CPU Version to setup h_ij, h_ij_inv, det_h
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_range_iter conf_iter;

  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *h_ij_d = gkyl_array_cfetch(h_ij, midx);
    const double *h_ij_inv_d = gkyl_array_cfetch(h_ij_inv, midx);
    const double *det_h_d = gkyl_array_cfetch(det_h, midx);
    double *h_ij_quad = gkyl_array_fetch(up->h_ij_quad, midx);
    double *h_ij_inv_quad = gkyl_array_fetch(up->h_ij_inv_quad, midx);
    double *det_h_quad = gkyl_array_fetch(up->det_h_quad, midx);

    // Sum over basis 
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

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


struct gkyl_vlasov_dist_proj_on_basis* 
gkyl_vlasov_dist_proj_on_basis_inew(const struct gkyl_vlasov_dist_proj_on_basis_inp *inp)
{
  gkyl_vlasov_dist_proj_on_basis *up = gkyl_malloc(sizeof(*up));

  up->dist_proj = gkyl_dist_proj_type_acquire(inp->dist_proj);
  up->num_comp_mom = gkyl_dist_proj_type_num_mom(up->dist_proj);
  up->phase_grid = *inp->phase_grid;
  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;

  up->cdim = up->conf_basis.ndim;
  up->pdim = up->phase_basis.ndim;
  int vdim = up->pdim - up->cdim;
  up->num_conf_basis = up->conf_basis.num_basis;
  up->num_phase_basis = up->phase_basis.num_basis;
  up->use_gpu = inp->use_gpu;
  up->vel_map = 0;
  if (inp->vel_map != 0) {
    up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);
  }

  int num_quad = up->conf_basis.poly_order+1;
  // initialize data needed for conf-space quadrature 
  up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, 
    inp->quad_type, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, false);

  // initialize data needed for phase-space quadrature 
  up->tot_quad = init_quad_values(up->cdim, &up->phase_basis, 
    inp->quad_type, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, false);

  up->fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad); // Only used in CPU implementation.

  // To avoid creating iterators over ranges in device kernel, we'll
  // create a map between phase-space and conf-space ordinates.
  int num_quad_v = num_quad;  // Hybrid basis have p=2 in velocity space.
  // hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (up->phase_basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    for (int d=0; d<vdim; d++) {
      is_vdim_p2[d] = true;
    }
  }
  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Number density ratio: num_ratio = n_target/n0 and bin_op memory to compute ratio
  // Used for fixing the density with simple rescaling
  if (up->use_gpu) {
    up->num_ratio = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->conf_basis.num_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_local_ncells, up->conf_basis.num_basis);
  }
  else {
    up->num_ratio = gkyl_array_new(GKYL_DOUBLE, up->conf_basis.num_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, up->conf_basis.num_basis);
    up->moms_dist_quad = gkyl_array_new(GKYL_DOUBLE, up->tot_conf_quad*up->num_comp_mom, conf_local_ext_ncells);
  }
    

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Allocate device copies of arrays needed for quadrature.

    int p2c_qidx_ho[up->phase_qrange.volume];
    up->p2c_qidx = (int*) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);

    // Allocate f_dist_quad at phase-space quadrature points
    // moms_dist_quad at configuration-space quadrature points.
    up->f_dist_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
      up->tot_quad, inp->conf_range_ext->volume*inp->vel_range->volume);
    up->moms_dist_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
      up->tot_conf_quad*up->num_comp_mom, inp->conf_range_ext->volume);

    // Allocate the memory for computing the specific phase nodal to modal calculation
    struct gkyl_mat_mm_array_mem *phase_nodal_to_modal_mem_ho;
    phase_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(up->num_phase_basis, up->tot_quad, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, false);

    // Compute the matrix A for the phase nodal to modal memory
    const double *phase_w = (const double*) up->weights->data;
    const double *phaseb_o = (const double*) up->basis_at_ords->data;
    for (int n=0; n<up->tot_quad; ++n){
      for (int k=0; k<up->num_phase_basis; ++k){
        gkyl_mat_set(phase_nodal_to_modal_mem_ho->A, k, n, phase_w[n]*phaseb_o[k+up->num_phase_basis*n]);
      }
    }
    
    // copy to device
    up->phase_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(up->num_phase_basis, up->tot_quad, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
    gkyl_mat_copy(up->phase_nodal_to_modal_mem->A, phase_nodal_to_modal_mem_ho->A);
    gkyl_mat_mm_array_mem_release(phase_nodal_to_modal_mem_ho);

    // initialize data needed for conf-space quadrature on device 
    up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, 
      inp->quad_type, num_quad,
      &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, up->use_gpu);

    // initialize data needed for phase-space quadrature on device 
    up->tot_quad = init_quad_values(up->cdim, &up->phase_basis, 
      inp->quad_type, num_quad,
      &up->ordinates, &up->weights, &up->basis_at_ords, up->use_gpu);

    int pidx[GKYL_MAX_DIM];
    for (int n=0; n<up->tot_quad; ++n) {
      gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
      int cqidx = gkyl_range_idx(&up->conf_qrange, pidx);
      p2c_qidx_ho[n] = cqidx;
    }
    gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);
  }
#endif

    // gkyl_array_clear(up->h_ij_quad, 0.0); 
    // gkyl_array_clear(up->h_ij_inv_quad, 0.0); 
    // gkyl_array_clear(up->det_h_quad, 0.0); 
    // gkyl_vlasov_dist_proj_on_basis_geom_quad_vars(up, inp->conf_range, inp->h_ij, inp->h_ij_inv, inp->det_h);
  
  return up;
}
  


static void
proj_on_basis(const gkyl_vlasov_dist_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_phase_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k) {
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
    }
  }
}


void 
gkyl_vlasov_dist_proj_on_basis_moments_advance(gkyl_vlasov_dist_proj_on_basis *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_dist)
  {
    int num_conf_basis = up->num_conf_basis;
    int tot_conf_quad = up->tot_conf_quad;
    int num_comp_mom = up->num_comp_mom;

    struct gkyl_range_iter conf_iter;
    gkyl_range_iter_init(&conf_iter, conf_range);

    while (gkyl_range_iter_next(&conf_iter)){
      long midx = gkyl_range_idx(conf_range,conf_iter.idx);

      const double *moms_dist_d = gkyl_array_cfetch(moms_dist, midx);
      double *moms_dist_quad = gkyl_array_fetch(up->moms_dist_quad, midx);

    
     //The input moments are modal.
     //Evaluate each moment at each configuration-space quadrature point.
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      for (int m=0; m<num_comp_mom; ++m) {
        double mom_quad = 0.0;

        for (int k=0; k<num_conf_basis; ++k) {
          mom_quad += moms_dist_d[m*num_conf_basis+k] * b_ord[k];
        }

        //Store all moments for quadrature point n contiguously.
        moms_dist_quad[n*num_comp_mom+m] = mom_quad;
      }
    }
  }
}


void
gkyl_vlasov_dist_proj_on_basis_advance(gkyl_vlasov_dist_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  struct gkyl_array *f_dist)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_vlasov_dist_proj_on_basis_advance_cu(up, phase_range, conf_range, f_dist);
#endif

  int pdim = up->pdim;
  int num_comp_mom = up->num_comp_mom;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *moms_dist_quad = gkyl_array_cfetch(up->moms_dist_quad, midx);
    
    // inner loop over velocity space
    gkyl_range_deflate(&vel_rng, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {

      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->phase_grid, pidx, xc);

      struct gkyl_range_iter qiter;
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {
        
        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        const double *xcomp_d = gkyl_array_cfetch(up->ordinates, pqidx);

        comp_to_phys(pdim, xcomp_d, up->phase_grid.dx, xc, xmu);
        
        // All distribution parameters corresponding to this configuration-space quadrature point.
        const double *moms = &moms_dist_quad[cqidx*num_comp_mom];

        double *fq = gkyl_array_fetch(up->fun_at_ords, pqidx);
        
        // Generic distribution evaluation.
        gkyl_dist_proj_type_calc(up->dist_proj, xmu, up->phase_grid.dx, pidx, moms, fq);
      }      
      
      // compute expansion coefficients of distributiuon function on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(f_dist, lidx));
    }
  }
}

void
gkyl_vlasov_dist_proj_on_basis_release(gkyl_vlasov_dist_proj_on_basis* up)
{
  if (up->vel_map != 0) {
    gkyl_velocity_map_release(up->vel_map);
  }

  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_array_release(up->fun_at_ords);
  gkyl_array_release(up->moms_dist_quad);
  gkyl_array_release(up->num_ratio);
  gkyl_dg_bin_op_mem_release(up->mem);

  if (up->is_canonical_pb) {
    gkyl_array_release(up->h_ij_quad);
    gkyl_array_release(up->h_ij_inv_quad);
    gkyl_array_release(up->det_h_quad);
  }

  if (up->use_gpu) {
    gkyl_cu_free(up->p2c_qidx);
    gkyl_array_release(up->f_dist_quad);
    gkyl_mat_mm_array_mem_release(up->phase_nodal_to_modal_mem);
  }

  gkyl_dist_proj_type_release(up->dist_proj);

  gkyl_free(up);
}

