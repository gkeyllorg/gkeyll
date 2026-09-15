#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mat.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_lte_proj_on_basis_priv.h>
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
gkyl_vlasov_lte_proj_on_basis_geom_quad_vars(gkyl_vlasov_lte_proj_on_basis *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array *h_ij,
  const struct gkyl_array *h_ij_inv, const struct gkyl_array *det_h,
  const struct gkyl_array *background_flows, const struct gkyl_array *effective_potential)
{
// Setup the intial geometric vars, on GPU 
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_vlasov_lte_proj_on_basis_geom_quad_vars_cu(up, conf_range, h_ij,
      h_ij_inv, det_h, background_flows, effective_potential);
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

    // Geometry for the extended can-pb when H has time-independant but spatially varying potentials
    if ( up->use_extended_hamil_def ) {
      const double *background_flows_d = gkyl_array_cfetch(background_flows, midx);
      const double *effective_potential_d = gkyl_array_cfetch(effective_potential, midx);
      double *background_flows_quad = gkyl_array_fetch(up->background_flows_quad, midx);
      double *effective_potential_quad = gkyl_array_fetch(up->effective_potential_quad, midx);

      // Sum over basis 
      for (int n=0; n<tot_conf_quad; ++n) {
        const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

        for (int k=0; k<num_conf_basis; ++k) {
          effective_potential_quad[n] += effective_potential_d[k]*b_ord[k];
          for (int j=0; j<vdim; ++j) {
            background_flows_quad[tot_conf_quad*j + n] += background_flows_d[num_conf_basis*j+k]*b_ord[k];
          }
        }
      }
    }
  }
}


struct gkyl_vlasov_lte_proj_on_basis* 
gkyl_vlasov_lte_proj_on_basis_inew(const struct gkyl_vlasov_lte_proj_on_basis_inp *inp)
{
  gkyl_vlasov_lte_proj_on_basis *up = gkyl_malloc(sizeof(*up));

  up->phase_grid = *inp->phase_grid;
  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;

  up->cdim = up->conf_basis.ndim;
  up->pdim = up->phase_basis.ndim;
  int vdim = up->pdim - up->cdim;
  up->num_conf_basis = up->conf_basis.num_basis;
  up->num_phase_basis = up->phase_basis.num_basis;
  up->use_gpu = inp->use_gpu;

  int num_quad = up->conf_basis.poly_order+1;
  // initialize data needed for conf-space quadrature 
  up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, 
    inp->quad_type, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, false);

  up->vel_map = 0;
  up->vmap = 0;
  up->jacob_vel_gauss = 0;
  up->use_vmap = false;
  if (inp->vel_map) {
    up->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
    if (inp->vel_map->is_mapped) {
      up->use_vmap = true;
      if (up->use_gpu) {
        up->vmap_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
        gkyl_cart_modal_tensor_cu_dev(up->vmap_basis_on_dev, 1, 3);
      }
      else {
        gkyl_cart_modal_tensor(&up->vmap_basis, 1, 3);
      }
      // Borrowed pointers; kept alive by the acquired vel_map.
      up->vmap = inp->vel_map->vmap;
      up->jacob_vel_gauss = inp->vel_map->jacob_vel_gauss;
    }
  }

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
  up->vel_qrange = get_qrange(vdim, vdim, num_quad_v, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);

  up->vel_range = *inp->vel_range;

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
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Allocate device copies of arrays needed for quadrature.

    int p2c_qidx_ho[up->phase_qrange.volume];
    up->p2c_qidx = (int*) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);
    int p2v_qidx_ho[up->phase_qrange.volume];
    up->p2v_qidx = (int*) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);

    // Allocate f_lte_quad at phase-space quadrature points
    // moms_lte_quad (n, V_drift, T/m) at configuration-space quadrature points.
    // expamp_quad, the exponential pre-factor in the LTE distribution, at quadrature points.
    up->f_lte_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
      up->tot_quad, inp->conf_range_ext->volume*inp->vel_range->volume);
    up->moms_lte_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
      up->tot_conf_quad*(vdim+2), inp->conf_range_ext->volume);
    up->expamp_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
      up->tot_conf_quad, inp->conf_range_ext->volume);

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

    int pidx[GKYL_MAX_DIM], qidx_vel[GKYL_MAX_DIM];
    for (int n=0; n<up->tot_quad; ++n) {
      gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
      // Fetch the configuration-space quadrature point corresponding to the phase-space quadrature point. 
      int cqidx = gkyl_range_idx(&up->conf_qrange, pidx);
      p2c_qidx_ho[n] = cqidx;

      // Fetch the velocity-space quadrature point corresponding to the phase-space quadrature point. 
      for (int i=0; i<vdim; ++i) {
        qidx_vel[i] = pidx[up->cdim+i];          
      }
      int vqidx = gkyl_range_idx(&up->vel_qrange, qidx_vel); 
      p2v_qidx_ho[n] = vqidx; 
    }
    gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);
    gkyl_cu_memcpy(up->p2v_qidx, p2v_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);
  }
#endif

  up->is_relativistic = false;
  if (inp->model_id == GKYL_MODEL_SR || inp->model_id == GKYL_MODEL_TRIAD_GR) {
    up->is_relativistic = true;
  }

  // Flag for using the extended hamiltonian defintion which includes background flows
  up->use_extended_hamil_def = inp->use_extended_hamil_def;

  up->is_canonical_pb = false;
  if (inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    up->is_canonical_pb = true;
    // Allocate and obtain geometric variables at quadrature points for canonical-pb
    // since these quantities are time-independent.
    if (up->use_gpu) { 
      up->h_ij_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
        up->tot_conf_quad*(vdim*(vdim+1)/2), inp->conf_range_ext->volume);
      up->h_ij_inv_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
        up->tot_conf_quad*(vdim*(vdim+1)/2), inp->conf_range_ext->volume);
      up->det_h_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
        up->tot_conf_quad, inp->conf_range_ext->volume);
      if (up->use_extended_hamil_def) {
        up->background_flows_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
          up->tot_conf_quad*vdim, inp->conf_range_ext->volume);
        up->effective_potential_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
          up->tot_conf_quad, inp->conf_range_ext->volume);
      }
    }
    else {
      up->h_ij_quad = gkyl_array_new(GKYL_DOUBLE, 
        up->tot_conf_quad*(vdim*(vdim+1)/2), inp->conf_range_ext->volume);
      up->h_ij_inv_quad = gkyl_array_new(GKYL_DOUBLE, 
        up->tot_conf_quad*(vdim*(vdim+1)/2), inp->conf_range_ext->volume);
      up->det_h_quad = gkyl_array_new(GKYL_DOUBLE, 
        up->tot_conf_quad, inp->conf_range_ext->volume);
      if (up->use_extended_hamil_def) {
        up->background_flows_quad = gkyl_array_new(GKYL_DOUBLE, 
          up->tot_conf_quad*vdim, inp->conf_range_ext->volume);
        up->effective_potential_quad = gkyl_array_new(GKYL_DOUBLE, 
          up->tot_conf_quad, inp->conf_range_ext->volume);
      }
    }
    gkyl_array_clear(up->h_ij_quad, 0.0); 
    gkyl_array_clear(up->h_ij_inv_quad, 0.0); 
    gkyl_array_clear(up->det_h_quad, 0.0); 
    if (up->use_extended_hamil_def) {
      gkyl_array_clear(up->background_flows_quad, 0.0); 
      gkyl_array_clear(up->effective_potential_quad, 0.0); 
    }
    gkyl_vlasov_lte_proj_on_basis_geom_quad_vars(up, inp->conf_range, inp->h_ij, inp->h_ij_inv, inp->det_h, 
      inp->use_extended_hamil_def ? inp->background_flows : 0,
      inp->use_extended_hamil_def ? inp->effective_potential : 0);
  }

  up->is_bimaxwellian = inp->is_bimaxwellian; 

  // Store a LTE moment calculation updater to compute and correct the density
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .vel_map = inp->vel_map,
    .hamil_range = inp->hamil_range,
    .hamil = inp->hamil,
    .model_id = inp->model_id,
    .hamil_id = inp->hamil_id,
    .gamma_inv = inp->gamma_inv,
    .h_ij = inp->h_ij,
    .h_ij_inv = inp->h_ij_inv,
    .det_h = inp->det_h,
    .use_extended_hamil_def = inp->use_extended_hamil_def,
    .effective_potential = inp->effective_potential,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_vlasov_lte_moments_inew( &inp_mom );

  return up;
}

static void
proj_on_basis(const gkyl_vlasov_lte_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
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
gkyl_vlasov_lte_proj_on_basis_advance(gkyl_vlasov_lte_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_lte, struct gkyl_array *f_lte)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_vlasov_lte_proj_on_basis_advance_cu(up, phase_range, conf_range, moms_lte, f_lte);
#endif

  double f_floor = 1.0e-40;  
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;  

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double n_quad[tot_conf_quad], V_drift_quad[tot_conf_quad][vdim], T_over_m_quad[tot_conf_quad];
  double Tperp_over_m_quad[tot_conf_quad];
  double expamp_quad[tot_conf_quad];

  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *moms_lte_d = gkyl_array_cfetch(moms_lte, midx);
    const double *n_d = moms_lte_d;
    const double *V_drift_d = &moms_lte_d[num_conf_basis];
    const double *T_over_m_d = &moms_lte_d[num_conf_basis*(vdim+1)];

    // Sum over basis for given LTE moments (n, V_drift, T/m) in the stationary frame
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      // Zero out quadrature values
      n_quad[n] = 0.0;
      for (int d=0; d<vdim; ++d) {
        V_drift_quad[n][d] = 0.0;
      }
      T_over_m_quad[n] = 0.0;

      // Compute the configuration-space quadrature
      for (int k=0; k<num_conf_basis; ++k) {
        n_quad[n] += n_d[k]*b_ord[k];
        for (int d=0; d<vdim; ++d) {
          V_drift_quad[n][d] += V_drift_d[num_conf_basis*d+k]*b_ord[k];
        }
        T_over_m_quad[n] += T_over_m_d[k]*b_ord[k];
      }   
      // Amplitude of the exponential.
      if ((n_quad[n] > 0.0) && (T_over_m_quad[n] > 0.0)) {
        if (up->is_relativistic) {
          if (up->is_bimaxwellian) {
            const double *Tperp_over_m_d = &moms_lte_d[num_conf_basis*(vdim+2)];
            Tperp_over_m_quad[n] = 0.0;
            for (int k=0; k<num_conf_basis; ++k) { 
              Tperp_over_m_quad[n] += Tperp_over_m_d[k]*b_ord[k];
            } 
            expamp_quad[n] = n_quad[n]*(1.0/(4.0*GKYL_PI*Tperp_over_m_quad[n]*(1 + (Tperp_over_m_quad[n] - T_over_m_quad[n]))))*(sqrt(2.0*T_over_m_quad[n]/GKYL_PI));
          }
          else {
            expamp_quad[n] = n_quad[n]*(1.0/(4.0*GKYL_PI*T_over_m_quad[n]))*(sqrt(2.0*T_over_m_quad[n]/GKYL_PI));
          }
        }
        else if (up->is_canonical_pb) { 
          const double *det_h_quad = gkyl_array_cfetch(up->det_h_quad, midx);
          expamp_quad[n] = (1.0/det_h_quad[n])*n_quad[n]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[n], vdim));
        }
        else {
          expamp_quad[n] = n_quad[n]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[n], vdim));
        }
      }
      else {
        expamp_quad[n] = 0.0;
      }      
    }

    // inner loop over velocity space
    double jacob_vel_qidx;
    int qidx_vel[GKYL_MAX_DIM];
    gkyl_range_deflate(&vel_rng, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->phase_grid, pidx, xc);

      struct gkyl_range_iter qiter;
      // compute LTE distribution function at phase-space quadrature nodes
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        if (up->use_vmap) {
          const double *xcomp_d = gkyl_array_cfetch(up->ordinates, pqidx);
          long loc_vel = gkyl_range_idx(&up->vel_range, vel_iter.idx);
          const double *vmap_d = gkyl_array_cfetch(up->vmap, loc_vel);
          const double *jacob_vel_quad_d = gkyl_array_cfetch(up->jacob_vel_gauss, loc_vel);
          double xcomp[1];
          for (int vd=0; vd<vdim; vd++) {
            xcomp[0] = xcomp_d[cdim+vd];
            xmu[cdim+vd] = up->vmap_basis.eval_expand(xcomp, vmap_d+vd*up->vmap_basis.num_basis);
          }
          for (int i=0; i<vdim; ++i) {
            qidx_vel[i] = qiter.idx[cdim+i];          
          }
          int vqidx = gkyl_range_idx(&up->vel_qrange, qidx_vel);
          jacob_vel_qidx = jacob_vel_quad_d[vqidx];
        }
        else {
          comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, pqidx),
            up->phase_grid.dx, xc, xmu);
          jacob_vel_qidx = 1.0;
        }

        double *fq = gkyl_array_fetch(up->fun_at_ords, pqidx);
        fq[0] = f_floor;
        if (T_over_m_quad[cqidx] > 0.0) {
          if (up->is_relativistic) {
            double vv = 0.0;
            double vu = 0.0;
            double uu = 0.0;
            // V_drift_quad is the spatial component of the four-velocity u_i = GammaV*V_drift
            for (int d=0; d<vdim; ++d) {
              vv += (V_drift_quad[cqidx][d]*V_drift_quad[cqidx][d]);
              vu += (V_drift_quad[cqidx][d]*xmu[cdim+d]);
              uu += (xmu[cdim+d]*xmu[cdim+d]);
            }
            double GammaV_quad = sqrt(1.0 + vv);
            double gamma_quad = sqrt(1.0 + uu);
            if (up->is_bimaxwellian) {
            // NOTE: THIS PROJECTION ROUTINE CURRENTLY ASSUMES Tpar = Txx AND THUS THE INITIAL MAGNETIC FIELD IS B = B0 x_hat (JJ: 10/23/24)
              fq[0] += jacob_vel_qidx*expamp_quad[cqidx]*exp((1.0/T_over_m_quad[cqidx]) 
                - (1.0/Tperp_over_m_quad[cqidx])*(GammaV_quad*gamma_quad - vu) 
                - ((1.0/T_over_m_quad[cqidx]) - (1.0/Tperp_over_m_quad[cqidx]))*sqrt(1.0 + 
                (1.0 + vv)*xmu[cdim]*xmu[cdim] - 2.0*GammaV_quad*gamma_quad*xmu[cdim]*V_drift_quad[cqidx][0] + (1.0 + uu)*V_drift_quad[cqidx][0]*V_drift_quad[cqidx][0]));
            }
            else {
              fq[0] += jacob_vel_qidx*expamp_quad[cqidx]*exp((1.0/T_over_m_quad[cqidx]) 
                - (1.0/T_over_m_quad[cqidx])*(GammaV_quad*gamma_quad - vu));
            }
          }
          else if (up->is_canonical_pb) {
            // Assumes a (particle) hamiltonian in canocial form: H = 1/2 g^{ij} p_i p_j
            const double *h_ij_inv_quad = gkyl_array_cfetch(up->h_ij_inv_quad, midx);
            const double *background_flows_quad;
            if (up->use_extended_hamil_def) {
              background_flows_quad = gkyl_array_cfetch(up->background_flows_quad, midx);
            }
            double efact = 0.0;
            for (int d0=0; d0<vdim; ++d0) {
              for (int d1=d0; d1<vdim; ++d1) {
                int sym_tensor_index = (d0*(2*vdim - d0 + 1))/2 + (d1-d0);
                // Grab the spatial metric component, the ctx includes geometry that isn't 
                // part of the canonical set of variables, like R on the surf of a sphere
                // q_can includes the canonical variables list
                double h_ij_inv_loc = h_ij_inv_quad[tot_conf_quad*sym_tensor_index + cqidx]; 
                // For off-diagnol components, we need to count these twice, due to symmetry
                int sym_fact = (d0 == d1) ? 1 : 2;
                if (up->use_extended_hamil_def) {
                  //printf(" (d0: %d, d1: %d) background_flows_quad[tot_conf_quad*d0 + cqidx]: %1.16e, background_flows_quad[tot_conf_quad*d1 + cqidx]: %1.16e\n",d0, d1, background_flows_quad[tot_conf_quad*d0 + cqidx], background_flows_quad[tot_conf_quad*d1 + cqidx]);
                  efact += sym_fact*h_ij_inv_loc*(xmu[cdim+d0]-V_drift_quad[cqidx][d0]-background_flows_quad[tot_conf_quad*d0 + cqidx]) 
                    * (xmu[cdim+d1]-V_drift_quad[cqidx][d1]-background_flows_quad[tot_conf_quad*d1 + cqidx]);
                }
                else {
                  efact += sym_fact*h_ij_inv_loc*(xmu[cdim+d0]-V_drift_quad[cqidx][d0]) * (xmu[cdim+d1]-V_drift_quad[cqidx][d1]);
                }
              }
            }
            fq[0] += expamp_quad[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
          }
          else {
            double efact = 0.0;
            for (int d=0; d<vdim; ++d) {
              efact += (xmu[cdim+d]-V_drift_quad[cqidx][d])*(xmu[cdim+d]-V_drift_quad[cqidx][d]);
            }
            // The velocity-space Jacobian is included since the projected
            // distribution is the evolved quantity Jf on mapped velocity grids.
            fq[0] += jacob_vel_qidx*expamp_quad[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
          }
        }
      }
      // compute expansion coefficients of LTE distributiuon function on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(f_lte, lidx));
    }
  }
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

void
gkyl_vlasov_lte_proj_on_basis_release(gkyl_vlasov_lte_proj_on_basis* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_array_release(up->fun_at_ords);

  if (up->vel_map) {
    gkyl_vlasov_velocity_map_release(up->vel_map);
  }

  gkyl_array_release(up->num_ratio);
  gkyl_dg_bin_op_mem_release(up->mem);

  if (up->is_canonical_pb) {
    gkyl_array_release(up->h_ij_quad);
    gkyl_array_release(up->h_ij_inv_quad);
    gkyl_array_release(up->det_h_quad);
    if (up->use_extended_hamil_def) {
      gkyl_array_release(up->background_flows_quad);
      gkyl_array_release(up->effective_potential_quad);
    }
  }

  if (up->use_gpu) {
    if (up->use_vmap) {
      gkyl_cu_free(up->vmap_basis_on_dev);
    }
    gkyl_cu_free(up->p2c_qidx);
    gkyl_cu_free(up->p2v_qidx);
    gkyl_array_release(up->f_lte_quad);
    gkyl_array_release(up->moms_lte_quad);
    gkyl_array_release(up->expamp_quad);
    gkyl_mat_mm_array_mem_release(up->phase_nodal_to_modal_mem);
  }

  gkyl_vlasov_lte_moments_release(up->moments_up);

  gkyl_free(up);
}
