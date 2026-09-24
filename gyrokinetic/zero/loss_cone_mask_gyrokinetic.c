#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mat.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>
#include <assert.h>

//
// mu_bound = (0.5*mass*pow(vpar,2)+charge*Delta_phi)/(bmag[0]*(Rm-1));
//          = 0.5*mass*pow(vpar,2)/(bmag[0]*(Rm-1)) + charge*Delta_phi/(bmag[0]*(Rm-1));
//          = 0.5*mass*pow(vpar,2)/(bmag_max-bmag[0]) + charge*(phi-phi_m)/(bmag_max-bmag[0]);
//

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
init_quad_values(int cdim, const struct gkyl_basis *basis, enum gkyl_quad_type qtype,
  int num_quad, struct gkyl_array **ordinates, struct gkyl_array **weights,
  struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  int num_quad_v = num_quad;
  // Hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[2] = {false};  // 2 is the max vdim for GK.
  if (num_quad > 1 && basis->b_type == GKYL_BASIS_MODAL_GKHYBRID) {
    num_quad_v = num_quad+1;
    is_vdim_p2[0] = true;  // only vpar is quadratic in GK hybrid.
  }

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];

  if (qtype == GKYL_GAUSS_QUAD) {
    if (num_quad <= gkyl_gauss_max) {
      // Use pre-computed values if possible (these are more accurate than computing them on the fly).
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
  else if (qtype == GKYL_GAUSS_LOBATTO_QUAD) {
    assert( (num_quad > 1) && (num_quad <= gkyl_gauss_max) );
    // Gauss-Lobatto quadrature
    memcpy(ordinates1, gkyl_gauss_lobatto_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_lobatto_weights[num_quad], sizeof(double[num_quad]));

    assert( (num_quad_v > 1) && (num_quad_v <= gkyl_gauss_max) );
    memcpy(ordinates1_v, gkyl_gauss_lobatto_ordinates[num_quad_v], sizeof(double[num_quad_v]));
    memcpy(weights1_v, gkyl_gauss_lobatto_weights[num_quad_v], sizeof(double[num_quad_v]));
  }
  else {
    fprintf(stderr, "Quadrature rule not available.\n");
    assert(false);
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
    for (int i=0; i<cdim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];

    for (int i=cdim; i<ndim; ++i)
      ord[i] = is_vdim_p2[i-cdim] ? 
        ordinates1_v[iter.idx[i]-qrange.lower[i]] : ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<cdim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];

    for (int i=cdim; i<ndim; ++i)
      wgt[0] *= is_vdim_p2[i-cdim] ? 
        weights1_v[iter.idx[i]-qrange.lower[i]] : weights1[iter.idx[i]-qrange.lower[i]];
  }

  // Pre-compute basis functions at ordinates.
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu)
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  else
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);

  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));

  // Copy host array to device array.
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

static void
gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad(gkyl_loss_cone_mask_gyrokinetic *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const double *bmag_max)
{
  // Get bmag_max-bmag at quadrature nodes.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu(up, conf_range, bmag, bmag_max);
#endif

  int cdim = up->cdim, pdim = up->pdim;

  int tot_quad_conf = up->tot_quad_conf;
  int num_basis_conf = up->num_basis_conf;

  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *bmag_d = gkyl_array_cfetch(bmag, linidx);
    double *Dbmag_quad = gkyl_array_fetch(up->Dbmag_quad, linidx);

    // Sum over basis 
    for (int n=0; n<tot_quad_conf; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->basis_at_ords_conf, n);
      for (int k=0; k<num_basis_conf; ++k)
        Dbmag_quad[n] += bmag_d[k]*b_ord[k];

      Dbmag_quad[n] = bmag_max[0] - Dbmag_quad[n];
    }
  }
}

struct gkyl_loss_cone_mask_gyrokinetic* 
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp)
{
  gkyl_loss_cone_mask_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->grid_phase = inp->phase_grid;
  up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);
  up->mass = inp->mass;
  up->charge = inp->charge;

  up->cdim = inp->conf_basis->ndim;
  up->pdim = inp->phase_basis->ndim;

  up->cellwise_trap_loss = inp->cellwise_trap_loss;
  int num_quad = inp->num_quad? inp->num_quad : inp->phase_basis->poly_order+1;
  up->norm_fac = 1;
  if (!up->cellwise_trap_loss)
    up->norm_fac = num_quad == 1? 1.0/pow(sqrt(2.0),up->pdim) : 1.0;

  if (num_quad == 1) {
    up->num_basis_conf = 1;
    up->num_basis_phase = 1;
  }
  else {
    up->num_basis_conf = inp->conf_basis->num_basis;
    up->num_basis_phase = inp->phase_basis->num_basis;
  }
  up->use_gpu = inp->use_gpu;

  // Initialize data needed for conf-space quadrature.
  up->tot_quad_conf = init_quad_values(up->cdim, inp->conf_basis, inp->qtype, num_quad,
    &up->ordinates_conf, &up->weights_conf, &up->basis_at_ords_conf, false);

  // Initialize data needed for phase-space quadrature.
  up->tot_quad_phase = init_quad_values(up->cdim, inp->phase_basis, inp->qtype, num_quad,
    &up->ordinates_phase, &up->weights_phase, &up->basis_at_ords_phase, false);

  up->fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad_phase); // Only used in CPU implementation.

  // To avoid creating iterators over ranges in device kernel, we'll
  // create a map between phase-space and conf-space ordinates.
  int num_quad_v = num_quad;  // Hybrid basis have p=2 in velocity space.
  // hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[2] = {false};  // 2 is the max vdim for GK.
  if (num_quad > 1 && inp->phase_basis->b_type == GKYL_BASIS_MODAL_GKHYBRID) {
    num_quad_v = num_quad+1;
    is_vdim_p2[0] = true;  // only vpar is quadratic in GK hybrid.
  }
  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Allocate device copies of arrays needed for quadrature.

    int p2c_qidx_ho[up->phase_qrange.volume];
    up->p2c_qidx = (int*) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);

    // Allocate mask_quad at phase-space quadrature points.
    // Dbmag_quad at configuration-space quadrature points.
    // qDphiDbmag_quad, the term proportional to (phi-phi_m)/(bmag_max-bmag), at quadrature points.
    up->mask_out_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_quad_phase,
      inp->conf_range_ext->volume*inp->vel_range->volume);
    up->qDphiDbmag_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_quad_conf, inp->conf_range_ext->volume);

    // Allocate the memory for computing the specific phase nodal to modal calculation
    struct gkyl_mat_mm_array_mem *phase_nodal_to_modal_mem_ho;
    phase_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(up->num_basis_phase, up->tot_quad_phase, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, false);

    // Compute the matrix A for the phase nodal to modal memory
    const double *phase_w = (const double*) up->weights_phase->data;
    const double *phaseb_o = (const double*) up->basis_at_ords_phase->data;
    for (int n=0; n<up->tot_quad_phase; ++n) {
      for (int k=0; k<up->num_basis_phase; ++k)
        gkyl_mat_set(phase_nodal_to_modal_mem_ho->A, k, n, phase_w[n]*phaseb_o[k+up->num_basis_phase*n]);
    }
    
    // Copy to device
    up->phase_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(up->num_basis_phase, up->tot_quad_phase, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
    gkyl_mat_copy(up->phase_nodal_to_modal_mem->A, phase_nodal_to_modal_mem_ho->A);
    gkyl_mat_mm_array_mem_release(phase_nodal_to_modal_mem_ho);

    // Initialize data needed for conf-space quadrature on device.
    up->tot_quad_conf = init_quad_values(up->cdim, inp->conf_basis, inp->qtype, num_quad,
      &up->ordinates_conf, &up->weights_conf, &up->basis_at_ords_conf, up->use_gpu);

    // Initialize data needed for phase-space quadrature on device.
    up->tot_quad_phase = init_quad_values(up->cdim, inp->phase_basis, inp->qtype, num_quad,
      &up->ordinates_phase, &up->weights_phase, &up->basis_at_ords_phase, up->use_gpu);

    int pidx[GKYL_MAX_DIM];
    for (int n=0; n<up->tot_quad_phase; ++n) {
      gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
      int cqidx = gkyl_range_idx(&up->conf_qrange, pidx);
      p2c_qidx_ho[n] = cqidx;
    }
    gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Allocate and obtain bmag_max-bmag at quadrature points.
  if (up->use_gpu) 
    up->Dbmag_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->tot_quad_conf, inp->conf_range_ext->volume);
  else
    up->Dbmag_quad = gkyl_array_new(GKYL_DOUBLE, up->tot_quad_conf, inp->conf_range_ext->volume);

  gkyl_array_clear(up->Dbmag_quad, 0.0); 
  gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad(up, inp->conf_range, inp->bmag, inp->bmag_max);

  // Save the location of bmag_max in this updater.
  if (up->use_gpu) {
    up->bmag_max_loc = gkyl_cu_malloc(sizeof(double)*up->cdim);
    gkyl_cu_memcpy(up->bmag_max_loc, inp->bmag_max_loc, sizeof(double)*up->cdim, GKYL_CU_MEMCPY_D2D);
  }
  else {
    up->bmag_max_loc = gkyl_malloc(sizeof(double)*up->cdim);
    memcpy(up->bmag_max_loc, inp->bmag_max_loc, sizeof(double)*up->cdim);
  }
    
  return up;
}

static void
proj_on_basis(const gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis_phase;
  int tot_quad = up->tot_quad_phase;

  const double* GKYL_RESTRICT weights = up->weights_phase->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords_phase->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

static void
nod_to_mod_reduce(const gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis_phase;
  int tot_quad = up->tot_quad_phase;

  const double* GKYL_RESTRICT weights = up->weights_phase->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords_phase->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  f[0] = 1.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    if (func_at_ords[imu] < 1e-14) {
      f[0] = 0.0;
      break;
    }
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_advance(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *phi, const double *phi_m, struct gkyl_array *mask_out)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_loss_cone_mask_gyrokinetic_advance_cu(up, phase_range, conf_range, 
      phi, phi_m, mask_out);
#endif

  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;

  int tot_quad_conf = up->tot_quad_conf;
  int num_basis_conf = up->num_basis_conf;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.0};
  double phi_quad[tot_quad_conf];
  double qDphiDbmag_quad[tot_quad_conf]; // charge*(phi-phi_m)/(bmag_max-bmag[0]).

  // Outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space.
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx_conf = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *phi_d = gkyl_array_cfetch(phi, linidx_conf);
    const double *Dbmag_quad = gkyl_array_cfetch(up->Dbmag_quad, linidx_conf);

    // Sum over basis for given potential phi.
    for (int n=0; n<tot_quad_conf; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->basis_at_ords_conf, n);

      // Compute the configuration-space quadrature
      phi_quad[n] = 0.0;
      for (int k=0; k<num_basis_conf; ++k)
        phi_quad[n] += phi_d[k]*b_ord[k];

      if (Dbmag_quad[n] > 0.0)
        qDphiDbmag_quad[n] = up->charge*(phi_quad[n]-phi_m[0])/Dbmag_quad[n];
      else
        qDphiDbmag_quad[n] = 0.0;
    }

    // Inner loop over velocity space.
    gkyl_range_deflate(&vel_rng, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      long linidx_phase = gkyl_range_idx(&vel_rng, vel_iter.idx);

      // Compute the mask function at phase-space quadrature nodes.
      struct gkyl_range_iter qiter;
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        const double *xcomp_d = gkyl_array_cfetch(up->ordinates_phase, pqidx);

        // Convert comp position coordinate to phys pos coord.
        gkyl_rect_grid_cell_center(up->grid_phase, pidx, xc);
        log_to_comp(up->cdim, xcomp_d, up->grid_phase->dx, xc, xmu);

        // Convert comp velocity coordinate to phys velocity coord.
        const struct gkyl_velocity_map *gvm = up->vel_map;
        long linidx_vel = gkyl_range_idx(&gvm->local_ext_vel, vel_iter.idx);
        const double *vmap_d = gkyl_array_cfetch(gvm->vmap, linidx_vel);
        double xcomp[1];
        for (int vd = 0; vd < vdim; vd++) {
          xcomp[0] = xcomp_d[cdim+vd];
          xmu[cdim+vd] = gvm->vmap_basis->eval_expand(xcomp, vmap_d+vd*gvm->vmap_basis->num_basis);
        }

        // KEparDbmag = 0.5*mass*pow(vpar,2)/(bmag_max-bmag[0]).
        double KEparDbmag = 0.0;
        if (Dbmag_quad[cqidx] > 0.0)
          KEparDbmag = 0.5*up->mass*pow(xmu[cdim], 2.0)/Dbmag_quad[cqidx];
        else
          KEparDbmag = 0.0;

	double mu_bound = GKYL_MAX2(0.0, KEparDbmag+qDphiDbmag_quad[cqidx]);

        double *fq = gkyl_array_fetch(up->fun_at_ords, pqidx);
	if (mu_bound < xmu[cdim+1] && fabs(xmu[cdim-1]) < fabs(up->bmag_max_loc[cdim-1])) 
          fq[0] = 1.0 * up->norm_fac;
        else
          fq[0] = 0.0;
      }
      // Compute DG expansion coefficients of the mask.
      if (up->cellwise_trap_loss)
        nod_to_mod_reduce(up, up->fun_at_ords, gkyl_array_fetch(mask_out, linidx_phase));
      else
        proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(mask_out, linidx_phase));
    }
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic* up)
{
  gkyl_velocity_map_release(up->vel_map);

  gkyl_array_release(up->ordinates_phase);
  gkyl_array_release(up->weights_phase);
  gkyl_array_release(up->basis_at_ords_phase);

  gkyl_array_release(up->ordinates_conf);
  gkyl_array_release(up->weights_conf);
  gkyl_array_release(up->basis_at_ords_conf);

  gkyl_array_release(up->fun_at_ords);
  gkyl_array_release(up->Dbmag_quad);

  if (up->use_gpu) {
    gkyl_cu_free(up->p2c_qidx);
    gkyl_array_release(up->mask_out_quad);
    gkyl_array_release(up->qDphiDbmag_quad);
    gkyl_mat_mm_array_mem_release(up->phase_nodal_to_modal_mem);
    gkyl_cu_free(up->bmag_max_loc);
  }
  else {
    gkyl_free(up->bmag_max_loc);
  }

  gkyl_free(up);
}
