#include <math.h>
#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_proj_maxwellian_pots_on_basis {
  struct gkyl_rect_grid grid;
  int cdim; // Configuration space dimension
  int pdim; // Phase space dimension
  int num_quad; // Number of 1D quadrature points
  
  bool use_gpu;

  struct gkyl_basis conf_basis;
  struct gkyl_basis surf_basis;

  int num_conf_basis; // number of configuration space basis functions
  int num_phase_basis; // number of phase space basis functions
  int num_surf_basis; // Number of surface basis functions
 
  int tot_conf_quad;
  struct gkyl_range conf_qrange;
  struct gkyl_array *conf_ordinates;
  struct gkyl_array *conf_weights;
  struct gkyl_array *conf_basis_at_ords;
  struct gkyl_array *conf_basis_at_nodes;

  int tot_quad;
  struct gkyl_range phase_qrange;
  struct gkyl_array *ordinates;
  struct gkyl_array *weights;
  struct gkyl_array *basis_at_ords;

  int tot_surf_quad;
  struct gkyl_range surf_qrange;
  struct gkyl_array *surf_ordinates;
  struct gkyl_array *surf_weights;
  struct gkyl_array *surf_basis_at_ords;

  // Potentials and derivatives evaluated at quadrature nodes
  struct gkyl_array *fpo_h_at_ords;
  struct gkyl_array *fpo_g_at_ords;
  struct gkyl_array *fpo_dhdv_at_ords;
  struct gkyl_array *fpo_d2gdv2_at_ords;

  // Potentials and derivatives evaluated at surface quadrature nodes
  struct gkyl_array *fpo_h_at_surf_ords;
  struct gkyl_array *fpo_g_at_surf_ords;
  struct gkyl_array *fpo_dhdv_at_surf_ords;
  struct gkyl_array *fpo_dgdv_at_surf_ords;
  struct gkyl_array *fpo_d2gdv2_at_surf_ords;

  struct gkyl_array *prim_moms_conf_quad; // Array of primitive moments evaluated at conf space quadrature nodes 
  struct gkyl_array *prim_moms_conf_nodes; // Array of primitive moments evaluated at conf space nodes 
  struct gkyl_array *pot_phase_quad; // Array of potential H or G evaluated at phase space quadrature nodes
  struct gkyl_array *pot_surf_quad; // Array of potentials (H, G) evaluated at phase space surface quadrature nodes in each velocity direction
  struct gkyl_array *pot_deriv_surf_quad; // Array of derivatives of potentials (dH/dv, d2G/dv2) evaluated at phase space surface quadrature nodes in each velocity direction 
  struct gkyl_array *sol_pot_surf_modal; // Array to store solutions of quad to modal conversion.

  struct gkyl_mat_mm_array_mem *phase_quad_nodal_to_modal_mem;
  struct gkyl_mat_mm_array_mem *surf_quad_nodal_to_modal_mem;
  struct gkyl_mat_mm_array_mem *dgdv_surf_quad_nodal_to_modal_mem;

  int *p2c_qidx;  // Mapping between conf-space and phase-space ordinates.
  int *surf2c_qidx;  // Mapping between conf-space and phase-space surface ordinates.
  int *surf2c_nidx;  // Mapping between conf-space and phase-space surface nodes.

  struct gkyl_basis *surf_basis_dev; // Basis that lives on device for nodal-modal conversion
};

GKYL_CU_DH
static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

GKYL_CU_DH
static inline void
surf_comp_to_phys(int dir, int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) {
    if (d < dir)
      xout[d] = 0.5*dx[d]*eta[d]+xc[d];
    else if (d > dir)
      xout[d] = 0.5*dx[d]*eta[d-1]+xc[d];
    else 
      xout[d] = 0.0;
  }
}

GKYL_CU_DH
static inline void
edge_idx_to_phase_idx(int ndim, int dir, const int *surf_idx, int edge_idx, int *phase_idx) {
  phase_idx[dir] = edge_idx;
  for (int i=0; i<ndim; ++i) {
    if (i < dir) phase_idx[i] = surf_idx[i];
    else if (i == dir) phase_idx[i] = edge_idx;
    else phase_idx[i] = surf_idx[i-1];
  }
}

static void
proj_on_basis(const gkyl_proj_maxwellian_pots_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_phase_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = (const double*)up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = (const double*)up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = (const double*)fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

static void
proj_on_surf_basis(const gkyl_proj_maxwellian_pots_on_basis *up, int offset, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_surf_basis = up->num_surf_basis;
  int tot_surf_quad = up->tot_surf_quad;

  const double* GKYL_RESTRICT weights = (const double*)up->surf_weights->data;
  const double* GKYL_RESTRICT surf_basis_at_ords = (const double*)up->surf_basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = (const double*)fun_at_ords->data;

  for (int k=0; k<num_surf_basis; ++k) f[offset*num_surf_basis+k] = 0.0;
  
  for (int imu=0; imu<tot_surf_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_surf_basis; ++k)
      f[offset*num_surf_basis+k] += tmp*surf_basis_at_ords[num_surf_basis*imu+k];
  }
}

static void
nod2mod(int num_ret_vals, const struct gkyl_basis *basis, const struct gkyl_array *fun_at_nodes, double *f) {
  const double *fao = (const double*)gkyl_array_cfetch(fun_at_nodes, 0);

  int num_basis = basis->num_basis;
  double fnodal[num_basis];
  for (int i=0; i<num_ret_vals; ++i) {
    for (int k=0; k<num_basis; ++k) {
      fnodal[k] = fao[num_ret_vals*k+i];
    }

    basis->nodal_to_modal(fnodal, &f[num_basis*i]);
  }
}

// create range to loop over quadrature points.
static inline struct gkyl_range get_qrange(int cdim, int dim, int num_quad, int num_quad_v, bool *is_vdim_p2) {
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  for (int i=cdim; i<dim; ++i) qshape[i] = is_vdim_p2[i-cdim]? num_quad_v : num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

// create range to loop over nodes.
static inline struct gkyl_range get_nrange(int cdim, int num_conf_basis) {
  int nshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) nshape[i] = num_conf_basis;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, cdim, nshape);
  return nrange;
}

// Sets ordinates, weights and basis functions at ords.
// Returns the total number of quadrature nodes
static int
init_quad_values(int cdim, const struct gkyl_basis *basis, int num_quad, struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  int vdim = ndim-cdim;
  int num_quad_v = num_quad;
  // hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    for (int d=0; d<vdim; d++) is_vdim_p2[d] = true;
  }

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];
  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }
  if (num_quad_v <= gkyl_gauss_max) {
    memcpy(ordinates1_v, gkyl_gauss_ordinates[num_quad_v], sizeof(double[num_quad_v]));
    memcpy(weights1_v, gkyl_gauss_weights[num_quad_v], sizeof(double[num_quad_v]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1_v, weights1_v, num_quad_v);
  }

  struct gkyl_range qrange = get_qrange(cdim, ndim, num_quad, num_quad_v, is_vdim_p2);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *ordinates = gkyl_array_cu_dev_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  } else {
    *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = (double *)gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<cdim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    for (int i=cdim; i<ndim; ++i)
      ord[i] = is_vdim_p2[i-cdim]? ordinates1_v[iter.idx[i]-qrange.lower[i]] :
                                   ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = (double *)gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<cdim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
    for (int i=cdim; i<ndim; ++i)
      wgt[0] *= is_vdim_p2[i-cdim]? weights1_v[iter.idx[i]-qrange.lower[i]] :
                                    weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu) 
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  else
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval((double *)gkyl_array_fetch(ordinates_ho, n), (double *)gkyl_array_fetch(basis_at_ords_ho, n));

  // copy host array to device array
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

// Note we're multiplying the input vtsq=T/m by sqrt(2)
GKYL_CU_D
static inline double eval_fpo_h(double den,
  double rel_speed, double vtsq) 
{
  double vth = sqrt(2.0*vtsq);
  return den/rel_speed * erf(rel_speed/vth);
}

GKYL_CU_D
static inline double eval_fpo_g(double den, double rel_speed,
  double vtsq) 
{
  double rel_speedsq = pow(rel_speed, 2);
  double vth = sqrt(2.0*vtsq);
  return  den*vth*
    (1.0/(sqrt(GKYL_PI))*exp(-rel_speedsq/pow(vth, 2)) + 
    erf(rel_speed/vth)*(vth/(2.0*rel_speed) + 
    rel_speed/vth));
}

GKYL_CU_D
static inline double eval_fpo_dhdv(double den,
  double rel_vel_in_dir, double vtsq, double rel_speed) 
{
  double rel_speedsq = pow(rel_speed, 2);
  double vth = sqrt(2.0*vtsq);
  double exp_term = exp(-pow(rel_speed/vth, 2));
  double erf_term = erf(rel_speed / vth);
  double term1 =  (2.0 / sqrt(GKYL_PI)) * exp_term * (rel_vel_in_dir / (vth * pow(rel_speed, 2)));
  double term2 = erf_term * ((-1.0) * rel_vel_in_dir / pow(rel_speed, 3)); 

  double dHdvi = den*(term1 + term2);
  return dHdvi;
}

GKYL_CU_D
static inline double eval_fpo_dgdv(double den, double rel_vel_in_dir,
  double vtsq, double rel_speedsq) {
  double rel_speed = sqrt(rel_speedsq);
  double vth = sqrt(2.0*vtsq);
  double dGdvi = den*rel_vel_in_dir/pow(rel_speed, 3)*(
      vth*exp(-rel_speedsq/pow(vth,2))*rel_speed/sqrt(GKYL_PI) +
      (2.0*rel_speedsq - pow(vth,2))*erf(rel_speed/vth)/2.0);
  return dGdvi;
}

GKYL_CU_D
static inline double eval_fpo_d2gdv2(double den,
  double rel_vel_in_dir, double vtsq, double rel_speed) {
  double rel_speedsq = pow(rel_speed,2);
  double vth = sqrt(2.0*vtsq);
  double exp_term = exp(-pow(rel_speed/vth, 2));
  double erf_term = erf(rel_speed / vth);

  double term1 = (1.0/sqrt(GKYL_PI))*((rel_speedsq - 3.0*pow(rel_vel_in_dir, 2))/(rel_speedsq*rel_speedsq))*exp_term;
  double term2 = erf_term * ((3.0*vth*pow(rel_vel_in_dir, 2) - rel_speedsq*vth)/(2*pow(rel_speed, 5)));
  double term3 = erf_term * ((rel_speedsq - pow(rel_vel_in_dir, 2))/(vth*pow(rel_speed, 3)));
  double d2Gdvi2 = den*vth*(term1 + term2 + term3);
  return d2Gdvi2;
}

GKYL_CU_D
static inline double eval_fpo_d2gdv2_cross(double den,
  double rel_vel_in_dir1, double rel_vel_in_dir2, double rel_speed,
  double vtsq) {
  double rel_speedsq = pow(rel_speed,2);
  double vth = sqrt(2.0*vtsq);
  double expTerm = exp(-pow(rel_speed/vth, 2));
  double erfTerm = erf(rel_speed / vth);

  double term1 = (2.0/sqrt(GKYL_PI))*((-1.5*rel_vel_in_dir1*rel_vel_in_dir2)/(rel_speedsq*rel_speedsq))*expTerm;
  double term2 = ((1.5 * vth * rel_vel_in_dir1*rel_vel_in_dir2)/pow(rel_speed, 5))*erfTerm;
  double term3 = ((-1.0 * rel_vel_in_dir1*rel_vel_in_dir2)/(vth * pow(rel_speed, 3)))*erfTerm;

  double d2Gdvidvj = den*vth*(term1 + term2 + term3);
  return d2Gdvidvj;
}

void 
gkyl_proj_maxwellian_pots_on_basis_advance_cu(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* prim_moms,
  struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
  struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
  struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf,
  struct gkyl_array *fpo_d2gdv2_surf);
