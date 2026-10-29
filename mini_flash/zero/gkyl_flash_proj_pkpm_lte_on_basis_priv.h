#pragma once

#include <string.h>

#include <gkyl_basis.h>
#include <gkyl_flash_proj_pkpm_lte_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct gkyl_flash_proj_pkpm_lte_on_basis {
  struct gkyl_rect_grid phase_grid;
  int cdim; // Configuration-space dimension
  int pdim; // Phase-space dimension
  int vdim; // Velocity-space dimension

  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // number of Configuration-space basis functions
  int num_phase_basis; // number of Phase-space basis functions

  struct gkyl_range conf_qrange; // Range of Configuration-space ordinates.
  struct gkyl_range phase_qrange; // Range of Phase-space ordinates.

  // for quadrature in Phase-space
  int tot_quad; // total number of Phase-space quadrature points
  struct gkyl_array *ordinates; // Phase-space ordinates for quadrature
  struct gkyl_array *weights; // Phase-space weights for quadrature
  struct gkyl_array *basis_at_ords; // Phase-space basis functions at ordinates

  // for quadrature in Configuration-space
  int tot_conf_quad; // total number of Configuration-space quadrature points
  struct gkyl_array *conf_basis_at_ords; // Configuration-space basis functions at ordinates

  struct gkyl_array *f0_at_ords; // function LTE distribution evaluated at
                                 // ordinates in a cell.
  struct gkyl_array *g_at_ords; // function LTE distribution evaluated at
                                 // ordinates in a cell.
};

GKYL_CU_DH
static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

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

// Number of quadrature points in velocity space, and which velocity directions
// use it. Hybrid basis have p=2 in velocity space, so they need one more point.
static inline int
vel_quad_shape(const struct gkyl_basis *phase_basis, int vdim, int num_quad, bool is_vdim_p2[3])
{
  bool is_hybrid = phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID;
  for (int d=0; d<3; ++d) is_vdim_p2[d] = is_hybrid && d < vdim;  // 3 is the max vdim.
  return is_hybrid ? num_quad+1 : num_quad;
}

// Sets ordinates, weights and basis functions at ords.
// Returns the total number of quadrature nodes
static int
init_quad_values(int cdim, const struct gkyl_basis *basis,
  enum gkyl_quad_type quad_type, int num_quad, int num_quad_v, bool *is_vdim_p2,
  struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords)
{
  int ndim = basis->ndim;

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];

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

  struct gkyl_range qrange = get_qrange(cdim, ndim, num_quad, num_quad_v, is_vdim_p2);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(*ordinates, node);
    for (int i=0; i<cdim; ++i) {
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    for (int i=cdim; i<ndim; ++i) {
      ord[i] = is_vdim_p2[i-cdim] ? 
        ordinates1_v[iter.idx[i]-qrange.lower[i]] : ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    
    // set weights
    double *wgt = gkyl_array_fetch(*weights, node);
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
  *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n) {
    basis->eval(gkyl_array_fetch(*ordinates, n), gkyl_array_fetch(*basis_at_ords, n));
  }

  return tot_quad;
}

static void
proj_on_basis(const gkyl_flash_proj_pkpm_lte_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
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