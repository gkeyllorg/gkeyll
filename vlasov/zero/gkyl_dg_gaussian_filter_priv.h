// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h>

GKYL_CU_DH static inline void
comp_to_phys(
  int ndim, const double *eta, const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xc,
  double *GKYL_RESTRICT xout
)
{
  for (int d = 0; d < ndim; ++d) {
    xout[d] = 0.5 * dx[d] * eta[d] + xc[d];
  }
}

GKYL_CU_DH static inline double
filter_norm(
  int cdim, int tot_quad, const struct gkyl_array *ordinates, const struct gkyl_array *weights,
  const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xmu, const double *GKYL_RESTRICT xcl,
  double *GKYL_RESTRICT xmul, const double *GKYL_RESTRICT xcc, double *GKYL_RESTRICT xmuc,
  const double *GKYL_RESTRICT xcr, double *GKYL_RESTRICT xmur
)
{
  double sum = 0.0;
  for (int j = 0; j < tot_quad; ++j) {
    // Fetch the quadrature point locations for all other quadrature points.
    comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, j), dx, xcl, xmul);
    comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, j), dx, xcc, xmuc);
    comp_to_phys(cdim, (const double *)gkyl_array_cfetch(ordinates, j), dx, xcr, xmur);
    const double *wq = (const double *)gkyl_array_cfetch(weights, j);
    // Sum the contributions from each cell in the three-cell filter.
    double expamp_l = -pow((xmu[0] - xmul[0]) / (sqrt(2.0) * dx[0]), 2.0);
    double expamp_c = -pow((xmu[0] - xmuc[0]) / (sqrt(2.0) * dx[0]), 2.0);
    double expamp_r = -pow((xmu[0] - xmur[0]) / (sqrt(2.0) * dx[0]), 2.0);
    sum += wq[0] * (exp(expamp_l) + exp(expamp_c) + exp(expamp_r));
  }
  return sum;
}

struct gkyl_dg_gaussian_filter {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid
  int cdim; // Configuration-space dimension

  struct gkyl_basis conf_basis; // Configuration-space basis
  int num_conf_basis; // number of configuration-space basis functions

  // for quadrature in configuration-space
  int tot_quad; // total number of configuration-space quadrature points
  struct gkyl_array *ordinates; // Configuration-space ordinates for quadrature
  struct gkyl_array *weights; // weights for configuration-space quadrature
  struct gkyl_array *basis_at_ords; // Configuration-space basis functions at ordinates

  struct gkyl_mat_mm_array_mem *n2m_mem; // structure of data which stores the info to
    // convert nodal to modal gkyl arrays
  struct gkyl_mat_mm_array_mem *m2n_mem; // structure of data which stores the info to
    // convert modal to nodal gkyl arrays

  struct gkyl_array *conf_arr_quad; // Configuration-space array at quadrature points.
  struct gkyl_array *conf_arr_filter; // Filtered configuration-space array at quadrature points.

  bool use_gpu; // Boolean if we are filtering on device.
  bool extend_filter; // Boolean if we are extending filter to 5 cells.
};
