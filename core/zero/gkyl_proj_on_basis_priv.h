// Private header for the proj_on_basis updater. Not for direct use in user code.
#pragma once

#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_proj_on_basis {
  struct gkyl_rect_grid grid;
  int num_quad; // number of quadrature points to use in each direction
  int num_ret_vals; // number of values returned by eval function
  evalf_t eval; // function to project
  void *ctx; // evaluation context

  int num_basis; // number of basis functions
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // ordinates for quadrature
  struct gkyl_array *weights; // weights for quadrature
  struct gkyl_array *basis_at_ords; // basis functions at ordinates

  proj_on_basis_c2p_t c2p; // Function transformin comp to phys coords.
  void *c2p_ctx; // Context for the c2p mapping.

  bool use_gpu; // Whether the projection runs on the GPU.
  // Device copies of the quadrature data (only allocated when use_gpu=true).
  struct gkyl_array *ordinates_cu;
  struct gkyl_array *weights_cu;
  struct gkyl_array *basis_at_ords_cu;
  struct gkyl_proj_on_basis *on_dev; // Device clone of this updater (points to itself on CPU).
};

GKYL_CU_DH
static inline void
proj_on_basis_log_to_comp(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  // Convert logical to computational coordinates.
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

#ifdef GKYL_HAVE_CUDA

/**
 * Create the device clone of the updater. The host-side updater must already
 * have its quadrature device mirrors (ordinates_cu, weights_cu,
 * basis_at_ords_cu) allocated and populated, and its eval/c2p members must be
 * device function pointers (c2p may be NULL, selecting an identity mapping
 * assigned on the device).
 *
 * @param up Host-side updater.
 * @return Pointer to the device clone.
 */
struct gkyl_proj_on_basis* gkyl_proj_on_basis_cu_dev_new(struct gkyl_proj_on_basis *up);

/**
 * Run the projection on the GPU.
 *
 * @param up Project on basis updater (created with use_gpu=true).
 * @param tm Time at which projection must be computed.
 * @param update_range Range on which to run projection.
 * @param arr Output array (must be a device array).
 */
void gkyl_proj_on_basis_advance_cu(const struct gkyl_proj_on_basis *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr);

#endif
