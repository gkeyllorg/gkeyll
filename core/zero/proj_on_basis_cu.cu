/* -*- c++ -*- */

extern "C" {
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_on_basis_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}

// Identity comp to phys coord mapping, for when user doesn't provide a map.
GKYL_CU_D static void
proj_on_basis_c2p_identity_cu(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_rect_grid *grid = (struct gkyl_rect_grid *) ctx;
  int ndim = grid->ndim;
  for (int d=0; d<ndim; d++) xphys[d] = xcomp[d];
}

__global__ static void
proj_on_basis_set_c2p_identity_cu_ker(struct gkyl_proj_on_basis *up)
{
  // Assigned in device code so the function pointer is a device address.
  up->c2p = proj_on_basis_c2p_identity_cu;
  up->c2p_ctx = &up->grid; // Use grid as the context since all we need is ndim.
}

struct gkyl_proj_on_basis*
gkyl_proj_on_basis_cu_dev_new(struct gkyl_proj_on_basis *up)
{
  struct gkyl_proj_on_basis *up_cu = (struct gkyl_proj_on_basis*)
    gkyl_cu_malloc(sizeof(struct gkyl_proj_on_basis));

  // Clone with device pointers to the quadrature data. The eval/c2p members
  // and their contexts are device pointers provided by the user.
  struct gkyl_proj_on_basis up_ho = *up;
  up_ho.ordinates = up->ordinates_cu->on_dev;
  up_ho.weights = up->weights_cu->on_dev;
  up_ho.basis_at_ords = up->basis_at_ords_cu->on_dev;
  up_ho.on_dev = up_cu;

  gkyl_cu_memcpy(up_cu, &up_ho, sizeof(struct gkyl_proj_on_basis), GKYL_CU_MEMCPY_H2D);

  // When the user does not provide a c2p mapping use the identity mapping,
  // which must be assigned on the device to obtain a device function pointer.
  if (up->c2p == 0)
    proj_on_basis_set_c2p_identity_cu_ker<<<1,1>>>(up_cu);

  return up_cu;
}

__global__ static void
gkyl_proj_on_basis_advance_cu_ker(const struct gkyl_proj_on_basis *up, double tm,
  struct gkyl_range update_range, struct gkyl_array* GKYL_RESTRICT fun_at_ords,
  struct gkyl_array* GKYL_RESTRICT arr)
{
  int idx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_basis = up->num_basis;
  int tot_quad = up->tot_quad;
  int num_ret_vals = up->num_ret_vals;
  int ndim = up->grid.ndim;

  const double* GKYL_RESTRICT weights = (const double*) up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = (const double*) up->basis_at_ords->data;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < update_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&update_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->grid, idx, xc);

    // Scratch space for the function evaluated at one quadrature node,
    // one slot per cell (i.e. per linc1).
    double *fq = (double*) gkyl_array_fetch(fun_at_ords, linc1);

    long lidx = gkyl_range_idx(&update_range, idx);
    double *f = (double*) gkyl_array_fetch(arr, lidx);

    // Arrangement of f is as:
    // c0[0], c0[1], ... c1[0], c1[1], ....
    // where c0, c1, ... are components of f (num_ret_vals).
    for (int k=0; k<num_ret_vals*num_basis; ++k) f[k] = 0.0;

    for (int imu=0; imu<tot_quad; ++imu) {
      proj_on_basis_log_to_comp(ndim, (const double*) gkyl_array_cfetch(up->ordinates, imu),
        up->grid.dx, xc, xmu);
      up->c2p(xmu, xmu, up->c2p_ctx);
      up->eval(tm, xmu, fq, up->ctx);

      long offset = 0;
      for (int n=0; n<num_ret_vals; ++n) {
        double tmp = weights[imu]*fq[n];
        for (int k=0; k<num_basis; ++k)
          f[offset+k] += tmp*basis_at_ords[k+num_basis*imu];
        offset += num_basis;
      }
    }
  }
}

void
gkyl_proj_on_basis_advance_cu(const struct gkyl_proj_on_basis *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr)
{
  assert(gkyl_array_is_cu_dev(arr));

  // Scratch storage for the function evaluated at a quadrature node, one slot
  // per cell (allocated per call, like the CPU path allocates fun_at_ords).
  struct gkyl_array *fun_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    up->num_ret_vals, update_range->volume);

  gkyl_proj_on_basis_advance_cu_ker<<<update_range->nblocks, update_range->nthreads>>>(
    up->on_dev, tm, *update_range, fun_at_ords->on_dev, arr->on_dev);

  gkyl_array_release(fun_at_ords);
}
