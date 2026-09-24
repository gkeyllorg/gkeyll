/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_priv.h>
}

// CUDA kernel to set device pointers to kernel that computes the reflected f.
__global__ static void
gkyl_bc_gksheath_set_cu_ker_ptrs(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int dim = basis->ndim;
  enum gkyl_basis_type b_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->reflectedf = ser_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gkyl_bc_gksheath_choose_reflectedf_kernel_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  gkyl_bc_gksheath_set_cu_ker_ptrs<<<1,1>>>(basis, edge, kers);
}

__global__ static void
gkyl_bc_sheath_gyrokinetic_advance_cu_ker(int cdim, int dir, const struct gkyl_range skin_r, const struct gkyl_range ghost_r,
  const struct gkyl_range conf_r, const struct gkyl_range vel_r, const struct gkyl_basis *basis,
  const struct gkyl_array *vmap, double q2Dm, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_bc_sheath_gyrokinetic_kernels *kers, struct gkyl_array *distf)
{
  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int pidx[GKYL_MAX_DIM];
  int vidx[2];

  int pdim = skin_r.ndim;
  int vpar_dir = cdim;
  int uplo = skin_r.upper[vpar_dir]+skin_r.lower[vpar_dir];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < skin_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange

    gkyl_sub_range_inv_idx(&skin_r, linc, pidx);

    gkyl_copy_int_arr(pdim, pidx, fidx);
    fidx[vpar_dir] = uplo - pidx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[dir] = ghost_r.lower[dir];

    long skin_loc = gkyl_range_idx(&skin_r, pidx);
    long ghost_loc = gkyl_range_idx(&ghost_r, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(distf, skin_loc);
    double *out = (double*) gkyl_array_fetch(distf, ghost_loc);

    for (int d=cdim; d<pdim; d++) vidx[d-cdim] = pidx[d]; 
    long conf_loc = gkyl_range_idx(&conf_r, pidx);
    long vel_loc = gkyl_range_idx(&vel_r, vidx);

    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *vmap_p = (const double*) gkyl_array_cfetch(vmap, vel_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[112];  // MF 2022/08/24: hardcoded to number of DG coeffs in 3x2v p2 for now.
    kers->reflectedf(vmap_p, q2Dm, phi_p, phi_wall_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(dir, basis, cdim, out, fhat);
  }
}

void
gkyl_bc_sheath_gyrokinetic_advance_cu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r)
{
  if (up->skin_r->volume > 0) {
    int nblocks = up->skin_r->nblocks, nthreads = up->skin_r->nthreads;

    gkyl_bc_sheath_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, *up->skin_r, *up->ghost_r,
      *conf_r, up->vel_map->local_vel, up->basis, up->vel_map->vmap->on_dev, up->q2Dm, phi->on_dev, phi_wall->on_dev,
      up->kernels_cu, distf->on_dev);
  }
}
