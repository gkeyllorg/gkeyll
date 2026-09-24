/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_ambi_bolt_potential_priv.h>
}

// CUDA kernel to set device pointers to l2g, RHS src and solution
// kernels. Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
ambi_bolt_potential_set_cu_ker_ptrs(struct gkyl_ambi_bolt_potential_kernels* kers,
  enum gkyl_basis_type b_type, int dim, int poly_order)
{
  const sheath_calc_kern_edge_list *sheath_calc_list;
  const phi_calc_kern_list *phi_calc_list;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      sheath_calc_list = ser_sheath_calc_list;
      phi_calc_list = ser_phi_calc_list;
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
  for (int k=0; k<2; k++)
    kers->sheath_calc[k] = CSHEATHK(sheath_calc_list, dim, poly_order, k);

  kers->phi_calc = CPHIK(phi_calc_list, dim, poly_order);;
}

__global__ static void
gkyl_ambi_bolt_potential_sheath_calc_cu_ker(double dz, double charge_e, double mass_e, double temp_e,
  struct gkyl_ambi_bolt_potential_kernels *kers, enum gkyl_edge_loc edge,
  struct gkyl_range skin_r, struct gkyl_range ghost_r,
  const struct gkyl_array *cmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *gammai, const struct gkyl_array *m0i, const struct gkyl_array *Jm0i,
  struct gkyl_array *sheath_vals)
{
  unsigned int keridx = (edge == GKYL_LOWER_EDGE) ? 0 : 1;

  int idx_g[GKYL_MAX_CDIM], idx_s[GKYL_MAX_CDIM]; // Ghost and skin indices.

  int ndim = ghost_r.ndim;

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < ghost_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc to idx
    // must use gkyl_sub_range_inv_idx so that linc=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&ghost_r, linc, idx_g);

    gkyl_copy_int_arr(ndim, idx_g, idx_s);
    // Assume only 1 ghost cell on either side along the field line.
    idx_s[ndim-1] = edge == GKYL_LOWER_EDGE? idx_g[ndim-1]+1 : idx_g[ndim-1]-1;

    long ghost_loc = gkyl_range_idx(&ghost_r, idx_g);
    long skin_loc = gkyl_range_idx(&skin_r, idx_s);

    const double *cmag_p = (const double*) gkyl_array_cfetch(cmag, skin_loc);
    const double *jactotinv_p = (const double*) gkyl_array_cfetch(jacobtot_inv, skin_loc);
    const double *m0i_p = (const double*) gkyl_array_cfetch(m0i, skin_loc);
    const double *Jm0i_p = (const double*) gkyl_array_cfetch(Jm0i, skin_loc);
    const double *gammai_p = (const double*) gkyl_array_cfetch(gammai, ghost_loc);
    double *out_p = (double*) gkyl_array_cfetch(sheath_vals, ghost_loc);

    kers->sheath_calc[keridx](dz, charge_e, mass_e, temp_e, cmag_p, jactotinv_p, gammai_p, m0i_p, Jm0i_p, out_p);
  }
}

__global__ static void
gkyl_ambi_bolt_potential_phi_calc_cu_ker(double charge_e, double temp_e,
  struct gkyl_ambi_bolt_potential_kernels *kers, struct gkyl_range local_r, struct gkyl_range extlocal_r,
  const struct gkyl_array *m0i, const struct gkyl_array *sheath_vals,
  struct gkyl_array *phi)
{
  int idx[GKYL_MAX_CDIM], idx_g[GKYL_MAX_CDIM]; // Volume and ghost indices.

  int ndim = local_r.ndim;

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < local_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc to idx
    // must use gkyl_sub_range_inv_idx so that linc=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&local_r, linc, idx);

    // We assume each MPI rank calls this over the local range and
    // that the sheath values are defined on the lower local ghost range.
    gkyl_copy_int_arr(ndim, idx, idx_g);
    idx_g[ndim-1] = extlocal_r.lower[ndim-1];

    long loc = gkyl_range_idx(&local_r, idx);
    long ghost_loc = gkyl_range_idx(&extlocal_r, idx_g);

    const double *m0i_p = (const double*) gkyl_array_cfetch(m0i, loc);
    const double *sheathvals_p = (const double*) gkyl_array_cfetch(sheath_vals, ghost_loc);
    double *phi_p = (double*) gkyl_array_cfetch(phi, loc);

    kers->phi_calc(charge_e, temp_e, m0i_p, sheathvals_p, phi_p);
  }
}

void
ambi_bolt_potential_choose_kernels_cu(const struct gkyl_basis *basis, struct gkyl_ambi_bolt_potential_kernels *kers)
{
  ambi_bolt_potential_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, basis->ndim, basis->poly_order);
}

void
gkyl_ambi_bolt_potential_sheath_calc_cu(struct gkyl_ambi_bolt_potential *up, enum gkyl_edge_loc edge,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_array *cmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *gammai, const struct gkyl_array *m0i, const struct gkyl_array *Jm0i,
  struct gkyl_array *sheath_vals)
{
  int nblocks = ghost_r->nblocks, nthreads = ghost_r->nthreads;

  gkyl_ambi_bolt_potential_sheath_calc_cu_ker<<<nblocks, nthreads>>>(up->dz, up->charge_e, up->mass_e, up->temp_e,
    up->kernels_cu, edge, *skin_r, *ghost_r, cmag->on_dev, jacobtot_inv->on_dev,
    gammai->on_dev, m0i->on_dev, Jm0i->on_dev, sheath_vals->on_dev);
}

void
gkyl_ambi_bolt_potential_phi_calc_cu(struct gkyl_ambi_bolt_potential *up, const struct gkyl_range *local_r,
  const struct gkyl_range *extlocal_r, const struct gkyl_array *m0i, const struct gkyl_array *sheath_vals,
  struct gkyl_array *phi)
{
  int nblocks = local_r->nblocks, nthreads = local_r->nthreads;

  gkyl_ambi_bolt_potential_phi_calc_cu_ker<<<nblocks, nthreads>>>(up->charge_e, up->temp_e, up->kernels_cu,
    *local_r, *extlocal_r, m0i->on_dev, sheath_vals->on_dev, phi->on_dev);
}
