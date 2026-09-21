/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_updown_tok_core.h>
#include <gkyl_bc_updown_tok_core_priv.h>
}

// One thread per buffer cell, the y average is a serial inner loop (no race).
__global__ static void
gkyl_bc_updown_tok_core_pack_cu_ker(struct bc_updown_tok_core_params par,
  const struct gkyl_basis *basis, const struct gkyl_range skin_r, const struct gkyl_range buff_r,
  const struct gkyl_array *distf, struct gkyl_array *buff_send)
{
  int bidx[GKYL_MAX_DIM];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < buff_r.volume; linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&buff_r, linc, bidx);

    long buff_loc = gkyl_range_idx(&buff_r, bidx);
    double *out = (double*) gkyl_array_fetch(buff_send, buff_loc);
    bc_updown_tok_core_pack_cell(&par, basis, &skin_r, bidx, distf, out);
  }
}

__global__ static void
gkyl_bc_updown_tok_core_unpack_cu_ker(struct bc_updown_tok_core_params par,
  const struct gkyl_range buff_r, const struct gkyl_range ghost_r,
  const struct gkyl_array *buff_recv, struct gkyl_array *distf)
{
  int gidx[GKYL_MAX_DIM];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < ghost_r.volume; linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&ghost_r, linc, gidx);

    bc_updown_tok_core_unpack_cell(&par, &buff_r, &ghost_r, gidx, buff_recv, distf);
  }
}

void
gkyl_bc_updown_tok_core_pack_cu(struct gkyl_bc_updown_tok_core *up, const struct gkyl_array *distf)
{
  if (up->buff_r.volume > 0) {
    int nblocks = up->buff_r.nblocks, nthreads = up->buff_r.nthreads;
    gkyl_bc_updown_tok_core_pack_cu_ker<<<nblocks, nthreads>>>(up->par, up->basis,
      *up->skin_r, up->buff_r, distf->on_dev, up->buff_send->on_dev);
  }
}

void
gkyl_bc_updown_tok_core_unpack_cu(struct gkyl_bc_updown_tok_core *up, struct gkyl_array *distf)
{
  if (up->ghost_r->volume > 0) {
    int nblocks = up->ghost_r->nblocks, nthreads = up->ghost_r->nthreads;
    gkyl_bc_updown_tok_core_unpack_cu_ker<<<nblocks, nthreads>>>(up->par, up->buff_r,
      *up->ghost_r, up->buff_recv->on_dev, distf->on_dev);
  }
}
