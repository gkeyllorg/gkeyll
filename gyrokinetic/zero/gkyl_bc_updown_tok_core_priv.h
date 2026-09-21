#pragma once

// Private header for bc_updown_tok_core updater, not for direct use in user code.

#include <gkyl_bc_updown_tok_core.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_util.h>
#include <assert.h>

// Max number of DG coefficients per cell supported (3x2v p2 serendipity).
#define BC_UPDOWN_TOK_CORE_MAX_NUM_BASIS 112

// Parameters of the mirror map (passed by value to device kernels).
struct bc_updown_tok_core_params {
  int dir; // Direction of the BC (x).
  int cdim; // Conf-space dimensionality.
  int zdir; // z direction (cdim-1).
  int ydir; // y direction (1), or -1 if there is no y direction.
  int uplo_z; // lower+upper of the local skin range in z.
  int ylo, yup; // Extents of the skin range in y.
  bool avg_y; // Average the donor over y.
  int num_basis; // Number of DG coefficients per cell.
};

// Primary struct in this updater.
struct gkyl_bc_updown_tok_core {
  struct bc_updown_tok_core_params par; // Mirror map parameters.
  const struct gkyl_basis *basis; // Phase-space basis (on device if use_gpu).
  const struct gkyl_range *skin_r, *ghost_r; // Skin and ghost ranges.
  struct gkyl_range buff_r; // Range of the send/recv buffers (skin with y collapsed).
  struct gkyl_array *buff_send, *buff_recv; // Send/recv buffers (distinct, needed for the same-rank path).
  struct gkyl_comm *comm; // Communicator.
  int partner; // Rank owning the mirrored z extent.
  struct gkyl_multib_comm_conn *mbcc_send, *mbcc_recv; // Connection to the partner rank.
  int local_blocks[1]; // Block ID list for the transfer (single block).
  bool use_gpu; // Whether to run on GPU.
};

// Fill the buffer cell bidx with the z-mirrored, x-reflected and y-averaged
// skin cell(s).
GKYL_CU_DH
static inline void
bc_updown_tok_core_pack_cell(const struct bc_updown_tok_core_params *par, const struct gkyl_basis *basis,
  const struct gkyl_range *skin_r, const int *bidx, const struct gkyl_array *distf, double *out)
{
  int pdim = skin_r->ndim;
  int sidx[GKYL_MAX_DIM];
  for (int d=0; d<pdim; d++) sidx[d] = bidx[d];
  sidx[par->zdir] = par->uplo_z - bidx[par->zdir];

  int ylo = 0, yup = 0;
  bool do_yavg = par->avg_y && (par->ydir >= 0);
  if (do_yavg) {
    ylo = par->ylo;
    yup = par->yup;
  }
  double fac = 1.0/(yup-ylo+1);

  for (int k=0; k<par->num_basis; k++) out[k] = 0.0;

  double tmp[BC_UPDOWN_TOK_CORE_MAX_NUM_BASIS], tmp_y[BC_UPDOWN_TOK_CORE_MAX_NUM_BASIS];
  for (int iy=ylo; iy<=yup; iy++) {
    if (do_yavg) sidx[par->ydir] = iy;

    long skin_loc = gkyl_range_idx(skin_r, sidx);
    const double *inp = (const double*) gkyl_array_cfetch(distf, skin_loc);

    // z -> -z.
    basis->flip_odd_sign(par->zdir, inp, tmp);
    // Reflect in x so the ghost trace at the face equals the donor trace.
    basis->flip_odd_sign(par->dir, tmp, tmp);

    if (do_yavg) {
      // Keep only the y-even part (in-cell y average), and average over cells.
      basis->flip_odd_sign(par->ydir, tmp, tmp_y);
      for (int k=0; k<par->num_basis; k++) out[k] += fac*0.5*(tmp[k] + tmp_y[k]);
    }
    else {
      for (int k=0; k<par->num_basis; k++) out[k] += tmp[k];
    }
  }
}

// Copy the buffer cell corresponding to the ghost cell gidx into the ghost.
GKYL_CU_DH
static inline void
bc_updown_tok_core_unpack_cell(const struct bc_updown_tok_core_params *par, const struct gkyl_range *buff_r,
  const struct gkyl_range *ghost_r, const int *gidx, const struct gkyl_array *buff, struct gkyl_array *distf)
{
  int pdim = ghost_r->ndim;
  int bidx[GKYL_MAX_DIM];
  for (int d=0; d<pdim; d++) bidx[d] = gidx[d];
  bidx[par->dir] = buff_r->lower[par->dir];
  if (par->avg_y && (par->ydir >= 0))
    bidx[par->ydir] = buff_r->lower[par->ydir];

  long buff_loc = gkyl_range_idx(buff_r, bidx);
  long ghost_loc = gkyl_range_idx(ghost_r, gidx);

  const double *inp = (const double*) gkyl_array_cfetch(buff, buff_loc);
  double *out = (double*) gkyl_array_fetch(distf, ghost_loc);
  for (int k=0; k<par->num_basis; k++) out[k] = inp[k];
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to pack the mirrored skin into the send buffer.
 *
 * @param up BC updater.
 * @param distf Distribution function.
 */
void gkyl_bc_updown_tok_core_pack_cu(struct gkyl_bc_updown_tok_core *up, const struct gkyl_array *distf);

/**
 * CUDA device function to unpack the received buffer into the ghost cells.
 *
 * @param up BC updater.
 * @param distf Distribution function.
 */
void gkyl_bc_updown_tok_core_unpack_cu(struct gkyl_bc_updown_tok_core *up, struct gkyl_array *distf);

#endif
