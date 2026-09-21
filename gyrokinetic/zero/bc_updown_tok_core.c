#include <gkyl_bc_updown_tok_core.h>
#include <gkyl_bc_updown_tok_core_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>
#include <math.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct gkyl_bc_updown_tok_core*
gkyl_bc_updown_tok_core_new(const struct gkyl_bc_updown_tok_core_inp *inp)
{
  assert(inp->dir == 0); // Only the lower radial boundary is supported.
  assert(inp->edge == GKYL_LOWER_EDGE);
  assert(inp->cdim >= 2); // Need a z direction distinct from x.

  int cdim = inp->cdim;
  int zdir = cdim-1;

  // The z grid must be symmetric about z=0.
  const struct gkyl_rect_grid *grid = inp->grid;
  assert(fabs(grid->lower[zdir]+grid->upper[zdir]) < 1e-10*(grid->upper[zdir]-grid->lower[zdir]));

  const struct gkyl_range *skin_r = inp->skin_r, *ghost_r = inp->ghost_r;
  int pdim = skin_r->ndim;
  assert(ghost_r->ndim == pdim);
  for (int d=0; d<pdim; d++)
    assert(gkyl_range_shape(skin_r, d) == gkyl_range_shape(ghost_r, d));

  struct gkyl_bc_updown_tok_core *up = gkyl_malloc(sizeof(*up));

  up->basis = inp->basis;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;
  up->use_gpu = inp->use_gpu;

  struct bc_updown_tok_core_params *par = &up->par;
  par->dir = inp->dir;
  par->cdim = cdim;
  par->zdir = zdir;
  par->ydir = cdim == 3 ? 1 : -1;
  par->uplo_z = skin_r->lower[zdir] + skin_r->upper[zdir];
  par->ylo = par->ydir >= 0 ? skin_r->lower[par->ydir] : 0;
  par->yup = par->ydir >= 0 ? skin_r->upper[par->ydir] : 0;
  par->avg_y = inp->avg_y && (par->ydir >= 0);

  // Number of basis functions (the basis may live on the device).
  int num_basis;
  if (inp->use_gpu) {
    struct gkyl_basis basis_ho;
    gkyl_cu_memcpy(&basis_ho, (void*) inp->basis, sizeof(struct gkyl_basis), GKYL_CU_MEMCPY_D2H);
    num_basis = basis_ho.num_basis;
  }
  else {
    num_basis = inp->basis->num_basis;
  }
  assert(num_basis <= BC_UPDOWN_TOK_CORE_MAX_NUM_BASIS);
  par->num_basis = num_basis;

  // Buffer range: the skin range with y collapsed to a single cell.
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d=0; d<pdim; d++) {
    lower[d] = skin_r->lower[d];
    upper[d] = skin_r->upper[d];
  }
  if (par->avg_y)
    upper[par->ydir] = lower[par->ydir];
  gkyl_range_init(&up->buff_r, pdim, lower, upper);

  up->buff_send = mkarr(up->use_gpu, num_basis, up->buff_r.volume);
  up->buff_recv = mkarr(up->use_gpu, num_basis, up->buff_r.volume);

  // Find the rank owning the mirrored z extent (it also has the same x and y
  // extents since only z may be decomposed). In serial, or for the middle
  // rank of an odd number of z cuts, the partner is the rank itself.
  up->comm = gkyl_comm_acquire(inp->comm);
  int my_rank;
  gkyl_comm_get_rank(up->comm, &my_rank);
  const struct gkyl_rect_decomp *decomp = inp->decomp;
  const struct gkyl_range *local_conf = &decomp->ranges[my_rank];
  int glo_z = decomp->parent_range.lower[zdir], gup_z = decomp->parent_range.upper[zdir];
  int mirror_lo_z = glo_z + gup_z - local_conf->upper[zdir];
  int mirror_up_z = glo_z + gup_z - local_conf->lower[zdir];
  up->partner = -1;
  for (int r=0; r<decomp->ndecomp; r++) {
    const struct gkyl_range *rr = &decomp->ranges[r];
    bool match = (rr->lower[zdir] == mirror_lo_z) && (rr->upper[zdir] == mirror_up_z);
    for (int d=0; d<cdim-1; d++)
      match = match && (rr->lower[d] == local_conf->lower[d]) && (rr->upper[d] == local_conf->upper[d]);
    if (match) {
      up->partner = r;
      break;
    }
  }
  assert(up->partner >= 0); // Requires a uniform decomposition in z (and no x/y cuts).

  struct gkyl_comm_conn cc_send = {
    .sr = GKYL_COMM_CONN_SEND,
    .block_id = 0,
    .rank = up->partner,
    .range = up->buff_r,
    .src_edge = GKYL_LOWER_POSITIVE,
    .tar_edge = GKYL_LOWER_POSITIVE,
  };
  struct gkyl_comm_conn cc_recv = cc_send;
  cc_recv.sr = GKYL_COMM_CONN_RECV;
  up->mbcc_send = gkyl_multib_comm_conn_new(1, &cc_send);
  up->mbcc_recv = gkyl_multib_comm_conn_new(1, &cc_recv);
  up->local_blocks[0] = 0;

  return up;
}

// Pack the mirrored, y-averaged skin cells into the send buffer.
static void
bc_updown_tok_core_pack(struct gkyl_bc_updown_tok_core *up, const struct gkyl_array *distf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_updown_tok_core_pack_cu(up, distf);
    return;
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->buff_r);
  while (gkyl_range_iter_next(&iter)) {
    long buff_loc = gkyl_range_idx(&up->buff_r, iter.idx);
    double *out = (double*) gkyl_array_fetch(up->buff_send, buff_loc);
    bc_updown_tok_core_pack_cell(&up->par, up->basis, up->skin_r, iter.idx, distf, out);
  }
}

// Unpack the received buffer into the ghost cells.
static void
bc_updown_tok_core_unpack(struct gkyl_bc_updown_tok_core *up, struct gkyl_array *distf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_updown_tok_core_unpack_cu(up, distf);
    return;
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    bc_updown_tok_core_unpack_cell(&up->par, &up->buff_r, up->ghost_r, iter.idx, up->buff_recv, distf);
  }
}

void
gkyl_bc_updown_tok_core_advance(struct gkyl_bc_updown_tok_core *up, struct gkyl_array *distf)
{
  bc_updown_tok_core_pack(up, distf);

  gkyl_multib_comm_conn_array_transfer(up->comm, 1, up->local_blocks,
    &up->mbcc_send, &up->mbcc_recv, &up->buff_send, &up->buff_recv);

  bc_updown_tok_core_unpack(up, distf);
}

void
gkyl_bc_updown_tok_core_release(struct gkyl_bc_updown_tok_core *up)
{
  gkyl_multib_comm_conn_release(up->mbcc_send);
  gkyl_multib_comm_conn_release(up->mbcc_recv);
  gkyl_comm_release(up->comm);
  gkyl_array_release(up->buff_send);
  gkyl_array_release(up->buff_recv);
  gkyl_free(up);
}
