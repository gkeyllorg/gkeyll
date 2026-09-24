#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_skin_surf_from_ghost_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_skin_surf_from_ghost*
gkyl_skin_surf_from_ghost_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis basis,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_skin_surf_from_ghost *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;

  // Choose the kernel that does the skin surf from ghost copy
  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_skin_surf_from_ghost_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_skin_surf_from_ghost_kernels));
#endif

  skin_surf_from_ghost_choose_kernel(basis, edge, up->dir, use_gpu, up->kernels);

  return up;
}

void
gkyl_skin_surf_from_ghost_advance(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    skin_surf_from_ghost_advance_cu(up, field);
    return;
  }
#endif

  int gidx[GKYL_MAX_DIM]; // ghost index.

  int ndim = up->skin_r->ndim; 

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->skin_r);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_copy_int_arr(ndim, iter.idx, gidx);
    // Get ghost cell corresponding to skin cell
    gidx[up->dir] = up->edge == GKYL_LOWER_EDGE? iter.idx[up->dir]-1 : iter.idx[up->dir]+1; 

    long ghost_linidx = gkyl_range_idx(up->ghost_r, gidx);
    long skin_linidx  = gkyl_range_idx(up->skin_r, iter.idx);

    const double *inp = (const double*) gkyl_array_cfetch(field, ghost_linidx);
    double *out = (double*) gkyl_array_fetch(field, skin_linidx);

    // Now call the kernel to ensure the nodal values to match the ghost nodal value 
    // while ensuring that the other nodal value remains unchanged
    up->kernels->ghost_to_skin(inp,out);

  }
}

void gkyl_skin_surf_from_ghost_release(struct gkyl_skin_surf_from_ghost *up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu)
    gkyl_free(up->kernels);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels);
#endif
  gkyl_free(up);
}