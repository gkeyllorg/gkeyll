#include <gkyl_translate_dim.h>
#include <gkyl_translate_dim_priv.h>
#include <gkyl_alloc.h>

struct gkyl_translate_dim*
gkyl_translate_dim_new(int cdim_do, struct gkyl_basis basis_do, int cdim_tar,
  struct gkyl_basis basis_tar, int dir, enum gkyl_edge_loc edge, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_translate_dim *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->cdim_do = cdim_do;
  up->cdim_tar = cdim_tar;
  up->vdim_do = basis_do.ndim - cdim_do;
  up->vdim_tar = basis_tar.ndim - cdim_tar;
  up->num_basis_do = basis_do.num_basis;
  up->num_basis_tar = basis_tar.num_basis;
  up->dir = dir;

  // Perform some basic checks:
  assert(up->vdim_do == up->vdim_tar);
  assert(basis_do.b_type == basis_tar.b_type);
  assert(basis_do.poly_order == basis_tar.poly_order);
  // Set pointer to the function performing basic checks on the ranges in advance method.
  if (up->vdim_do == 0) {
    up->range_check_func = cdim_do > cdim_tar? translate_dim_range_check_conf_deflate
                                             : translate_dim_range_check_conf_inflate;
  }
  else {
    up->range_check_func = cdim_do > cdim_tar? translate_dim_range_check_phase_deflate
                                             : translate_dim_range_check_phase_inflate;
  }

  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_translate_dim_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_translate_dim_kernels));
#endif

  // Choose kernels that translates the DG coefficients.
  trans_dim_choose_kernel(up->kernels, cdim_do, basis_do, cdim_tar, basis_tar, dir, edge, use_gpu);

  return up;
}

void
gkyl_translate_dim_advance(gkyl_translate_dim* up,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, int ncomp,
  struct gkyl_array *GKYL_RESTRICT ftar)
{
  // Perform some basic checks.
  up->range_check_func(up->dir, up->cdim_do, up->cdim_tar, up->vdim_do, rng_do, rng_tar);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_translate_dim_advance_cu(up, rng_do, rng_tar, fdo, ncomp, ftar);
    return;
  }
#endif

  int idx_do[GKYL_MAX_DIM] = {0};

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, rng_tar);
  while (gkyl_range_iter_next(&iter)) {

    // Translate the target idx to the donor idx:
    up->kernels->get_idx_do(up->cdim_tar, up->vdim_do, iter.idx, rng_do, up->cdim_do, idx_do, up->dir);

    long linidx_do = gkyl_range_idx(rng_do, idx_do);
    long linidx_tar = gkyl_range_idx(rng_tar, iter.idx);

    const double *fdo_c = gkyl_array_cfetch(fdo, linidx_do);
    double *ftar_c = gkyl_array_fetch(ftar, linidx_tar);

    for (int n=0; n<ncomp; n++) {
      up->kernels->translate(&fdo_c[n*up->num_basis_do], &ftar_c[n*up->num_basis_tar]);
    }
  }
}

void
gkyl_translate_dim_release(gkyl_translate_dim* up)
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
