#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_deflate_zsurf_priv.h>

struct gkyl_deflate_zsurf* 
gkyl_deflate_zsurf_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis, 
  int edge, bool use_gpu) 
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_deflate_zsurf_cu_dev_new(cbasis, deflated_cbasis, edge);
  } 
#endif 

  gkyl_deflate_zsurf *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;
  up->num_basis = cbasis->num_basis;
  up->num_deflated_basis = deflated_cbasis->num_basis;
  up->cdim = cbasis->ndim;
  up->kernel = deflate_zsurf_choose_kernel(cbasis->b_type, cbasis->ndim, edge, cbasis->poly_order); // edge = 0,1 = lo, up

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}


void 
gkyl_deflate_zsurf_advance(const gkyl_deflate_zsurf *up, int zidx, 
  const struct gkyl_range *range, const struct gkyl_range *deflated_range, 
  const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp) 
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(deflated_field)) {
    return gkyl_deflate_zsurf_advance_cu(up, zidx, range, deflated_range, 
      field, deflated_field, ncomp);
  }
#endif
  int do_idx[3];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, deflated_range);

  while (gkyl_range_iter_next(&iter)) {
    for(int i = 0; i < up->cdim-1; i++)
      do_idx[i] = iter.idx[i];
    do_idx[up->cdim-1] = zidx;

    long loc = gkyl_range_idx(range, do_idx);
    const double *fld = gkyl_array_cfetch(field, loc);

    long loc_deflated = gkyl_range_idx(deflated_range, iter.idx);
    double *fld_deflated = gkyl_array_fetch(deflated_field, loc_deflated);
    for(int c = 0; c<ncomp; c++)
      up->kernel(&fld[c*up->num_basis], &fld_deflated[c*up->num_deflated_basis]);
  }
}

void gkyl_deflate_zsurf_release(gkyl_deflate_zsurf* up) 
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  gkyl_free(up);
}
