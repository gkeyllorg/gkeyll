#include <gkyl_rz_calc_derived_geo.h>
#include <gkyl_rz_calc_derived_geo_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>

gkyl_rz_calc_derived_geo*
gkyl_rz_calc_derived_geo_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, int node_type, bool use_gpu)
{
  gkyl_rz_calc_derived_geo *up = gkyl_malloc(sizeof(gkyl_rz_calc_derived_geo));
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->cbasis = *cbasis;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->kernel = rz_derived_geo_choose_kernel(up->cdim, cbasis->b_type, node_type, up->poly_order);
  return up;
}


void
gkyl_rz_calc_derived_geo_advance(const gkyl_rz_calc_derived_geo *up, const struct gkyl_range *crange,
  struct gkyl_array *g_ij, struct gkyl_array *bmag, struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_inv,
  struct gkyl_array *gij, struct gkyl_array *b_i, struct gkyl_array *cmag, struct gkyl_array *jacobtot, struct gkyl_array *jacobtot_inv,
  struct gkyl_array *gxxj,  struct gkyl_array *gxyj, struct gkyl_array *gyyj, struct gkyl_array *gxzj, struct gkyl_array *eps2)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    const double *g_ij_i = gkyl_array_cfetch(g_ij, loc);
    const double *bmag_i = gkyl_array_cfetch(bmag, loc);
    const double *jacobgeo_i = gkyl_array_cfetch(jacobgeo, loc);
    double *jacobgeo_inv_i = gkyl_array_fetch(jacobgeo_inv, loc);
    double *gij_i = gkyl_array_fetch(gij, loc);
    double *bi_i = gkyl_array_fetch(b_i, loc);
    double *cmag_i = gkyl_array_fetch(cmag, loc);
    double *jacobtot_i = gkyl_array_fetch(jacobtot, loc);
    double *jacobtot_inv_i = gkyl_array_fetch(jacobtot_inv, loc);
    double *gxxj_i= gkyl_array_fetch(gxxj, loc);
    double *gxyj_i= gkyl_array_fetch(gxyj, loc);
    double *gyyj_i= gkyl_array_fetch(gyyj, loc);
    double *gxzj_i= gkyl_array_fetch(gxzj, loc);
    double *eps2_i= gkyl_array_fetch(eps2, loc);
    up->kernel(g_ij_i, bmag_i, jacobgeo_i, jacobgeo_inv_i, gij_i, bi_i, cmag_i, jacobtot_i, jacobtot_inv_i,
      gxxj_i, gxyj_i, gyyj_i, gxzj_i, eps2_i);
  }
}

void
gkyl_rz_calc_derived_geo_release(gkyl_rz_calc_derived_geo* up)
{
  gkyl_free(up);
}
