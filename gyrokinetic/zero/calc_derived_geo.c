#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_derived_geo_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>

gkyl_calc_derived_geo*
gkyl_calc_derived_geo_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, int node_type, bool use_gpu)
{
  gkyl_calc_derived_geo *up = gkyl_malloc(sizeof(*up));
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->cbasis = *cbasis;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->kernel = derived_geo_choose_kernel(up->cdim, cbasis->b_type, node_type, up->poly_order);
  return up;
}


void
gkyl_calc_derived_geo_advance(const gkyl_calc_derived_geo *up, const struct gkyl_range *crange, struct gkyl_array *gFld,
  struct gkyl_array *bmagFld, struct gkyl_array *jFld, struct gkyl_array *jinvFld, struct gkyl_array *grFld,
  struct gkyl_array *biFld, struct gkyl_array *cmagFld, struct gkyl_array *jtotFld, struct gkyl_array *jtotinvFld,
  struct gkyl_array *gxxJFld,  struct gkyl_array *gxyJFld, struct gkyl_array *gyyJFld, struct gkyl_array *gxzJFld,
  struct gkyl_array *eps2Fld)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    const double *gij = gkyl_array_cfetch(gFld, loc);
    const double *bmag_i = gkyl_array_cfetch(bmagFld, loc);
    double *j_i = gkyl_array_fetch(jFld, loc);
    double *jinv_i = gkyl_array_fetch(jinvFld, loc);
    double *grij = gkyl_array_fetch(grFld, loc);
    double *bi_i = gkyl_array_fetch(biFld, loc);
    double *cmag_i = gkyl_array_fetch(cmagFld, loc);
    double *jtot_i = gkyl_array_fetch(jtotFld, loc);
    double *jtotinv_i = gkyl_array_fetch(jtotinvFld, loc);
    double *gxxJ_i= gkyl_array_fetch(gxxJFld, loc);
    double *gxyJ_i= gkyl_array_fetch(gxyJFld, loc);
    double *gyyJ_i= gkyl_array_fetch(gyyJFld, loc);
    double *gxzJ_i= gkyl_array_fetch(gxzJFld, loc);
    double *eps2_i= gkyl_array_fetch(eps2Fld, loc);
    up->kernel(gij, bmag_i, j_i, jinv_i, grij, bi_i, cmag_i, jtot_i, jtotinv_i, gxxJ_i, gxyJ_i, gyyJ_i, gxzJ_i, eps2_i);
  }
}

void
gkyl_calc_derived_geo_release(gkyl_calc_derived_geo* up)
{
  gkyl_free(up);
}
