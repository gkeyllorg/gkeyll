#include <gkyl_calc_bmag.h>
#include <gkyl_calc_bmag_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_on_basis.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_bmag*
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, bool use_gpu)
{
  gkyl_calc_bmag *up = gkyl_malloc(sizeof(gkyl_calc_bmag));
  up->cbasis = cbasis;
  up->pbasis = pbasis;
  up->cgrid = cgrid;
  up->pgrid = pgrid;
  up->use_gpu = use_gpu;
  return up;
}

void gkyl_calc_bmag_global(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_bmag_ctx *gc = (struct gkyl_bmag_ctx*) ctx;
  // Need a crude an manual deflated coordinate because this works on deflated geometry due to the allgather
  double xpt[GKYL_MAX_CDIM];
  if (gc->cgrid->ndim == 1)
    xpt[0] = xn[2];
  else if (gc->cgrid->ndim == 2){
    xpt[0] = xn[0];
    xpt[1] = xn[2];
  }
  else{
    xpt[0] = xn[0];
    xpt[1] = xn[1];
    xpt[2] = xn[2];
  }
  int cidx[GKYL_MAX_CDIM];
  for(int i = 0; i < gc->cgrid->ndim; i++){
    int idxtemp = gc->crange_global->lower[i] + (int) floor((xpt[i] - (gc->cgrid->lower[i]) )/gc->cgrid->dx[i]);
    idxtemp = GKYL_MIN2(idxtemp, gc->crange_global->upper[i]);
    idxtemp = GKYL_MAX2(idxtemp, gc->crange_global->lower[i]);
    cidx[i] = idxtemp;
  }
  long lidx = gkyl_range_idx(gc->crange_global, cidx);
  const double *mcoeffs = gkyl_array_cfetch(gc->bmag, lidx);
  double cxc[gc->cgrid->ndim];
  double xyz[gc->cgrid->ndim];
  gkyl_rect_grid_cell_center(gc->cgrid, cidx, cxc);
  for(int i = 0; i < gc->cgrid->ndim; i++)
    xyz[i] = (xpt[i]-cxc[i])/(gc->cgrid->dx[i]*0.5);
  fout[0] = gc->cbasis->eval_expand(xyz, mcoeffs);
}

static inline void bmag_comp(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_bmag_ctx *gc = (struct gkyl_bmag_ctx*) ctx;
  double RZPHI[gc->cgrid->ndim];

  int cidx[GKYL_MAX_CDIM];
  for(int i = 0; i < gc->cgrid->ndim; i++){
    int idxtemp = gc->crange_global->lower[i] + (int) floor((xn[i] - (gc->cgrid->lower[i]) )/gc->cgrid->dx[i]);
    idxtemp = GKYL_MIN2(idxtemp, gc->crange->upper[i]);
    idxtemp = GKYL_MAX2(idxtemp, gc->crange->lower[i]);
    cidx[i] = idxtemp;
  }

  long lidx = gkyl_range_idx(gc->crange, cidx);
  const double *mcoeffs = gkyl_array_cfetch(gc->mapc2p, lidx);
  
  double cxc[gc->cgrid->ndim];
  double xyz[gc->cgrid->ndim];
  gkyl_rect_grid_cell_center(gc->cgrid, cidx, cxc);
  for(int i = 0; i < gc->cgrid->ndim; i++)
    xyz[i] = (xn[i]-cxc[i])/(gc->cgrid->dx[i]*0.5);
  for(int i = 0; i < gc->cgrid->ndim; i++){
    RZPHI[i] = gc->cbasis->eval_expand(xyz, &mcoeffs[i*gc->cbasis->num_basis]);
  }

  double R = RZPHI[0];
  double Z = RZPHI[1];

  int rzidx[2];
  int idxtemp = gc->range->lower[0] + (int) floor((R - gc->grid->lower[0])/gc->grid->dx[0]);
  idxtemp = GKYL_MIN2(idxtemp, gc->range->upper[0]);
  idxtemp = GKYL_MAX2(idxtemp, gc->range->lower[0]);
  rzidx[0] = idxtemp;
  idxtemp = gc->range->lower[1] + (int) floor((Z - gc->grid->lower[1])/gc->grid->dx[1]);
  idxtemp = GKYL_MIN2(idxtemp, gc->range->upper[1]);
  idxtemp = GKYL_MAX2(idxtemp, gc->range->lower[1]);
  rzidx[1] = idxtemp;


  long loc = gkyl_range_idx(gc->range, rzidx);
  const double *coeffs = gkyl_array_cfetch(gc->bmagdg,loc);

  double xc[2];
  gkyl_rect_grid_cell_center(gc->grid, rzidx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(gc->grid->dx[0]*0.5);
  xy[1] = (Z-xc[1])/(gc->grid->dx[1]*0.5);
  fout[0] = gc->basis->eval_expand(xy, coeffs);
}

void gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, 
  const struct gkyl_range *crange, const struct gkyl_range *crange_ext, const struct gkyl_range *crange_global, 
  const struct gkyl_range *prange, const struct gkyl_range *prange_ext, 
  const struct gkyl_array *bmagrz, struct gkyl_array* bmag_compdg, struct gkyl_array* mapc2p, bool use_quad)
{
  // Convert bmag into computational coordinates
  struct gkyl_bmag_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->grid = up->pgrid;
  ctx->cgrid = up->cgrid;
  ctx->range = prange;
  ctx->crange = crange;
  ctx->crange_global = crange_global;
  ctx->bmagdg = bmagrz;
  ctx->basis = up->pbasis;
  ctx->cbasis = up->cbasis;
  ctx->mapc2p = mapc2p;
  if (use_quad) {
    gkyl_proj_on_basis *eval_bmag_comp = gkyl_proj_on_basis_new(up->cgrid, up->cbasis, 2, 1, bmag_comp, ctx);
    gkyl_proj_on_basis_advance(eval_bmag_comp, 0.0, crange, bmag_compdg); //on ghosts with ext_range
    gkyl_proj_on_basis_release(eval_bmag_comp);
  }
  else {
    gkyl_eval_on_nodes *eval_bmag_comp = gkyl_eval_on_nodes_new(up->cgrid, up->cbasis, 1, bmag_comp, ctx);
    gkyl_eval_on_nodes_advance(eval_bmag_comp, 0.0, crange, bmag_compdg); //on ghosts with ext_range
    gkyl_eval_on_nodes_release(eval_bmag_comp);
  }

  // Free allocated memory
  gkyl_free(ctx);
}

void
gkyl_calc_bmag_release(gkyl_calc_bmag* up)
{
  gkyl_free(up);
}
