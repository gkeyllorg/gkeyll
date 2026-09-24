#include <math.h>

#include <gkyl_fv_proj.h>
#include <gkyl_array_ops.h>

gkyl_fv_proj*
gkyl_fv_proj_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  // This updater is just a wrapper around more general
  // gkyl_proj_on_basis updater, however specialized to poly_order=0
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, grid->ndim, 0);
  return gkyl_proj_on_basis_new(grid, &basis, num_quad, num_ret_vals, eval, ctx);
}

void
gkyl_fv_proj_advance(const gkyl_fv_proj *pob,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  gkyl_proj_on_basis_advance(pob, tm, update_rng, out);

  // from projections, compute cell average
  double denorm = 1.0/sqrt(pow(2, update_rng->ndim));
  gkyl_array_scale_range(out, denorm, update_rng);
}

void
gkyl_fv_proj_release(gkyl_fv_proj* pob)
{
  gkyl_proj_on_basis_release(pob);
}
