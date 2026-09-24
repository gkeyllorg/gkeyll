#include <gkyl_bgk_collisions.h>
#include <gkyl_bgk_collisions_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_array_ops_priv.h>

gkyl_bgk_collisions*
gkyl_bgk_collisions_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  bool use_gpu)
{
  gkyl_bgk_collisions *up = gkyl_malloc(sizeof(gkyl_bgk_collisions));

  up->cdim = cbasis->ndim;
  up->vdim = pbasis->ndim-up->cdim;
  up->cnum_basis = cbasis->num_basis;
  up->pnum_basis = pbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->use_gpu = use_gpu;
  assert(up->pnum_basis > up->cnum_basis);
  assert(up->pnum_basis <= 160); // MF 2022/11/18: hardcode to 3x3v p=1 hybrid.

  int poly_order = cbasis->poly_order;
  up->pb_type = pbasis->b_type;
  if (!up->use_gpu)
    up->mul_op = choose_mul_conf_phase_kern(up->pb_type, up->cdim, up->vdim, poly_order);

  up->cellav_fac = 1./sqrt(pow(2,up->cdim));

  return up;
}

void
gkyl_bgk_collisions_advance(const gkyl_bgk_collisions *up,
  const struct gkyl_range *crange, const struct gkyl_range *prange,
  const struct gkyl_array *nu, const struct gkyl_array *nufM, const struct gkyl_array *fin,
  bool implicit_step, double dt, struct gkyl_array *out, struct gkyl_array *cflfreq)
{
  // Compute nu*f_M - nu*f, and its contribution to the CFL rate.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_bgk_collisions_advance_cu(up, crange, prange, nu, nufM, fin, implicit_step, dt, out, cflfreq);
#endif

  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, prange);

  if (implicit_step) {

    const double cellav_fac_dt = up->cellav_fac*dt;

    while (gkyl_range_iter_next(&piter)) {
      long ploc = gkyl_range_idx(prange, piter.idx);
      long cloc = gkyl_range_idx(crange, piter.idx);

      const double *nu_d = gkyl_array_cfetch(nu, cloc);
      double *out_d = gkyl_array_fetch(out, ploc);

      // Add nu*f_M.
      array_acc1(up->pnum_basis, out_d, 1./(1.0 + nu_d[0]*cellav_fac_dt), gkyl_array_cfetch(nufM, ploc));

      // Calculate and add -nu*f.
      double incr[160]; // mul_op assigns, but need increment, so use a buffer.
      up->mul_op(nu_d, gkyl_array_cfetch(fin, ploc), incr);
      array_acc1(up->pnum_basis, out_d, -1.0/(1.0 + nu_d[0]*cellav_fac_dt), incr);

      // No CFL contribution in the implicit case.
    }

  }
  else {

    while (gkyl_range_iter_next(&piter)) {
      long ploc = gkyl_range_idx(prange, piter.idx);
      long cloc = gkyl_range_idx(crange, piter.idx);
  
      const double *nu_d = gkyl_array_cfetch(nu, cloc);
      double *out_d = gkyl_array_fetch(out, ploc);
  
      // Add nu*f_M.
      array_acc1(up->pnum_basis, out_d, 1., gkyl_array_cfetch(nufM, ploc));
  
      // Calculate and add -nu*f.
      double incr[160]; // mul_op assigns, but need increment, so use a buffer.
      up->mul_op(nu_d, gkyl_array_cfetch(fin, ploc), incr);
      array_acc1(up->pnum_basis, out_d, -1., incr);
  
      // Add contribution to CFL frequency.
      double *cflfreq_d = gkyl_array_fetch(cflfreq, ploc);
      cflfreq_d[0] += nu_d[0]*up->cellav_fac;
    }

  }
}

void
gkyl_bgk_collisions_release(gkyl_bgk_collisions* up)
{
  gkyl_free(up);
}
