#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>

gkyl_dg_cx*
gkyl_dg_cx_new(struct gkyl_dg_cx_inp *inp, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_dg_cx_cu_dev_new(inp);
  } 
#endif  
  gkyl_dg_cx *up = gkyl_malloc(sizeof(struct gkyl_dg_cx));

  up->cbasis = inp->cbasis;
  up->conf_rng = inp->conf_rng;
  up->type_ion = inp->type_ion;
  up->vt_sq_ion_min = inp->vt_sq_ion_min;
  up->vt_sq_neut_min = inp->vt_sq_neut_min;

  // Get parameters for fitting function.
  fit_param(up->type_ion, &up->a, &up->b);

  up->react_rate = choose_kern(*up->cbasis);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up;
  
  return up;
}

void gkyl_dg_cx_coll(const struct gkyl_dg_cx *up, 
  struct gkyl_array *maxwellian_moms_ion, struct gkyl_array *maxwellian_moms_neut,
  struct gkyl_array *upar_b_i, struct gkyl_array *coef_cx, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(coef_cx)) {
    return gkyl_dg_cx_coll_cu(up, maxwellian_moms_ion, maxwellian_moms_neut, 
      upar_b_i, coef_cx, cflrate);
  }
#endif
  
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *maxwellian_moms_ion_d = gkyl_array_cfetch(maxwellian_moms_ion, linidx);
    const double *maxwellian_moms_neut_d = gkyl_array_cfetch(maxwellian_moms_neut, linidx);
    const double *upar_b_i_d = gkyl_array_cfetch(upar_b_i, linidx);

    double *coef_cx_d = gkyl_array_fetch(coef_cx, linidx);
    
    double cflr = up->react_rate(up->a, up->b, up->vt_sq_ion_min, up->vt_sq_neut_min,
      maxwellian_moms_ion_d, maxwellian_moms_neut_d, upar_b_i_d, coef_cx_d);
  }
}

void
gkyl_dg_cx_release(gkyl_dg_cx* cx)
{
  free(cx);
}
