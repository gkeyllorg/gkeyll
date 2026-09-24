#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk_priv.h>

gkyl_gyrokinetic_cross_prim_moms_bgk* 
gkyl_gyrokinetic_cross_prim_moms_bgk_new(const struct gkyl_basis *phase_basis,
  const struct gkyl_basis *conf_basis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_gyrokinetic_cross_prim_moms_bgk_cu_dev_new(phase_basis, conf_basis); 
  }
#endif
  gkyl_gyrokinetic_cross_prim_moms_bgk *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;

  // Select the kernel.
  up->cross_prim_moms_calc = choose_gyrokinetic_cross_prim_moms_bgk_kern(conf_basis->ndim,
    phase_basis->ndim-conf_basis->ndim, phase_basis->poly_order); 
  
  up->on_dev = up; // host-side points to itself
  return up;   
}

void gkyl_gyrokinetic_cross_prim_moms_bgk_advance(gkyl_gyrokinetic_cross_prim_moms_bgk *up,
  const struct gkyl_range *conf_rng, double delta_sr, double betap1,
  double m_self, const struct gkyl_array *prim_moms_self, double m_other, const struct gkyl_array *prim_moms_other,
  struct gkyl_array *prim_moms_cross)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)  {
    return gkyl_gyrokinetic_cross_prim_moms_bgk_advance_cu(up, conf_rng, delta_sr, betap1, 
      m_self, prim_moms_self, m_other, prim_moms_other, prim_moms_cross);
  }
#endif
  struct gkyl_range_iter conf_iter;  
  
  // loop over configuration space cells
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);
    const double *prim_moms_self_d = gkyl_array_cfetch(prim_moms_self, midx);
    const double *prim_moms_other_d = gkyl_array_cfetch(prim_moms_other, midx);
    double *out_d = gkyl_array_fetch(prim_moms_cross, midx);

    up->cross_prim_moms_calc(delta_sr, betap1, m_self, prim_moms_self_d,
      m_other, prim_moms_other_d, out_d);
  }
}

void gkyl_gyrokinetic_cross_prim_moms_bgk_release(gkyl_gyrokinetic_cross_prim_moms_bgk *up)
{
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
