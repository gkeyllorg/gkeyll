#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_vlasov_cross_prim_moms_bgk.h>
#include <gkyl_vlasov_cross_prim_moms_bgk_priv.h>

gkyl_vlasov_cross_prim_moms_bgk* 
gkyl_vlasov_cross_prim_moms_bgk_new(const struct gkyl_basis *pbasis,
  const struct gkyl_basis *cbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_vlasov_cross_prim_moms_bgk_cu_dev_new(pbasis, cbasis); 
  }
#endif
  gkyl_vlasov_cross_prim_moms_bgk *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = use_gpu;

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  // Select the kernel.
  // choose kernel tables based on basis-function type
  const gkyl_vlasov_cross_prim_moms_bgk_kern_list *cross_prim_moms_bgk_kernels;

  switch (gkyl_basis_phase_kernel_type(cbasis, pbasis)) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      cross_prim_moms_bgk_kernels = ser_vlasov_cross_prim_moms_bgk_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      cross_prim_moms_bgk_kernels = ten_vlasov_cross_prim_moms_bgk_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != cross_prim_moms_bgk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
  up->cross_prim_moms_calc = cross_prim_moms_bgk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  
  up->on_dev = up; // host-side points to itself
  return up;   
}

void gkyl_vlasov_cross_prim_moms_bgk_advance(gkyl_vlasov_cross_prim_moms_bgk *up,
  const struct gkyl_range *conf_rng, double delta_sr, double betap1,
  double m_self, const struct gkyl_array *prim_moms_self, double m_other, const struct gkyl_array *prim_moms_other,
  struct gkyl_array *prim_moms_cross)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)  {
    return gkyl_vlasov_cross_prim_moms_bgk_advance_cu(up, conf_rng, delta_sr, betap1, 
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

void gkyl_vlasov_cross_prim_moms_bgk_release(gkyl_vlasov_cross_prim_moms_bgk *up)
{
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
