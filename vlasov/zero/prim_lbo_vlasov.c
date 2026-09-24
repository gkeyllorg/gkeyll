#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
#include <gkyl_util.h>

void
prim_lbo_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_prim_lbo_type *prim_ty = container_of(ref, struct gkyl_prim_lbo_type, ref_count);
  if (GKYL_IS_CU_ALLOC(prim_ty->flag))
    gkyl_cu_free(prim_ty->on_dev);

  struct prim_lbo_type_vlasov *vlasov = container_of(prim_ty, struct prim_lbo_type_vlasov, prim);
  gkyl_free(vlasov);
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_prim_lbo_vlasov_cu_dev_new(cbasis, pbasis);
  } 
#endif  
  struct prim_lbo_type_vlasov *prim_vlasov = gkyl_malloc(sizeof(struct prim_lbo_type_vlasov));

  int cdim = prim_vlasov->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov->prim.pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = prim_vlasov->prim.poly_order = cbasis->poly_order;
  prim_vlasov->prim.num_config = cbasis->num_basis;
  prim_vlasov->prim.num_phase = pbasis->num_basis;
  prim_vlasov->prim.self_prim = self_prim;
  prim_vlasov->prim.cross_prim = cross_prim;
  prim_vlasov->prim.udim = vdim;

  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_self_kern_list *self_prim_kernels;
  const gkyl_prim_lbo_vlasov_cross_kern_list *cross_prim_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      cross_prim_kernels = ser_cross_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
  assert(NULL != cross_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
  prim_vlasov->self_prim = self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  prim_vlasov->cross_prim = cross_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

  prim_vlasov->prim.flag = 0;
  GKYL_CLEAR_CU_ALLOC(prim_vlasov->prim.flag);
  prim_vlasov->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_free);

  prim_vlasov->prim.on_dev = &prim_vlasov->prim;
    
  return &prim_vlasov->prim;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(false);
  return 0;
}

#endif
