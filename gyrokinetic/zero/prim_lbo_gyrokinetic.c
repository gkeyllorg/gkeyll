#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_gyrokinetic_priv.h>

void
prim_lbo_gyrokinetic_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_prim_lbo_type *prim_ty = container_of(ref, struct gkyl_prim_lbo_type, ref_count);
  if (GKYL_IS_CU_ALLOC(prim_ty->flag))
    gkyl_cu_free(prim_ty->on_dev);

  struct prim_lbo_type_gyrokinetic *gk = container_of(prim_ty, struct prim_lbo_type_gyrokinetic, prim);
  gkyl_free(gk);
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_gyrokinetic_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_prim_lbo_gyrokinetic_cu_dev_new(cbasis, pbasis);
  } 
#endif     
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic = gkyl_malloc(sizeof(struct prim_lbo_type_gyrokinetic));
  int cdim = prim_gyrokinetic->prim.cdim = cbasis->ndim;
  int pdim = prim_gyrokinetic->prim.pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = prim_gyrokinetic->prim.poly_order = cbasis->poly_order;
  prim_gyrokinetic->prim.num_config = cbasis->num_basis;
  prim_gyrokinetic->prim.num_phase = pbasis->num_basis;
  prim_gyrokinetic->prim.self_prim = self_prim;
  prim_gyrokinetic->prim.cross_prim = cross_prim;
  prim_gyrokinetic->prim.udim = 1;

  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_gyrokinetic_kern_list *self_prim_kernels;
  const gkyl_prim_lbo_gyrokinetic_cross_kern_list *cross_prim_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      cross_prim_kernels = ser_cross_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  assert(NULL != self_prim_kernels[pdim-2].kernels[poly_order]);
  assert(NULL != cross_prim_kernels[pdim-2].kernels[poly_order]);
    
  prim_gyrokinetic->self_prim = self_prim_kernels[pdim-2].kernels[poly_order];
  prim_gyrokinetic->cross_prim = cross_prim_kernels[pdim-2].kernels[poly_order];

  prim_gyrokinetic->prim.flag = 0;
  GKYL_CLEAR_CU_ALLOC(prim_gyrokinetic->prim.flag);
  prim_gyrokinetic->prim.ref_count = gkyl_ref_count_init(prim_lbo_gyrokinetic_free);

  prim_gyrokinetic->prim.on_dev = &prim_gyrokinetic->prim;
    
  return &prim_gyrokinetic->prim;
}
