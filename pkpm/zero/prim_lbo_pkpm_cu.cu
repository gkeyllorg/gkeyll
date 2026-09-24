/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_pkpm.h>
#include <gkyl_prim_lbo_pkpm_priv.h>
}

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__ static void
gkyl_prim_lbo_pkpm_set_cu_dev_ptrs(struct prim_lbo_type_pkpm *prim_pkpm, 
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  prim_pkpm->prim.self_prim = self_prim;
  
  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_pkpm_self_kern_list *self_prim_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      self_prim_kernels = ten_self_prim_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }

  prim_pkpm->self_prim = CK(self_prim_kernels, cdim, poly_order);
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_pkpm_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_pkpm *prim_pkpm =
    (struct prim_lbo_type_pkpm*) gkyl_malloc(sizeof(struct prim_lbo_type_pkpm));
  
  int cdim = prim_pkpm->prim.cdim = cbasis->ndim;
  int pdim = prim_pkpm->prim.pdim = pbasis->ndim;

  int poly_order = prim_pkpm->prim.poly_order = cbasis->poly_order;
  prim_pkpm->prim.num_config = cbasis->num_basis;
  prim_pkpm->prim.num_phase = pbasis->num_basis;
  prim_pkpm->prim.udim = 1;

  prim_pkpm->conf_range = *conf_range;

  prim_pkpm->prim.flag = 0;
  GKYL_SET_CU_ALLOC(prim_pkpm->prim.flag);
  prim_pkpm->prim.ref_count = gkyl_ref_count_init(prim_lbo_pkpm_free);
  
  // copy the host struct to device struct
  struct prim_lbo_type_pkpm *prim_pkpm_cu = (struct prim_lbo_type_pkpm*)
    gkyl_cu_malloc(sizeof(struct prim_lbo_type_pkpm));
  gkyl_cu_memcpy(prim_pkpm_cu, prim_pkpm, sizeof(struct prim_lbo_type_pkpm), GKYL_CU_MEMCPY_H2D);
  
  gkyl_prim_lbo_pkpm_set_cu_dev_ptrs<<<1,1>>>(prim_pkpm_cu, cbasis->b_type, cdim, poly_order);

  prim_pkpm->prim.on_dev = &prim_pkpm_cu->prim;
    
  return &prim_pkpm->prim;
}
