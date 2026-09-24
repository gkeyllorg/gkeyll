/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
}

__global__ static void
gkyl_prim_lbo_vlasov_set_cu_dev_ptrs(struct prim_lbo_type_vlasov *prim_vlasov, int cdim, int vdim, int poly_order, enum gkyl_basis_type b_type, int tblidx)
{
  prim_vlasov->prim.self_prim = self_prim;
  prim_vlasov->prim.cross_prim = cross_prim;
  
  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_self_kern_list *self_prim_kernels;
  const gkyl_prim_lbo_vlasov_cross_kern_list *cross_prim_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      cross_prim_kernels = ser_cross_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  prim_vlasov->self_prim = self_prim_kernels[tblidx].kernels[poly_order];
  prim_vlasov->cross_prim = cross_prim_kernels[tblidx].kernels[poly_order];
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_vlasov *prim_vlasov =
    (struct prim_lbo_type_vlasov*) gkyl_malloc(sizeof(struct prim_lbo_type_vlasov));
  
  int cdim = prim_vlasov->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov->prim.pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = prim_vlasov->prim.poly_order = cbasis->poly_order;
  prim_vlasov->prim.num_config = cbasis->num_basis;
  prim_vlasov->prim.num_phase = pbasis->num_basis;
  prim_vlasov->prim.udim = vdim;

  prim_vlasov->prim.flag = 0;
  GKYL_SET_CU_ALLOC(prim_vlasov->prim.flag);
  prim_vlasov->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_free);
  
  // copy the host struct to device struct
  struct prim_lbo_type_vlasov *prim_vlasov_cu = (struct prim_lbo_type_vlasov*)
    gkyl_cu_malloc(sizeof(struct prim_lbo_type_vlasov));
  gkyl_cu_memcpy(prim_vlasov_cu, prim_vlasov, sizeof(struct prim_lbo_type_vlasov), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);
  
  gkyl_prim_lbo_vlasov_set_cu_dev_ptrs<<<1,1>>>(prim_vlasov_cu, cdim, vdim, poly_order,
    cbasis->b_type, cv_index[cdim].vdim[vdim]);

  prim_vlasov->prim.on_dev = &prim_vlasov_cu->prim;
    
  return &prim_vlasov->prim;
}
