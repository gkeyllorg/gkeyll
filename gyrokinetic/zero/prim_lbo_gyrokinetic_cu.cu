/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_gyrokinetic_priv.h>
}

__global__ static void
gkyl_prim_lbo_gyrokinetic_set_cu_dev_ptrs(struct prim_lbo_type_gyrokinetic *prim_gyrokinetic, int cdim, int vdim, int poly_order, enum gkyl_basis_type b_type)
{
  int pdim = cdim+vdim;

  prim_gyrokinetic->prim.self_prim = self_prim;
  prim_gyrokinetic->prim.cross_prim = cross_prim;
  
  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_gyrokinetic_kern_list *self_prim_kernels;
  const gkyl_prim_lbo_gyrokinetic_cross_kern_list *cross_prim_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      cross_prim_kernels = ser_cross_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  prim_gyrokinetic->self_prim = self_prim_kernels[pdim-2].kernels[poly_order];
  prim_gyrokinetic->cross_prim = cross_prim_kernels[pdim-2].kernels[poly_order];
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic =
    (struct prim_lbo_type_gyrokinetic*) gkyl_malloc(sizeof(struct prim_lbo_type_gyrokinetic));
  
  int cdim = prim_gyrokinetic->prim.cdim = cbasis->ndim;
  int pdim = prim_gyrokinetic->prim.pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = prim_gyrokinetic->prim.poly_order = cbasis->poly_order;
  prim_gyrokinetic->prim.num_config = cbasis->num_basis;
  prim_gyrokinetic->prim.num_phase = pbasis->num_basis;
  prim_gyrokinetic->prim.udim = 1;

  prim_gyrokinetic->prim.flag = 0;
  GKYL_SET_CU_ALLOC(prim_gyrokinetic->prim.flag);
  prim_gyrokinetic->prim.ref_count = gkyl_ref_count_init(prim_lbo_gyrokinetic_free);
  
  // copy the host struct to device struct
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic_cu = (struct prim_lbo_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct prim_lbo_type_gyrokinetic));
  gkyl_cu_memcpy(prim_gyrokinetic_cu, prim_gyrokinetic, sizeof(struct prim_lbo_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  gkyl_prim_lbo_gyrokinetic_set_cu_dev_ptrs<<<1,1>>>(prim_gyrokinetic_cu, cdim, vdim, poly_order,
    cbasis->b_type);

  prim_gyrokinetic->prim.on_dev = &prim_gyrokinetic_cu->prim;
    
  return &prim_gyrokinetic->prim;
}
