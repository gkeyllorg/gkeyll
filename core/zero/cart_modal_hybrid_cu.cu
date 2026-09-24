/* -*- c++ -*- */

#include <gkyl_util.h>
extern "C" {
#include <assert.h>
#include <string.h>    
    
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>

#include <gkyl_cart_modal_hybrid_priv.h>
}

__global__ void static
gkyl_cart_modal_hybrid_cu_dev_kern(struct gkyl_basis *basis, int cdim, int vdim)
{
  assert(ev_list[cdim].ev[vdim]);

  basis->ndim = cdim+vdim;
  basis->poly_order = 1;
  basis->num_basis = num_basis_list[cdim].count[vdim];
  basis->num_quad = num_quad_list[cdim].count[vdim];
  basis->b_type = GKYL_BASIS_MODAL_HYBRID;  
  
  // function pointers
  basis->eval = ev_list[cdim].ev[vdim];
  basis->eval_expand = eve_list[cdim].ev[vdim];
  basis->eval_grad_expand = eveg_list[cdim].ev[vdim];
  basis->flip_odd_sign = fos_list[cdim].fs[vdim];
  basis->flip_even_sign = fes_list[cdim].fs[vdim];
  basis->node_list = nl_list[cdim].nl[vdim];
  basis->nodal_to_modal = n2m_list[cdim].n2m[vdim];
  basis->quad_nodal_to_modal = qn2m_list[cdim].n2m[vdim];
  basis->modal_to_quad_nodal = m2qn_list[cdim].n2m[vdim];
}

void
gkyl_cart_modal_hybrid_cu_dev(struct gkyl_basis *basis, int cdim, int vdim)
{
  int ndim = cdim+vdim;
  assert(ndim>1 && ndim<7);
  assert(cdim<4 && vdim>0 && vdim<4);

  struct gkyl_basis ho_basis;

  strcpy(ho_basis.id, "hybrid");
  // this copy needs to be done here as the strcpy needed in the
  // "type" field can't be done on the device
  gkyl_cu_memcpy(basis, &ho_basis, sizeof(struct gkyl_basis),
    GKYL_CU_MEMCPY_H2D);
  
  gkyl_cart_modal_hybrid_cu_dev_kern<<<1,1>>>(basis, cdim, vdim);
}

struct gkyl_basis *
gkyl_cart_modal_hybrid_cu_dev_new(int cdim, int vdim)
{
  struct gkyl_basis *basis = (struct gkyl_basis *) gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_hybrid_cu_dev(basis, cdim, vdim);
  return basis;
}
