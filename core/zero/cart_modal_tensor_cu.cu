/* -*- c++ -*- */

extern "C" {
#include <assert.h>
#include <string.h>    
    
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_cart_modal_tensor_priv.h>
#include <gkyl_util.h>
}

__global__ void static
gkyl_cart_modal_tensor_cu_dev_kern(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ev_list[ndim].ev[poly_order]);

  basis->ndim = ndim;
  basis->poly_order = poly_order;
  basis->num_basis = pow(poly_order+1, ndim);
  basis->num_quad = num_quad_list[ndim].count[poly_order];
  basis->b_type = GKYL_BASIS_MODAL_TENSOR;  
  
  // function pointers
  basis->eval = ev_list[ndim].ev[poly_order];
  basis->eval_expand = eve_list[ndim].ev[poly_order];
  basis->eval_grad_expand = eveg_list[ndim].ev[poly_order];
  basis->flip_odd_sign = fos_list[ndim].fs[poly_order];
  basis->flip_even_sign = fes_list[ndim].fs[poly_order];
  basis->node_list = nl_list[ndim].nl[poly_order];
  basis->nodal_to_modal = n2m_list[ndim].n2m[poly_order];
  basis->quad_nodal_to_modal = qn2m_list[ndim].n2m[poly_order];
  basis->modal_to_quad_nodal = m2qn_list[ndim].n2m[poly_order];
}

void
gkyl_cart_modal_tensor_cu_dev(struct gkyl_basis *basis, int ndim, int poly_order)
{
  assert(ndim>0 && ndim<=6);

  struct gkyl_basis ho_basis;

  strcpy(ho_basis.id, "tensor");
  // this copy needs to be done here as the strcpy needed in the
  // "type" field can't be done on the device
  gkyl_cu_memcpy(basis, &ho_basis, sizeof(struct gkyl_basis),
    GKYL_CU_MEMCPY_H2D);
  
  gkyl_cart_modal_tensor_cu_dev_kern<<<1,1>>>(basis, ndim, poly_order);
}

struct gkyl_basis *
gkyl_cart_modal_tensor_cu_dev_new(int ndim, int poly_order)
{
  struct gkyl_basis *basis = (struct gkyl_basis *) gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_tensor_cu_dev(basis, ndim, poly_order);
  return basis;
}
