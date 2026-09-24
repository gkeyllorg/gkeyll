#include <assert.h>
#include <string.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>

#include <gkyl_cart_modal_hybrid_priv.h>

void
gkyl_cart_modal_hybrid(struct gkyl_basis *basis, int cdim, int vdim)
{
  int ndim = cdim + vdim;
  assert(ndim > 1 && ndim < 7);
  assert(cdim < 4 && vdim > 0 && vdim < 4);

  basis->ndim = ndim;
  basis->poly_order = 1;
  basis->num_basis = num_basis_list[cdim].count[vdim];
  basis->num_quad = num_quad_list[ndim].count[1];
  strcpy(basis->id, "hybrid");
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

struct gkyl_basis *
gkyl_cart_modal_hybrid_new(int cdim, int vdim)
{
  struct gkyl_basis *basis = gkyl_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_hybrid(basis, cdim, vdim);
  return basis;
}

#ifndef GKYL_HAVE_CUDA
void
gkyl_cart_modal_hybrid_cu_dev(struct gkyl_basis *basis, int cdim, int vdim)
{
  assert(false);
}

struct gkyl_basis *
gkyl_cart_modal_hybrid_cu_dev_new(int cdim, int vdim)
{
  assert(false);
}
#endif
