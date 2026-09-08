#include <assert.h>
#include <string.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>

#include <gkyl_cart_modal_gkhybrid_vel_priv.h>

void gkyl_cart_modal_gkhybrid_vel(struct gkyl_basis *basis, int vdim)
{
  int ndim = vdim;
  assert(ndim == 1 || ndim == 2);

  basis->ndim = ndim;
  basis->poly_order = 1;
  basis->num_basis = num_basis_list[ndim].count[1];
  basis->num_quad = num_quad_list[ndim].count[1];
  strcpy(basis->id, "gkhybrid_vel");
  basis->b_type = GKYL_BASIS_MODAL_GKHYBRID_VEL;

  // function pointers
  basis->eval = ev_list[ndim].ev[1];
  basis->eval_expand = eve_list[ndim].ev[1];
  basis->eval_grad_expand = eveg_list[ndim].ev[1];
  basis->flip_odd_sign = fos_list[ndim].fs[1];
  basis->flip_even_sign = fes_list[ndim].fs[1];
  basis->node_list = nl_list[ndim].nl[1];
  basis->nodal_to_modal = n2m_list[ndim].n2m[1];
  basis->quad_nodal_to_modal = qn2m_list[ndim].n2m[1];
  basis->modal_to_quad_nodal = m2qn_list[ndim].n2m[1];
}

struct gkyl_basis *gkyl_cart_modal_gkhybrid_vel_new(int vdim)
{
  struct gkyl_basis *basis = gkyl_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_gkhybrid_vel(basis, vdim);
  return basis;
}

#ifndef GKYL_HAVE_CUDA
void gkyl_cart_modal_gkhybrid_vel_cu_dev(struct gkyl_basis *basis, int vdim)
{
  assert(false);
}

struct gkyl_basis *gkyl_cart_modal_gkhybrid_vel_cu_dev_new(int vdim)
{
  assert(false);
}
#endif
