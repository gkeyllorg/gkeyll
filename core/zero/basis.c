#include <gkyl_alloc.h>
#include <gkyl_basis.h>

void
gkyl_cart_modal_basis_release(struct gkyl_basis *basis)
{
  gkyl_free(basis);
}

void
gkyl_cart_modal_basis_release_cu(struct gkyl_basis *basis)
{
  gkyl_cu_free(basis);
}

unsigned
gkyl_cart_modal_basis_get_ndim(const struct gkyl_basis *basis)
{
  return basis->ndim;
}

unsigned
gkyl_cart_modal_basis_get_poly_order(const struct gkyl_basis *basis)
{
  return basis->poly_order;
}

unsigned
gkyl_cart_modal_basis_get_num_basis(const struct gkyl_basis *basis)
{
  return basis->num_basis;
}

enum gkyl_basis_type
gkyl_cart_modal_basis_get_type(const struct gkyl_basis *basis)
{
  return basis->b_type;
}
