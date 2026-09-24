#include <gkyl_proj_powsqrt_on_basis.h>

struct gkyl_proj_powsqrt_on_basis {
  int num_quad; // number of quadrature points to use in each direction
  int ndim; // Configuration-space dimension

  int num_basis; // number of conf-space basis functions

  bool use_gpu;

  // for quadrature in conf space
  int tot_quad; // total number of quadrature points
  struct gkyl_array *weights; // weights for conf-space quadrature
  struct gkyl_array *basis_at_ords; // conf-space basis functions at ordinates

  struct gkyl_array *fun_at_ords; // function (Maxwellian) evaluated at
                                  // ordinates in a cell.
};

void
gkyl_proj_powsqrt_on_basis_advance_cu(const gkyl_proj_powsqrt_on_basis *up,
  const struct gkyl_range *range, double expIn, const struct gkyl_array *fIn,
  struct gkyl_array *fOut);
