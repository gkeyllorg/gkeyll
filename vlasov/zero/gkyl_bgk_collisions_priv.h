#include <gkyl_bgk_collisions.h>
#include <gkyl_dg_bin_ops_priv.h>

struct gkyl_bgk_collisions {
  unsigned cdim; // Configuration-space dimension.
  unsigned vdim; // Velocity-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned pnum_basis; // Number of phase-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  enum gkyl_basis_type pb_type; // Phase basis type.

  double cellav_fac; // Multiply 0-th DG coeff by this to get the cell avg.

  mul_op_t mul_op; // Conf*phase bin_op multiplication kernel.

  bool use_gpu;
};

void
gkyl_bgk_collisions_advance_cu(const gkyl_bgk_collisions *up,
  const struct gkyl_range *crange, const struct gkyl_range *prange,
  const struct gkyl_array *nu, const struct gkyl_array *nufM, const struct gkyl_array *fin,
  bool implicit_step, double dt, struct gkyl_array *out, struct gkyl_array *cflfreq);
