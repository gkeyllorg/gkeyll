#include <gkyl_spitzer_coll_freq.h>

struct gkyl_spitzer_coll_freq {
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

  // Time independent factors that can be precomputed.
  double hbar_fac, r4pieps0_fac, nufraceps0_fac, cellav_fac;
  double eps0;
};

void
gkyl_spitzer_coll_freq_advance_normnu_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *momsSelf, double vtSqMinSelf,
  const struct gkyl_array *momsOther, double vtSqMinOther,
  double normNu, struct gkyl_array *nuOut);

void
gkyl_spitzer_coll_freq_advance_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *bmag,
  double qSelf, double mSelf, const struct gkyl_array *momsSelf, double vtSqMinSelf,
  double qOther, double mOther, const struct gkyl_array *momsOther, double vtSqMinOther,
  struct gkyl_array *nuOut);
