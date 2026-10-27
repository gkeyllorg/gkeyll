#pragma once
#include <gkyl_gpython_core.h>

#define GPYTHON_VLASOV_API_VERSION 2
int gpython_vlasov_api_version(void);

/* samples contains n rows of [ns, nr, Ts, Tr, |B|]. Temperatures are energies.
 * Return 1 for invalid physical parameters, otherwise 0. */
int gpython_spitzer_coll_freq_calc_coulomb_log(
  size_t n, const double *samples, double ms, double mr, double qs, double qr, double eps0,
  double hbar, double *out
);

/* Maxwellian moment blocks are [n, u..., v_t^2]. Output is one modal field.
 * Return 1 for incompatible shapes or invalid parameters, otherwise 0. */
int gpython_spitzer_coll_freq(
  const gpython_basis *basis, const gpython_array *moms_self, const gpython_array *moms_other,
  const gpython_array *bmag, double ms, double mr, double qs, double qr, double eps0, double hbar,
  double nufrac, double vtmin_self, double vtmin_other, int num_quad, gpython_array *out
);
