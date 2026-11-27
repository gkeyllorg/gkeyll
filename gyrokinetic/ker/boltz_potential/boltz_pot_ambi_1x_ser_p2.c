#include <gkyl_boltz_potential_kernels.h>

GKYL_CU_DH void boltz_pot_ambi_reference_1x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out)
{
  const int nq = 1;
  double nref[1], phiref[1];
  double Te_eV = -T_e/q_e;
  double vte_sided = sqrt(T_e/(2.0*M_PI*m_e));
  for (int q=0; q<nq; ++q) {
    double density = lower[2*nq+q], flux = lower[q], area = lower[nq+q];
    double phi = 0.0;
    if (!gkyl_is_dfinite(density)) density = 0.0;
    if (density>0.0 && gkyl_is_dfinite(flux) && flux>0.0 && gkyl_is_dfinite(area) && area>0.0) {
      phi = lower[3*nq+q]+Te_eV*(log(density)+log(vte_sided)+log(area)-log(flux));
      if (!gkyl_is_dfinite(phi)) phi = 0.0;
    }
    nref[q] = density;
    phiref[q] = phi;
  }
  out[0] = 1.4142135623730951*nref[0];
  out[3] = 1.4142135623730951*phiref[0];
  out[1] = 0.0;
  out[4] = 0.0;
  out[2] = 0.0;
  out[5] = 0.0;
}

