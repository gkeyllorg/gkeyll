#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_spitzer_coll_freq_priv.h>
#include <gkyl_range.h>

// Calculate the plasma frequency
double plasma_frequency(double n, double m, double eps0, double eV)
{
  return sqrt(n*eV*eV/m/eps0);
}

// Calculate the Coulomb Logarithm
double coulomb_log(double ns, double nr, double ms, double mr, double Ts, double Tr,
  double qs, double qr, double bmag, double eps0, double hbar, double eV)
{
  double vts = sqrt(Ts/ms); // Thermal velocity for species s
  double vtr = sqrt(Tr/mr);  // Thermal velocity for species r
  double wps = plasma_frequency(ns, ms, eps0, eV); // Plasma Frequency for species s
  double wpr = plasma_frequency(nr, mr, eps0, eV); // Plasma frequency for species r
  double wcs = qs*bmag/ms; // Cyclotron frequency for species s
  double wcr = qr*bmag/mr; // Cyclotron frequency for species r
  double inner1 = (wps*wps + wcs*wcs)/(Ts/ms + 3*Ts/ms) + (wpr*wpr + wcr*wcr)/(Tr/mr + 3*Ts/ms);
  double u = 3*(vts*vts + vtr*vtr); // Relative velocity
  double msr = ms*mr/(ms+mr); // Reduced mass
  double inner2 = fmax(fabs(qs*qr)/(4*M_PI*eps0*msr*u*u), hbar/(2*sqrt(eV)*msr*u));
  double inner = (1/inner1)*(1/inner2/inner2) + 1;
  return 0.5*log(inner);
}

double gkyl_calc_Morse_alpha_E_const(double ns, double nr, double ms, double mr, double qs, double qr,
  double Ts, double Tr, double bmag, double eps0, double hbar, double eV)
{
  double clog = 0.5*( coulomb_log(ns, nr, ms, mr, Ts, Tr, qs, qr, bmag, eps0, hbar, eV)
                     +coulomb_log(nr, ns, mr, ms, Tr, Ts, qr, qs, bmag, eps0, hbar, eV) );
  return 2.0*pow(qs*qr,2)*clog/(3.0*pow(2.0*M_PI,1.5)*pow(eps0,2)*ms*mr);
}

double gkyl_calc_Morse_alpha_E(double ns, double nr, double ms, double mr, double qs, double qr,
  double Ts, double Tr, double bmag, double eps0, double hbar, double eV)
{
  double alpha_E_fac = gkyl_calc_Morse_alpha_E_const(ns, nr, ms, mr, qs, qr, Ts, Tr, bmag, eps0, hbar, eV);
  return alpha_E_fac * ns * nr / pow(sqrt(Ts/ms+Tr/mr),3.0);
}

// create range to loop over quadrature points.
static inline struct gkyl_range get_qrange(int dim, int num_quad) {
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<dim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

// Sets weights and basis functions at ords. Returns total
// number of quadrature nodes.
static int
init_quad_values(const struct gkyl_basis *basis, int num_quad,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  double ordinates1[num_quad], weights1[num_quad];

  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }

  struct gkyl_range qrange = get_qrange(ndim, num_quad);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  } else {
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);

    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];

    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));

  if (use_gpu)
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  else
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);

  // copy host array to device array
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

gkyl_spitzer_coll_freq*
gkyl_spitzer_coll_freq_new(const struct gkyl_basis *basis, int num_quad,
  double nufrac, double eps0, double hbar, bool use_gpu)
{
  struct gkyl_spitzer_coll_freq *up = gkyl_malloc(sizeof(struct gkyl_spitzer_coll_freq));

  up->ndim = basis->ndim;
  up->num_quad = num_quad;
  up->num_basis = basis->num_basis;
  up->use_gpu = use_gpu;

  // initialize data needed for quadrature
  up->tot_quad = init_quad_values(basis, num_quad, &up->weights,
    &up->basis_at_ords, use_gpu);

  if (up->use_gpu)
    up->fun_at_ords = NULL;
  else
    up->fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad); // Only used in CPU implementation.

  up->eps0 = eps0;
  // Pre-compute time-independent factors for the case in which we calculate
  // the collision frequency from scratch (instead of normnu).
  up->hbar_fac = hbar/(2.0*exp(0.5));
  up->r4pieps0_fac = 1./(4.*M_PI*eps0);
  up->nufraceps0_fac = nufrac/(3.*sqrt(pow(2.*M_PI,3))*pow(eps0,2));
  up->cellav_fac = 1./pow(sqrt(2.),up->ndim);
  
  return up;
}

static void
proj_on_basis(const gkyl_spitzer_coll_freq *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;

  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

void
calc_nu(const gkyl_spitzer_coll_freq *up, struct gkyl_range qrange, const double *vtSqSelf_d,
  double vtSqMinSelf, const double *m0Other_d, const double *vtSqOther_d, double vtSqMinOther,
  double normNu, long linidx, struct gkyl_array *nuOut)
{
  // Perform the multiplication of normNu*n_r/(v_ts^2+v_tr^2)^(3/2) via
  // quadrature in one cell.
  struct gkyl_range_iter qiter;
  gkyl_range_iter_init(&qiter, &qrange);
  while (gkyl_range_iter_next(&qiter)) {
    
    int qidx = gkyl_range_idx(&qrange, qiter.idx);

    // Evaluate densities and thermal speeds (squared) at quad point.
    const double *b_ord = gkyl_array_cfetch(up->basis_at_ords, qidx);
    double vtSqSelf_q=0., m0Other_q=0., vtSqOther_q=0.;
    for (int k=0; k<up->num_basis; ++k) {
      vtSqSelf_q += vtSqSelf_d[k]*b_ord[k];
      m0Other_q += m0Other_d[k]*b_ord[k];
      vtSqOther_q += vtSqOther_d[k]*b_ord[k];
    }

    double *fq = gkyl_array_fetch(up->fun_at_ords, qidx);

    if (m0Other_q<0.) {
      fq[0] = 0.;
    } else if ((vtSqSelf_q < vtSqMinSelf) && (vtSqOther_q < vtSqMinOther)) {
      fq[0] = normNu*m0Other_q/pow(sqrt(vtSqMinSelf+vtSqMinOther),3);
    } else if (vtSqSelf_q < vtSqMinSelf) {
      fq[0] = normNu*m0Other_q/pow(sqrt(vtSqMinSelf+vtSqOther_q),3);
    } else if (vtSqOther_q < vtSqMinOther) {
      fq[0] = normNu*m0Other_q/pow(sqrt(vtSqSelf_q+vtSqMinOther),3);
    } else {
      fq[0] = normNu*m0Other_q/pow(sqrt(vtSqSelf_q+vtSqOther_q),3);
    }
  }

  // compute expansion coefficients of Maxwellian on basis
  proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(nuOut, linidx));
}

void
gkyl_spitzer_coll_freq_advance_normnu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *momsSelf, double vtSqMinSelf,
  const struct gkyl_array *momsOther, double vtSqMinOther,
  double normNu, struct gkyl_array *nuOut)
{
  // Scale project normNu*n_r/(v_ts^2+v_tr^2)^(3/2) onto the basis using
  // quadrature.

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_spitzer_coll_freq_advance_normnu_cu(up, range, momsSelf, 
      vtSqMinSelf, momsOther, vtSqMinOther, normNu, nuOut);
#endif

  int vtsq_idx = momsSelf->ncomp-up->num_basis;

  // Create range to loop over quadrature points.
  struct gkyl_range qrange = get_qrange(up->ndim, up->num_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(range, iter.idx);

    const double *momsSelf_d = gkyl_array_cfetch(momsSelf, linidx);
    const double *momsOther_d = gkyl_array_cfetch(momsOther, linidx);

    const double *vtSqSelf_d = &momsSelf_d[vtsq_idx];
    const double *m0Other_d = momsOther_d;
    const double *vtSqOther_d = &momsOther_d[vtsq_idx];

    calc_nu(up, qrange, vtSqSelf_d, vtSqMinSelf, m0Other_d, vtSqOther_d, vtSqMinOther, normNu, linidx, nuOut);
  }

}

void
gkyl_spitzer_coll_freq_advance(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *bmag,  
  double qSelf, double mSelf, const struct gkyl_array *momsSelf, double vtSqMinSelf,
  double qOther, double mOther, const struct gkyl_array *momsOther, double vtSqMinOther,
  struct gkyl_array *nuOut)
{
  // Compute the Spitzer-like collision frequency
  //   nu_sr = nu_frac * (n_r/m_s)*(1/m_s+1/m_r)
  //     * (q_s^2*q_r^2*log(Lambda_sr)/(3*(2*pi)^(3/2)*epsilon_0^2))
  //     * (1/(v_ts^2+v_tr^2)^(3/2))
  // where log(Lambda_sr) is the Coulomb logarithm (see Gkeyll docs).

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_spitzer_coll_freq_advance_cu(up, range, bmag, qSelf, mSelf, momsSelf,
      vtSqMinSelf, qOther, mOther, momsOther, vtSqMinOther, nuOut);
#endif

  int vtsq_idx = momsSelf->ncomp-up->num_basis;
  double mReduced = 1./(1./mSelf+1./mOther);
  double timeConstFac = up->nufraceps0_fac*pow(qSelf*qOther,2)/(mSelf*mReduced);

  // Create range to loop over quadrature points.
  struct gkyl_range qrange = get_qrange(up->ndim, up->num_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(range, iter.idx);

    const double *bmag_d = gkyl_array_cfetch(bmag, linidx);
    const double *momsSelf_d = gkyl_array_cfetch(momsSelf, linidx);
    const double *momsOther_d = gkyl_array_cfetch(momsOther, linidx);

    const double *m0Self_d = momsSelf_d;
    const double *vtSqSelf_d = &momsSelf_d[vtsq_idx];
    const double *m0Other_d = momsOther_d;
    const double *vtSqOther_d = &momsOther_d[vtsq_idx];

    // Compute the Coulomb logarithm using cell-average values.
    double bmagAv = bmag_d[0]*up->cellav_fac;
    double m0SelfAv    =    m0Self_d[0] < 0.? 1.e-14 : m0Self_d[0]*up->cellav_fac;
    double vtSqSelfAv  =  vtSqSelf_d[0] < vtSqMinSelf?  vtSqMinSelf*up->cellav_fac  : vtSqSelf_d[0]*up->cellav_fac;
    double m0OtherAv   =   m0Other_d[0] < 0.? 1.e-14 : m0Other_d[0]*up->cellav_fac;
    double vtSqOtherAv = vtSqOther_d[0] < vtSqMinOther? vtSqMinOther*up->cellav_fac : vtSqOther_d[0]*up->cellav_fac;

    double omegaSqSumSelf  = m0SelfAv*pow(qSelf,2)/(up->eps0*mSelf)+pow(qSelf*bmagAv/mSelf,2);
    double omegaSqSumOther = m0OtherAv*pow(qOther,2)/(up->eps0*mOther)+pow(qOther*bmagAv/mOther,2);

    double rmaxSumSelf = omegaSqSumSelf/(vtSqSelfAv+3.*vtSqSelfAv)+omegaSqSumOther/(vtSqOtherAv+3.*vtSqSelfAv);
    double rmaxSumOther = omegaSqSumSelf/(vtSqSelfAv+3.*vtSqOtherAv)+omegaSqSumOther/(vtSqOtherAv+3.*vtSqOtherAv);

    double rmaxSelf = 1./sqrt(rmaxSumSelf);
    double rmaxOther = 1./sqrt(rmaxSumOther);

    double uRelSq = 3.*(vtSqOtherAv+vtSqSelfAv);

    double rMin = GKYL_MAX2(fabs(qSelf*qOther)*up->r4pieps0_fac/(mReduced*uRelSq), up->hbar_fac/(mReduced*sqrt(uRelSq)));

    double logLambda = 0.5*(0.5*log(1.+pow(rmaxSelf/rMin,2))+0.5*log(1.+pow(rmaxOther/rMin,2)));

    // Normalized nu (nu missing density and temperature factors).
    double normNu = timeConstFac*logLambda;

    calc_nu(up, qrange, vtSqSelf_d, vtSqMinSelf, m0Other_d, vtSqOther_d, vtSqMinOther, normNu, linidx, nuOut);
  }

}

void
gkyl_spitzer_coll_freq_release(gkyl_spitzer_coll_freq* up)
{
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  if (!up->use_gpu)
    gkyl_array_release(up->fun_at_ords);
  gkyl_free(up);
}
