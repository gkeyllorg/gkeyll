/* -*- c++ -*- */

extern "C" {
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_spitzer_coll_freq_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

GKYL_CU_D void
calc_nu_cu(const struct gkyl_array* GKYL_RESTRICT basis_at_ords, const struct gkyl_array* GKYL_RESTRICT weights,
  const double* vtSqSelf_d, double vtSqMinSelf, const double* GKYL_RESTRICT m0Other_d,
  const double* vtSqOther_d, double vtSqMinOther, double normNu, long linidx,
  struct gkyl_array* GKYL_RESTRICT nuOut)
{
  // Perform the multiplication of normNu*n_r/(v_ts^2+v_tr^2)^(3/2) via
  // quadrature in one cell.
  int num_basis = basis_at_ords->ncomp;
  int tot_quad = basis_at_ords->size;

  double *nuOut_d = (double *) gkyl_array_fetch(nuOut, linidx);
  for (int k=0; k<num_basis; ++k) nuOut_d[k] = 0.0;

  // Compute expansion coefficients of nuOut using quadrature.
  const double *w_d = (const double *) weights->data;
  const double *bo_d = (const double *) basis_at_ords->data;

  for (int n=0; n<tot_quad; ++n) {

    const double *b_ord = (const double *) gkyl_array_cfetch(basis_at_ords, n);

    // Evaluate densities and thermal speeds (squared) at quad point.
    double vtSqSelf_q=0., m0Other_q=0., vtSqOther_q=0.;
    for (int k=0; k<num_basis; ++k) {
      vtSqSelf_q += vtSqSelf_d[k]*b_ord[k];
      m0Other_q += m0Other_d[k]*b_ord[k];
      vtSqOther_q += vtSqOther_d[k]*b_ord[k];
    }

    double nu_o;
    if (m0Other_q<0.) {
      nu_o = 0.;
    } else if ((vtSqSelf_q < vtSqMinSelf) && (vtSqOther_q < vtSqMinOther)) {
      nu_o = normNu*m0Other_q/pow(sqrt(vtSqMinSelf+vtSqMinOther),3);
    } else if (vtSqSelf_q < vtSqMinSelf) {
      nu_o = normNu*m0Other_q/pow(sqrt(vtSqMinSelf+vtSqOther_q),3);
    } else if (vtSqOther_q < vtSqMinOther) {
      nu_o = normNu*m0Other_q/pow(sqrt(vtSqSelf_q+vtSqMinOther),3);
    } else {
      nu_o = normNu*m0Other_q/pow(sqrt(vtSqSelf_q+vtSqOther_q),3);
    }

    double tmp = w_d[n]*nu_o;
    for (int k=0; k<num_basis; ++k)
      nuOut_d[k] += tmp*bo_d[k+num_basis*n];
  }
}


__global__ static void
gkyl_spitzer_coll_freq_advance_normnu_cu_ker(const struct gkyl_range range,
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, const struct gkyl_array* GKYL_RESTRICT weights,
  const struct gkyl_array* momsSelf, double vtSqMinSelf,
  const struct gkyl_array* GKYL_RESTRICT momsOther, double vtSqMinOther,
  double normNu, struct gkyl_array* GKYL_RESTRICT nuOut)
{
  int vtsq_idx = momsSelf->ncomp-basis_at_ords->ncomp;
  int idx[3];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, tid, idx);
    long linidx = gkyl_range_idx(&range, idx);

    const double *momsSelf_d = (const double *) gkyl_array_cfetch(momsSelf, linidx);
    const double *momsOther_d = (const double *) gkyl_array_cfetch(momsOther, linidx);

    const double *vtSqSelf_d = &momsSelf_d[vtsq_idx];
    const double *m0Other_d = momsOther_d;
    const double *vtSqOther_d = &momsOther_d[vtsq_idx];

    calc_nu_cu(basis_at_ords, weights, vtSqSelf_d, vtSqMinSelf, m0Other_d, vtSqOther_d, vtSqMinOther, normNu, linidx, nuOut);
  }
}

__global__ static void
gkyl_spitzer_coll_freq_advance_cu_ker(const struct gkyl_range range,
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, const struct gkyl_array* GKYL_RESTRICT weights,
  double nufraceps0_fac, double cellav_fac, double r4pieps0_fac, double hbar_fac, double eps0,
  const struct gkyl_array* GKYL_RESTRICT bmag,
  double qSelf, double mSelf, const struct gkyl_array* momsSelf, double vtSqMinSelf,
  double qOther, double mOther, const struct gkyl_array* momsOther, double vtSqMinOther,
  struct gkyl_array* GKYL_RESTRICT nuOut)
{
  int vtsq_idx = momsSelf->ncomp-basis_at_ords->ncomp;
  double mReduced = 1./(1./mSelf+1./mOther);
  double timeConstFac = nufraceps0_fac*pow(qSelf*qOther,2)/(mSelf*mReduced);

  int idx[3];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, tid, idx);
    long linidx = gkyl_range_idx(&range, idx);

    const double *bmag_d = (const double *) gkyl_array_cfetch(bmag, linidx);
    const double *momsSelf_d = (const double *) gkyl_array_cfetch(momsSelf, linidx);
    const double *momsOther_d = (const double *) gkyl_array_cfetch(momsOther, linidx);

    const double *m0Self_d = momsSelf_d;
    const double *vtSqSelf_d = &momsSelf_d[vtsq_idx];
    const double *m0Other_d = momsOther_d;
    const double *vtSqOther_d = &momsOther_d[vtsq_idx];

    // Compute the Coulomb logarithm using cell-average values.
    double bmagAv = bmag_d[0]*cellav_fac;
    double m0SelfAv    =    m0Self_d[0] < 0.? 1.e-14 : m0Self_d[0]*cellav_fac;
    double vtSqSelfAv  =  vtSqSelf_d[0] < vtSqMinSelf?  vtSqMinSelf*cellav_fac  : vtSqSelf_d[0]*cellav_fac;
    double m0OtherAv   =   m0Other_d[0] < 0.? 1.e-14 : m0Other_d[0]*cellav_fac;
    double vtSqOtherAv = vtSqOther_d[0] < vtSqMinOther? vtSqMinOther*cellav_fac : vtSqOther_d[0]*cellav_fac;

    double omegaSqSumSelf  = m0SelfAv*pow(qSelf,2)/(eps0*mSelf)+pow(qSelf*bmagAv/mSelf,2);
    double omegaSqSumOther = m0OtherAv*pow(qOther,2)/(eps0*mOther)+pow(qOther*bmagAv/mOther,2);

    double rmaxSumSelf = omegaSqSumSelf/(vtSqSelfAv+3.*vtSqSelfAv)+omegaSqSumOther/(vtSqOtherAv+3.*vtSqSelfAv);
    double rmaxSumOther = omegaSqSumSelf/(vtSqSelfAv+3.*vtSqOtherAv)+omegaSqSumOther/(vtSqOtherAv+3.*vtSqOtherAv);

    double rmaxSelf = 1./sqrt(rmaxSumSelf);
    double rmaxOther = 1./sqrt(rmaxSumOther);

    double uRelSq = 3.*(vtSqOtherAv+vtSqSelfAv);

    double rMin = GKYL_MAX2(fabs(qSelf*qOther)*r4pieps0_fac/(mReduced*uRelSq), hbar_fac/(mReduced*sqrt(uRelSq)));

    double logLambda = 0.5*(0.5*log(1.+pow(rmaxSelf/rMin,2))+0.5*log(1.+pow(rmaxOther/rMin,2)));

    // Normalized nu (nu missing density and temperature factors).
    double normNu = timeConstFac*logLambda;

    calc_nu_cu(basis_at_ords, weights, vtSqSelf_d, vtSqMinSelf, m0Other_d, vtSqOther_d, vtSqMinOther, normNu, linidx, nuOut);
  }
}

void
gkyl_spitzer_coll_freq_advance_normnu_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *momsSelf, double vtSqMinSelf,
  const struct gkyl_array *momsOther, double vtSqMinOther,
  double normNu, struct gkyl_array *nuOut)
{
  int nblocks = range->nblocks, nthreads = range->nthreads;
  gkyl_spitzer_coll_freq_advance_normnu_cu_ker<<<nblocks, nthreads>>>
    (*range, up->basis_at_ords->on_dev, up->weights->on_dev, momsSelf->on_dev, vtSqMinSelf,
     momsOther->on_dev, vtSqMinOther, normNu, nuOut->on_dev);
}

void
gkyl_spitzer_coll_freq_advance_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *bmag,
  double qSelf, double mSelf, const struct gkyl_array *momsSelf, double vtSqMinSelf,
  double qOther, double mOther, const struct gkyl_array *momsOther, double vtSqMinOther,
  struct gkyl_array *nuOut)
{
  int nblocks = range->nblocks, nthreads = range->nthreads;
  gkyl_spitzer_coll_freq_advance_cu_ker<<<nblocks, nthreads>>>
    (*range, up->basis_at_ords->on_dev, up->weights->on_dev,
     up->nufraceps0_fac, up->cellav_fac, up->r4pieps0_fac, up->hbar_fac, up->eps0, bmag->on_dev,
     qSelf, mSelf, momsSelf->on_dev, vtSqMinSelf,
     qOther, mOther, momsOther->on_dev, vtSqMinOther, nuOut->on_dev);
}
