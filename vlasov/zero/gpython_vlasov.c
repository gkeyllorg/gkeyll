#include <gkyl_gpython_vlasov.h>
#include <gkyl_gpython_core_priv.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_range.h>
#include <limits.h>
#include <math.h>

int
gpython_vlasov_api_version(void)
{
  return GPYTHON_VLASOV_API_VERSION;
}

static bool
physical_parameters(double ms, double mr, double qs, double qr, double eps0, double hbar)
{
  return isfinite(ms) && ms > 0 && isfinite(mr) && mr > 0 && isfinite(qs) && qs != 0 &&
         isfinite(qr) && qr != 0 && isfinite(eps0) && eps0 > 0 && isfinite(hbar) && hbar >= 0;
}

int
gpython_spitzer_coll_freq_calc_coulomb_log(
  size_t n, const double *samples, double ms, double mr, double qs, double qr, double eps0,
  double hbar, double *out
)
{
  if (!physical_parameters(ms, mr, qs, qr, eps0, hbar)) {
    return 1;
  }
  for (size_t i = 0; i < n; ++i) {
    const double *p = samples + 5 * i;
    for (int j = 0; j < 5; ++j) {
      if (!isfinite(p[j]) || (j < 4 ? p[j] <= 0 : p[j] < 0)) {
        return 1;
      }
    }
    out[i] = gkyl_spitzer_coll_freq_calc_coulomb_log(
      p[0], p[1], ms, mr, p[2], p[3], qs, qr, p[4], eps0, hbar
    );
  }
  return 0;
}

int
gpython_spitzer_coll_freq(
  const gpython_basis *basis, const gpython_array *moms_self, const gpython_array *moms_other,
  const gpython_array *bmag, double ms, double mr, double qs, double qr, double eps0, double hbar,
  double nufrac, double vtmin_self, double vtmin_other, int num_quad, gpython_array *out
)
{
  const struct gkyl_basis *b = CBAS(basis);
  size_t nb = b->num_basis, size = CARR(moms_self)->size;
  if (!physical_parameters(ms, mr, qs, qr, eps0, hbar) || !isfinite(nufrac) || nufrac < 0 ||
      !isfinite(vtmin_self) || vtmin_self < 0 || !isfinite(vtmin_other) || vtmin_other < 0 ||
      num_quad < 1 || num_quad > 64 || b->ndim < 1 || b->ndim > 3 || size > INT_MAX || !size ||
      CARR(moms_self)->ncomp < 3 * nb || CARR(moms_self)->ncomp % nb ||
      CARR(moms_other)->ncomp != CARR(moms_self)->ncomp || CARR(moms_other)->size != size ||
      CARR(bmag)->ncomp != nb || CARR(bmag)->size != size || CARR(out)->ncomp != nb ||
      CARR(out)->size != size) {
    return 1;
  }
  int cells = (int)size;
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, &cells);
  gkyl_spitzer_coll_freq *up = gkyl_spitzer_coll_freq_new(b, num_quad, nufrac, eps0, hbar, false);
  gkyl_spitzer_coll_freq_advance(
    up, &range, CARR(bmag), qs, ms, CARR(moms_self), vtmin_self, qr, mr, CARR(moms_other),
    vtmin_other, ARR(out)
  );
  gkyl_spitzer_coll_freq_release(up);
  return 0;
}
