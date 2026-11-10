// Tests for the cross-species primitive moments of the Vlasov BGK operator: the pair of
// species must exchange no net momentum or energy, for a relative drift along any velocity
// direction (Maxwellian targets at (u_sr, vtsq_sr) with rho_s nu_sr = rho_r nu_rs).
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_cross_prim_moms_bgk.h>

struct cross_species {
  double mass, den, vtsq, flow[GKYL_MAX_VDIM];
};

static void
set_prim_moms(
  struct gkyl_array *arr, const struct gkyl_range *range, int nc, int vdim,
  const struct cross_species *sp
)
{
  double sq = pow(sqrt(2.0), range->ndim);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
    d[0] = sq * sp->den;
    for (int v = 0; v < vdim; ++v) {
      d[(1 + v) * nc] = sq * sp->flow[v];
    }
    d[(1 + vdim) * nc] = sq * sp->vtsq;
  }
}

// Rate of change of the momentum (per direction, in mom) and energy of one species relaxing at
// unit rho nu to the target held in cross, starting from the state sp.
static double
species_rates(
  const struct cross_species *sp, const double *cross, int cdim, int vdim, int nc, double *mom
)
{
  double sq = pow(sqrt(2.0), cdim), energy = 0.0;
  for (int v = 0; v < vdim; ++v) {
    double u_sr = cross[(1 + v) * nc] / sq;
    mom[v] += u_sr - sp->flow[v];
    energy += 0.5 * (u_sr * u_sr - sp->flow[v] * sp->flow[v]);
  }
  energy += 0.5 * vdim * (cross[(1 + vdim) * nc] / sq - sp->vtsq);
  return energy;
}

static void
test_pair_conservation(int vdim, const double *drift)
{
  struct gkyl_basis cbasis, pbasis;
  int cdim = 1, poly_order = 1;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
  int nc = cbasis.num_basis, ncomp = (vdim + 2) * nc;

  int lower[] = {1}, upper[] = {2};
  struct gkyl_range range;
  gkyl_range_init(&range, cdim, lower, upper);

  struct cross_species elc = {.mass = 1.0, .den = 1.0, .vtsq = 1.0};
  struct cross_species ion = {.mass = 100.0, .den = 1.0, .vtsq = 0.005};
  for (int v = 0; v < vdim; ++v) {
    ion.flow[v] = drift[v];
  }

  struct gkyl_array *pe = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *pi = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *cei = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *cie = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);
  set_prim_moms(pe, &range, nc, vdim, &elc);
  set_prim_moms(pi, &range, nc, vdim, &ion);

  gkyl_vlasov_cross_prim_moms_bgk *up =
    gkyl_vlasov_cross_prim_moms_bgk_new(&pbasis, &cbasis, false);
  double delta_sr = 2.0, betap1 = 1.0;
  gkyl_vlasov_cross_prim_moms_bgk_advance(
    up, &range, delta_sr, betap1, elc.mass, pe, ion.mass, pi, cei
  );
  gkyl_vlasov_cross_prim_moms_bgk_advance(
    up, &range, delta_sr, betap1, ion.mass, pi, elc.mass, pe, cie
  );

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double mom[GKYL_MAX_VDIM] = {0.0};
    double energy = species_rates(&elc, gkyl_array_cfetch(cei, loc), cdim, vdim, nc, mom);
    energy += species_rates(&ion, gkyl_array_cfetch(cie, loc), cdim, vdim, nc, mom);
    TEST_CHECK(gkyl_compare_double(energy, 0.0, 1e-13));
    TEST_MSG("vdim %d: pair energy rate %g", vdim, energy);
    for (int v = 0; v < vdim; ++v) {
      TEST_CHECK(gkyl_compare_double(mom[v], 0.0, 1e-13));
    }
  }

  gkyl_vlasov_cross_prim_moms_bgk_release(up);
  gkyl_array_release(pe);
  gkyl_array_release(pi);
  gkyl_array_release(cei);
  gkyl_array_release(cie);
}

// One relative drift along each velocity direction, and one along the diagonal.
static void
test_all_directions(int vdim)
{
  double a = 0.02;
  for (int d = 0; d < vdim; ++d) {
    double drift[GKYL_MAX_VDIM] = {0.0};
    drift[d] = a;
    test_pair_conservation(vdim, drift);
  }
  double diag[GKYL_MAX_VDIM] = {0.0};
  for (int d = 0; d < vdim; ++d) {
    diag[d] = a / sqrt(vdim);
  }
  test_pair_conservation(vdim, diag);
}

void
test_1x1v()
{
  test_all_directions(1);
}

void
test_1x2v()
{
  test_all_directions(2);
}

void
test_1x3v()
{
  test_all_directions(3);
}

TEST_LIST = {
  {"cross_prim_moms_bgk_1x1v", test_1x1v},
  {"cross_prim_moms_bgk_1x2v", test_1x2v},
  {"cross_prim_moms_bgk_1x3v", test_1x3v},
  {NULL, NULL}
};
