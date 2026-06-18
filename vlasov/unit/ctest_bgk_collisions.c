// Tests for the BGK collision operator updater: out = nu*f_M - nu*f.
//
// We use fin = 0 to isolate exact, analytically-known behavior:
//   explicit: out = nufM,  and cflfreq += nu_d[0]*cellav_fac
//   implicit: out = nufM / (1 + nu_d[0]*cellav_fac*dt)
// where cellav_fac = 1/sqrt(2^cdim).
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void
test_bgk_explicit_1x1v()
{
  int poly_order = 1;
  int cdim = 1, vdim = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);

  // phase grid 1x1v
  double lower[] = {0.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {3, 4};
  double clower[] = {lower[0]}, cupper[] = {upper[0]};
  int ccells[] = {cells[0]};

  struct gkyl_rect_grid grid, cgrid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_rect_grid_init(&cgrid, 1, clower, cupper, ccells);

  struct gkyl_range prange, prange_ext, crange, crange_ext;
  int pghost[] = {0, 0};
  int cghost[] = {0};
  gkyl_create_grid_ranges(&grid, pghost, &prange_ext, &prange);
  gkyl_create_grid_ranges(&cgrid, cghost, &crange_ext, &crange);

  gkyl_bgk_collisions *up = gkyl_bgk_collisions_new(&cbasis, &pbasis, false);

  struct gkyl_array *nu = mkarr(cbasis.num_basis, crange_ext.volume);
  struct gkyl_array *nufM = mkarr(pbasis.num_basis, prange_ext.volume);
  struct gkyl_array *fin = mkarr(pbasis.num_basis, prange_ext.volume);
  struct gkyl_array *out = mkarr(pbasis.num_basis, prange_ext.volume);
  struct gkyl_array *cfl = mkarr(1, prange_ext.volume);

  // nu = constant: cell-average coefficient = nu0_phys * sqrt(2)^cdim.
  // Set the [0] component directly; pick a value.
  double nu0 = 2.5;
  gkyl_array_clear(nu, 0.0);
  gkyl_array_shiftc(nu, nu0, 0); // sets component 0 in every cell to nu0

  // nufM = known nonzero values; fin = 0.
  gkyl_array_clear(fin, 0.0);
  gkyl_array_clear(out, 0.0);
  gkyl_array_clear(cfl, 0.0);
  // fill nufM deterministically
  gkyl_array_clear(nufM, 0.0);
  gkyl_array_shiftc(nufM, 3.0, 0);
  if (pbasis.num_basis > 1) gkyl_array_shiftc(nufM, -1.5, 1);

  gkyl_bgk_collisions_advance(up, &crange, &prange, nu, nufM, fin,
    false, 0.0, out, cfl);

  // out should equal nufM exactly (since fin = 0).
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &prange);
  double cellav_fac = 1.0/sqrt(pow(2.0, cdim));
  while (gkyl_range_iter_next(&it)) {
    long ploc = gkyl_range_idx(&prange, it.idx);
    double *o = gkyl_array_fetch(out, ploc);
    double *fm = gkyl_array_fetch(nufM, ploc);
    for (int k=0; k<pbasis.num_basis; ++k)
      TEST_CHECK( gkyl_compare(fm[k], o[k], 1e-14) );
    // cfl frequency increment = nu0 * cellav_fac
    double *c = gkyl_array_fetch(cfl, ploc);
    TEST_CHECK( gkyl_compare(nu0*cellav_fac, c[0], 1e-14) );
  }

  gkyl_array_release(nu); gkyl_array_release(nufM); gkyl_array_release(fin);
  gkyl_array_release(out); gkyl_array_release(cfl);
  gkyl_bgk_collisions_release(up);
}

void
test_bgk_implicit_1x1v()
{
  int poly_order = 1;
  int cdim = 1, vdim = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);

  double lower[] = {0.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 4};
  double clower[] = {lower[0]}, cupper[] = {upper[0]};
  int ccells[] = {cells[0]};

  struct gkyl_rect_grid grid, cgrid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_rect_grid_init(&cgrid, 1, clower, cupper, ccells);

  struct gkyl_range prange, prange_ext, crange, crange_ext;
  int pghost[] = {0, 0};
  int cghost[] = {0};
  gkyl_create_grid_ranges(&grid, pghost, &prange_ext, &prange);
  gkyl_create_grid_ranges(&cgrid, cghost, &crange_ext, &crange);

  gkyl_bgk_collisions *up = gkyl_bgk_collisions_new(&cbasis, &pbasis, false);

  struct gkyl_array *nu = mkarr(cbasis.num_basis, crange_ext.volume);
  struct gkyl_array *nufM = mkarr(pbasis.num_basis, prange_ext.volume);
  struct gkyl_array *fin = mkarr(pbasis.num_basis, prange_ext.volume);
  struct gkyl_array *out = mkarr(pbasis.num_basis, prange_ext.volume);
  struct gkyl_array *cfl = mkarr(1, prange_ext.volume);

  double nu0 = 4.0, dt = 0.1;
  gkyl_array_clear(nu, 0.0);
  gkyl_array_shiftc(nu, nu0, 0);

  gkyl_array_clear(fin, 0.0);
  gkyl_array_clear(out, 0.0);
  gkyl_array_clear(cfl, 0.0);
  gkyl_array_clear(nufM, 0.0);
  gkyl_array_shiftc(nufM, 5.0, 0);

  gkyl_bgk_collisions_advance(up, &crange, &prange, nu, nufM, fin,
    true, dt, out, cfl);

  double cellav_fac = 1.0/sqrt(pow(2.0, cdim));
  double factor = 1.0/(1.0 + nu0*cellav_fac*dt);

  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &prange);
  while (gkyl_range_iter_next(&it)) {
    long ploc = gkyl_range_idx(&prange, it.idx);
    double *o = gkyl_array_fetch(out, ploc);
    double *fm = gkyl_array_fetch(nufM, ploc);
    for (int k=0; k<pbasis.num_basis; ++k)
      TEST_CHECK( gkyl_compare(factor*fm[k], o[k], 1e-13) );
  }

  gkyl_array_release(nu); gkyl_array_release(nufM); gkyl_array_release(fin);
  gkyl_array_release(out); gkyl_array_release(cfl);
  gkyl_bgk_collisions_release(up);
}

TEST_LIST = {
  { "bgk_explicit_1x1v", test_bgk_explicit_1x1v },
  { "bgk_implicit_1x1v", test_bgk_implicit_1x1v },
  { NULL, NULL },
};
