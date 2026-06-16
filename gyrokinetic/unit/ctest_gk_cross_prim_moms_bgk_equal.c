// Test the gyrokinetic BGK cross primitive-moment calculator in the special
// case of two identical species.  The cross primitive moments (n, upar, vt^2)
// of a species relaxing against an identical species must be identical to its
// own primitive moments: the density is unchanged, and the cross drift speed
// and thermal speed equal the (common) input values.  Uses spatially-constant
// fields so the expected DG coefficients are exact.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk.h>
#include <math.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
run_equal_species(int cdim, int vdim, int poly_order,
  double n0, double upar0, double vtsq0)
{
  int ndim = cdim + vdim;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++) { lower[d] = -0.5; upper[d] = 0.5; cells[d] = 2; }
  for (int d=0; d<vdim; d++) {
    lower[cdim+d] = (d==0)? -6.0 : 0.0;
    upper[cdim+d] = (d==0)?  6.0 : 6.0;
    cells[cdim+d] = 8;
  }

  double confLower[GKYL_MAX_CDIM], confUpper[GKYL_MAX_CDIM];
  int confCells[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d]; confUpper[d] = upper[d]; confCells[d] = cells[d];
  }

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) confGhost[d] = 1;
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int nb = confBasis.num_basis;
  double fac = pow(sqrt(2.0), cdim); // 0th coeff scale of a constant

  // prim_moms layout: [den | upar | vtsq], each nb coefficients.
  struct gkyl_array *prim_moms = mkarr(3*nb, confLocal_ext.volume);
  gkyl_array_clear(prim_moms, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&confLocal_ext, iter.idx);
    double *d = gkyl_array_fetch(prim_moms, lidx);
    d[0*nb + 0] = n0*fac;
    d[1*nb + 0] = upar0*fac;
    d[2*nb + 0] = vtsq0*fac;
  }

  double mass = 1.5;
  double delta_sr = 1.0, betap1 = 1.0;

  gkyl_gyrokinetic_cross_prim_moms_bgk *calc =
    gkyl_gyrokinetic_cross_prim_moms_bgk_new(&basis, &confBasis, false);

  struct gkyl_array *out = mkarr(3*nb, confLocal_ext.volume);
  gkyl_array_clear(out, 0.0);

  // Self and other are identical.
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance(calc, &confLocal, delta_sr, betap1,
    mass, prim_moms, mass, prim_moms, out);

  gkyl_gyrokinetic_cross_prim_moms_bgk_release(calc);

  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&confLocal, iter.idx);
    const double *o = gkyl_array_cfetch(out, lidx);
    // Recover pointwise values from 0th coefficients.
    double n_out    = o[0*nb + 0]/fac;
    double upar_out = o[1*nb + 0]/fac;
    double vtsq_out = o[2*nb + 0]/fac;
    TEST_CHECK( gkyl_compare(n_out, n0, 1e-10*fabs(n0)) );
    TEST_CHECK( gkyl_compare(upar_out, upar0, 1e-9*fmax(fabs(upar0), 1.0)) );
    TEST_CHECK( gkyl_compare(vtsq_out, vtsq0, 1e-9*fabs(vtsq0)) );
    TEST_MSG("n=%g upar=%g vtsq=%g", n_out, upar_out, vtsq_out);
  }

  gkyl_array_release(prim_moms);
  gkyl_array_release(out);
}

void test_equal_1x1v_p1() { run_equal_species(1, 1, 1, 2.0e19, 1.0e4, 5.0e11); }
void test_equal_1x2v_p1() { run_equal_species(1, 2, 1, 3.0e19, -2.0e4, 8.0e11); }
void test_equal_1x1v_zeroflow_p1() { run_equal_species(1, 1, 1, 1.0e19, 0.0, 1.0e12); }

TEST_LIST = {
  { "equal_1x1v_p1", test_equal_1x1v_p1 },
  { "equal_1x2v_p1", test_equal_1x2v_p1 },
  { "equal_1x1v_zeroflow_p1", test_equal_1x1v_zeroflow_p1 },
  { NULL, NULL },
};
