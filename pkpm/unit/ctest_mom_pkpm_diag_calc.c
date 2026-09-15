// Test computation of PKPM diagnostic moments (8 components) via gkyl_mom_calc.
//
// For the diagnostic kernel the (mass-weighted) moment block layout is:
//   0: rho            = \int F_0 dvpar
//   1: M1             = \int vpar F_0 dvpar
//   2: p_par          = \int vpar^2 F_0 dvpar
//   3: \int G_1 dvpar
//   4: q_par          = \int vpar^3 F_0 dvpar
//   5: \int vpar G_1 dvpar
//   6: r_parpar       = \int vpar^4 F_0 dvpar
//   7: \int vpar^2 G_1 dvpar
//
// With F_0 = a and G_1 = b constant over a symmetric velocity domain [-V,V]
// the odd-in-v moments (M1, q_par, \int vpar G_1) all vanish.
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_pkpm.h>
#include <gkyl_mom_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static const double AVAL = 1.25;
static const double BVAL = 0.60;

static void
eval_distf(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = AVAL;
  fout[1] = BVAL;
}

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void
test_mom_pkpm_diag_calc_1x1v_p1()
{
  int poly_order = 1;
  double mass = 2.0;
  double lower[] = {-1.0, -3.0}, upper[] = {1.0, 3.0};
  int cells[] = {2, 8};
  int cdim = 1, vdim = 1;
  int pdim = cdim+vdim;

  double V = upper[1]; // symmetric domain [-V,V]
  double Lv = upper[1] - lower[1];

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 2, eval_distf, NULL);
  struct gkyl_array *distf = mkarr(2*basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, distf);

  struct gkyl_mom_type *mt = gkyl_mom_pkpm_new(&confBasis, &basis, mass, true, false);
  TEST_CHECK( mt->num_mom == 8 );
  gkyl_mom_calc *mcalc = gkyl_mom_calc_new(&grid, mt, false);

  struct gkyl_array *mom = mkarr(mt->num_mom*confBasis.num_basis, confLocal_ext.volume);
  gkyl_mom_calc_advance(mcalc, &local, &confLocal, distf, mom);

  double sqrt2 = sqrt(2.0);
  double rho_exp     = mass * AVAL * Lv;
  double ppar_exp    = mass * AVAL * (2.0*V*V*V)/3.0;       // \int_{-V}^{V} v^2 dv = 2V^3/3
  double intG_exp    = mass * BVAL * Lv;
  double rparpar_exp = mass * AVAL * (2.0*V*V*V*V*V)/5.0;   // \int v^4 = 2V^5/5
  double intv2G_exp  = mass * BVAL * (2.0*V*V*V)/3.0;       // \int v^2 b dv

  for (int i=1; i<=cells[0]; ++i) {
    int cidx[] = {i};
    long lidx = gkyl_range_idx(&confLocal, cidx);
    double *m = gkyl_array_fetch(mom, lidx);
    double rho     = m[0]/sqrt2;
    double m1      = m[2]/sqrt2;
    double ppar    = m[4]/sqrt2;
    double intG    = m[6]/sqrt2;
    double qpar    = m[8]/sqrt2;
    double intvG   = m[10]/sqrt2;
    double rparpar = m[12]/sqrt2;
    double intv2G  = m[14]/sqrt2;

    TEST_CHECK( gkyl_compare(rho_exp,     rho,     1e-12) );
    TEST_CHECK( gkyl_compare(0.0,         m1,      1e-12) );  // odd moment
    TEST_CHECK( gkyl_compare(ppar_exp,    ppar,    1e-12) );
    TEST_CHECK( gkyl_compare(intG_exp,    intG,    1e-12) );
    TEST_CHECK( gkyl_compare(0.0,         qpar,    1e-12) );  // odd moment
    TEST_CHECK( gkyl_compare(0.0,         intvG,   1e-12) );  // odd moment
    TEST_CHECK( gkyl_compare(rparpar_exp, rparpar, 1e-11) );
    TEST_CHECK( gkyl_compare(intv2G_exp,  intv2G,  1e-12) );
  }

  gkyl_array_release(mom);
  gkyl_array_release(distf);
  gkyl_mom_calc_release(mcalc);
  gkyl_mom_type_release(mt);
  gkyl_proj_on_basis_release(proj);
}

TEST_LIST = {
  { "mom_pkpm_diag_calc_1x1v_p1", test_mom_pkpm_diag_calc_1x1v_p1 },
  { NULL, NULL },
};
