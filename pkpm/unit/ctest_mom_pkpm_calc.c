// Test actual computation of PKPM moments via gkyl_mom_calc.
//
// The PKPM distribution f is a 2-component vector [F_0, G_1] (each with
// num_basis components).  The non-diagnostic moment kernel produces
//   out = [rho, p_par, p_perp, M1]   (mass weighted)
// where
//   rho     = mass * \int F_0 dvpar
//   p_par   = mass * \int vpar^2 F_0 dvpar
//   p_perp  = mass * \int G_1   dvpar
//   M1      = mass * \int vpar  F_0 dvpar
//
// We project a distribution that is constant in velocity:
//   F_0 = a,  G_1 = b
// over a symmetric velocity domain [-V, V] (so M1 = 0), and check the
// resulting moments against the exact analytic integrals.
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

static const double AVAL = 0.75; // F_0 amplitude
static const double BVAL = 0.40; // G_1 amplitude

static void
eval_distf(double t, const double *xn, double *restrict fout, void *ctx)
{
  // 2-component distribution: F_0 (constant), G_1 (constant).
  fout[0] = AVAL;
  fout[1] = BVAL;
}

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void
test_mom_pkpm_calc_1x1v_p1()
{
  int poly_order = 1;
  double mass = 1.5;
  double lower[] = {-2.0, -3.0}, upper[] = {2.0, 3.0};
  int cells[] = {4, 8};
  int cdim = 1, vdim = 1;
  int pdim = cdim+vdim;

  double Vlo = lower[1], Vhi = upper[1];
  double Lv = Vhi - Vlo;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // PKPM p1 uses hybrid phase basis (p=2 in velocity space).
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Distribution function: 2 components per basis function (F_0, G_1).
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 2, eval_distf, NULL);
  struct gkyl_array *distf = mkarr(2*basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, distf);

  struct gkyl_mom_type *mt = gkyl_mom_pkpm_new(&confBasis, &basis, mass, false, false);
  TEST_CHECK( mt->num_mom == 4 );
  gkyl_mom_calc *mcalc = gkyl_mom_calc_new(&grid, mt, false);

  struct gkyl_array *mom = mkarr(mt->num_mom*confBasis.num_basis, confLocal_ext.volume);
  gkyl_mom_calc_advance(mcalc, &local, &confLocal, distf, mom);

  // Exact expectations (cell-average component, index 0 of each moment block;
  // recall the orthonormal cell-average basis is 1/sqrt(2^cdim) = 1/sqrt(2)).
  double sqrt2 = sqrt(2.0);
  double rho_exp    = mass * AVAL * Lv;            // \int a dv
  double ppar_exp   = mass * AVAL * (Vhi*Vhi*Vhi - Vlo*Vlo*Vlo)/3.0; // \int v^2 a dv
  double pperp_exp  = mass * BVAL * Lv;            // \int b dv
  double m1_exp     = 0.0;                         // symmetric => 0

  // moment block layout: [rho(2), p_par(2), p_perp(2), M1(2)] per conf cell.
  for (int i=1; i<=cells[0]; ++i) {
    int cidx[] = {i};
    long lidx = gkyl_range_idx(&confLocal, cidx);
    double *m = gkyl_array_fetch(mom, lidx);
    // cell-average physical value = m[block*2 + 0]/sqrt(2) for p1 1D conf basis.
    double rho   = m[0]/sqrt2;
    double ppar  = m[2]/sqrt2;
    double pperp = m[4]/sqrt2;
    double m1    = m[6]/sqrt2;
    TEST_CHECK( gkyl_compare(rho_exp,   rho,   1e-12) );
    TEST_CHECK( gkyl_compare(ppar_exp,  ppar,  1e-12) );
    TEST_CHECK( gkyl_compare(pperp_exp, pperp, 1e-12) );
    TEST_CHECK( gkyl_compare(m1_exp,    m1,    1e-12) );
    // x-slope (component 1 of each block) should be ~0 for uniform profile.
    TEST_CHECK( gkyl_compare(0.0, m[1], 1e-12) );
    TEST_CHECK( gkyl_compare(0.0, m[5], 1e-12) );
  }

  gkyl_array_release(mom);
  gkyl_array_release(distf);
  gkyl_mom_calc_release(mcalc);
  gkyl_mom_type_release(mt);
  gkyl_proj_on_basis_release(proj);
}

TEST_LIST = {
  { "mom_pkpm_calc_1x1v_p1", test_mom_pkpm_calc_1x1v_p1 },
  { NULL, NULL },
};
