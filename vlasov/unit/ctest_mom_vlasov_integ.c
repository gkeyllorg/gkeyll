// Verify Vlasov integrated moments (M0, M1i, M2) of analytically-known
// distribution functions, summed over the whole domain.
//
// For f(x,vx) = A (constant):
//   int f dx dvx               = A * Lx * Lv
//   int vx f dx dvx            = A * Lx * 0       (symmetric vel domain) = 0
//   int vx^2 f dx dvx          = A * Lx * (vhi^3 - vlo^3)/3
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static const double AMP = 0.75;

static void
eval_const(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = AMP;
}

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void
test_integ_1x1v_const()
{
  int poly_order = 2;
  // config x in [-1,3] (Lx=4), velocity vx in [-3,3] (Lv=6).
  double lower[] = {-1.0, -3.0}, upper[] = {3.0, 3.0};
  int cells[] = {6, 8};
  int cdim = 1, vdim = 1, ndim = 2;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_range local, local_ext, confLocal, confLocal_ext;
  int ghost[] = {0, 0};
  int confGhost[] = {0};
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_const, NULL);

  struct gkyl_array *distf = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, distf);

  // integrated moments: vdim+2 components = [M0, M1x, M2]
  struct gkyl_mom_type *int_t = gkyl_int_mom_vlasov_new(&confBasis, &basis,
    GKYL_F_MOMENT_M0M1M2, false);
  gkyl_mom_calc *intcalc = gkyl_mom_calc_new(&grid, int_t, false);

  struct gkyl_array *int_mom = mkarr(vdim+2, confLocal_ext.volume);
  gkyl_mom_calc_advance(intcalc, &local, &confLocal, distf, int_mom);

  double red[vdim+2];
  gkyl_array_reduce_range(red, int_mom, GKYL_SUM, &confLocal);

  double Lx = 4.0, vlo = -3.0, vhi = 3.0, Lv = vhi - vlo;
  double exp_M0 = AMP*Lx*Lv;                       // 0.75*4*6 = 18
  double exp_M1 = 0.0;                             // symmetric vel domain
  double exp_M2 = AMP*Lx*(vhi*vhi*vhi - vlo*vlo*vlo)/3.0; // 0.75*4*18 = 54

  TEST_CHECK( gkyl_compare( exp_M0, red[0], 1e-11) );
  TEST_MSG("M0: expected %g got %g", exp_M0, red[0]);
  TEST_CHECK( gkyl_compare( exp_M1, red[1], 1e-11) );
  TEST_MSG("M1: expected %g got %g", exp_M1, red[1]);
  TEST_CHECK( gkyl_compare( exp_M2, red[2], 1e-11) );
  TEST_MSG("M2: expected %g got %g", exp_M2, red[2]);

  gkyl_array_release(int_mom);
  gkyl_mom_calc_release(intcalc);
  gkyl_mom_type_release(int_t);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_release(distf);
}

void
test_integ_1x2v_const()
{
  int poly_order = 2;
  // x in [0,2] (Lx=2), vx,vy in [-2,2] (Lv=4 each).
  double lower[] = {0.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {4, 6, 6};
  int cdim = 1, vdim = 2, ndim = 3;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_range local, local_ext, confLocal, confLocal_ext;
  int ghost[] = {0, 0, 0};
  int confGhost[] = {0};
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_const, NULL);

  struct gkyl_array *distf = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, distf);

  struct gkyl_mom_type *int_t = gkyl_int_mom_vlasov_new(&confBasis, &basis,
    GKYL_F_MOMENT_M0M1M2, false);
  gkyl_mom_calc *intcalc = gkyl_mom_calc_new(&grid, int_t, false);

  struct gkyl_array *int_mom = mkarr(vdim+2, confLocal_ext.volume);
  gkyl_mom_calc_advance(intcalc, &local, &confLocal, distf, int_mom);

  double red[vdim+2];
  gkyl_array_reduce_range(red, int_mom, GKYL_SUM, &confLocal);

  double Lx = 2.0, Lvx = 4.0, Lvy = 4.0;
  double phaseVol = AMP*Lx*Lvx*Lvy;     // 0.75*2*4*4 = 24
  // M0 = integral of f over velocity space (per unit) summed
  double exp_M0 = phaseVol;
  double exp_M1x = 0.0, exp_M1y = 0.0;  // symmetric velocity domains
  // M2 = int (vx^2+vy^2) f. Each: AMP*Lx * (v^3lim) * Lother
  // int vx^2 dvx = (2^3 - (-2)^3)/3 = 16/3; times Lvy=4; times Lx=2; times AMP
  double iv2 = (2.0*2.0*2.0 - (-2.0)*(-2.0)*(-2.0))/3.0; // 16/3
  double exp_M2 = AMP*Lx*(iv2*Lvy + iv2*Lvx); // both directions

  TEST_CHECK( gkyl_compare( exp_M0, red[0], 1e-11) );
  TEST_MSG("M0: expected %g got %g", exp_M0, red[0]);
  TEST_CHECK( gkyl_compare( exp_M1x, red[1], 1e-11) );
  TEST_CHECK( gkyl_compare( exp_M1y, red[2], 1e-11) );
  TEST_CHECK( gkyl_compare( exp_M2, red[3], 1e-11) );
  TEST_MSG("M2: expected %g got %g", exp_M2, red[3]);

  gkyl_array_release(int_mom);
  gkyl_mom_calc_release(intcalc);
  gkyl_mom_type_release(int_t);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_release(distf);
}

TEST_LIST = {
  { "integ_1x1v_const", test_integ_1x1v_const },
  { "integ_1x2v_const", test_integ_1x2v_const },
  { NULL, NULL },
};
