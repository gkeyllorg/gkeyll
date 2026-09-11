#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <math.h>
#include <stdio.h>

// Allocate cu_dev array
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Create ghost ranges
struct skin_ghost_ranges
{
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
                       const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d = 0; d < ndim; ++d){
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
                           d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
                           d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_den_e(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_den_i(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_vtsq_e(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double eV = 1.602e-19;
  double me = 9.11e-31;
  double Te = 30.0*eV;
  fout[0] = Te/me;
}
void eval_vtsq_i(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double eV = 1.602e-19;
  double mi = 1.67e-27;
  double Ti = 10.0*eV;
  fout[0] = Ti/mi;
}
void eval_upar_e(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq[1];
  eval_vtsq_e(t, xn, vtsq, ctx);
  double vt = sqrt(vtsq[0]); 
  fout[0] = 0.01*vt;
}
void eval_upar_i(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq[1];
  eval_vtsq_i(t, xn, vtsq, ctx);
  double vt = sqrt(vtsq[0]); 
  fout[0] = 0.01*vt;
}
void eval_nu_ei(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double eV = 1.602e-19;
  double me = 9.11e-31;
  double Te = 30.0*eV;
  double logLambdaElc, nuElc;
  double den[1];
  eval_den_e(t, xn, den, ctx); 
  logLambdaElc = 6.6 - 0.5*log(den[0]/1.0e20) + 1.5*log(Te/eV);
  nuElc = logLambdaElc*pow(eV,4)*den[0]/(6.0*sqrt(2.0)*pow(M_PI,3.0/2.0)*pow(8.85e-12,2)*sqrt(me)*pow(Te,3.0/2.0));
  fout[0] = nuElc/1.96;
}
void eval_nu_ie(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double me = 9.11e-31;
  double mi = 1.67e-27;
  double nu_ei[1];
  eval_nu_ei(t, xn, nu_ei, ctx);  
  fout[0] = nu_ei[0]*me/mi;
}

void test_1x1v(int poly_order, bool use_gpu)
{
  double eV = 1.602e-19;
  double me = 9.11e-31;
  double mi = 1.67e-27;
  double Te = 30.0*eV;
  double delta_sr = 1.0;
  double betaGreenep1 = 1.0;
  double vt = sqrt(Te/me);

  double lower[] = {-0.5, -5.0*vt}, upper[] = {0.5, 5.0*vt};
  int cells[] = {2, 32};
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order==1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create moment arrays
  // (1) electron
  struct gkyl_array *den_e_ho, *upar_e_ho, *vtsq_e_ho;
  den_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  upar_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  vtsq_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_den_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_den_e, NULL);
  gkyl_proj_on_basis *proj_upar_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar_e, NULL);
  gkyl_proj_on_basis *proj_vtsq_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq_e, NULL);
  gkyl_proj_on_basis_advance(proj_den_e, 0.0, &confLocal, den_e_ho);
  gkyl_proj_on_basis_advance(proj_upar_e, 0.0, &confLocal, upar_e_ho);
  gkyl_proj_on_basis_advance(proj_vtsq_e, 0.0, &confLocal, vtsq_e_ho);
  struct gkyl_array *den_e, *upar_e, *vtsq_e;
  if (use_gpu) {
    den_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    upar_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    vtsq_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(den_e, den_e_ho);
    gkyl_array_copy(upar_e, upar_e_ho);
    gkyl_array_copy(vtsq_e, vtsq_e_ho);
  } else {
    den_e = den_e_ho;
    upar_e = upar_e_ho;
    vtsq_e = vtsq_e_ho;
  }
  struct gkyl_array *prim_moms_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(prim_moms_e, 1., den_e, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_e, 1., upar_e, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_e, 1., vtsq_e, 2*confBasis.num_basis);
  // (2) ion
  struct gkyl_array *den_i_ho, *upar_i_ho, *vtsq_i_ho;
  den_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  upar_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  vtsq_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_den_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_den_i, NULL);
  gkyl_proj_on_basis *proj_upar_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar_i, NULL);
  gkyl_proj_on_basis *proj_vtsq_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq_i, NULL);
  gkyl_proj_on_basis_advance(proj_den_i, 0.0, &confLocal, den_i_ho);
  gkyl_proj_on_basis_advance(proj_upar_i, 0.0, &confLocal, upar_i_ho);
  gkyl_proj_on_basis_advance(proj_vtsq_i, 0.0, &confLocal, vtsq_i_ho);
  struct gkyl_array *den_i, *upar_i, *vtsq_i;
  if (use_gpu) {
    den_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    upar_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    vtsq_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(den_i, den_i_ho);
    gkyl_array_copy(upar_i, upar_i_ho);
    gkyl_array_copy(vtsq_i, vtsq_i_ho);
  } else {
    den_i = den_i_ho;
    upar_i = upar_i_ho;
    vtsq_i = vtsq_i_ho;
  }
  struct gkyl_array *prim_moms_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(prim_moms_i, 1., den_i, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_i, 1., upar_i, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_i, 1., vtsq_i, 2*confBasis.num_basis);
  
  // Create collisionality arrays
  struct gkyl_array *nu_ei = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_proj_on_basis *proj_nu_ei = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_nu_ei, NULL);
  struct gkyl_array *nu_ie = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_proj_on_basis *proj_nu_ie = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_nu_ie, NULL);
  gkyl_proj_on_basis_advance(proj_nu_ei, 0.0, &confLocal, nu_ei);
  gkyl_proj_on_basis_advance(proj_nu_ie, 0.0, &confLocal, nu_ie);

  // Calculate the cross primitive moments
  gkyl_gyrokinetic_cross_prim_moms_bgk *crossPrimMomsCalc = gkyl_gyrokinetic_cross_prim_moms_bgk_new(&basis, &confBasis, use_gpu);
  struct gkyl_array *prim_moms_cross_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  struct gkyl_array *prim_moms_cross_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance(crossPrimMomsCalc, &confLocal, delta_sr, betaGreenep1,
     me, prim_moms_e, mi, prim_moms_i, prim_moms_cross_e);
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance(crossPrimMomsCalc, &confLocal, delta_sr, betaGreenep1,
     mi, prim_moms_i, me, prim_moms_e, prim_moms_cross_i);
  gkyl_gyrokinetic_cross_prim_moms_bgk_release(crossPrimMomsCalc);

  // Write out on host
  //char fname[1024];
  //sprintf(fname, "ctest_gyrokinetic_cross_prim_moms_bgk_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  //gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, prim_moms_cross, fname);

  // Compare with the expected cross moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *primMomsCross_e = gkyl_array_cfetch(prim_moms_cross_e, linidx);
    const double *primMomsCross_i = gkyl_array_cfetch(prim_moms_cross_i, linidx);
    TEST_CHECK( gkyl_compare(1.0e19, primMomsCross_e[0*confBasis.num_basis]/sqrt(2), 1e-12*1.0e19) );
    TEST_CHECK( gkyl_compare(1.16391130e4, primMomsCross_e[1*confBasis.num_basis]/sqrt(2), 1e-12*1.16391130e4) );
    TEST_CHECK( gkyl_compare(5.27398867e12, primMomsCross_e[2*confBasis.num_basis]/sqrt(2), 1e-12*5.27398867e12) );
    TEST_CHECK( gkyl_compare(1.0e19, primMomsCross_i[0*confBasis.num_basis]/sqrt(2), 1e-12*1.0e19) );
    TEST_CHECK( gkyl_compare(1.16391130e4, primMomsCross_i[1*confBasis.num_basis]/sqrt(2), 1e-12*1.16391130e4) );
    TEST_CHECK( gkyl_compare(2.74816328e9, primMomsCross_i[2*confBasis.num_basis]/sqrt(2), 1e-12*2.74816328e9) );
    TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e", primMomsCross_e[0*confBasis.num_basis]/sqrt(2), primMomsCross_e[1*confBasis.num_basis]/sqrt(2), primMomsCross_e[2*confBasis.num_basis]/sqrt(2));
  } // The hard coded numbers are expected values. The basis is looked up in maxima with: load("basis-precalc/basisSer1x"); polyOrder:1$ basis:basisC[polyOrder];

  // Release memory for moment data object
  gkyl_array_release(den_e_ho);
  gkyl_array_release(upar_e_ho);
  gkyl_array_release(vtsq_e_ho);
  gkyl_array_release(den_i_ho);
  gkyl_array_release(upar_i_ho);
  gkyl_array_release(vtsq_i_ho);
  if (use_gpu) {
    gkyl_array_release(den_e);
    gkyl_array_release(upar_e);
    gkyl_array_release(vtsq_e);
    gkyl_array_release(den_i);
    gkyl_array_release(upar_i);
    gkyl_array_release(vtsq_i);
  }  
  gkyl_array_release(prim_moms_e);
  gkyl_array_release(prim_moms_i);
  gkyl_array_release(prim_moms_cross_e);
  gkyl_array_release(prim_moms_cross_i);
  gkyl_array_release(nu_ei);
  gkyl_array_release(nu_ie);
  gkyl_proj_on_basis_release(proj_nu_ei);
  gkyl_proj_on_basis_release(proj_nu_ie);
  gkyl_proj_on_basis_release(proj_den_e);
  gkyl_proj_on_basis_release(proj_upar_e);
  gkyl_proj_on_basis_release(proj_vtsq_e);
  gkyl_proj_on_basis_release(proj_den_i);
  gkyl_proj_on_basis_release(proj_upar_i);
  gkyl_proj_on_basis_release(proj_vtsq_i);
}

void test_1x2v(int poly_order, bool use_gpu)
{
  double eV = 1.602e-19;
  double me = 9.11e-31;
  double mi = 1.67e-27;
  double Te = 30.0*eV;
  double delta_sr = 1.0;
  double betaGreenep1 = 1.0;
  double vt = sqrt(Te/me);

  double lower[] = {-0.5, -5.0*vt, 0.0}, upper[] = {0.5, 5.0*vt, 5.0*vt};
  int cells[] = {2, 32, 32};
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order==1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create moment arrays
  // (1) electron
  struct gkyl_array *den_e_ho, *upar_e_ho, *vtsq_e_ho;
  den_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  upar_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  vtsq_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_den_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_den_e, NULL);
  gkyl_proj_on_basis *proj_upar_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar_e, NULL);
  gkyl_proj_on_basis *proj_vtsq_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq_e, NULL);
  gkyl_proj_on_basis_advance(proj_den_e, 0.0, &confLocal, den_e_ho);
  gkyl_proj_on_basis_advance(proj_upar_e, 0.0, &confLocal, upar_e_ho);
  gkyl_proj_on_basis_advance(proj_vtsq_e, 0.0, &confLocal, vtsq_e_ho);
  struct gkyl_array *den_e, *upar_e, *vtsq_e;
  if (use_gpu) {
    den_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    upar_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    vtsq_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(den_e, den_e_ho);
    gkyl_array_copy(upar_e, upar_e_ho);
    gkyl_array_copy(vtsq_e, vtsq_e_ho);
  } else {
    den_e = den_e_ho;
    upar_e = upar_e_ho;
    vtsq_e = vtsq_e_ho;
  }
  struct gkyl_array *prim_moms_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(prim_moms_e, 1., den_e, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_e, 1., upar_e, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_e, 1., vtsq_e, 2*confBasis.num_basis);
  // (2) ion
  struct gkyl_array *den_i_ho, *upar_i_ho, *vtsq_i_ho;
  den_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  upar_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  vtsq_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_den_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_den_i, NULL);
  gkyl_proj_on_basis *proj_upar_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar_i, NULL);
  gkyl_proj_on_basis *proj_vtsq_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq_i, NULL);
  gkyl_proj_on_basis_advance(proj_den_i, 0.0, &confLocal, den_i_ho);
  gkyl_proj_on_basis_advance(proj_upar_i, 0.0, &confLocal, upar_i_ho);
  gkyl_proj_on_basis_advance(proj_vtsq_i, 0.0, &confLocal, vtsq_i_ho);
  struct gkyl_array *den_i, *upar_i, *vtsq_i;
  if (use_gpu) {
    den_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    upar_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    vtsq_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(den_i, den_i_ho);
    gkyl_array_copy(upar_i, upar_i_ho);
    gkyl_array_copy(vtsq_i, vtsq_i_ho);
  } else {
    den_i = den_i_ho;
    upar_i = upar_i_ho;
    vtsq_i = vtsq_i_ho;
  }
  struct gkyl_array *prim_moms_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(prim_moms_i, 1., den_i, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_i, 1., upar_i, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_i, 1., vtsq_i, 2*confBasis.num_basis);
  
  // Create collisionality arrays
  struct gkyl_array *nu_ei = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_proj_on_basis *proj_nu_ei = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_nu_ei, NULL);
  struct gkyl_array *nu_ie = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_proj_on_basis *proj_nu_ie = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_nu_ie, NULL);
  gkyl_proj_on_basis_advance(proj_nu_ei, 0.0, &confLocal, nu_ei);
  gkyl_proj_on_basis_advance(proj_nu_ie, 0.0, &confLocal, nu_ie);

  // Calculate the cross primitive moments
  gkyl_gyrokinetic_cross_prim_moms_bgk *crossPrimMomsCalc = gkyl_gyrokinetic_cross_prim_moms_bgk_new(&basis, &confBasis, use_gpu);
  struct gkyl_array *prim_moms_cross_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  struct gkyl_array *prim_moms_cross_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance(crossPrimMomsCalc, &confLocal, delta_sr, betaGreenep1,
    me, prim_moms_e, mi, prim_moms_i, prim_moms_cross_e);
  gkyl_gyrokinetic_cross_prim_moms_bgk_advance(crossPrimMomsCalc, &confLocal, delta_sr, betaGreenep1,
    mi, prim_moms_i, me, prim_moms_e, prim_moms_cross_i);
  gkyl_gyrokinetic_cross_prim_moms_bgk_release(crossPrimMomsCalc);

  // Write out on host
  //char fname[1024];
  //sprintf(fname, "ctest_gyrokinetic_cross_prim_moms_bgk_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  //gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, prim_moms_cross, fname);

  // Compare with the expected cross moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *primMomsCross_e = gkyl_array_cfetch(prim_moms_cross_e, linidx);
    const double *primMomsCross_i = gkyl_array_cfetch(prim_moms_cross_i, linidx);
    TEST_CHECK( gkyl_compare(1.0e19       , primMomsCross_e[0*confBasis.num_basis]/sqrt(2), 1e-12*1.0e19) );
    TEST_CHECK( gkyl_compare(1.16391130e4 , primMomsCross_e[1*confBasis.num_basis]/sqrt(2), 1e-12*1.16391130e4) );
    TEST_CHECK( gkyl_compare(5.27373215e12, primMomsCross_e[2*confBasis.num_basis]/sqrt(2), 1e-12*5.27373215e12) );
    TEST_CHECK( gkyl_compare(1.0e19       , primMomsCross_i[0*confBasis.num_basis]/sqrt(2), 1e-12*1.0e19) );
    TEST_CHECK( gkyl_compare(1.16391130e4 , primMomsCross_i[1*confBasis.num_basis]/sqrt(2), 1e-12*1.16391130e4) );
    TEST_CHECK( gkyl_compare(2.86262992e9 , primMomsCross_i[2*confBasis.num_basis]/sqrt(2), 1e-12*2.86262992e9) );
    TEST_MSG( "Expeced: %.9e | Got: %.9e\n", 2.86262992e9 , primMomsCross_i[2*confBasis.num_basis]/sqrt(2));
  } // The hard coded numbers are expected values. The basis is looked up in maxima with: load("basis-precalc/basisSer1x"); polyOrder:1$ basis:basisC[polyOrder];

  // Release memory for moment data object
  gkyl_array_release(den_e_ho);
  gkyl_array_release(upar_e_ho);
  gkyl_array_release(vtsq_e_ho);
  gkyl_array_release(den_i_ho);
  gkyl_array_release(upar_i_ho);
  gkyl_array_release(vtsq_i_ho);
  if (use_gpu) {
    gkyl_array_release(den_e);
    gkyl_array_release(upar_e);
    gkyl_array_release(vtsq_e);
    gkyl_array_release(den_i);
    gkyl_array_release(upar_i);
    gkyl_array_release(vtsq_i);
  }  
  gkyl_array_release(prim_moms_e);
  gkyl_array_release(prim_moms_i);
  gkyl_array_release(prim_moms_cross_e);
  gkyl_array_release(prim_moms_cross_i);
  gkyl_array_release(nu_ei);
  gkyl_array_release(nu_ie);
  gkyl_proj_on_basis_release(proj_nu_ei);
  gkyl_proj_on_basis_release(proj_nu_ie);
  gkyl_proj_on_basis_release(proj_den_e);
  gkyl_proj_on_basis_release(proj_upar_e);
  gkyl_proj_on_basis_release(proj_vtsq_e);
  gkyl_proj_on_basis_release(proj_den_i);
  gkyl_proj_on_basis_release(proj_upar_i);
  gkyl_proj_on_basis_release(proj_vtsq_i);
}

// Test the gyrokinetic BGK cross primitive-moment calculator in the special
// case of two identical species.  The cross primitive moments (n, upar, vt^2)
// of a species relaxing against an identical species must be identical to its
// own primitive moments: the density is unchanged, and the cross drift speed
// and thermal speed equal the (common) input values.  Uses spatially-constant
// fields so the expected DG coefficients are exact.
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
  struct gkyl_array *prim_moms = mkarr(3*nb, confLocal_ext.volume, false);
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

  struct gkyl_array *out = mkarr(3*nb, confLocal_ext.volume, false);
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

void test_cross_prim_moms_bgk_1x1v_p1_ho() {test_1x1v(1, false);}
void test_cross_prim_moms_bgk_1x2v_p1_ho() {test_1x2v(1, false);}

TEST_LIST = {
  {"test_cross_prim_moms_bgk_1x1v_p1_ho", test_cross_prim_moms_bgk_1x1v_p1_ho},
  {"test_cross_prim_moms_bgk_1x2v_p1_ho", test_cross_prim_moms_bgk_1x2v_p1_ho},
  { "equal_1x1v_p1", test_equal_1x1v_p1 },
  { "equal_1x2v_p1", test_equal_1x2v_p1 },
  { "equal_1x1v_zeroflow_p1", test_equal_1x1v_zeroflow_p1 },
  {NULL, NULL},
};
