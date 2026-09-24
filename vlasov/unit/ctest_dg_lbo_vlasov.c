// A test for the vlasov LBO collision updater
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_lbo_vlasov_drag.h>
#include <gkyl_dg_lbo_vlasov_diff.h>
#include <gkyl_dg_updater_lbo_vlasov.h>
#include <gkyl_array_rio.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

void nu_prof(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double vy  = xn[2];
  fout[0] = 1.0;
}

void maxwellian1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double vy  = xn[2];
  fout[0] = 1.0/(2*M_PI)*exp(-(pow(vx, 2) + pow(vy, 2))/2);
}

void maxwellian1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  fout[0] = 1.0/sqrt(2*M_PI)*exp(-(pow(vx, 2))/2);
}

void
test_1x1v_p2()
{
  // initialize grid and ranges
  int cdim = 1, vdim = 1;
  int pdim = cdim+vdim;

  int cells[] = {2, 4};
  int ghost[] = {0, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  gkyl_proj_on_basis *projF = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, maxwellian1x1v, NULL);
  gkyl_proj_on_basis *projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, nu_prof, NULL);

  struct gkyl_array *cflrate, *rhs, *fin, *nuSum, *nuUSum, *nuVtSqSum, *nuPrimMomsSum;
  cflrate = mkarr(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  nuSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  nuUSum = mkarr(vdim*confBasis.num_basis, confRange_ext.volume);
  nuVtSqSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  nuPrimMomsSum = mkarr((vdim+1)*confBasis.num_basis, confRange_ext.volume);

  gkyl_proj_on_basis_advance(projF, 0.0, &phaseRange_ext, fin);
  gkyl_proj_on_basis_advance(projNu, 0.0, &confRange_ext, nuSum);
  gkyl_proj_on_basis_release(projF);
  gkyl_proj_on_basis_release(projNu);
  
  gkyl_array_clear(nuUSum, 0.0);
  gkyl_array_clear(nuVtSqSum, 1.0);
  gkyl_array_set_offset(nuPrimMomsSum, 1.0, nuUSum, 0);
  gkyl_array_set_offset(nuPrimMomsSum, 1.0, nuVtSqSum, vdim*confBasis.num_basis);

  // initialize hyper_dg slvr
  int up_dirs[] = {1};
  int zero_flux_flags[] = {1};

  gkyl_dg_updater_collisions *slvr;
  // LBO updater
  struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = { .nuSum = nuSum, .nuPrimMomsSum = nuPrimMomsSum };
  struct gkyl_dg_lbo_vlasov_diff_auxfields diff_inp = { .nuSum = nuSum, .nuPrimMomsSum = nuPrimMomsSum };
  slvr = gkyl_dg_updater_lbo_vlasov_new(&phaseGrid, &confBasis, &basis, &confRange, &drag_inp, &diff_inp, false);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_updater_lbo_vlasov_advance(slvr, &phaseRange, fin, cflrate, rhs);
  }

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d1, *rhs_d2;
  int i = 0;
  while(val==0) {
    rhs_d1 = gkyl_array_fetch(rhs, i);
    val = rhs_d1[0];
    if(val==0) i++;
    }
  TEST_CHECK(i == linl);

  // get linear index of some other cell
  int idx1[] = {1, 1};
  int idx2[] = {2, 3};
  int linl1 = gkyl_range_idx(&phaseRange, idx1);
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d1 = gkyl_array_fetch(rhs, linl1);
  rhs_d2 = gkyl_array_fetch(rhs, linl2);
  TEST_CHECK( gkyl_compare_double(rhs_d1[0], -0.20752566554314, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[1], 0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[2], 0.52331159647886, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[3], -1.3061806047047, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[4], 0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[5], -0.59363481492976, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[6], -1.3061806047047, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[7], 1.5729635651345, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[0], 0.20752566554313, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[1], -0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[2], -0.021201659682202, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[3], 0.05635142305445, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[4], -0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[5], -0.00090746168862665, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[6], 0.056351423054451, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[7], 0.015080609728081, 1e-12) );

  // release memory for moment data object
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(fin);
  gkyl_array_release(nuSum);
  gkyl_array_release(nuUSum);
  gkyl_array_release(nuVtSqSum);
  gkyl_array_release(nuPrimMomsSum);
  gkyl_dg_updater_lbo_vlasov_release(slvr);
}

void
test_1x2v_p2()
{
  // initialize grid and ranges
  int cdim = 1, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {24, 12, 12};
  int ghost[] = {0, 0, 0};
  double lower[] = {0., -1., -1.};
  double upper[] = {1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  gkyl_proj_on_basis *projF = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, maxwellian1x2v, NULL);
  gkyl_proj_on_basis *projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, nu_prof, NULL);

  struct gkyl_array *cflrate, *rhs, *fin, *nuSum, *nuUSum, *nuVtSqSum, *nuPrimMomsSum;
  cflrate = mkarr(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  nuSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  nuUSum = mkarr(vdim*confBasis.num_basis, confRange_ext.volume);
  nuVtSqSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  nuPrimMomsSum = mkarr((vdim+1)*confBasis.num_basis, confRange_ext.volume);

  gkyl_proj_on_basis_advance(projF, 0.0, &phaseRange_ext, fin);
  gkyl_proj_on_basis_advance(projNu, 0.0, &confRange_ext, nuSum);
  gkyl_proj_on_basis_release(projF);
  gkyl_proj_on_basis_release(projNu);
  
  gkyl_array_clear(nuUSum, 0.0);
  gkyl_array_clear(nuVtSqSum, 1.0);
  gkyl_array_set_offset(nuPrimMomsSum, 1.0, nuUSum, 0);
  gkyl_array_set_offset(nuPrimMomsSum, 1.0, nuVtSqSum, vdim*confBasis.num_basis);

  // initialize hyper_dg slvr
  int up_dirs[] = {1, 2};
  int zero_flux_flags[] = {1, 1};

  gkyl_dg_updater_collisions *slvr;
  // LBO updater
  struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = { .nuSum = nuSum, .nuPrimMomsSum = nuPrimMomsSum };
  struct gkyl_dg_lbo_vlasov_diff_auxfields diff_inp = { .nuSum = nuSum, .nuPrimMomsSum = nuPrimMomsSum };
  slvr = gkyl_dg_updater_lbo_vlasov_new(&phaseGrid, &confBasis, &basis, &confRange, &drag_inp, &diff_inp, false);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    
    gkyl_dg_updater_lbo_vlasov_advance(slvr, &phaseRange, fin, cflrate, rhs);
  }

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);
  
  // check that ghost cells are empty
  double val = 0;
  double *rhs_d1, *rhs_d2;
  int i = 0;
  while(val==0) {
    rhs_d1 = gkyl_array_fetch(rhs, i);
    val = rhs_d1[0];
    if(val==0) i++;
    }
  TEST_CHECK(i == linl);

  // get linear index of some other cell
  int idx1[] = {1, 1, 1};
  int idx2[] = {3, 3, 3};
  int linl1 = gkyl_range_idx(&phaseRange, idx1);
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d1 = gkyl_array_fetch(rhs, linl1);
  rhs_d2 = gkyl_array_fetch(rhs, linl2);
  TEST_CHECK( gkyl_compare_double(rhs_d1[0], -0.61212875932713, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[1], 1.4774325329977, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[2], 0.53637935257539, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[3], 0.53637935257539, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[4], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[5], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[6], 0.048437402462207, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[7], 1.4774325329977, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[8], -0.6997889965392, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[9], -0.69978899653921, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[10], -0.11735479971929, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[11], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[12], 1.7060158928832, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[13], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[14], -0.030917104954987, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[15], 1.7060158928832, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[16], -0.030917104954985, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[17], -0.11735479971929, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[18], 0.075366767497302, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[19], 0.075366767497302, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[0], 0.12351087529995, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[1], -0.29806419452944, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[2], 0.008464559638675, 1e-11) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[3], 0.0084645596386717, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[4], -0.021026986912361, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[5], -0.02102698691237, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[6], 0.00037753552357966, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[7], -0.29806419452945, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[8], -0.000056168510255361, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[9], -0.000056168510261349, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[10], -0.00094472239586452, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[11], -0.02102698691235, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[12], 0.00038921843594153, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[13], -0.021026986912361, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[14], -0.000010512115344419, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[15], 0.00038921843590858, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[16], -0.000010512115345234, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[17], -0.00094472239587751, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[18], 0.000010945446315969, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[19], 0.000010945446318925, 1e-8) );

  // release memory for moment data object
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(fin);
  gkyl_array_release(nuSum);
  gkyl_array_release(nuUSum);
  gkyl_array_release(nuVtSqSum);
  gkyl_array_release(nuPrimMomsSum);
  gkyl_dg_updater_lbo_vlasov_release(slvr);
}


#ifdef GKYL_HAVE_CUDA

void
test_1x1v_p2_cu()
{
  // initialize grid and ranges
  int cdim = 1, vdim = 1;
  int pdim = cdim+vdim;

  int cells[] = {2, 4};
  int ghost[] = {0, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  gkyl_proj_on_basis *projF = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, maxwellian1x1v, NULL);
  gkyl_proj_on_basis *projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, nu_prof, NULL);

  struct gkyl_array *cflrate_cu, *rhs, *rhs_cu, *fin, *fin_cu, *nuSum, *nuSum_cu,
                    *nuUSum_cu, *nuVtSqSum_cu, *nuPrimMomsSum_cu;
  cflrate_cu = mkarr_cu(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  rhs_cu = mkarr_cu(basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  fin_cu = mkarr_cu(basis.num_basis, phaseRange_ext.volume);
  nuSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  nuSum_cu = mkarr_cu(confBasis.num_basis, confRange_ext.volume);
  nuUSum_cu = mkarr_cu(vdim*confBasis.num_basis, confRange_ext.volume);
  nuVtSqSum_cu = mkarr_cu(confBasis.num_basis, confRange_ext.volume);
  nuPrimMomsSum_cu = mkarr_cu((vdim+1)*confBasis.num_basis, confRange_ext.volume);

  gkyl_proj_on_basis_advance(projF, 0.0, &phaseRange_ext, fin);
  gkyl_proj_on_basis_advance(projNu, 0.0, &confRange_ext, nuSum);
  gkyl_proj_on_basis_release(projF);
  gkyl_proj_on_basis_release(projNu);
  
  gkyl_array_clear(nuUSum_cu, 0.0);
  gkyl_array_clear(nuVtSqSum_cu, 1.0);
  gkyl_array_set_offset(nuPrimMomsSum_cu, 1.0, nuUSum_cu, 0);
  gkyl_array_set_offset(nuPrimMomsSum_cu, 1.0, nuVtSqSum_cu, vdim*confBasis.num_basis);

  gkyl_array_copy(fin_cu, fin);
  gkyl_array_copy(nuSum_cu, nuSum);

  // initialize hyper_dg slvr
  int up_dirs[] = {1};
  int zero_flux_flags[] = {1};

  gkyl_dg_updater_collisions *slvr;
  // LBO updater
  struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = { .nuSum = nuSum_cu, .nuPrimMomsSum = nuPrimMomsSum_cu };
  struct gkyl_dg_lbo_vlasov_diff_auxfields diff_inp = { .nuSum = nuSum_cu, .nuPrimMomsSum = nuPrimMomsSum_cu };
  slvr = gkyl_dg_updater_lbo_vlasov_new(&phaseGrid, &confBasis, &basis, &confRange, &drag_inp, &diff_inp, true);
  
  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs_cu, 0.0);
    gkyl_array_clear(cflrate_cu, 0.0);
    gkyl_dg_updater_lbo_vlasov_advance(slvr, &phaseRange,
      fin_cu, cflrate_cu, rhs_cu);
  }
  gkyl_array_copy(rhs, rhs_cu);

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d1, *rhs_d2;
  int i = 0;
  while(val==0) {
    rhs_d1 = gkyl_array_fetch(rhs, i);
    val = rhs_d1[0];
    if(val==0) i++;
    }
  TEST_CHECK(i == linl);

  // get linear index of some other cell
  int idx1[] = {1, 1};
  int idx2[] = {2, 3};
  int linl1 = gkyl_range_idx(&phaseRange, idx1);
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d1 = gkyl_array_fetch(rhs, linl1);
  rhs_d2 = gkyl_array_fetch(rhs, linl2);
  TEST_CHECK( gkyl_compare_double(rhs_d1[0], -0.20752566554314, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[1], 0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[2], 0.52331159647886, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[3], -1.3061806047047, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[4], 0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[5], -0.59363481492976, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[6], -1.3061806047047, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[7], 1.5729635651345, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[0], 0.20752566554313, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[1], -0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[2], -0.021201659682202, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[3], 0.05635142305445, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[4], -0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[5], -0.00090746168862665, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[6], 0.056351423054451, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[7], 0.015080609728081, 1e-12) );

  // release memory for moment data object
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_cu);
  gkyl_array_release(cflrate_cu);
  gkyl_array_release(fin);
  gkyl_array_release(fin_cu);
  gkyl_array_release(nuSum);
  gkyl_array_release(nuSum_cu);
  gkyl_array_release(nuUSum_cu);
  gkyl_array_release(nuVtSqSum_cu);
  gkyl_array_release(nuPrimMomsSum_cu);
  gkyl_dg_updater_lbo_vlasov_release(slvr);
}

void
test_1x2v_p2_cu()
{
  // initialize grid and ranges
  int cdim = 1, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {24, 12, 12};
  int ghost[] = {0, 0, 0};
  double lower[] = {0., -1., -1.};
  double upper[] = {1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  gkyl_proj_on_basis *projF = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, maxwellian1x2v, NULL);
  gkyl_proj_on_basis *projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, nu_prof, NULL);

  struct gkyl_array *cflrate_cu, *rhs, *rhs_cu, *fin, *nuSum, *fin_cu,
                    *nuSum_cu, *nuUSum_cu, *nuVtSqSum_cu, *nuPrimMomsSum_cu;
  cflrate_cu = mkarr_cu(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  rhs_cu = mkarr_cu(basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  nuSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  fin_cu = mkarr_cu(basis.num_basis, phaseRange_ext.volume);
  nuSum_cu = mkarr_cu(confBasis.num_basis, confRange_ext.volume);
  nuUSum_cu = mkarr_cu(vdim*confBasis.num_basis, confRange_ext.volume);
  nuVtSqSum_cu = mkarr_cu(confBasis.num_basis, confRange_ext.volume);
  nuPrimMomsSum_cu = mkarr_cu((vdim+1)*confBasis.num_basis, confRange_ext.volume);

  gkyl_proj_on_basis_advance(projF, 0.0, &phaseRange_ext, fin);
  gkyl_proj_on_basis_advance(projNu, 0.0, &confRange_ext, nuSum);
  gkyl_proj_on_basis_release(projF);
  gkyl_proj_on_basis_release(projNu);
  
  gkyl_array_copy(fin_cu, fin);
  gkyl_array_copy(nuSum_cu, nuSum);
  
  gkyl_array_clear(nuUSum_cu, 0.0);
  gkyl_array_clear(nuVtSqSum_cu, 1.0);
  gkyl_array_set_offset(nuPrimMomsSum_cu, 1.0, nuUSum_cu, 0);
  gkyl_array_set_offset(nuPrimMomsSum_cu, 1.0, nuVtSqSum_cu, vdim*confBasis.num_basis);

  // initialize hyper_dg slvr
  int up_dirs[] = {1, 2};
  int zero_flux_flags[] = {1, 1};

  gkyl_dg_updater_collisions *slvr;
  // LBO updater
  struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = { .nuSum = nuSum_cu, .nuPrimMomsSum = nuPrimMomsSum_cu };
  struct gkyl_dg_lbo_vlasov_diff_auxfields diff_inp = { .nuSum = nuSum_cu, .nuPrimMomsSum = nuPrimMomsSum_cu };
  slvr = gkyl_dg_updater_lbo_vlasov_new(&phaseGrid, &confBasis, &basis, &confRange, &drag_inp, &diff_inp, true);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs_cu, 0.0);
    gkyl_array_clear(cflrate_cu, 0.0);
    
    gkyl_dg_updater_lbo_vlasov_advance(slvr, &phaseRange,
      fin_cu, cflrate_cu, rhs_cu);
  }
  gkyl_array_copy(rhs, rhs_cu);
    
  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d1, *rhs_d2;
  int i = 0;
  while(val==0) {
    rhs_d1 = gkyl_array_fetch(rhs, i);
    val = rhs_d1[0];
    if(val==0) i++;
    }
  TEST_CHECK(i == linl);
  
  // get linear index of some other cell
  int idx1[] = {1, 1, 1};
  int idx2[] = {3, 3, 3};
  int linl1 = gkyl_range_idx(&phaseRange, idx1);
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d1 = gkyl_array_fetch(rhs, linl1);
  rhs_d2 = gkyl_array_fetch(rhs, linl2);
  TEST_CHECK( gkyl_compare_double(rhs_d1[0], -0.61212875932713, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[1], 1.4774325329977, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[2], 0.53637935257539, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[3], 0.53637935257539, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[4], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[5], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[6], 0.048437402462207, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[7], 1.4774325329977, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[8], -0.6997889965392, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[9], -0.69978899653921, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[10], -0.11735479971929, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[11], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[12], 1.7060158928832, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[13], -1.2996714533493, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[14], -0.030917104954987, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[15], 1.7060158928832, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[16], -0.030917104954985, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[17], -0.11735479971929, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[18], 0.075366767497302, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[19], 0.075366767497302, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[0], 0.12351087529995, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[1], -0.29806419452944, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[2], 0.008464559638675, 1e-11) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[3], 0.0084645596386717, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[4], -0.021026986912361, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[5], -0.02102698691237, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[6], 0.00037753552357966, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[7], -0.29806419452945, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[8], -0.000056168510255361, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[9], -0.000056168510261349, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[10], -0.00094472239586452, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[11], -0.02102698691235, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[12], 0.00038921843594153, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[13], -0.021026986912361, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[14], -0.000010512115344419, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[15], 0.00038921843590858, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[16], -0.000010512115345234, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[17], -0.00094472239587751, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[18], 0.000010945446315969, 1e-8) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[19], 0.000010945446318925, 1e-8) );

  // release memory for moment data object
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_cu);
  gkyl_array_release(cflrate_cu);
  gkyl_array_release(fin);
  gkyl_array_release(nuSum);
  gkyl_array_release(fin_cu);
  gkyl_array_release(nuSum_cu);
  gkyl_array_release(nuUSum_cu);
  gkyl_array_release(nuVtSqSum_cu);
  gkyl_array_release(nuPrimMomsSum_cu);
  gkyl_dg_updater_lbo_vlasov_release(slvr);
}
#endif

TEST_LIST = {
  { "test_1x1v_p2", test_1x1v_p2 },
  { "test_1x2v_p2", test_1x2v_p2 },
  #ifdef GKYL_HAVE_CUDA
  { "test_1x1v_p2_cu", test_1x1v_p2_cu },
  { "test_1x2v_p2_cu", test_1x2v_p2_cu },
  #endif
  { NULL, NULL },
};
