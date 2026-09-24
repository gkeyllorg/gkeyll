#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_M0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.25;
}

void eval_M1i_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_M2_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double n = 1.0, vth2 = 1.0, ux = 0.5;
  double x = xn[0];
  fout[0] = n*vth2 + n*ux*ux;
}

void eval_udrift_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_vtsq_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void eval_M1i_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void eval_M2_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double n = 1.0, vth2 = 1.0, ux = 0.5, uy = 0.25;
  double x = xn[0];
  fout[0] = 2*n*vth2 + n*(ux*ux+uy*uy);
}

void eval_udrift_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void eval_vtsq_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void
test_1x1v(int poly_order, bool use_gpu)
{
  double lower[] = {0.1, -6.0}, upper[] = {1.0, 6.0};
  int cells[] = {2, 32};
  int vdim = 1, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms, *moms_diag;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume); 
  moms_diag = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume); 

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_1v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // Projection routine expects all the moments in a single array.
  gkyl_array_set_offset(moms, 1., m0, 0);
  gkyl_array_set_offset(moms, 1., m1i, confBasis.num_basis);
  gkyl_array_set_offset(moms, 1., m2, (vdim+1)*confBasis.num_basis);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .model_id = GKYL_MODEL_DEFAULT,
    .mass = 1.0,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // Compute the moments of our corrected distribution function
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .model_id = GKYL_MODEL_DEFAULT,
    .mass = 1.0,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );

  // correction updater
  struct gkyl_vlasov_lte_correct_inp inp = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .model_id = GKYL_MODEL_DEFAULT,
    .use_gpu = false,
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_lte = gkyl_vlasov_lte_correct_inew( &inp );

 // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_correct_maxwellian_test_1x1v_p%d_corr_m0.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Moments computed from *only* density-corrected LTE distribution function
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms_diag);
  struct gkyl_array *m0_n_corr_only, *m1i_n_corr_only, *m2_n_corr_only;
  m0_n_corr_only = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i_n_corr_only = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_n_corr_only = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset_range(m0_n_corr_only, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i_n_corr_only, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2_n_corr_only, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  struct gkyl_vlasov_lte_correct_status stat_corr = gkyl_vlasov_lte_correct_all_moments(corr_lte, 
    distf, moms, &local, &confLocal);

  // Moments computed from all-moment-corrected LTE distribution function 
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr;
  m0_corr = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_corr = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset_range(m0_corr, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i_corr, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2_corr, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

   // write distribution function to file
  sprintf(fname, "ctest_correct_maxwellian_test_1x1v_p%d_corr_all_moms.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);
  
  // Compare m0 to the computed m0 (density correction only)
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *n0 = gkyl_array_cfetch(m0, gkyl_range_idx(&confLocal, iter.idx));
    const double *nr = gkyl_array_cfetch(m0_n_corr_only, gkyl_range_idx(&confLocal, iter.idx));
    
    for (int k=0; k<confBasis.num_basis; ++k)
      TEST_CHECK( gkyl_compare_double(n0[k], nr[k], 1e-14) );
  }

  // Compare m0, m1i, m2 to the computed m0, m1i, m2 (all corrections)
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *n0 = gkyl_array_cfetch(m0, gkyl_range_idx(&confLocal, iter.idx));
    const double *vb0 = gkyl_array_cfetch(m1i, gkyl_range_idx(&confLocal, iter.idx));
    const double *T0 = gkyl_array_cfetch(m2, gkyl_range_idx(&confLocal, iter.idx));
    const double *nr = gkyl_array_cfetch(m0_corr, gkyl_range_idx(&confLocal, iter.idx));
    const double *vbr = gkyl_array_cfetch(m1i_corr, gkyl_range_idx(&confLocal, iter.idx));
    const double *Tr = gkyl_array_cfetch(m2_corr, gkyl_range_idx(&confLocal, iter.idx));
    
    for (int k=0; k<confBasis.num_basis; ++k){
      TEST_CHECK( gkyl_compare_double(n0[k], nr[k], 1e-14) );
      TEST_CHECK( gkyl_compare_double(vb0[k], vbr[k], 1e-14) );
      TEST_CHECK( gkyl_compare_double(T0[k], Tr[k], 1e-14) );
    }
  }
  
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(m0_n_corr_only);
  gkyl_array_release(m1i_n_corr_only);
  gkyl_array_release(m2_n_corr_only);
  gkyl_array_release(moms_diag);

  gkyl_array_release(m0); 
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);

  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_vlasov_lte_correct_release(corr_lte);
  gkyl_vlasov_lte_moments_release(lte_moms);
}

void test_1x1v_p1() { test_1x1v(1, false); }
void test_1x1v_p2() { test_1x1v(2, false); }

TEST_LIST = {
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x1v_p2", test_1x1v_p2 },
  { NULL, NULL },
};
