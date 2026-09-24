#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void 
eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void 
eval_M2(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0; 
}

void 
eval_M1i_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
  fout[1] = 0.5;
}

void 
eval_M1i_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.5; 
  fout[1] = 1.0; 
  fout[2] = 0.5;
}

void
eval_M0_x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.1 + sin(xn[0]);
}

void 
eval_M1i_1v_x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5*(sin(xn[0])); // 0.5;
}

void 
eval_M2_x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T*(1.1 + sin(xn[0]));
}

void 
test_1x1v(int poly_order)
{
  double lower[] = {0.1, -10.0}, upper[] = {1.0, 10.0};
  int cells[] = {2, 32}; 
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
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

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Correct the distribution function
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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  gkyl_array_set_offset_range(m0, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {2.7845923966306263e-01, 4.6204926597763572e-17, 
    6.5920727847853647e-02, -1.4768663476574924e-19, -3.3457863588588799e-17, 
    2.5460491687342435e-03, -4.0874057255360184e-17, -1.0795797155980702e-17};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 16}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
      // printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}


void 
test_1x1v_spatially_varied(int poly_order)
{
  double lower[] = {0.0, -10.0}, upper[] = {10.0, 10.0};
  int cells[] = {10, 32}; // 1001
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
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

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0_x, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v_x, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_x, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Correct the distribution function
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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  gkyl_array_set_offset_range(m0, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  // Write the output (moments)
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_n_corr.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m0_corr,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_vb_corr.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m1i_corr,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_T_corr.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m2_corr,fname);


  // Write the output (moments)
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_n.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m0,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_vb.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m1i,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_T.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m2,fname);


  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double m0_vals[] = {2.2057459659974716e+00,3.4938172762332936e-01,-2.4953670224477906e-02};
  double m1i_vals[] = {3.2505552369353341e-01,1.7469086381166479e-01,-1.2476835112239079e-02};
  double m2_vals[] = {2.2057459659974712e+00,3.4938172762332981e-01,-2.4953670224477753e-02};


  //const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 16}));
  const double *m0_fixed = gkyl_array_cfetch(m0, gkyl_range_idx(&confLocal_ext, (int[1]){1}));
  const double *m1i_fixed = gkyl_array_cfetch(m1i, gkyl_range_idx(&confLocal_ext, (int[1]){1}));
  const double *m2_fixed = gkyl_array_cfetch(m2, gkyl_range_idx(&confLocal_ext, (int[1]){1}));

  if (poly_order == 2)
    for (int i = 0; i < confBasis.num_basis; ++i){
      TEST_CHECK(gkyl_compare_double(m0_vals[i], m0_fixed[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(m1i_vals[i], m1i_fixed[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(m2_vals[i], m2_fixed[i], 1e-12));
    }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}


void 
test_1x2v(int poly_order)
{
  double lower[] = {0.1, -10.0, -10.0}, upper[] = {1.0, 10.0, 10.0};
  int cells[] = {2, 32, 32};
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Correct the distribution function
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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  gkyl_array_set_offset_range(m0, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {7.1154795741657104e-02, -1.0154349268741242e-17, 1.7043680613403746e-02, 
    1.0985977050895525e-02, -7.6488419515851888e-18, -7.8668793275235743e-18, 2.8173732360475563e-03, 
    -1.2268169102706324e-17, 6.6921571942179570e-04, -4.4309884559864081e-04, -4.2230625985903307e-18, 
    -1.0393016114178111e-17, -5.8390877219590108e-19, -1.0792789352242940e-17, 1.0442989813709548e-04, 
    5.9309279220864833e-19, -1.1960215598429369e-04, -2.6085561563408531e-18, 2.0845874350775922e-18, 
    -3.4160776312939708e-19};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]){1, 16, 16}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
      // printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void 
test_1x3v(int poly_order)
{
  double lower[] = {0.1, -10.0, -10.0, -10.0}, upper[] = {1.0, 10.0, 10.0, 10.0};
  int cells[] = {2, 16, 16, 16};
  int vdim = 3, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Correct the distribution function
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
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  gkyl_array_set_offset_range(m0, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x3v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {2.3681169627173117e-03, -1.4446897638564469e-19, 1.7212114821089531e-03, 
    1.4278819468396819e-03, 1.1009127860830206e-03, -1.3416643175438523e-20, -9.1242971227460240e-21, 
    1.0665480449733901e-03, 9.1027098303971133e-20, 8.3350273031015259e-04, 7.0324247191191436e-04, 
    -4.4239730499907700e-19, 4.8462190085462665e-04, 2.7916039109363608e-04, 9.5955866471644934e-05, 
    2.8545136484013421e-20, 3.5931632283878487e-21, 7.0377414655311627e-20, 5.4868903058922895e-04, 
    -2.4033547244859593e-19, -4.3147003734719234e-20, -1.1369323846038000e-19, 3.1070439565220937e-04, 
    -3.6244873179843119e-21, 2.1725375592137358e-04, -2.5466173336449186e-19, 2.4678946252732544e-04, 
    1.4954398636806508e-04, 1.3237463325222323e-20, 8.1278508031050744e-05, 7.1527130914181915e-05, 
    -1.6410396660305442e-20, -6.0143678114279295e-20, -2.9118247539011877e-20, -2.7750766817365425e-20, 
    -1.0216466400671430e-19, -9.9249885235718083e-21, -1.7576662987358121e-19, 1.6898546616575218e-04, 
    9.7560025734548306e-21, 1.2228452029340925e-04, -1.1579931030235590e-20, -1.1437237301704353e-20, 
    6.2517215908454443e-05, -3.7022575562915509e-19, -2.1532599752429503e-20, -2.0849216060814057e-20, 
    1.6541329748965388e-20};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[4]){1, 8, 8, 8}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
      // printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }
void test_1x1v_p2_spatially_varied() { test_1x1v_spatially_varied(2); }
void test_1x2v_p2() { test_1x2v(2); }
void test_1x3v_p2() { test_1x3v(2); }

TEST_LIST = {
  {"test_1x1v_p2", test_1x1v_p2},
  {"test_1x1v_p2_spatially_varied", test_1x1v_p2_spatially_varied},
  {"test_1x2v_p2", test_1x2v_p2},
  {"test_1x3v_p2", test_1x3v_p2},
  {NULL, NULL},
};
