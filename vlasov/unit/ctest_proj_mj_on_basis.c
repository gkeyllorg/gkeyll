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
eval_M1i_1v_no_drift(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
}

void 
eval_M2_1v_no_drift(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; 
}

void 
eval_M2_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
  fout[1] = 0.25;
}

void 
eval_M2_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
  fout[1] = 0.25;
  fout[2] = -0.5;
}

void 
eval_M2_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
test_1x1v_no_drift(int poly_order)
{
  double lower[] = {0.1, -15.0}, upper[] = {1.0, 15.0};
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

  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

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

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v_no_drift, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_1v_no_drift, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

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
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p1_vals[] = {5.3918752026566863e-01, -1.0910243387206232e-17, -6.0196985297046972e-02,
    5.0006050167249552e-18};
  double p2_vals[] = {5.3922143701031633e-01, -9.6625223288531320e-18, -5.7898881215132203e-02,
    7.8842251929957589e-18, 1.9166441863144966e-17, -1.0173903909543560e-02,
    1.7916734900988946e-17, 1.4245174569363429e-18};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 17}));

  if (poly_order == 1)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p1_vals[i], fv[i], 1e-12));

  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x1v_p%d_no_drift.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void test_1x1v_no_drift_p2() { test_1x1v_no_drift(2); }

void 
test_1x1v(int poly_order)
{
  double lower[] = {0.1, -15.0}, upper[] = {1.0, 15.0}; // +/- 15 on velocity
  int cells[] = {2, 32};                                // default {2, 32}
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

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_1v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

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
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // test accuracy of the projection:
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
  double p2_vals[] = {5.9297594654488650e-01, -4.1867292592431142e-18, 7.1286851491369927e-03, 
    6.9143629007589357e-18, 1.0932899315657513e-17, -1.6063381083048084e-02, 
    -5.5123762524241300e-18, -3.3195546731973899e-18};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 17}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
      // printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }
  
  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }

void 
test_1x2v(int poly_order)
{
  double lower[] = {0.1, -15.0, -15.0}, upper[] = {1.0, 15.0, 15.0};
  int cells[] = {2, 16, 16};
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

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_2v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

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
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p2_vals[] = {1.6408879023240103e-01, -1.3576896184824113e-17, -9.8151809581183153e-03, 
    -2.9682559802577287e-02, 1.2829937358239629e-18, 3.8215290052990313e-19, 
    8.6596866231148772e-03, -2.9507573862413946e-18, -9.7332674540940786e-03, 
    -7.2599184095971433e-03, 6.4445077190172000e-19, 1.4972406793044089e-17, 
    1.6202061433993536e-20, 1.0653781194370923e-17, 1.7475225696783907e-03, 
    -5.7723245596542925e-19, -4.1851388419497057e-04, -4.9153245858838362e-19, 
    5.5404332026569672e-19, 5.9845689181369784e-19};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]){1, 9, 9}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
      // printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void test_1x2v_p2() { test_1x2v(2); }

void 
test_1x3v(int poly_order)
{
  double lower[] = {0.1, -15.0, -15.0, -15.0}, upper[] = {1.0, 15.0, 15.0, 15.0};
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

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_3v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

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
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // values to compare  at index (1, 9, 9, 9) [remember, lower-left index is (1,1,1,1)]
  double p2_vals[] = {2.2004825355599965e-02, 7.1324343268494727e-19, -1.2968205139278119e-03, 
    -3.9683045523914544e-03, -1.1492140196787010e-02, 2.0065506128885464e-19, 
    2.6863348725598927e-19, 1.0399887120789773e-03, 6.2237804549528725e-19, 
    1.2878059167286534e-03, 2.6451147610928889e-03, -6.3077058326192259e-18, 
    -1.3820826308538342e-03, -1.0478012443518811e-03, 1.6782192015094876e-03, 
    7.3099489316101678e-20, -2.2315380074151361e-19, -1.5210474522584246e-19, 
    -8.4301245320312257e-04, 2.8388294826914462e-18, -5.5972719240227798e-19, 
    2.7769888555314303e-18, 2.6587058568191109e-04, -2.6096614467447103e-19, 
    -2.0844243229831450e-05, 6.3431871295172754e-18, 7.3781907692053105e-04, 
    4.9347115341968450e-04, 2.2046700923304946e-20, -3.7605267297870676e-04, 
    -5.6100948433461376e-04, 3.5408403867716293e-20, -3.8294647162266746e-19, 
    -3.0565669788894797e-20, -2.3269941088952652e-20, -1.0021845103459830e-18, 
    2.4164727245648497e-19, -1.7030071204748659e-18, -1.4254568013549241e-04, 
    1.2495807446563311e-19, 3.8839940617297632e-05, -5.3617945589348144e-20, 
    -4.8911545428669850e-20, 2.5014295050500632e-04, 2.0619708529882749e-19, 
    -6.1276211233163564e-20, -6.5891739087463357e-20, -1.0199525829477060e-19};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[4]){1, 9, 9, 9}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
      // printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x3v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void test_1x3v_p2() { test_1x3v(2); }

TEST_LIST = {
  {"test_1x1v_no_drift_p2", test_1x1v_no_drift_p2},
  {"test_1x1v_p2", test_1x1v_p2},
  {"test_1x2v_p2", test_1x2v_p2},
  {"test_1x3v_p2", test_1x3v_p2},
  {NULL, NULL},
};
