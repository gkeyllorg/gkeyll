  // Write m0 to file after the positivity shift.
#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_positivity_shift_vlasov.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>

#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_ctx {
  double n0; // Density.
  double ux; // Mean speed along vx.
  double uy; // Mean speed along vy.
  double temp; // Temperature.
  double mass; // Species mass.
  int vdim; // Number of velocity space dimensions.
  double vx_max; // Maximum vx of the grid.
  double vy_max; // Maximum vy of the grid.
};

void eval_distf_1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vy = xn[1], vx = xn[2];

  struct test_ctx *tctx = ctx;
  double mass = tctx->mass;
  double n0 = tctx->n0;
  double ux = tctx->ux;
  double uy = tctx->uy;
  double vtsq = tctx->temp/mass;
  int vdim = tctx->vdim;

  fout[0] = (n0/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vx-ux,2)+pow(vy-uy,2))/(2.0*vtsq));

  // Intentionally set some places to be negative.
  if (fabs(vx) > 0.8*tctx->vx_max || fabs(vy) > 0.8*tctx->vy_max)
    fout[0] = -0.2 * (n0/pow(2.0*M_PI*vtsq,vdim/2.0));
}

void
test_1x2v(int poly_order, bool use_gpu)
{
  int cdim = 1;
  double vx_max = 6.0, vy_max = 6.0;
  double lower[] = {0.0, -vx_max, -vy_max}, upper[] = {1.0, vx_max, vy_max};
  int cells[] = {2, 12, 8};

  int ndim = sizeof(cells)/sizeof(cells[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .ux = 0, // Mean speed along vx.
    .uy = 0, // Mean speed along vy.
    .temp = 2.75, // Temperature.
    .mass = 1.0, // Species mass.
    .vdim = vdim, // Number of velocity space dimensions.
    .vx_max = vx_max, // Maximum vx of the grid.
    .vy_max = vy_max, // Maximum vy of the grid.
  };

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_tensor(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = {0};
  for (int d=0; d<cdim; d++) confGhost[d] = 1;
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = {0};
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create distribution function arrays.
  struct gkyl_array *distf_ho, *distf;
  distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  distf_ho = distf;
  if (use_gpu)
    distf_ho = mkarr(false, distf->ncomp, distf->size);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x2v, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(use_gpu, velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(use_gpu, velBasis.num_basis, velLocal.volume);
  // Identity velocity map: used to build the Hamiltonian and by the moment/LTE updaters.
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velLocal, &velBasis, inp_vmap, false, use_gpu);
  gkyl_dg_vlasov_calc_hamil(&velGrid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, use_gpu);

  // Compute M0 of the original f.
  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .hamil_id = gkyl_hamil_id_from_model_id(GKYL_MODEL_DEFAULT),
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = use_gpu,
    .vel_map = vel_map,
  };


  struct gkyl_array *m0_pre = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_mom_type *m0_mom_up = gkyl_mom_vlasov_inew(&inp_mom);
  struct gkyl_mom_calc *mom_calc = gkyl_mom_calc_new(&grid, m0_mom_up, use_gpu);
  gkyl_mom_calc_advance(mom_calc, &local, &confLocal, distf, m0_pre);


//  // Write m0 to file.
//  char fname0M0[1024];
//  sprintf(fname0M0, "ctest_positivity_shift_vlasov_1x2v_p%d_m0_pre.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, m0_pre, fname0M0);

  inp_mom.mom_type = GKYL_F_MOMENT_M0M1M2;
  struct gkyl_mom_type *int_mom_up = gkyl_int_mom_vlasov_inew(&inp_mom);
  struct gkyl_mom_calc *int_mom_calc = gkyl_mom_calc_new(&grid, int_mom_up, use_gpu);

  // Compute the integrated moments of the original f.
  //struct gkyl_dg_updater_moment *int_mom_up = gkyl_dg_updater_moment_new(
  //  &grid, &confBasis, &basis, &confLocal, 0, &local, 0, 0, GKYL_F_MOMENT_M0M1M2, true, use_gpu);

  int num_mom = int_mom_up->num_mom;
  //int num_mom = gkyl_dg_updater_moment_num_mom(int_mom_up);
  struct gkyl_array *intmom_grid = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  double *red_intmom;
  if (use_gpu)
    red_intmom = gkyl_cu_malloc(sizeof(double[num_mom]));
  else
    red_intmom = gkyl_malloc(sizeof(double[num_mom]));

  gkyl_mom_calc_advance(int_mom_calc, &local, &confLocal, distf, intmom_grid);
  //gkyl_dg_updater_moment_advance(int_mom_up, &local, &confLocal, distf, intmom_grid);
  gkyl_array_reduce_range(red_intmom, intmom_grid, GKYL_SUM, &confLocal);
  double intmom_pre[num_mom];
  if (use_gpu)
    gkyl_cu_memcpy(intmom_pre, red_intmom, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intmom_pre, red_intmom, sizeof(double[num_mom]));

//  printf("\nintmom_pre = %16.14e %16.14e %16.14e %16.14e\n",intmom_pre[0],intmom_pre[1],intmom_pre[2],intmom_pre[3]);
//  // Write distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_positivity_shift_vlasov_1x2v_p%d_pre.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf, fname0);

  // Run the positivity shift. First time it sets ffloor in the pos_shift updater.
  struct gkyl_array *m0 = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *ps_delta_m0 = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *deltaf;
  deltaf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_array_set(deltaf, -1.0, distf);

  struct gkyl_positivity_shift_vlasov* pos_shift = gkyl_positivity_shift_vlasov_new(confBasis,
    basis, grid, &confLocal_ext, use_gpu);
  gkyl_positivity_shift_vlasov_advance(pos_shift, &confLocal, &local, distf, m0, ps_delta_m0);

//  // Commenting this out as we are now using f=0 as the floor.
//  // Project distf and apply the positivity shift again (using new ffloor).
//  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
//  gkyl_array_copy(distf, distf_ho);
//  gkyl_positivity_shift_vlasov_advance(pos_shift, &confLocal, &local, distf, m0, ps_delta_m0);

  // Compute delta f:
  gkyl_array_accumulate(deltaf, 1.0, distf);

  // Compute M0 after the positivity shift.
  struct gkyl_array *m0_post = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  gkyl_mom_calc_advance(mom_calc, &local, &confLocal, distf, m0_post);
  //gkyl_dg_updater_moment_advance(m0_mom_up, &local, &confLocal, distf, m0_post);

//  // Write m0 to file after the positivity shift.
//  char fname1M0[1024];
//  sprintf(fname1M0, "ctest_positivity_shift_vlasov_1x2v_p%d_m0_post.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, m0_post, fname1M0);

  // Compute the integrated moments after the positivity shift.
  gkyl_mom_calc_advance(int_mom_calc, &local, &confLocal, distf, intmom_grid);
  //gkyl_dg_updater_moment_advance(int_mom_up, &local, &confLocal, distf, intmom_grid);
  gkyl_array_reduce_range(red_intmom, intmom_grid, GKYL_SUM, &confLocal);
  double intmom_post[num_mom];
  if (use_gpu)
    gkyl_cu_memcpy(intmom_post, red_intmom, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intmom_post, red_intmom, sizeof(double[num_mom]));

  // Compute the integrated moments of the shift.
  struct gkyl_array *ps_intmom_grid = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  gkyl_mom_calc_advance(int_mom_calc, &local, &confLocal, deltaf, ps_intmom_grid);
  //gkyl_dg_updater_moment_advance(int_mom_up, &local, &confLocal, deltaf, ps_intmom_grid);
  gkyl_array_reduce_range(red_intmom, ps_intmom_grid, GKYL_SUM, &confLocal);
  double intmom_shift[num_mom];
  if (use_gpu)
    gkyl_cu_memcpy(intmom_shift, red_intmom, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intmom_shift, red_intmom, sizeof(double[num_mom]));

//  printf("intmom_post = %16.14e %16.14e %16.14e %16.14e\n",intmom_post[0],intmom_post[1],intmom_post[2],intmom_post[3]);
//  printf("intmom_shift = %16.14e %16.14e %16.14e %16.14e\n",intmom_shift[0],intmom_shift[1],intmom_shift[2],intmom_shift[3]);
//  // Write distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_positivity_shift_vlasov_1x2v_p%d_post.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf, fname1);
 
  // Check the integrated moments.
  TEST_CHECK( gkyl_compare( intmom_shift[0],-4.22405166400353e-16, 1e-10));
  TEST_MSG("intmom_shift[0]: produced: %.14e | expected: %.14e", intmom_shift[0],-4.22405166400353e-16);
  TEST_CHECK( gkyl_compare( intmom_shift[1], 4.16333634234434e-16, 1e-10));
  TEST_MSG("intmom_shift[1]: produced: %.14e | expected: %.14e", intmom_shift[1], 4.16333634234434e-16);
  TEST_CHECK( gkyl_compare( intmom_shift[2], 1.38777878078145e-17, 1e-10));
  TEST_MSG("intmom_shift[2]: produced: %.14e | expected: %.14e", intmom_shift[2], 1.38777878078145e-17);
  // note: half the original value because this is now <H> instead of <T>.
  TEST_CHECK( gkyl_compare( intmom_shift[3], 1.04952716460250e+01, 1e-10));
  TEST_MSG("intmom_shift[3]: produced: %.14e | expected: %.14e", intmom_shift[3], 1.04952716460250e+01);

  gkyl_array_release(distf);
  gkyl_array_release(deltaf);
  gkyl_array_release(intmom_grid);
  gkyl_array_release(ps_intmom_grid);
  gkyl_array_release(ps_delta_m0);
  gkyl_array_release(m0);
  gkyl_array_release(m0_pre);
  gkyl_array_release(m0_post);
  if (use_gpu) {
    gkyl_array_release(distf_ho);
    gkyl_cu_free(red_intmom);
  } else{
    gkyl_free(red_intmom);
  }
  gkyl_proj_on_basis_release(proj_distf);
  gkyl_mom_type_release(m0_mom_up);
  gkyl_mom_type_release(int_mom_up);
  gkyl_positivity_shift_vlasov_release(pos_shift);
  gkyl_vlasov_velocity_map_release(vel_map);
}

void test_1x2v_ho()
{
  test_1x2v(1, false);
}

void test_1x2v_dev()
{
  test_1x2v(1, true);
}

TEST_LIST = {
  { "test_1x2v_ho", test_1x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_dev", test_1x2v_dev },
#endif
  { NULL, NULL },
};
