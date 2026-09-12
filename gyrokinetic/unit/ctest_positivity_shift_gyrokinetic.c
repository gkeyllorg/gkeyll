#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>

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
  double upar; // Parallel flow speed.
  double temp; // Temperature.
  double mass; // Species mass.
  double B0; // Magnetic field.
  int vdim; // Number of velocity space dimensions.
  double vpar_max; // Maximum vpar of the grid.
  double mu_max; // Maximum mu of the grid.
};

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void eval_bfield_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
}

void eval_bmag_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct test_ctx *tctx = ctx;
  fout[0] = tctx->B0;
}

void eval_distf_1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1], mu = xn[2];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;
  double mass = tctx->mass;
  double n0 = tctx->n0;
  double upar = tctx->upar;
  double vtsq = tctx->temp/mass;
  int vdim = tctx->vdim;

  fout[0] = (n0/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2)+2.0*mu*B0/mass)/(2.0*vtsq));

  // Intentionally set some places to be negative.
  if (fabs(vpar) > 0.8*tctx->vpar_max || mu > 0.8*tctx->mu_max)
    fout[0] = -0.2 * (n0/pow(2.0*M_PI*vtsq,vdim/2.0));
}

void
test_1x2v(int poly_order, bool use_gpu)
{
  int cdim = 1;
  double vpar_max = 6.0;
  double mu_max = 36.0;
  double lower[] = {0.1, -vpar_max, 0.0}, upper[] = {1.0, vpar_max, mu_max};
  int cells[] = {2, 12, 8};

  int ndim = sizeof(cells)/sizeof(cells[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar = 0, // Parallel flow speed.
    .temp = 2.75, // Temperature.
    .mass = 1.0, // Species mass.
    .B0 = 1.0, // Magnetic field.
    .vdim = vdim, // Number of velocity space dimensions.
    .vpar_max = vpar_max, // Maximum vpar of the grid.
    .mu_max = mu_max, // Maximum mu of the grid.
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
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag_ho, *bmag;
  bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  bmag_ho = bmag;
  if (use_gpu)
    bmag_ho = mkarr(false, bmag->ncomp, bmag->size);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

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

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0},  .mapc2p = mapc2p,  .c2p_ctx = 0,
    .bfield_func = eval_bfield_1x,  .bfield_ctx = &proj_ctx,
    .basis = confBasis,  .grid = confGrid,
    .local = confLocal,  .local_ext = confLocal_ext,
    .global = confLocal, .global_ext = confLocal_ext,
    .position_map = pmap,
  };
  int geo_ghost[3] = {1, 1, 1};
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, geo_ghost, &geometry_input.geo_global_ext, &geometry_input.geo_global);
  memcpy(&geometry_input.geo_local, &geometry_input.geo_global, sizeof(struct gkyl_range));
  memcpy(&geometry_input.geo_local_ext, &geometry_input.geo_global_ext, sizeof(struct gkyl_range));
  // Deflate geometry.
  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);
  // If we are on the gpu, copy from host.
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Compute the integrated moments of the original f.
  struct gkyl_dg_updater_moment *int_mom_up = gkyl_dg_updater_moment_gyrokinetic_new(
    &grid, &confBasis, &basis, &confLocal, proj_ctx.mass, 0, gvm, gk_geom, NULL, GKYL_F_MOMENT_M0M1M2PARM2PERP, true, use_gpu);

  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(int_mom_up);
  struct gkyl_array *intmom_grid = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  double *red_intmom;
  if (use_gpu)
    red_intmom = gkyl_cu_malloc(sizeof(double[2+vdim]));
  else
    red_intmom = gkyl_malloc(sizeof(double[2+vdim]));

  gkyl_dg_updater_moment_gyrokinetic_advance(int_mom_up, &local, &confLocal, distf, intmom_grid);
  gkyl_array_reduce_range(red_intmom, intmom_grid, GKYL_SUM, &confLocal);
  double intmom_pre[2+vdim];
  if (use_gpu)
    gkyl_cu_memcpy(intmom_pre, red_intmom, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intmom_pre, red_intmom, sizeof(double[2+vdim]));

//  printf("\nintmom_pre = %16.14e %16.14e %16.14e %16.14e\n",intmom_pre[0],intmom_pre[1],intmom_pre[2],intmom_pre[3]);
//  // Write distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_positivity_shift_gyrokinetic_1x2v_p%d_pre.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf, fname0);

  // Run the positivity shift. First time it sets ffloor in the pos_shift updater.
  struct gkyl_array *m0 = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *ps_delta_m0 = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *deltaf;
  deltaf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_array_set(deltaf, -1.0, distf);


  struct gkyl_positivity_shift_gyrokinetic* pos_shift = gkyl_positivity_shift_gyrokinetic_new(confBasis,
    basis, grid, proj_ctx.mass, gk_geom, gvm, &confLocal_ext, use_gpu);
  gkyl_positivity_shift_gyrokinetic_advance(pos_shift, &confLocal, &local, distf, m0, ps_delta_m0);

  // Project distf and apply the positivity shift again (using new ffloor).
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  gkyl_positivity_shift_gyrokinetic_advance(pos_shift, &confLocal, &local, distf, m0, ps_delta_m0);

  // Compute delta f:
  gkyl_array_accumulate(deltaf, 1.0, distf);

  // Compute the integrated moments after the positivity shift.
  gkyl_dg_updater_moment_gyrokinetic_advance(int_mom_up, &local, &confLocal, distf, intmom_grid);
  gkyl_array_reduce_range(red_intmom, intmom_grid, GKYL_SUM, &confLocal);
  double intmom_post[2+vdim];
  if (use_gpu)
    gkyl_cu_memcpy(intmom_post, red_intmom, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intmom_post, red_intmom, sizeof(double[2+vdim]));

  // Compute the integrated moments of the shift.
  struct gkyl_array *ps_intmom_grid = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  gkyl_dg_updater_moment_gyrokinetic_advance(int_mom_up, &local, &confLocal, deltaf, ps_intmom_grid);
  gkyl_array_reduce_range(red_intmom, ps_intmom_grid, GKYL_SUM, &confLocal);
  double intmom_shift[2+vdim];
  if (use_gpu)
    gkyl_cu_memcpy(intmom_shift, red_intmom, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intmom_shift, red_intmom, sizeof(double[2+vdim]));

//  printf("intmom_post = %16.14e %16.14e %16.14e %16.14e\n",intmom_post[0],intmom_post[1],intmom_post[2],intmom_post[3]);
//  printf("intmom_shift = %16.14e %16.14e %16.14e %16.14e\n",intmom_shift[0],intmom_shift[1],intmom_shift[2],intmom_shift[3]);
//  // Write distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_positivity_shift_gyrokinetic_1x2v_p%d_post.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf, fname1);
 
  // Check the integrated moments.
  TEST_CHECK( gkyl_compare( intmom_shift[0], 9.13090909090910e+00, 1e-10));
  TEST_MSG("intmom_shift[0]: produced: %.14e | expected: %.14e", intmom_shift[0], 9.13090909090910e+00);
  TEST_CHECK( gkyl_compare( intmom_shift[1], 8.61055942680578e-16, 1e-10));
  TEST_MSG("intmom_shift[1]: produced: %.14e | expected: %.14e", intmom_shift[1], 8.61055942680578e-16);
  TEST_CHECK( gkyl_compare( intmom_shift[2], 1.79770909090909e+02, 1e-10));
  TEST_MSG("intmom_shift[2]: produced: %.14e | expected: %.14e", intmom_shift[2], 1.79770909090909e+02);
  TEST_CHECK( gkyl_compare( intmom_shift[3], 4.58457166783993e+02, 1e-10));
  TEST_MSG("intmom_shift[3]: produced: %.14e | expected: %.14e", intmom_shift[3], 4.58457166783993e+02);

  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(deltaf);
  gkyl_array_release(intmom_grid);
  gkyl_array_release(ps_intmom_grid);
  gkyl_array_release(ps_delta_m0);
  gkyl_array_release(m0);
  if (use_gpu) {
    gkyl_array_release(bmag_ho);
    gkyl_array_release(distf_ho);
    gkyl_cu_free(red_intmom);
  } else{
    gkyl_free(red_intmom);
  }
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
  gkyl_dg_updater_moment_gyrokinetic_release(int_mom_up);
  gkyl_positivity_shift_gyrokinetic_release(pos_shift);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_velocity_map_release(gvm);
  gkyl_position_map_release(pmap);
}

void test_positivity_shift_1x2v_ho()
{
  test_1x2v(1, false);
}

void test_positivity_shift_1x2v_dev()
{
  test_1x2v(1, true);
}

TEST_LIST = {
  { "test_positivity_shift_1x2v_ho", test_positivity_shift_1x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_positivity_shift_1x2v_dev", test_positivity_shift_1x2v_dev },
#endif
  { NULL, NULL },
};
