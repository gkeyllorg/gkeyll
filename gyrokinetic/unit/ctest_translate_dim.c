#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_translate_dim.h>
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

void create_lower_dim_objects(int cdim_tar, struct gkyl_rect_grid grid_tar, int poly_order,
  struct gkyl_rect_grid *grid, struct gkyl_rect_grid *confGrid,
  struct gkyl_basis *basis, struct gkyl_basis *confBasis,
  struct gkyl_range *confLocal, struct gkyl_range *confLocal_ext,
  struct gkyl_range *local, struct gkyl_range *local_ext)
{
  // Create lower dimensional grid, basis and range based on the target
  // dimensionality and grid.
  const int ndim = grid_tar.ndim-1;
  const int cdim = cdim_tar-1;
  const int vdim = ndim-cdim;

  double confLower[GKYL_MAX_CDIM] = {0.0}, confUpper[GKYL_MAX_CDIM] = {0.0};
  int confCells[GKYL_MAX_CDIM] = {0};
  confLower[cdim-1] = grid_tar.lower[cdim_tar-1];
  confUpper[cdim-1] = grid_tar.upper[cdim_tar-1];
  confCells[cdim-1] = grid_tar.cells[cdim_tar-1];
  if (cdim_tar == 3) {
    confLower[0] = grid_tar.lower[0];
    confUpper[0] = grid_tar.upper[0];
    confCells[0] = grid_tar.cells[0];
  }

  double lower[GKYL_MAX_DIM] = {0.0}, upper[GKYL_MAX_DIM] = {0.0};
  int cells[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) {
    lower[d] = confLower[d];
    upper[d] = confUpper[d];
    cells[d] = confCells[d];
  }
  for (int d=0; d<vdim; d++) {
    lower[cdim+d] = grid_tar.lower[cdim_tar+d];
    upper[cdim+d] = grid_tar.upper[cdim_tar+d];
    cells[cdim+d] = grid_tar.cells[cdim_tar+d];
  }

  // Grids.
  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions.
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(basis, ndim, poly_order);
  gkyl_cart_modal_serendip(confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(confGrid, confGhost, confLocal_ext, confLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  gkyl_create_grid_ranges(grid, ghost, local_ext, local);
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

double den_profile_1x(double n0, double z)
{
  return n0*(1.0+0.3*cos(2.0*M_PI*z));
}

void eval_distf_2x2v_low(double t, const double *xn, double* restrict fout, void *ctx)
{
  // This projects a low-dim (1x2v) distribution for the 2x2v test.
  double x = xn[0], vpar = xn[1], mu = xn[2];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;
  double mass = tctx->mass;
  double n0 = tctx->n0;
  double upar = tctx->upar;
  double vtsq = tctx->temp/mass;
  int vdim = tctx->vdim;

  double den = den_profile_1x(n0, x);

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2)+2.0*mu*B0/mass)/(2.0*vtsq));
}

void
test_2x2v(int poly_order, bool use_gpu)
{
  const int cdim = 2;
  double vpar_max = 6.0;
  double mu_max = 36.0;
  double lower[] = {0.1, -M_PI, -vpar_max, 0.0}, upper[] = {1.0, M_PI, vpar_max, mu_max};
  int cells[] = {2, 4, 6, 4};

  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int vdim = ndim-cdim;

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

  double confLower[GKYL_MAX_CDIM] = {0.0}, confUpper[GKYL_MAX_CDIM] = {0.0};
  int confCells[GKYL_MAX_CDIM] = {0};
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

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

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create lower dimensional grid, basis and ranges.
  struct gkyl_rect_grid grid_low, confGrid_low;
  struct gkyl_basis basis_low, confBasis_low;
  struct gkyl_range confLocal_low, confLocal_ext_low;
  struct gkyl_range local_low, local_ext_low;
  create_lower_dim_objects(cdim, grid, poly_order, &grid_low,  &confGrid_low, &basis_low,
    &confBasis_low, &confLocal_low, &confLocal_ext_low, &local_low, &local_ext_low);

  // Create donor distribution function arrays.
  struct gkyl_array *distf_low_ho, *distf_low;
  distf_low = mkarr(use_gpu, basis_low.num_basis, local_ext_low.volume);
  distf_low_ho = use_gpu? mkarr(false, distf_low->ncomp, distf_low->size)
                        : gkyl_array_acquire(distf_low);

  // Project the donor distribution.
  gkyl_proj_on_basis *proj_distf_low = gkyl_proj_on_basis_new(&grid_low, &basis_low,
    poly_order+1, 1, eval_distf_2x2v_low, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf_low, 0.0, &local_low, distf_low_ho);
  gkyl_array_copy(distf_low, distf_low_ho);

//  // Write distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_translate_dim_2x2v_p%d_low.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid_low, &local_low, NULL, distf_low, fname0);

  // Create target distribution function arrays.
  struct gkyl_array *distf_ho, *distf;
  distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                    : gkyl_array_acquire(distf);

  // Translate the DG coefficients.
  int cdim_do = confGrid_low.ndim;
  int vdim_do = grid_low.ndim - cdim_do;
  struct gkyl_translate_dim* trans_dim_upd = gkyl_translate_dim_new(cdim_do,
    basis_low, cdim, basis, 0, GKYL_NO_EDGE, use_gpu);
  gkyl_translate_dim_advance(trans_dim_upd, &local_low, &local, distf_low, 1, distf);
  gkyl_array_copy(distf_ho, distf);

//  // Write distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_translate_dim_2x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname1);
 
  // How DG coefficients of the higher dim field are mapped to those of the
  // lower dim field. If <0, its amplitude is 0.
  int dg_map[] = {0,-1,1,2,3,-1,-1,4,-1,5,6,-1,-1,-1,7,-1,8,-1,9,10,-1,-1,11,-1};

  // Check coefficients of the higher dimensional field.
  int pidx_do[GKYL_MAX_DIM] = {-1};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {

    // Translate the target idx to the donor idx:
    for (int d=0; d<cdim_do-1; d++) pidx_do[d] = iter.idx[d]; 
    pidx_do[cdim_do-1] = iter.idx[cdim-1]; 
    for (int d=0; d<vdim_do; d++) pidx_do[cdim_do+d] = iter.idx[cdim+d]; 

    long plinidx_tar = gkyl_range_idx(&local, iter.idx);
    long plinidx_do = gkyl_range_idx(&local_low, pidx_do);

    const double *f_c = gkyl_array_cfetch(distf_ho, plinidx_tar);
    const double *flow_c = gkyl_array_cfetch(distf_low_ho, plinidx_do);

    for (int k=0; k<basis.num_basis; k++) {
      if (dg_map[k] < 0)
        TEST_CHECK( gkyl_compare( f_c[k], 0.0, 1e-16));
      else
        TEST_CHECK( gkyl_compare( f_c[k], 1.4142135623730951*flow_c[dg_map[k]], 1e-14));
    }
  }

  gkyl_translate_dim_release(trans_dim_upd);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_distf_low);
  gkyl_array_release(distf_low);
  gkyl_array_release(distf_low_ho);
}

double den_profile_2x(double n0, double x, double z)
{
  return n0*(1.0-x)*(1.0+0.3*cos(2.0*M_PI*z));
}

void eval_distf_3x2v_low(double t, const double *xn, double* restrict fout, void *ctx)
{
  // This projects a low-dim (2x2v) distribution for the 3x2v test.
  double x = xn[0], y = xn[1], vpar = xn[2], mu = xn[3];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;
  double mass = tctx->mass;
  double n0 = tctx->n0;
  double upar = tctx->upar;
  double vtsq = tctx->temp/mass;
  int vdim = tctx->vdim;

  double den = den_profile_2x(n0, x, y);

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2)+2.0*mu*B0/mass)/(2.0*vtsq));
}

void
test_3x2v(int poly_order, bool use_gpu)
{
  const int cdim = 3;
  double vpar_max = 6.0;
  double mu_max = 36.0;
  double lower[] = {0.1, -0.3, -M_PI, -vpar_max, 0.0}, upper[] = {1.0, 0.3, M_PI, vpar_max, mu_max};
  int cells[] = {2, 2, 4, 6, 4};

  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int vdim = ndim-cdim;

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

  double confLower[GKYL_MAX_CDIM] = {0.0}, confUpper[GKYL_MAX_CDIM] = {0.0};
  int confCells[GKYL_MAX_CDIM] = {0};
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

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

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create lower dimensional grid, basis and ranges.
  struct gkyl_rect_grid grid_low, confGrid_low;
  struct gkyl_basis basis_low, confBasis_low;
  struct gkyl_range confLocal_low, confLocal_ext_low;
  struct gkyl_range local_low, local_ext_low;
  create_lower_dim_objects(cdim, grid, poly_order, &grid_low,  &confGrid_low, &basis_low,
    &confBasis_low, &confLocal_low, &confLocal_ext_low, &local_low, &local_ext_low);

  // Create donor distribution function arrays.
  struct gkyl_array *distf_low_ho, *distf_low;
  distf_low = mkarr(use_gpu, basis_low.num_basis, local_ext_low.volume);
  distf_low_ho = use_gpu? mkarr(false, distf_low->ncomp, distf_low->size)
                        : gkyl_array_acquire(distf_low);

  // Project the donor distribution.
  gkyl_proj_on_basis *proj_distf_low = gkyl_proj_on_basis_new(&grid_low, &basis_low,
    poly_order+1, 1, eval_distf_3x2v_low, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf_low, 0.0, &local_low, distf_low_ho);
  gkyl_array_copy(distf_low, distf_low_ho);

//  // Write distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_translate_dim_3x2v_p%d_low.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid_low, &local_low, NULL, distf_low, fname0);

  // Create target distribution function arrays.
  struct gkyl_array *distf_ho, *distf;
  distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                    : gkyl_array_acquire(distf);

  // Translate the DG coefficients.
  int cdim_do = confGrid_low.ndim;
  int vdim_do = grid_low.ndim - cdim_do;
  struct gkyl_translate_dim* trans_dim_upd = gkyl_translate_dim_new(cdim_do,
    basis_low, cdim, basis, 0, GKYL_NO_EDGE, use_gpu);
  gkyl_translate_dim_advance(trans_dim_upd, &local_low, &local, distf_low, 1, distf);
  gkyl_array_copy(distf_ho, distf);

//  // Write distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_translate_dim_3x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname1);
 
  // How DG coefficients of the higher dim field are mapped to those of the
  // lower dim field. If <0, its amplitude is 0.
  int dg_map[] = {
    0,1,-1,2,3,4,-1,5,-1,6,-1,7,8,-1,9,10,-1,-1,11,-1,-1,12,-1,13,-1,14,
    -1,-1,-1,15,-1,-1,16,17,-1,18,19,-1,20,-1,21,-1,22,-1,-1,23,-1,-1
  };

  // Check coefficients of the higher dimensional field.
  int pidx_do[GKYL_MAX_DIM] = {-1};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {

    // Translate the target idx to the donor idx:
    for (int d=0; d<cdim_do-1; d++) pidx_do[d] = iter.idx[d]; 
    pidx_do[cdim_do-1] = iter.idx[cdim-1]; 
    for (int d=0; d<vdim_do; d++) pidx_do[cdim_do+d] = iter.idx[cdim+d]; 

    long plinidx_tar = gkyl_range_idx(&local, iter.idx);
    long plinidx_do = gkyl_range_idx(&local_low, pidx_do);

    const double *f_c = gkyl_array_cfetch(distf_ho, plinidx_tar);
    const double *flow_c = gkyl_array_cfetch(distf_low_ho, plinidx_do);

    for (int k=0; k<basis.num_basis; k++) {
      if (dg_map[k] < 0)
        TEST_CHECK( gkyl_compare( f_c[k], 0.0, 1e-16));
      else
        TEST_CHECK( gkyl_compare( f_c[k], 1.4142135623730951*flow_c[dg_map[k]], 1e-14));
    }
  }

  gkyl_translate_dim_release(trans_dim_upd);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_distf_low);
  gkyl_array_release(distf_low);
  gkyl_array_release(distf_low_ho);
}

void test_2x2v_ho()
{
  test_2x2v(1, false);
}

void test_2x2v_dev()
{
  test_2x2v(1, true);
}

void test_3x2v_ho()
{
  test_3x2v(1, false);
}

void test_3x2v_dev()
{
  test_3x2v(1, true);
}

TEST_LIST = {
  { "test_2x2v_ho", test_2x2v_ho },
  { "test_3x2v_ho", test_3x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_2x2v_dev", test_2x2v_dev },
  { "test_3x2v_dev", test_3x2v_dev },
#endif
  { NULL, NULL },
};
