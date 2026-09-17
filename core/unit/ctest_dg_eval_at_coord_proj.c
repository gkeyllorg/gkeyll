#include <math.h>
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_dg_eval_at_coord_proj.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // Allocate array (filled with zeros).
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
    : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
eval_f_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0 + 2.0*x;
}

static void
test_1x_to_scalar(int poly_order, bool use_gpu)
{
  // Evaluate a 1D DG field at a fixed x0, reducing to a scalar.
  // Verified by comparing against the analytic value f(x0) * (1/sqrt(2)).

  // 1D donor grid.
  double lower_do[] = {0.0}, upper_do[] = {1.0};
  int cells_do[] = {2};

  double x0 = 0.75; // lies in the second x-cell [0.5, 1.0]

  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // Donor basis and range.
  struct gkyl_basis basis_do;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);

  int ghost_do[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &local_ext_do, &local_do);

  // Project f(x) onto the donor basis.
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_1x, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Scalar target: 1D range with a single cell and a 1-component array.
  int lower_tar[] = {1}, upper_tar[] = {1};
  struct gkyl_range local_tar;
  gkyl_range_init(&local_tar, 1, lower_tar, upper_tar);

  struct gkyl_array *ftar = mkarr(use_gpu, 1, local_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, 1, 1)
                                      : gkyl_array_acquire(ftar);

  // Create the updater: evaluate in x (dir 0).
  int eval_dirs[] = {0};
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    basis_do.ndim, &basis_do, num_eval_dirs, eval_dirs, use_gpu);

  // Apply updater at x = x0.
  double eval_coords[] = {x0};
  bool pick_lower[] = {false};
  int known_index[] = {-1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Reference: f(x0) / sqrt(2).
  double fref;
  eval_f_1x(0.0, (const double[]){x0}, &fref, NULL);
  fref /= sqrt(2.0);

  gkyl_array_copy(ftar_ho, ftar);
  const double *ftar_c = gkyl_array_cfetch(ftar_ho, 0);
  TEST_CHECK(gkyl_compare(fref, ftar_c[0], 1e-14));
  TEST_MSG(" Expected %.6e | Got %.6e\n", fref, ftar_c[0]);

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
}

struct dg_evproj_tst_ctx {
  int ndim_do;
  int num_eval_dirs;
  int eval_dirs[GKYL_MAX_DIM];
  double eval_coords[GKYL_MAX_DIM];
  int ndim_tar;
  int dirs_tar[GKYL_MAX_DIM];
};

static void
eval_f_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];

  fout[0] = 1.0 + 2.0*x + 3.0*y + 4.0*x*y;
}

static void
eval_f_2x_at_coord(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct dg_evproj_tst_ctx *params = (struct dg_evproj_tst_ctx *) ctx;

  int ndim_do = params->ndim_do;
  int num_eval_dirs = params->num_eval_dirs;
  int *eval_dirs = params->eval_dirs;
  double *eval_coords = params->eval_coords;
  int ndim_tar = params->ndim_tar;
  int *dirs_tar = params->dirs_tar;

  double xn_p[GKYL_MAX_DIM];
  for (int d=0; d<num_eval_dirs; d++)
    xn_p[eval_dirs[d]] = eval_coords[d];

  for (int d=0; d<ndim_tar; d++)
    xn_p[dirs_tar[d]] = xn[d];

  eval_f_2x(t, xn_p, fout, ctx);
}

static void
test_2x_ev_at_1dcoord(int ndim_do, const double *lower_do, const double *upper_do, const int *cells_do,
  int num_eval_dirs, const int *eval_dirs, const double *eval_coords, int poly_order, bool use_gpu)
{
  // Project a 2D DG field onto 1D by evaluating in x at a fixed
  // computational coordinate x0. Verified by comparing against a direct 1D
  // projection of f at eval_coords.

  // Donor grid.
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // Target grid.
  int ndim_tar = ndim_do - num_eval_dirs;
  int dirs_tar[ndim_tar];
  int c = 0;
  for (int d=0; d<ndim_do; d++) {
    bool in_eval_dirs = false;
    for (int i=0; i<num_eval_dirs; i++) {
      if (d == eval_dirs[i]) {
        in_eval_dirs = true;
        break;
      }
    }

    if (!in_eval_dirs)
      dirs_tar[c++] = d;
  }
  double lower_tar[ndim_tar], upper_tar[ndim_tar];
  int cells_tar[ndim_tar];
  for (int d=0; d<ndim_tar; d++) {
    lower_tar[d] = lower_do[dirs_tar[d]];
    upper_tar[d] = upper_do[dirs_tar[d]];
    cells_tar[d] = cells_do[dirs_tar[d]];
  }
  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim_tar, lower_tar, upper_tar, cells_tar);

  // Basis.
  struct gkyl_basis basis_do, basis_tar;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);
  gkyl_cart_modal_serendip(&basis_tar, ndim_tar, poly_order);

  // Ranges.
  int num_ghost[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, num_ghost, &local_ext_do, &local_do);

  struct gkyl_range local_tar, local_ext_tar;
  gkyl_create_grid_ranges(&grid_tar, num_ghost, &local_ext_tar, &local_tar);

  // Project f(x,y).
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_2x, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Context for projecting reference solution.
  struct dg_evproj_tst_ctx eval_params;
  eval_params.ndim_do = ndim_do;
  eval_params.num_eval_dirs = num_eval_dirs;
  for (int d=0; d<num_eval_dirs; d++) {
    eval_params.eval_dirs[d] = eval_dirs[d];
    eval_params.eval_coords[d] = eval_coords[d];
  }
  eval_params.ndim_tar = ndim_tar;
  for (int d=0; d<ndim_tar; d++)
    eval_params.dirs_tar[d] = dirs_tar[d];

  // Project f at eval_coords[0] directly onto the target basis.
  struct gkyl_array *fref_ho = mkarr(false, basis_tar.num_basis, local_ext_tar.volume);
  gkyl_proj_on_basis *proj_tar = gkyl_proj_on_basis_new(&grid_tar, &basis_tar,
    poly_order+1, 1, eval_f_2x_at_coord, &eval_params);
  gkyl_proj_on_basis_advance(proj_tar, 0.0, &local_tar, fref_ho);

  // Evaluate eval_coords and project.
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    basis_do.ndim, &basis_do, num_eval_dirs, eval_dirs, use_gpu);

  struct gkyl_array *ftar = mkarr(use_gpu, basis_tar.num_basis, local_ext_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, ftar->ncomp, ftar->size)
                                      : gkyl_array_acquire(ftar);
  bool pick_lower[] = {false, false, false, false, false, false};
  int known_index[] = {-1, -1, -1, -1, -1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Check answer.
  gkyl_array_copy(ftar_ho, ftar);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_tar);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local_tar, iter.idx);
    const double *ftar_c = gkyl_array_cfetch(ftar_ho, linidx);
    const double *fref_c = gkyl_array_cfetch(fref_ho, linidx);
    for (int k = 0; k < basis_tar.num_basis; k++) {
      TEST_CHECK(gkyl_compare(fref_c[k], ftar_c[k], 1e-14));
      TEST_MSG(" k:%d | Expected %.6e | Got %6.e\n", k, fref_c[k], ftar_c[k]);
    }
  }

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_proj_on_basis_release(proj_tar);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
  gkyl_array_release(fref_ho);
}

static void
test_2x_to_1x(int poly_order, bool use_gpu)
{
  // Project a 2D DG field onto 1D by evaluating in x at a fixed
  // computational coordinate x0. Verified by comparing against a direct 1D
  // projection of f at eval_coords.

  double lower_do[] = {0.0, 0.0}, upper_do[] = {1.0, 1.0};
  int cells_do[] = {2, 2};
  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);

  int eval_dirs[1];
  double eval_coords[1];
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);

  // Evaluate at x.
  eval_dirs[0] = 0;
  eval_coords[0] = 0.75; // lies in the second x-cell [0.5, 1.0]
  test_2x_ev_at_1dcoord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);

  // Evaluate at y.
  eval_dirs[0] = 1;
  eval_coords[0] = 0.25; // lies in the second x-cell [0.5, 1.0]
  test_2x_ev_at_1dcoord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);
}

static void
test_2x_to_scalar(int poly_order, bool use_gpu)
{
  // Evaluate a 2D DG field at a fixed (x0, y0), reducing to a scalar.
  // Verified by comparing against the analytic value f(x0, y0).

  // 2D donor grid.
  double lower_do[] = {0.0, 0.0}, upper_do[] = {1.0, 1.0};
  int cells_do[] = {2, 2};

  double x0 = 0.75, y0 = 0.25; // x0 in [0.5, 1.0], y0 in [0.0, 0.5]

  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // Donor basis and range.
  struct gkyl_basis basis_do;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);

  int ghost_do[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &local_ext_do, &local_do);

  // Project f(x,y) onto the donor basis.
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_2x, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Target range with 1 cell.
  int lower_tar[] = {1}, upper_tar[] = {1};
  struct gkyl_range local_tar;
  gkyl_range_init(&local_tar, 1, lower_tar, upper_tar);

  struct gkyl_array *ftar = mkarr(use_gpu, 1, local_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, ftar->ncomp, ftar->size)
                                      : gkyl_array_acquire(ftar);

  // Evaluate and project field.
  int eval_dirs[] = {0, 1};
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    basis_do.ndim, &basis_do, num_eval_dirs, eval_dirs, use_gpu);

  // Apply updater at (x0, y0).
  double eval_coords[] = {x0, y0};
  bool pick_lower[] = {false, false};
  int known_index[] = {-1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Check the results.
  double fref;
  eval_f_2x(0.0, (const double[2]){x0,y0}, &fref, 0);
  // Apply normalization for p0.
  fref *= 1.0/sqrt(2.0);

  gkyl_array_copy(ftar_ho, ftar);
  const double *ftar_c = gkyl_array_cfetch(ftar_ho, 0);
  TEST_CHECK(gkyl_compare(fref, ftar_c[0], 1e-14));
  TEST_MSG(" Expected %.6e | Got %.6e\n", fref, ftar_c[0]);

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
}

static void
eval_f_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  fout[0] = 1.0 + 2.0*x + 3.0*y + 4.0*z
    + 5.0*x*y + 6.0*x*z + 7.0*y*z + 8.0*x*y*z;
}

static void
eval_f_3x_at_coord(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct dg_evproj_tst_ctx *params = (struct dg_evproj_tst_ctx *) ctx;

  int num_eval_dirs = params->num_eval_dirs;
  int *eval_dirs = params->eval_dirs;
  double *eval_coords = params->eval_coords;
  int ndim_tar = params->ndim_tar;
  int *dirs_tar = params->dirs_tar;

  double xn_p[GKYL_MAX_DIM];
  for (int d=0; d<num_eval_dirs; d++)
    xn_p[eval_dirs[d]] = eval_coords[d];

  for (int d=0; d<ndim_tar; d++)
    xn_p[dirs_tar[d]] = xn[d];

  eval_f_3x(t, xn_p, fout, NULL);
}

static void
test_3x_ev_at_coord(int ndim_do, const double *lower_do, const double *upper_do, const int *cells_do,
  int num_eval_dirs, const int *eval_dirs, const double *eval_coords, int poly_order, bool use_gpu)
{
  // Project a 3D DG field onto a lower-dimensional target by evaluating at
  // fixed coordinates in eval_dirs. Verified by comparing against a direct
  // projection of f at eval_coords onto the target basis.

  // Donor grid.
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // Target grid: directions not in eval_dirs.
  int ndim_tar = ndim_do - num_eval_dirs;
  int dirs_tar[GKYL_MAX_DIM];
  int c = 0;
  for (int d=0; d<ndim_do; d++) {
    bool in_eval_dirs = false;
    for (int i=0; i<num_eval_dirs; i++) {
      if (d == eval_dirs[i]) {
        in_eval_dirs = true;
        break;
      }
    }
    if (!in_eval_dirs)
      dirs_tar[c++] = d;
  }
  double lower_tar[GKYL_MAX_DIM], upper_tar[GKYL_MAX_DIM];
  int cells_tar[GKYL_MAX_DIM];
  for (int d=0; d<ndim_tar; d++) {
    lower_tar[d] = lower_do[dirs_tar[d]];
    upper_tar[d] = upper_do[dirs_tar[d]];
    cells_tar[d] = cells_do[dirs_tar[d]];
  }
  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim_tar, lower_tar, upper_tar, cells_tar);

  // Basis.
  struct gkyl_basis basis_do, basis_tar;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);
  gkyl_cart_modal_serendip(&basis_tar, ndim_tar, poly_order);

  // Ranges.
  int num_ghost[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, num_ghost, &local_ext_do, &local_do);

  struct gkyl_range local_tar, local_ext_tar;
  gkyl_create_grid_ranges(&grid_tar, num_ghost, &local_ext_tar, &local_tar);

  // Project f(x,y,z).
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_3x, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Context for projecting reference solution.
  struct dg_evproj_tst_ctx eval_params;
  eval_params.ndim_do = ndim_do;
  eval_params.num_eval_dirs = num_eval_dirs;
  for (int d=0; d<num_eval_dirs; d++) {
    eval_params.eval_dirs[d] = eval_dirs[d];
    eval_params.eval_coords[d] = eval_coords[d];
  }
  eval_params.ndim_tar = ndim_tar;
  for (int d=0; d<ndim_tar; d++)
    eval_params.dirs_tar[d] = dirs_tar[d];

  // Project f at eval_coords directly onto the target basis.
  struct gkyl_array *fref_ho = mkarr(false, basis_tar.num_basis, local_ext_tar.volume);
  gkyl_proj_on_basis *proj_tar = gkyl_proj_on_basis_new(&grid_tar, &basis_tar,
    poly_order+1, 1, eval_f_3x_at_coord, &eval_params);
  gkyl_proj_on_basis_advance(proj_tar, 0.0, &local_tar, fref_ho);

  // Evaluate eval_coords and project.
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    basis_do.ndim, &basis_do, num_eval_dirs, eval_dirs, use_gpu);

  struct gkyl_array *ftar = mkarr(use_gpu, basis_tar.num_basis, local_ext_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, ftar->ncomp, ftar->size)
                                      : gkyl_array_acquire(ftar);
  bool pick_lower[] = {false, false, false, false, false, false};
  int known_index[] = {-1, -1, -1, -1, -1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Check answer.
  gkyl_array_copy(ftar_ho, ftar);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_tar);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local_tar, iter.idx);
    const double *ftar_c = gkyl_array_cfetch(ftar_ho, linidx);
    const double *fref_c = gkyl_array_cfetch(fref_ho, linidx);
    for (int k = 0; k < basis_tar.num_basis; k++) {
      TEST_CHECK(gkyl_compare(fref_c[k], ftar_c[k], 1e-14));
      TEST_MSG(" k:%d | Expected %.6e | Got %.6e\n", k, fref_c[k], ftar_c[k]);
    }
  }

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_proj_on_basis_release(proj_tar);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
  gkyl_array_release(fref_ho);
}

static void
test_3x_to_2x(int poly_order, bool use_gpu)
{
  double lower_do[] = {0.0, 0.0, 0.0}, upper_do[] = {1.0, 1.0, 1.0};
  int cells_do[] = {2, 2, 2};
  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);

  int eval_dirs[1];
  double eval_coords[1];
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);

  // Evaluate at x.
  eval_dirs[0] = 0;
  eval_coords[0] = 0.75;
  test_3x_ev_at_coord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);

  // Evaluate at y.
  eval_dirs[0] = 1;
  eval_coords[0] = 0.25;
  test_3x_ev_at_coord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);

  // Evaluate at z.
  eval_dirs[0] = 2;
  eval_coords[0] = 0.6;
  test_3x_ev_at_coord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);
}

static void
test_3x_to_1x(int poly_order, bool use_gpu)
{
  double lower_do[] = {0.0, 0.0, 0.0}, upper_do[] = {1.0, 1.0, 1.0};
  int cells_do[] = {2, 2, 2};
  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);

  int eval_dirs[2];
  double eval_coords[2];
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);

  // Evaluate at (x, y).
  eval_dirs[0] = 0;
  eval_dirs[1] = 1;
  eval_coords[0] = 0.75;
  eval_coords[1] = 0.25;
  test_3x_ev_at_coord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);

  // Evaluate at (x, z).
  eval_dirs[0] = 0;
  eval_dirs[1] = 2;
  eval_coords[0] = 0.75;
  eval_coords[1] = 0.6;
  test_3x_ev_at_coord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);

  // Evaluate at (y, z).
  eval_dirs[0] = 1;
  eval_dirs[1] = 2;
  eval_coords[0] = 0.25;
  eval_coords[1] = 0.6;
  test_3x_ev_at_coord(ndim_do, lower_do, upper_do, cells_do, num_eval_dirs,
    eval_dirs, eval_coords, poly_order, use_gpu);
}

static void
test_3x_to_scalar(int poly_order, bool use_gpu)
{
  // Evaluate a 3D DG field at a fixed (x0, y0, z0), reducing to a scalar.
  // Verified by comparing against the analytic value f(x0, y0, z0) / sqrt(2).

  // 3D donor grid.
  double lower_do[] = {0.0, 0.0, 0.0}, upper_do[] = {1.0, 1.0, 1.0};
  int cells_do[] = {2, 2, 2};

  // x0 in [0.5,1.0], y0 in [0.0,0.5], z0 in [0.5,1.0]
  double x0 = 0.75, y0 = 0.25, z0 = 0.6;

  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // Donor basis and range.
  struct gkyl_basis basis_do;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);

  int ghost_do[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &local_ext_do, &local_do);

  // Project f(x,y,z) onto the donor basis.
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_3x, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Scalar target: 1D range with a single cell and a 1-component array.
  int lower_tar[] = {1}, upper_tar[] = {1};
  struct gkyl_range local_tar;
  gkyl_range_init(&local_tar, 1, lower_tar, upper_tar);

  struct gkyl_array *ftar = mkarr(use_gpu, 1, local_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, 1, 1)
                                      : gkyl_array_acquire(ftar);

  // Create the updater: evaluate in x (dir 0), y (dir 1), and z (dir 2).
  int eval_dirs[] = {0, 1, 2};
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    basis_do.ndim, &basis_do, num_eval_dirs, eval_dirs, use_gpu);

  // Apply updater at (x0, y0, z0).
  double eval_coords[] = {x0, y0, z0};
  bool pick_lower[] = {false, false, false};
  int known_index[] = {-1, -1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Reference: f(x0, y0, z0) / sqrt(2).
  double fref;
  eval_f_3x(0.0, (const double[]){x0, y0, z0}, &fref, NULL);
  fref /= sqrt(2.0);

  gkyl_array_copy(ftar_ho, ftar);
  const double *ftar_c = gkyl_array_cfetch(ftar_ho, 0);
  TEST_CHECK(gkyl_compare(fref, ftar_c[0], 1e-14));
  TEST_MSG(" Expected %.6e | Got %.6e\n", fref, ftar_c[0]);

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
}

void test_dg_evproj_1x_to_scalar_p1_ho(void) {
  test_1x_to_scalar(1, false);
}

void test_dg_evproj_1x_to_scalar_p2_ho(void) {
  test_1x_to_scalar(2, false);
}

void test_dg_evproj_2x_to_1x_p1_ho(void) {
  test_2x_to_1x(1, false);
}

void test_dg_evproj_2x_to_1x_p2_ho(void) {
  test_2x_to_1x(2, false);
}

void test_dg_evproj_2x_to_scalar_p1_ho(void) {
  test_2x_to_scalar(1, false);
}

void test_dg_evproj_2x_to_scalar_p2_ho(void) {
  test_2x_to_scalar(2, false);
}

void test_dg_evproj_3x_to_2x_p1_ho(void) {
  test_3x_to_2x(1, false);
}

void test_dg_evproj_3x_to_1x_p1_ho(void) {
  test_3x_to_1x(1, false);
}

void test_dg_evproj_3x_to_scalar_p1_ho(void) {
  test_3x_to_scalar(1, false);
}

#ifdef GKYL_HAVE_CUDA
void test_dg_evproj_1x_to_scalar_p1_dev(void) {
  test_1x_to_scalar(1, true);
}

void test_dg_evproj_1x_to_scalar_p2_dev(void) {
  test_1x_to_scalar(2, true);
}

void test_dg_evproj_2x_to_1x_p1_dev(void) {
  test_2x_to_1x(1, true);
}

void test_dg_evproj_2x_to_1x_p2_dev(void) {
  test_2x_to_1x(2, true);
}

void test_dg_evproj_2x_to_scalar_p1_dev(void) {
  test_2x_to_scalar(1, true);
}

void test_dg_evproj_2x_to_scalar_p2_dev(void) {
  test_2x_to_scalar(2, true);
}

void test_dg_evproj_3x_to_2x_p1_dev(void) {
  test_3x_to_2x(1, true);
}

void test_dg_evproj_3x_to_1x_p1_dev(void) {
  test_3x_to_1x(1, true);
}

void test_dg_evproj_3x_to_scalar_p1_dev(void) {
  test_3x_to_scalar(1, true);
}
#endif

TEST_LIST = {
  { "test_dg_evproj_1x_to_scalar_p1_ho", test_dg_evproj_1x_to_scalar_p1_ho },
  { "test_dg_evproj_1x_to_scalar_p2_ho", test_dg_evproj_1x_to_scalar_p2_ho },
  { "test_dg_evproj_2x_to_1x_p1_ho", test_dg_evproj_2x_to_1x_p1_ho },
  { "test_dg_evproj_2x_to_1x_p2_ho", test_dg_evproj_2x_to_1x_p2_ho },
  { "test_dg_evproj_2x_to_scalar_p1_ho", test_dg_evproj_2x_to_scalar_p1_ho },
  { "test_dg_evproj_2x_to_scalar_p2_ho", test_dg_evproj_2x_to_scalar_p2_ho },
  { "test_dg_evproj_3x_to_2x_p1_ho", test_dg_evproj_3x_to_2x_p1_ho },
  { "test_dg_evproj_3x_to_1x_p1_ho", test_dg_evproj_3x_to_1x_p1_ho },
  { "test_dg_evproj_3x_to_scalar_p1_ho", test_dg_evproj_3x_to_scalar_p1_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_dg_evproj_1x_to_scalar_p1_dev", test_dg_evproj_1x_to_scalar_p1_dev },
  { "test_dg_evproj_1x_to_scalar_p2_dev", test_dg_evproj_1x_to_scalar_p2_dev },
  { "test_dg_evproj_2x_to_1x_p1_dev", test_dg_evproj_2x_to_1x_p1_dev },
  { "test_dg_evproj_2x_to_1x_p2_dev", test_dg_evproj_2x_to_1x_p2_dev },
  { "test_dg_evproj_2x_to_scalar_p1_dev", test_dg_evproj_2x_to_scalar_p1_dev },
  { "test_dg_evproj_2x_to_scalar_p2_dev", test_dg_evproj_2x_to_scalar_p2_dev },
  { "test_dg_evproj_3x_to_2x_p1_dev", test_dg_evproj_3x_to_2x_p1_dev },
  { "test_dg_evproj_3x_to_1x_p1_dev", test_dg_evproj_3x_to_1x_p1_dev },
  { "test_dg_evproj_3x_to_scalar_p1_dev", test_dg_evproj_3x_to_scalar_p1_dev },
#endif
  { NULL, NULL },
};
