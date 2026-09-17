/* Test gkyl_array_average updater on decomposed domains.
Each test computes a (weighted) average on a decomposed domain and compares, on every rank,
the local part of the result to the average computed without decomposition on the global range.
*/

#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <mpi.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static void
eval_func(double t, const double *xn, double* restrict fout, void *ctx)
{
  int ndim = *(int *) ctx;
  double x = xn[0], y = xn[1], z = ndim > 2 ? xn[2] : 0.0;
  fout[0] = cos(0.5*x)*y + x*y*y + sin(z);
}

static void
eval_weight(double t, const double *xn, double* restrict fout, void *ctx)
{
  int ndim = *(int *) ctx;
  double x = xn[0], y = xn[1], z = ndim > 2 ? xn[2] : 0.0;
  fout[0] = 1.0 + 0.1*x*x*(y+4.0) + 0.05*z*z;
}

static struct gkyl_array*
project(struct gkyl_rect_grid *grid, struct gkyl_basis *basis, struct gkyl_range *range_ext,
  evalf_t func, int *ndim)
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, range_ext->volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(grid, basis, basis->poly_order+1, 1, func, ndim);
  gkyl_proj_on_basis_advance(proj, 0.0, range_ext, arr);
  gkyl_proj_on_basis_release(proj);
  return arr;
}

// Build the reduced range made of the non-averaged dimensions of rng (a single cell if none remain).
static void
reduced_range(const struct gkyl_range *rng, const int *avg_dim, int ndim, int nghost,
  struct gkyl_range *red, struct gkyl_range *red_ext)
{
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM], lo_ext[GKYL_MAX_DIM], up_ext[GKYL_MAX_DIM];
  int k = 0;
  for (int d=0; d<ndim; ++d) {
    if (!avg_dim[d]) {
      lo[k] = rng->lower[d];  up[k] = rng->upper[d];
      lo_ext[k] = lo[k]-nghost;  up_ext[k] = up[k]+nghost;
      k++;
    }
  }
  if (k == 0) {
    lo[0] = up[0] = lo_ext[0] = up_ext[0] = 1;
    k = 1;
  }
  gkyl_range_init(red, k, lo, up);
  gkyl_range_init(red_ext, k, lo_ext, up_ext);
}

// Returns the max relative difference between the decomposed and global averages
// on this rank's cells. If use_comm is false, the updater is used without communicator.
static double
run_avg(int ndim, double *lower, double *upper, int *cells, int *cuts, int *avg_dim,
  int poly_order, bool weighted, bool use_comm)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int nremain = 0;
  for (int d=0; d<ndim; ++d) nremain += 1-avg_dim[d];
  struct gkyl_basis basis_avg;
  gkyl_cart_modal_serendip(&basis_avg, nremain > 0 ? nremain : 1, poly_order);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; ++d) ghost[d] = 1;

  struct gkyl_range global, global_ext;
  gkyl_create_grid_ranges(&grid, ghost, &global_ext, &global);

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(ndim, cuts, &global);
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], ghost, &local_ext, &local);

  struct gkyl_range global_avg, global_avg_ext, local_avg, local_avg_ext;
  reduced_range(&global, avg_dim, ndim, 1, &global_avg, &global_avg_ext);
  reduced_range(&local, avg_dim, ndim, 1, &local_avg, &local_avg_ext);

  // reference: average computed on the undecomposed global range
  struct gkyl_array *f_glob = project(&grid, &basis, &global_ext, eval_func, &ndim);
  struct gkyl_array *w_glob = project(&grid, &basis, &global_ext, eval_weight, &ndim);
  struct gkyl_array *avg_ref = gkyl_array_new(GKYL_DOUBLE, basis_avg.num_basis, global_avg_ext.volume);
  struct gkyl_array_average *up_ref = gkyl_array_average_inew( &(struct gkyl_array_average_inp) {
      .grid = &grid, .basis = basis, .basis_avg = basis_avg,
      .local = &global, .local_avg = &global_avg, .local_avg_ext = &global_avg_ext,
      .weight = weighted ? w_glob : NULL, .avg_dim = avg_dim, .use_gpu = false,
    }
  );
  gkyl_array_average_advance(up_ref, f_glob, avg_ref);

  // decomposed average
  struct gkyl_array *f_loc = project(&grid, &basis, &local_ext, eval_func, &ndim);
  struct gkyl_array *w_loc = project(&grid, &basis, &local_ext, eval_weight, &ndim);
  struct gkyl_array *avg = gkyl_array_new(GKYL_DOUBLE, basis_avg.num_basis, local_avg_ext.volume);
  struct gkyl_array_average *up = gkyl_array_average_inew( &(struct gkyl_array_average_inp) {
      .grid = &grid, .basis = basis, .basis_avg = basis_avg,
      .local = &local, .local_avg = &local_avg, .local_avg_ext = &local_avg_ext,
      .weight = weighted ? w_loc : NULL, .avg_dim = avg_dim, .use_gpu = false,
      .comm = use_comm ? comm : NULL, .global_avg = &global_avg,
    }
  );
  gkyl_array_average_advance(up, f_loc, avg);

  double max_err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_avg);
  while (gkyl_range_iter_next(&iter)) {
    const double *a = gkyl_array_cfetch(avg, gkyl_range_idx(&local_avg, iter.idx));
    const double *r = gkyl_array_cfetch(avg_ref, gkyl_range_idx(&global_avg, iter.idx));
    for (int k=0; k<basis_avg.num_basis; ++k)
      max_err = fmax(max_err, fabs(a[k]-r[k])/fmax(fabs(r[k]), 1.0));
  }

  gkyl_array_average_release(up);
  gkyl_array_average_release(up_ref);
  gkyl_array_release(f_loc);
  gkyl_array_release(w_loc);
  gkyl_array_release(avg);
  gkyl_array_release(f_glob);
  gkyl_array_release(w_glob);
  gkyl_array_release(avg_ref);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);

  return max_err;
}

static void
get_cuts_2x(int *cuts)
{
  int sz;
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  cuts[0] = sz == 4 ? 2 : sz;
  cuts[1] = sz == 4 ? 2 : 1;
}

static void
test_2x()
{
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {16, 8};
  int cuts[2];
  get_cuts_2x(cuts);

  int avg_x[] = {1, 0}, avg_y[] = {0, 1}, avg_xy[] = {1, 1};
  for (int p=1; p<=2; ++p) {
    for (int w=0; w<2; ++w) {
      double err;
      err = run_avg(2, lower, upper, cells, cuts, avg_x, p, w, true);
      TEST_CHECK( err < 1e-12 );
      TEST_MSG("avg_x p=%d weighted=%d err=%g", p, w, err);
      err = run_avg(2, lower, upper, cells, cuts, avg_y, p, w, true);
      TEST_CHECK( err < 1e-12 );
      TEST_MSG("avg_y p=%d weighted=%d err=%g", p, w, err);
      err = run_avg(2, lower, upper, cells, cuts, avg_xy, p, w, true);
      TEST_CHECK( err < 1e-12 );
      TEST_MSG("avg_xy p=%d weighted=%d err=%g", p, w, err);
    }
  }
}

static void
test_3x()
{
  double lower[] = {-2.0, -3.0, -1.0}, upper[] = {3.0, 5.0, 2.0};
  int cells[] = {8, 8, 4};
  int sz;
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  int cuts[] = {sz == 4 ? 2 : 1, sz == 4 ? 2 : sz, 1};

  int avg_yz[] = {0, 1, 1}, avg_x[] = {1, 0, 0};
  for (int p=1; p<=2; ++p) {
    for (int w=0; w<2; ++w) {
      double err;
      err = run_avg(3, lower, upper, cells, cuts, avg_yz, p, w, true);
      TEST_CHECK( err < 1e-12 );
      TEST_MSG("avg_yz p=%d weighted=%d err=%g", p, w, err);
      err = run_avg(3, lower, upper, cells, cuts, avg_x, p, w, true);
      TEST_CHECK( err < 1e-12 );
      TEST_MSG("avg_x p=%d weighted=%d err=%g", p, w, err);
    }
  }
}

// Without communicator, the weighted average over a decomposed dimension is incorrect.
static void
test_2x_no_comm()
{
  int sz;
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  if (sz == 1) return;

  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {16, 8};
  int cuts[2];
  get_cuts_2x(cuts);

  int avg_x[] = {1, 0};
  double err = run_avg(2, lower, upper, cells, cuts, avg_x, 1, true, false);
  TEST_CHECK( err > 1e-3 );
  TEST_MSG("weighted avg_x without comm err=%g", err);
}

TEST_LIST = {
  { "test_2x", test_2x },
  { "test_3x", test_3x },
  { "test_2x_no_comm", test_2x_no_comm },
  { NULL, NULL },
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  { NULL, NULL },
};

#endif
