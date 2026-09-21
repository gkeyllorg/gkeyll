// Test the up-down tokamak core BC updater (GKYL_BC_GK_SPECIES_UPDOWN_TOK_CORE),
// which fills the lower radial ghost cells with the up-down mirror image of the
// skin cells:
//   f_ghost(x_g, y, z, vpar, mu) = < f_skin(x_s, y', -z, vpar, mu) >_y'
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_updown_tok_core.h>
#include <gkyl_null_comm.h>
#include <gkyl_proj_on_basis.h>
#include <math.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_ctx {
  int cdim; // Number of conf-space dimensions.
  bool even; // Use a function even in z and vpar.
  // Parameters of the expected (ghost) function.
  double x_face; // Location of the boundary.
  double y_lo, y_up; // y extents (y average).
};

// y-dependence of the test function (cubic, so the quadrature is exact).
static inline double
fy(double y)
{
  return 1.0 + 0.5*y + 0.3*y*y - 0.2*y*y*y;
}

// Its y average over [a,b].
static inline double
fy_avg(double a, double b)
{
  // Antiderivative of fy.
  double Fb = b + 0.25*b*b + 0.1*b*b*b - 0.05*b*b*b*b;
  double Fa = a + 0.25*a*a + 0.1*a*a*a - 0.05*a*a*a*a;
  return (Fb-Fa)/(b-a);
}

// Part of the test function that does not depend on y.
static inline double
fxzvm(bool even, double x, double z, double vpar, double mu)
{
  if (even) {
    return (1.0 + 0.3*x + 0.2*x*x) * (1.0 + 0.4*cos(z)) * exp(-0.5*vpar*vpar) * (1.0 + 0.1*mu);
  }
  else {
    return (1.0 + 0.3*x + 0.2*x*x) * (1.0 + 0.4*z + 0.2*sin(2.0*z))
      * exp(-0.5*pow(vpar-0.5,2)) * (1.0 + 0.1*mu + 0.05*mu*mu)
      + 0.5*x*z*vpar + 0.1*z*mu;
  }
}

// Test function f.
void
eval_f(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  int cdim = tc->cdim;
  double x = xn[0], z = xn[cdim-1], vpar = xn[cdim], mu = xn[cdim+1];
  double yfac = cdim == 3 ? fy(xn[1]) : 1.0;
  fout[0] = fxzvm(tc->even, x, z, vpar, mu) * yfac;
}

// Expected ghost function.
void
eval_ghost(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  int cdim = tc->cdim;
  double x = xn[0], z = xn[cdim-1], vpar = xn[cdim], mu = xn[cdim+1];

  // Donor x location (reflected about the boundary).
  double xd = 2.0*tc->x_face - x;
  double yfac = cdim == 3 ? fy_avg(tc->y_lo, tc->y_up) : 1.0;
  fout[0] = fxzvm(tc->even, xd, -z, vpar, mu) * yfac;
}

void
test_bc_updown_tok_core(int cdim, bool even, bool use_gpu)
{
  enum gkyl_edge_loc edge = GKYL_LOWER_EDGE;
  int vdim = 2;
  int ndim = cdim+vdim;
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];

  // x.
  lower[0] = 0.5, upper[0] = 1.5, cells[0] = 4;
  // y.
  if (cdim == 3) {
    lower[1] = 0.0, upper[1] = 1.0, cells[1] = 4;
  }
  // z.
  lower[cdim-1] = -M_PI, upper[cdim-1] = M_PI, cells[cdim-1] = 6;
  // vpar, mu.
  lower[cdim] = -3.0, upper[cdim] = 3.0, cells[cdim] = 6;
  lower[cdim+1] = 0.0, upper[cdim+1] = 2.0, cells[cdim+1] = 3;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis *basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                                    : gkyl_cart_modal_gkhybrid_new(cdim, vdim);
  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int dir = 0;
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  struct test_ctx tc = {
    .cdim = cdim,
    .even = even,
    .x_face = lower[0],
    .y_lo = cdim == 3 ? lower[1] : 0.0,
    .y_up = cdim == 3 ? upper[1] : 1.0,
  };

  // Project f in the interior.
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  gkyl_array_clear(distf_ho, 0.0);
  gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis_ho, 2, 1, eval_f, &tc);
  gkyl_proj_on_basis_advance(proj_f, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Keep a copy to check that the interior is untouched.
  struct gkyl_array *distf0_ho = mkarr(false, distf->ncomp, distf->size);
  gkyl_array_copy(distf0_ho, distf_ho);

  // Expected ghost.
  struct gkyl_array *distf_ex = mkarr(false, basis_ho.num_basis, local_ext.volume);
  gkyl_array_clear(distf_ex, 0.0);
  gkyl_proj_on_basis *proj_g = gkyl_proj_on_basis_new(&grid, &basis_ho, 2, 1, eval_ghost, &tc);
  gkyl_proj_on_basis_advance(proj_g, 0.0, &ghost_r, distf_ex);

  // Serial: a single-rank conf-space decomposition and a null comm.
  struct gkyl_range global_conf;
  gkyl_create_global_range(cdim, cells, &global_conf);
  int cuts[GKYL_MAX_CDIM] = {1, 1, 1};
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &global_conf);
  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = use_gpu
    }
  );

  // Create and apply the updater.
  struct gkyl_bc_updown_tok_core_inp inp = {
    .dir = dir,
    .edge = edge,
    .cdim = cdim,
    .grid = &grid,
    .basis = basis,
    .skin_r = &skin_r,
    .ghost_r = &ghost_r,
    .decomp = decomp,
    .comm = comm,
    .avg_y = true,
    .use_gpu = use_gpu,
  };
  struct gkyl_bc_updown_tok_core *up = gkyl_bc_updown_tok_core_new(&inp);

  gkyl_bc_updown_tok_core_advance(up, distf);

  gkyl_array_copy(distf_ho, distf);

  // Check the ghost against the expected function.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ghost_r, iter.idx);
    const double *f_p = gkyl_array_cfetch(distf_ho, loc);
    const double *fex_p = gkyl_array_cfetch(distf_ex, loc);
    // Tolerance relative to the cell average (quadrature roundoff).
    double tol = 1e-12*(1.0+fabs(fex_p[0]));
    for (int k=0; k<basis_ho.num_basis; k++) {
      TEST_CHECK( fabs(fex_p[k] - f_p[k]) < tol );
      TEST_MSG("idx=(%d,%d,%d,%d,%d) k=%d | Expected: %.14e | Got: %.14e", iter.idx[0], iter.idx[1],
        iter.idx[2], iter.idx[3], ndim>4? iter.idx[4] : 0, k, fex_p[k], f_p[k]);
    }

    // Ghost must be independent of y.
    if (cdim == 3) {
      int idx0[GKYL_MAX_DIM];
      for (int d=0; d<ndim; d++) idx0[d] = iter.idx[d];
      idx0[1] = ghost_r.lower[1];
      const double *f0_p = gkyl_array_cfetch(distf_ho, gkyl_range_idx(&ghost_r, idx0));
      for (int k=0; k<basis_ho.num_basis; k++)
        TEST_CHECK( f0_p[k] == f_p[k] );
    }
  }

  // Check the interior (incl. the skin) is untouched.
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *f_p = gkyl_array_cfetch(distf_ho, loc);
    const double *f0_p = gkyl_array_cfetch(distf0_ho, loc);
    for (int k=0; k<basis_ho.num_basis; k++)
      TEST_CHECK( f0_p[k] == f_p[k] );
  }

  gkyl_bc_updown_tok_core_release(up);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_proj_on_basis_release(proj_f);
  gkyl_proj_on_basis_release(proj_g);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf0_ho);
  gkyl_array_release(distf_ex);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);
}

void
test_all(int cdim, bool use_gpu)
{
  // Generic f.
  test_bc_updown_tok_core(cdim, false, use_gpu);
  // f even in z and vpar.
  test_bc_updown_tok_core(cdim, true, use_gpu);
}

void test_2x2v_ho() { test_all(2, false); }
void test_3x2v_ho() { test_all(3, false); }

#ifdef GKYL_HAVE_CUDA
void test_2x2v_dev() { test_all(2, true); }
void test_3x2v_dev() { test_all(3, true); }
#endif

TEST_LIST = {
  { "test_2x2v_ho", test_2x2v_ho },
  { "test_3x2v_ho", test_3x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_2x2v_dev", test_2x2v_dev },
  { "test_3x2v_dev", test_3x2v_dev },
#endif
  { NULL, NULL },
};
