// Test the perpendicular FEM Helmholtz/Poisson solver.
// The following cases are omitted:
// 2D Px–Px, p=1 — solver runs, solution not verified
// 2D any BC, p=2 — no code exists for it
// 3D Px–Px / Py–Py, p=1 — solver runs, check block is a no-op comment
// 3D p=2, all BCs — check code exists but tests are disabled in TEST_LIST
//
#include <acutest.h>
#include <assert.h>
#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
//#include <gkyl_fem_parproj.h>

static double error_L2norm(struct gkyl_rect_grid grid, struct gkyl_range range,
  struct gkyl_basis basis, struct gkyl_array* field1, struct gkyl_array* field2)
{
  // Compute the L2 norm of the difference between 2 fields.
  assert(field1->ncomp == field2->ncomp);
  assert(field1->size == field2->size);

  struct gkyl_array *diff = gkyl_array_new(GKYL_DOUBLE, field1->ncomp, field1->size);
  gkyl_array_copy(diff, field1);
  gkyl_array_accumulate(diff, -1.0, field2);

  struct gkyl_array *l2_cell = gkyl_array_new(GKYL_DOUBLE, 1, field1->size);
  gkyl_dg_calc_l2_range(&basis, 0, l2_cell, 0, diff, range);
  gkyl_array_scale_range(l2_cell, grid.cellVolume, &range);

  double l2[1];
  gkyl_array_reduce_range(l2, l2_cell, GKYL_SUM, &range);

  gkyl_array_release(diff);
  gkyl_array_release(l2_cell);
  return sqrt(l2[0]);
}

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // allocate array (filled with zeros)
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

double poly_test_func_1x(double x, double a, double *c)
{
  // Function that can be used to produce homogeneous Dirichlet or Neumann
  // boundary values depending on the choice of a and c. It assumes x \in [0,1].
  return pow(x,2)/2.-a*pow(x,4)/12.+c[0]*x+c[1];
}

void evalFunc_consteps_periodicx_sol_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  // These values have to match those in the test below.
  double gxx = 1.0;
  double amn[] = {0., 10., 10.};
  double bmn[] = {0., 10., 10.};
  fout[0] = 0.;
  for (int m=1; m<4; m++) {
    double a = amn[m-1];
    double b = bmn[m-1];
    double t1 = a*cos(m*x);
    double t2 = b*sin(m*x);
    fout[0] += t1+t2;
  }
  double kz = 1.;
  fout[0] *= (1.+kz*z);
//  fout[0] *= (1.+kz*z+0.5*pow(z,2));
}
void evalFunc_consteps_periodicx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  // These values have to match those in the test below.
  double gxx = 1.0;
  double amn[] = {0., 10., 10.};
  double bmn[] = {0., 10., 10.};
  fout[0] = 0.;
  for (int m=1; m<4; m++) {
    double a = amn[m-1];
    double b = bmn[m-1];
    double t1 = a*gxx*pow(m,2)*cos(m*x);
    double t2 = b*gxx*pow(m,2)*sin(m*x);
    fout[0] += t1+t2;
  }
  double kz = 1.;
  fout[0] *= (1.+kz*z);
//  fout[0] *= (1.+kz*z+0.5*pow(z,2));
}

void evalFunc_consteps_dirichletx_sol_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double kz = 1.;
  double xp = x, zp = z;
  fout[0] = poly_test_func_1x(xp, a, c)
//           *sin(kz*z);
           *(1.+kz*z);
//           *(1.+kz*z+0.5*pow(z,2));
}
void evalFunc_consteps_dirichletx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double kz = 1.;
  double xp = x, zp = z;
  fout[0] = -(1 - a*pow(xp,2))
//              *sin(kz*z);
              *(1.+kz*z);
//              *(1.+kz*z+0.5*pow(z,2));
}

void evalFunc_consteps_neumannx_dirichletx_sol_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  fout[0] = poly_test_func_1x(x, a, c);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_neumannx_dirichletx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  fout[0] = -(1 - a*pow(x,2));

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_dirichletx_neumannx_sol_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  fout[0] = poly_test_func_1x(x-1., a, c);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_dirichletx_neumannx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  fout[0] = -(1 - a*pow(x-1.,2));

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

// Persistent objects created by consteps test.
struct fem_poisson_perp_consteps_objs {
  struct gkyl_fem_poisson_perp *poisson;
  struct gkyl_array *eps;
  struct gkyl_array *rho;
  struct gkyl_array *phi;
  struct gkyl_array *phisol_ho;
  struct gkyl_array *rho_ho;
  struct gkyl_array *phi_ho;
  struct gkyl_range localRange, localRange_ext;
  struct gkyl_basis basis;
};

void
fem_poisson_perp_consteps_objs_release(struct fem_poisson_perp_consteps_objs *objs)
{
  gkyl_fem_poisson_perp_release(objs->poisson);
  gkyl_array_release(objs->eps);
  gkyl_array_release(objs->rho);
  gkyl_array_release(objs->phi);
  gkyl_array_release(objs->phisol_ho);
  gkyl_array_release(objs->rho_ho);
  gkyl_array_release(objs->phi_ho);
  gkyl_free(objs);
}


// Helper: verify phi_ho against reference solution sol over a range.
// pass filter_dim < 0 to check all cells; otherwise only checks cells
// where iter.idx[filter_dim] == filter_idx.
// Solution; checked convergence but note that p=2 serendipity doesn't
// converge as p+1. One must use tensor basis and a RHS in the space of
// the basis for that.
static void
check_sol(const struct gkyl_range *range, const struct gkyl_array *phi_ho,
  int num_basis, const double *sol, double scale_fac,
  int filter_dim, int filter_idx, double tol)
{
  long i = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    if (filter_dim >= 0 && iter.idx[filter_dim] != filter_idx) continue;
    long loc = gkyl_range_idx(range, iter.idx);
    const double *phi_p = gkyl_array_cfetch(phi_ho, loc);
    for (int m = 0; m < num_basis; m++) {
      TEST_CHECK( gkyl_compare(sol[i]*scale_fac, phi_p[m], tol) );
      TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[i]*scale_fac,
        iter.idx[0], iter.idx[1], iter.idx[2]);
      TEST_MSG("Produced: %.13e", phi_p[m]);
      i += 1;
    }
  }
}
// ============================================================
// Accessors for solution arrays (defined at end of file).
// Naming: get_sol_<ndim>x_[bias[B]_]p<order>_<bcx>[_<bcy>]
//   D=Dirichlet  N=Neumann  P=Periodic
// ============================================================

static const double *get_sol_2x_p1_DxDx(void);
static const double *get_sol_2x_p1_NxDx(void);
static const double *get_sol_2x_p1_DxNx(void);
static const double *get_sol_2x_bias_p1_DxDx(void);
static const double *get_sol_2x_bias_p1_NxDx(void);
static const double *get_sol_2x_bias_p1_DxNx(void);
static const double *get_sol_3x_bias_p1_DxDx_DyDy(void);
static const double *get_sol_3x_bias_p1_DxDx_PyPy(void);
static const double *get_sol_3x_bias_p1_PxPx_DyDy(void);
static const double *get_sol_3x_bias_p1_DxDx_NyDy(void);
static const double *get_sol_3x_bias_p1_DxDx_DyNy(void);
static const double *get_sol_3x_bias_p1_NxDx_DyDy(void);
static const double *get_sol_3x_bias_p1_DxNx_DyDy(void);
static const double *get_sol_3x_bias_p1_NxDx_PyPy(void);
static const double *get_sol_3x_bias_p2_PxPx_PyPy(void);
static const double *get_sol_3x_bias_p2_DxDx_DyDy(void);
static const double *get_sol_3x_bias_p2_DxDx_PyPy(void);
static const double *get_sol_3x_bias_p2_PxPx_DyDy(void);
static const double *get_sol_3x_bias_p2_DxDx_NyDy(void);
static const double *get_sol_3x_bias_p2_DxDx_DyNy(void);
static const double *get_sol_3x_bias_p2_NxDx_DyDy(void);
static const double *get_sol_3x_bias_p2_DxNx_DyDy(void);
static const double *get_sol_3x_bias_p2_NxDx_PyPy(void);
static const double *get_sol_3x_biasB_p1_DxDx_DyDy(void);
static const double *get_sol_3x_biasB_p1_DxDx_PyPy(void);

void
fem_poisson_perp_consteps_2x_check(struct fem_poisson_perp_consteps_objs *objs,
  int poly_order, struct gkyl_poisson_bc bcs, double scale_fac)
{
  // Check results in fem_poisson_perp_consteps_2x test.

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution not yet verified for periodic 2x.
    } else if (bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) {
      check_sol(&objs->localRange, objs->phi_ho, objs->basis.num_basis,
        get_sol_2x_p1_DxDx(), scale_fac, -1, -1, 1e-10);
    } else if (bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) {
      check_sol(&objs->localRange, objs->phi_ho, objs->basis.num_basis,
        get_sol_2x_p1_NxDx(), scale_fac, -1, -1, 1e-10);
    } else if (bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_NEUMANN) {
      check_sol(&objs->localRange, objs->phi_ho, objs->basis.num_basis,
        get_sol_2x_p1_DxNx(), scale_fac, -1, -1, 1e-10);
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
      TEST_MSG("This BC combination is not available");
    }
  } else {
    TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
    TEST_MSG("This poly_order is not available");
  }

}

struct fem_poisson_perp_consteps_objs*
test_fem_poisson_perp_consteps_2x_objs(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  if (   (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET)
      || (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET)
      || (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) )
  {
    lower[0] = 0.;  upper[0] = 1.;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);
  int dim_perp = dim-1; 

  struct fem_poisson_perp_consteps_objs *objs = gkyl_malloc(sizeof(*objs));

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  gkyl_cart_modal_serendip(&objs->basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, ghost, &objs->localRange_ext, &objs->localRange);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if (bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &objs->basis,
      poly_order+1, 1, evalFunc_consteps_periodicx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &objs->basis,
      2*(poly_order+1), 1, evalFunc_consteps_periodicx_sol_2x, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &objs->basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &objs->basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_sol_2x, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &objs->basis,
      poly_order+1, 1, evalFunc_consteps_neumannx_dirichletx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &objs->basis,
      2*(poly_order+1), 1, evalFunc_consteps_neumannx_dirichletx_sol_2x, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) {
    projob = gkyl_proj_on_basis_new(&grid, &objs->basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_neumannx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &objs->basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_neumannx_sol_2x, NULL);
  }

  // Create DG field we wish to make continuous.
  objs->rho = mkarr(use_gpu, objs->basis.num_basis, objs->localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  objs->phi = mkarr(use_gpu, objs->basis.num_basis, objs->localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim_perp+ceil((pow(3.,dim_perp-1)-dim_perp)/2);
  objs->eps = mkarr(use_gpu, epsnum*objs->basis.num_basis, objs->localRange_ext.volume);
  // Analytic solution.
  objs->phisol_ho = mkarr(false, objs->basis.num_basis, objs->localRange_ext.volume);
  // Device copies:
  if (use_gpu) {
    objs->rho_ho = mkarr(false, objs->rho->ncomp, objs->rho->size);
    objs->phi_ho = mkarr(false, objs->phi->ncomp, objs->phi->size);
  }
  else {
    objs->rho_ho = gkyl_array_acquire(objs->rho);
    objs->phi_ho = gkyl_array_acquire(objs->phi);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &objs->localRange, objs->rho_ho);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_copy(objs->rho, objs->rho_ho);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(objs->eps, epsilon_0*dg0norm, 0*objs->basis.num_basis);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &objs->localRange, objs->phisol_ho);
  gkyl_proj_on_basis_release(projob_sol);

  // FEM poisson solver.
  objs->poisson = gkyl_fem_poisson_perp_new(&objs->localRange, &grid, objs->basis, &bcs, NULL, objs->eps, NULL, use_gpu);

//  struct gkyl_fem_parproj* smooth_op = gkyl_fem_parproj_new(&objs->localRange, &grid,
//    &objs->localRange_ext, &objs->basis, GKYL_FEM_PARPROJ_DIRICHLET, NULL, NULL, NULL, use_gpu);
//  gkyl_fem_parproj_set_rhs(smooth_op, objs->rho, objs->rho);
//  gkyl_fem_parproj_solve  (smooth_op, objs->rho);

  // Check that the solver can be called multiple times without issue.
  for (int i=0; i<2; i++) {
    // Set the RHS source.
    gkyl_fem_poisson_perp_set_rhs(objs->poisson, objs->rho);

    // Solve the problem.
    gkyl_fem_poisson_perp_solve(objs->poisson, objs->phi);
  }
  gkyl_array_copy(objs->phi_ho, objs->phi);
//  gkyl_fem_parproj_set_rhs(smooth_op, phi, phi);
//  gkyl_fem_parproj_solve  (smooth_op, phi);
//  gkyl_fem_parproj_release(smooth_op);

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    // Subtract the volume averaged sol from the numerical and analytic solutions.
    // This is not strictly necessary, as the potential is only known up to 
    // constant shift, but it makes unit testing more robust across CPU/GPU.
    struct gkyl_array *sol_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, objs->localRange_ext.volume);
    double sol_avg[1];
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim); // /perpRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_array_clear(sol_cellavg, 0.0);
    gkyl_dg_calc_average_range(&objs->basis, 0, sol_cellavg, 0, objs->phi_ho, objs->localRange);
    for (int kIdx=0; kIdx<cells[1]; kIdx++) {
      struct gkyl_range perp_range;
      gkyl_range_deflate(&perp_range, &objs->localRange, (int[]){0,0,1}, (int[]){0,0,kIdx+1});
      gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &perp_range);
      gkyl_array_shiftc_range(objs->phi_ho, mavgfac*sol_avg[0]/perp_range.volume, 0, &perp_range);
    }
    // Now do the same to the analytic solution.
    gkyl_array_clear(sol_cellavg, 0.0);
    gkyl_dg_calc_average_range(&objs->basis, 0, sol_cellavg, 0, objs->phisol_ho, objs->localRange);
    for (int kIdx=0; kIdx<cells[1]; kIdx++) {
      struct gkyl_range perp_range;
      gkyl_range_deflate(&perp_range, &objs->localRange, (int[]){0,0,1}, (int[]){0,0,kIdx+1});
      gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &perp_range);
      gkyl_array_shiftc_range(objs->phisol_ho, mavgfac*sol_avg[0]/perp_range.volume, 0, &perp_range);
    }
    gkyl_array_release(sol_cellavg);
  }

//  double errL2 = error_L2norm(grid, objs->localRange, objs->basis, phi, phisol);
//  printf("error L2 norm = %g\n",errL2);

  // Check results for correctness.
  fem_poisson_perp_consteps_2x_check(objs, poly_order, bcs, 1.0);

//  gkyl_grid_sub_array_write(&grid, &localRange, 0, objs->rho_ho, "ctest_fem_poisson_perp_2x_rho_1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, objs->phi_ho, "ctest_fem_poisson_perp_2x_phi_8x8_p1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, objs->phisol_ho, "ctest_fem_poisson_perp_2x_phisol_8x8_p1.gkyl");

  return objs;
}

void
test_fem_poisson_perp_consteps_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  struct fem_poisson_perp_consteps_objs *objs = test_fem_poisson_perp_consteps_2x_objs(poly_order, cells, bcs, use_gpu);
  fem_poisson_perp_consteps_objs_release(objs);
}

void
test_fem_poisson_perp_consteps_2x_update(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  // Run the first test.
  struct fem_poisson_perp_consteps_objs *objs = test_fem_poisson_perp_consteps_2x_objs(poly_order, cells, bcs, use_gpu);

  double global_fac = 1.0;
  // Check that the solver can be called multiple times without issue.
  for (int i=0; i<2; i++) {
    // Now update the LHS matrix. Multiply it by a constant so the solution should be the same as before but divided by that constant.
    double prob_fac = 1.3;
    gkyl_array_scale(objs->eps, prob_fac);
    gkyl_fem_poisson_perp_update_lhs(objs->poisson, objs->eps, NULL);

    // Set the RHS source.
    gkyl_fem_poisson_perp_set_rhs(objs->poisson, objs->rho);

    // Solve the problem.
    gkyl_fem_poisson_perp_solve(objs->poisson, objs->phi);

    global_fac *= prob_fac;
  }
  gkyl_array_copy(objs->phi_ho, objs->phi);

  // Check results for correctness.
  fem_poisson_perp_consteps_2x_check(objs, poly_order, bcs, 1.0/global_fac);

  // Release persistent objects. 
  fem_poisson_perp_consteps_objs_release(objs);
}

void
test_fem_poisson_perp_consteps_2x_bias(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  if (   (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET)
      || (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET)
      || (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) )
  {
    lower[0] = 0.;  upper[0] = 1.;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);
  int dim_perp = dim-1; 

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_sol_2x, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_neumannx_dirichletx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_neumannx_dirichletx_sol_2x, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_neumannx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_neumannx_sol_2x, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim_perp+ceil((pow(3.,dim_perp-1)-dim_perp)/2);
  struct gkyl_array *eps = mkarr(use_gpu, epsnum*basis.num_basis, localRange_ext.volume);
  // Analytic solution.
  struct gkyl_array *phisol_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, rho->ncomp, rho->size);
    phi_ho = mkarr(false, phi->ncomp, phi->size);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 0*basis.num_basis);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol_ho);

  // Specify the bias:
  struct gkyl_poisson_bias_line bias = {
    .perp_dirs = {0, 1},
    .perp_coords = {0.5, -M_PI}, // Location of the plane in the 'dir' dimension.
    .val = 0., // Biasing value.
  };
  struct gkyl_poisson_bias_line_list bll = {
    .num_bias_line = 1, // Number of bias lines.
    .bl = &bias,
  };

  // FEM poisson solver.
  struct gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&localRange, &grid, basis, &bcs, &bll, eps, NULL, use_gpu);

  // Check that the solver can be called multiple times without issue.
  for (int i=0; i<2; i++) {
    // Set the RHS source.
    gkyl_fem_poisson_perp_set_rhs(poisson, rho);

    // Solve the problem.
    gkyl_fem_poisson_perp_solve(poisson, phi);

    gkyl_array_copy(phi_ho, phi);

  //  double errL2 = error_L2norm(grid, localRange, basis, phi, phisol);
  //  printf("error L2 norm = %g\n",errL2);

    if (poly_order == 1) {
      if (bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_2x_bias_p1_DxDx(), 1.0, -1, -1, 1e-10);
      } else if (bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_2x_bias_p1_NxDx(), 1.0, -1, -1, 1e-10);
      } else if (bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_NEUMANN) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_2x_bias_p1_DxNx(), 1.0, -1, -1, 1e-10);
      } else {
        TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
        TEST_MSG("This BC combination is not available");
      }
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
      TEST_MSG("This poly_order is not available");
    }
  }

//  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_poisson_perp_2x_rho_1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_poisson_perp_2x_phi_8x8_p1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phisol_ho, "ctest_fem_poisson_perp_2x_phisol_8x8_p1.gkyl");

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  gkyl_array_release(phisol_ho);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
}

void evalFunc_consteps_periodicx_periodicy_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  // These values have to match those in the test below.
  double gxx = 1.0,  gxy = 0.,  gyy = 1.0;
  double amn[] = {0., 10., 0., 10., 0., 0., 10., 0., 0.};
  double bmn[] = {0., 10., 0., 10., 0., 0., 10., 0., 0.};
  fout[0] = 0.;
  for (int m=1; m<4; m++) {
    for (int n=1; n<4; n++) {
      double a = amn[(m-1)*3+(n-1)];
      double b = bmn[(m-1)*3+(n-1)];
      double t1 = a*cos(m*x)*cos(n*y);
      double t2 = b*sin(m*x)*sin(n*y);
      fout[0] += t1+t2;
    }
  }
  double kz = 1.;
  fout[0] *= (1.+kz*z);
//  fout[0] *= (1.+kz*z+0.5*pow(z,2));
}
void evalFunc_consteps_periodicx_periodicy_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  // These values have to match those in the test below.
  double gxx = 1.0,  gxy = 0.,  gyy = 1.0;
  double amn[] = {0., 10., 0., 10., 0., 0., 10., 0., 0.};
  double bmn[] = {0., 10., 0., 10., 0., 0., 10., 0., 0.};
  fout[0] = 0.;
  for (int m=1; m<4; m++) {
    for (int n=1; n<4; n++) {
      double a = amn[(m-1)*3+(n-1)];
      double b = bmn[(m-1)*3+(n-1)];
      double t1 = (a*gxx*pow(m,2) - 2*b*gxy*m*n + a*gyy*pow(n,2))*cos(m*x)*cos(n*y);
      double t2 = (b*gxx*pow(m,2) - 2*a*gxy*m*n + b*gyy*pow(n,2))*sin(m*x)*sin(n*y);
      fout[0] += t1+t2;
    }
  }
  double kz = 1.;
  fout[0] *= (1.+kz*z);
//  fout[0] *= (1.+kz*z+0.5*pow(z,2));
}

void evalFunc_consteps_dirichletx_dirichlety_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  double kz = 1.;
  double xp = x, yp = y, zp = z;
  fout[0] = poly_test_func_1x(xp, a, c)*poly_test_func_1x(yp, b, d)
//           *sin(kz*z);
           *(1.+kz*z);
//           *(1.+kz*z+0.5*pow(z,2));
}
void evalFunc_consteps_dirichletx_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  double kz = 1.;
  double xp = x, yp = y, zp = z;
  fout[0] = -( (1 - b*pow(yp,2))*poly_test_func_1x(xp, a, c)
              +(1 - a*pow(xp,2))*poly_test_func_1x(yp, b, d) )
//              *sin(kz*z);
              *(1.+kz*z);
//              *(1.+kz*z+0.5*pow(z,2));
}

void evalFunc_consteps_dirichletx_periodicy_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double n = 2.;
  fout[0] = poly_test_func_1x(x, a, c)*sin(n*y);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_dirichletx_periodicy_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double n = 2.;
  fout[0] = -( (1-a*pow(x,2))*sin(n*y)
              -pow(n,2)*poly_test_func_1x(x, a, c)*sin(n*y) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_periodicx_dirichlety_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  double m = 2.;
  fout[0] = sin(m*x)*poly_test_func_1x(y, b, d);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_periodicx_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  double m = 2.;
  fout[0] = -( (1-b*pow(y,2))*sin(m*x)
              -pow(m,2)*poly_test_func_1x(y, b, d)*sin(m*x) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_dirichletx_neumanny_dirichlety_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 5.;
  double d[] = {0., b/12.-1./2.};
  fout[0] = poly_test_func_1x(x, a, c)*poly_test_func_1x(y, b, d);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_dirichletx_neumanny_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 5.;
  double d[] = {0., b/12.-1./2.};
  fout[0] = -( (1 - b*pow(y,2))*poly_test_func_1x(x, a, c)
              +(1 - a*pow(x,2))*poly_test_func_1x(y, b, d) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_dirichletx_dirichlety_neumanny_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 5.;
  double d[] = {0., b/12.-1./2.};
  fout[0] = poly_test_func_1x(x, a, c)*poly_test_func_1x(y-1., b, d);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_dirichletx_dirichlety_neumanny_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 5.;
  double d[] = {0., b/12.-1./2.};
  fout[0] = -( (1 - b*pow(y-1.,2))*poly_test_func_1x(x   , a, c)
              +(1 - a*pow(x   ,2))*poly_test_func_1x(y-1., b, d) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_neumannx_dirichletx_dirichlety_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  fout[0] = poly_test_func_1x(x, a, c)*poly_test_func_1x(y, b, d);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_neumannx_dirichletx_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  fout[0] = -( (1 - b*pow(y,2))*poly_test_func_1x(x, a, c)
              +(1 - a*pow(x,2))*poly_test_func_1x(y, b, d) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_dirichletx_neumannx_dirichlety_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  fout[0] = poly_test_func_1x(x-1., a, c)*poly_test_func_1x(y, b, d);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_dirichletx_neumannx_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  fout[0] = -( (1 - b*pow(y   ,2))*poly_test_func_1x(x-1., a, c)
              +(1 - a*pow(x-1.,2))*poly_test_func_1x(y   , b, d) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void evalFunc_consteps_neumannx_dirichletx_periodicy_sol_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  double n = 2.;
  fout[0] = poly_test_func_1x(x, a, c)*sin(n*y);

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}
void evalFunc_consteps_neumannx_dirichletx_periodicy_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 5.;
  double c[] = {0., a/12.-1./2.};
  double n = 2.;
  fout[0] = -( (1 - a*pow(x,2))*sin(n*y)
              -pow(n,2)*sin(n*y)*poly_test_func_1x(x, a, c) );

  double kz = 1.;
  fout[0] *= (1.+kz*z);
}

void
test_fem_poisson_perp_consteps_3x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI,-M_PI}, upper[] = {M_PI,M_PI,M_PI};
  if (   ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) 
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) 
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_NEUMANN))
      || ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) )
  {
    lower[0] = 0.;  upper[0] = 1.;
    lower[1] = 0.;  upper[1] = 1.;
  } else if ( ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) 
           || ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) )
  {
    lower[0] = 0.;  upper[0] = 1.;
  } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    lower[1] = 0.;  upper[1] = 1.;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);
  int dim_perp = dim-1; 

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
      (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_periodicx_periodicy_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_periodicx_periodicy_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_dirichlety_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_periodicy_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_periodicy_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_periodicx_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_periodicx_dirichlety_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_neumanny_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_neumanny_dirichlety_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_NEUMANN)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_dirichlety_neumanny_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_dirichlety_neumanny_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_neumannx_dirichletx_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_neumannx_dirichletx_dirichlety_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_neumannx_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_neumannx_dirichlety_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_neumannx_dirichletx_periodicy_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_neumannx_dirichletx_periodicy_sol_3x, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim_perp+ceil((pow(3.,dim_perp-1)-dim_perp)/2);
  struct gkyl_array *eps = mkarr(use_gpu, epsnum*basis.num_basis, localRange_ext.volume);
  // Analytic solution.
  struct gkyl_array *phisol_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, rho->ncomp, rho->size);
    phi_ho = mkarr(false, phi->ncomp, phi->size);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps,        0.*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 2*basis.num_basis);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol_ho);

  // FEM poisson solver.
  struct gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&localRange, &grid, basis, &bcs, NULL, eps, NULL, use_gpu);

//  struct gkyl_fem_parproj* smooth_op = gkyl_fem_parproj_new(&localRange, &grid, &localRange_ext,
//    &basis, GKYL_FEM_PARPROJ_DIRICHLET, NULL, NULL, NULL, use_gpu);
//  gkyl_fem_parproj_set_rhs(smooth_op, rho, rho);
//  gkyl_fem_parproj_solve  (smooth_op, rho);

  // Check that the solver can be called multiple times without issue.
  for (int i=0; i<2; i++) {
    // Set the RHS source.
    gkyl_fem_poisson_perp_set_rhs(poisson, rho);

    // Solve the problem.
    gkyl_fem_poisson_perp_solve(poisson, phi);
    gkyl_array_copy(phi_ho, phi);

  //  gkyl_fem_parproj_set_rhs(smooth_op, phi, phi);
  //  gkyl_fem_parproj_solve  (smooth_op, phi);
  //  gkyl_fem_parproj_release(smooth_op);

    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Subtract the volume averaged sol from the numerical and analytic solutions.
      // This is not strictly necessary, as the potential is only known up to 
      // constant shift, but it makes unit testing more robust across CPU/GPU.
      struct gkyl_array *sol_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, localRange_ext.volume);
      double sol_avg[1];
      // Factor accounting for normalization when subtracting a constant from a
      // DG field and the 1/N to properly compute the volume averaged RHS.
      double mavgfac = -pow(sqrt(2.),dim); // /perpRange.volume;
      // Subtract the volume averaged sol from the sol.
      gkyl_array_clear(sol_cellavg, 0.0);
      gkyl_dg_calc_average_range(&basis, 0, sol_cellavg, 0, phi_ho, localRange);
      for (int kIdx=0; kIdx<cells[2]; kIdx++) {
        struct gkyl_range perp_range;
        gkyl_range_deflate(&perp_range, &localRange, (int[]){0,0,1}, (int[]){0,0,kIdx+1});
        gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &perp_range);
        gkyl_array_shiftc_range(phi_ho, mavgfac*sol_avg[0]/perp_range.volume, 0, &perp_range);
      }
      // Now do the same to the analytic solution.
      gkyl_array_clear(sol_cellavg, 0.0);
      gkyl_dg_calc_average_range(&basis, 0, sol_cellavg, 0, phisol_ho, localRange);
      for (int kIdx=0; kIdx<cells[2]; kIdx++) {
        struct gkyl_range perp_range;
        gkyl_range_deflate(&perp_range, &localRange, (int[]){0,0,1}, (int[]){0,0,kIdx+1});
        gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &perp_range);
        gkyl_array_shiftc_range(phisol_ho, mavgfac*sol_avg[0]/perp_range.volume, 0, &perp_range);
      }
      gkyl_array_release(sol_cellavg);
    }

  //  double errL2 = error_L2norm(grid, localRange, basis, phi, phisol);
  //  printf("error L2 norm = %g\n",errL2);

    if (poly_order == 1) {
      if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
          (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
        // check_sol(&localRange, phi_ho, basis.num_basis,
          // get_sol_3x_bias_p1_PxPx_PyPy(), 1.0, 2, 3, 1e-10);      
        } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_DxDx_PyPy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_PxPx_DyDy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_NEUMANN && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_DxDx_NyDy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_NEUMANN)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_DxDx_DyNy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_NxDx_DyDy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_NEUMANN) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_DxNx_DyDy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_NxDx_PyPy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis,
          get_sol_3x_bias_p1_DxDx_DyDy(), 1.0, 2, 3, 1e-10);
      } else {
        TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
        TEST_MSG("This BC combination is not available");
      }
    } else if (poly_order == 2) {
      if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
          (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_PxPx_PyPy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_DxDx_DyDy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_DxDx_PyPy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_PxPx_DyDy(), 1.0, 2, 3, 1e-10);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_NEUMANN && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_DxDx_NyDy(), 1.0, 2, 3, 5e-9);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_NEUMANN)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_DxDx_DyNy(), 1.0, 2, 3, 5e-9);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_NxDx_DyDy(), 1.0, 2, 3, 5e-9);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_NEUMANN) &&
                (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_DxNx_DyDy(), 1.0, 2, 3, 5e-9);
      } else if ((bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
                (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
        check_sol(&localRange, phi_ho, basis.num_basis/2,
          get_sol_3x_bias_p2_NxDx_PyPy(), 1.0, 2, 3, 5e-9);
      } else {
        TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
        TEST_MSG("This BC combination is not available");
      }
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
      TEST_MSG("This poly_order is not available");
    }
  }

//  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_poisson_perp_3x_rho_1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_poisson_perp_3x_phi_8x8x8_p1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phisol_ho, "ctest_fem_poisson_perp_3x_phisol_8x8x8_p1.gkyl");

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  gkyl_array_release(phisol_ho);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
}

void
test_fem_poisson_perp_consteps_3x_bias(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI,-M_PI}, upper[] = {M_PI,M_PI,M_PI};
  if (   ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) 
     )
  {
    lower[0] = 0.;  upper[0] = 1.;
    lower[1] = 0.;  upper[1] = 1.;
  } else if ( ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) 
            )
  {
    lower[0] = 0.;  upper[0] = 1.;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);
  int dim_perp = dim-1; 

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
      (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_dirichlety_sol_3x, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_consteps_dirichletx_periodicy_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_consteps_dirichletx_periodicy_sol_3x, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim_perp+ceil((pow(3.,dim_perp-1)-dim_perp)/2);
  struct gkyl_array *eps = mkarr(use_gpu, epsnum*basis.num_basis, localRange_ext.volume);
  // Analytic solution.
  struct gkyl_array *phisol_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, rho->ncomp, rho->size);
    phi_ho = mkarr(false, phi->ncomp, phi->size);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps,        0.*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 2*basis.num_basis);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol_ho);

  // Specify the bias:
  struct gkyl_poisson_bias_line bias = {
    .perp_dirs = {0, 2},
    .perp_coords = {0.5, -M_PI}, // Location of the plane in the 'dir' dimension.
    .val = 0., // Biasing value.
  };
  struct gkyl_poisson_bias_line_list bll = {
    .num_bias_line = 1, // Number of bias lines.
    .bl = &bias,
  };

  // FEM poisson solver.
  struct gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&localRange, &grid, basis, &bcs, &bll, eps, NULL, use_gpu);

//  struct gkyl_fem_parproj* smooth_op = gkyl_fem_parproj_new(&localRange, &grid, &localRange_ext,
//    &basis, GKYL_FEM_PARPROJ_DIRICHLET, NULL, NULL, NULL, use_gpu);
//  gkyl_fem_parproj_set_rhs(smooth_op, rho, rho);
//  gkyl_fem_parproj_solve  (smooth_op, rho);

  // Check that the solver can be called multiple times without issue.
  for (int i=0; i<2; i++) {
    // Set the RHS source.
    gkyl_fem_poisson_perp_set_rhs(poisson, rho);

    // Solve the problem.
    gkyl_fem_poisson_perp_solve(poisson, phi);
  }
  gkyl_array_copy(phi_ho, phi);

//  gkyl_fem_parproj_set_rhs(smooth_op, phi, phi);
//  gkyl_fem_parproj_solve  (smooth_op, phi);
//  gkyl_fem_parproj_release(smooth_op);

//  double errL2 = error_L2norm(grid, localRange, basis, phi, phisol);
//  printf("error L2 norm = %g\n",errL2);

  if (poly_order == 1) {
    if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution; checked convergence:
      check_sol(&localRange, phi_ho, basis.num_basis,
        get_sol_3x_biasB_p1_DxDx_DyDy(), 1.0, 1, 3, 1e-10);
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution; checked convergence:
      check_sol(&localRange, phi_ho, basis.num_basis,
        get_sol_3x_biasB_p1_DxDx_PyPy(), 1.0, 1, 3, 1e-10);
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
      TEST_MSG("This BC combination is not available");
    }
  } else {
    TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
    TEST_MSG("This poly_order is not available");
  }

//  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_poisson_perp_3x_rho_1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_poisson_perp_3x_phi_8x8x8_p1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phisol_ho, "ctest_fem_poisson_perp_3x_phisol_8x8x8_p1.gkyl");

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  gkyl_array_release(phisol_ho);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
}

void test_2x_p1_periodic_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, false);
}

void test_2x_p1_dirichletx_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, false);
}

void test_2x_p1_neumannx_dirichletx_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, false);
}

void test_2x_p1_dirichletx_neumannx_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, false);
}

void test_2x_p1_neumannx_dirichletx_consteps_update() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_update(1, cells, bc_tv, false);
}

void test_2x_p1_dirichletx_consteps_bias() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_bias(1, cells, bc_tv, false);
}

void test_2x_p1_neumannx_dirichletx_consteps_bias() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_bias(1, cells, bc_tv, false);
}

void test_2x_p1_dirichletx_neumannx_consteps_bias() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_bias(1, cells, bc_tv, false);
}

void test_3x_p1_periodicx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_periodicx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_neumanny_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_dirichlety_neumanny_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_neumannx_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_neumannx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_neumannx_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_dirichlety_consteps_bias() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x_bias(1, cells, bc_tv, false);
}

void test_3x_p1_dirichletx_periodicy_consteps_bias() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x_bias(1, cells, bc_tv, false);
}

void test_3x_p2_periodicx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_periodicx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_dirichletx_neumanny_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_dirichletx_dirichlety_neumanny_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_neumannx_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_dirichletx_neumannx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

void test_3x_p2_neumannx_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, false);
}

#ifdef GKYL_HAVE_CUDA
void gpu_test_2x_p1_periodic_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_dirichletx_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_neumannx_dirichletx_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_dirichletx_neumannx_consteps() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_neumannx_dirichletx_consteps_update() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_update(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_dirichletx_consteps_bias() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_bias(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_neumannx_dirichletx_consteps_bias() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_bias(1, cells, bc_tv, true);
}

void gpu_test_2x_p1_dirichletx_neumannx_consteps_bias() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_2x_bias(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_periodicx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_periodicx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_neumanny_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_dirichlety_neumanny_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_neumannx_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_neumannx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_neumannx_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_dirichlety_consteps_bias() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x_bias(1, cells, bc_tv, true);
}

void gpu_test_3x_p1_dirichletx_periodicy_consteps_bias() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x_bias(1, cells, bc_tv, true);
}

void gpu_test_3x_p2_periodicx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_periodicx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_dirichletx_neumanny_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_dirichletx_dirichlety_neumanny_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_neumannx_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_dirichletx_neumannx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

void gpu_test_3x_p2_neumannx_dirichletx_periodicy_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps_3x(2, cells, bc_tv, true);
}

#endif

TEST_LIST = {
  { "test_2x_p1_periodicx", test_2x_p1_periodic_consteps },
  { "test_2x_p1_dirichletx", test_2x_p1_dirichletx_consteps },
  { "test_2x_p1_neumannx_dirichletx", test_2x_p1_neumannx_dirichletx_consteps },
  { "test_2x_p1_dirichletx_neumannx", test_2x_p1_dirichletx_neumannx_consteps },
  { "test_2x_p1_neumannx_dirichletx_update", test_2x_p1_neumannx_dirichletx_consteps_update },
  { "test_2x_p1_dirichletx_bias", test_2x_p1_dirichletx_consteps_bias },
  { "test_2x_p1_neumannx_dirichletx_bias", test_2x_p1_neumannx_dirichletx_consteps_bias },
  { "test_2x_p1_dirichletx_neumannx_bias", test_2x_p1_dirichletx_neumannx_consteps_bias },

  { "test_3x_p1_periodicx_periodicy", test_3x_p1_periodicx_periodicy_consteps },
  { "test_3x_p1_dirichletx_dirichlety", test_3x_p1_dirichletx_dirichlety_consteps },
  { "test_3x_p1_dirichletx_periodicy", test_3x_p1_dirichletx_periodicy_consteps },
  { "test_3x_p1_periodicx_dirichlety", test_3x_p1_periodicx_dirichlety_consteps },
  { "test_3x_p1_dirichletx_neumanny_dirichlety", test_3x_p1_dirichletx_neumanny_dirichlety_consteps },
  { "test_3x_p1_dirichletx_dirichlety_neumanny", test_3x_p1_dirichletx_dirichlety_neumanny_consteps },
  { "test_3x_p1_neumannx_dirichletx_dirichlety", test_3x_p1_neumannx_dirichletx_dirichlety_consteps },
  { "test_3x_p1_dirichletx_neumannx_dirichlety", test_3x_p1_dirichletx_neumannx_dirichlety_consteps },
  { "test_3x_p1_neumannx_dirichletx_periodicy", test_3x_p1_neumannx_dirichletx_periodicy_consteps },
  { "test_3x_p1_dirichletx_dirichlety_bias", test_3x_p1_dirichletx_dirichlety_consteps_bias },
  { "test_3x_p1_dirichletx_periodicy_bias", test_3x_p1_dirichletx_periodicy_consteps_bias },
//  { "test_3x_p2_periodicx_periodicy", test_3x_p2_periodicx_periodicy_consteps },
//  { "test_3x_p2_dirichletx_dirichlety", test_3x_p2_dirichletx_dirichlety_consteps },
//  { "test_3x_p2_dirichletx_periodicy", test_3x_p2_dirichletx_periodicy_consteps },
//  { "test_3x_p2_periodicx_dirichlety", test_3x_p2_periodicx_dirichlety_consteps },
//  { "test_3x_p2_dirichletx_neumanny_dirichlety", test_3x_p2_dirichletx_neumanny_dirichlety_consteps },
//  { "test_3x_p2_dirichletx_dirichlety_neumanny", test_3x_p2_dirichletx_dirichlety_neumanny_consteps },
//  { "test_3x_p2_neumannx_dirichletx_dirichlety", test_3x_p2_neumannx_dirichletx_dirichlety_consteps },
//  { "test_3x_p2_dirichletx_neumannx_dirichlety", test_3x_p2_dirichletx_neumannx_dirichlety_consteps },
//  { "test_3x_p2_neumannx_dirichletx_periodicy", test_3x_p2_neumannx_dirichletx_periodicy_consteps },
#ifdef GKYL_HAVE_CUDA
  { "gpu_test_2x_p1_periodicx", gpu_test_2x_p1_periodic_consteps },
  { "gpu_test_2x_p1_dirichletx", gpu_test_2x_p1_dirichletx_consteps },
  { "gpu_test_2x_p1_neumannx_dirichletx", gpu_test_2x_p1_neumannx_dirichletx_consteps },
  { "gpu_test_2x_p1_dirichletx_neumannx", gpu_test_2x_p1_dirichletx_neumannx_consteps },
  { "gpu_test_2x_p1_neumannx_dirichletx_update", gpu_test_2x_p1_neumannx_dirichletx_consteps_update },
  { "gpu_test_2x_p1_dirichletx_bias", gpu_test_2x_p1_dirichletx_consteps_bias },
  { "gpu_test_2x_p1_neumannx_dirichletx_bias", gpu_test_2x_p1_neumannx_dirichletx_consteps_bias },
  { "gpu_test_2x_p1_dirichletx_neumannx_bias", gpu_test_2x_p1_dirichletx_neumannx_consteps_bias },

  { "gpu_test_3x_p1_periodicx_periodicy", gpu_test_3x_p1_periodicx_periodicy_consteps },
  { "gpu_test_3x_p1_dirichletx_dirichlety", gpu_test_3x_p1_dirichletx_dirichlety_consteps },
  { "gpu_test_3x_p1_dirichletx_periodicy", gpu_test_3x_p1_dirichletx_periodicy_consteps },
  { "gpu_test_3x_p1_periodicx_dirichlety", gpu_test_3x_p1_periodicx_dirichlety_consteps },
  { "gpu_test_3x_p1_dirichletx_neumanny_dirichlety", gpu_test_3x_p1_dirichletx_neumanny_dirichlety_consteps },
  { "gpu_test_3x_p1_dirichletx_dirichlety_neumanny", gpu_test_3x_p1_dirichletx_dirichlety_neumanny_consteps },
  { "gpu_test_3x_p1_neumannx_dirichletx_dirichlety", gpu_test_3x_p1_neumannx_dirichletx_dirichlety_consteps },
  { "gpu_test_3x_p1_dirichletx_neumannx_dirichlety", gpu_test_3x_p1_dirichletx_neumannx_dirichlety_consteps },
  { "gpu_test_3x_p1_neumannx_dirichletx_periodicy", gpu_test_3x_p1_neumannx_dirichletx_periodicy_consteps },
  { "gpu_test_3x_p1_dirichletx_dirichlety_bias", gpu_test_3x_p1_dirichletx_dirichlety_consteps_bias },
  { "gpu_test_3x_p1_dirichletx_periodicy_bias", gpu_test_3x_p1_dirichletx_periodicy_consteps_bias },
//  { "gpu_test_3x_p2_periodicx_periodicy", gpu_test_3x_p2_periodicx_periodicy_consteps },
//  { "gpu_test_3x_p2_dirichletx_dirichlety", gpu_test_3x_p2_dirichletx_dirichlety_consteps },
//  { "gpu_test_3x_p2_dirichletx_periodicy", gpu_test_3x_p2_dirichletx_periodicy_consteps },
//  { "gpu_test_3x_p2_periodicx_dirichlety", gpu_test_3x_p2_periodicx_dirichlety_consteps },
//  { "gpu_test_3x_p2_dirichletx_neumanny_dirichlety", gpu_test_3x_p2_dirichletx_neumanny_dirichlety_consteps },
//  { "gpu_test_3x_p2_dirichletx_dirichlety_neumanny", gpu_test_3x_p2_dirichletx_dirichlety_neumanny_consteps },
//  { "gpu_test_3x_p2_neumannx_dirichletx_dirichlety", gpu_test_3x_p2_neumannx_dirichletx_dirichlety_consteps },
//  { "gpu_test_3x_p2_dirichletx_neumannx_dirichlety", gpu_test_3x_p2_dirichletx_neumannx_dirichlety_consteps },
//  { "gpu_test_3x_p2_neumannx_dirichletx_periodicy", gpu_test_3x_p2_neumannx_dirichletx_periodicy_consteps },
#endif
  { NULL, NULL },
};


// ============================================================
// Solution array definitions (forward-declared above).
// ============================================================

static const double *
get_sol_2x_p1_DxDx(void)
{
  static const double sol[256] = {
         0.0592784971266791,  0.0342244562732445, -0.0076848087077185, -0.004436826376072 ,
         0.0326575388702468,  0.0188548388578077, -0.0076848087077184, -0.004436826376072 ,
         0.0060365806138147,  0.0034852214423708, -0.0076848087077184, -0.004436826376072 ,
        -0.0205843776426173, -0.011884395973066 , -0.0076848087077184, -0.004436826376072 ,
        -0.0472053358990496, -0.0272540133885028, -0.0076848087077183, -0.004436826376072 ,
        -0.0738262941554814, -0.0426236308039397, -0.0076848087077184, -0.004436826376072 ,
        -0.1004472524119136, -0.0579932482193765, -0.0076848087077186, -0.004436826376072 ,
        -0.1270682106683453, -0.0733628656348134, -0.0076848087077185, -0.004436826376072 ,
        
         0.1515053065818721,  0.0190227169922115, -0.0196410056887546, -0.0024660871694134,
         0.0834668670525276,  0.0104799404455762, -0.0196410056887545, -0.0024660871694134,
         0.0154284275231832,  0.0019371638989408, -0.0196410056887545, -0.0024660871694134,
        -0.0526100120061611, -0.0066056126476946, -0.0196410056887545, -0.0024660871694134,
        -0.1206484515355057, -0.0151483891943299, -0.0196410056887544, -0.0024660871694134,
        -0.1886868910648499, -0.0236911657409653, -0.0196410056887545, -0.0024660871694134,
        -0.2567253305941943, -0.0322339422876006, -0.0196410056887547, -0.0024660871694134,
        -0.3247637701235385, -0.040776718834236 , -0.0196410056887545, -0.0024660871694134,
        
         0.1936336022589356,  0.0053000658574412, -0.0251024783837958, -0.0006870955612404,
         0.1066760663456682,  0.0029198970139942, -0.0251024783837957, -0.0006870955612404,
         0.0197185304324009,  0.0005397281705472, -0.0251024783837957, -0.0006870955612404,
        -0.0672390054808664, -0.0018404406728998, -0.0251024783837957, -0.0006870955612404,
        -0.1541965413941337, -0.0042206095163468, -0.0251024783837957, -0.0006870955612405,
        -0.2411540773074009, -0.0066007783597938, -0.0251024783837957, -0.0006870955612404,
        -0.3281116132206683, -0.0089809472032408, -0.0251024783837958, -0.0006870955612404,
        -0.4150691491339354, -0.0113611160466878, -0.0251024783837957, -0.0006870955612404,
        
         0.192494999673069 , -0.0059574383668913, -0.0249548710136595,  0.000772316716123 ,
         0.1060487906890969, -0.0032820547831718, -0.0249548710136595,  0.000772316716123 ,
         0.0196025817051248, -0.0006066711994523, -0.0249548710136594,  0.000772316716123 ,
        -0.0668436272788473,  0.0020687123842672, -0.0249548710136594,  0.000772316716123 ,
        -0.1532898362628194,  0.0047440959679867, -0.0249548710136594,  0.000772316716123 ,
        -0.2397360452467914,  0.0074194795517062, -0.0249548710136595,  0.000772316716123 ,
        -0.3261822542307636,  0.0100948631354257, -0.0249548710136595,  0.000772316716123 ,
        -0.4126284632147356,  0.0127702467191452, -0.0249548710136594,  0.000772316716123 ,
        
         0.1583369220970716, -0.013763736916611 , -0.0205266499095721,  0.0017843179303531,
         0.0872305209919558, -0.0075826782921555, -0.020526649909572,  0.0017843179303531,
         0.0161241198868402, -0.0014016196677002, -0.020526649909572,  0.0017843179303531,
        -0.0549822812182755,  0.0047794389567552, -0.020526649909572,  0.0017843179303531,
        -0.1260886823233912,  0.0109604975812106, -0.020526649909572,  0.0017843179303531,
        -0.1971950834285068,  0.017141556205666 , -0.020526649909572,  0.0017843179303531,
        -0.2683014845336226,  0.0233226148301214, -0.020526649909572,  0.0017843179303531,
        -0.3394078856387381,  0.0295036734545768, -0.020526649909572,  0.0017843179303531,
        
         0.1048226005613422, -0.0171327710275426, -0.0135891035131684,  0.0022210763491261,
         0.0577485651331015, -0.0094387368591906, -0.0135891035131684,  0.0022210763491261,
         0.0106745297048609, -0.0017447026908387, -0.0135891035131683,  0.0022210763491261,
        -0.0363995057233797,  0.0059493314775132, -0.0135891035131683,  0.0022210763491261,
        -0.0834735411516203,  0.0136433656458652, -0.0135891035131683,  0.0022210763491261,
        -0.1305475765798609,  0.0213373998142171, -0.0135891035131683,  0.0022210763491261,
        -0.1776216120081016,  0.0290314339825691, -0.0135891035131684,  0.0022210763491261,
        -0.2246956474363422,  0.036725468150921 , -0.0135891035131683,  0.0022210763491261,
        
         0.0490310738538797, -0.0150784819355111, -0.0063563423764922,  0.0019547602401182,
         0.0270120579611045, -0.0083069938305107, -0.0063563423764922,  0.0019547602401182,
         0.0049930420683294, -0.0015355057255103, -0.0063563423764922,  0.0019547602401182,
        -0.0170259738244458,  0.0052359823794901, -0.0063563423764922,  0.0019547602401182,
        -0.0390449897172209,  0.0120074704844905, -0.0063563423764922,  0.0019547602401182,
        -0.0610640056099961,  0.0187789585894908, -0.0063563423764922,  0.0019547602401182,
        -0.0830830215027712,  0.0255504466944912, -0.0063563423764922,  0.0019547602401182,
        -0.1051020373955464,  0.0323219347994916, -0.0063563423764922,  0.0019547602401182,
        
         0.0114571885202825, -0.0066148108763414, -0.001485299161996,  0.0008575378710055,
         0.0063119612942494, -0.0036442125523494, -0.001485299161996,  0.0008575378710055,
         0.0011667340682163, -0.0006736142283574, -0.001485299161996,  0.0008575378710055,
        -0.0039784931578168,  0.0022969840956346, -0.001485299161996,  0.0008575378710055,
        -0.0091237203838499,  0.0052675824196266, -0.001485299161996,  0.0008575378710055,
        -0.014268947609883 ,  0.0082381807436186, -0.001485299161996,  0.0008575378710055,
        -0.0194141748359161,  0.0112087790676106, -0.001485299161996,  0.0008575378710055,
        -0.0245594020619492,  0.0141793773916026, -0.001485299161996,  0.0008575378710055,
  };
  return sol;
}

static const double *
get_sol_2x_p1_NxDx(void)
{
  static const double sol[256] = {
         0.2779969376054874, -0.0077857556587993, -0.0360392619648285,  0.00100933805314  ,
         0.153153272024753 , -0.0042893060786969, -0.0360392619648284,  0.00100933805314  ,
         0.0283096064440185, -0.0007928564985946, -0.0360392619648283,  0.00100933805314  ,
        -0.0965340591367159,  0.0027035930815078, -0.0360392619648281,  0.00100933805314  ,
        -0.2213777247174505,  0.0062000426616102, -0.0360392619648281,  0.00100933805314  ,
        -0.346221390298185 ,  0.0096964922417126, -0.0360392619648282,  0.00100933805314  ,
        -0.4710650558789194,  0.0131929418218149, -0.0360392619648277,  0.00100933805314  ,
        -0.5959087214596539,  0.0166893914019174, -0.0360392619648275,  0.00100933805314  ,
        
         0.226190519948558 , -0.0221246935211791, -0.0293231266236292,  0.0028682244940153,
         0.1246122296507558, -0.0121888724188829, -0.0293231266236291,  0.0028682244940153,
         0.0230339393529535, -0.0022530513165867, -0.029323126623629,  0.0028682244940153,
        -0.0785443509448487,  0.0076827697857095, -0.0293231266236288,  0.0028682244940153,
        -0.180122641242651 ,  0.0176185908880057, -0.0293231266236288,  0.0028682244940153,
        -0.2817009315404533,  0.027554411990302 , -0.0293231266236289,  0.0028682244940153,
        -0.3832792218382556,  0.0374902330925981, -0.0293231266236285,  0.0028682244940153,
        -0.4848575121360578,  0.0474260541948944, -0.0293231266236283,  0.0028682244940153,
        
         0.1311172040286984, -0.0327659110179022, -0.0169979112172524,  0.0042477419386764,
         0.0722347123270466, -0.0180513013074449, -0.0169979112172524,  0.0042477419386764,
         0.0133522206253946, -0.0033366915969877, -0.0169979112172523,  0.0042477419386764,
        -0.0455302710762573,  0.0113779181134696, -0.0169979112172521,  0.0042477419386764,
        -0.1044127627779092,  0.0260925278239268, -0.0169979112172521,  0.0042477419386764,
        -0.1632952544795612,  0.040807137534384 , -0.0169979112172522,  0.0042477419386764,
        -0.2221777461812131,  0.0555217472448413, -0.0169979112172518,  0.0042477419386764,
        -0.2810602378828649,  0.0702363569552985, -0.0169979112172516,  0.0042477419386764,
        
         0.0098560286339075, -0.0372442612385308, -0.001277726297742,  0.0048283110563137,
         0.0054298549021959, -0.0205185011099672, -0.001277726297742,  0.0048283110563137,
         0.0010036811704842, -0.0037927409814036, -0.0012777262977419,  0.0048283110563137,
        -0.0034224925612274,  0.0129330191471601, -0.0012777262977417,  0.0048283110563136,
        -0.0078486662929391,  0.0296587792757237, -0.0012777262977417,  0.0048283110563136,
        -0.0122748400246508,  0.0463845394042873, -0.0012777262977418,  0.0048283110563137,
        -0.0167010137563625,  0.0631102995328509, -0.0012777262977415,  0.0048283110563136,
        -0.0211271874880741,  0.0798360596614145, -0.0012777262977413,  0.0048283110563136,
        
        -0.1119744480538166, -0.0330945972726272,  0.0145162623068366,  0.0042903525161177,
        -0.0616886403509404, -0.0182323801920337,  0.0145162623068366,  0.0042903525161177,
        -0.0114028326480642, -0.0033701631114402,  0.0145162623068367,  0.0042903525161177,
         0.038882975054812 ,  0.0114920539691533,  0.0145162623068368,  0.0042903525161176,
         0.0891687827576881,  0.0263542710497468,  0.0145162623068368,  0.0042903525161176,
         0.1394545904605643,  0.0412164881303402,  0.0145162623068367,  0.0042903525161176,
         0.1897403981634405,  0.0560787052109337,  0.014516262306837,  0.0042903525161176,
         0.2400262058663167,  0.0709409222915272,  0.0145162623068371,  0.0042903525161176,
        
        -0.2002161484584767, -0.0178517722097536,  0.0259558334923958,  0.0023142869872789,
        -0.1103025037352214, -0.0098348469192286,  0.0259558334923958,  0.0023142869872789,
        -0.0203888590119662, -0.0018179216287036,  0.0259558334923959,  0.0023142869872789,
         0.0695247857112891,  0.0061990036618213,  0.025955833492396,  0.0023142869872788,
         0.1594384304345443,  0.0142159289523463,  0.025955833492396,  0.0023142869872788,
         0.2493520751577995,  0.0222328542428713,  0.0259558334923959,  0.0023142869872789,
         0.3392657198810548,  0.0302497795333963,  0.0259558334923961,  0.0023142869872788,
         0.42917936460431  ,  0.0382667048239213,  0.0259558334923963,  0.0023142869872788,
        
        -0.2121714756100758,  0.0109493608605278,  0.0275057108788264, -0.0014194648610122,
        -0.1168888981292207,  0.006032190342864 ,  0.0275057108788265, -0.0014194648610122,
        -0.0216063206483658,  0.0011150198252003,  0.0275057108788265, -0.0014194648610122,
         0.0736762568324892, -0.0038021506924635,  0.0275057108788265, -0.0014194648610122,
         0.1689588343133442, -0.0087193212101273,  0.0275057108788266, -0.0014194648610122,
         0.2642414117941991, -0.013636491727791 ,  0.0275057108788265, -0.0014194648610122,
         0.3595239892750541, -0.0185536622454548,  0.0275057108788266, -0.0014194648610123,
         0.4548065667559091, -0.0234708327631186,  0.0275057108788267, -0.0014194648610123,
        
        -0.0966033131446177,  0.0557739488486548,  0.0125235628099973, -0.0072304823595651,
        -0.053220418987227 ,  0.0307268232286602,  0.0125235628099973, -0.0072304823595651,
        -0.0098375248298362,  0.0056796976086655,  0.0125235628099973, -0.0072304823595651,
         0.0335453693275546, -0.0193674280113291,  0.0125235628099973, -0.0072304823595651,
         0.0769282634849454, -0.0444145536313237,  0.0125235628099973, -0.0072304823595651,
         0.1203111576423362, -0.0694616792513183,  0.0125235628099973, -0.0072304823595651,
         0.163694051799727 , -0.0945088048713129,  0.0125235628099974, -0.0072304823595652,
         0.2070769459571177, -0.1195559304913075,  0.0125235628099974, -0.0072304823595652,
  };
  return sol;
}

static const double *
get_sol_2x_p1_DxNx(void)
{
  static const double sol[256] = {
        -0.0966033131446173, -0.0557739488486548,  0.012523562809997 , 0.0072304823595651,
        -0.053220418987227 , -0.0307268232286602,  0.0125235628099973,  0.0072304823595651,
        -0.0098375248298362, -0.0056796976086655,  0.0125235628099973,  0.0072304823595651,
         0.0335453693275545,  0.0193674280113291,  0.0125235628099973,  0.0072304823595651,
         0.0769282634849452,  0.0444145536313237,  0.0125235628099971,  0.0072304823595651,
         0.1203111576423362,  0.0694616792513183,  0.0125235628099974,  0.0072304823595651,
         0.1636940517997271,  0.0945088048713129,  0.0125235628099971,  0.0072304823595652,
         0.2070769459571174,  0.1195559304913075,  0.0125235628099976,  0.0072304823595651,
        
        -0.2121714756100753, -0.0109493608605278,  0.0275057108788262,  0.0014194648610122,
        -0.1168888981292208, -0.0060321903428641,  0.0275057108788265,  0.0014194648610122,
        -0.0216063206483658, -0.0011150198252003,  0.0275057108788265,  0.0014194648610122,
         0.0736762568324891,  0.0038021506924635,  0.0275057108788265,  0.0014194648610122,
         0.168958834313344 ,  0.0087193212101273,  0.0275057108788263,  0.0014194648610122,
         0.2642414117941991,  0.013636491727791 ,  0.0275057108788266,  0.0014194648610122,
         0.3595239892750542,  0.0185536622454548,  0.0275057108788264,  0.0014194648610122,
         0.4548065667559088,  0.0234708327631186,  0.0275057108788267,  0.0014194648610122,
        
        -0.2002161484584764,  0.0178517722097536,  0.0259558334923957, -0.0023142869872788,
        -0.1103025037352215,  0.0098348469192286,  0.0259558334923959, -0.0023142869872789,
        -0.0203888590119662,  0.0018179216287036,  0.0259558334923959, -0.0023142869872789,
         0.069524785711289 , -0.0061990036618213,  0.0259558334923959, -0.0023142869872789,
         0.1594384304345441, -0.0142159289523463,  0.0259558334923958, -0.0023142869872788,
         0.2493520751577995, -0.0222328542428713,  0.0259558334923959, -0.0023142869872789,
         0.3392657198810548, -0.0302497795333963,  0.0259558334923958, -0.0023142869872789,
         0.4291793646043098, -0.0382667048239212,  0.0259558334923961, -0.0023142869872789,
        
        -0.1119744480538163,  0.0330945972726272,  0.0145162623068365, -0.0042903525161176,
        -0.0616886403509405,  0.0182323801920337,  0.0145162623068366, -0.0042903525161177,
        -0.0114028326480643,  0.0033701631114402,  0.0145162623068367, -0.0042903525161177,
         0.0388829750548119, -0.0114920539691533,  0.0145162623068367, -0.0042903525161177,
         0.0891687827576879, -0.0263542710497468,  0.0145162623068366, -0.0042903525161176,
         0.1394545904605643, -0.0412164881303402,  0.0145162623068367, -0.0042903525161177,
         0.1897403981634405, -0.0560787052109337,  0.0145162623068366, -0.0042903525161176,
         0.2400262058663166, -0.0709409222915272,  0.0145162623068367, -0.0042903525161177,
        
         0.0098560286339078,  0.0372442612385308, -0.0012777262977421, -0.0048283110563137,
         0.0054298549021959,  0.0205185011099672, -0.0012777262977419, -0.0048283110563137,
         0.0010036811704842,  0.0037927409814036, -0.0012777262977419, -0.0048283110563137,
        -0.0034224925612275, -0.0129330191471601, -0.0012777262977419, -0.0048283110563137,
        -0.0078486662929393, -0.0296587792757237, -0.001277726297742 , -0.0048283110563136,
        -0.0122748400246508, -0.0463845394042873, -0.0012777262977419, -0.0048283110563137,
        -0.0167010137563626, -0.063110299532851 , -0.0012777262977419, -0.0048283110563136,
        -0.0211271874880741, -0.0798360596614145, -0.0012777262977419, -0.0048283110563137,
        
         0.1311172040286987,  0.0327659110179021, -0.0169979112172525, -0.0042477419386764,
         0.0722347123270465,  0.0180513013074449, -0.0169979112172524, -0.0042477419386764,
         0.0133522206253946,  0.0033366915969877, -0.0169979112172523, -0.0042477419386764,
        -0.0455302710762574, -0.0113779181134696, -0.0169979112172523, -0.0042477419386764,
        -0.1044127627779095, -0.0260925278239268, -0.0169979112172524, -0.0042477419386764,
        -0.1632952544795612, -0.040807137534384 , -0.0169979112172523, -0.0042477419386764,
        -0.2221777461812133, -0.0555217472448413, -0.0169979112172523, -0.0042477419386764,
        -0.2810602378828649, -0.0702363569552985, -0.0169979112172524, -0.0042477419386764,
        
         0.2261905199485582,  0.0221246935211791, -0.0293231266236291, -0.0028682244940153,
         0.1246122296507557,  0.0121888724188829, -0.0293231266236291, -0.0028682244940153,
         0.0230339393529534,  0.0022530513165867, -0.029323126623629 , -0.0028682244940153,
        -0.0785443509448489, -0.0076827697857095, -0.029323126623629 , -0.0028682244940153,
        -0.1801226412426513, -0.0176185908880057, -0.029323126623629 , -0.0028682244940153,
        -0.2817009315404534, -0.0275544119903019, -0.0293231266236291, -0.0028682244940153,
        -0.3832792218382559, -0.0374902330925982, -0.029323126623629 , -0.0028682244940153,
        -0.4848575121360577, -0.0474260541948943, -0.0293231266236291, -0.0028682244940153,
        
         0.2779969376054876,  0.0077857556587993, -0.0360392619648284, -0.00100933805314  ,
         0.153153272024753 ,  0.0042893060786969, -0.0360392619648284, -0.00100933805314  ,
         0.0283096064440184,  0.0007928564985946, -0.0360392619648283, -0.00100933805314  ,
        -0.0965340591367161, -0.0027035930815078, -0.0360392619648283, -0.00100933805314  ,
        -0.2213777247174507, -0.0062000426616102, -0.0360392619648283, -0.00100933805314  ,
        -0.346221390298185 , -0.0096964922417126, -0.0360392619648284, -0.00100933805314  ,
        -0.4710650558789197, -0.013192941821815 , -0.0360392619648282, -0.00100933805314  ,
        -0.5959087214596538, -0.0166893914019173, -0.0360392619648285, -0.00100933805314  ,
  };
  return sol;
}

static const double *
get_sol_2x_bias_p1_DxDx(void)
{
  static const double sol[256] = {
         0.0523071668741079,  0.0301995568753129,  0.0043898894860765,  0.0025345038764989,
         0.0326575388702468,  0.0188548388578077, -0.0076848087077184, -0.004436826376072 ,
         0.0060365806138147,  0.0034852214423708, -0.0076848087077184, -0.004436826376072 ,
        -0.0205843776426173, -0.011884395973066 , -0.0076848087077184, -0.004436826376072 ,
        -0.0472053358990496, -0.0272540133885028, -0.0076848087077183, -0.004436826376072 ,
        -0.0738262941554814, -0.0426236308039397, -0.0076848087077184, -0.004436826376072 ,
        -0.1004472524119136, -0.0579932482193765, -0.0076848087077186, -0.004436826376072 ,
        -0.1270682106683453, -0.0733628656348134, -0.0076848087077185, -0.004436826376072 ,
                                                                       
         0.1305913158241592,  0.0149978175942799,  0.0165830888926301,  0.0045052430831576,
         0.0834668670525276,  0.0104799404455762, -0.0196410056887545, -0.0024660871694134,
         0.0154284275231832,  0.0019371638989408, -0.0196410056887545, -0.0024660871694134,
        -0.0526100120061611, -0.0066056126476946, -0.0196410056887545, -0.0024660871694134,
        -0.1206484515355057, -0.0151483891943299, -0.0196410056887544, -0.0024660871694134,
        -0.1886868910648499, -0.0236911657409653, -0.0196410056887545, -0.0024660871694134,
        -0.2567253305941943, -0.0322339422876006, -0.0196410056887547, -0.0024660871694134,
        -0.3247637701235385, -0.040776718834236 , -0.0196410056887545, -0.0024660871694134,
                                                                       
         0.1587769509960808,  0.0012751664595096,  0.0352710125851786,  0.0062842346913305,
         0.1066760663456682,  0.0029198970139942, -0.0251024783837957, -0.0006870955612404,
         0.0197185304324009,  0.0005397281705472, -0.0251024783837957, -0.0006870955612404,
        -0.0672390054808664, -0.0018404406728998, -0.0251024783837957, -0.0006870955612404,
        -0.1541965413941337, -0.0042206095163468, -0.0251024783837957, -0.0006870955612405,
        -0.2411540773074009, -0.0066007783597938, -0.0251024783837957, -0.0006870955612404,
        -0.3281116132206683, -0.0089809472032408, -0.0251024783837958, -0.0006870955612404,
        -0.4150691491339354, -0.0113611160466878, -0.0251024783837957, -0.0006870955612404,
                                                                       
         0.1436956879050723, -0.0099823377648229,  0.0595680163429045,  0.0077436469686939,
         0.1060487906890969, -0.0032820547831718, -0.0249548710136595,  0.000772316716123 ,
         0.0196025817051248, -0.0006066711994523, -0.0249548710136594,  0.000772316716123 ,
        -0.0668436272788473,  0.0020687123842672, -0.0249548710136594,  0.000772316716123 ,
        -0.1532898362628194,  0.0047440959679867, -0.0249548710136594,  0.000772316716123 ,
        -0.2397360452467914,  0.0074194795517062, -0.0249548710136595,  0.000772316716123 ,
        -0.3261822542307636,  0.0100948631354257, -0.0249548710136595,  0.000772316716123 ,
        -0.4126284632147356,  0.0127702467191452, -0.0249548710136594,  0.000772316716123 ,
                                                                       
         0.1095376103290749, -0.0097388375186793,  0.0639962374469919, -0.0051870123222179,
         0.0872305209919558, -0.0075826782921555, -0.020526649909572 ,  0.0017843179303531,
         0.0161241198868402, -0.0014016196677002, -0.020526649909572 ,  0.0017843179303531,
        -0.0549822812182755,  0.0047794389567552, -0.020526649909572 ,  0.0017843179303531,
        -0.1260886823233912,  0.0109604975812106, -0.020526649909572 ,  0.0017843179303531,
        -0.1971950834285068,  0.017141556205666 , -0.020526649909572 ,  0.0017843179303531,
        -0.2683014845336226,  0.0233226148301214, -0.020526649909572 ,  0.0017843179303531,
        -0.3394078856387381,  0.0295036734545768, -0.020526649909572 ,  0.0017843179303531,
                                                                       
         0.0699659492984875, -0.013107871629611 ,  0.0467843874558059, -0.0047502539034449,
         0.0577485651331015, -0.0094387368591906, -0.0135891035131684,  0.0022210763491261,
         0.0106745297048609, -0.0017447026908387, -0.0135891035131683,  0.0022210763491261,
        -0.0363995057233797,  0.0059493314775132, -0.0135891035131683,  0.0022210763491261,
        -0.0834735411516203,  0.0136433656458652, -0.0135891035131683,  0.0022210763491261,
        -0.1305475765798609,  0.0213373998142171, -0.0135891035131683,  0.0022210763491261,
        -0.1776216120081016,  0.0290314339825691, -0.0135891035131684,  0.0022210763491261,
        -0.2246956474363422,  0.036725468150921 , -0.0135891035131683,  0.0022210763491261,
                                                                       
         0.0281170830961668, -0.0110535825375795,  0.0298677522048923, -0.0050165700124528,
         0.0270120579611045, -0.0083069938305107, -0.0063563423764922,  0.0019547602401182,
         0.0049930420683294, -0.0015355057255103, -0.0063563423764922,  0.0019547602401182,
        -0.0170259738244458,  0.0052359823794901, -0.0063563423764922,  0.0019547602401182,
        -0.0390449897172209,  0.0120074704844905, -0.0063563423764922,  0.0019547602401182,
        -0.0610640056099961,  0.0187789585894908, -0.0063563423764922,  0.0019547602401182,
        -0.0830830215027712,  0.0255504466944912, -0.0063563423764922,  0.0019547602401182,
        -0.1051020373955464,  0.0323219347994916, -0.0063563423764922,  0.0019547602401182,
                                                                       
         0.0044858582677115, -0.0025899114784098,  0.0105893990317988, -0.0061137923815654,
         0.0063119612942494, -0.0036442125523494, -0.001485299161996 ,  0.0008575378710055,
         0.0011667340682163, -0.0006736142283574, -0.001485299161996 ,  0.0008575378710055,
        -0.0039784931578168,  0.0022969840956346, -0.001485299161996 ,  0.0008575378710055,
        -0.0091237203838499,  0.0052675824196266, -0.001485299161996 ,  0.0008575378710055,
        -0.014268947609883 ,  0.0082381807436186, -0.001485299161996 ,  0.0008575378710055,
        -0.0194141748359161,  0.0112087790676106, -0.001485299161996 ,  0.0008575378710055,
        -0.0245594020619492,  0.0141793773916026, -0.001485299161996 ,  0.0008575378710055,
  };
  return sol;
}

static const double *
get_sol_2x_bias_p1_NxDx(void)
{
  static const double sol[256] = {
         0.2947281302116576, -0.0077857556587993, -0.0650185376299358,  0.00100933805314  ,
         0.153153272024753 , -0.0042893060786969, -0.0360392619648284,  0.00100933805314  ,
         0.0283096064440185, -0.0007928564985946, -0.0360392619648283,  0.00100933805314  ,
        -0.0965340591367159,  0.0027035930815078, -0.0360392619648281,  0.00100933805314  ,
        -0.2213777247174505,  0.0062000426616102, -0.0360392619648281,  0.00100933805314  ,
        -0.346221390298185 ,  0.0096964922417126, -0.0360392619648282,  0.00100933805314  ,
        -0.4710650558789194,  0.0131929418218149, -0.0360392619648277,  0.00100933805314  ,
        -0.5959087214596539,  0.0166893914019174, -0.0360392619648275,  0.00100933805314  ,
                                                                                           
         0.2429217125547281, -0.0221246935211791, -0.0583024022887365,  0.0028682244940153,
         0.1246122296507558, -0.0121888724188829, -0.0293231266236291,  0.0028682244940153,
         0.0230339393529535, -0.0022530513165867, -0.029323126623629 ,  0.0028682244940153,
        -0.0785443509448487,  0.0076827697857095, -0.0293231266236288,  0.0028682244940153,
        -0.180122641242651 ,  0.0176185908880057, -0.0293231266236288,  0.0028682244940153,
        -0.2817009315404533,  0.027554411990302 , -0.0293231266236289,  0.0028682244940153,
        -0.3832792218382556,  0.0374902330925981, -0.0293231266236285,  0.0028682244940153,
        -0.4848575121360578,  0.0474260541948944, -0.0293231266236283,  0.0028682244940153,
                                                                                           
         0.1478483966348686, -0.0327659110179021, -0.0459771868823598,  0.0042477419386764,
         0.0722347123270466, -0.0180513013074449, -0.0169979112172524,  0.0042477419386764,
         0.0133522206253946, -0.0033366915969877, -0.0169979112172523,  0.0042477419386764,
        -0.0455302710762573,  0.0113779181134696, -0.0169979112172521,  0.0042477419386764,
        -0.1044127627779092,  0.0260925278239268, -0.0169979112172521,  0.0042477419386764,
        -0.1632952544795612,  0.040807137534384 , -0.0169979112172522,  0.0042477419386764,
        -0.2221777461812131,  0.0555217472448413, -0.0169979112172518,  0.0042477419386764,
        -0.2810602378828649,  0.0702363569552985, -0.0169979112172516,  0.0042477419386764,
                                                                                           
         0.0265872212400777, -0.0372442612385308, -0.0302570019628494,  0.0048283110563137,
         0.0054298549021959, -0.0205185011099672, -0.001277726297742 ,  0.0048283110563137,
         0.0010036811704842, -0.0037927409814036, -0.0012777262977419,  0.0048283110563137,
        -0.0034224925612274,  0.0129330191471601, -0.0012777262977417,  0.0048283110563136,
        -0.0078486662929391,  0.0296587792757237, -0.0012777262977417,  0.0048283110563136,
        -0.0122748400246508,  0.0463845394042873, -0.0012777262977418,  0.0048283110563137,
        -0.0167010137563625,  0.0631102995328509, -0.0012777262977415,  0.0048283110563136,
        -0.0211271874880741,  0.0798360596614145, -0.0012777262977413,  0.0048283110563136,
                                                                                           
        -0.0973346545234177, -0.0343020670920067, -0.0108406039001324,  0.0063817515918889,
        -0.0616886403509404, -0.0182323801920337,  0.0145162623068366,  0.0042903525161177,
        -0.0114028326480642, -0.0033701631114402,  0.0145162623068367,  0.0042903525161177,
         0.038882975054812 ,  0.0114920539691533,  0.0145162623068368,  0.0042903525161176,
         0.0891687827576881,  0.0263542710497468,  0.0145162623068368,  0.0042903525161176,
         0.1394545904605643,  0.0412164881303402,  0.0145162623068367,  0.0042903525161176,
         0.1897403981634405,  0.0560787052109337,  0.014516262306837 ,  0.0042903525161176,
         0.2400262058663167,  0.0709409222915272,  0.0145162623068371,  0.0042903525161176,
                                                                                           
        -0.1897591530796203, -0.0190592420291331,  0.0078437862017036,  0.0044056860630501,
        -0.1103025037352214, -0.0098348469192286,  0.0259558334923958,  0.0023142869872789,
        -0.0203888590119662, -0.0018179216287036,  0.0259558334923959,  0.0023142869872789,
         0.0695247857112891,  0.0061990036618213,  0.025955833492396 ,  0.0023142869872788,
         0.1594384304345443,  0.0142159289523463,  0.025955833492396 ,  0.0023142869872788,
         0.2493520751577995,  0.0222328542428713,  0.0259558334923959,  0.0023142869872789,
         0.3392657198810548,  0.0302497795333963,  0.0259558334923961,  0.0023142869872788,
         0.42917936460431  ,  0.0382667048239213,  0.0259558334923963,  0.0023142869872788,
                                                                                           
        -0.2058972783827619,  0.0097418910411483,  0.0166384825044111,  0.000671934214759 ,
        -0.1168888981292207,  0.006032190342864 ,  0.0275057108788265, -0.0014194648610122,
        -0.0216063206483658,  0.0011150198252003,  0.0275057108788265, -0.0014194648610122,
         0.0736762568324892, -0.0038021506924635,  0.0275057108788265, -0.0014194648610122,
         0.1689588343133442, -0.0087193212101273,  0.0275057108788266, -0.0014194648610122,
         0.2642414117941991, -0.013636491727791 ,  0.0275057108788265, -0.0014194648610122,
         0.3595239892750541, -0.0185536622454548,  0.0275057108788266, -0.0014194648610123,
         0.4548065667559091, -0.0234708327631186,  0.0275057108788267, -0.0014194648610123,
                                                                                           
        -0.0945119140688464,  0.0545664790292753,  0.0089011533518588, -0.0051390832837939,
        -0.053220418987227 ,  0.0307268232286602,  0.0125235628099973, -0.0072304823595651,
        -0.0098375248298362,  0.0056796976086655,  0.0125235628099973, -0.0072304823595651,
         0.0335453693275546, -0.0193674280113291,  0.0125235628099973, -0.0072304823595651,
         0.0769282634849454, -0.0444145536313237,  0.0125235628099973, -0.0072304823595651,
         0.1203111576423362, -0.0694616792513183,  0.0125235628099973, -0.0072304823595651,
         0.163694051799727 , -0.0945088048713129,  0.0125235628099974, -0.0072304823595652,
         0.2070769459571177, -0.1195559304913075,  0.0125235628099974, -0.0072304823595652,
  };
  return sol;
}

static const double *
get_sol_2x_bias_p1_DxNx(void)
{
  static const double sol[256] = {
        -0.0945119140688464, -0.0545664790292753,  0.0089011533518588,  0.0051390832837939,
        -0.053220418987227 , -0.0307268232286602,  0.0125235628099973,  0.0072304823595651,
        -0.0098375248298362, -0.0056796976086655,  0.0125235628099973,  0.0072304823595651,
         0.0335453693275545,  0.0193674280113291,  0.0125235628099973,  0.0072304823595651,
         0.0769282634849452,  0.0444145536313237,  0.0125235628099971,  0.0072304823595651,
         0.1203111576423362,  0.0694616792513183,  0.0125235628099974,  0.0072304823595651,
         0.1636940517997271,  0.0945088048713129,  0.0125235628099971,  0.0072304823595652,
         0.2070769459571174,  0.1195559304913075,  0.0125235628099976,  0.0072304823595651,
                                                                       
        -0.2058972783827619, -0.0097418910411483,  0.0166384825044111, -0.000671934214759 ,
        -0.1168888981292208, -0.0060321903428641,  0.0275057108788265,  0.0014194648610122,
        -0.0216063206483658, -0.0011150198252003,  0.0275057108788265,  0.0014194648610122,
         0.0736762568324891,  0.0038021506924635,  0.0275057108788265,  0.0014194648610122,
         0.168958834313344 ,  0.0087193212101273,  0.0275057108788263,  0.0014194648610122,
         0.2642414117941991,  0.013636491727791 ,  0.0275057108788266,  0.0014194648610122,
         0.3595239892750542,  0.0185536622454548,  0.0275057108788264,  0.0014194648610122,
         0.4548065667559088,  0.0234708327631186,  0.0275057108788267,  0.0014194648610122,
                                                                       
        -0.1897591530796203,  0.0190592420291331,  0.0078437862017036, -0.0044056860630501,
        -0.1103025037352215,  0.0098348469192286,  0.0259558334923959, -0.0023142869872789,
        -0.0203888590119662,  0.0018179216287036,  0.0259558334923959, -0.0023142869872789,
         0.069524785711289 , -0.0061990036618213,  0.0259558334923959, -0.0023142869872789,
         0.1594384304345441, -0.0142159289523463,  0.0259558334923958, -0.0023142869872788,
         0.2493520751577995, -0.0222328542428713,  0.0259558334923959, -0.0023142869872789,
         0.3392657198810548, -0.0302497795333963,  0.0259558334923958, -0.0023142869872789,
         0.4291793646043098, -0.0382667048239212,  0.0259558334923961, -0.0023142869872789,
                                                                       
        -0.0973346545234177,  0.0343020670920067, -0.0108406039001325, -0.0063817515918889,
        -0.0616886403509405,  0.0182323801920337,  0.0145162623068366, -0.0042903525161177,
        -0.0114028326480643,  0.0033701631114402,  0.0145162623068367, -0.0042903525161177,
         0.0388829750548119, -0.0114920539691533,  0.0145162623068367, -0.0042903525161177,
         0.0891687827576879, -0.0263542710497468,  0.0145162623068366, -0.0042903525161176,
         0.1394545904605643, -0.0412164881303402,  0.0145162623068367, -0.0042903525161177,
         0.1897403981634405, -0.0560787052109337,  0.0145162623068366, -0.0042903525161176,
         0.2400262058663166, -0.0709409222915272,  0.0145162623068367, -0.0042903525161177,
                                                                       
         0.0265872212400778,  0.0372442612385308, -0.0302570019628495, -0.0048283110563137,
         0.0054298549021959,  0.0205185011099672, -0.0012777262977419, -0.0048283110563137,
         0.0010036811704842,  0.0037927409814036, -0.0012777262977419, -0.0048283110563137,
        -0.0034224925612275, -0.0129330191471601, -0.0012777262977419, -0.0048283110563137,
        -0.0078486662929393, -0.0296587792757237, -0.001277726297742 , -0.0048283110563136,
        -0.0122748400246508, -0.0463845394042873, -0.0012777262977419, -0.0048283110563137,
        -0.0167010137563626, -0.063110299532851 , -0.0012777262977419, -0.0048283110563136,
        -0.0211271874880741, -0.0798360596614145, -0.0012777262977419, -0.0048283110563137,
                                                                       
         0.1478483966348688,  0.0327659110179022, -0.04597718688236  , -0.0042477419386764,
         0.0722347123270465,  0.0180513013074449, -0.0169979112172524, -0.0042477419386764,
         0.0133522206253946,  0.0033366915969877, -0.0169979112172523, -0.0042477419386764,
        -0.0455302710762574, -0.0113779181134696, -0.0169979112172523, -0.0042477419386764,
        -0.1044127627779095, -0.0260925278239268, -0.0169979112172524, -0.0042477419386764,
        -0.1632952544795612, -0.040807137534384 , -0.0169979112172523, -0.0042477419386764,
        -0.2221777461812133, -0.0555217472448413, -0.0169979112172523, -0.0042477419386764,
        -0.2810602378828649, -0.0702363569552985, -0.0169979112172524, -0.0042477419386764,
                                                                       
         0.2429217125547283,  0.0221246935211791, -0.0583024022887367, -0.0028682244940153,
         0.1246122296507557,  0.0121888724188829, -0.0293231266236291, -0.0028682244940153,
         0.0230339393529534,  0.0022530513165867, -0.029323126623629 , -0.0028682244940153,
        -0.0785443509448489, -0.0076827697857095, -0.029323126623629 , -0.0028682244940153,
        -0.1801226412426513, -0.0176185908880057, -0.029323126623629 , -0.0028682244940153,
        -0.2817009315404534, -0.0275544119903019, -0.0293231266236291, -0.0028682244940153,
        -0.3832792218382559, -0.0374902330925982, -0.029323126623629 , -0.0028682244940153,
        -0.4848575121360577, -0.0474260541948943, -0.0293231266236291, -0.0028682244940153,
                                                                       
         0.2947281302116578,  0.0077857556587993, -0.065018537629936 , -0.00100933805314  ,
         0.153153272024753 ,  0.0042893060786969, -0.0360392619648284, -0.00100933805314  ,
         0.0283096064440184,  0.0007928564985946, -0.0360392619648283, -0.00100933805314  ,
        -0.0965340591367161, -0.0027035930815078, -0.0360392619648283, -0.00100933805314  ,
        -0.2213777247174507, -0.0062000426616102, -0.0360392619648283, -0.00100933805314  ,
        -0.346221390298185 , -0.0096964922417126, -0.0360392619648284, -0.00100933805314  ,
        -0.4710650558789197, -0.013192941821815 , -0.0360392619648282, -0.00100933805314  ,
        -0.5959087214596538, -0.0166893914019173, -0.0360392619648285, -0.00100933805314  ,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_DxDx_DyDy(void)
{
  static const double sol[512] = {
        -1.4857902425065184e-04, -8.5782139646996372e-05, -8.5782139646865672e-05, 1.8914704406197892e-04, -4.9526341416772864e-05, 1.0920409680558287e-04,
        1.0920409680550554e-04, 6.3049014687256998e-05, -3.7869881831633913e-04, -2.1864186469671895e-04, -4.7077585402820190e-05, 4.8209875139245301e-04,
        -2.7180256605144123e-05, 2.7833984389239624e-04, 5.9931650281304105e-05, 3.4601554422903200e-05, -4.8238508035628331e-04, -2.7850515599677118e-04,
        -1.2785705897207671e-05, 6.1409551253425998e-04, -7.3818307415179943e-06, 3.5454820946987623e-04, 1.6276715296076932e-05, 9.3973659577562421e-06,
        -4.7821812488875947e-04, -2.7609936313588445e-04, 1.5191498758083683e-05, 6.0879081146090802e-04, 8.7708158973654205e-06, 3.5148553887719326e-04,
        -1.9339385888747259e-05, -1.1165599648880918e-05, -3.9208040614287438e-04, -2.2636772803058295e-04, 3.4540136347233381e-05, 4.9913404823207755e-04,
        1.9941757017920681e-05, 2.8817517710849280e-04, -4.3970975879881029e-05, -2.5386654760770837e-05, -2.5807022372215224e-04, -1.4899691313582957e-04,
        4.2830678547521361e-05, 3.2853372287024223e-04, 2.4728303788985857e-05, 1.8967903333702248e-04, -5.4525167891612326e-05, -3.1480120359827638e-05,
        -1.1891747993225567e-04, -6.8657039050248661e-05, 3.7509195538050490e-05, 1.5138671107812069e-04, 2.1655944140984915e-05, 8.7403158392697900e-05,
        -4.7750707052706831e-05, -2.7568883570882659e-05, -2.6974823754671981e-05, -1.5573921756102582e-05, 1.5573921756086665e-05, 3.4340030182754069e-05,
        8.9916079182148039e-06, 1.9826225669992933e-05, -1.9826225669985302e-05, -1.1446676727580283e-05, -3.7869881831637274e-04, -4.7077585402708903e-05,
        -2.1864186469670288e-04, 4.8209875139248608e-04, -2.7180256605210324e-05, 5.9931650281245863e-05, 2.7833984389239038e-04, 3.4601554422944427e-05,
        -9.6593915317034899e-04, -1.2040150071024966e-04, -1.2040150071024630e-04, 1.2296791992509439e-03, -1.5153325635635661e-05, 1.5327592891960440e-04,
        1.5327592891959123e-04, 1.9290789976223312e-05, -1.2313170233906195e-03, -1.5389090291881229e-04, -3.2814484098396489e-05, 1.5675158485681691e-03,
        -4.1817897444776702e-06, 1.9590927819016461e-04, 4.1774151505829872e-05, 5.3235857015899117e-06, -1.2211748055640191e-03, -1.5284687644820770e-04,
        3.8670096290763775e-05, 1.5546043993793627e-03, 4.7845587083223804e-06, 1.9458019071073231e-04, -4.9228580170615438e-05, -6.0909347156166817e-06,
        -1.0016098602570024e-03, -1.2554426638121464e-04, 8.8095783986847693e-05, 1.2750894369280114e-03, 1.0978610563435414e-05, 1.5982287543424525e-04,
        -1.1214945864320036e-04, -1.3976210615571825e-05, -6.5979684339352812e-04, -8.2940058872077652e-05, 1.0925005331179741e-04, 8.3994778696934870e-04,
        1.3618940110575512e-05, 1.0558601423795712e-04, -1.3907969009605144e-04, -1.7337455796105482e-05, -3.0491150826769746e-04, -3.8726663276881513e-05,
        9.5643090454551556e-05, 3.8816455270337567e-04, 1.1907675738098372e-05, 4.9300592207781632e-05, -1.2175748182276777e-04, -1.5158947764468026e-05,
        -6.9626408103776196e-05, -9.0509819492167675e-06, 4.0198825461396472e-05, 8.8637203992388735e-05, 5.2255868648144249e-06, 1.1522262245235528e-05,
        -5.1174713585175428e-05, -6.6523812089424644e-06, -4.8238508035620151e-04, -1.2785705897252288e-05, -2.7850515599680664e-04, 6.1409551253430368e-04,
        -7.3818307414815423e-06, 1.6276715296141225e-05, 3.5454820946987097e-04, 9.3973659577153728e-06, -1.2313170233906072e-03, -3.2814484098386040e-05,
        -1.5389090291881698e-04, 1.5675158485681897e-03, -4.1817897444823730e-06, 4.1774151505828591e-05, 1.9590927819015664e-04, 5.3235857015929212e-06,
        -1.5708705777981958e-03, -4.2150433122696279e-05, -4.2150433122698495e-05, 1.9997811123959371e-03, -1.2083229378437410e-06, 5.3659188242121851e-05,
        5.3659188242117866e-05, 1.5382434574345833e-06, -1.5587645637327890e-03, -4.2060661306192221e-05, 4.9139843602174565e-05, 1.9843696720032254e-03,
        1.2601527202679231e-06, 5.3544905127007844e-05, -6.2556987501449741e-05, -1.6042248447095148e-06, -1.2792161864866266e-03, -3.4731820796201924e-05,
        1.1225748725242052e-04, 1.6284934001327763e-03, 2.9711553210228856e-06, 4.4214997854706335e-05, -1.4290827386115453e-04, -3.7823996304679761e-06,
        -8.4357750896041524e-04, -2.3165757864813304e-05, 1.3925862048675021e-04, 1.0739079292105724e-03, 3.7065142258784857e-06, 2.9490936864656485e-05,
        -1.7728179706443526e-04, -4.7185409457017731e-06, -3.9114247992314450e-04, -1.1058811420873844e-05, 1.2195486531873429e-04, 4.9794003061819580e-04,
        3.2834345625943677e-06, 1.4078309516759583e-05, -1.5525342423245637e-04, -4.1799436025258761e-06, -8.9955228480428413e-05, -2.6858679675509360e-06,
        5.1935675378201352e-05, 1.1451660589930813e-04, 1.5506865940911167e-06, 3.4192174121908953e-06, -6.6116193242601623e-05, -1.9740860933638614e-06,
        -4.7821812488872613e-04, 1.5191498758100369e-05, -2.7609936313588413e-04, 6.0879081146098728e-04, 8.7708158973496454e-06, -1.9339385888790925e-05,
        3.5148553887716583e-04, -1.1165599648852878e-05, -1.2211748055639972e-03, 3.8670096290759126e-05, -1.5284687644821456e-04, 1.5546043993793848e-03,
        4.7845587083257541e-06, -4.9228580170613290e-05, 1.9458019071072680e-04, -6.0909347156182750e-06, -1.5587645637327831e-03, 4.9139843602175452e-05,
        -4.2060661306194789e-05, 1.9843696720032327e-03, 1.2601527202678208e-06, -6.2556987501449768e-05, 5.3544905127004523e-05, -1.6042248447091459e-06,
        -1.5474982642617466e-03, 4.8565262338569481e-05, 4.8565262338568776e-05, 1.9700272218948593e-03, -1.5918873674820159e-06, -6.1825522557911871e-05,
        -6.1825522557912508e-05, 2.0265363267635107e-06, -1.2706869384245629e-03, 3.9656184460819470e-05, 1.1125183114826664e-04, 1.6176353260059105e-03,
        -3.5517711434683191e-06, -5.0483909874733382e-05, -1.4162803339381479e-04, 4.5215468089140504e-06, -8.3879554668975795e-04, 2.5926625069031862e-05,
        1.3810078013050519e-04, 1.0678202998640151e-03, -4.3749936672351580e-06, -3.3005631306618052e-05, -1.7580781995370960e-04, 5.5695420273585997e-06,
        -3.8988499721200347e-04, 1.1784819402652368e-05, 1.2107784645253936e-04, 4.9633920480195821e-04, -3.7897816410764171e-06, -1.5002546732690018e-05,
        -1.5413694411715089e-04, 4.8245436976423641e-06, -9.0086007742574192e-05, 2.6103625253466520e-06, 5.2011180820413507e-05, 1.1468309313383139e-04,
        -1.5070935067046508e-06, -3.3230959625222529e-06, -6.6212314692308015e-05, 1.9185903485009741e-06, -3.9208040614284722e-04, 3.4540136347213093e-05,
        -2.2636772803059100e-04, 4.9913404823208253e-04, 1.9941757017931554e-05, -4.3970975879880263e-05, 2.8817517710850375e-04, -2.5386654760776736e-05,
        -1.0016098602569918e-03, 8.8095783986845972e-05, -1.2554426638121602e-04, 1.2750894369280275e-03, 1.0978610563435211e-05, -1.1214945864320599e-04,
        1.5982287543424081e-04, -1.3976210615569642e-05, -1.2792161864866210e-03, 1.1225748725241937e-04, -3.4731820796203604e-05, 1.6284934001327808e-03,
        2.9711553210233968e-06, -1.4290827386115594e-04, 4.4214997854704112e-05, -3.7823996304677249e-06, -1.2706869384245618e-03, 1.1125183114826646e-04,
        3.9656184460818718e-05, 1.6176353260059101e-03, -3.5517711434682683e-06, -1.4162803339381555e-04, -5.0483909874733890e-05, 4.5215468089141537e-06,
        -1.0441173998824859e-03, 9.1153799606617992e-05, 9.1153799606618073e-05, 1.3292032360397348e-03, -8.0518327772840401e-06, -1.1604243491016443e-04,
        -1.1604243491016380e-04, 1.0250305362999503e-05, -6.9005433895489503e-04, 5.9949151256281410e-05, 1.1326460393002962e-04, 8.7846678973585207e-04,
        -9.9641793477503907e-06, -7.6317668737877270e-05, -1.4419037369695810e-04, 1.2684799080065410e-05, -3.2161392308227690e-04, 2.7631503624010198e-05,
        9.9454569354364179e-05, 4.0942739519954688e-04, -8.6944232123168641e-06, -3.5176009937019723e-05, -1.2660964699911679e-04, 1.1068348703517426e-05,
        -7.4676777957819192e-05, 6.2861604388864368e-06, 4.3114657856159156e-05, 9.5066526934358737e-05, -3.6293164215592856e-06, -8.0025338133193964e-06,
        -5.4886684916474010e-05, 4.6202650513188232e-06, -2.5807022372216449e-04, 4.2830678547518969e-05, -1.4899691313582071e-04, 3.2853372287025562e-04,
        2.4728303788984804e-05, -5.4525167891608260e-05, 1.8967903333701682e-04, -3.1480120359831297e-05, -6.5979684339352498e-04, 1.0925005331179476e-04,
        -8.2940058872077706e-05, 8.3994778696934924e-04, 1.3618940110576470e-05, -1.3907969009605478e-04, 1.0558601423795538e-04, -1.7337455796106072e-05,
        -8.4357750896041372e-04, 1.3925862048674929e-04, -2.3165757864814100e-05, 1.0739079292105683e-03, 3.7065142258784857e-06, -1.7728179706443878e-04,
        2.9490936864655540e-05, -4.7185409457013005e-06, -8.3879554668975860e-04, 1.3810078013050441e-04, 2.5926625069031462e-05, 1.0678202998640096e-03,
        -4.3749936672350818e-06, -1.7580781995371169e-04, -3.3005631306617849e-05, 5.5695420273589242e-06, -6.9005433895489578e-04, 1.1326460393002904e-04,
        5.9949151256281458e-05, 8.7846678973584904e-04, -9.9641793477503653e-06, -1.4419037369695927e-04, -7.6317668737876267e-05, 1.2684799080065638e-05,
        -4.5700790975607651e-04, 7.4600267375337498e-05, 7.4600267375337959e-05, 5.8178935875592588e-04, -1.2358685770139086e-05, -9.4969125900869786e-05,
        -9.4969125900869000e-05, 1.5733101584855904e-05, -2.1405750025111815e-04, 3.4466226050632516e-05, 6.5667216952079626e-05, 2.7250376448509594e-04,
        -1.0812713792346902e-05, -4.3876885114389393e-05, -8.3596995207867956e-05, 1.3765017386719724e-05, -5.0159272049234077e-05, 7.8690281972934216e-06,
        2.8959469219980895e-05, 6.3854760712488500e-05, -4.5431855479681580e-06, -1.0017587816761052e-05, -3.6866563286394175e-05, 5.7836570226375672e-06,
        -1.1891747993226251e-04, 3.7509195538055992e-05, -6.8657039050244378e-05, 1.5138671107812595e-04, 2.1655944140983332e-05, -4.7750707052715593e-05,
        8.7403158392689037e-05, -2.7568883570880843e-05, -3.0491150826769691e-04, 9.5643090454552613e-05, -3.8726663276881581e-05, 3.8816455270336027e-04,
        1.1907675738097432e-05, -1.2175748182277376e-04, 4.9300592207778536e-05, -1.5158947764468221e-05, -3.9114247992314591e-04, 1.2195486531873341e-04,
        -1.1058811420874908e-05, 4.9794003061817552e-04, 3.2834345625941889e-06, -1.5525342423246220e-04, 1.4078309516759917e-05, -4.1799436025256102e-06,
        -3.8988499721200818e-04, 1.2107784645253757e-04, 1.1784819402651616e-05, 4.9633920480194390e-04, -3.7897816410767364e-06, -1.5413694411715395e-04,
        -1.5002546732686913e-05, 4.8245436976436923e-06, -3.2161392308228053e-04, 9.9454569354363149e-05, 2.7631503624011527e-05, 4.0942739519954103e-04,
        -8.6944232123161492e-06, -1.2660964699911720e-04, -3.5176009937017941e-05, 1.1068348703517651e-05, -2.1405750025111896e-04, 6.5667216952079626e-05,
        3.4466226050632868e-05, 2.7250376448509437e-04, -1.0812713792346967e-05, -8.3596995207868051e-05, -4.3876885114388695e-05, 1.3765017386719678e-05,
        -1.0163877852160976e-04, 3.0438753201852812e-05, 3.0438753201852917e-05, 1.2939023268193514e-04, -9.5264492369834572e-06, -3.8749750996843554e-05,
        -3.8749750996843378e-05, 1.2127551130928198e-05, -2.4458655728475410e-05, 6.9692295538358886e-06, 1.4121211468851588e-05, 3.1136847587379800e-05,
        -4.0236865589514017e-06, -8.8721081333432448e-06, -1.7976867336290102e-05, 5.1223140190652031e-06, -2.6974823754661329e-05, 1.5573921756091266e-05,
        -1.5573921756096677e-05, 3.4340030182743030e-05, 8.9916079182173179e-06, -1.9826225669985950e-05, 1.9826225669986620e-05, -1.1446676727580623e-05,
        -6.9626408103756491e-05, 4.0198825461406474e-05, -9.0509819492174502e-06, 8.8637203992336435e-05, 5.2255868648150289e-06, -5.1174713585190749e-05,
        1.1522262245218026e-05, -6.6523812089505900e-06, -8.9955228480439079e-05, 5.1935675378196955e-05, -2.6858679675677970e-06, 1.1451660589925285e-04,
        1.5506865940821917e-06, -6.6116193242616016e-05, 3.4192174122066818e-06, -1.9740860933552026e-06, -9.0086007742595225e-05, 5.2011180820405816e-05,
        2.6103625253575236e-06, 1.1468309313381435e-04, -1.5070935066976264e-06, -6.6212314692306511e-05, -3.3230959625159697e-06, 1.9185903485014870e-06,
        -7.4676777957820601e-05, 4.3114657856161548e-05, 6.2861604388869069e-06, 9.5066526934354780e-05, -3.6293164215604896e-06, -5.4886684916472533e-05,
        -8.0025338133181208e-06, 4.6202650513183065e-06, -5.0159272049234396e-05, 2.8959469219981146e-05, 7.8690281972935876e-06, 6.3854760712487497e-05,
        -4.5431855479681952e-06, -3.6866563286393782e-05, -1.0017587816760623e-05, 5.7836570226374503e-06, -2.4458655728475423e-05, 1.4121211468851677e-05,
        6.9692295538358971e-06, 3.1136847587379671e-05, -4.0236865589514525e-06, -1.7976867336290004e-05, -8.8721081333431685e-06, 5.1223140190651497e-06,
        -6.1937980258105384e-06, 3.5759909575078880e-06, 3.5759909575078880e-06, 7.8849527650921135e-06, -2.0645993419368465e-06, -4.5523796014734159e-06,
        -4.5523796014734159e-06, 2.6283175883640376e-06,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_DxDx_PyPy(void)
{
  static const double sol[512] = {
        5.0226946551362277e-03, 2.8998541112002828e-03, 2.8998541112001328e-03, -6.3940913062055585e-03, 1.6742315517121347e-03, -3.6916303368608990e-03,
        -3.6916303368608040e-03, -2.1313637687352333e-03, 5.0226946551361054e-03, 2.8998541112002867e-03, -2.8998541112002039e-03, -6.3940913062055281e-03,
        -1.6742315517121328e-03, -3.6916303368608964e-03, 3.6916303368608209e-03, 2.1313637687352342e-03, -5.0226946551365859e-03, -2.8998541112002637e-03,
        -2.8998541112003253e-03, 6.3940913062057849e-03, -1.6742315517121109e-03, 3.6916303368608938e-03, 3.6916303368609415e-03, 2.1313637687352177e-03,
        -5.0226946551364645e-03, -2.8998541112002676e-03, 2.8998541112003947e-03, 6.3940913062057719e-03, 1.6742315517121085e-03, 3.6916303368608890e-03,
        -3.6916303368609489e-03, -2.1313637687352207e-03, 5.0226946551364515e-03, 2.8998541112002685e-03, 2.8998541112002637e-03, -6.3940913062056999e-03,
        1.6742315517121265e-03, -3.6916303368608890e-03, -3.6916303368609046e-03, -2.1313637687352246e-03, 5.0226946551364463e-03, 2.8998541112002746e-03,
        -2.8998541112002672e-03, -6.3940913062057441e-03, -1.6742315517121232e-03, -3.6916303368608881e-03, 3.6916303368608790e-03, 2.1313637687352246e-03,
        -5.0226946551362259e-03, -2.8998541112002832e-03, -2.8998541112002503e-03, 6.3940913062056105e-03, -1.6742315517121248e-03, 3.6916303368609046e-03,
        3.6916303368609068e-03, 2.1313637687352290e-03, -5.0226946551362190e-03, -2.8998541112002893e-03, 2.8998541112002538e-03, 6.3940913062056392e-03,
        1.6742315517121209e-03, 3.6916303368609068e-03, -3.6916303368608916e-03, -2.1313637687352281e-03, 1.2790705487348994e-02, 1.5850090338456867e-03,
        7.3847172562461283e-03, -1.6283079974463467e-02, 9.1510539235880120e-04, -2.0177799327709654e-03, -9.4010406064926836e-03, -1.1649657873507447e-03,
        1.2790705487348885e-02, 1.5850090338456915e-03, -7.3847172562461908e-03, -1.6283079974463436e-02, -9.1510539235879838e-04, -2.0177799327709650e-03,
        9.4010406064927009e-03, 1.1649657873507447e-03, -1.2790705487349287e-02, -1.5850090338456703e-03, -7.3847172562462498e-03, 1.6283079974463682e-02,
        -9.1510539235878569e-04, 2.0177799327709607e-03, 9.4010406064927773e-03, 1.1649657873507358e-03, -1.2790705487349180e-02, -1.5850090338456742e-03,
        7.3847172562463122e-03, 1.6283079974463655e-02, 9.1510539235878320e-04, 2.0177799327709585e-03, -9.4010406064927911e-03, -1.1649657873507374e-03,
        1.2790705487349174e-02, 1.5850090338456761e-03, 7.3847172562462315e-03, -1.6283079974463582e-02, 9.1510539235879448e-04, -2.0177799327709585e-03,
        -9.4010406064927565e-03, -1.1649657873507384e-03, 1.2790705487349183e-02, 1.5850090338456796e-03, -7.3847172562462263e-03, -1.6283079974463623e-02,
        -9.1510539235879263e-04, -2.0177799327709572e-03, 9.4010406064927322e-03, 1.1649657873507393e-03, -1.2790705487348996e-02, -1.5850090338456880e-03,
        -7.3847172562462159e-03, 1.6283079974463543e-02, -9.1510539235879599e-04, 2.0177799327709711e-03, 9.4010406064927738e-03, 1.1649657873507428e-03,
        -1.2790705487349006e-02, -1.5850090338456926e-03, 7.3847172562462116e-03, 1.6283079974463575e-02, 9.1510539235879350e-04, 2.0177799327709715e-03,
        -9.4010406064927547e-03, -1.1649657873507419e-03, 1.6258847328107631e-02, 4.1732359150411819e-04, 9.3870498815959488e-03, -2.0698163334150768e-02,
        2.4094188789408792e-04, -5.3126963343913351e-04, -1.1950090172702790e-02, -3.0672866587835884e-04, 1.6258847328107544e-02, 4.1732359150412318e-04,
        -9.3870498815960008e-03, -2.0698163334150737e-02, -2.4094188789408486e-04, -5.3126963343913417e-04, 1.1950090172702807e-02, 3.0672866587835884e-04,
        -1.6258847328107867e-02, -4.1732359150410052e-04, -9.3870498815960286e-03, 2.0698163334150962e-02, -2.4094188789407700e-04, 5.3126963343912647e-04,
        1.1950090172702856e-02, 3.0672866587835233e-04, -1.6258847328107773e-02, -4.1732359150410545e-04, 9.3870498815960806e-03, 2.0698163334150928e-02,
        2.4094188789407423e-04, 5.3126963343912430e-04, -1.1950090172702876e-02, -3.0672866587835314e-04, 1.6258847328107780e-02, 4.1732359150410756e-04,
        9.3870498815960286e-03, -2.0698163334150858e-02, 2.4094188789407990e-04, -5.3126963343912723e-04, -1.1950090172702843e-02, -3.0672866587835341e-04,
        1.6258847328107794e-02, 4.1732359150410756e-04, -9.3870498815960164e-03, -2.0698163334150893e-02, -2.4094188789407995e-04, -5.3126963343912517e-04,
        1.1950090172702823e-02, 3.0672866587835455e-04, -1.6258847328107638e-02, -4.1732359150411965e-04, -9.3870498815960199e-03, 2.0698163334150855e-02,
        -2.4094188789408291e-04, 5.3126963343913568e-04, 1.1950090172702871e-02, 3.0672866587835650e-04, -1.6258847328107656e-02, -4.1732359150412036e-04,
        9.3870498815960095e-03, 2.0698163334150886e-02, 2.4094188789408158e-04, 5.3126963343913492e-04, -1.1950090172702854e-02, -3.0672866587835694e-04,
        1.6064983708873892e-02, -5.2925080425479339e-04, 9.2751226688452867e-03, -2.0451367188367861e-02, -3.0556309430532816e-04, 6.7375745464185800e-04,
        -1.1807602351500069e-02, 3.8899404780599303e-04, 1.6064983708873822e-02, -5.2925080425478700e-04, -9.2751226688453266e-03, -2.0451367188367833e-02,
        3.0556309430533130e-04, 6.7375745464185615e-04, 1.1807602351500086e-02, -3.8899404780599412e-04, -1.6064983708874072e-02, 5.2925080425480684e-04,
        -9.2751226688453318e-03, 2.0451367188368024e-02, 3.0556309430533526e-04, -6.7375745464186689e-04, 1.1807602351500116e-02, -3.8899404780599748e-04,
        -1.6064983708873996e-02, 5.2925080425480402e-04, 9.2751226688453751e-03, 2.0451367188367989e-02, -3.0556309430533781e-04, -6.7375745464186602e-04,
        -1.1807602351500135e-02, 3.8899404780599791e-04, 1.6064983708874003e-02, -5.2925080425480196e-04, 9.2751226688453387e-03, -2.0451367188367930e-02,
        -3.0556309430533493e-04, 6.7375745464186396e-04, -1.1807602351500104e-02, 3.8899404780599742e-04, 1.6064983708874016e-02, -5.2925080425480337e-04,
        -9.2751226688453300e-03, -2.0451367188367958e-02, 3.0556309430533374e-04, 6.7375745464186678e-04, 1.1807602351500088e-02, -3.8899404780599602e-04,
        -1.6064983708873902e-02, 5.2925080425479133e-04, -9.2751226688453405e-03, 2.0451367188367951e-02, 3.0556309430533249e-04, -6.7375745464185865e-04,
        1.1807602351500142e-02, -3.8899404780599656e-04, -1.6064983708873923e-02, 5.2925080425479058e-04, 9.2751226688453266e-03, 2.0451367188367976e-02,
        -3.0556309430533211e-04, -6.7375745464186027e-04, -1.1807602351500125e-02, 3.8899404780599553e-04, 1.3108851408420495e-02, -1.1774729751721242e-03,
        7.5683988894183754e-03, -1.6688092464312451e-02, -6.7981433917913405e-04, 1.4989702203260478e-03, -9.6348746765321652e-03, 8.6543086021247605e-04,
        1.3108851408420441e-02, -1.1774729751721192e-03, -7.5683988894184049e-03, -1.6688092464312427e-02, 6.7981433917913730e-04, 1.4989702203260455e-03,
        9.6348746765321774e-03, -8.6543086021247735e-04, -1.3108851408420628e-02, 1.1774729751721368e-03, -7.5683988894183997e-03, 1.6688092464312576e-02,
        6.7981433917913838e-04, -1.4989702203260585e-03, 9.6348746765321982e-03, -8.6543086021248039e-04, -1.3108851408420566e-02, 1.1774729751721311e-03,
        7.5683988894184361e-03, 1.6688092464312545e-02, -6.7981433917914142e-04, -1.4989702203260570e-03, -9.6348746765322138e-03, 8.6543086021248104e-04,
        1.3108851408420576e-02, -1.1774729751721318e-03, 7.5683988894184092e-03, -1.6688092464312500e-02, -6.7981433917913817e-04, 1.4989702203260522e-03,
        -9.6348746765321843e-03, 8.6543086021247909e-04, 1.3108851408420587e-02, -1.1774729751721333e-03, -7.5683988894184040e-03, -1.6688092464312517e-02,
        6.7981433917913730e-04, 1.4989702203260546e-03, 9.6348746765321774e-03, -8.6543086021247779e-04, -1.3108851408420507e-02, 1.1774729751721248e-03,
        -7.5683988894184135e-03, 1.6688092464312531e-02, 6.7981433917913882e-04, -1.4989702203260522e-03, 9.6348746765322225e-03, -8.6543086021248115e-04,
        -1.3108851408420528e-02, 1.1774729751721261e-03, 7.5683988894184023e-03, 1.6688092464312552e-02, -6.7981433917913817e-04, -1.4989702203260537e-03,
        -9.6348746765322086e-03, 8.6543086021248017e-04, 8.5642422695473425e-03, -1.4463583345179229e-03, 4.9445675797283326e-03, -1.0902623153480380e-02,
        -8.3505537377857994e-04, 1.8412720436711647e-03, -6.2946324125349529e-03, 1.0630589100648789e-03, 8.5642422695473095e-03, -1.4463583345179173e-03,
        -4.9445675797283508e-03, -1.0902623153480369e-02, 8.3505537377858319e-04, 1.8412720436711617e-03, 6.2946324125349607e-03, -1.0630589100648806e-03,
        -8.5642422695474396e-03, 1.4463583345179336e-03, -4.9445675797283465e-03, 1.0902623153480473e-02, 8.3505537377858200e-04, -1.8412720436711751e-03,
        6.2946324125349737e-03, -1.0630589100648819e-03, -8.5642422695473945e-03, 1.4463583345179286e-03, 4.9445675797283725e-03, 1.0902623153480451e-02,
        -8.3505537377858482e-04, -1.8412720436711725e-03, -6.2946324125349859e-03, 1.0630589100648832e-03, 8.5642422695474049e-03, -1.4463583345179279e-03,
        4.9445675797283543e-03, -1.0902623153480420e-02, -8.3505537377858363e-04, 1.8412720436711679e-03, -6.2946324125349659e-03, 1.0630589100648815e-03,
        8.5642422695474084e-03, -1.4463583345179294e-03, -4.9445675797283508e-03, -1.0902623153480428e-02, 8.3505537377858276e-04, 1.8412720436711697e-03,
        6.2946324125349616e-03, -1.0630589100648802e-03, -8.5642422695473546e-03, 1.4463583345179234e-03, -4.9445675797283560e-03, 1.0902623153480446e-02,
        8.3505537377858406e-04, -1.8412720436711708e-03, 6.2946324125349911e-03, -1.0630589100648848e-03, -8.5642422695473720e-03, 1.4463583345179244e-03,
        4.9445675797283465e-03, 1.0902623153480463e-02, -8.3505537377858352e-04, -1.8412720436711725e-03, -6.2946324125349815e-03, 1.0630589100648841e-03,
        3.8932555050064832e-03, -1.2504371313709222e-03, 2.2477721138394924e-03, -4.9562700674910002e-03, -7.2194021440170409e-04, 1.5918565112213662e-03,
        -2.8615038576424224e-03, 9.1905878526491434e-04, 3.8932555050064706e-03, -1.2504371313709168e-03, -2.2477721138395002e-03, -4.9562700674909959e-03,
        7.2194021440170734e-04, 1.5918565112213634e-03, 2.8615038576424237e-03, -9.1905878526491597e-04, -3.8932555050065431e-03, 1.2504371313709329e-03,
        -2.2477721138394993e-03, 4.9562700674910557e-03, 7.2194021440170637e-04, -1.5918565112213766e-03, 2.8615038576424324e-03, -9.1905878526491694e-04,
        -3.8932555050065162e-03, 1.2504371313709277e-03, 2.2477721138395149e-03, 4.9562700674910419e-03, -7.2194021440170908e-04, -1.5918565112213743e-03,
        -2.8615038576424411e-03, 9.1905878526491824e-04, 3.8932555050065266e-03, -1.2504371313709268e-03, 2.2477721138395011e-03, -4.9562700674910297e-03,
        -7.2194021440170745e-04, 1.5918565112213688e-03, -2.8615038576424263e-03, 9.1905878526491640e-04, 3.8932555050065266e-03, -1.2504371313709283e-03,
        -2.2477721138395011e-03, -4.9562700674910315e-03, 7.2194021440170680e-04, 1.5918565112213706e-03, 2.8615038576424263e-03, -9.1905878526491532e-04,
        -3.8932555050064949e-03, 1.2504371313709220e-03, -2.2477721138395028e-03, 4.9562700674910410e-03, 7.2194021440170723e-04, -1.5918565112213730e-03,
        2.8615038576424424e-03, -9.1905878526491922e-04, -3.8932555050065088e-03, 1.2504371313709229e-03, 2.2477721138394937e-03, 4.9562700674910531e-03,
        -7.2194021440170680e-04, -1.5918565112213738e-03, -2.8615038576424354e-03, 9.1905878526491857e-04, 8.6371743090067675e-04, -4.9866749123428583e-04,
        4.9866749123428833e-04, -1.0995468558481246e-03, -2.8790581030022907e-04, 6.3482367321052548e-04, -6.3482367321052939e-04, 3.6651561861604439e-04,
        8.6371743090068336e-04, -4.9866749123428031e-04, -4.9866749123428453e-04, -1.0995468558481305e-03, 2.8790581030023222e-04, 6.3482367321052299e-04,
        6.3482367321052603e-04, -3.6651561861604585e-04, -8.6371743090070082e-04, 4.9866749123429613e-04, -4.9866749123428887e-04, 1.0995468558481442e-03,
        2.8790581030023097e-04, -6.3482367321053622e-04, 6.3482367321053134e-04, -3.6651561861604666e-04, -8.6371743090069041e-04, 4.9866749123429136e-04,
        4.9866749123429473e-04, 1.0995468558481370e-03, -2.8790581030023368e-04, -6.3482367321053416e-04, -6.3482367321053546e-04, 3.6651561861604786e-04,
        8.6371743090070527e-04, -4.9866749123429028e-04, 4.9866749123428529e-04, -1.0995468558481448e-03, -2.8790581030023249e-04, 6.3482367321052863e-04,
        -6.3482367321052613e-04, 3.6651561861604677e-04, 8.6371743090070006e-04, -4.9866749123429147e-04, -4.9866749123428833e-04, -1.0995468558481398e-03,
        2.8790581030023173e-04, 6.3482367321053036e-04, 6.3482367321052884e-04, -3.6651561861604590e-04, -8.6371743090068976e-04, 4.9866749123428540e-04,
        -4.9866749123428844e-04, 1.0995468558481435e-03, 2.8790581030023178e-04, -6.3482367321053210e-04, 6.3482367321053361e-04, -3.6651561861604856e-04,
        -8.6371743090070180e-04, 4.9866749123428594e-04, 4.9866749123428161e-04, 1.0995468558481518e-03, -2.8790581030023140e-04, -6.3482367321053329e-04,
        -6.3482367321052884e-04, 3.6651561861604791e-04,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_PxPx_DyDy(void)
{
  static const double sol[512] = {
        5.0226946551374047e-03, 2.8998541112002477e-03, 2.8998541112002247e-03, -6.3940913062064493e-03, 1.6742315517121274e-03, -3.6916303368605963e-03,
        -3.6916303368608513e-03, -2.1313637687352671e-03, 1.2790705487349972e-02, 7.3847172562462203e-03, 1.5850090338456310e-03, -1.6283079974464206e-02,
        9.1510539235879567e-04, -9.4010406064925674e-03, -2.0177799327709238e-03, -1.1649657873507653e-03, 1.6258847328108415e-02, 9.3870498815960199e-03,
        4.1732359150406154e-04, -2.0698163334151354e-02, 2.4094188789408060e-04, -1.1950090172702729e-02, -5.3126963343909015e-04, -3.0672866587836957e-04,
        1.6064983708874485e-02, 9.2751226688453300e-03, -5.2925080425484793e-04, -2.0451367188368302e-02, -3.0556309430533428e-04, -1.1807602351500038e-02,
        6.7375745464189811e-04, 3.8899404780598740e-04, 1.3108851408420902e-02, 7.5683988894183979e-03, -1.1774729751721765e-03, -1.6688092464312750e-02,
        -6.7981433917913892e-04, -9.6348746765321513e-03, 1.4989702203260884e-03, 8.6543086021247096e-04, 8.5642422695475732e-03, 4.9445675797283361e-03,
        -1.4463583345179739e-03, -1.0902623153480544e-02, -8.3505537377858644e-04, -6.2946324125349503e-03, 1.8412720436712046e-03, 1.0630589100648776e-03,
        3.8932555050065370e-03, 2.2477721138394720e-03, -1.2504371313709736e-03, -4.9562700674910228e-03, -7.2194021440171125e-04, -2.8615038576424259e-03,
        1.5918565112214074e-03, 9.1905878526491250e-04, 8.6371743090054871e-04, 4.9866749123423770e-04, -4.9866749123433972e-04, -1.0995468558480034e-03,
        -2.8790581030023926e-04, -6.3482367321054012e-04, 6.3482367321056777e-04, 3.6651561861604206e-04, 5.0226946551373657e-03, -2.8998541112002702e-03,
        2.8998541112002286e-03, -6.3940913062062463e-03, -1.6742315517121252e-03, 3.6916303368607138e-03, -3.6916303368608734e-03, 2.1313637687352550e-03,
        1.2790705487349948e-02, -7.3847172562462315e-03, 1.5850090338456367e-03, -1.6283079974464067e-02, -9.1510539235879220e-04, 9.4010406064926489e-03,
        -2.0177799327709398e-03, 1.1649657873507560e-03, 1.6258847328108415e-02, -9.3870498815960199e-03, 4.1732359150406864e-04, -2.0698163334151271e-02,
        -2.4094188789407586e-04, 1.1950090172702777e-02, -5.3126963343910468e-04, 3.0672866587836117e-04, 1.6064983708874509e-02, -9.2751226688453179e-03,
        -5.2925080425484088e-04, -2.0451367188368264e-02, 3.0556309430533862e-04, 1.1807602351500059e-02, 6.7375745464188727e-04, -3.8899404780599363e-04,
        1.3108851408420951e-02, -7.5683988894183702e-03, -1.1774729751721673e-03, -1.6688092464312750e-02, 6.7981433917914467e-04, 9.6348746765321531e-03,
        1.4989702203260782e-03, -8.6543086021247681e-04, 8.5642422695476547e-03, -4.9445675797282892e-03, -1.4463583345179648e-03, -1.0902623153480574e-02,
        8.3505537377859176e-04, 6.2946324125349338e-03, 1.8412720436711953e-03, -1.0630589100648830e-03, 3.8932555050066524e-03, -2.2477721138394052e-03,
        -1.2504371313709630e-03, -4.9562700674910896e-03, 7.2194021440171732e-04, 2.8615038576423877e-03, 1.5918565112213962e-03, -9.1905878526491900e-04,
        8.6371743090070375e-04, -4.9866749123414814e-04, -4.9866749123432703e-04, -1.0995468558481086e-03, 2.8790581030024653e-04, 6.3482367321047940e-04,
        6.3482367321055606e-04, -3.6651561861604878e-04, -5.0226946551355815e-03, -2.8998541112004064e-03, -2.8998541112002945e-03, 6.3940913062051005e-03,
        -1.6742315517121022e-03, 3.6916303368610673e-03, 3.6916303368609120e-03, 2.1313637687351951e-03, -1.2790705487348395e-02, -7.3847172562463070e-03,
        -1.5850090338457032e-03, 1.6283079974463061e-02, -9.1510539235877886e-04, 9.4010406064928363e-03, 2.0177799327709784e-03, 1.1649657873507200e-03,
        -1.6258847328107093e-02, -9.3870498815960667e-03, -4.1732359150413803e-04, 2.0698163334150407e-02, -2.4094188789407412e-04, 1.1950090172702878e-02,
        5.3126963343915065e-04, 3.0672866587834696e-04, -1.6064983708873427e-02, -9.2751226688453613e-03, 5.2925080425477149e-04, 2.0451367188367559e-02,
        3.0556309430533705e-04, 1.1807602351500126e-02, -6.7375745464184141e-04, -3.8899404780600122e-04, -1.3108851408420103e-02, -7.5683988894184265e-03,
        1.1774729751721036e-03, 1.6688092464312205e-02, 6.7981433917913914e-04, 9.6348746765321912e-03, -1.4989702203260322e-03, -8.6543086021248473e-04,
        -8.5642422695470250e-03, -4.9445675797283699e-03, 1.4463583345179017e-03, 1.0902623153480194e-02, 8.3505537377858319e-04, 6.2946324125349572e-03,
        -1.8412720436711478e-03, -1.0630589100648835e-03, -3.8932555050062343e-03, -2.2477721138395206e-03, 1.2504371313709038e-03, 4.9562700674908649e-03,
        7.2194021440170626e-04, 2.8615038576424111e-03, -1.5918565112213521e-03, -9.1905878526491868e-04, -8.6371743090048723e-04, -4.9866749123431348e-04,
        4.9866749123427011e-04, 1.0995468558480334e-03, 2.8790581030022875e-04, 6.3482367321050141e-04, -6.3482367321051421e-04, -3.6651561861605013e-04,
        -5.0226946551355824e-03, 2.8998541112004064e-03, -2.8998541112002976e-03, 6.3940913062052853e-03, 1.6742315517121005e-03, -3.6916303368609619e-03,
        3.6916303368608946e-03, -2.1313637687352047e-03, -1.2790705487348402e-02, 7.3847172562463009e-03, -1.5850090338457047e-03, 1.6283079974463197e-02,
        9.1510539235877832e-04, -9.4010406064927565e-03, 2.0177799327709676e-03, -1.1649657873507263e-03, -1.6258847328107107e-02, 9.3870498815960580e-03,
        -4.1732359150413803e-04, 2.0698163334150515e-02, 2.4094188789407342e-04, -1.1950090172702817e-02, 5.3126963343914230e-04, -3.0672866587835125e-04,
        -1.6064983708873441e-02, 9.2751226688453526e-03, 5.2925080425477008e-04, 2.0451367188367642e-02, -3.0556309430533781e-04, -1.1807602351500076e-02,
        -6.7375745464184553e-04, 3.8899404780599835e-04, -1.3108851408420125e-02, 7.5683988894184127e-03, 1.1774729751720999e-03, 1.6688092464312271e-02,
        -6.7981433917914098e-04, -9.6348746765321531e-03, -1.4989702203260372e-03, 8.6543086021248158e-04, -8.5642422695470614e-03, 4.9445675797283499e-03,
        1.4463583345178976e-03, 1.0902623153480243e-02, -8.3505537377858601e-04, -6.2946324125349286e-03, -1.8412720436711528e-03, 1.0630589100648804e-03,
        -3.8932555050062894e-03, 2.2477721138394889e-03, 1.2504371313708973e-03, 4.9562700674909014e-03, -7.2194021440171005e-04, -2.8615038576423912e-03,
        -1.5918565112213552e-03, 9.1905878526491705e-04, -8.6371743090057050e-04, 4.9866749123426534e-04, 4.9866749123426046e-04, 1.0995468558480613e-03,
        -2.8790581030023438e-04, -6.3482367321048537e-04, -6.3482367321051594e-04, 3.6651561861604910e-04, 5.0226946551371540e-03, 2.8998541112001484e-03,
        2.8998541112002308e-03, -6.3940913062060711e-03, 1.6742315517121297e-03, -3.6916303368608261e-03, -3.6916303368608634e-03, -2.1313637687352298e-03,
        1.2790705487349752e-02, 7.3847172562461301e-03, 1.5850090338456432e-03, -1.6283079974463873e-02, 9.1510539235879914e-04, -9.4010406064926923e-03,
        -2.0177799327709372e-03, -1.1649657873507423e-03, 1.6258847328108252e-02, 9.3870498815959418e-03, 4.1732359150407997e-04, -2.0698163334151087e-02,
        2.4094188789408404e-04, -1.1950090172702791e-02, -5.3126963343911292e-04, -3.0672866587835694e-04, 1.6064983708874391e-02, 9.2751226688452693e-03,
        -5.2925080425482386e-04, -2.0451367188368114e-02, -3.0556309430532870e-04, -1.1807602351500069e-02, 6.7375745464187588e-04, 3.8899404780599265e-04,
        1.3108851408420885e-02, 7.5683988894183589e-03, -1.1774729751721552e-03, -1.6688092464312639e-02, -6.7981433917913383e-04, -9.6348746765321687e-03,
        1.4989702203260650e-03, 8.6543086021247508e-04, 8.5642422695476304e-03, 4.9445675797283161e-03, -1.4463583345179526e-03, -1.0902623153480515e-02,
        -8.3505537377857994e-04, -6.2946324125349555e-03, 1.8412720436711812e-03, 1.0630589100648806e-03, 3.8932555050066667e-03, 2.2477721138394777e-03,
        -1.2504371313709515e-03, -4.9562700674910783e-03, -7.2194021440170301e-04, -2.8615038576424228e-03, 1.5918565112213816e-03, 9.1905878526491402e-04,
        8.6371743090075666e-04, 4.9866749123427510e-04, -4.9866749123431652e-04, -1.0995468558481457e-03, -2.8790581030022923e-04, -6.3482367321053286e-04,
        6.3482367321054283e-04, 3.6651561861604319e-04, 5.0226946551366995e-03, -2.8998541112004112e-03, 2.8998541112002733e-03, -6.3940913062058456e-03,
        -1.6742315517121052e-03, 3.6916303368609562e-03, -3.6916303368608825e-03, 2.1313637687352185e-03, 1.2790705487349422e-02, -7.3847172562463209e-03,
        1.5850090338456729e-03, -1.6283079974463713e-02, -9.1510539235878201e-04, 9.4010406064927877e-03, -2.0177799327709533e-03, 1.1649657873507328e-03,
        1.6258847328108013e-02, -9.3870498815960806e-03, 4.1732359150410263e-04, -2.0698163334150969e-02, -2.4094188789407177e-04, 1.1950090172702859e-02,
        -5.3126963343912354e-04, 3.0672866587835081e-04, 1.6064983708874218e-02, -9.2751226688453682e-03, -5.2925080425480825e-04, -2.0451367188368035e-02,
        3.0556309430533862e-04, 1.1807602351500116e-02, 6.7375745464186602e-04, -3.8899404780599835e-04, 1.3108851408420767e-02, -7.5683988894184265e-03,
        -1.1774729751721411e-03, -1.6688092464312590e-02, 6.7981433917914185e-04, 9.6348746765321965e-03, 1.4989702203260583e-03, -8.6543086021247898e-04,
        8.5642422695475610e-03, -4.9445675797283551e-03, -1.4463583345179385e-03, -1.0902623153480489e-02, 8.3505537377858807e-04, 6.2946324125349711e-03,
        1.8412720436711760e-03, -1.0630589100648837e-03, 3.8932555050066476e-03, -2.2477721138394889e-03, -1.2504371313709372e-03, -4.9562700674910679e-03,
        7.2194021440171125e-04, 2.8615038576424285e-03, 1.5918565112213771e-03, -9.1905878526491651e-04, 8.6371743090078821e-04, -4.9866749123425689e-04,
        -4.9866749123430134e-04, -1.0995468558481546e-03, 2.8790581030023796e-04, 6.3482367321052776e-04, 6.3482367321053687e-04, -3.6651561861604666e-04,
        -5.0226946551359909e-03, -2.8998541112001167e-03, -2.8998541112002759e-03, 6.3940913062055568e-03, -1.6742315517121373e-03, 3.6916303368608578e-03,
        3.6916303368608851e-03, 2.1313637687352211e-03, -1.2790705487348746e-02, -7.3847172562461188e-03, -1.5850090338456876e-03, 1.6283079974463436e-02,
        -9.1510539235880228e-04, 9.4010406064927061e-03, 2.0177799327709607e-03, 1.1649657873507406e-03, -1.6258847328107395e-02, -9.3870498815959453e-03,
        -4.1732359150412247e-04, 2.0698163334150727e-02, -2.4094188789409087e-04, 1.1950090172702807e-02, 5.3126963343913221e-04, 3.0672866587835851e-04,
        -1.6064983708873673e-02, -9.2751226688452884e-03, 5.2925080425478559e-04, 2.0451367188367823e-02, 3.0556309430532767e-04, 1.1807602351500090e-02,
        -6.7375745464185247e-04, -3.8899404780598913e-04, -1.3108851408420302e-02, -7.5683988894183780e-03, 1.1774729751721164e-03, 1.6688092464312423e-02,
        6.7981433917913426e-04, 9.6348746765321878e-03, -1.4989702203260481e-03, -8.6543086021247952e-04, -8.5642422695471846e-03, -4.9445675797283413e-03,
        1.4463583345179110e-03, 1.0902623153480366e-02, 8.3505537377857669e-04, 6.2946324125349694e-03, -1.8412720436711595e-03, -1.0630589100648785e-03,
        -3.8932555050063661e-03, -2.2477721138395106e-03, 1.2504371313709108e-03, 4.9562700674909976e-03, 7.2194021440170138e-04, 2.8615038576424359e-03,
        -1.5918565112213645e-03, -9.1905878526491716e-04, -8.6371743090059934e-04, -4.9866749123431847e-04, 4.9866749123427478e-04, 1.0995468558481266e-03,
        2.8790581030022511e-04, 6.3482367321053329e-04, -6.3482367321052462e-04, -3.6651561861604710e-04, -5.0226946551354965e-03, 2.8998541112004025e-03,
        -2.8998541112003192e-03, 6.3940913062049453e-03, 1.6742315517121124e-03, -3.6916303368612121e-03, 3.6916303368609424e-03, -2.1313637687351882e-03,
        -1.2790705487348383e-02, 7.3847172562463270e-03, -1.5850090338457210e-03, 1.6283079974462999e-02, 9.1510539235878320e-04, -9.4010406064929612e-03,
        2.0177799327710032e-03, -1.1649657873507161e-03, -1.6258847328107139e-02, 9.3870498815960927e-03, -4.1732359150415077e-04, 2.0698163334150421e-02,
        2.4094188789407423e-04, -1.1950090172702984e-02, 5.3126963343916604e-04, -3.0672866587833851e-04, -1.6064983708873506e-02, 9.2751226688453838e-03,
        5.2925080425476434e-04, 2.0451367188367625e-02, -3.0556309430534030e-04, -1.1807602351500206e-02, -6.7375745464182590e-04, 3.8899404780600453e-04,
        -1.3108851408420209e-02, 7.5683988894184318e-03, 1.1774729751720971e-03, 1.6688092464312309e-02, -6.7981433917914597e-04, -9.6348746765322554e-03,
        -1.4989702203260255e-03, 8.6543086021249242e-04, -8.5642422695471621e-03, 4.9445675797283543e-03, 1.4463583345178917e-03, 1.0902623153480321e-02,
        -8.3505537377858807e-04, -6.2946324125349954e-03, -1.8412720436711400e-03, 1.0630589100648895e-03, -3.8932555050064065e-03, 2.2477721138394872e-03,
        1.2504371313708923e-03, 4.9562700674910176e-03, -7.2194021440171179e-04, -2.8615038576424237e-03, -1.5918565112213459e-03, 9.1905878526492789e-04,
        -8.6371743090070266e-04, 4.9866749123425895e-04, 4.9866749123425667e-04, 1.0995468558482129e-03, -2.8790581030023547e-04, -6.3482367321048342e-04,
        -6.3482367321050521e-04, 3.6651561861605837e-04,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_DxDx_NyDy(void)
{
  static const double sol[512] = {
        -6.9079330390629801e-04, -3.9882969996476189e-04, 1.9369084511113992e-05, 8.7940752169275444e-04, 1.1182746156426905e-05, 5.0772616937670115e-04,
        -2.4657619741027554e-05, -1.4236083395042380e-05, -5.6185672699465936e-04, -3.2438813257638433e-04, 5.5072482877209447e-05, 7.1526609919163396e-04,
        3.1796112814150550e-05, 4.1295907491048809e-04, -7.0109474725162422e-05, -4.0477724105344443e-05, -3.2502275688893505e-04, -1.8765197618263087e-04,
        8.1663673516578359e-05, 4.1376697705834621e-04, 4.7148543887737761e-05, 2.3888847558643323e-04, -1.0396112459892713e-04, -6.0021983272398073e-05,
        -2.2431527513896648e-05, -1.2950848448637168e-05, 9.3037454217519903e-05, 2.8556232243850333e-05, 5.3715199237214681e-05, 1.6486948373111031e-05,
        -1.1844040261445826e-04, -6.8381598332397376e-05, 2.8251817425953193e-04, 1.6311194395994056e-04, 8.3025338190665209e-05, -3.5965694232108491e-04,
        4.7934701354154866e-05, -2.0764803246514594e-04, -1.0569457822356623e-04, -6.1022793189381905e-05, 5.0481002617883098e-04, 2.9145220450425613e-04,
        4.5314922353602512e-05, -6.4264336602178908e-04, 2.6162582618972440e-05, -3.7103032036580117e-04, -5.7687709677044810e-05, -3.3306014710878311e-05,
        5.3693018033440596e-04, 3.0999678415207708e-04, -2.6770342705469010e-05, -6.8353360772284585e-04, -1.5455864567294547e-05, -3.9463831241894835e-04,
        3.4079717623710872e-05, 1.9675934143967217e-05, 2.4528129331624464e-04, 1.4161322072328672e-04, -1.4161322072332120e-04, -3.1225290264544673e-04,
        -8.1760431105425965e-05, -1.8027929739759170e-04, 1.8027929739763954e-04, 1.0408430088183608e-04, -1.7637269696838631e-03, -2.2062884079452740e-04,
        4.9468484637507319e-05, 2.2452950174551758e-03, 6.1951506089924524e-06, 2.8086934398447670e-04, -6.2975360691777362e-05, -7.8866746576113936e-06,
        -1.4344421799159486e-03, -1.7939931355867353e-04, 1.4064417737504758e-04, 1.8261023019737138e-03, 1.7608728039188745e-05, 2.2838250579132761e-04,
        -1.7904566643377488e-04, -2.2416615502115621e-05, -8.2967964857643035e-04, -1.0371181608154435e-04, 2.0851563288962292e-04, 1.0562153967439481e-03,
        2.6089469003520896e-05, 1.3202929245949860e-04, -2.6544874554619840e-04, -3.3212938152301273e-05, -5.7211901773873569e-05, -7.1296099928764446e-06,
        2.3746882866713706e-04, 7.2833040600970948e-05, 2.9672293680699401e-05, 9.0762788507300187e-06, -3.0230732248927766e-04, -3.7774017352396973e-05,
        7.2082558778282157e-04, 8.9944959225539376e-05, 2.1173132536804161e-04, -9.1763982097132471e-04, 2.6373734989019922e-05, -1.1450353272677223e-04,
        -2.6954245076450602e-04, -3.3574820128277174e-05, 1.2870709745399684e-03, 1.6018636460857875e-04, 1.1519060110358107e-04, -1.6384928596764488e-03,
        1.4180159317135439e-05, -2.0392365286806921e-04, -1.4664224517805719e-04, -1.8051910306269282e-05, 1.3669669930979470e-03, 1.6922519313425260e-04,
        -6.9062613281952698e-05, -1.7402036887708442e-03, -8.9615892346788474e-06, -2.1543044331873715e-04, 8.7919470620908099e-05, 1.1408461742067528e-05,
        6.2367351899503742e-04, 7.6851632631626000e-05, -3.6007807407825444e-04, -7.9396134934042945e-04, -4.4370310787516813e-05, -9.7835203972673063e-05,
        4.5839379876786319e-04, 5.6485181349801060e-05, -2.2508116665430275e-03, -6.0589640055258542e-05, 6.3151521941853120e-05, 2.8653733298781617e-03,
        1.7047546620032634e-06, 7.7133036611583063e-05, -8.0394414780745167e-05, -2.1702209096828419e-06, -1.8304802004777839e-03, -4.9253344222028937e-05,
        1.7952696313978540e-04, 2.3302745517463802e-03, 4.8402587889251854e-06, 6.2701478332860363e-05, -2.2854500881667886e-04, -6.1618431473634482e-06,
        -1.0586408004009461e-03, -2.8478966558306632e-05, 2.6609472226573651e-04, 1.3476921061319536e-03, 7.1538337474718335e-06, 3.6254864168172705e-05,
        -3.3874923065981108e-04, -9.1071166597731378e-06, -7.3130563309161035e-05, -2.0610335296603170e-06, 3.0288987840833758e-04, 9.3098133806493509e-05,
        8.0985669980506337e-06, 2.6237781666312881e-06, -3.8559093698597503e-04, -1.0309799896356255e-05, 9.1873938193641036e-04, 2.4320623105375322e-05,
        2.6976650148029252e-04, -1.1695920015167591e-03, 7.1328895618554546e-06, -3.0961126534047540e-05, -3.4342355254598974e-04, -9.0804538732886023e-06,
        1.6392548325407505e-03, 4.3147080612413811e-05, 1.4622328788140799e-04, -2.0868370054480294e-03, 3.7365707477203179e-06, -5.4927960382745753e-05,
        -1.8614809738656867e-04, -4.7568041008803349e-06, 1.7385437019493174e-03, 4.5304719745681212e-05, -8.8898832400838445e-05, -2.2132357097847683e-03,
        -2.4908572133142122e-06, -5.7674721348964573e-05, 1.1317177141253271e-04, 3.1709609176343877e-06, 7.9228320374875242e-04, 2.0495214248909324e-05,
        -4.5742492095876726e-04, -1.0086082258578600e-03, -1.1832917463708485e-05, -2.6091227965377788e-05, 5.8232023070588939e-04, 1.5063777489328235e-05,
        -2.2340461539592346e-03, 7.0269213258613311e-05, 6.2692883129639407e-05, 2.8440301613965667e-03, -1.9695499036957175e-06, -8.9455520679791617e-05,
        -7.9810549217891458e-05, 2.5073158495669902e-06, -1.8168031392969803e-03, 5.7149799176488832e-05, 1.7820248370310410e-04, 2.3128631055018322e-03,
        -5.6049473482291679e-06, -7.2753981509130874e-05, -2.2685889348759451e-04, 7.1353222864950245e-06, -1.0507771170938870e-03, 3.3019066232459204e-05,
        2.6406284642211032e-04, 1.3376813225743581e-03, -8.3269378127493067e-06, -4.2034592749949409e-05, -3.3616257139440161e-04, 1.0600525082959169e-05,
        -7.3041740691808134e-05, 2.1123152916991204e-06, 3.0043293633984995e-04, 9.2985059059888983e-05, -9.5170831626752805e-06, -2.6890619020219609e-06,
        -3.8246315140499850e-04, 1.2115627743500475e-05, 9.1036543165672220e-04, -2.9155325553533829e-05, 2.6733745933858288e-04, -1.1589316276819181e-03,
        -8.5352976962444854e-06, 3.7115896237255145e-05, -3.4033128468836796e-04, 1.0865775553292953e-05, 1.6238014988128172e-03, -5.2069067000113288e-05,
        1.4456504608357674e-04, -2.0671642931638563e-03, -4.6939570960796462e-06, 6.6286006115703452e-05, -1.8403708921444040e-04, 5.9755952373201247e-06,
        1.7210227955432630e-03, -5.5420419775662968e-05, -8.8434304245289404e-05, -2.1909308947362749e-03, 2.7590526689666105e-06, 7.0552412321437028e-05,
        1.1258040847991929e-04, -3.5123844659675123e-06, 7.8392504372920872e-04, -2.5320800186122937e-05, -4.5259933502154824e-04, -9.9796795365608205e-04,
        1.4618970803559581e-05, 3.2234391985325116e-05, 5.7617706668597769e-04, -1.8610534889900269e-05, -1.8339090424506843e-03, 1.6075005578360761e-04,
        5.1471540987780427e-05, 2.3346396048014714e-03, -4.5090950025754530e-06, -2.0464125429304014e-04, -6.5525331588775609e-05, 5.7402583940354572e-06,
        -1.4913852692680597e-03, 1.3073029568382339e-04, 1.4628465166305414e-04, 1.8985931336037614e-03, -1.2822821572042931e-05, -1.6642489828339903e-04,
        -1.8622621593640859e-04, 1.6323965035577423e-05, -8.6268379122426020e-04, 7.5576666091141834e-05, 2.1669631592544916e-04, 1.0982309911733820e-03,
        -1.9020141320076695e-05, -9.6212120541949218e-05, -2.7586308244495619e-04, 2.4213401094000003e-05, -6.0575372069967455e-05, 5.0851459879377340e-06,
        2.4640119579254699e-04, 7.7114872895307326e-05, -2.1678156787094102e-05, -6.4735943521887308e-06, -3.1367858331675602e-04, 2.7597161158338475e-05,
        7.4559480183159898e-04, -6.5975041930545276e-05, 2.1904137112216987e-04, -9.4917202172906947e-04, -1.9348461836307453e-05, 8.3988868724731795e-05,
        -2.7884843156041284e-04, 2.4631366250678167e-05, 1.3297080252287668e-03, -1.1772587914055043e-04, 1.1819655564340515e-04, -1.6927715315547805e-03,
        -1.0529897957677521e-05, 1.4986975558194376e-04, -1.5046894560673518e-04, 1.3404981510784434e-05, 1.4085793923499623e-03, -1.2496886316451210e-04,
        -7.2660150604629122e-05, -1.7931779383632163e-03, 6.3481591817070706e-06, 1.5909036411151771e-04, 9.2499279607521743e-05, -8.0814606941427644e-06,
        6.4136415990857299e-04, -5.6986764463632038e-05, -3.7029177037178836e-04, -8.1648224320973327e-04, 3.2901323803321743e-05, 7.2546431794944845e-05,
        4.7139624290569047e-04, -4.1884701925559310e-05, -1.2100112667458897e-03, 1.9945749296636440e-04, 3.3972790274495772e-05, 1.5403927676947034e-03,
        -5.5938134322215696e-06, -2.5391737091356098e-04, -4.3248721623872660e-05, 7.1211483658335299e-06, -9.8397636034453450e-04, 1.6222237447238830e-04,
        9.6528523782578516e-05, 1.2526412858398512e-03, -1.5903892253583135e-05, -2.0651557490663052e-04, -1.2288467388474000e-04, 2.0246291318839486e-05,
        -5.6925297375807257e-04, 9.3835695361493563e-05, 1.4291213510234853e-04, 7.2468181732210975e-04, -2.3579175340076496e-05, -1.1945659553664517e-04,
        -1.8193286738518963e-04, 3.0017233855795637e-05, -4.0541836830676359e-05, 6.4812209752868387e-06, 1.6233938212623677e-04, 5.1611380785723758e-05,
        -2.6854953961718241e-05, -8.2508536825548366e-06, -2.0666453033267172e-04, 3.4187431139095443e-05, 4.9013450248755377e-04, -8.1515030663001578e-05,
        1.4404674523170864e-04, -6.2396083704236598e-04, -2.3949705609325615e-05, 1.0377189629761350e-04, -1.8337727148740832e-04, 3.0488933717316898e-05,
        8.7363231171246156e-04, -1.4558955682893914e-04, 7.7365818157923587e-05, -1.1121688959189492e-03, -1.3043739320776170e-05, 1.8534133239460972e-04,
        -9.8489782725531659e-05, 1.6605202170929833e-05, 9.2429031121327328e-04, -1.5463556819535182e-04, -4.8118408509521975e-05, -1.1766585566365616e-03,
        7.8210222232853527e-06, 1.9685726688898814e-04, 6.1256659750270513e-05, -9.9564742906303526e-06, 4.2047339144771385e-04, -7.0544580168748378e-05,
        -2.4276042573941192e-04, -5.3527945482356677e-04, 4.0728932350295730e-05, 8.9806073776681691e-05, 3.0904373733406230e-04, -5.1849560869830381e-05,
        -5.6182935935282811e-04, 1.7477050575086486e-04, 1.5794234573042703e-05, 7.1523125908826493e-04, -4.9015805954909655e-06, -2.2248984820528401e-04,
        -2.0106692702968120e-05, 6.2399082612448420e-06, -4.5678597105276730e-04, 1.4215113869939708e-04, 4.4852593938587782e-05, 5.8150682190454017e-04,
        -1.3931219752469142e-05, -1.8096408850885645e-04, -5.7099147102286593e-05, 1.7735000277791322e-05, -2.6422489214980252e-04, 8.2272349665436255e-05,
        6.6322596801483795e-05, 3.3636886209089138e-04, -2.0639815215036942e-05, -1.0473599369591789e-04, -8.4431319984719147e-05, 2.6275310782271061e-05,
        -1.9172276725498251e-05, 5.8565003039016137e-06, 7.5158596679360178e-05, 2.4407075553624133e-05, -2.3478895984198962e-05, -7.4555592663146275e-06,
        -9.5679901449448646e-05, 2.9889574222544112e-05, 2.2595287291712833e-04, -7.1010304284437612e-05, 6.6364474451957597e-05, -2.8764704994641428e-04,
        -2.0900174336625055e-05, 9.0398959214416647e-05, -8.4484631909731762e-05, 2.6606758363727290e-05, 4.0170043536066267e-04, -1.2688043902193714e-04,
        3.5103428034237590e-05, -5.1138073042365968e-04, -1.1356462993730947e-05, 1.6152387668003697e-04, -4.4688068740591611e-05, 1.4457231881138533e-05,
        4.2273028355949746e-04, -1.3494024878531186e-04, -2.2961839515623369e-05, -5.3815256880357897e-04, 6.7031296572293737e-06, 1.7178433706557214e-04,
        2.9231340645244033e-05, -8.5333523155403361e-06, 1.9147960544159756e-04, -6.1665043824634419e-05, -1.1055080174602965e-04, -2.4376120081633150e-04,
        3.5602329651742721e-05, 7.8502068648085319e-05, 1.4073559490929551e-04, -4.5323190465914605e-05, -1.2955898186390926e-04, 7.4800913055063687e-05,
        3.6522239721318641e-06, 1.6493376891421171e-04, -2.1086124934482099e-06, -9.5224555881080988e-05, -4.6494272799621577e-06, 2.6843480916631399e-06,
        -1.0528648823581542e-04, 6.0787182324984421e-05, 1.0361506757951604e-05, 1.3403391313090328e-04, -5.9822187159139861e-06, -7.7384515826667597e-05,
        -1.3190612774451585e-05, 7.6156038361061387e-06, -6.0862501235595707e-05, 3.5138981471927000e-05, 1.5286694095101679e-05, 7.7480399814173790e-05,
        -8.8257769508283041e-06, -4.4733329689633267e-05, -1.9460573362581266e-05, 1.1235567269471202e-05, -4.5142603222996201e-06, 2.6063094122710723e-06,
        1.7245977964551712e-05, 5.7468340527632004e-06, -9.9569700202718593e-06, -3.3179361873508551e-06, -2.1954820139700859e-05, 1.2675621317665824e-05,
        5.1479709017778204e-05, -2.9721823859217862e-05, 1.5082155306937880e-05, -6.5535729817849114e-05, -8.7076864264203733e-06, 3.7837071251873586e-05,
        -1.9200187299523851e-05, 1.1085233305871462e-05, 9.0968534244011416e-05, -5.2520707733565707e-05, 7.7167285674099497e-06, -1.1580658468917291e-04,
        -4.4552553156573706e-06, 6.6860962844225174e-05, -9.8237042928275584e-06, 5.6717183179032468e-06, 9.4503458338676474e-05, -5.4561597111185429e-05,
        -5.6758391897901815e-06, -1.2030668453073529e-04, 3.2769472841023753e-06, 6.9459096699131382e-05, 7.2255704379213280e-06, -4.1716850373824735e-06,
        4.2336308243184655e-05, -2.4442878960697636e-05, -2.4442878960697636e-05, -5.3895814709293923e-05, 1.4112102747728219e-05, 3.1116763130605032e-05,
        3.1116763130605032e-05, -1.7965271569764636e-05,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_DxDx_DyNy(void)
{
  static const double sol[512] = {
        2.4528129331630134e-04, 1.4161322072330483e-04, 1.4161322072328683e-04, -3.1225290264548392e-04, 8.1760431105401652e-05, -1.8027929739760661e-04,
        -1.8027929739761000e-04, -1.0408430088182000e-04, 5.3693018033431348e-04, 3.0999678415214295e-04, 2.6770342705417267e-05, -6.8353360772280454e-04,
        1.5455864567346426e-05, -3.9463831241897844e-04, -3.4079717623695104e-05, -1.9675934143992066e-05, 5.0481002617899188e-04, 2.9145220450423542e-04,
        -4.5314922353404482e-05, -6.4264336602190574e-04, -2.6162582619074331e-05, -3.7103032036578127e-04, 5.7687709676937853e-05, 3.3306014710932019e-05,
        2.8251817425993943e-04, 1.6311194395979506e-04, -8.3025338190720856e-05, -3.5965694232132722e-04, -4.7934701354125037e-05, -2.0764803246505980e-04,
        1.0569457822360065e-04, 6.1022793189366422e-05, -2.2431527513860707e-05, -1.2950848448608950e-05, -9.3037454217678806e-05, 2.8556232243852698e-05,
        -5.3715199237144208e-05, 1.6486948373083493e-05, 1.1844040261456511e-04, 6.8381598332347218e-05, -3.2502275688907155e-04, -1.8765197618256113e-04,
        -8.1663673516519012e-05, 4.1376697705843729e-04, -4.7148543887784260e-05, 2.3888847558638216e-04, 1.0396112459887155e-04, 6.0021983272434631e-05,
        -5.6185672699466782e-04, -3.2438813257639734e-04, -5.5072482877194797e-05, 7.1526609919159482e-04, -3.1796112814151838e-05, 4.1295907491050522e-04,
        7.0109474725142730e-05, 4.0477724105347290e-05, -6.9079330390630332e-04, -3.9882969996476124e-04, -1.9369084511126962e-05, 8.7940752169269742e-04,
        -1.1182746156417725e-05, 5.0772616937670993e-04, 2.4657619741037014e-05, 1.4236083395034633e-05, 6.2367351899500370e-04, 7.6851632631555648e-05,
        3.6007807407826783e-04, -7.9396134934036906e-04, 4.4370310787568746e-05, -9.7835203972601709e-05, -4.5839379876788942e-04, -5.6485181349849416e-05,
        1.3669669930979761e-03, 1.6922519313425672e-04, 6.9062613281975358e-05, -1.7402036887708588e-03, 8.9615892346698773e-06, -2.1543044331873924e-04,
        -8.7919470620924958e-05, -1.1408461742061553e-05, 1.2870709745400632e-03, 1.6018636460856113e-04, -1.1519060110356575e-04, -1.6384928596765102e-03,
        -1.4180159317139120e-05, -2.0392365286805704e-04, 1.4664224517804710e-04, 1.8051910306271542e-05, 7.2082558778293704e-04, 8.9944959225516229e-05,
        -2.1173132536804489e-04, -9.1763982097139551e-04, -2.6373734989019563e-05, -1.1450353272675932e-04, 2.6954245076451052e-04, 3.3574820128275405e-05,
        -5.7211901773804499e-05, -7.1296099928855400e-06, -2.3746882866716045e-04, 7.2833040600929667e-05, -2.9672293680691639e-05, 9.0762788507323819e-06,
        3.0230732248929024e-04, 3.7774017352392656e-05, -8.2967964857642753e-04, -1.0371181608153377e-04, -2.0851563288963774e-04, 1.0562153967439373e-03,
        -2.6089469003517304e-05, 1.3202929245949090e-04, 2.6544874554620365e-04, 3.3212938152299769e-05, -1.4344421799159632e-03, -1.7939931355866409e-04,
        -1.4064417737504262e-04, 1.8261023019736982e-03, -1.7608728039193037e-05, 2.2838250579132398e-04, 1.7904566643376686e-04, 2.2416615502119518e-05,
        -1.7637269696838659e-03, -2.2062884079452653e-04, -4.9468484637505455e-05, 2.2452950174551419e-03, -6.1951506089930910e-06, 2.8086934398448092e-04,
        6.2975360691774624e-05, 7.8866746576120576e-06, 7.9228320374871762e-04, 2.0495214248979126e-05, 4.5742492095879333e-04, -1.0086082258578147e-03,
        1.1832917463663955e-05, -2.6091227965457842e-05, -5.8232023070591639e-04, -1.5063777489280264e-05, 1.7385437019493434e-03, 4.5304719745675459e-05,
        8.8898832400847294e-05, -2.2132357097847822e-03, 2.4908572133151579e-06, -5.7674721348962351e-05, -1.1317177141254004e-04, -3.1709609176348154e-06,
        1.6392548325408021e-03, 4.3147080612406641e-05, -1.4622328788140213e-04, -2.0868370054480641e-03, -3.7365707477220556e-06, -5.4927960382742277e-05,
        1.8614809738656446e-04, 4.7568041008815300e-06, 9.1873938193646955e-04, 2.4320623105365984e-05, -2.6976650148029420e-04, -1.1695920015167992e-03,
        -7.1328895618548160e-06, -3.0961126534042701e-05, 3.4342355254599055e-04, 9.0804538732882195e-06, -7.3130563309116366e-05, -2.0610335296653078e-06,
        -3.0288987840834431e-04, 9.3098133806460183e-05, -8.0985669980487482e-06, 2.6237781666334921e-06, 3.8559093698597817e-04, 1.0309799896355122e-05,
        -1.0586408004009239e-03, -2.8478966558305924e-05, -2.6609472226574274e-04, 1.3476921061319287e-03, -7.1538337474704410e-06, 3.6254864168172258e-05,
        3.3874923065981298e-04, 9.1071166597727465e-06, -1.8304802004777767e-03, -4.9253344222025752e-05, -1.7952696313978798e-04, 2.3302745517463576e-03,
        -4.8402587889251591e-06, 6.2701478332859929e-05, 2.2854500881667826e-04, 6.1618431473638616e-06, -2.2508116665430249e-03, -6.0589640055256245e-05,
        -6.3151521941853120e-05, 2.8653733298781361e-03, -1.7047546620036213e-06, 7.7133036611583822e-05, 8.0394414780744056e-05, 2.1702209096831663e-06,
        7.8392504372925523e-04, -2.5320800186145820e-05, 4.5259933502152905e-04, -9.9796795365614537e-04, -1.4618970803541214e-05, 3.2234391985342402e-05,
        -5.7617706668595004e-04, 1.8610534889883806e-05, 1.7210227955432836e-03, -5.5420419775660136e-05, 8.8434304245293659e-05, -2.1909308947362935e-03,
        -2.7590526689701884e-06, 7.0552412321432068e-05, -1.1258040847992118e-04, 3.5123844659711566e-06, 1.6238014988128491e-03, -5.2069067000117360e-05,
        -1.4456504608357444e-04, -2.0671642931638810e-03, 4.6939570960792625e-06, 6.6286006115705404e-05, 1.8403708921443886e-04, -5.9755952373196969e-06,
        9.1036543165675754e-04, -2.9155325553538169e-05, -2.6733745933858315e-04, -1.1589316276819452e-03, 8.5352976962446396e-06, 3.7115896237257653e-05,
        3.4033128468836791e-04, -1.0865775553293072e-05, -7.3041740691777587e-05, 2.1123152916959775e-06, -3.0043293633985239e-04, 9.2985059059862732e-05,
        9.5170831626758565e-06, -2.6890619020200572e-06, 3.8246315140499926e-04, -1.2115627743500696e-05, -1.0507771170938653e-03, 3.3019066232458140e-05,
        -2.6406284642211314e-04, 1.3376813225743343e-03, 8.3269378127498946e-06, -4.2034592749948440e-05, 3.3616257139440221e-04, -1.0600525082959479e-05,
        -1.8168031392969675e-03, 5.7149799176489184e-05, -1.7820248370310621e-04, 2.3128631055018096e-03, 5.6049473482294999e-06, -7.2753981509130305e-05,
        2.2685889348759484e-04, -7.1353222864949212e-06, -2.2340461539592260e-03, 7.0269213258614450e-05, -6.2692883129639935e-05, 2.8440301613965437e-03,
        1.9695499036957175e-06, -8.9455520679790750e-05, 7.9810549217890943e-05, -2.5073158495669610e-06, 6.4136415990857722e-04, -5.6986764463633543e-05,
        3.7029177037179703e-04, -8.1648224320975094e-04, -3.2901323803324068e-05, 7.2546431794953979e-05, -4.7139624290569486e-04, 4.1884701925557338e-05,
        1.4085793923499818e-03, -1.2496886316451560e-04, 7.2660150604629298e-05, -1.7931779383632382e-03, -6.3481591817058695e-06, 1.5909036411152085e-04,
        -9.2499279607519873e-05, 8.0814606941412736e-06, 1.3297080252287884e-03, -1.1772587914055242e-04, -1.1819655564340408e-04, -1.6927715315548004e-03,
        1.0529897957677164e-05, 1.4986975558194462e-04, 1.5046894560673474e-04, -1.3404981510784272e-05, 7.4559480183162229e-04, -6.5975041930547919e-05,
        -2.1904137112216997e-04, -9.4917202172908996e-04, 1.9348461836307433e-05, 8.3988868724733218e-05, 2.7884843156041278e-04, -2.4631366250678018e-05,
        -6.0575372069946429e-05, 5.0851459879353555e-06, -2.4640119579254818e-04, 7.7114872895287390e-05, 2.1678156787094268e-05, -6.4735943521869927e-06,
        3.1367858331675646e-04, -2.7597161158338441e-05, -8.6268379122424329e-04, 7.5576666091140113e-05, -2.1669631592545041e-04, 1.0982309911733632e-03,
        1.9020141320076898e-05, -9.6212120541947429e-05, 2.7586308244495630e-04, -2.4213401094000017e-05, -1.4913852692680471e-03, 1.3073029568382250e-04,
        -1.4628465166305546e-04, 1.8985931336037425e-03, 1.2822821572043187e-05, -1.6642489828339754e-04, 1.8622621593640876e-04, -1.6323965035577555e-05,
        -1.8339090424506746e-03, 1.6075005578360718e-04, -5.1471540987780780e-05, 2.3346396048014523e-03, 4.5090950025754530e-06, -2.0464125429303859e-04,
        6.5525331588775406e-05, -5.7402583940352802e-06, 4.2047339144771884e-04, -7.0544580168746372e-05, 2.4276042573941360e-04, -5.3527945482357089e-04,
        -4.0728932350297417e-05, 8.9806073776680376e-05, -3.0904373733406609e-04, 5.1849560869832665e-05, 9.2429031121328423e-04, -1.5463556819535337e-04,
        4.8118408509523703e-05, -1.1766585566365740e-03, -7.8210222232856864e-06, 1.9685726688899052e-04, -6.1256659750271475e-05, 9.9564742906302052e-06,
        8.7363231171247609e-04, -1.4558955682894122e-04, -7.7365818157923235e-05, -1.1121688959189638e-03, 1.3043739320776170e-05, 1.8534133239461173e-04,
        9.8489782725531320e-05, -1.6605202170929907e-05, 4.9013450248756884e-04, -8.1515030663003706e-05, -1.4404674523170886e-04, -6.2396083704238073e-04,
        2.3949705609325615e-05, 1.0377189629761538e-04, 1.8337727148740848e-04, -3.0488933717316904e-05, -4.0541836830662847e-05, 6.4812209752848744e-06,
        -1.6233938212623750e-04, 5.1611380785709785e-05, 2.6854953961718350e-05, -8.2508536825531239e-06, 2.0666453033267208e-04, -3.4187431139095545e-05,
        -5.6925297375806162e-04, 9.3835695361491842e-05, -1.4291213510234937e-04, 7.2468181732209653e-04, 2.3579175340076526e-05, -1.1945659553664353e-04,
        1.8193286738518982e-04, -3.0017233855795586e-05, -9.8397636034452582e-04, 1.6222237447238684e-04, -9.6528523782579045e-05, 1.2526412858398384e-03,
        1.5903892253583209e-05, -2.0651557490662865e-04, 1.2288467388474008e-04, -2.0246291318839398e-05, -1.2100112667458827e-03, 1.9945749296636326e-04,
        -3.3972790274496219e-05, 1.5403927676946908e-03, 5.5938134322216204e-06, -2.5391737091355897e-04, 4.3248721623872789e-05, -7.1211483658335443e-06,
        1.9147960544160306e-04, -6.1665043824636221e-05, 1.1055080174602976e-04, -2.4376120081633551e-04, -3.5602329651741962e-05, 7.8502068648086647e-05,
        -1.4073559490929594e-04, 4.5323190465914280e-05, 4.2273028355950445e-04, -1.3494024878531264e-04, 2.2961839515624121e-05, -5.3815256880358536e-04,
        -6.7031296572295651e-06, 1.7178433706557325e-04, -2.9231340645244964e-05, 8.5333523155405428e-06, 4.0170043536067135e-04, -1.2688043902193852e-04,
        -3.5103428034237394e-05, -5.1138073042366814e-04, 1.1356462993730807e-05, 1.6152387668003857e-04, 4.4688068740591320e-05, -1.4457231881138467e-05,
        2.2595287291713638e-04, -7.1010304284439536e-05, -6.6364474451958058e-05, -2.8764704994642274e-04, 2.0900174336624889e-05, 9.0398959214418450e-05,
        8.4484631909732033e-05, -2.6606758363727233e-05, -1.9172276725491939e-05, 5.8565003038994394e-06, -7.5158596679360720e-05, 2.4407075553616452e-05,
        2.3478895984198975e-05, -7.4555592663126946e-06, 9.5679901449448836e-05, -2.9889574222544105e-05, -2.6422489214979829e-04, 8.2272349665434154e-05,
        -6.6322596801484445e-05, 3.3636886209088439e-04, 2.0639815215036979e-05, -1.0473599369591595e-04, 8.4431319984719323e-05, -2.6275310782271078e-05,
        -4.5678597105276454e-04, 1.4215113869939523e-04, -4.4852593938587911e-05, 5.8150682190453355e-04, 1.3931219752469295e-05, -1.8096408850885461e-04,
        5.7099147102286688e-05, -1.7735000277791359e-05, -5.6182935935282605e-04, 1.7477050575086326e-04, -1.5794234573042971e-05, 7.1523125908825875e-04,
        4.9015805954909655e-06, -2.2248984820528222e-04, 2.0106692702968208e-05, -6.2399082612448200e-06, 4.2336308243184337e-05, -2.4442878960699218e-05,
        2.4442878960699808e-05, -5.3895814709293103e-05, -1.4112102747727776e-05, 3.1116763130606509e-05, -3.1116763130606143e-05, 1.7965271569764575e-05,
        9.4503458338682613e-05, -5.4561597111185097e-05, 5.6758391897917417e-06, -1.2030668453073874e-04, -3.2769472841017108e-06, 6.9459096699131965e-05,
        -7.2255704379226824e-06, 4.1716850373820128e-06, 9.0968534244016471e-05, -5.2520707733566358e-05, -7.7167285674121299e-06, -1.1580658468917736e-04,
        4.4552553156561441e-06, 6.6860962844225892e-05, 9.8237042928283342e-06, -5.6717183179027064e-06, 5.1479709017779133e-05, -2.9721823859220088e-05,
        -1.5082155306938090e-05, -6.5535729817851269e-05, 8.7076864264206867e-06, 3.7837071251875416e-05, 1.9200187299524397e-05, -1.1085233305871360e-05,
        -4.5142603223004349e-06, 2.6063094122691225e-06, -1.7245977964552505e-05, 5.7468340527625220e-06, 9.9569700202717051e-06, -3.3179361873487524e-06,
        2.1954820139701161e-05, -1.2675621317665762e-05, -6.0862501235599170e-05, 3.5138981471924683e-05, -1.5286694095102407e-05, 7.7480399814173464e-05,
        8.8257769508282516e-06, -4.4733329689631329e-05, 1.9460573362581161e-05, -1.1235567269471365e-05, -1.0528648823582079e-04, 6.0787182324981514e-05,
        -1.0361506757951987e-05, 1.3403391313090323e-04, 5.9822187159136922e-06, -7.7384515826665673e-05, 1.3190612774451854e-05, -7.6156038361059896e-06,
        -1.2955898186391113e-04, 7.4800913055062928e-05, -3.6522239721294407e-06, 1.6493376891421133e-04, 2.1086124934497494e-06, -9.5224555881079430e-05,
        4.6494272799617046e-06, -2.6843480916634906e-06,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_NxDx_DyDy(void)
{
  static const double sol[512] = {
        -6.9079330390640719e-04, 1.9369084511107887e-05, -3.9882969996474807e-04, 8.7940752169279673e-04, 1.1182746156434002e-05, -2.4657619741025694e-05,
        5.0772616937669475e-04, -1.4236083395045803e-05, -1.7637269696839362e-03, 4.9468484637509975e-05, -2.2062884079452052e-04, 2.2452950174552027e-03,
        6.1951506089905610e-06, -6.2975360691778988e-05, 2.8086934398447404e-04, -7.8866746576099621e-06, -2.2508116665430774e-03, 6.3151521941853120e-05,
        -6.0589640055252084e-05, 2.8653733298781795e-03, 1.7047546620036724e-06, -8.0394414780744964e-05, 7.7133036611580502e-05, -2.1702209096831959e-06,
        -2.2340461539592654e-03, 6.2692883129640016e-05, 7.0269213258618177e-05, 2.8440301613965775e-03, -1.9695499036957175e-06, -7.9810549217891811e-05,
        -8.9455520679792837e-05, 2.5073158495670495e-06, -1.8339090424507008e-03, 5.1471540987781308e-05, 1.6075005578361081e-04, 2.3346396048014796e-03,
        -4.5090950025754014e-06, -6.5525331588775718e-05, -2.0464125429304049e-04, 5.7402583940354860e-06, -1.2100112667458984e-03, 3.3972790274496971e-05,
        1.9945749296636572e-04, 1.5403927676947125e-03, -5.5938134322214671e-06, -4.3248721623873019e-05, -2.5391737091356016e-04, 7.1211483658333825e-06,
        -5.6182935935283765e-04, 1.5794234573048036e-05, 1.7477050575086299e-04, 7.1523125908827686e-04, -4.9015805954886396e-06, -2.0106692702968323e-05,
        -2.2248984820528309e-04, 6.2399082612450707e-06, -1.2955898186394569e-04, 3.6522239721253135e-06, 7.4800913055049985e-05, 1.6493376891422626e-04,
        -2.1086124934573583e-06, -4.6494272799621899e-06, -9.5224555881080460e-05, 2.6843480916630107e-06, -5.6185672699473244e-04, 5.5072482877236356e-05,
        -3.2438813257639447e-04, 7.1526609919165174e-04, 3.1796112814129558e-05, -7.0109474725178509e-05, 4.1295907491049807e-04, -4.0477724105331453e-05,
        -1.4344421799160239e-03, 1.4064417737504346e-04, -1.7939931355866488e-04, 1.8261023019737433e-03, 1.7608728039191811e-05, -1.7904566643377166e-04,
        2.2838250579132406e-04, -2.2416615502117464e-05, -1.8304802004778312e-03, 1.7952696313978692e-04, -4.9253344222021144e-05, 2.3302745517463975e-03,
        4.8402587889253895e-06, -2.2854500881667913e-04, 6.2701478332856839e-05, -6.1618431473636253e-06, -1.8168031392970061e-03, 1.7820248370310613e-04,
        5.7149799176493609e-05, 2.3128631055018408e-03, -5.6049473482291162e-06, -2.2685889348759538e-04, -7.2753981509132270e-05, 7.1353222864948475e-06,
        -1.4913852692680716e-03, 1.4628465166305601e-04, 1.3073029568382662e-04, 1.8985931336037672e-03, -1.2822821572042982e-05, -1.8622621593640965e-04,
        -1.6642489828339919e-04, 1.6323965035577511e-05, -9.8397636034453667e-04, 9.6528523782581091e-05, 1.6222237447239060e-04, 1.2526412858398577e-03,
        -1.5903892253582650e-05, -1.2288467388474103e-04, -2.0651557490663017e-04, 2.0246291318839398e-05, -4.5678597105276139e-04, 4.4852593938591367e-05,
        1.4215113869939943e-04, 5.8150682190455003e-04, -1.3931219752469052e-05, -5.7099147102287596e-05, -1.8096408850885493e-04, 1.7735000277791423e-05,
        -1.0528648823579432e-04, 1.0361506757991369e-05, 6.0787182324990825e-05, 1.3403391313092093e-04, -5.9822187158932228e-06, -1.3190612774449762e-05,
        -7.7384515826664616e-05, 7.6156038361076735e-06, -3.2502275688907800e-04, 8.1663673516511111e-05, -1.8765197618261133e-04, 4.1376697705841322e-04,
        4.7148543887775912e-05, -1.0396112459888250e-04, 2.3888847558642490e-04, -6.0021983272421675e-05, -8.2967964857650755e-04, 2.0851563288962598e-04,
        -1.0371181608152599e-04, 1.0562153967439828e-03, 2.6089469003523400e-05, -2.6544874554619845e-04, 1.3202929245948819e-04, -3.3212938152303442e-05,
        -1.0586408004009799e-03, 2.6609472226574263e-04, -2.8478966558299951e-05, 1.3476921061319647e-03, 7.1538337474710669e-06, -3.3874923065981428e-04,
        3.6254864168169568e-05, -9.1071166597727990e-06, -1.0507771170939024e-03, 2.6406284642211417e-04, 3.3019066232463053e-05, 1.3376813225743620e-03,
        -8.3269378127498438e-06, -3.3616257139440362e-04, -4.2034592749950561e-05, 1.0600525082959492e-05, -8.6268379122426443e-04, 2.1669631592545168e-04,
        7.5576666091144490e-05, 1.0982309911733842e-03, -1.9020141320076871e-05, -2.7586308244495760e-04, -9.6212120541949218e-05, 2.4213401094000003e-05,
        -5.6925297375806660e-04, 1.4291213510235064e-04, 9.3835695361496762e-05, 7.2468181732211192e-04, -2.3579175340076520e-05, -1.8193286738519131e-04,
        -1.1945659553664515e-04, 3.0017233855795525e-05, -2.6422489214978181e-04, 6.6322596801488701e-05, 8.2272349665441635e-05, 3.3636886209089366e-04,
        -2.0639815215035295e-05, -8.4431319984722521e-05, -1.0473599369591775e-04, 2.6275310782270177e-05, -6.0862501235539220e-05, 1.5286694095082349e-05,
        3.5138981471942260e-05, 7.7480399814187993e-05, -8.8257769508439589e-06, -1.9460573362585078e-05, -4.4733329689626524e-05, 1.1235567269471844e-05,
        -2.2431527514113289e-05, 9.3037454217544569e-05, -1.2950848448533081e-05, 2.8556232243984154e-05, 5.3715199237225361e-05, -1.1844040261446435e-04,
        1.6486948373042730e-05, -6.8381598332408421e-05, -5.7211901773897239e-05, 2.3746882866716487e-04, -7.1296099928691076e-06, 7.2833040600977603e-05,
        2.9672293680690542e-05, -3.0230732248929393e-04, 9.0762788507249009e-06, -3.7774017352391789e-05, -7.3130563309168069e-05, 3.0288987840834702e-04,
        -2.0610335296580487e-06, 9.3098133806489877e-05, 8.0985669980488566e-06, -3.8559093698598045e-04, 2.6237781666304576e-06, -1.0309799896355221e-05,
        -7.3041740691809096e-05, 3.0043293633985466e-04, 2.1123152917003706e-06, 9.2985059059884443e-05, -9.5170831626762140e-06, -3.8246315140500132e-04,
        -2.6890619020216670e-06, 1.2115627743500962e-05, -6.0575372069963200e-05, 2.4640119579254937e-04, 5.0851459879394492e-06, 7.7114872895303545e-05,
        -2.1678156787094499e-05, -3.1367858331675787e-04, -6.4735943521885902e-06, 2.7597161158338549e-05, -4.0541836830664609e-05, 1.6233938212623799e-04,
        6.4812209752894441e-06, 5.1611380785719807e-05, -2.6854953961718529e-05, -2.0666453033267354e-04, -8.2508536825550670e-06, 3.4187431139095416e-05,
        -1.9172276725473528e-05, 7.5158596679357562e-05, 5.8565003039065181e-06, 2.4407075553617848e-05, -2.3478895984200893e-05, -9.5679901449450205e-05,
        -7.4555592663157456e-06, 2.9889574222544275e-05, -4.5142603222569322e-06, 1.7245977964563073e-05, 2.6063094122765399e-06, 5.7468340527402128e-06,
        -9.9569700202618558e-06, -2.1954820139718512e-05, -3.3179361873593818e-06, 1.2675621317656365e-05, 2.8251817425973917e-04, 8.3025338190885302e-05,
        1.6311194395986483e-04, -3.5965694232120829e-04, 4.7934701354040375e-05, -1.0569457822370867e-04, -2.0764803246510629e-04, -6.1022793189308512e-05,
        7.2082558778287632e-04, 2.1173132536805907e-04, 8.9944959225527071e-05, -9.1763982097136667e-04, 2.6373734989017391e-05, -2.6954245076451784e-04,
        -1.1450353272676479e-04, -3.3574820128275134e-05, 9.1873938193643518e-04, 2.6976650148030157e-04, 2.4320623105370497e-05, -1.1695920015167823e-03,
        7.1328895618531744e-06, -3.4342355254599554e-04, -3.0961126534044084e-05, -9.0804538732871234e-06, 9.1036543165673597e-04, 2.6733745933858662e-04,
        -2.9155325553535425e-05, -1.1589316276819316e-03, -8.5352976962452274e-06, -3.4033128468837051e-04, 3.7115896237257009e-05, 1.0865775553293366e-05,
        7.4559480183160993e-04, 2.1904137112217147e-04, -6.5975041930545168e-05, -9.4917202172907879e-04, -1.9348461836308022e-05, -2.7884843156041446e-04,
        8.3988868724732270e-05, 2.4631366250678289e-05, 4.9013450248756711e-04, 1.4404674523170834e-04, -8.1515030663000454e-05, -6.2396083704237433e-04,
        -2.3949705609326177e-05, -1.8337727148740908e-04, 1.0377189629761366e-04, 3.0488933717317135e-05, 2.2595287291714611e-04, 6.6364474451956215e-05,
        -7.1010304284436148e-05, -2.8764704994642436e-04, -2.0900174336625099e-05, -8.4484631909732398e-05, 9.0398959214415590e-05, 2.6606758363727161e-05,
        5.1479709017806705e-05, 1.5082155306918322e-05, -2.9721823859213156e-05, -6.5535729817872018e-05, -8.7076864264308121e-06, -1.9200187299506148e-05,
        3.7837071251867237e-05, 1.1085233305882173e-05, 5.0481002617914204e-04, 4.5314922353442369e-05, 2.9145220450415503e-04, -6.4264336602200712e-04,
        2.6162582619072271e-05, -5.7687709676956983e-05, -3.7103032036573172e-04, -3.3306014710934500e-05, 1.2870709745400660e-03, 1.1519060110358824e-04,
        1.6018636460855650e-04, -1.6384928596765208e-03, 1.4180159317132296e-05, -1.4664224517806261e-04, -2.0392365286805419e-04, -1.8051910306266937e-05,
        1.6392548325407947e-03, 1.4622328788140983e-04, 4.3147080612405137e-05, -2.0868370054480650e-03, 3.7365707477203691e-06, -1.8614809738656999e-04,
        -5.4927960382739824e-05, -4.7568041008804086e-06, 1.6238014988128394e-03, 1.4456504608357777e-04, -5.2069067000117178e-05, -2.0671642931638762e-03,
        -4.6939570960801832e-06, -1.8403708921444116e-04, 6.6286006115706569e-05, 5.9755952373204787e-06, 1.3297080252287802e-03, 1.1819655564340479e-04,
        -1.1772587914055171e-04, -1.6927715315547928e-03, -1.0529897957677752e-05, -1.5046894560673534e-04, 1.4986975558194465e-04, 1.3404981510784509e-05,
        8.7363231171247240e-04, 7.7365818157922530e-05, -1.4558955682893938e-04, -1.1121688959189588e-03, -1.3043739320776373e-05, -9.8489782725531469e-05,
        1.8534133239461051e-04, 1.6605202170929968e-05, 4.0170043536067303e-04, 3.5103428034234649e-05, -1.2688043902193719e-04, -5.1138073042366749e-04,
        -1.1356462993731803e-05, -4.4688068740589646e-05, 1.6152387668003727e-04, 1.4457231881139441e-05, 9.0968534244010386e-05, 7.7167285674124671e-06,
        -5.2520707733572158e-05, -1.1580658468916999e-04, -4.4552553156533683e-06, -9.8237042928303434e-06, 6.6860962844231083e-05, 5.6717183178996131e-06,
        5.3693018033442353e-04, -2.6770342705478280e-05, 3.0999678415212110e-04, -6.8353360772288336e-04, -1.5455864567310600e-05, 3.4079717623727183e-05,
        -3.9463831241897117e-04, 1.9675934143970172e-05, 1.3669669930980236e-03, -6.9062613281972092e-05, 1.6922519313424257e-04, -1.7402036887709063e-03,
        -8.9615892346685745e-06, 8.7919470620919008e-05, -2.1543044331872856e-04, 1.1408461742061493e-05, 1.7385437019493595e-03, -8.8898832400841277e-05,
        4.5304719745671474e-05, -2.2132357097848017e-03, -2.4908572133149280e-06, 1.1317177141253521e-04, -5.7674721348956699e-05, 3.1709609176355087e-06,
        1.7210227955432816e-03, -8.8434304245292331e-05, -5.5420419775666512e-05, -2.1909308947362927e-03, 2.7590526689673258e-06, 1.1258040847992115e-04,
        7.0552412321438289e-05, -3.5123844659690170e-06, 1.4085793923499720e-03, -7.2660150604630884e-05, -1.2496886316451367e-04, -1.7931779383632289e-03,
        6.3481591817070960e-06, 9.2499279607521838e-05, 1.5909036411151944e-04, -8.0814606941423358e-06, 9.2429031121327957e-04, -4.8118408509523615e-05,
        -1.5463556819535234e-04, -1.1766585566365692e-03, 7.8210222232854560e-06, 6.1256659750271475e-05, 1.9685726688898922e-04, -9.9564742906303086e-06,
        4.2273028355950131e-04, -2.2961839515624145e-05, -1.3494024878531275e-04, -5.3815256880358276e-04, 6.7031296572297566e-06, 2.9231340645244466e-05,
        1.7178433706557334e-04, -8.5333523155407037e-06, 9.4503458338678696e-05, -5.6758391897908312e-06, -5.4561597111185503e-05, -1.2030668453073674e-04,
        3.2769472841020526e-06, 7.2255704379215940e-06, 6.9459096699131477e-05, -4.1716850373821949e-06, 2.4528129331627364e-04, -1.4161322072330534e-04,
        1.4161322072332695e-04, -3.1225290264551026e-04, -8.1760431105412073e-05, 1.8027929739760823e-04, -1.8027929739762065e-04, 1.0408430088182956e-04,
        6.2367351899509228e-04, -3.6007807407824750e-04, 7.6851632631600670e-05, -7.9396134934046816e-04, -4.4370310787535901e-05, 4.5839379876786568e-04,
        -9.7835203972629722e-05, 5.6485181349827101e-05, 7.9228320374876609e-04, -4.5742492095878081e-04, 2.0495214248910872e-05, -1.0086082258578591e-03,
        -1.1832917463701227e-05, 5.8232023070590674e-04, -2.6091227965398181e-05, 1.5063777489310737e-05, 7.8392504372921544e-04, -4.5259933502155209e-04,
        -2.5320800186128338e-05, -9.9796795365610460e-04, 1.4618970803557920e-05, 5.7617706668597292e-04, 3.2234391985331892e-05, -1.8610534889895576e-05,
        6.4136415990857201e-04, -3.7029177037179280e-04, -5.6986764463631171e-05, -8.1648224320974032e-04, 3.2901323803323126e-05, 4.7139624290569340e-04,
        7.2546431794947067e-05, -4.1884701925559473e-05, 4.2047339144771472e-04, -2.4276042573941338e-04, -7.0544580168748093e-05, -5.3527945482356937e-04,
        4.0728932350296049e-05, 3.0904373733406420e-04, 8.9806073776681975e-05, -5.1849560869830896e-05, 1.9147960544159881e-04, -1.1055080174603042e-04,
        -6.1665043824634527e-05, -2.4376120081633304e-04, 3.5602329651742815e-05, 1.4073559490929629e-04, 7.8502068648085658e-05, -4.5323190465914727e-05,
        4.2336308243185197e-05, -2.4442878960697954e-05, -2.4442878960697954e-05, -5.3895814709294411e-05, 1.4112102747728401e-05, 3.1116763130605317e-05,
        3.1116763130605317e-05, -1.7965271569764799e-05,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_DxNx_DyDy(void)
{
  static const double sol[512] = {
        2.4528129331641009e-04, 1.4161322072339013e-04, 1.4161322072323059e-04, -3.1225290264558366e-04, 8.1760431105386229e-05, -1.8027929739766388e-04,
        -1.8027929739757061e-04, -1.0408430088180949e-04, 6.2367351899508393e-04, 3.6007807407826723e-04, 7.6851632631613504e-05, -7.9396134934045298e-04,
        4.4370310787524178e-05, -4.5839379876787668e-04, -9.7835203972628665e-05, -5.6485181349821416e-05, 7.9228320374878549e-04, 4.5742492095878005e-04,
        2.0495214248914016e-05, -1.0086082258578244e-03, 1.1832917463701158e-05, -5.8232023070591693e-04, -2.6091227965387976e-05, -1.5063777489315933e-05,
        7.8392504372923799e-04, 4.5259933502155545e-04, -2.5320800186129666e-05, -9.9796795365606730e-04, -1.4618970803555413e-05, -5.7617706668598300e-04,
        3.2234391985323307e-05, 1.8610534889900770e-05, 6.4136415990862026e-04, 3.7029177037178126e-04, -5.6986764463615016e-05, -8.1648224320972546e-04,
        -3.2901323803334279e-05, -4.7139624290569432e-04, 7.2546431794942595e-05, 4.1884701925559656e-05, 4.2047339144775083e-04, 2.4276042573940094e-04,
        -7.0544580168771268e-05, -5.3527945482357338e-04, -4.0728932350285410e-05, -3.0904373733406024e-04, 8.9806073776675660e-05, 5.1849560869833525e-05,
        1.9147960544158748e-04, 1.1055080174602805e-04, -6.1665043824638729e-05, -2.4376120081634332e-04, -3.5602329651747633e-05, -1.4073559490929369e-04,
        7.8502068648088301e-05, 4.5323190465911326e-05, 4.2336308243161155e-05, 2.4442878960684070e-05, -2.4442878960701085e-05, -5.3895814709309427e-05,
        -1.4112102747730210e-05, -3.1116763130613984e-05, 3.1116763130599923e-05, 1.7965271569761685e-05, 5.3693018033450050e-04, 2.6770342705359106e-05,
        3.0999678415209730e-04, -6.8353360772295210e-04, 1.5455864567378393e-05, -3.4079717623668880e-05, -3.9463831241894352e-04, -1.9675934144003216e-05,
        1.3669669930980539e-03, 6.9062613281974694e-05, 1.6922519313423946e-04, -1.7402036887709158e-03, 8.9615892346710530e-06, -8.7919470620922207e-05,
        -2.1543044331872175e-04, -1.1408461742063897e-05, 1.7385437019493855e-03, 8.8898832400845966e-05, 4.5304719745672186e-05, -2.2132357097847956e-03,
        2.4908572133136756e-06, -1.1317177141254158e-04, -5.7674721348954423e-05, -3.1709609176349484e-06, 1.7210227955433098e-03, 8.8434304245292155e-05,
        -5.5420419775666248e-05, -2.1909308947362814e-03, -2.7590526689688848e-06, -1.1258040847992607e-04, 7.0552412321439102e-05, 3.5123844659692677e-06,
        1.4085793923499944e-03, 7.2660150604627523e-05, -1.2496886316451722e-04, -1.7931779383632191e-03, -6.3481591817073264e-06, -9.2499279607523736e-05,
        1.5909036411151779e-04, 8.0814606941438113e-06, 9.2429031121328542e-04, 4.8118408509518702e-05, -1.5463556819535819e-04, -1.1766585566365631e-03,
        -7.8210222232860947e-06, -6.1256659750269686e-05, 1.9685726688898863e-04, 9.9564742906309574e-06, 4.2273028355947908e-04, 2.2961839515620249e-05,
        -1.3494024878532319e-04, -5.3815256880358146e-04, -6.7031296572285555e-06, -2.9231340645240418e-05, 1.7178433706557119e-04, 8.5333523155413458e-06,
        9.4503458338600580e-05, 5.6758391897734915e-06, -5.4561597111207316e-05, -1.2030668453075738e-04, -3.2769472841110353e-06, -7.2255704379161823e-06,
        6.9459096699120947e-05, 4.1716850373823499e-06, 5.0481002617904804e-04, -4.5314922353421830e-05, 2.9145220450423634e-04, -6.4264336602196082e-04,
        -2.6162582619079407e-05, 5.7687709676965182e-05, -3.7103032036577476e-04, 3.3306014710926822e-05, 1.2870709745401126e-03, -1.1519060110358146e-04,
        1.6018636460855631e-04, -1.6384928596765438e-03, -1.4180159317133011e-05, 1.4664224517805789e-04, -2.0392365286805118e-04, 1.8051910306267174e-05,
        1.6392548325408359e-03, -1.4622328788140585e-04, 4.3147080612402216e-05, -2.0868370054480776e-03, -3.7365707477212382e-06, 1.8614809738656579e-04,
        -5.4927960382736503e-05, 4.7568041008804382e-06, 1.6238014988128691e-03, -1.4456504608357661e-04, -5.2069067000120898e-05, -2.0671642931638788e-03,
        4.6939570960793650e-06, 1.8403708921443785e-04, 6.6286006115708751e-05, -5.9755952373199925e-06, 1.3297080252287945e-03, -1.1819655564340603e-04,
        -1.1772587914055666e-04, -1.6927715315547887e-03, 1.0529897957677164e-05, 1.5046894560673396e-04, 1.4986975558194660e-04, -1.3404981510783889e-05,
        8.7363231171246579e-04, -7.7365818157924915e-05, -1.4558955682894646e-04, -1.1121688959189475e-03, 1.3043739320776322e-05, 9.8489782725532608e-05,
        1.8534133239461251e-04, -1.6605202170929098e-05, 4.0170043536063898e-04, -3.5103428034237590e-05, -1.2688043902194603e-04, -5.1138073042364775e-04,
        1.1356462993731508e-05, 4.4688068740596321e-05, 1.6152387668004025e-04, -1.4457231881137132e-05, 9.0968534243944562e-05, -7.7167285673880403e-06,
        -5.2520707733581692e-05, -1.1580658468914775e-04, 4.4552553156694357e-06, 9.8237042928496913e-06, 6.6860962844229497e-05, -5.6717183178946029e-06,
        2.8251817425989937e-04, -8.3025338190759114e-05, 1.6311194395982403e-04, -3.5965694232128895e-04, -4.7934701354103719e-05, 1.0569457822362715e-04,
        -2.0764803246508764e-04, 6.1022793189351798e-05, 7.2082558778295428e-04, -2.1173132536804774e-04, 8.9944959225520430e-05, -9.1763982097140733e-04,
        -2.6373734989020383e-05, 2.6954245076451236e-04, -1.1450353272676037e-04, 3.3574820128275771e-05, 9.1873938193649102e-04, -2.6976650148029723e-04,
        2.4320623105364212e-05, -1.1695920015168078e-03, -7.1328895618542281e-06, 3.4342355254599196e-04, -3.0961126534039862e-05, 9.0804538732875994e-06,
        9.1036543165677055e-04, -2.6733745933858483e-04, -2.9155325553541354e-05, -1.1589316276819444e-03, 8.5352976962447920e-06, 3.4033128468836813e-04,
        3.7115896237260282e-05, -1.0865775553293173e-05, 7.4559480183162316e-04, -2.1904137112217098e-04, -6.5975041930551768e-05, -9.4917202172907988e-04,
        1.9348461836307690e-05, 2.7884843156041278e-04, 8.3988868724735793e-05, -2.4631366250678018e-05, 4.9013450248755551e-04, -1.4404674523170886e-04,
        -8.1515030663008111e-05, -6.2396083704236186e-04, 2.3949705609325855e-05, 1.8337727148740878e-04, 1.0377189629761795e-04, -3.0488933717316698e-05,
        2.2595287291710662e-04, -6.6364474451956418e-05, -7.1010304284444564e-05, -2.8764704994639070e-04, 2.0900174336625617e-05, 8.4484631909733713e-05,
        9.0398959214423478e-05, -2.6606758363726633e-05, 5.1479709017752949e-05, -1.5082155306935768e-05, -2.9721823859212993e-05, -6.5535729817781379e-05,
        8.7076864264203513e-06, 1.9200187299526287e-05, 3.7837071251892241e-05, -1.1085233305871834e-05, -2.2431527513872884e-05, -9.3037454217624434e-05,
        -1.2950848448593547e-05, 2.8556232243847948e-05, -5.3715199237173366e-05, 1.1844040261451382e-04, 1.6486948373078488e-05, 6.8381598332375028e-05,
        -5.7211901773788798e-05, -2.3746882866715864e-04, -7.1296099928848429e-06, 7.2833040600918297e-05, -2.9672293680692828e-05, 3.0230732248928867e-04,
        9.0762788507335424e-06, 3.7774017352393571e-05, -7.3130563309101350e-05, -3.0288987840834501e-04, -2.0610335296664254e-06, 9.3098133806454356e-05,
        -8.0985669980490108e-06, 3.8559093698597833e-04, 2.6237781666355622e-06, 1.0309799896355250e-05, -7.3041740691768927e-05, -3.0043293633985325e-04,
        2.1123152916934326e-06, 9.2985059059864345e-05, 9.5170831626760479e-06, 3.8246315140499947e-04, -2.6890619020178464e-06, -1.2115627743500844e-05,
        -6.0575372069947465e-05, -2.4640119579254834e-04, 5.0851459879323011e-06, 7.7114872895296986e-05, 2.1678156787094434e-05, 3.1367858331675619e-04,
        -6.4735943521846032e-06, -2.7597161158338549e-05, -4.0541836830675247e-05, -1.6233938212623694e-04, 6.4812209752813668e-06, 5.1611380785727993e-05,
        2.6854953961718624e-05, 2.0666453033267132e-04, -8.2508536825505439e-06, -3.4187431139095714e-05, -1.9172276725516337e-05, -7.5158596679359270e-05,
        5.8565003038960250e-06, 2.4407075553646800e-05, 2.3478895984199182e-05, 9.5679901449446152e-05, -7.4555592663082866e-06, -2.9889574222545057e-05,
        -4.5142603223459884e-06, -1.7245977964566013e-05, 2.6063094122603197e-06, 5.7468340528147339e-06, 9.9569700202628604e-06, 2.1954820139689062e-05,
        -3.3179361873405256e-06, -1.2675621317670252e-05, -3.2502275688898297e-04, -8.1663673516515218e-05, -1.8765197618259547e-04, 4.1376697705835895e-04,
        -4.7148543887783813e-05, 1.0396112459888037e-04, 2.3888847558641401e-04, 6.0021983272428126e-05, -8.2967964857640498e-04, -2.0851563288963555e-04,
        -1.0371181608153745e-04, 1.0562153967439215e-03, -2.6089469003518687e-05, 2.6544874554620268e-04, 1.3202929245949519e-04, 3.3212938152300663e-05,
        -1.0586408004009109e-03, -2.6609472226574323e-04, -2.8478966558307740e-05, 1.3476921061319242e-03, -7.1538337474705816e-06, 3.3874923065981352e-04,
        3.6254864168174548e-05, 9.1071166597727244e-06, -1.0507771170938588e-03, -2.6406284642211346e-04, 3.3019066232456236e-05, 1.3376813225743371e-03,
        8.3269378127500487e-06, 3.3616257139440253e-04, -4.2034592749946726e-05, -1.0600525082959604e-05, -8.6268379122424405e-04, -2.1669631592545016e-04,
        7.5576666091137728e-05, 1.0982309911733714e-03, 1.9020141320077115e-05, 2.7586308244495586e-04, -9.6212120541945694e-05, -2.4213401094000305e-05,
        -5.6925297375807105e-04, -1.4291213510234823e-04, 9.3835695361489321e-05, 7.2468181732211106e-04, 2.3579175340076791e-05, 1.8193286738518824e-04,
        -1.1945659553664183e-04, -3.0017233855795931e-05, -2.6422489214981732e-04, -6.6322596801482765e-05, 8.2272349665431173e-05, 3.3636886209090379e-04,
        2.0639815215037020e-05, 8.4431319984715745e-05, -1.0473599369591468e-04, -2.6275310782271935e-05, -6.0862501235637246e-05, -1.5286694095084582e-05,
        3.5138981471916619e-05, 7.7480399814190974e-05, 8.8257769508375181e-06, 1.9460573362573226e-05, -4.4733329689633741e-05, -1.1235567269473010e-05,
        -5.6185672699463572e-04, -5.5072482877231199e-05, -3.2438813257640027e-04, 7.1526609919158365e-04, -3.1796112814134138e-05, 7.0109474725172762e-05,
        4.1295907491050273e-04, 4.0477724105333954e-05, -1.4344421799159428e-03, -1.4064417737504625e-04, -1.7939931355866808e-04, 1.8261023019736887e-03,
        -1.7608728039191811e-05, 1.7904566643377122e-04, 2.2838250579132729e-04, 2.2416615502118027e-05, -1.8304802004777663e-03, -1.7952696313978887e-04,
        -4.9253344222027344e-05, 2.3302745517463559e-03, -4.8402587889248271e-06, 2.2854500881667938e-04, 6.2701478332861487e-05, 6.1618431473634482e-06,
        -1.8168031392969621e-03, -1.7820248370310657e-04, 5.7149799176487768e-05, 2.3128631055018131e-03, 5.6049473482295761e-06, 2.2685889348759498e-04,
        -7.2753981509129112e-05, -7.1353222864950982e-06, -1.4913852692680471e-03, -1.4628465166305517e-04, 1.3073029568382095e-04, 1.8985931336037497e-03,
        1.2822821572043416e-05, 1.8622621593640838e-04, -1.6642489828339656e-04, -1.6323965035577717e-05, -9.8397636034453146e-04, -9.6528523782577987e-05,
        1.6222237447238529e-04, 1.2526412858398478e-03, 1.5903892253583568e-05, 1.2288467388473875e-04, -2.0651557490662816e-04, -2.0246291318839724e-05,
        -4.5678597105277500e-04, -4.4852593938584658e-05, 1.4215113869939401e-04, 5.8150682190454331e-04, 1.3931219752470214e-05, 5.7099147102284574e-05,
        -1.8096408850885512e-04, -1.7735000277791535e-05, -1.0528648823582545e-04, -1.0361506757950515e-05, 6.0787182324986102e-05, 1.3403391313090627e-04,
        5.9822187159117466e-06, 1.3190612774451422e-05, -7.7384515826668979e-05, -7.6156038361048664e-06, -6.9079330390632197e-04, -1.9369084511119769e-05,
        -3.9882969996474915e-04, 8.7940752169272311e-04, -1.1182746156426749e-05, 2.4657619741028168e-05, 5.0772616937670007e-04, 1.4236083395043767e-05,
        -1.7637269696838601e-03, -4.9468484637510063e-05, -2.2062884079452458e-04, 2.2452950174551459e-03, -6.1951506089909447e-06, 6.2975360691778500e-05,
        2.8086934398447859e-04, 7.8866746576102568e-06, -2.2508116665430171e-03, -6.3151521941853567e-05, -6.0589640055257214e-05, 2.8653733298781379e-03,
        -1.7047546620036213e-06, 8.0394414780744815e-05, 7.7133036611584486e-05, 2.1702209096830778e-06, -2.2340461539592216e-03, -6.2692883129640193e-05,
        7.0269213258613216e-05, 2.8440301613965480e-03, 1.9695499036957683e-06, 7.9810549217891404e-05, -8.9455520679789923e-05, -2.5073158495670495e-06,
        -1.8339090424506739e-03, -5.1471540987780868e-05, 1.6075005578360620e-04, 2.3346396048014584e-03, 4.5090950025756064e-06, 6.5525331588775406e-05,
        -2.0464125429303824e-04, -5.7402583940354860e-06, -1.2100112667458858e-03, -3.3972790274495860e-05, 1.9945749296636228e-04, 1.5403927676946978e-03,
        5.5938134322217229e-06, 4.3248721623872484e-05, -2.5391737091355897e-04, -7.1211483658334858e-06, -5.6182935935283006e-04, -1.5794234573042483e-05,
        1.7477050575086367e-04, 7.1523125908826428e-04, 4.9015805954910036e-06, 2.0106692702967991e-05, -2.2248984820528303e-04, -6.2399082612448268e-06,
        -1.2955898186391279e-04, -3.6522239721291976e-06, 7.4800913055063917e-05, 1.6493376891421342e-04, 2.1086124934495961e-06, 4.6494272799615834e-06,
        -9.5224555881080636e-05, -2.6843480916634093e-06,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p1_NxDx_PyPy(void)
{
  static const double sol[512] = {
        2.3004015067766781e-02, -6.5620995429994662e-04, 1.3281374291817376e-02, -2.9285031811004718e-02, -3.7886299375998559e-04, 8.3538153360449507e-04,
        -1.6907720999310375e-02, 4.8230775330259821e-04, 2.3004015067766743e-02, -6.5620995429994662e-04, -1.3281374291817399e-02, -2.9285031811004721e-02,
        3.7886299375998613e-04, 8.3538153360449409e-04, 1.6907720999310372e-02, -4.8230775330259897e-04, -2.3004015067766920e-02, 6.5620995429994662e-04,
        -1.3281374291817378e-02, 2.9285031811004874e-02, 3.7886299375998532e-04, -8.3538153360449409e-04, 1.6907720999310375e-02, -4.8230775330259572e-04,
        -2.3004015067766864e-02, 6.5620995429994662e-04, 1.3281374291817409e-02, 2.9285031811004846e-02, -3.7886299375998613e-04, -8.3538153360449333e-04,
        -1.6907720999310392e-02, 4.8230775330259610e-04, 2.3004015067766864e-02, -6.5620995429994803e-04, 1.3281374291817413e-02, -2.9285031811004836e-02,
        -3.7886299375998738e-04, 8.3538153360449452e-04, -1.6907720999310399e-02, 4.8230775330259767e-04, 2.3004015067766868e-02, -6.5620995429994738e-04,
        -1.3281374291817413e-02, -2.9285031811004839e-02, 3.7886299375998700e-04, 8.3538153360449496e-04, 1.6907720999310406e-02, -4.8230775330259751e-04,
        -2.3004015067766875e-02, 6.5620995429994738e-04, -1.3281374291817416e-02, 2.9285031811004916e-02, 3.7886299375998608e-04, -8.3538153360449333e-04,
        1.6907720999310434e-02, -4.8230775330259653e-04, -2.3004015067766889e-02, 6.5620995429994803e-04, 1.3281374291817407e-02, 2.9285031811004950e-02,
        -3.7886299375998700e-04, -8.3538153360449496e-04, -1.6907720999310406e-02, 4.8230775330259659e-04, 1.8634512108127495e-02, -1.8665237556726594e-03,
        1.0758640581844770e-02, -2.3722479674155590e-02, -1.0776379927864414e-03, 2.3761594399866787e-03, -1.3696180025719200e-02, 1.3718762923137797e-03,
        1.8634512108127461e-02, -1.8665237556726600e-03, -1.0758640581844791e-02, -2.3722479674155590e-02, 1.0776379927864405e-03, 2.3761594399866791e-03,
        1.3696180025719193e-02, -1.3718762923137797e-03, -1.8634512108127634e-02, 1.8665237556726600e-03, -1.0758640581844775e-02, 2.3722479674155736e-02,
        1.0776379927864412e-03, -2.3761594399866822e-03, 1.3696180025719200e-02, -1.3718762923137810e-03, -1.8634512108127582e-02, 1.8665237556726615e-03,
        1.0758640581844805e-02, 2.3722479674155708e-02, -1.0776379927864405e-03, -2.3761594399866826e-03, -1.3696180025719215e-02, 1.3718762923137808e-03,
        1.8634512108127579e-02, -1.8665237556726637e-03, 1.0758640581844798e-02, -2.3722479674155701e-02, -1.0776379927864448e-03, 2.3761594399866843e-03,
        -1.3696180025719224e-02, 1.3718762923137829e-03, 1.8634512108127579e-02, -1.8665237556726650e-03, -1.0758640581844798e-02, -2.3722479674155698e-02,
        1.0776379927864440e-03, 2.3761594399866848e-03, 1.3696180025719226e-02, -1.3718762923137823e-03, -1.8634512108127582e-02, 1.8665237556726615e-03,
        -1.0758640581844805e-02, 2.3722479674155774e-02, 1.0776379927864433e-03, -2.3761594399866839e-03, 1.3696180025719250e-02, -1.3718762923137827e-03,
        -1.8634512108127600e-02, 1.8665237556726628e-03, 1.0758640581844798e-02, 2.3722479674155805e-02, -1.0776379927864422e-03, -2.3761594399866848e-03,
        -1.3696180025719233e-02, 1.3718762923137823e-03, 1.0604640693740158e-02, -2.7695246669819525e-03, 6.1225921591901529e-03, -1.3500132005025856e-02,
        -1.5989858120093413e-03, 3.5257157385353115e-03, -7.7943048471972058e-03, 2.0355729307294634e-03, 1.0604640693740116e-02, -2.7695246669819560e-03,
        -6.1225921591901780e-03, -1.3500132005025865e-02, 1.5989858120093389e-03, 3.5257157385353119e-03, 7.7943048471972023e-03, -2.0355729307294629e-03,
        -1.0604640693740295e-02, 2.7695246669819538e-03, -6.1225921591901581e-03, 1.3500132005025995e-02, 1.5989858120093409e-03, -3.5257157385353158e-03,
        7.7943048471972006e-03, -2.0355729307294647e-03, -1.0604640693740237e-02, 2.7695246669819542e-03, 6.1225921591901902e-03, 1.3500132005025971e-02,
        -1.5989858120093407e-03, -3.5257157385353158e-03, -7.7943048471972188e-03, 2.0355729307294647e-03, 1.0604640693740227e-02, -2.7695246669819551e-03,
        6.1225921591901780e-03, -1.3500132005025953e-02, -1.5989858120093394e-03, 3.5257157385353167e-03, -7.7943048471972205e-03, 2.0355729307294655e-03,
        1.0604640693740225e-02, -2.7695246669819551e-03, -6.1225921591901798e-03, -1.3500132005025953e-02, 1.5989858120093394e-03, 3.5257157385353167e-03,
        7.7943048471972223e-03, -2.0355729307294655e-03, -1.0604640693740237e-02, 2.7695246669819555e-03, -6.1225921591901867e-03, 1.3500132005026026e-02,
        1.5989858120093404e-03, -3.5257157385353175e-03, 7.7943048471972483e-03, -2.0355729307294651e-03, -1.0604640693740246e-02, 2.7695246669819573e-03,
        6.1225921591901815e-03, 1.3500132005026061e-02, -1.5989858120093394e-03, -3.5257157385353184e-03, -7.7943048471972318e-03, 2.0355729307294647e-03,
        3.3969299063247823e-04, -3.1569456526246991e-03, 1.9612183958349583e-04, -4.3244277172230813e-04, -1.8226634223598899e-03, 4.0189181579988681e-03,
        -2.4967095066302401e-04, 2.3203234803717216e-03, 3.3969299063242467e-04, -3.1569456526247021e-03, -1.9612183958352686e-04, -4.3244277172231269e-04,
        1.8226634223598884e-03, 4.0189181579988673e-03, 2.4967095066302146e-04, -2.3203234803717211e-03, -3.3969299063260508e-04, 3.1569456526247017e-03,
        -1.9612183958349588e-04, 4.3244277172242685e-04, 1.8226634223598920e-03, -4.0189181579988733e-03, 2.4967095066301458e-04, -2.3203234803717242e-03,
        -3.3969299063254366e-04, 3.1569456526247043e-03, 1.9612183958353144e-04, 4.3244277172240067e-04, -1.8226634223598910e-03, -4.0189181579988725e-03,
        -2.4967095066302965e-04, 2.3203234803717246e-03, 3.3969299063253238e-04, -3.1569456526247051e-03, 1.9612183958351835e-04, -4.3244277172238121e-04,
        -1.8226634223598925e-03, 4.0189181579988742e-03, -2.4967095066303089e-04, 2.3203234803717246e-03, 3.3969299063252832e-04, -3.1569456526247056e-03,
        -1.9612183958352063e-04, -4.3244277172238012e-04, 1.8226634223598920e-03, 4.0189181579988733e-03, 2.4967095066303165e-04, -2.3203234803717246e-03,
        -3.3969299063254079e-04, 3.1569456526247051e-03, -1.9612183958352521e-04, 4.3244277172244989e-04, 1.8226634223598923e-03, -4.0189181579988759e-03,
        2.4967095066305648e-04, -2.3203234803717259e-03, -3.3969299063254426e-04, 3.1569456526247064e-03, 1.9612183958352312e-04, 4.3244277172247927e-04,
        -1.8226634223598914e-03, -4.0189181579988768e-03, -2.4967095066303968e-04, 2.3203234803717255e-03, -1.0004294319619630e-02, -2.8151582054434295e-03,
        -5.7759820184846399e-03, 1.2735867044672292e-02, -1.6253323477241504e-03, 3.5838090592689460e-03, 7.3530562666047940e-03, 2.0691131250931472e-03,
        -1.0004294319619697e-02, -2.8151582054434343e-03, 5.7759820184846017e-03, 1.2735867044672292e-02, 1.6253323477241478e-03, 3.5838090592689481e-03,
        -7.3530562666047940e-03, -2.0691131250931464e-03, 1.0004294319619515e-02, 2.8151582054434343e-03, 5.7759820184846485e-03, -1.2735867044672193e-02,
        1.6253323477241535e-03, -3.5838090592689512e-03, -7.3530562666048105e-03, -2.0691131250931490e-03, 1.0004294319619585e-02, 2.8151582054434378e-03,
        -5.7759820184846095e-03, -1.2735867044672217e-02, -1.6253323477241519e-03, -3.5838090592689503e-03, 7.3530562666047966e-03, 2.0691131250931494e-03,
        -1.0004294319619599e-02, -2.8151582054434373e-03, -5.7759820184846242e-03, 1.2735867044672245e-02, -1.6253323477241524e-03, 3.5838090592689529e-03,
        7.3530562666047940e-03, 2.0691131250931481e-03, -1.0004294319619604e-02, -2.8151582054434373e-03, 5.7759820184846208e-03, 1.2735867044672245e-02,
        1.6253323477241519e-03, 3.5838090592689520e-03, -7.3530562666047948e-03, -2.0691131250931490e-03, 1.0004294319619588e-02, 2.8151582054434360e-03,
        5.7759820184846173e-03, -1.2735867044672182e-02, 1.6253323477241522e-03, -3.5838090592689546e-03, -7.3530562666047732e-03, -2.0691131250931498e-03,
        1.0004294319619592e-02, 2.8151582054434386e-03, -5.7759820184846147e-03, -1.2735867044672156e-02, -1.6253323477241513e-03, -3.5838090592689559e-03,
        7.3530562666047896e-03, 2.0691131250931490e-03, -1.7517266800300840e-02, -1.5224584786921185e-03, -1.0113598702620196e-02, 2.2300181684695344e-02,
        -8.7899181250292451e-04, 1.9381505727626407e-03, 1.2875015898636386e-02, 1.1189917549145393e-03, -1.7517266800300923e-02, -1.5224584786921235e-03,
        1.0113598702620148e-02, 2.2300181684695351e-02, 8.7899181250292147e-04, 1.9381505727626422e-03, -1.2875015898636382e-02, -1.1189917549145384e-03,
        1.7517266800300735e-02, 1.5224584786921213e-03, 1.0113598702620219e-02, -2.2300181684695265e-02, 8.7899181250292884e-04, -1.9381505727626449e-03,
        -1.2875015898636408e-02, -1.1189917549145410e-03, 1.7517266800300819e-02, 1.5224584786921235e-03, -1.0113598702620170e-02, -2.2300181684695282e-02,
        -8.7899181250292635e-04, -1.9381505727626449e-03, 1.2875015898636398e-02, 1.1189917549145421e-03, -1.7517266800300829e-02, -1.5224584786921235e-03,
        -1.0113598702620184e-02, 2.2300181684695317e-02, -8.7899181250292581e-04, 1.9381505727626440e-03, 1.2875015898636386e-02, 1.1189917549145388e-03,
        -1.7517266800300836e-02, -1.5224584786921235e-03, 1.0113598702620180e-02, 2.2300181684695306e-02, 8.7899181250292548e-04, 1.9381505727626431e-03,
        -1.2875015898636391e-02, -1.1189917549145393e-03, 1.7517266800300819e-02, 1.5224584786921243e-03, 1.0113598702620179e-02, -2.2300181684695261e-02,
        8.7899181250292613e-04, -1.9381505727626485e-03, -1.2875015898636375e-02, -1.1189917549145419e-03, 1.7517266800300829e-02, 1.5224584786921256e-03,
        -1.0113598702620170e-02, -2.2300181684695240e-02, -8.7899181250292472e-04, -1.9381505727626490e-03, 1.2875015898636389e-02, 1.1189917549145416e-03,
        -1.8501401394323238e-02, 9.5426810591446164e-04, -1.0681789075397860e-02, 2.3553024408339691e-02, 5.5094694782878678e-04, -1.2148214890143143e-03,
        1.3598344982384715e-02, -7.0137751369975460e-04, -1.8501401394323349e-02, 9.5426810591445394e-04, 1.0681789075397798e-02, 2.3553024408339705e-02,
        -5.5094694782879036e-04, -1.2148214890143111e-03, -1.3598344982384707e-02, 7.0137751369975688e-04, 1.8501401394323141e-02, -9.5426810591446240e-04,
        1.0681789075397893e-02, -2.3553024408339618e-02, -5.5094694782878277e-04, 1.2148214890143139e-03, -1.3598344982384747e-02, 7.0137751369975341e-04,
        1.8501401394323241e-02, -9.5426810591445600e-04, -1.0681789075397836e-02, -2.3553024408339642e-02, 5.5094694782878548e-04, 1.2148214890143119e-03,
        1.3598344982384733e-02, -7.0137751369975449e-04, -1.8501401394323245e-02, 9.5426810591445676e-04, -1.0681789075397853e-02, 2.3553024408339673e-02,
        5.5094694782878439e-04, -1.2148214890143098e-03, 1.3598344982384717e-02, -7.0137751369975406e-04, -1.8501401394323259e-02, 9.5426810591445600e-04,
        1.0681789075397847e-02, 2.3553024408339663e-02, -5.5094694782878548e-04, -1.2148214890143111e-03, -1.3598344982384724e-02, 7.0137751369975406e-04,
        1.8501401394323241e-02, -9.5426810591445318e-04, 1.0681789075397848e-02, -2.3553024408339635e-02, -5.5094694782878374e-04, 1.2148214890143083e-03,
        -1.3598344982384712e-02, 7.0137751369975319e-04, 1.8501401394323262e-02, -9.5426810591445112e-04, -1.0681789075397840e-02, -2.3553024408339611e-02,
        5.5094694782878461e-04, 1.2148214890143072e-03, 1.3598344982384726e-02, -7.0137751369975406e-04, -8.4242802754184182e-03, 4.8637604847417016e-03,
        -4.8637604847417077e-03, 1.0724445933620218e-02, 2.8080934251394820e-03, -6.1917617466852098e-03, 6.1917617466851959e-03, -3.5748153112067321e-03,
        -8.4242802754185605e-03, 4.8637604847416903e-03, 4.8637604847416253e-03, 1.0724445933620250e-02, -2.8080934251394894e-03, -6.1917617466852046e-03,
        -6.1917617466851795e-03, 3.5748153112067356e-03, 8.4242802754183228e-03, -4.8637604847417016e-03, 4.8637604847417597e-03, -1.0724445933620146e-02,
        -2.8080934251394755e-03, 6.1917617466852089e-03, -6.1917617466852376e-03, 3.5748153112067278e-03, 8.4242802754184443e-03, -4.8637604847416921e-03,
        -4.8637604847416895e-03, -1.0724445933620177e-02, 2.8080934251394803e-03, 6.1917617466852072e-03, 6.1917617466852185e-03, -3.5748153112067295e-03,
        -8.4242802754184443e-03, 4.8637604847416973e-03, -4.8637604847417016e-03, 1.0724445933620217e-02, 2.8080934251394829e-03, -6.1917617466852055e-03,
        6.1917617466851968e-03, -3.5748153112067325e-03, -8.4242802754184512e-03, 4.8637604847416973e-03, 4.8637604847416964e-03, 1.0724445933620199e-02,
        -2.8080934251394829e-03, -6.1917617466852089e-03, -6.1917617466852072e-03, 3.5748153112067312e-03, 8.4242802754184477e-03, -4.8637604847416964e-03,
        4.8637604847416973e-03, -1.0724445933620185e-02, -2.8080934251394833e-03, 6.1917617466852011e-03, -6.1917617466852029e-03, 3.5748153112067291e-03,
        8.4242802754184738e-03, -4.8637604847416938e-03, -4.8637604847416825e-03, -1.0724445933620168e-02, 2.8080934251394850e-03, 6.1917617466852011e-03,
        6.1917617466852133e-03, -3.5748153112067286e-03,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_PxPx_PyPy(void)
{
  static const double sol[640] = {
        5.4502158943918211e+00, 1.4782758980221811e+00, -1.0143961216287423e+00, -7.0788264442624271e+00, 4.1786262474674174e-01, -1.8819043793655346e+00,
        1.2913668593654497e+00, 1.4300160023058239e-01, -4.3919645391728357e-01, 1.5665697539380616e+01, -1.4624872597200635e+00, 1.6993190682438275e+00,
        -2.2418353145413441e+00, 1.7213213200231410e+00, -3.4316478965842806e-01, -2.1633011812923613e+00, 2.8539460745428178e+00, 5.2659329213503925e-01,
        1.2999235352023769e-01, 1.5665697539380627e+01, -5.1231786808679800e+00, -3.9077797649552437e-01, 4.6941444639961460e-01, 6.3815280014559725e+00,
        -6.5917730762051951e-01, 4.9747600316724139e-01, -5.9758337641761972e-01, 5.4963863363049092e-01, 3.9419202503096262e-01, 1.5665697539380638e+01,
        7.1054273576010019e-15, -1.7069809072513831e+00, 1.8181896694868822e+00, -1.4048349091629908e-01, 6.5229796081267557e-16, 2.1730550089788037e+00,
        -2.3146282139230303e+00, 1.2499260451994659e-01, -1.2499260451994597e-01, 1.5665697539380631e+01, 2.7278721498501426e+00, -7.0755542166833685e-01,
        -5.2355818804680432e-01, -3.6131739081926675e+00, 4.1518369083043022e-01, 9.0074636843026545e-01, 6.6651052638833630e-01, -3.2079725539104875e-01,
        -3.4211688329996676e-01, 1.5665697539380625e+01, -1.9249877224862288e+00, 1.6136475856827914e-01, -1.4711148381875034e+00, 2.3101029185750042e+00,
        -1.8942923385713109e-03, -2.0542379553857743e-01, 1.8727880636075556e+00, -4.5294792043656273e-01, 1.9363822721812923e-01, 1.5665697539380639e+01,
        -3.0549093633740085e+00, -3.7994249985831274e-01, 1.0685398632759302e+00, 3.7485383873339733e+00, -1.7386900795665511e-01, 4.8368200776802284e-01,
        -1.3602940093361655e+00, -3.7184297847002318e-01, 3.8712131218628787e-01, 1.5665697539380631e+01, 3.3874749822062711e+00, -1.5370291956071291e-01,
        1.8947604832419591e+00, -4.4528747113469986e+00, 3.4505908199700336e-01, 1.9566996785213275e-01, -2.4121059242273448e+00, -1.9863797621842202e-01,
        -1.9863797621842175e-01, 1.5665697539380615e+01, 6.2531003217557561e+00, -2.2418353145413406e+00, 1.6993190682438222e+00, -8.1009304520475087e+00,
        8.3115202323859816e-01, 2.8539460745428187e+00, -2.1633011812923626e+00, -8.0750981210547657e-01, -4.1090887349067584e-01, 1.5665697539380639e+01,
        6.1153471320564314e+00, 9.3575965172465736e-01, -1.5569123679262598e+00, -7.9255651286233890e+00, 7.6024277282743713e-01, -1.1912594861150771e+00,
        1.9820117526159073e+00, -8.2119407263786293e-01, -2.3899601848999669e-01, 1.5665697539380625e+01, 1.4210854715202004e-14, 2.2494971535489068e+00,
        -1.2756734231893674e+00, -1.4048349091629930e-01, -9.5208296129717008e-16, -2.8636999022292637e+00, 1.6239833206725696e+00, -3.0175884102839090e-01,
        3.0175884102839184e-01, 1.5665697539380625e+01, -3.3256561883227675e-01, 1.7005280561147451e+00, 8.4033563321960258e-01, 2.8288585126417820e-01,
        -1.7119007404034631e-01, -2.1648402712362214e+00, -1.0697808916513580e+00, 2.6872211366005361e-01, 1.1327550506052576e-01, 1.5665697539380623e+01,
        1.9249877224862217e+00, 1.4711148381875025e+00, -1.6136475856828625e-01, -2.5910699004075894e+00, 1.8942923385745344e-03, -1.8727880636075520e+00,
        2.0542379553857337e-01, 6.2971415694500832e-01, -3.7040446372657576e-01, 1.5665697539380616e+01, -2.7278721498501426e+00, 9.2492417508745106e-01,
        -2.1560377848025780e+00, 3.3322069263600800e+00, -4.1518369083043161e-01, -1.1774654907158588e+00, 2.7447223855344651e+00, 7.4754870093938741e-01,
        -8.4634562248372333e-02, 1.5665697539380632e+01, -8.1780880442419672e+00, -1.4787766771950537e+00, -2.6228088648617565e-01, 1.0270549879706238e+01,
        -8.3304631557717246e-01, 1.8825418912939917e+00, 3.3389406508122027e-01, 4.7955449618885981e-01, 4.7955449618885876e-01, 1.5665697539380645e+01,
        -3.0549093633739872e+00, -3.5612118829268411e+00, 2.8726145195092379e+00, 3.7485383873339710e+00, -1.7386900795665752e-01, 4.5335652480671502e+00,
        -3.6569532464990155e+00, -1.9507674196157704e-01, 2.1035507567784320e-01, 1.5665697539380650e+01, -5.1231786808679587e+00, -2.0232575732512998e+00,
        2.1018940431553883e+00, 6.3815280014559841e+00, -6.5917730762052029e-01, 2.5756878623133690e+00, -2.6757952355637431e+00, 9.7639007917882614e-01,
        -3.2559420517375662e-02, 1.5665697539380657e+01, 7.1054273576010019e-15, -3.0167309868706003e+00, 5.0843958986766680e-01, -1.4048349091628809e-01,
        -4.6926773889305484e-15, 3.8404192770477850e+00, -6.4726394585405389e-01, 3.0175884102839057e-01, -3.0175884102839112e-01, 1.5665697539380634e+01,
        -6.5960283235613559e-01, -2.7880259776413983e+00, -2.9535395799047315e-01, 6.9921731223805828e-01, 7.0124608833432911e-02, 3.5492686474345190e+00,
        3.7599740870352377e-01, -4.2391812020101916e-01, 1.5827993394684697e-01, 1.5665697539380623e+01, 1.4624872597200849e+00, -2.2418353145413441e+00,
        1.6993190682438222e+00, -2.0022883018557298e+00, 3.4316478965842812e-01, 2.8539460745428231e+00, -2.1633011812923608e+00, -9.5334473768337713e-01,
        2.9675909202810047e-01, 1.5665697539380618e+01, 8.5106536630742653e+00, 1.6257374643920475e-01, 1.6110561095734441e+00, -1.0974886203719290e+01,
        1.0042363896175210e+00, -2.0696288548243347e-01, -2.0509389025866294e+00, -1.0500354508773038e+00, -2.9107116022099289e-01, 1.5665697539380604e+01,
        8.1780880442419708e+00, 3.7874514632244436e+00, -2.0463938995432054e+00, -1.0551516861538820e+01, 8.3304631557717457e-01, -4.8215772879830547e+00,
        2.6051413316078422e+00, -4.7955449618885798e-01, -4.7955449618885859e-01, 1.5665697539380607e+01, -2.7278721498501213e+00, 4.6487098044534978e+00,
        -3.4175961947383620e+00, 3.3322069263600782e+00, -4.1518369083043488e-01, -5.9179936242654509e+00, 4.3507367294468446e+00, 4.9756349189949461e-01,
        1.6535064679152256e-01, 1.5665697539380638e+01, -9.6405753039620095e+00, 1.4711148381874999e+00, -1.6136475856828092e-01, 1.2132354690645652e+01,
        -1.1762111052355997e+00, -1.8727880636075553e+00, 2.0542379553857698e-01, 1.1311403928438426e+00, 4.8455424518915058e-01, 1.5665697539380664e+01,
        2.1316282072803006e-14, 3.5592472331681240e+00, 3.4076656429851582e-02, -1.4048349091629220e-01, 3.8695282355514626e-15, -4.5310641702982455e+00,
        -4.3380947396400749e-02, -1.2499260451994629e-01, 1.2499260451994701e-01, 1.5665697539380647e+01, -3.3256561883228386e-01, 3.3330076528705277e+00,
        -7.9214396353617822e-01, 2.8288585126419385e-01, -1.7119007404034725e-01, -4.2430521303823481e+00, 1.0084309674947685e+00, 6.9547355920839171e-01,
        -3.1347594048781219e-01, 1.5665697539380625e+01, -6.2531003217557455e+00, 1.6993190682438204e+00, -2.2418353145413406e+00, 7.8199634702149483e+00,
        -8.3115202323860138e-01, -2.1633011812923604e+00, 2.8539460745428218e+00, 1.2342612576538134e+00, -1.5842572057662735e-02, 1.5665697539380620e+01,
        -1.0905960194092088e+01, -2.2455097313438745e+00, 2.4716228830705056e-01, 1.3743240296982613e+01, -1.2482300064076051e+00, 2.8586237541840607e+00,
        -3.1464748454693747e-01, 1.1021105926083008e+00, 5.1991253846043550e-01, 1.5665697539380623e+01, -3.3874749822062640e+00, -5.4199310599802200e+00,
        3.6788734962989853e+00, 4.1719077295144302e+00, -3.4505908199700158e-01, 6.8997891471291757e+00, -4.6833531907539676e+00, 1.9863797621842311e-01,
        1.9863797621842255e-01, 1.5665697539380625e+01, 8.5106536630742688e+00, -4.1037281292243684e+00, 2.3300982732117133e+00, -1.0974886203719295e+01,
        1.0042363896175179e+00, 5.2242101413176201e+00, -2.9663083532485550e+00, -8.7326921436885940e-01, -4.6783739672943780e-01, 1.5665697539380618e+01,
        9.6405753039620379e+00, 1.6136475856828092e-01, -1.4711148381874981e+00, -1.2413321672478254e+01, 1.1762111052356046e+00, -2.0542379553857670e-01,
        1.8727880636075478e+00, -1.3079066293522887e+00, -3.0778800868070544e-01, 1.5665697539380631e+01, 2.7278721498501639e+00, 3.0162302076977188e+00,
        -1.7851165979825812e+00, -3.6131739081926835e+00, 4.1518369083042983e-01, -3.8397817651193273e+00, 2.2725248703007228e+00, -9.2431493744783222e-01,
        2.6140079875681571e-01, 1.5665697539380655e+01, 2.7278721498501568e+00, -3.0162302076977152e+00, 1.7851165979825847e+00, -3.6131739081926679e+00,
        4.1518369083043077e-01, 3.8397817651193265e+00, -2.2725248703007153e+00, -9.2431493744783277e-01, 2.6140079875681599e-01, 1.5665697539380639e+01,
        9.6405753039620414e+00, -1.6136475856828447e-01, 1.4711148381874963e+00, -1.2413321672478233e+01, 1.1762111052356037e+00, 2.0542379553857706e-01,
        -1.8727880636075500e+00, -1.3079066293522899e+00, -3.0778800868070444e-01, 1.5665697539380607e+01, 8.5106536630742688e+00, 4.1037281292243657e+00,
        -2.3300982732117159e+00, -1.0974886203719272e+01, 1.0042363896175219e+00, -5.2242101413176147e+00, 2.9663083532485537e+00, -8.7326921436885985e-01,
        -4.6783739672943597e-01, 1.5665697539380588e+01, -3.3874749822062782e+00, 5.4199310599802200e+00, -3.6788734962989835e+00, 4.1719077295144462e+00,
        -3.4505908199700419e-01, -6.8997891471291819e+00, 4.6833531907539623e+00, 1.9863797621842127e-01, 1.9863797621842422e-01, 1.5665697539380609e+01,
        -1.0905960194092074e+01, 2.2455097313438710e+00, -2.4716228830704523e-01, 1.3743240296982599e+01, -1.2482300064076046e+00, -2.8586237541840540e+00,
        3.1464748454693053e-01, 1.1021105926083008e+00, 5.1991253846043417e-01, 1.5665697539380636e+01, -6.2531003217557384e+00, -1.6993190682438204e+00,
        2.2418353145413406e+00, 7.8199634702149314e+00, -8.3115202323859982e-01, 2.1633011812923595e+00, -2.8539460745428191e+00, 1.2342612576538141e+00,
        -1.5842572057662315e-02, 1.5665697539380639e+01, -3.3256561883227675e-01, -3.3330076528705250e+00, 7.9214396353617644e-01, 2.8288585126417620e-01,
        -1.7119007404034572e-01, 4.2430521303823490e+00, -1.0084309674947725e+00, 6.9547355920839105e-01, -3.1347594048781230e-01, 1.5665697539380634e+01,
        2.8421709430404007e-14, -3.5592472331681222e+00, -3.4076656429846253e-02, -1.4048349091630685e-01, -8.1004038685081991e-16, 4.5310641702982473e+00,
        4.3380947396408673e-02, -1.2499260451994648e-01, 1.2499260451994688e-01, 1.5665697539380645e+01, -9.6405753039620237e+00, -1.4711148381875070e+00,
        1.6136475856828447e-01, 1.2132354690645675e+01, -1.1762111052356012e+00, 1.8727880636075578e+00, -2.0542379553857232e-01, 1.1311403928438450e+00,
        4.8455424518914980e-01, 1.5665697539380643e+01, -2.7278721498501284e+00, -4.6487098044534987e+00, 3.4175961947383620e+00, 3.3322069263601044e+00,
        -4.1518369083043150e-01, 5.9179936242654536e+00, -4.3507367294468464e+00, 4.9756349189949367e-01, 1.6535064679152167e-01, 1.5665697539380613e+01,
        8.1780880442419601e+00, -3.7874514632244480e+00, 2.0463938995432098e+00, -1.0551516861538797e+01, 8.3304631557717412e-01, 4.8215772879830503e+00,
        -2.6051413316078431e+00, -4.7955449618885987e-01, -4.7955449618885732e-01, 1.5665697539380586e+01, 8.5106536630742582e+00, -1.6257374643921096e-01,
        -1.6110561095734477e+00, -1.0974886203719283e+01, 1.0042363896175217e+00, 2.0696288548243486e-01, 2.0509389025866258e+00, -1.0500354508773055e+00,
        -2.9107116022099189e-01, 1.5665697539380597e+01, 1.4624872597200707e+00, 2.2418353145413406e+00, -1.6993190682438186e+00, -2.0022883018557285e+00,
        3.4316478965842906e-01, -2.8539460745428187e+00, 2.1633011812923577e+00, -9.5334473768337691e-01, 2.9675909202809997e-01, 1.5665697539380623e+01,
        -6.5960283235613559e-01, 2.7880259776413947e+00, 2.9535395799047315e-01, 6.9921731223805195e-01, 7.0124608833429983e-02, -3.5492686474345136e+00,
        -3.7599740870352399e-01, -4.2391812020101871e-01, 1.5827993394684667e-01, 1.5665697539380631e+01, 0.0000000000000000e+00, 3.0167309868705985e+00,
        -5.0843958986766680e-01, -1.4048349091628975e-01, -1.2979111161188702e-15, -3.8404192770477770e+00, 6.4726394585405544e-01, 3.0175884102839184e-01,
        -3.0175884102839107e-01, 1.5665697539380629e+01, -5.1231786808679658e+00, 2.0232575732513016e+00, -2.1018940431553936e+00, 6.3815280014559876e+00,
        -6.5917730762052162e-01, -2.5756878623133637e+00, 2.6757952355637480e+00, 9.7639007917883069e-01, -3.2559420517376217e-02, 1.5665697539380636e+01,
        -3.0549093633740156e+00, 3.5612118829268429e+00, -2.8726145195092361e+00, 3.7485383873339750e+00, -1.7386900795665677e-01, -4.5335652480671609e+00,
        3.6569532464990164e+00, -1.9507674196157757e-01, 2.1035507567784251e-01, 1.5665697539380638e+01, -8.1780880442419814e+00, 1.4787766771950537e+00,
        2.6228088648617742e-01, 1.0270549879706252e+01, -8.3304631557717312e-01, -1.8825418912939991e+00, -3.3389406508121872e-01, 4.7955449618885904e-01,
        4.7955449618885837e-01, 1.5665697539380629e+01, -2.7278721498501497e+00, -9.2492417508744751e-01, 2.1560377848025762e+00, 3.3322069263600911e+00,
        -4.1518369083043205e-01, 1.1774654907158566e+00, -2.7447223855344633e+00, 7.4754870093938719e-01, -8.4634562248372056e-02, 1.5665697539380616e+01,
        1.9249877224862075e+00, -1.4711148381874990e+00, 1.6136475856828802e-01, -2.5910699004075779e+00, 1.8942923385743744e-03, 1.8727880636075500e+00,
        -2.0542379553857354e-01, 6.2971415694500765e-01, -3.7040446372657621e-01, 1.5665697539380616e+01, -3.3256561883229097e-01, -1.7005280561147487e+00,
        -8.4033563321960258e-01, 2.8288585126418925e-01, -1.7119007404034650e-01, 2.1648402712362236e+00, 1.0697808916513587e+00, 2.6872211366005372e-01,
        1.1327550506052518e-01, 1.5665697539380611e+01, 0.0000000000000000e+00, -2.2494971535489121e+00, 1.2756734231893638e+00, -1.4048349091628354e-01,
        -1.9482695447842184e-16, 2.8636999022292722e+00, -1.6239833206725702e+00, -3.0175884102839207e-01, 3.0175884102839212e-01, 1.5665697539380615e+01,
        6.1153471320564172e+00, -9.3575965172466091e-01, 1.5569123679262589e+00, -7.9255651286233739e+00, 7.6024277282743702e-01, 1.1912594861150809e+00,
        -1.9820117526159078e+00, -8.2119407263786348e-01, -2.3899601848999721e-01, 1.5665697539380615e+01, 6.2531003217557384e+00, 2.2418353145413388e+00,
        -1.6993190682438248e+00, -8.1009304520474963e+00, 8.3115202323859849e-01, -2.8539460745428178e+00, 2.1633011812923590e+00, -8.0750981210547501e-01,
        -4.1090887349067573e-01, 1.5665697539380620e+01, 3.3874749822062569e+00, 1.5370291956072180e-01, -1.8947604832419573e+00, -4.4528747113470075e+00,
        3.4505908199700430e-01, -1.9566996785213062e-01, 2.4121059242273430e+00, -1.9863797621842222e-01, -1.9863797621842241e-01, 1.5665697539380620e+01,
        -3.0549093633740227e+00, 3.7994249985831807e-01, -1.0685398632759249e+00, 3.7485383873339648e+00, -1.7386900795665666e-01, -4.8368200776803066e-01,
        1.3602940093361668e+00, -3.7184297847002262e-01, 3.8712131218628748e-01, 1.5665697539380634e+01, -1.9249877224862217e+00, -1.6136475856827737e-01,
        1.4711148381874981e+00, 2.3101029185750024e+00, -1.8942923385708280e-03, 2.0542379553857140e-01, -1.8727880636075536e+00, -4.5294792043656368e-01,
        1.9363822721812987e-01, 1.5665697539380638e+01, 2.7278721498501355e+00, 7.0755542166834040e-01, 5.2355818804680254e-01, -3.6131739081926684e+00,
        4.1518369083042905e-01, -9.0074636843026978e-01, -6.6651052638833552e-01, -3.2079725539105020e-01, -3.4211688329996731e-01, 1.5665697539380632e+01,
        0.0000000000000000e+00, 1.7069809072513848e+00, -1.8181896694868858e+00, -1.4048349091628820e-01, -4.7951458347135731e-16, -2.1730550089788068e+00,
        2.3146282139230339e+00, 1.2499260451994659e-01, -1.2499260451994658e-01, 1.5665697539380618e+01, -5.1231786808679871e+00, 3.9077797649552615e-01,
        -4.6941444639961638e-01, 6.3815280014559947e+00, -6.5917730762052062e-01, -4.9747600316724699e-01, 5.9758337641762205e-01, 5.4963863363049015e-01,
        3.9419202503096307e-01, 1.5665697539380618e+01, -1.4624872597200849e+00, -1.6993190682438257e+00, 2.2418353145413406e+00, 1.7213213200231645e+00,
        -3.4316478965842862e-01, 2.1633011812923599e+00, -2.8539460745428191e+00, 5.2659329213503836e-01, 1.2999235352023908e-01, 1.5665697539380609e+01,
        5.4502158943918175e+00, -1.4782758980221766e+00, 1.0143961216287396e+00, -7.0788264442624165e+00, 4.1786262474674263e-01, 1.8819043793655328e+00,
        -1.2913668593654577e+00, 1.4300160023058200e-01, -4.3919645391728296e-01, 1.5665697539380606e+01,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_DxDx_DyDy(void)
{
  static const double sol[640] = {
        -1.6736386918814842e-04, -8.9993566755474939e-05, -8.9993566755553327e-05, 2.1306090344382802e-04, -4.8127663587508573e-05, 1.1456541205790448e-04,
        1.1456541205785753e-04, 5.1386806341961409e-06, 5.1386806342542322e-06, 1.3021735293941976e-16, -4.0904809409466058e-04, -2.1943276528524800e-04,
        -4.9944248969825226e-05, 5.2073459404767034e-04, -2.6835819803993821e-05, 2.7934669199395864e-04, 6.3581027727126449e-05, 1.2959979779996389e-05,
        4.8277858685092419e-06, 2.5963965559206471e-18, -5.1818273509544963e-04, -2.7802320360133129e-04, -1.3907438669093784e-05, 6.5966735965373763e-04,
        -7.4779609092899315e-06, 3.5393466478226404e-04, 1.7704726007174729e-05, 1.6382516606644185e-05, 4.1749764935150647e-06, -6.3547404895491045e-17,
        -5.1295357548712922e-04, -2.7521036368747691e-04, 1.5668491170627079e-05, 6.5301043019852262e-04, 8.3756457236890595e-06, 3.5035380700139911e-04,
        -1.9946616319957822e-05, 1.6222781764510503e-05, 3.2005306256819650e-06, 2.2995172754667680e-16, -4.2024623592899952e-04, -2.2547403683496499e-04,
        3.6177829984709971e-05, 5.3499027675668620e-04, 1.9370676929407923e-05, 2.8703747245083336e-04, -4.6055825421175033e-05, 1.3288392231881120e-05,
        1.9005356520037570e-06, 3.1362767201675851e-16, -2.7595668815097143e-04, -1.4805299908269591e-04, 4.5030778835311037e-05, 3.5130390791152642e-04,
        2.4117676294458697e-05, 1.8847739297140733e-04, -5.7325983606757273e-05, 8.7302227497018262e-06, 2.7620610285165812e-07, -2.5080839944164045e-17,
        -1.2491297408360295e-04, -6.7053537941976415e-05, 3.9653588273064960e-05, 1.5901921507471145e-04, 2.1192021856326956e-05, 8.5361837309139718e-05,
        -5.0480604823932580e-05, 3.9233331137339112e-06, -1.6763676439127430e-06, -1.9950606076705151e-17, -2.1186748034872204e-05, -1.1173638897301938e-05,
        1.7314566130892644e-05, 2.6971578150237605e-05, 9.3854234968676341e-06, 1.4224489489855110e-05, -2.2042135620260495e-05, 8.1993829229793671e-07,
        -3.9368034886668810e-06, 1.3826372781016619e-17, -4.0904809409451714e-04, -4.9944248969747733e-05, -2.1943276528526315e-04, 5.2073459404797934e-04,
        -2.6835819804019046e-05, 6.3581027727065517e-05, 2.7934669199371020e-04, 4.8277858685077774e-06, 1.2959979779994376e-05, 6.9723122910651338e-17,
        -9.9858541584180540e-04, -1.2197991056199567e-04, -1.2197991056206587e-04, 1.2712391981478468e-03, -1.4892793976961156e-05, 1.5528530782951420e-04,
        1.5528530782953203e-04, 1.2151963571855302e-05, 1.2151963571851188e-05, 5.5315123985650645e-17, -1.2650670719583249e-03, -1.5452803962138343e-04,
        -3.3989146299085161e-05, 1.6104810110849606e-03, -4.1472845964592906e-06, 1.9672037871097946e-04, 4.3269543497749215e-05, 1.5346594041114546e-05,
        1.0513071189119953e-05, 8.1224921361588226e-17, -1.2522942054227102e-03, -1.5296725544019501e-04, 3.8195274286687170e-05, 1.5942206408098643e-03,
        4.6718088449734621e-06, 1.9473343798517319e-04, -4.8624112756910573e-05, 1.5201005557162338e-05, 8.0589243504866540e-06, 1.2646430323484297e-16,
        -1.0259673331826800e-03, -1.2532171759625012e-04, 8.8249173983067660e-05, 1.3060974747577253e-03, 1.0787713458074866e-05, 1.5953956192457048e-04,
        -1.1234473024707882e-04, 1.2449879038406603e-05, 4.7859185315923841e-06, 1.2492617262243828e-16, -6.7370284092416264e-04, -8.2291875715507218e-05,
        1.0985003716686473e-04, 8.5765067834920466e-04, 1.3428120325720441e-05, 1.0476085113908371e-04, -1.3984349355511364e-04, 8.1833370321413966e-06,
        6.9543751779662684e-07, 7.0997481676593796e-17, -3.0500488785042019e-04, -3.7252512702599006e-05, 9.6677053695450602e-05, 3.8828342864919080e-04,
        1.1825246932944356e-05, 4.7423939530637341e-05, -1.2307375840803705e-04, 3.6677659788390534e-06, -4.2161030397036606e-06, -3.8675472809755950e-18,
        -5.1280526705659661e-05, -6.2839193715145735e-06, 4.2430283014355643e-05, 6.5282162763118009e-05, 5.1630088156727532e-06, 7.9996808448424854e-06,
        -5.4015448353894996e-05, 7.5571876572633797e-07, -9.9330071707124047e-06, -8.6830981711185261e-17, -5.1818273509495415e-04, -1.3907438669079623e-05,
        -2.7802320360142870e-04, 6.5966735965313611e-04, -7.4779609092294263e-06, 1.7704726006891892e-05, 3.5393466478249959e-04, 4.1749764934271104e-06,
        1.6382516606623020e-05, 4.0749186165041001e-17, -1.2650670719581497e-03, -3.3989146299039705e-05, -1.5452803962143685e-04, 1.6104810110848080e-03,
        -4.1472845964594609e-06, 4.3269543497699613e-05, 1.9672037871103638e-04, 1.0513071189135437e-05, 1.5346594041119834e-05, 7.9240261984326072e-17,
        -1.6026670347201000e-03, -4.3053493382510712e-05, -4.3053493382558044e-05, 2.0402592745641132e-03, -1.1564832204205967e-06, 5.4808820093699731e-05,
        5.4808820093733822e-05, 1.3279889860017377e-05, 1.3279889860022863e-05, 9.2097956645809326e-17, -1.5864825433925698e-03, -4.2619862890565045e-05,
        4.8395012432343505e-05, 2.0196557693944101e-03, 1.3017448917312349e-06, 5.4256791123399684e-05, -6.1608787613830381e-05, 1.3153443349271750e-05,
        1.0179488166753925e-05, 1.0205130438863518e-16, -1.2997582939416842e-03, -3.4916529755849637e-05, 1.1180845736502144e-04, 1.6546443250261163e-03,
        3.0057340765579082e-06, 4.4450139752483907e-05, -1.4233664084413466e-04, 1.0772906330799025e-05, 6.0455797587637271e-06, 9.1330888093378114e-17,
        -8.5348426805717222e-04, -2.2928719752622003e-05, 1.3917704869122568e-04, 1.0865196300130122e-03, 3.7402936569764150e-06, 2.9189177861487167e-05,
        -1.7717795290418068e-04, 7.0805126791898829e-06, 8.7830508376714211e-07, 5.6615583707764242e-17, -3.8639859219410038e-04, -1.0375616284744785e-05,
        1.2249026701637098e-04, 4.9190087168677211e-04, 3.2966811714911069e-06, 1.3208574766734493e-05, -1.5593501202052922e-04, 3.1755470736365147e-06,
        -5.3221282290210776e-06, 5.8555555451602192e-18, -6.4987942613438211e-05, -1.7584195977303762e-06, 5.3733951100266865e-05, 8.2732251790000952e-05,
        1.4372749293382423e-06, 2.2385385205262093e-06, -6.8405470204461853e-05, 6.5629254036916556e-07, -1.2558648067133677e-05, -6.2643782381705725e-17,
        -5.1295357548637873e-04, 1.5668491170882496e-05, -2.7521036368765222e-04, 6.5301043019802346e-04, 8.3756457234651734e-06, -1.9946616319747704e-05,
        3.5035380700151241e-04, 3.2005306256888069e-06, 1.6222781764533938e-05, 1.1021652978936597e-16, -1.2522942054223754e-03, 3.8195274286697429e-05,
        -1.5296725544028156e-04, 1.5942206408096325e-03, 4.6718088449811464e-06, -4.8624112756891437e-05, 1.9473343798523865e-04, 8.0589243504738045e-06,
        1.5201005557168734e-05, 8.5918816888484507e-17, -1.5864825433924558e-03, 4.8395012432365846e-05, -4.2619862890609910e-05, 2.0196557693943381e-03,
        1.3017448917293096e-06, -6.1608787613840504e-05, 5.4256791123429228e-05, 1.0179488166751226e-05, 1.3153443349275177e-05, 8.8436997644373772e-17,
        -1.5704614601383636e-03, 4.7905128339743241e-05, 4.7905128339721455e-05, 1.9992602892418077e-03, -1.4574932682125792e-06, -6.0985145558572417e-05,
        -6.0985145558559136e-05, 1.0082546116344673e-05, 1.0082546116347507e-05, 8.7079524843457519e-17, -1.2866320417064825e-03, 3.9248066346110315e-05,
        1.1067785968105715e-04, 1.6379340806129848e-03, -3.3712632304338978e-06, -4.9964359181647862e-05, -1.4089734474548839e-04, 8.2578616887689264e-06,
        5.9880257100210153e-06, 7.5321941371283178e-17, -8.4486437881703771e-04, 2.5771383674393462e-05, 1.3776937488854907e-04, 1.0755461660389958e-03,
        -4.1977162667747768e-06, -3.2808002798411447e-05, -1.7538592781772406e-04, 5.4274990248750375e-06, 8.6992119349358339e-07, 5.0348534537109036e-17,
        -3.8249295195446459e-04, 1.1673090903164162e-05, 1.2125163449742982e-04, 4.8692883535647754e-04, -3.6878033555764420e-06, -1.4860311881417390e-05,
        -1.5435818325335691e-04, 2.4339111771927973e-06, -5.2720559362171523e-06, 1.5157875640171339e-17, -6.4327335572091798e-05, 1.9416918503506198e-06,
        5.3193484611823988e-05, 8.1891272588506745e-05, -1.6349233391438523e-06, -2.4718514327569146e-06, -6.7717434735062581e-05, 5.0282216664009615e-07,
        -1.2435436975576558e-05, -9.9002448715467669e-18, -4.2024623592846105e-04, 3.6177829984429400e-05, -2.2547403683507238e-04, 5.3499027675656065e-04,
        1.9370676929591888e-05, -4.6055825421170601e-05, 2.8703747245078457e-04, 1.9005356520859239e-06, 1.3288392231919035e-05, 9.9271504070519187e-17,
        -1.0259673331823702e-03, 8.8249173983027626e-05, -1.2532171759632019e-04, 1.3060974747575818e-03, 1.0787713458090848e-05, -1.1234473024706425e-04,
        1.5953956192458536e-04, 4.7859185315675474e-06, 1.2449879038409287e-05, 7.5617121957013773e-17, -1.2997582939415506e-03, 1.1180845736501496e-04,
        -3.4916529755881803e-05, 1.6546443250260369e-03, 3.0057340765646887e-06, -1.4233664084413027e-04, 4.4450139752502529e-05, 6.0455797587646961e-06,
        1.0772906330801502e-05, 7.4237365136813508e-17, -1.2866320417064361e-03, 1.1067785968106318e-04, 3.9248066346090854e-05, 1.6379340806129577e-03,
        -3.3712632304321939e-06, -1.4089734474549226e-04, -4.9964359181636593e-05, 5.9880257100190332e-06, 8.2578616887703291e-06, 7.2093480822636435e-17,
        -1.0540980928973197e-03, 9.0676106148561848e-05, 9.0676106148553960e-05, 1.3419090576788013e-03, -7.7938401140048707e-06, -1.1543431202056919e-04,
        -1.1543431202056450e-04, 4.9043496343493622e-06, 4.9043496343501398e-06, 6.3037991309102385e-17, -6.9217011202356304e-04, 5.9541218839865003e-05,
        1.1287158611772119e-04, 8.8116025352627286e-04, -9.7032178227066747e-06, -7.5798354446155174e-05, -1.4369004629315670e-04, 3.2234047118027443e-06,
        7.1260045700449238e-07, 4.4875646352635977e-17, -3.1336327392571109e-04, 2.6963244897811973e-05, 9.9338449335989494e-05, 3.9892398862896702e-04,
        -8.5312831858724501e-06, -3.4325289834583248e-05, -1.2646182156855073e-04, 1.4457563333310716e-06, -4.3177809666449573e-06, 2.0591662977301158e-17,
        -5.2701904096224544e-05, 4.5069561575177094e-06, 4.3578076238989909e-05, 6.7091633065451405e-05, -3.7645201112378417e-06, -5.7375355586665948e-06,
        -5.5476635063979589e-05, 2.9887242508465092e-07, -1.0186424686089527e-05, -4.2994398398959147e-18, -2.7595668815092784e-04, 4.5030778835202332e-05,
        -1.4805299908260343e-04, 3.5130390791139420e-04, 2.4117676294505365e-05, -5.7325983606700216e-05, 1.8847739297142972e-04, 2.7620610285355521e-07,
        8.7302227496579397e-06, 2.2068417882700666e-17, -6.7370284092403904e-04, 1.0985003716684141e-04, -8.2291875715501757e-05, 8.5765067834911381e-04,
        1.3428120325718364e-05, -1.3984349355510245e-04, 1.0476085113909073e-04, 6.9543751780616676e-07, 8.1833370321377933e-06, 4.8437492265353730e-17,
        -8.5348426805707486e-04, 1.3917704869120990e-04, -2.2928719752636328e-05, 1.0865196300129493e-03, 3.7402936569831785e-06, -1.7717795290417499e-04,
        2.9189177861496532e-05, 8.7830508376708133e-07, 7.0805126791911475e-06, 5.2844436006943742e-17, -8.4486437881698816e-04, 1.3776937488854695e-04,
        2.5771383674380513e-05, 1.0755461660389629e-03, -4.1977162667738222e-06, -1.7538592781772452e-04, -3.2808002798403451e-05, 8.6992119349322436e-07,
        5.4274990248764156e-06, 5.3857020936747595e-17, -6.9217011202354765e-04, 1.1287158611772247e-04, 5.9541218839857820e-05, 8.8116025352626148e-04,
        -9.7032178227062241e-06, -1.4369004629315852e-04, -7.5798354446150322e-05, 7.1260045700381200e-07, 3.2234047118037552e-06, 4.7224205907133811e-17,
        -4.5451002284198453e-04, 7.4115147977006358e-05, 7.4115147977004203e-05, 5.7860944874782883e-04, -1.2080644800994201e-05, -9.4351549492118655e-05,
        -9.4351549492116771e-05, 4.6835828133488412e-07, 4.6835828133561553e-07, 3.4282674743364552e-17, -2.0576429737189132e-04, 3.3563630657697868e-05,
        6.5229065397315528e-05, 2.6194618532257558e-04, -1.0621583438880923e-05, -4.2727845056957630e-05, -8.3039210743648229e-05, 2.0994736147339859e-07,
        -2.8385484135676012e-06, 1.7422616860709244e-17, -3.4602943899175216e-05, 5.6125950401316104e-06, 2.8616778662176850e-05, 4.4050932407094513e-05,
        -4.6829279478430336e-06, -7.1450581043514569e-06, -3.6430304491689317e-05, 4.3342012469426337e-08, -6.6915544801017071e-06, -3.5904178472779480e-18,
        -1.2491297408348596e-04, 3.9653588273158940e-05, -6.7053537942000783e-05, 1.5901921507466678e-04, 2.1192021856246867e-05, -5.0480604823934877e-05,
        8.5361837309134270e-05, -1.6763676439549742e-06, 3.9233331137395939e-06, 4.4671544762584355e-18, -3.0500488785035785e-04, 9.6677053695436562e-05,
        -3.7252512702608323e-05, 3.8828342864913893e-04, 1.1825246932935438e-05, -1.2307375840802916e-04, 4.7423939530636514e-05, -4.2161030396958425e-06,
        3.6677659788430429e-06, 2.4128934571038136e-17, -3.8639859219403587e-04, 1.2249026701636486e-04, -1.0375616284750664e-05, 4.9190087168672289e-04,
        3.2966811714913877e-06, -1.5593501202052686e-04, 1.3208574766736104e-05, -5.3221282290233434e-06, 3.1755470736278715e-06, 3.1432284870142483e-17,
        -3.8249295195442838e-04, 1.2125163449742528e-04, 1.1673090903159696e-05, 4.8692883535644642e-04, -3.6878033555725236e-06, -1.5435818325335593e-04,
        -1.4860311881412347e-05, -5.2720559362166873e-06, 2.4339111771886909e-06, 3.2077275460478534e-17, -3.1336327392569852e-04, 9.9338449335986784e-05,
        2.6963244897806884e-05, 3.9892398862895548e-04, -8.5312831858720622e-06, -1.2646182156854913e-04, -3.4325289834580172e-05, -4.3177809666453419e-06,
        1.4457563333301612e-06, 2.8122834349512110e-17, -2.0576429737189064e-04, 6.5229065397313983e-05, 3.3563630657696838e-05, 2.6194618532257119e-04,
        -1.0621583438879969e-05, -8.3039210743649097e-05, -4.2727845056955679e-05, -2.8385484135681285e-06, 2.0994736147304273e-07, 2.0500639862063439e-17,
        -9.3151209302059127e-05, 2.9542138267396792e-05, 2.9542138267397802e-05, 1.1858521738957481e-04, -9.3312153718748475e-06, -3.7608324302396471e-05,
        -3.7608324302396125e-05, -1.2706382352886107e-06, -1.2706382352886850e-06, 1.0852179600489556e-17, -1.5684011184944028e-05, 4.9183748803106635e-06,
        1.2939562203297638e-05, 1.9966373918734288e-05, -4.1420646084805540e-06, -6.2612880579344652e-06, -1.6472580531169222e-05, -2.6020768542684102e-07,
        -3.0088387502362375e-06, 2.9995937327780627e-18, -2.1186748034736431e-05, 1.7314566130796099e-05, -1.1173638897437015e-05, 2.6971578150202090e-05,
        9.3854234968672513e-06, -2.2042135620263812e-05, 1.4224489489858336e-05, -3.9368034887195114e-06, 8.1993829231012340e-07, 7.6416996781065047e-18,
        -5.1280526705648921e-05, 4.2430283014323821e-05, -6.2839193714807337e-06, 6.5282162763124067e-05, 5.1630088157216431e-06, -5.4015448353865011e-05,
        7.9996808448244386e-06, -9.9330071707169702e-06, 7.5571876571604143e-07, 9.7641800636517703e-18, -6.4987942613452184e-05, 5.3733951100238621e-05,
        -1.7584195977080992e-06, 8.2732251789935588e-05, 1.4372749293546847e-06, -6.8405470204479200e-05, 2.2385385205287847e-06, -1.2558648067127453e-05,
        6.5629254041339423e-07, 9.0893247680704442e-18, -6.4327335572037439e-05, 5.3193484611841084e-05, 1.9416918503227109e-06, 8.1891272588485915e-05,
        -1.6349233391667857e-06, -6.7717434735056211e-05, -2.4718514327320817e-06, -1.2435436975574417e-05, 5.0282216664938090e-07, 9.9415504832613775e-18,
        -5.2701904096233008e-05, 4.3578076238977095e-05, 4.5069561575094119e-06, 6.7091633065443219e-05, -3.7645201112378235e-06, -5.5476635063977041e-05,
        -5.7375355586687267e-06, -1.0186424686092530e-05, 2.9887242509399168e-07, 8.8524613742982187e-18, -3.4602943899185231e-05, 2.8616778662172665e-05,
        5.6125950401316646e-06, 4.4050932407092101e-05, -4.6829279478433995e-06, -3.6430304491688470e-05, -7.1450581043495511e-06, -6.6915544801019002e-06,
        4.3342012473073687e-08, 6.6378527613512372e-18, -1.5684011184945674e-05, 1.2939562203297628e-05, 4.9183748803111208e-06, 1.9966373918733990e-05,
        -4.1420646084811334e-06, -1.6472580531168789e-05, -6.2612880579338876e-06, -3.0088387502362185e-06, -2.6020768542654058e-07, 3.8580014055187301e-18,
        -2.4631643688340385e-06, 2.2363599383471548e-06, 2.2363599383472056e-06, 3.1357068183325297e-06, -1.7612712352079330e-06, -2.8469757015205160e-06,
        -2.8469757015205410e-06, -6.3071636577514211e-07, -6.3071636577518499e-07, 1.2372631589442725e-18,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_DxDx_PyPy(void)
{
  static const double sol[640] = {
        5.8616184533469527e-03, 3.1565760734898127e-03, 2.7094777984668392e-03, -7.4620748753877917e-03, 1.4328949642727182e-03, -4.0184476689694810e-03,
        -3.4492736718164069e-03, -1.7632215123710278e-04, -5.2264298589365711e-04, 3.5131649393126775e-17, 5.8616184533454548e-03, 3.1565760734899480e-03,
        -2.7094777984674428e-03, -7.4620748753869296e-03, -1.4328949642726360e-03, -4.0184476689695782e-03, 3.4492736718167269e-03, -1.7632215123710747e-04,
        -5.2264298589345328e-04, 1.3113006726376057e-16, -5.8616184533450645e-03, -3.1565760734899063e-03, -2.7094777984659453e-03, 7.4620748753872505e-03,
        -1.4328949642727893e-03, 4.0184476689695140e-03, 3.4492736718158296e-03, 1.7632215123710446e-04, 5.2264298589332133e-04, -6.7007510729665259e-16,
        -5.8616184533434972e-03, -3.1565760734900798e-03, 2.7094777984662758e-03, 7.4620748753870840e-03, 1.4328949642727672e-03, 4.0184476689695556e-03,
        -3.4492736718152593e-03, 1.7632215123710763e-04, 5.2264298589287594e-04, -1.2842098510653892e-15, 5.8616184533454565e-03, 3.1565760734899358e-03,
        2.7094777984669749e-03, -7.4620748753874231e-03, 1.4328949642725803e-03, -4.0184476689694975e-03, -3.4492736718174520e-03, -1.7632215123710847e-04,
        -5.2264298589330311e-04, 2.1405198315438111e-16, 5.8616184533454088e-03, 3.1565760734899601e-03, -2.7094777984659605e-03, -7.4620748753875914e-03,
        -1.4328949642728730e-03, -4.0184476689694984e-03, 3.4492736718165261e-03, -1.7632215123710143e-04, -5.2264298589249668e-04, 4.4941804705996507e-16,
        -5.8616184533449344e-03, -3.1565760734899059e-03, -2.7094777984672976e-03, 7.4620748753882150e-03, -1.4328949642725777e-03, 4.0184476689694298e-03,
        3.4492736718163448e-03, 1.7632215123710584e-04, 5.2264298589430025e-04, -1.1539732792029041e-15, -5.8616184533448580e-03, -3.1565760734899163e-03,
        2.7094777984665547e-03, 7.4620748753871950e-03, 1.4328949642728101e-03, 4.0184476689695513e-03, -3.4492736718163097e-03, 1.7632215123710075e-04,
        5.2264298589368996e-04, -8.7123700173506059e-16, 1.4341205890542634e-02, 1.7536458636411641e-03, 6.5546632008048515e-03, -1.8256929039364878e-02,
        7.9551267430263116e-04, -2.2324613660129437e-03, -8.3443485748260857e-03, -1.6506737964353752e-04, -1.3363619667712365e-03, -1.1204417257716884e-16,
        1.4341205890541529e-02, 1.7536458636412541e-03, -6.5546632008052566e-03, -1.8256929039364305e-02, -7.9551267430259690e-04, -2.2324613660130101e-03,
        8.3443485748262783e-03, -1.6506737964354181e-04, -1.3363619667710529e-03, -2.9814579716794004e-17, -1.4341205890541050e-02, -1.7536458636412450e-03,
        -6.5546632008041819e-03, 1.8256929039364444e-02, -7.9551267430269512e-04, 2.2324613660129719e-03, 8.3443485748256590e-03, 1.6506737964353847e-04,
        1.3363619667711248e-03, -4.3980537973069959e-16, -1.4341205890539973e-02, -1.7536458636413575e-03, 6.5546632008044474e-03, 1.8256929039364378e-02,
        7.9551267430268016e-04, 2.2324613660129871e-03, -8.3443485748254404e-03, 1.6506737964354083e-04, 1.3363619667708508e-03, -8.3643101508846570e-16,
        1.4341205890541506e-02, 1.7536458636412552e-03, 6.5546632008047102e-03, -1.8256929039364576e-02, 7.9551267430260947e-04, -2.2324613660129641e-03,
        -8.3443485748265888e-03, -1.6506737964354324e-04, -1.3363619667710507e-03, -4.1976905387266643e-17, 1.4341205890541499e-02, 1.7536458636412508e-03,
        -6.5546632008043762e-03, -1.8256929039364728e-02, -7.9551267430271193e-04, -2.2324613660129507e-03, 8.3443485748262540e-03, -1.6506737964353831e-04,
        -1.3363619667707905e-03, 1.2001911151799687e-16, -1.4341205890540923e-02, -1.7536458636412493e-03, -6.5546632008049869e-03, 1.8256929039365152e-02,
        -7.9551267430258573e-04, 2.2324613660129086e-03, 8.3443485748260250e-03, 1.6506737964353893e-04, 1.3363619667713777e-03, -7.5128207929077665e-16,
        -1.4341205890540861e-02, -1.7536458636412447e-03, 6.5546632008047874e-03, 1.8256929039364479e-02, 7.9551267430266900e-04, 2.2324613660129867e-03,
        -8.3443485748261013e-03, 1.6506737964353663e-04, 1.3363619667711953e-03, -5.9938370361699392e-16, 1.8173569455656034e-02, 4.8776519957604707e-04,
        8.2864356215894778e-03, -2.3135681230452648e-02, 2.2095135802113725e-04, -6.2094461961560558e-04, -1.0548964172698889e-02, -1.4276292489187592e-04,
        -1.7088219069646235e-03, -2.2100584181581480e-16, 1.8173569455655188e-02, 4.8776519957610518e-04, -8.2864356215897883e-03, -2.3135681230452277e-02,
        -2.2095135802111752e-04, -6.2094461961565469e-04, 1.0548964172699016e-02, -1.4276292489187839e-04, -1.7088219069644811e-03, -1.4876109397575244e-16,
        -1.8173569455654716e-02, -4.8776519957612248e-04, -8.2864356215890025e-03, 2.3135681230452315e-02, -2.2095135802118989e-04, 6.2094461961563778e-04,
        1.0548964172698589e-02, 1.4276292489187435e-04, 1.7088219069645867e-03, -2.4290390950766229e-16, -1.8173569455653963e-02, -4.8776519957620059e-04,
        8.2864356215892210e-03, 2.3135681230452270e-02, 2.2095135802117612e-04, 6.2094461961563290e-04, -1.0548964172698514e-02, 1.4276292489187487e-04,
        1.7088219069644197e-03, -4.9858079137137815e-16, 1.8173569455655174e-02, 4.8776519957611424e-04, 8.2864356215893130e-03, -2.3135681230452416e-02,
        2.2095135802114696e-04, -6.2094461961562466e-04, -1.0548964172699123e-02, -1.4276292489187687e-04, -1.7088219069645082e-03, -2.1645988308624400e-16,
        1.8173569455655150e-02, 4.8776519957610638e-04, -8.2864356215892158e-03, -2.3135681230452523e-02, -2.2095135802118319e-04, -6.2094461961560764e-04,
        1.0548964172699019e-02, -1.4276292489187549e-04, -1.7088219069644247e-03, -9.9140963264315763e-17, -1.8173569455654598e-02, -4.8776519957612372e-04,
        -8.2864356215895125e-03, 2.3135681230452829e-02, -2.2095135802112378e-04, 6.2094461961558758e-04, 1.0548964172698830e-02, 1.4276292489187435e-04,
        1.7088219069646298e-03, -4.4848074315116613e-16, -1.8173569455654511e-02, -4.8776519957611250e-04, 8.2864356215895073e-03, 2.3135681230452385e-02,
        2.2095135802115403e-04, 6.2094461961564071e-04, -1.0548964172698927e-02, 1.4276292489187402e-04, 1.7088219069645871e-03, -3.8396631432671612e-16,
        1.7990573879785168e-02, -5.5033150535898810e-04, 8.1929460267051273e-03, -2.2902720538814127e-02, -2.5005174898973937e-04, 7.0059403080543526e-04,
        -1.0429948176920740e-02, -1.0938846751228960e-04, -1.6994005263250107e-03, -2.9381282230482898e-16, 1.7990573879784484e-02, -5.5033150535894278e-04,
        -8.1929460267053771e-03, -2.2902720538813877e-02, 2.5005174898975357e-04, 7.0059403080541163e-04, 1.0429948176920835e-02, -1.0938846751228980e-04,
        -1.6994005263248999e-03, -2.1887573342775259e-16, -1.7990573879784102e-02, 5.5033150535891903e-04, -8.1929460267048202e-03, 2.2902720538813891e-02,
        2.5005174898969345e-04, -7.0059403080540827e-04, 1.0429948176920542e-02, 1.0938846751228834e-04, 1.6994005263249875e-03, -9.5982628332163999e-17,
        -1.7990573879783599e-02, 5.5033150535885561e-04, 8.1929460267049902e-03, 2.2902720538813839e-02, -2.5005174898970787e-04, -7.0059403080540989e-04,
        -1.0429948176920533e-02, 1.0938846751229020e-04, 1.6994005263248945e-03, -2.6117596237755658e-16, 1.7990573879784501e-02, -5.5033150535894180e-04,
        8.1929460267050076e-03, -2.2902720538813939e-02, -2.5005174898971990e-04, 7.0059403080542258e-04, -1.0429948176920834e-02, -1.0938846751229017e-04,
        -1.6994005263249294e-03, -2.9830977257397701e-16, 1.7990573879784477e-02, -5.5033150535894039e-04, -8.1929460267049954e-03, -2.2902720538814009e-02,
        2.5005174898970760e-04, 7.0059403080542919e-04, 1.0429948176920820e-02, -1.0938846751228923e-04, -1.6994005263249057e-03, -2.2477885568782595e-16,
        -1.7990573879783998e-02, 5.5033150535891339e-04, -8.1929460267051394e-03, 2.2902720538814262e-02, 2.5005174898973828e-04, -7.0059403080544437e-04,
        1.0429948176920719e-02, 1.0938846751228730e-04, 1.6994005263249799e-03, -2.4211931635057277e-16, -1.7990573879783870e-02, 5.5033150535892814e-04,
        8.1929460267052071e-03, 2.2902720538813964e-02, -2.5005174898972549e-04, -7.0059403080541304e-04, -1.0429948176920808e-02, 1.0938846751228806e-04,
        1.6994005263249764e-03, -2.2479850252782248e-16, 1.4738458963226696e-02, -1.2698766175711731e-03, 6.7039542837525066e-03, -1.8762648099116208e-02,
        -5.7647737886046312e-04, 1.6166037551301662e-03, -8.5344020981064049e-03, -6.4925612047121364e-05, -1.3983785661310436e-03, -3.1436539032981156e-16,
        1.4738458963226148e-02, -1.2698766175711388e-03, -6.7039542837527217e-03, -1.8762648099116017e-02, 5.7647737886047179e-04, 1.6166037551301536e-03,
        8.5344020981064830e-03, -6.4925612047123370e-05, -1.3983785661309629e-03, -2.2651529080434071e-16, -1.4738458963225858e-02, 1.2698766175711083e-03,
        -6.7039542837523470e-03, 1.8762648099116059e-02, 5.7647737886042745e-04, -1.6166037551301447e-03, 8.5344020981062939e-03, 6.4925612047122083e-05,
        1.3983785661310230e-03, -2.6426688005693668e-18, -1.4738458963225537e-02, 1.2698766175710658e-03, 6.7039542837524684e-03, 1.8762648099116003e-02,
        -5.7647737886044100e-04, -1.6166037551301432e-03, -8.5344020981063112e-03, 6.4925612047124061e-05, 1.3983785661309737e-03, -1.0693961517654768e-16,
        1.4738458963226163e-02, -1.2698766175711421e-03, 6.7039542837524580e-03, -1.8762648099116073e-02, -5.7647737886044490e-04, 1.6166037551301495e-03,
        -8.5344020981064517e-03, -6.4925612047124346e-05, -1.3983785661309850e-03, -2.9001981455784633e-16, 1.4738458963226148e-02, -1.2698766175711347e-03,
        -6.7039542837524701e-03, -1.8762648099116111e-02, 5.7647737886044295e-04, 1.6166037551301597e-03, 8.5344020981064691e-03, -6.4925612047121879e-05,
        -1.3983785661309876e-03, -2.6297640124286923e-16, -1.4738458963225773e-02, 1.2698766175711076e-03, -6.7039542837525265e-03, 1.8762648099116309e-02,
        5.7647737886046138e-04, -1.6166037551301768e-03, 8.5344020981064049e-03, 6.4925612047123655e-05, 1.3983785661310136e-03, -1.1897205361844484e-16,
        -1.4738458963225598e-02, 1.2698766175711174e-03, 6.7039542837526315e-03, 1.8762648099116094e-02, -5.7647737886045455e-04, -1.6166037551301579e-03,
        -8.5344020981064795e-03, 6.4925612047122815e-05, 1.3983785661310164e-03, -1.2970455818161604e-16, 9.6780483124899331e-03, -1.5800101497750775e-03,
        4.3916291488594571e-03, -1.2320542821102249e-02, -7.1712350226570088e-04, 2.0114161532917741e-03, -5.5907196612254642e-03, -9.3539679266368939e-06,
        -9.2641347200032116e-04, -2.8000043592001863e-16, 9.6780483124894682e-03, -1.5800101497750636e-03, -4.3916291488596462e-03, -1.2320542821102098e-02,
        7.1712350226570555e-04, 2.0114161532917641e-03, 5.5907196612255249e-03, -9.3539679266382915e-06, -9.2641347200025893e-04, -1.7598141661156963e-16,
        -9.6780483124893086e-03, 1.5800101497750179e-03, -4.3916291488593964e-03, 1.2320542821102162e-02, 7.1712350226567508e-04, -2.0114161532917589e-03,
        5.5907196612253905e-03, 9.3539679266371463e-06, 9.2641347200030392e-04, 4.0721693196094477e-17, -9.6780483124891143e-03, 1.5800101497749840e-03,
        4.3916291488594753e-03, 1.2320542821102124e-02, -7.1712350226568712e-04, -2.0114161532917519e-03, -5.5907196612254156e-03, 9.3539679266377680e-06,
        9.2641347200027888e-04, -2.4654046496217829e-17, 9.6780483124894769e-03, -1.5800101497750604e-03, 4.3916291488594718e-03, -1.2320542821102155e-02,
        -7.1712350226568169e-04, 2.0114161532917675e-03, -5.5907196612254920e-03, -9.3539679266380458e-06, -9.2641347200028029e-04, -2.2870225702025344e-16,
        9.6780483124894959e-03, -1.5800101497750504e-03, -4.3916291488594770e-03, -1.2320542821102181e-02, 7.1712350226568690e-04, 2.0114161532917693e-03,
        5.5907196612255145e-03, -9.3539679266376545e-06, -9.2641347200029232e-04, -2.2829444272414866e-16, -9.6780483124892132e-03, 1.5800101497750205e-03,
        -4.3916291488594935e-03, 1.2320542821102315e-02, 7.1712350226569264e-04, -2.0114161532917840e-03, 5.5907196612254538e-03, 9.3539679266373140e-06,
        9.2641347200030078e-04, -5.5723315262335439e-17, -9.6780483124890137e-03, 1.5800101497750283e-03, 4.3916291488596115e-03, 1.2320542821102166e-02,
        -7.1712350226569026e-04, -2.0114161532917680e-03, -5.5907196612255119e-03, 9.3539679266372666e-06, 9.2641347200030143e-04, -9.4706440125226137e-17,
        4.3851076902750950e-03, -1.3897533935028055e-03, 1.9712094061431666e-03, -5.5824175834556997e-03, -6.3058672814617642e-04, 1.7692116884069039e-03,
        -2.5094284625967966e-03, 5.7352089444010442e-05, -4.3418753644622294e-04, -2.1350469860381366e-16, 4.3851076902746621e-03, -1.3897533935027994e-03,
        -1.9712094061433475e-03, -5.5824175834555817e-03, 6.3058672814617761e-04, 1.7692116884068941e-03, 2.5094284625968382e-03, 5.7352089444009608e-05,
        -4.3418753644616927e-04, -9.6273135977081884e-17, -4.3851076902746717e-03, 1.3897533935027498e-03, -1.9712094061431966e-03, 5.5824175834556684e-03,
        6.3058672814615073e-04, -1.7692116884068868e-03, 2.5094284625967606e-03, -5.7352089444010137e-05, 4.3418753644620695e-04, 5.1581442514916895e-17,
        -4.3851076902745945e-03, 1.3897533935027168e-03, 1.9712094061432365e-03, 5.5824175834556606e-03, -6.3058672814616135e-04, -1.7692116884068772e-03,
        -2.5094284625967792e-03, -5.7352089444009873e-05, 4.3418753644620283e-04, 1.1137537455112442e-17, 4.3851076902746916e-03, -1.3897533935027931e-03,
        1.9712094061432530e-03, -5.5824175834556251e-03, -6.3058672814615398e-04, 1.7692116884068991e-03, -2.5094284625968122e-03, 5.7352089444009642e-05,
        -4.3418753644619464e-04, -1.4813473259956328e-16, 4.3851076902747445e-03, -1.3897533935027821e-03, -1.9712094061432421e-03, -5.5824175834556502e-03,
        6.3058672814615831e-04, 1.7692116884068967e-03, 2.5094284625968312e-03, 5.7352089444009920e-05, -4.3418753644621064e-04, -1.5910457447133395e-16,
        -4.3851076902745728e-03, 1.3897533935027500e-03, -1.9712094061432430e-03, 5.5824175834557430e-03, 6.3058672814616276e-04, -1.7692116884069080e-03,
        2.5094284625967874e-03, -5.7352089444010144e-05, 4.3418753644620494e-04, -3.1243208329125511e-17, -4.3851076902743412e-03, 1.3897533935027602e-03,
        1.9712094061433735e-03, 5.5824175834556468e-03, -6.3058672814615756e-04, -1.7692116884068922e-03, -2.5094284625968291e-03, -5.7352089444009412e-05,
        4.3418753644620240e-04, -9.6295615860319939e-17, 7.5073464735969588e-04, -6.0801547049924465e-04, 3.0194257199911358e-04, -9.5571525078467329e-04,
        -2.7511963833464409e-04, 7.7402802696399620e-04, -3.8438497801566524e-04, 1.3522801738399896e-04, -1.0185503033916061e-04, -1.4242012451345626e-16,
        7.5073464735926740e-04, -6.0801547049924812e-04, -3.0194257199929502e-04, -9.5571525078458190e-04, 2.7511963833464257e-04, 7.7402802696399035e-04,
        3.8438497801568947e-04, 1.3522801738399787e-04, -1.0185503033910946e-04, -3.8640659624400708e-19, -7.5073464735945486e-04, 6.0801547049919554e-04,
        -3.0194257199923441e-04, 9.5571525078469509e-04, 2.7511963833461666e-04, -7.7402802696398200e-04, 3.8438497801567277e-04, -1.3522801738399831e-04,
        1.0185503033916149e-04, 4.8896644747051199e-17, -7.5073464735949226e-04, 6.0801547049916149e-04, 3.0194257199923653e-04, 9.5571525078472295e-04,
        -2.7511963833462766e-04, -7.7402802696397040e-04, -3.8438497801568844e-04, -1.3522801738399839e-04, 1.0185503033917997e-04, 2.9084605296716279e-17,
        7.5073464735932194e-04, -6.0801547049923880e-04, 3.0194257199927599e-04, -9.5571525078461019e-04, -2.7511963833462154e-04, 7.7402802696399393e-04,
        -3.8438497801566952e-04, 1.3522801738399825e-04, -1.0185503033914513e-04, -4.9539149096075250e-17, 7.5073464735941865e-04, -6.0801547049922568e-04,
        -3.0194257199925669e-04, -9.5571525078464576e-04, 2.7511963833462230e-04, 7.7402802696398948e-04, 3.8438497801569158e-04, 1.3522801738399858e-04,
        -1.0185503033917348e-04, -7.4179190831097889e-17, -7.5073464735936021e-04, 6.0801547049919218e-04, -3.0194257199924238e-04, 9.5571525078470029e-04,
        2.7511963833462598e-04, -7.7402802696400087e-04, 3.8438497801564849e-04, -1.3522801738399882e-04, 1.0185503033914628e-04, -2.0403375001758317e-17,
        -7.5073464735908807e-04, 6.0801547049920519e-04, 3.0194257199940230e-04, 9.5571525078465302e-04, -2.7511963833461422e-04, -7.7402802696398753e-04,
        -3.8438497801567906e-04, -1.3522801738399798e-04, 1.0185503033914848e-04, -1.1793080162498876e-16,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_PxPx_DyDy(void)
{
  static const double sol[640] = {
        5.8616184533438415e-03, 2.7094777984664276e-03, 3.1565760734900668e-03, -7.4620748753882350e-03, 1.4328949642728596e-03, -3.4492736718159259e-03,
        -4.0184476689694194e-03, -5.2264298589377572e-04, -1.7632215123710571e-04, 7.5979082301536452e-16, 1.4341205890540316e-02, 6.5546632008047770e-03,
        1.7536458636413647e-03, -1.8256929039365127e-02, 7.9551267430269046e-04, -8.3443485748257561e-03, -2.2324613660128874e-03, -1.3363619667711227e-03,
        -1.6506737964354078e-04, 3.6150316530706261e-16, 1.8173569455654341e-02, 8.2864356215895402e-03, 4.8776519957621312e-04, -2.3135681230452718e-02,
        2.2095135802116113e-04, -1.0548964172698679e-02, -6.2094461961555993e-04, -1.7088219069645069e-03, -1.4276292489187270e-04, 5.1918134700881396e-17,
        1.7990573879784012e-02, 8.1929460267052539e-03, -5.5033150535884162e-04, -2.2902720538814019e-02, -2.5005174898972505e-04, -1.0429948176920631e-02,
        7.0059403080549435e-04, -1.6994005263249402e-03, -1.0938846751228856e-04, -1.2584677568311243e-16, 1.4738458963226040e-02, 6.7039542837526826e-03,
        -1.2698766175710315e-03, -1.8762648099115885e-02, -5.7647737886045108e-04, -8.5344020981063303e-03, 1.6166037551302282e-03, -1.3983785661310184e-03,
        -6.4925612047121907e-05, -1.8098485043007565e-16, 9.6780483124897562e-03, 4.3916291488596626e-03, -1.5800101497749428e-03, -1.2320542821101696e-02,
        -7.1712350226569319e-04, -5.5907196612254087e-03, 2.0114161532918473e-03, -9.2641347200033764e-04, -9.3539679266380255e-06, -1.6937631708429027e-16,
        4.3851076902753881e-03, 1.9712094061434021e-03, -1.3897533935026676e-03, -5.5824175834548635e-03, -6.3058672814616796e-04, -2.5094284625967645e-03,
        1.7692116884069950e-03, -4.3418753644630751e-04, 5.7352089444009778e-05, -1.4852136885022670e-16, 7.5073464736049515e-04, 3.0194257199939721e-04,
        -6.0801547049908864e-04, -9.5571525078346354e-04, -2.7511963833462213e-04, -3.8438497801564561e-04, 7.7402802696411872e-04, -1.0185503033937700e-04,
        1.3522801738399966e-04, -1.4286296276697706e-16, 5.8616184533469779e-03, -2.7094777984656123e-03, 3.1565760734897567e-03, -7.4620748753869183e-03,
        -1.4328949642729051e-03, 3.4492736718160564e-03, -4.0184476689695886e-03, -5.2264298589454821e-04, -1.7632215123711094e-04, -6.7813185122706456e-16,
        1.4341205890542550e-02, -6.5546632008040961e-03, 1.7536458636411602e-03, -1.8256929039364284e-02, -7.9551267430271930e-04, 8.3443485748258672e-03,
        -2.2324613660130005e-03, -1.3363619667715958e-03, -1.6506737964354276e-04, -7.8197593865398220e-16, 1.8173569455655954e-02, -8.2864356215889539e-03,
        4.8776519957606008e-04, -2.3135681230452221e-02, -2.2095135802118867e-04, 1.0548964172698731e-02, -6.2094461961564851e-04, -1.7088219069647786e-03,
        -1.4276292489187592e-04, -8.4118925304500103e-16, 1.7990573879785136e-02, -8.1929460267047751e-03, -5.5033150535897292e-04, -2.2902720538813804e-02,
        2.5005174898969182e-04, 1.0429948176920627e-02, 7.0059403080542312e-04, -1.6994005263250738e-03, -1.0938846751228902e-04, -8.4256328087254971e-16,
        1.4738458963226729e-02, -6.7039542837523261e-03, -1.2698766175711523e-03, -1.8762648099115882e-02, 5.7647737886041476e-04, 8.5344020981062939e-03,
        1.6166037551301748e-03, -1.3983785661310505e-03, -6.4925612047123343e-05, -7.5173252314282128e-16, 9.6780483124900302e-03, -4.3916291488594302e-03,
        -1.5800101497750632e-03, -1.2320542821101879e-02, 7.1712350226565871e-04, 5.5907196612253506e-03, 2.0114161532917923e-03, -9.2641347200027953e-04,
        -9.3539679266398297e-06, -5.7142268612241973e-16, 4.3851076902752207e-03, -1.9712094061432820e-03, -1.3897533935028036e-03, -5.5824175834552616e-03,
        6.3058672814613696e-04, 2.5094284625966955e-03, 1.7692116884069258e-03, -4.3418753644613930e-04, 5.7352089444007216e-05, -3.7186439493645241e-16,
        7.5073464735979530e-04, -3.0194257199933616e-04, -6.0801547049926167e-04, -9.5571525078413867e-04, 2.7511963833461807e-04, 3.8438497801561016e-04,
        7.7402802696402645e-04, -1.0185503033901673e-04, 1.3522801738399568e-04, -1.7344565127837030e-16, -5.8616184533437782e-03, -2.7094777984653143e-03,
        -3.1565760734898899e-03, 7.4620748753859260e-03, -1.4328949642730510e-03, 3.4492736718157403e-03, 4.0184476689695582e-03, 5.2264298589424073e-04,
        1.7632215123711454e-04, -2.6639418630864231e-15, -1.4341205890539730e-02, -6.5546632008041325e-03, -1.7536458636412510e-03, 1.8256929039363289e-02,
        -7.9551267430277178e-04, 8.3443485748254734e-03, 2.2324613660130274e-03, 1.3363619667713857e-03, 1.6506737964354411e-04, -1.9292192738161615e-15,
        -1.8173569455653408e-02, -8.2864356215890805e-03, -4.8776519957612161e-04, 2.3135681230451340e-02, -2.2095135802118347e-04, 1.0548964172698379e-02,
        6.2094461961568668e-04, 1.7088219069646359e-03, 1.4276292489187820e-04, -1.3481263701032338e-15, -1.7990573879782777e-02, -8.1929460267048185e-03,
        5.5033150535891741e-04, 2.2902720538813103e-02, 2.5005174898972982e-04, 1.0429948176920324e-02, -7.0059403080535579e-04, 1.6994005263249862e-03,
        1.0938846751228772e-04, -8.9840064440486445e-16, -1.4738458963224539e-02, -6.7039542837522108e-03, 1.2698766175711124e-03, 1.8762648099115476e-02,
        5.7647737886046832e-04, 8.5344020981060979e-03, -1.6166037551300773e-03, 1.3983785661310171e-03, 6.4925612047124915e-05, -6.1187347955024573e-16,
        -9.6780483124879711e-03, -4.3916291488591110e-03, 1.5800101497750261e-03, 1.2320542821101815e-02, 7.1712350226572485e-04, 5.5907196612252595e-03,
        -2.0114161532916912e-03, 9.2641347200032853e-04, 9.3539679266387489e-06, -4.6292186796819506e-16, -4.3851076902733134e-03, -1.9712094061426948e-03,
        1.3897533935027539e-03, 5.5824175834555583e-03, 6.3058672814622597e-04, 2.5094284625967163e-03, -1.7692116884068161e-03, 4.3418753644631916e-04,
        -5.7352089444008754e-05, -3.9198879635261493e-16, -7.5073464735810850e-04, -3.0194257199840505e-04, 6.0801547049918589e-04, 9.5571525078483321e-04,
        2.7511963833473397e-04, 3.8438497801580814e-04, -7.7402802696390958e-04, 1.0185503033949897e-04, -1.3522801738399619e-04, -3.6578027336901036e-16,
        -5.8616184533451425e-03, 2.7094777984647540e-03, -3.1565760734898244e-03, 7.4620748753849216e-03, 1.4328949642729910e-03, -3.4492736718154675e-03,
        4.0184476689696406e-03, 5.2264298589441637e-04, 1.7632215123710636e-04, -1.2133063787372501e-15, -1.4341205890540856e-02, 6.5546632008035184e-03,
        -1.7536458636411736e-03, 1.8256929039362498e-02, 7.9551267430280159e-04, -8.3443485748255965e-03, 2.2324613660130695e-03, 1.3363619667714104e-03,
        1.6506737964353888e-04, -9.9280213755089207e-16, -1.8173569455654268e-02, 8.2864356215885965e-03, -4.8776519957604241e-04, 2.3135681230450685e-02,
        2.2095135802123198e-04, -1.0548964172698592e-02, 6.2094461961572202e-04, 1.7088219069646433e-03, 1.4276292489187519e-04, -7.4002985514707971e-16,
        -1.7990573879783364e-02, 8.1929460267045028e-03, 5.5033150535899775e-04, 2.2902720538812545e-02, -2.5005174898968104e-04, -1.0429948176920546e-02,
        -7.0059403080533129e-04, 1.6994005263250003e-03, 1.0938846751228570e-04, -4.9004278444642846e-16, -1.4738458963224845e-02, 6.7039542837520451e-03,
        1.2698766175711876e-03, 1.8762648099114980e-02, -5.7647737886043146e-04, -8.5344020981063216e-03, -1.6166037551300651e-03, 1.3983785661310227e-03,
        6.4925612047122652e-05, -3.2296552980969570e-16, -9.6780483124880180e-03, 4.3916291488590511e-03, 1.5800101497751003e-03, 1.2320542821101354e-02,
        -7.1712350226569850e-04, -5.5907196612254790e-03, -2.0114161532916825e-03, 9.2641347200030295e-04, 9.3539679266380831e-06, -2.4057760356606480e-16,
        -4.3851076902730957e-03, 1.9712094061426987e-03, 1.3897533935028320e-03, 5.5824175834551185e-03, -6.3058672814621599e-04, -2.5094284625969167e-03,
        -1.7692116884068106e-03, 4.3418753644622441e-04, -5.7352089444008971e-05, -2.0305749738373731e-16, -7.5073464735759318e-04, 3.0194257199835984e-04,
        6.0801547049927826e-04, 9.5571525078440365e-04, -2.7511963833477229e-04, -3.8438497801593683e-04, -7.7402802696390968e-04, 1.0185503033923353e-04,
        -1.3522801738399603e-04, -1.8989785385171275e-16, 5.8616184533443853e-03, 2.7094777984674944e-03, 3.1565760734900872e-03, -7.4620748753883720e-03,
        1.4328949642726720e-03, -3.4492736718159966e-03, -4.0184476689694845e-03, -5.2264298589279788e-04, -1.7632215123711129e-04, -5.4265573850165279e-16,
        1.4341205890540911e-02, 6.5546632008053659e-03, 1.7536458636413798e-03, -1.8256929039365485e-02, 7.9551267430259494e-04, -8.3443485748257423e-03,
        -2.2324613660129611e-03, -1.3363619667708334e-03, -1.6506737964354484e-04, -5.2931364482778452e-16, 1.8173569455654976e-02, 8.2864356215898716e-03,
        4.8776519957622092e-04, -2.3135681230453297e-02, 2.2095135802110061e-04, -1.0548964172698645e-02, -6.2094461961561544e-04, -1.7088219069643828e-03,
        -1.4276292489187880e-04, -4.8820771494737806e-16, 1.7990573879784651e-02, 8.1929460267054152e-03, -5.5033150535884672e-04, -2.2902720538814807e-02,
        -2.5005174898976051e-04, -1.0429948176920573e-02, 7.0059403080542843e-04, -1.6994005263248472e-03, -1.0938846751229052e-04, -4.5741231311481576e-16,
        1.4738458963226618e-02, 6.7039542837527503e-03, -1.2698766175710536e-03, -1.8762648099116912e-02, -5.7647737886047234e-04, -8.5344020981063320e-03,
        1.6166037551301538e-03, -1.3983785661309375e-03, -6.4925612047122964e-05, -3.8795519996012429e-16, 9.6780483124902384e-03, 4.3916291488596557e-03,
        -1.5800101497749808e-03, -1.2320542821102991e-02, -7.1712350226571595e-04, -5.5907196612254781e-03, 2.0114161532917654e-03, -9.2641347200029004e-04,
        -9.3539679266394197e-06, -2.8306978253837860e-16, 4.3851076902757255e-03, 1.9712094061433358e-03, -1.3897533935027134e-03, -5.5824175834564777e-03,
        -6.3058672814617761e-04, -2.5094284625968811e-03, 1.7692116884068887e-03, -4.3418753644631651e-04, 5.7352089444008633e-05, -1.7232291504532142e-16,
        7.5073464736067936e-04, 3.0194257199930120e-04, -6.0801547049913504e-04, -9.5571525078550705e-04, -2.7511963833463054e-04, -3.8438497801588398e-04,
        7.7402802696397853e-04, -1.0185503033950893e-04, 1.3522801738399557e-04, -4.0183722521239453e-17, 5.8616184533470351e-03, -2.7094777984668556e-03,
        3.1565760734897884e-03, -7.4620748753884292e-03, -1.4328949642727011e-03, 3.4492736718147151e-03, -4.0184476689695669e-03, -5.2264298589348862e-04,
        -1.7632215123709768e-04, -1.7670040387231012e-16, 1.4341205890542679e-02, -6.5546632008048776e-03, 1.7536458636411584e-03, -1.8256929039365721e-02,
        -7.9551267430265718e-04, 8.3443485748250831e-03, -2.2324613660129741e-03, -1.3363619667712435e-03, -1.6506737964353709e-04, -1.4159583548387092e-16,
        1.8173569455656079e-02, -8.2864356215895836e-03, 4.8776519957605444e-04, -2.3135681230453547e-02, -2.2095135802115623e-04, 1.0548964172698277e-02,
        -6.2094461961560113e-04, -1.7088219069646528e-03, -1.4276292489187457e-04, -1.4370112997316304e-16, 1.7990573879785251e-02, -8.1929460267053007e-03,
        -5.5033150535897248e-04, -2.2902720538814959e-02, 2.5005174898971768e-04, 1.0429948176920377e-02, 7.0059403080546291e-04, -1.6994005263250254e-03,
        -1.0938846751228595e-04, -1.5154734466435251e-16, 1.4738458963226835e-02, -6.7039542837527511e-03, -1.2698766175711551e-03, -1.8762648099116933e-02,
        5.7647737886044566e-04, 8.5344020981062575e-03, 1.6166037551302000e-03, -1.3983785661310325e-03, -6.4925612047119887e-05, -1.3440270325811173e-16,
        9.6780483124901291e-03, -4.3916291488597260e-03, -1.5800101497750667e-03, -1.2320542821102842e-02, 7.1712350226570121e-04, 5.5907196612255084e-03,
        2.0114161532918192e-03, -9.2641347200029612e-04, -9.3539679266370565e-06, -1.4142510264724841e-16, 4.3851076902753161e-03, -1.9712094061434342e-03,
        -1.3897533935028018e-03, -5.5824175834561316e-03, 6.3058672814617610e-04, 2.5094284625970021e-03, 1.7692116884069483e-03, -4.3418753644620977e-04,
        5.7352089444010801e-05, -1.9072614114096155e-16, 7.5073464735991131e-04, -3.0194257199937617e-04, -6.0801547049925137e-04, -9.5571525078493903e-04,
        2.7511963833464702e-04, 3.8438497801606439e-04, 7.7402802696404770e-04, -1.0185503033922346e-04, 1.3522801738399928e-04, -2.7159400245815774e-16,
        -5.8616184533466023e-03, -2.7094777984677395e-03, -3.1565760734897845e-03, 7.4620748753841835e-03, -1.4328949642725887e-03, 3.4492736718177599e-03,
        4.0184476689697820e-03, 5.2264298589374536e-04, 1.7632215123710281e-04, 6.9624920624588576e-16, -1.4341205890542202e-02, -6.5546632008054067e-03,
        -1.7536458636411407e-03, 1.8256929039362189e-02, -7.9551267430255700e-04, 8.3443485748269514e-03, 2.2324613660131649e-03, 1.3363619667713061e-03,
        1.6506737964353991e-04, 7.9619770587187039e-16, -1.8173569455655531e-02, -8.2864356215898473e-03, -4.8776519957603753e-04, 2.3135681230450622e-02,
        -2.2095135802110795e-04, 1.0548964172699432e-02, 6.2094461961576712e-04, 1.7088219069646905e-03, 1.4276292489187359e-04, 8.8902897800480822e-16,
        -1.7990573879784658e-02, -8.1929460267054447e-03, 5.5033150535898246e-04, 2.2902720538812607e-02, 2.5005174898973703e-04, 1.0429948176921063e-02,
        -7.0059403080529117e-04, 1.6994005263250337e-03, 1.0938846751228540e-04, 8.3212135948128831e-16, -1.4738458963226208e-02, -6.7039542837528223e-03,
        1.2698766175711698e-03, 1.8762648099115160e-02, 5.7647737886047179e-04, 8.5344020981066252e-03, -1.6166037551300441e-03, 1.3983785661310264e-03,
        6.4925612047122625e-05, 6.5153847757252344e-16, -9.6780483124894370e-03, -4.3916291488597390e-03, 1.5800101497750823e-03, 1.2320542821101609e-02,
        7.1712350226571194e-04, 5.5907196612255735e-03, -2.0114161532916600e-03, 9.2641347200029536e-04, 9.3539679266377731e-06, 4.4198077128583854e-16,
        -4.3851076902745763e-03, -1.9712094061434095e-03, 1.3897533935028131e-03, 5.5824175834554360e-03, 6.3058672814618846e-04, 2.5094284625968039e-03,
        -1.7692116884067968e-03, 4.3418753644618364e-04, -5.7352089444011329e-05, 2.2471737429018154e-16, -7.5073464735914738e-04, -3.0194257199931464e-04,
        6.0801547049925528e-04, 9.5571525078477911e-04, 2.7511963833465596e-04, 3.8438497801555562e-04, -7.7402802696389049e-04, 1.0185503033907287e-04,
        -1.3522801738399904e-04, -4.6964796235091913e-18, -5.8616184533469675e-03, 2.7094777984668465e-03, -3.1565760734898287e-03, 7.4620748753871542e-03,
        1.4328949642727234e-03, -3.4492736718168822e-03, 4.0184476689694862e-03, 5.2264298589321757e-04, 1.7632215123710289e-04, 3.1776391741962064e-16,
        -1.4341205890542656e-02, 6.5546632008048533e-03, -1.7536458636411461e-03, 1.8256929039364295e-02, 7.9551267430261793e-04, -8.3443485748262783e-03,
        2.2324613660129611e-03, 1.3363619667710806e-03, 1.6506737964354023e-04, 6.1540592612503074e-16, -1.8173569455655979e-02, 8.2864356215894587e-03,
        -4.8776519957602799e-04, 2.3135681230452138e-02, 2.2095135802114259e-04, -1.0548964172698903e-02, 6.2094461961562705e-04, 1.7088219069645884e-03,
        1.4276292489187310e-04, 8.2893778027973621e-16, -1.7990573879785053e-02, 8.1929460267051672e-03, 5.5033150535899850e-04, 2.2902720538813707e-02,
        -2.5005174898971009e-04, -1.0429948176920639e-02, -7.0059403080539537e-04, 1.6994005263249936e-03, 1.0938846751228491e-04, 8.3447909816312473e-16,
        -1.4738458963226533e-02, 6.7039542837526367e-03, 1.2698766175711900e-03, 1.8762648099115934e-02, -5.7647737886044523e-04, -8.5344020981062904e-03,
        -1.6166037551301263e-03, 1.3983785661310295e-03, 6.4925612047121229e-05, 6.7923021830212132e-16, -9.6780483124897076e-03, 4.3916291488596357e-03,
        1.5800101497751005e-03, 1.2320542821102105e-02, -7.1712350226568917e-04, -5.5907196612253263e-03, -2.0114161532917402e-03, 9.2641347200033612e-04,
        9.3539679266379476e-06, 4.8083758871934514e-16, -4.3851076902747940e-03, 1.9712094061433835e-03, 1.3897533935028239e-03, 5.5824175834556494e-03,
        -6.3058672814616601e-04, -2.5094284625966556e-03, -1.7692116884068809e-03, 4.3418753644626132e-04, -5.7352089444010286e-05, 2.7793202384748856e-16,
        -7.5073464735935424e-04, 3.0194257199937368e-04, 6.0801547049924899e-04, 9.5571525078468316e-04, -2.7511963833462994e-04, -3.8438497801557200e-04,
        -7.7402802696398536e-04, 1.0185503033921120e-04, -1.3522801738399890e-04, 8.1559707371039833e-17,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_DxDx_NyDy(void)
{
  static const double sol[640] = {
        -7.4185325279775777e-04, -3.9802597843886317e-04, 2.0456784067311337e-05, 9.4440888007033875e-04, 1.0964684083884220e-05, 5.0670299970896452e-04,
        -2.6042304806390264e-05, 2.3457263346557451e-05, 5.0764812234692848e-06, 4.6249790389122405e-16, -6.0391518322723753e-04, -3.2401990289661228e-04,
        5.8133025438359082e-05, 7.6880819717104528e-04, 3.1157234701922163e-05, 4.1249030379146855e-04, -7.4005667889377098e-05, 1.9094342915952928e-05,
        4.2641058539285365e-06, 2.2220549807193982e-16, -3.5046776690493651e-04, -1.8804143086081278e-04, 8.6097341540881974e-05, 4.4615949312709391e-04,
        4.6138941825176955e-05, 2.3938426697811351e-04, -1.0960536143111664e-04, 1.1077684120019632e-05, 2.6393373958545555e-06, -6.0886009644883065e-17,
        -2.6369173463173775e-05, -1.4156115092826504e-05, 9.7874810518121645e-05, 3.3569013123064332e-05, 5.2437299075732224e-05, 1.8021301046369540e-05,
        -1.2459855077791895e-04, 8.2737327434932787e-07, 2.0225634437670926e-07, 6.6278316245383176e-17, 3.0109676195020615e-04, 1.6152938865401258e-04,
        8.6992163559215918e-05, -3.8330822797399874e-04, 4.6572459248672606e-05, -2.0563337622609749e-04, -1.1074450567150355e-04, -9.5344390692561088e-06,
        -3.0476196625405792e-06, 2.2276701672177206e-17, 5.4219548657417969e-04, 2.9091105641263198e-04, 4.6964511263907207e-05, -6.9023655328536414e-04,
        2.5099762187585532e-05, -3.7034141718852477e-04, -5.9787702377023312e-05, -1.7138457671756282e-05, -7.1078530314263283e-06, -8.0749235629294875e-17,
        5.8489869202431376e-04, 3.1367048785931479e-04, -2.8617271620076177e-05, -7.4459944282300513e-04, -1.5601216384924636e-05, -3.9931508426118849e-04,
        3.6430932046668254e-05, -1.8606531928361270e-05, -1.1993596362267388e-05, 1.4335740025532410e-16, 2.9366410440273647e-04, 1.5833417110391160e-04,
        -1.4679523967225448e-04, -3.7384615745770508e-04, -7.8278512648658480e-05, -2.0156573641113469e-04, 1.8687621490631206e-04, -8.6854584349767026e-06,
        -1.7623476288397682e-05, 1.0848651979579840e-16, -1.8111297630617701e-03, -2.2122899042749104e-04, 4.9921555573584552e-05, 2.3056406705032101e-03,
        6.1001783803918497e-06, 2.8163335848549242e-04, -6.3552138125501126e-05, 2.1979419790257315e-05, 1.2783176280283245e-05, 2.7857679634059921e-16,
        -1.4743782859745810e-03, -1.8009465973339172e-04, 1.4186131475810724e-04, 1.8769425632445258e-03, 1.7335155668372967e-05, 2.2926770929982310e-04,
        -1.8059513103274984e-04, 1.7891815716313254e-05, 1.0737709998699541e-05, 1.9230809006369143e-16, -8.5563118982650362e-04, -1.0451452877111586e-04,
        2.1009107982633018e-04, 1.0892527473457480e-03, 2.5673892807397301e-05, 1.3305117783823689e-04, -2.6745435254661471e-04, 1.0380762902201684e-05,
        6.6467490018292703e-06, 1.0350256818363629e-16, -6.4401076442146066e-05, -7.8652359084521495e-06, 2.3880385491964199e-04, 8.1985147666923674e-05,
        2.9185828660963689e-05, 1.0012760081344271e-05, -3.0400686433692739e-04, 7.7809221149302604e-07, 5.1040642097443571e-07, 4.0017417234473824e-17,
        7.3504616021202588e-04, 8.9786913323456805e-05, 2.1219385883249455e-04, -9.3574317878674154e-04, 2.5939534318854238e-05, -1.1430233396837355e-04,
        -2.7013127437547352e-04, -8.9337392931330060e-06, -7.6719397877145017e-06, 7.5198342416082037e-18, 1.3236741350855302e-03, 1.6169038423836875e-04,
        1.1444951102739500e-04, -1.6850901479238580e-03, 1.4004905074148330e-05, -2.0583833004827264e-04, -1.4569880785232482e-04, -1.6042765456802835e-05,
        -1.7897506519468700e-05, -3.8723351581988297e-18, 1.4277325347643744e-03, 1.7441241606392551e-04, -7.0228009414156255e-05, -1.8175606551730900e-03,
        -8.5314378017897350e-06, -2.2203398570291407e-04, 8.9403066536777049e-05, -1.7465893456537019e-05, -3.0180009723344792e-05, -1.0928219871388014e-17,
        7.1871323474241155e-04, 8.7710214994633017e-05, -3.5755576442373147e-04, -9.1495071101386225e-04, -4.3774882151247583e-05, -1.1165861388560425e-04,
        4.5518279763373005e-04, -8.1880463208147463e-06, -4.4456825652039233e-05, 7.1223580357852086e-18, -2.2944506126900003e-03, -6.1638761923244069e-05,
        6.3246271865104631e-05, 2.9209274547706780e-03, 1.6996066367278695e-06, 7.8468610735826896e-05, -8.0515035224969260e-05, 1.9018655925888491e-05,
        1.6147340377872063e-05, 2.0553230239201092e-16, -1.8678324615435337e-03, -5.0177806705329094e-05, 1.7972636255985107e-04, 2.3778254749349001e-03,
        4.8298278531543210e-06, 6.3878356071487048e-05, -2.2879885225845367e-04, 1.5481586715521866e-05, 1.3563696440515726e-05, 1.6862058727884781e-16,
        -1.0839644414274879e-03, -2.9119379734342048e-05, 2.6616909306346494e-04, 1.3799301146204916e-03, 7.1531646979130580e-06, 3.7070135770886574e-05,
        -3.3884390766164969e-04, 8.9821168968253373e-06, 8.3963955227197131e-06, 1.1080275815901326e-16, -8.1583353672125855e-05, -2.1908826904612249e-06,
        3.0254949112894263e-04, 1.0385887422230145e-04, 8.1313422416431917e-06, 2.7890813449652836e-06, -3.8515761035684532e-04, 6.7280421817372825e-07,
        6.4546625765658320e-07, 5.5833352067559234e-17, 9.3120558827591447e-04, 2.5018632896846171e-05, 2.6884245930155251e-04, -1.1854619810895742e-03,
        7.2279537638850403e-06, -3.1849720933209709e-05, -3.4224720987189211e-04, -7.7309086961042199e-06, -9.6892163186684390e-06, 1.2233071025946061e-17,
        1.6769291755919060e-03, 4.5047880146268317e-05, 1.4502400361684611e-04, -2.1347979518942318e-03, 3.8980585320111308e-06, -5.7347754260091483e-05,
        -1.8462136052193066e-04, -1.3884451253037142e-05, -2.2607359364033955e-05, -2.0946695208114841e-17, 1.8087408305902113e-03, 4.8609262873996871e-05,
        -8.8933756155094482e-05, -2.3025994638613619e-03, -2.3689535246092704e-06, -6.1881537000419012e-05, 1.1321623074936586e-04, -1.5106133674163060e-05,
        -3.8108004663823115e-05, -4.9996761330612133e-17, 9.1041892369570464e-04, 2.4407192087710259e-05, -4.5307301340336866e-04, -1.1589997251880580e-03,
        -1.2201146101856527e-05, -3.1071332312954140e-05, 5.7678007822296613e-04, -7.0756666654389634e-06, -5.6202873171474166e-05, -1.3729881301442140e-16,
        -2.2712814698433652e-03, 6.9283273198997747e-05, 6.2607182589646940e-05, 2.8914322086887731e-03, -1.9084463650281187e-06, -8.8200379526165773e-05,
        -7.9701449000710122e-05, 1.4578513220970450e-05, 1.5993584902369550e-05, 1.6206959329244378e-16, -1.8489709303060869e-03, 5.6401394648714923e-05,
        1.7791020436132239e-04, 2.3538139908238103e-03, -5.4230258469318264e-06, -7.1801232593806587e-05, -2.2648681018836513e-04, 1.1867224005866736e-05,
        1.3434536925262127e-05, 1.3884044675756950e-16, -1.0730176634537302e-03, 3.2732351109901716e-05, 2.6347918201770733e-04, 1.3659944281656225e-03,
        -8.0305360856003998e-06, -4.1669592924460280e-05, -3.3541954324910466e-04, 6.8851645536146903e-06, 8.3164350440197039e-06, 9.8989004851509373e-17,
        -8.0757424528281682e-05, 2.4649783637425219e-06, 2.9949136733266980e-04, 1.0280743336824143e-04, -9.1266848261765632e-06, -3.1380161064322520e-06,
        -3.8126449637687740e-04, 5.1575107170413136e-07, 6.3927833037311136e-07, 5.4320200009267459e-17, 9.2180607376450324e-04, -2.8115105526507709e-05,
        2.6612401553146117e-04, -1.1734960229442387e-03, -8.1052313173914392e-06, 3.5791654512843503e-05, -3.3878652215938284e-04, -5.9258354008351692e-06,
        -9.5969568198567602e-06, 1.2354355107121336e-17, 1.6599996084286698e-03, -5.0637183330130939e-05, 1.4355365040744227e-04, -2.1132459353674824e-03,
        -4.3684028831008672e-06, 6.4463160899271821e-05, -1.8274954204225148e-04, -1.0642774281309240e-05, -2.2392186132941561e-05, -2.4439344978395069e-17,
        1.7904898493183452e-03, -5.4594815714699535e-05, -8.8038739116137082e-05, -2.2793652342907554e-03, 2.7188239063177206e-06, 6.9501385310809908e-05,
        1.1207683823986927e-04, -1.1580436012146539e-05, -3.7747530559387761e-05, -5.5342821763728846e-17, 9.0124227629169657e-04, -2.7574015762423465e-05,
        -4.4848460248174822e-04, -1.1473174857899373e-03, 1.3619750416073549e-05, 3.5102825588544167e-05, 5.7093884748957700e-04, -5.4247562557651676e-06,
        -5.5653119508561458e-05, -9.2909246370792702e-17, -1.8607939734754306e-03, 1.6006863502989764e-04, 5.1292513034897540e-05, 2.3688651979412874e-03,
        -4.4105638650402665e-06, -2.0377377840281942e-04, -6.5297421839984664e-05, 8.6581482607711247e-06, 1.3099218422382047e-05, 1.2130732581462570e-16,
        -1.5148068254142147e-03, 1.3030661548299124e-04, 1.4575747047953651e-04, 1.9284097119175174e-03, -1.2533216551743824e-05, -1.6588553643187131e-04,
        -1.8555509319178660e-04, 7.0478959388719097e-06, 1.1003322845519252e-05, 1.0671497867529702e-16, -8.7908999816598422e-04, 7.5622040169426828e-05,
        2.1586222088963472e-04, 1.1191167491928862e-03, -1.8560279471199331e-05, -9.6269883559477111e-05, -2.7480124608354376e-04, 4.0890266760258649e-06,
        6.8115232787057308e-06, 7.9717774256069590e-17, -6.6158684055542149e-05, 5.6930861126356883e-06, 2.4536651900581864e-04, 8.4222652499154474e-05,
        -2.1095228673489870e-05, -7.2475264609192327e-06, -3.1236139835906653e-04, 3.0620542597938747e-07, 5.2383308233409013e-07, 4.5367438840226148e-17,
        7.5521454677788689e-04, -6.4959947881067446e-05, 2.1803021136887151e-04, -9.6141834203176198e-04, -1.8739033755341927e-05, 8.2696613375289981e-05,
        -2.7756118472743952e-04, -3.5194659963031186e-06, -7.8598788404717967e-06, 1.0409110540003064e-17, 1.3599971624203907e-03, -1.1698967736129318e-04,
        1.1761250922651231e-04, -1.7313308153829561e-03, -1.0103627331035097e-05, 1.4893254125386500e-04, -1.4972543114425156e-04, -6.3208097101481016e-06,
        -1.8339153449490702e-05, -1.8592739562088397e-17, 1.4669018971550508e-03, -1.2615623684791310e-04, -7.2126947445431660e-05, -1.8674248210697011e-03,
        6.2473383435715822e-06, 1.6060193832964713e-04, 9.1820490646870895e-05, -6.8766677416233056e-06, -3.0915439591944671e-05, -3.4169910043710598e-17,
        7.3835886561778822e-04, -6.3621467426451799e-05, -3.6743879842920073e-04, -9.3996038534386532e-04, 3.1529332913731504e-05, 8.0992674189995619e-05,
        4.6776429544571370e-04, -3.2207980489607478e-06, -4.5587253648864800e-05, -1.3033565763243572e-17, -1.2218883130104699e-03, 1.9925006776766980e-04,
        3.3680823229324827e-05, 1.5555127229134420e-03, -5.4906542718873293e-06, -2.5365331033436574e-04, -4.2877035890732717e-05, 1.2579210421270085e-06,
        8.6095141065868681e-06, 8.5843190723295419e-17, -9.9469578895192406e-04, 1.6220276314314548e-04, 9.5710452873525549e-05, 1.2662875474526807e-03,
        -1.5602511542291027e-05, -2.0649060889963317e-04, -1.2184323687799138e-04, 1.0239596690626961e-06, 7.2319800498269018e-06, 7.6128518192850284e-17,
        -5.7725205425157491e-04, 9.4132221087062096e-05, 1.4174402118255613e-04, 7.3486496691657513e-04, -2.3105816616276255e-05, -1.1983408465235449e-04,
        -1.8044581161690831e-04, 5.9407987857311818e-07, 4.4769112358823013e-06, 5.8602639297672241e-17, -4.3439872336690766e-05, 7.0856641229987313e-06,
        1.6111739776339482e-04, 5.5300695965309908e-05, -2.6262233432286444e-05, -9.0203339996338009e-06, -2.0510889533447221e-04, 4.4421321296519763e-08,
        3.4429239733250568e-07, 3.4926715594934423e-17, 4.9591678847464517e-04, -8.0861991798899719e-05, 1.4316667347179717e-04, -6.3132191851335946e-04,
        -2.3329868400283431e-05, 1.0294055168874697e-04, -1.8225690491621555e-04, -5.1117197886643198e-07, -5.1659110589382652e-06, 1.0102921131673717e-17,
        8.9304807164371051e-04, -1.4562879398958924e-04, 7.7225405099948114e-05, -1.1368859353378222e-03, -1.2582857886102611e-05, 1.8539128287040646e-04,
        -9.8311031283336125e-05, -9.1822109815642119e-07, -1.2053450101606588e-05, -8.9611256229312162e-18, 9.6325796444312359e-04, -1.5703540377679862e-04,
        -4.7364249466624545e-05, -1.2262659386990057e-03, 7.7720809383173407e-06, 1.9991235362653694e-04, 6.0296585106922121e-05, -9.9955532671637544e-07,
        -2.0321917006824889e-05, -1.4644761409477663e-17, 4.8485609369285158e-04, -7.9184942178627942e-05, -2.4126875678009328e-04, -6.1724121140272088e-04,
        3.9263369648519171e-05, 1.0080560040594159e-04, 3.0714478305167290e-04, -4.6828117363387576e-07, -2.9948261566553504e-05, -1.2605868310791414e-17,
        -5.5318357068331902e-04, 1.7536157032520592e-04, 1.5250285570937922e-05, 7.0422482410403353e-04, -4.8299508614752260e-06, -2.2324229706304245e-04,
        -1.9414223854234676e-05, -7.6230039760333720e-06, 3.8612474615733436e-06, 5.1737179513716968e-17, -4.5032633391021244e-04, 1.4275601491743441e-04,
        4.3336922314421925e-05, 5.7328344530477055e-04, -1.3724632897095447e-05, -1.8173412014178895e-04, -5.5169636466947756e-05, -6.2054078588158429e-06,
        3.2435488692399899e-06, 4.8028491712169638e-17, -2.6133624469453566e-04, 8.2846618815280416e-05, 6.4181688197263934e-05, 3.3269149827546476e-04,
        -2.0323581324053096e-05, -1.0546706130613400e-04, -8.1705857651543159e-05, -3.6005452741340406e-06, 2.0081691410636277e-06, 3.8520388310352902e-17,
        -1.9661496835333231e-05, 6.2364545034745673e-06, 7.2956872225951117e-05, 2.5029872332176844e-05, -2.3096613004447442e-05, -7.9392561682747264e-06,
        -9.2877018106376599e-05, -2.7032276650755302e-07, 1.5500833779823407e-07, 2.3954220237922715e-17, 2.2452284971681333e-04, -7.1167813182502039e-05,
        6.4833668634053107e-05, -2.8582657318236793e-04, -2.0511575929727710e-05, 9.0599474345224277e-05, -8.2535854840078805e-05, 3.0981066058928703e-06,
        -2.3154572454531680e-06, 8.4437989845047043e-18, 4.0432786569248309e-04, -1.2816820392756137e-04, 3.4987563703451466e-05, -5.1472555438688986e-04,
        -1.1047540869100629e-05, 1.6316325294175941e-04, -4.4540568810560594e-05, 5.5618567638296239e-06, -5.4055108406485347e-06, -3.3020213440190191e-18,
        4.3610842151560496e-04, -1.3819686293535001e-04, -2.1413925090722030e-05, -5.5518347382006619e-04, 6.8888362637683960e-06, 1.7593013720955338e-04,
        2.7260783634202321e-05, 6.0629566842151169e-06, -9.1041122183861232e-06, -7.2371742466117540e-18, 2.1943348004293205e-04, -6.9776478809210302e-05,
        -1.0931001556628712e-04, -2.7934760190890623e-04, 3.4451757069420710e-05, 8.8828250006272910e-05, 1.3915602444573250e-04, 2.8472996610849330e-06,
        -1.3462461610689810e-05, -1.4223434379638855e-18, -9.3040091887961748e-05, 7.6930184467915048e-05, 2.5570471988835521e-06, 1.1844376047444965e-04,
        -2.1248537374219816e-06, -9.7935203603906131e-05, -3.2552234182110448e-06, -1.7981070634123119e-05, 7.9801219628408649e-07, 1.3255963992479744e-17,
        -7.5742210399364090e-05, 6.2625581143333398e-05, 7.2652034204423278e-06, 9.6422865071444226e-05, -6.0388313855885586e-06, -7.9724871095860508e-05,
        -9.2488947105183267e-06, -1.4636620319132153e-05, 6.7044703542109376e-07, 1.5454990145859717e-17, -4.3958996718337654e-05, 3.6342108184726977e-05,
        1.0755488950520903e-05, 5.5961562078806149e-05, -8.9457858332929588e-06, -4.6264958144601579e-05, -1.3692167872881542e-05, -8.4914149761173461e-06,
        4.1533005655934768e-07, 2.3331565406350569e-17, -3.3149552045454746e-06, 2.7320505889872851e-06, 1.2216240033802679e-05, 4.2200706411281198e-06,
        -1.0173710041968676e-05, -3.4780097375416221e-06, -1.5551762461788650e-05, -6.3374425050565725e-07, 3.2586873501894366e-08, 1.4287004275686782e-17,
        3.7749886887661497e-05, -3.1230076485256521e-05, 1.0834213338389489e-05, -4.8057116773325226e-05, -9.0542379286305008e-06, 3.9757137205703429e-05,
        -1.3792387169260664e-05, 7.3084505808951379e-06, -4.7736440058582091e-07, -5.1426902537239170e-19, 6.8002926475318857e-05, -5.6225462188814054e-05,
        5.8069646962412217e-06, -8.6570446906972983e-05, -4.9002968245301933e-06, 7.1577263531492970e-05, -7.3924984553271851e-06, 1.3140222484471872e-05,
        -1.1167339315051476e-06, -1.0014636178582984e-18, 7.3270287126041972e-05, -6.0708847522397337e-05, -3.7405038611585289e-06, -9.3276007817187916e-05,
        2.8745282594752142e-06, 7.7284792487988441e-05, 4.7618111116923385e-06, 1.4257402535386491e-05, -1.8716955469369723e-06, -1.7259023988470547e-18,
        3.7604074979957093e-05, -3.0294674009570444e-05, -1.8072883786921887e-05, -4.7871492379369695e-05, 1.5390330853863949e-05, 3.8566332418962707e-05,
        2.3007504344675715e-05, 6.6491000063331675e-06, -2.8178579922985129e-06, -7.5910951338837699e-19,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_DxDx_DyNy(void)
{
  static const double sol[640] = {
        2.9366410440246786e-04, 1.5833417110396627e-04, 1.4679523967276470e-04, -3.7384615745761769e-04, 7.8278512648382849e-05, -2.0156573641117145e-04,
        -1.8687621490639760e-04, -8.6854584348496205e-06, -1.7623476288500814e-05, 5.2300407698711697e-17, 5.8489869202463327e-04, 3.1367048785918631e-04,
        2.8617271619829673e-05, -7.4459944282316928e-04, 1.5601216385170638e-05, -3.9931508426113612e-04, -3.6430932046669887e-05, -1.8606531928366969e-05,
        -1.1993596362429246e-05, -2.0892545508593297e-17, 5.4219548657418501e-04, 2.9091105641273601e-04, -4.6964511263591833e-05, -6.9023655328552677e-04,
        -2.5099762187887353e-05, -3.7034141718842486e-04, 5.9787702377095255e-05, -1.7138457671773087e-05, -7.1078530313942427e-06, 3.9115154576600269e-17,
        3.0109676195039589e-04, 1.6152938865406988e-04, -8.6992163559288790e-05, -3.8330822797356349e-04, -4.6572459248522587e-05, -2.0563337622632095e-04,
        1.1074450567142924e-04, -9.5344390692464763e-06, -3.0476196624031540e-06, 5.6810779361459721e-17, -2.6369173462459848e-05, -1.4156115092992755e-05,
        -9.7874810517917558e-05, 3.3569013122980706e-05, -5.2437299075891127e-05, 1.8021301046486428e-05, 1.2459855077798471e-04, 8.2737327436896982e-07,
        2.0225634438136138e-07, -1.3373254239006441e-16, -3.5046776690432491e-04, -1.8804143086091797e-04, -8.6097341541172933e-05, 4.4615949312728348e-04,
        -4.6138941824949177e-05, 2.3938426697800533e-04, 1.0960536143108869e-04, 1.1077684120011821e-05, 2.6393373958376242e-06, -1.6032627550668000e-16,
        -6.0391518322698437e-04, -3.2401990289661098e-04, -5.8133025438105893e-05, 7.6880819717087235e-04, -3.1157234702228416e-05, 4.1249030379158797e-04,
        7.4005667889267377e-05, 1.9094342915948621e-05, 4.2641058540425827e-06, 2.0275386120011121e-16, -7.4185325279781252e-04, -3.9802597843876847e-04,
        -2.0456784067544166e-05, 9.4440888007044500e-04, -1.0964684083691834e-05, 5.0670299970890988e-04, 2.6042304806449675e-05, 2.3457263346559995e-05,
        5.0764812237368312e-06, 2.3746653068095925e-16, 7.1871323474278386e-04, 8.7710214994792733e-05, 3.5755576442357139e-04, -9.1495071101398823e-04,
        4.3774882151079233e-05, -1.1165861388563981e-04, -4.5518279763369802e-04, -8.1880463208082919e-06, -4.4456825652008327e-05, 3.4242442663156875e-18,
        1.4277325347644839e-03, 1.7441241606391770e-04, 7.0228009414156716e-05, -1.8175606551731503e-03, 8.5314378018246768e-06, -2.2203398570289667e-04,
        -8.9403066536790140e-05, -1.7465893456554210e-05, -3.0180009723319835e-05, 8.4266422106792768e-18, 1.3236741350856989e-03, 1.6169038423838374e-04,
        -1.1444951102738993e-04, -1.6850901479238550e-03, -1.4004905074141413e-05, -2.0583833004829191e-04, 1.4569880785236580e-04, -1.6042765456800721e-05,
        -1.7897506519466129e-05, 7.4908691137896685e-18, 7.3504616021230983e-04, 8.9786913323435284e-05, -2.1219385883246752e-04, -9.3574317878671607e-04,
        -2.5939534318855813e-05, -1.1430233396838016e-04, 2.7013127437548387e-04, -8.9337392931378460e-06, -7.6719397877386337e-06, -3.8697968593087614e-18,
        -6.4401076441772531e-05, -7.8652359084996697e-06, -2.3880385491960462e-04, 8.1985147666986571e-05, -2.9185828660990960e-05, 1.0012760081324107e-05,
        3.0400686433691541e-04, 7.7809221149930562e-07, 5.1040642096020408e-07, -1.3552017760322912e-17, -8.5563118982616470e-04, -1.0451452877114972e-04,
        -2.1009107982637403e-04, 1.0892527473457853e-03, -2.5673892807371934e-05, 1.3305117783822941e-04, 2.6745435254661032e-04, 1.0380762902208026e-05,
        6.6467490018253722e-06, 2.7973538917408030e-17, -1.4743782859743448e-03, -1.8009465973339782e-04, -1.4186131475814963e-04, 1.8769425632445393e-03,
        -1.7335155668355575e-05, 2.2926770929982388e-04, 1.8059513103275228e-04, 1.7891815716312892e-05, 1.0737709998674950e-05, 1.0215351539542659e-16,
        -1.8111297630616155e-03, -2.2122899042747670e-04, -4.9921555573624397e-05, 2.3056406705032179e-03, -6.1001783803642330e-06, 2.8163335848549524e-04,
        6.3552138125521157e-05, 2.1979419790251000e-05, 1.2783176280231583e-05, 1.6952945329051220e-16, 9.1041892369545115e-04, 2.4407192087632311e-05,
        4.5307301340360138e-04, -1.1589997251880938e-03, 1.2201146101923685e-05, -3.1071332313057810e-05, -5.7678007822300462e-04, -7.0756666650892133e-06,
        -5.6202873171458269e-05, 3.4656927675756581e-18, 1.8087408305903975e-03, 4.8609262873996356e-05, 8.8933756155082204e-05, -2.3025994638614634e-03,
        2.3689535246087084e-06, -6.1881537000423172e-05, -1.1321623074934100e-04, -1.5106133674220888e-05, -3.8108004663833130e-05, -2.0277669099632116e-18,
        1.6769291755920760e-03, 4.5047880146264644e-05, -1.4502400361682781e-04, -2.1347979518942664e-03, -3.8980585320159521e-06, -5.7347754260095698e-05,
        1.8462136052194050e-04, -1.3884451253027081e-05, -2.2607359364032397e-05, -3.6018740075868448e-18, 9.3120558827613370e-04, 2.5018632896834594e-05,
        -2.6884245930153939e-04, -1.1854619810895682e-03, -7.2279537638890951e-06, -3.1849720933219562e-05, 3.4224720987189965e-04, -7.7309086961055159e-06,
        -9.6892163186640667e-06, 2.9794417583722753e-18, -8.1583353671867056e-05, -2.1908826904851967e-06, -3.0254949112893797e-04, 1.0385887422231885e-04,
        -8.1313422416429613e-06, 2.7890813449587394e-06, 3.8515761035684765e-04, 6.7280421817594080e-07, 6.4546625765690147e-07, 2.4525852988621017e-17,
        -1.0839644414272339e-03, -2.9119379734364328e-05, -2.6616909306347378e-04, 1.3799301146205072e-03, -7.1531646979110311e-06, 3.7070135770882644e-05,
        3.3884390766164801e-04, 8.9821168968263622e-06, 8.3963955227189118e-06, 6.1062881766822862e-17, -1.8678324615433244e-03, -5.0177806705337605e-05,
        -1.7972636255985985e-04, 2.3778254749349070e-03, -4.8298278531550706e-06, 6.3878356071485530e-05, 2.2879885225844988e-04, 1.5481586715522154e-05,
        1.3563696440521318e-05, 1.0170969510065392e-16, -2.2944506126898303e-03, -6.1638761923241793e-05, -6.3246271865108493e-05, 2.9209274547706884e-03,
        -1.6996066367282102e-06, 7.8468610735822411e-05, 8.0515035224968325e-05, 1.9018655925887671e-05, 1.6147340377885233e-05, 1.3017529747857753e-16,
        9.0124227629200860e-04, -2.7574015762490868e-05, 4.4848460248164451e-04, -1.1473174857901441e-03, -1.3619750415973458e-05, 3.5102825588677369e-05,
        -5.7093884748951151e-04, -5.4247562557807852e-06, -5.5653119508569121e-05, 7.1998202120845336e-18, 1.7904898493184785e-03, -5.4594815714664041e-05,
        8.8038739116150133e-05, -2.2793652342908430e-03, -2.7188239063540939e-06, 6.9501385310805707e-05, -1.1207683823985288e-04, -1.1580436012153722e-05,
        -3.7747530559385986e-05, -1.5499910534328783e-17, 1.6599996084288332e-03, -5.0637183330141720e-05, -1.4355365040743972e-04, -2.1132459353675240e-03,
        4.3684028831044450e-06, 6.4463160899274558e-05, 1.8274954204226132e-04, -1.0642774281306924e-05, -2.2392186132941347e-05, -9.7248139365350449e-18,
        9.2180607376468311e-04, -2.8115105526518080e-05, -2.6612401553145607e-04, -1.1734960229442541e-03, 8.1052313173892082e-06, 3.5791654512840400e-05,
        3.3878652215938977e-04, -5.9258354008358425e-06, -9.5969568198578901e-06, 6.9752705364097240e-18, -8.0757424528087745e-05, 2.4649783637273693e-06,
        -2.9949136733266676e-04, 1.0280743336824149e-04, 9.1266848261757754e-06, -3.1380161064356613e-06, 3.8126449637687908e-04, 5.1575107170505907e-07,
        6.3927833037202208e-07, 3.5192119948421573e-17, -1.0730176634535381e-03, 3.2732351109887696e-05, -2.6347918201771086e-04, 1.3659944281656262e-03,
        8.0305360856018143e-06, -4.1669592924463010e-05, 3.3541954324910477e-04, 6.8851645536152629e-06, 8.3164350440190415e-06, 6.7623191256853138e-17,
        -1.8489709303059082e-03, 5.6401394648705165e-05, -1.7791020436132767e-04, 2.3538139908238099e-03, 5.4230258469341092e-06, -7.1801232593809148e-05,
        2.2648681018836388e-04, 1.1867224005866285e-05, 1.3434536925260396e-05, 9.6858086594756248e-17, -2.2712814698432004e-03, 6.9283273198992435e-05,
        -6.2607182589649867e-05, 2.8914322086887683e-03, 1.9084463650295836e-06, -8.8200379526170002e-05, 7.9701449000709159e-05, 1.4578513220969622e-05,
        1.5993584902367352e-05, 1.1458621488120533e-16, 7.3835886561804734e-04, -6.3621467426478321e-05, 3.6743879842914538e-04, -9.3996038534395249e-04,
        -3.1529332913732161e-05, 8.0992674190009836e-05, -4.6776429544570817e-04, -3.2207980490257788e-06, -4.5587253648868073e-05, -5.2983440882930716e-17,
        1.4669018971551885e-03, -1.2615623684792459e-04, 7.2126947445427228e-05, -1.8674248210697438e-03, -6.2473383435674419e-06, 1.6060193832965298e-04,
        -9.1820490646871234e-05, -6.8766677416135300e-06, -3.0915439591939562e-05, -2.2073630823871677e-17, 1.3599971624205301e-03, -1.1698967736130279e-04,
        -1.1761250922651303e-04, -1.7313308153829919e-03, 1.0103627331035624e-05, 1.4893254125386966e-04, 1.4972543114425736e-04, -6.3208097101507841e-06,
        -1.8339153449490641e-05, -1.3511233360745342e-17, 7.5521454677802871e-04, -6.4959947881077407e-05, -2.1803021136886910e-04, -9.6141834203178074e-04,
        1.8739033755341243e-05, 8.2696613375290469e-05, 2.7756118472744331e-04, -3.5194659963027014e-06, -7.8598788404715155e-06, 6.3647442170826156e-18,
        -6.6158684055394061e-05, 5.6930861126234088e-06, -2.4536651900581766e-04, 8.4222652499145448e-05, 2.1095228673489402e-05, -7.2475264609209123e-06,
        3.1236139835906821e-04, 3.0620542597954841e-07, 5.2383308233434244e-07, 3.2945633218545722e-17, -8.7908999816583612e-04, 7.5622040169414726e-05,
        -2.1586222088963581e-04, 1.1191167491928801e-03, 1.8560279471199578e-05, -9.6269883559480283e-05, 2.7480124608354371e-04, 4.0890266760259386e-06,
        6.8115232787059417e-06, 6.0328404987001239e-17, -1.5148068254140711e-03, 1.3030661548298094e-04, -1.4575747047953794e-04, 1.9284097119175055e-03,
        1.2533216551744164e-05, -1.6588553643187519e-04, 1.8555509319178416e-04, 7.0478959388720808e-06, 1.1003322845519606e-05, 8.2874361330714665e-17,
        -1.8607939734752923e-03, 1.6006863502988851e-04, -5.1292513034897953e-05, 2.3688651979412710e-03, 4.4105638650403004e-06, -2.0377377840282172e-04,
        6.5297421839983865e-05, 8.6581482607713805e-06, 1.3099218422382941e-05, 9.7474000767547502e-17, 4.8485609369288524e-04, -7.9184942178639150e-05,
        2.4126875678013267e-04, -6.1724121140272294e-04, -3.9263369648531429e-05, 1.0080560040594443e-04, -3.0714478305168537e-04, -4.6828117362726122e-07,
        -2.9948261566555140e-05, -1.5921508748764577e-17, 9.6325796444322377e-04, -1.5703540377681838e-04, 4.7364249466624342e-05, -1.2262659386990328e-03,
        -7.7720809383103560e-06, 1.9991235362654542e-04, -6.0296585106920495e-05, -9.9955532671419899e-07, -2.0321917006825774e-05, -2.2758358789180729e-17,
        8.9304807164381145e-04, -1.4562879398959960e-04, -7.7225405099946433e-05, -1.1368859353378421e-03, 1.2582857886101947e-05, 1.8539128287041017e-04,
        9.8311031283338307e-05, -9.1822109815711544e-07, -1.2053450101606825e-05, -1.3072496749888570e-17, 4.9591678847475055e-04, -8.0861991798911103e-05,
        -1.4316667347179625e-04, -6.3132191851337312e-04, 2.3329868400283174e-05, 1.0294055168874902e-04, 1.8225690491621642e-04, -5.1117197886615828e-07,
        -5.1659110589384939e-06, 2.0803388797406004e-18, -4.3439872336584223e-05, 7.0856641229869017e-06, -1.6111739776339564e-04, 5.5300695965298171e-05,
        2.6262233432285939e-05, -9.0203339996339110e-06, 2.0510889533447205e-04, 4.4421321296593752e-08, 3.4429239733182826e-07, 2.2562023417140470e-17,
        -5.7725205425146931e-04, 9.4132221087049370e-05, -1.4174402118255600e-04, 7.3486496691656049e-04, 2.3105816616276611e-05, -1.1983408465235623e-04,
        1.8044581161690771e-04, 5.9407987857315407e-07, 4.4769112358815145e-06, 4.4436566150092569e-17, -9.9469578895181997e-04, 1.6220276314313309e-04,
        -9.5710452873525644e-05, 1.2662875474526614e-03, 1.5602511542291521e-05, -2.0649060889963382e-04, 1.2184323687798982e-04, 1.0239596690627556e-06,
        7.2319800498266867e-06, 6.1369779999155662e-17, -1.2218883130103663e-03, 1.9925006776765812e-04, -3.3680823229324854e-05, 1.5555127229134207e-03,
        5.4906542718873632e-06, -2.5365331033436612e-04, 4.2877035890732656e-05, 1.2579210421270104e-06, 8.6095141065869613e-06, 7.1629786937454728e-17,
        2.1943348004296780e-04, -6.9776478809209231e-05, 1.0931001556629078e-04, -2.7934760190890233e-04, -3.4451757069426652e-05, 8.8828250006276136e-05,
        -1.3915602444573700e-04, 2.8472996610827392e-06, -1.3462461610689957e-05, -1.3423725475374167e-17, 4.3610842151565722e-04, -1.3819686293535999e-04,
        2.1413925090727119e-05, -5.5518347382007280e-04, -6.8888362637695031e-06, 1.7593013720955733e-04, -2.7260783634204059e-05, 6.0629566842156252e-06,
        -9.1041122183868957e-06, -1.3232981739875361e-17, 4.0432786569254695e-04, -1.2816820392757143e-04, -3.4987563703448329e-05, -5.1472555438689680e-04,
        1.1047540869101181e-05, 1.6316325294176266e-04, 4.4540568810561190e-05, 5.5618567638297611e-06, -5.4055108406481688e-06, -1.1067173402525196e-17,
        2.2452284971688618e-04, -7.1167813182509290e-05, -6.4833668634051413e-05, -2.8582657318237232e-04, 2.0511575929728835e-05, 9.0599474345227855e-05,
        8.2535854840079563e-05, 3.0981066058932078e-06, -2.3154572454531503e-06, -2.0250909035767044e-18, -1.9661496835268277e-05, 6.2364545034624598e-06,
        -7.2956872225952743e-05, 2.5029872332164925e-05, 2.3096613004447052e-05, -7.9392561682745486e-06, 9.2877018106373739e-05, -2.7032276650741835e-07,
        1.5500833780181726e-07, 1.3277001112846085e-17, -2.6133624469447495e-04, 8.2846618815267040e-05, -6.4181688197265032e-05, 3.3269149827544773e-04,
        2.0323581324052313e-05, -1.0546706130613412e-04, 8.1705857651540774e-05, -3.6005452741341418e-06, 2.0081691410669324e-06, 2.8571060802081716e-17,
        -4.5032633391014668e-04, 1.4275601491742435e-04, -4.3336922314422691e-05, 5.7328344530475027e-04, 1.3724632897094390e-05, -1.8173412014178887e-04,
        5.5169636466947789e-05, -6.2054078588159285e-06, 3.2435488692397214e-06, 3.8969589344185024e-17, -5.5318357068325397e-04, 1.7536157032519553e-04,
        -1.5250285570937746e-05, 7.0422482410401348e-04, 4.8299508614755072e-06, -2.2324229706304123e-04, 1.9414223854234778e-05, -7.6230039760334219e-06,
        3.8612474615728133e-06, 4.3284701132573816e-17, 3.7604074979975477e-05, -3.0294674009582567e-05, 1.8072883786917177e-05, -4.7871492379365758e-05,
        -1.5390330853862085e-05, 3.8566332418961975e-05, -2.3007504344675498e-05, 6.6491000063301673e-06, -2.8178579922956300e-06, -7.9095918696594176e-18,
        7.3270287126054955e-05, -6.0708847522411039e-05, 3.7405038611583371e-06, -9.3276007817184013e-05, -2.8745282594779060e-06, 7.7284792487990650e-05,
        -4.7618111116928247e-06, 1.4257402535386187e-05, -1.8716955469354936e-06, -3.9505250537844865e-18, 6.8002926475353674e-05, -5.6225462188821847e-05,
        -5.8069646962339915e-06, -8.6570446906970746e-05, 4.9002968245327234e-06, 7.1577263531495111e-05, 7.3924984553274739e-06, 1.3140222484471143e-05,
        -1.1167339315079637e-06, -1.0990652265370899e-17, 3.7749886887707968e-05, -3.1230076485260194e-05, -1.0834213338389640e-05, -4.8057116773323850e-05,
        9.0542379286307905e-06, 3.9757137205704310e-05, 1.3792387169260654e-05, 7.3084505808988420e-06, -4.7736440058837683e-07, -1.5484276100280180e-17,
        -3.3149552045122865e-06, 2.7320505889780681e-06, -1.2216240033828867e-05, 4.2200706411434401e-06, 1.0173710041948193e-05, -3.4780097375270807e-06,
        1.5551762461777675e-05, -6.3374425050783931e-07, 3.2586873484871532e-08, 7.1248334377215630e-18, -4.3958996718323173e-05, 3.6342108184710239e-05,
        -1.0755488950499986e-05, 5.5961562078764271e-05, 8.9457858333129674e-06, -4.6264958144617056e-05, 1.3692167872886584e-05, -8.4914149761200091e-06,
        4.1533005654660920e-07, 2.1226076432594638e-17, -7.5742210399340238e-05, 6.2625581143322420e-05, -7.2652034204399960e-06, 9.6422865071429752e-05,
        6.0388313855885739e-06, -7.9724871095856144e-05, 9.2488947105202630e-06, -1.4636620319129859e-05, 6.7044703542218178e-07, 1.4353850908910043e-17,
        -9.3040091887932651e-05, 7.6930184467904680e-05, -2.5570471988813422e-06, 1.1844376047443472e-04, 2.1248537374239493e-06, -9.7935203603904451e-05,
        3.2552234182113150e-06, -1.7981070634123197e-05, 7.9801219628634659e-07, 1.4947070801485985e-17,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_NxDx_DyDy(void)
{
  static const double sol[640] = {
        -7.4185325279800725e-04, 2.0456784067503705e-05, -3.9802597843864769e-04, 9.4440888007044066e-04, 1.0964684083799404e-05, -2.6042304806460913e-05,
        5.0670299970892744e-04, 5.0764812236798420e-06, 2.3457263346544559e-05, 1.1444387639127215e-16, -1.8111297630616463e-03, 4.9921555573646745e-05,
        -2.2122899042748763e-04, 2.3056406705032361e-03, 6.1001783803440449e-06, -6.3552138125513161e-05, 2.8163335848548543e-04, 1.2783176280256957e-05,
        2.1979419790247300e-05, 1.2089115540443858e-16, -2.2944506126898988e-03, 6.3246271865103235e-05, -6.1638761923246536e-05, 2.9209274547706849e-03,
        1.6996066367289599e-06, -8.0515035224970196e-05, 7.8468610735821001e-05, 1.6147340377878863e-05, 1.9018655925888982e-05, 1.3207856945143187e-16,
        -2.2712814698432750e-03, 6.2607182589647211e-05, 6.9283273198993004e-05, 2.8914322086887648e-03, -1.9084463650291406e-06, -7.9701449000710786e-05,
        -8.8200379526168931e-05, 1.5993584902368012e-05, 1.4578513220970464e-05, 1.2605949657091669e-16, -1.8607939734753515e-03, 5.1292513034896070e-05,
        1.6006863502989648e-04, 2.3688651979412701e-03, -4.4105638650403004e-06, -6.5297421839985193e-05, -2.0377377840282183e-04, 1.3099218422382167e-05,
        8.6581482607717041e-06, 1.0239515161122573e-16, -1.2218883130103891e-03, 3.3680823229323207e-05, 1.9925006776767170e-04, 1.5555127229134171e-03,
        -5.4906542718873810e-06, -4.2877035890733429e-05, -2.5365331033436693e-04, 8.6095141065858618e-06, 1.2579210421273283e-06, 6.3066757853797817e-17,
        -5.5318357068322437e-04, 1.5250285570931653e-05, 1.7536157032521370e-04, 7.0422482410400979e-04, -4.8299508614777137e-06, -1.9414223854238400e-05,
        -2.2324229706304001e-04, 3.8612474615743837e-06, -7.6230039760317829e-06, 1.7937266001530596e-18, -9.3040091887798182e-05, 2.5570471989134180e-06,
        7.6930184467945161e-05, 1.1844376047445425e-04, -2.1248537373841836e-06, -3.2552234182275937e-06, -9.7935203603889923e-05, 7.9801219626825619e-07,
        -1.7981070634123007e-05, -1.0128036117885212e-16, -6.0391518322645214e-04, 5.8133025438251359e-05, -3.2401990289692670e-04, 7.6880819717100853e-04,
        3.1157234702082137e-05, -7.4005667889303007e-05, 4.1249030379150915e-04, 4.2641058537420046e-06, 1.9094342915989936e-05, 5.7972050663881779e-17,
        -1.4743782859743769e-03, 1.4186131475814385e-04, -1.8009465973344265e-04, 1.8769425632445367e-03, 1.7335155668334535e-05, -1.8059513103275561e-04,
        2.2926770929982291e-04, 1.0737709998713822e-05, 1.7891815716327078e-05, 1.1195061398383187e-16, -1.8678324615434192e-03, 1.7972636255984167e-04,
        -5.0177806705345845e-05, 2.3778254749349014e-03, 4.8298278531591253e-06, -2.2879885225845310e-04, 6.3878356071482792e-05, 1.3563696440508467e-05,
        1.5481586715523658e-05, 1.2519746461652156e-16, -1.8489709303060073e-03, 1.7791020436131827e-04, 5.6401394648710322e-05, 2.3538139908237995e-03,
        -5.4230258469329165e-06, -2.2648681018836580e-04, -7.1801232593808891e-05, 1.3434536925262129e-05, 1.1867224005867195e-05, 1.1843790739414284e-16,
        -1.5148068254141455e-03, 1.4575747047953139e-04, 1.3030661548298953e-04, 1.9284097119174970e-03, -1.2533216551743824e-05, -1.8555509319178703e-04,
        -1.6588553643187407e-04, 1.1003322845518524e-05, 7.0478959388723687e-06, 9.8852723025181421e-17, -9.9469578895185944e-04, 9.5710452873520142e-05,
        1.6220276314314437e-04, 1.2662875474526525e-03, -1.5602511542291162e-05, -1.2184323687799192e-04, -2.0649060889963542e-04, 7.2319800498275634e-06,
        1.0239596690629926e-06, 6.8625507719155729e-17, -4.5032633391014229e-04, 4.3336922314407004e-05, 1.4275601491743747e-04, 5.7328344530472913e-04,
        -1.3724632897102337e-05, -5.5169636466954247e-05, -1.8173412014179445e-04, 3.2435488692356264e-06, -6.2054078588162292e-06, 2.9730292178322782e-17,
        -7.5742210399298401e-05, 7.2652034204017456e-06, 6.2625581143333262e-05, 9.6422865071383131e-05, -6.0388313856066335e-06, -9.2488947105253537e-06,
        -7.9724871095872353e-05, 6.7044703544074461e-07, -1.4636620319128287e-05, -2.7309579880584324e-17, -3.5046776690446905e-04, 8.6097341540724331e-05,
        -1.8804143086101184e-04, 4.4615949312715598e-04, 4.6138941825169657e-05, -1.0960536143110376e-04, 2.3938426697808543e-04, 2.6393373956046143e-06,
        1.1077684120010877e-05, 6.7512769569450989e-17, -8.5563118982642501e-04, 2.1009107982625130e-04, -1.0451452877112072e-04, 1.0892527473457733e-03,
        2.5673892807441163e-05, -2.6745435254661563e-04, 1.3305117783822521e-04, 6.6467490018670767e-06, 1.0380762902208861e-05, 1.2439980280761793e-16,
        -1.0839644414274159e-03, 2.6616909306344830e-04, -2.9119379734350206e-05, 1.3799301146204899e-03, 7.1531646979159201e-06, -3.3884390766165056e-04,
        3.7070135770882529e-05, 8.3963955227112072e-06, 8.9821168968255152e-06, 1.1281238326745632e-16, -1.0730176634536760e-03, 2.6347918201769812e-04,
        3.2732351109899887e-05, 1.3659944281656090e-03, -8.0305360855985600e-06, -3.3541954324910591e-04, -4.1669592924462712e-05, 8.3164350440211320e-06,
        6.8851645536151122e-06, 1.0152614645118684e-16, -8.7908999816593467e-04, 2.1586222088962838e-04, 7.5622040169425608e-05, 1.1191167491928647e-03,
        -1.8560279471199128e-05, -2.7480124608354468e-04, -9.6269883559479470e-05, 6.8115232787051463e-06, 4.0890266760257836e-06, 8.8210735506913166e-17,
        -5.7725205425153165e-04, 1.4174402118254949e-04, 9.4132221087060022e-05, 7.3486496691654532e-04, -2.3105816616276570e-05, -1.8044581161690717e-04,
        -1.1983408465235718e-04, 4.4769112358833051e-06, 5.9407987857323031e-07, 7.2862062998879852e-17, -2.6133624469449972e-04, 6.4181688197260221e-05,
        8.2846618815276906e-05, 3.3269149827542507e-04, -2.0323581324049829e-05, -8.1705857651543105e-05, -1.0546706130613814e-04, 2.0081691410600528e-06,
        -3.6005452741349643e-06, 6.1693126814669362e-17, -4.3958996718367131e-05, 1.0755488950523258e-05, 3.6342108184685594e-05, 5.5961562078727293e-05,
        -8.9457858332897807e-06, -1.3692167872848242e-05, -4.6264958144620241e-05, 4.1533005659194891e-07, -8.4914149761237141e-06, 8.2812149854959705e-17,
        -2.6369173463935267e-05, 9.7874810518100015e-05, -1.4156115092418241e-05, 3.3569013123087527e-05, 5.2437299075712193e-05, -1.2459855077796993e-04,
        1.8021301046369524e-05, 2.0225634453753787e-07, 8.2737327431625208e-07, 2.8177397077573526e-16, -6.4401076442197511e-05, 2.3880385491959882e-04,
        -7.8652359084217054e-06, 8.1985147666929623e-05, 2.9185828660984950e-05, -3.0400686433693097e-04, 1.0012760081337085e-05, 5.1040642097596895e-07,
        7.7809221148224797e-07, 1.2817421272705128e-16, -8.1583353672124066e-05, 3.0254949112893249e-04, -2.1908826904486396e-06, 1.0385887422229175e-04,
        8.1313422416425124e-06, -3.8515761035684960e-04, 2.7890813449622491e-06, 6.4546625765873201e-07, 6.7280421817243716e-07, 9.2197833111684565e-17,
        -8.0757424528256949e-05, 2.9949136733266134e-04, 2.4649783637441816e-06, 1.0280743336822311e-04, -9.1266848261750555e-06, -3.8126449637687881e-04,
        -3.1380161064341044e-06, 6.3927833037391234e-07, 5.1575107170366369e-07, 7.5874538131131677e-17, -6.6158684055516007e-05, 2.4536651900581257e-04,
        5.6930861126352351e-06, 8.4222652499130987e-05, -2.1095228673489812e-05, -3.1236139835906734e-04, -7.2475264609201958e-06, 5.2383308233429554e-07,
        3.0620542597919858e-07, 6.8193629092451981e-17, -4.3439872336668323e-05, 1.6111739776338878e-04, 7.0856641229971160e-06, 5.5300695965283921e-05,
        -2.6262233432286572e-05, -2.0510889533447113e-04, -9.0203339996343701e-06, 3.4429239733307118e-07, 4.4421321296388115e-08, 6.6096051000889153e-17,
        -1.9661496835318848e-05, 7.2956872225944395e-05, 6.2364545034732409e-06, 2.5029872332157129e-05, -2.3096613004447269e-05, -9.2877018106369294e-05,
        -7.9392561682712807e-06, 1.5500833779624817e-07, -2.7032276650636750e-07, 7.6653401639516628e-17, -3.3149552044579930e-06, 1.2216240033808677e-05,
        2.7320505890367992e-06, 4.2200706410946925e-06, -1.0173710041956921e-05, -1.5551762461784462e-05, -3.4780097375448968e-06, 3.2586873488663109e-08,
        -6.3374425049983369e-07, 1.1208704012626897e-16, 3.0109676195024578e-04, 8.6992163559396031e-05, 1.6152938865388649e-04, -3.8330822797395835e-04,
        4.6572459248561327e-05, -1.1074450567149586e-04, -2.0563337622614048e-04, -3.0476196626152295e-06, -9.5344390692388496e-06, 1.9092368042348901e-16,
        7.3504616021194998e-04, 2.1219385883252588e-04, 8.9786913323488897e-05, -9.3574317878677841e-04, 2.5939534318833506e-05, -2.7013127437548854e-04,
        -1.1430233396836451e-04, -7.6719397877112796e-06, -8.9337392931369159e-06, 9.5643625673969897e-17, 9.3120558827590547e-04, 2.6884245930155592e-04,
        2.5018632896854560e-05, -1.1854619810895955e-03, 7.2279537638824679e-06, -3.4224720987189417e-04, -3.1849720933210353e-05, -9.6892163186664671e-06,
        -7.7309086961064527e-06, 5.6028473148456790e-17, 9.2180607376450964e-04, 2.6612401553145846e-04, -2.8115105526505151e-05, -1.1734960229442622e-03,
        -8.1052313173919847e-06, -3.3878652215938403e-04, 3.5791654512842331e-05, -9.5969568198564807e-06, -5.9258354008360331e-06, 4.5038336270219261e-17,
        7.5521454677789578e-04, 2.1803021136886734e-04, -6.4959947881067839e-05, -9.6141834203178659e-04, -1.8739033755342080e-05, -2.7756118472743892e-04,
        8.2696613375290293e-05, -7.8598788404716561e-06, -3.5194659963033452e-06, 4.4233911651540140e-17, 4.9591678847464625e-04, 1.4316667347179192e-04,
        -8.0861991798903798e-05, -6.3132191851337973e-04, -2.3329868400283943e-05, -1.8225690491621515e-04, 1.0294055168874891e-04, -5.1659110589369388e-06,
        -5.1117197886668249e-07, 4.9241705811745897e-17, 2.2452284971679964e-04, 6.4833668634040463e-05, -7.1167813182507325e-05, -2.8582657318238186e-04,
        -2.0511575929733324e-05, -8.2535854840072923e-05, 9.0599474345226703e-05, -2.3154572454575781e-06, 3.0981066058919801e-06, 6.5702071410912399e-17,
        3.7749886887585332e-05, 1.0834213338346070e-05, -3.1230076485287651e-05, -4.8057116773306991e-05, -9.0542379286509143e-06, -1.3792387169287603e-05,
        3.9757137205718499e-05, -4.7736440055484852e-07, 7.3084505808939707e-06, 1.0203989573196756e-16, 5.4219548657402584e-04, 4.6964511263801802e-05,
        2.9091105641275715e-04, -6.9023655328555940e-04, 2.5099762187725797e-05, -5.9787702377005517e-05, -3.7034141718845218e-04, -7.1078530313566328e-06,
        -1.7138457671768140e-05, 1.9161885033206326e-18, 1.3236741350855599e-03, 1.1444951102741669e-04, 1.6169038423836617e-04, -1.6850901479239001e-03,
        1.4004905074152640e-05, -1.4569880785233498e-04, -2.0583833004826857e-04, -1.7897506519471960e-05, -1.6042765456802161e-05, 1.4300776630829959e-17,
        1.6769291755919155e-03, 1.4502400361685275e-04, 4.5047880146265376e-05, -2.1347979518942612e-03, 3.8980585320075360e-06, -1.8462136052192925e-04,
        -5.7347754260088115e-05, -2.2607359364032996e-05, -1.3884451253031794e-05, 1.5312228686440909e-17, 1.6599996084286732e-03, 1.4355365040744238e-04,
        -5.0637183330136170e-05, -2.1132459353675054e-03, -4.3684028831049219e-06, -1.8274954204225104e-04, 6.4463160899273976e-05, -2.2392186132941893e-05,
        -1.0642774281307195e-05, 1.6238941250461819e-17, 1.3599971624203883e-03, 1.1761250922650970e-04, -1.1698967736129490e-04, -1.7313308153829791e-03,
        -1.0103627331035387e-05, -1.4972543114425037e-04, 1.4893254125386701e-04, -1.8339153449490943e-05, -6.3208097101488461e-06, 1.9476238450815159e-17,
        8.9304807164369880e-04, 7.7225405099946121e-05, -1.4562879398959312e-04, -1.1368859353378387e-03, -1.2582857886101920e-05, -9.8311031283335407e-05,
        1.8539128287040863e-04, -1.2053450101605257e-05, -9.1822109815714625e-07, 2.5781257654841598e-17, 4.0432786569244986e-04, 3.4987563703455363e-05,
        -1.2816820392757135e-04, -5.1472555438690146e-04, -1.1047540869096394e-05, -4.4540568810559333e-05, 1.6316325294176082e-04, -5.4055108406509953e-06,
        5.5618567638278036e-06, 3.8946897289047353e-17, 6.8002926475210410e-05, 5.8069646962579202e-06, -5.6225462188851432e-05, -8.6570446906967778e-05,
        -4.9002968245208387e-06, -7.3924984553276374e-06, 7.1577263531497699e-05, -1.1167339314804364e-06, 1.3140222484467003e-05, 6.6267225159319282e-17,
        5.8489869202466081e-04, -2.8617271619844676e-05, 3.1367048785920317e-04, -7.4459944282312168e-04, -1.5601216385088299e-05, 3.6430932046641732e-05,
        -3.9931508426116090e-04, -1.1993596362323995e-05, -1.8606531928345532e-05, -1.6210055668937958e-16, 1.4277325347644965e-03, -7.0228009414114608e-05,
        1.7441241606388127e-04, -1.8175606551731254e-03, -8.5314378018309973e-06, 8.9403066536796741e-05, -2.2203398570289293e-04, -3.0180009723339873e-05,
        -1.7465893456541061e-05, -4.3318933013182567e-17, 1.8087408305902700e-03, -8.8933756155081838e-05, 4.8609262873974651e-05, -2.3025994638613892e-03,
        -2.3689535246149096e-06, 1.1321623074936765e-04, -6.1881537000414255e-05, -3.8108004663828739e-05, -1.5106133674190260e-05, -1.2329061675074945e-17,
        1.7904898493183250e-03, -8.8038739116152721e-05, -5.4594815714690414e-05, -2.2793652342907728e-03, 2.7188239063321850e-06, 1.1207683823987290e-04,
        6.9501385310799283e-05, -3.7747530559389502e-05, -1.1580436012148740e-05, -4.7213831098352963e-18, 1.4669018971550354e-03, -7.2126947445435468e-05,
        -1.2615623684791142e-04, -1.8674248210697169e-03, 6.2473383435735753e-06, 9.1820490646871939e-05, 1.6060193832964490e-04, -3.0915439591943777e-05,
        -6.8766677416190290e-06, 2.3884772252174143e-18, 9.6325796444311643e-04, -4.7364249466621374e-05, -1.5703540377680101e-04, -1.2262659386990189e-03,
        7.7720809383160633e-06, 6.0296585106922731e-05, 1.9991235362653838e-04, -2.0321917006824412e-05, -9.9955532671635744e-07, 8.8559730946373446e-18,
        4.3610842151558675e-04, -2.1413925090714942e-05, -1.3819686293535334e-04, -5.5518347382007356e-04, 6.8888362637715224e-06, 2.7260783634202982e-05,
        1.7593013720955538e-04, -9.1041122183867737e-06, 6.0629566842156548e-06, 1.5573035261013121e-17, 7.3270287126041931e-05, -3.7405038611519128e-06,
        -6.0708847522380898e-05, -9.3276007817195505e-05, 2.8745282594694069e-06, 4.7618111116912238e-06, 7.7284792487987154e-05, -1.8716955469427051e-06,
        1.4257402535389059e-05, 8.5076182105527532e-18, 2.9366410440316701e-04, -1.4679523967258041e-04, 1.5833417110366465e-04, -3.7384615745768692e-04,
        -7.8278512648504361e-05, 1.8687621490632084e-04, -2.0156573641119286e-04, -1.7623476288564741e-05, -8.6854584349435904e-06, -1.2727645769649134e-16,
        7.1871323474263348e-04, -3.5755576442369596e-04, 8.7710214994732844e-05, -9.1495071101393652e-04, -4.3774882151098579e-05, 4.5518279763370506e-04,
        -1.1165861388559535e-04, -4.4456825652019548e-05, -8.1880463208024440e-06, -1.3893386733482794e-17, 9.1041892369559275e-04, -4.5307301340345680e-04,
        2.4407192087609418e-05, -1.1589997251879149e-03, -1.2201146101918565e-05, 5.7678007822305362e-04, -3.1071332312990467e-05, -5.6202873171461935e-05,
        -7.0756666652780254e-06, -3.2435527423937891e-17, 9.0124227629166816e-04, -4.4848460248175293e-04, -2.7574015762470031e-05, -1.1473174857899882e-03,
        1.3619750416026979e-05, 5.7093884748954675e-04, 3.5102825588581247e-05, -5.5653119508575166e-05, -5.4247562557556097e-06, -7.4236581082921709e-18,
        7.3835886561779212e-04, -3.6743879842918988e-04, -6.3621467426429139e-05, -9.3996038534390695e-04, 3.1529332913748146e-05, 4.6776429544570465e-04,
        8.0992674190014512e-05, -4.5587253648867070e-05, -3.2207980489842500e-06, 3.0495007042036106e-18, 4.8485609369284887e-04, -2.4126875678009098e-04,
        -7.9184942178625679e-05, -6.1724121140272674e-04, 3.9263369648520594e-05, 3.0714478305167810e-04, 1.0080560040594569e-04, -2.9948261566550803e-05,
        -4.6828117363520094e-07, 2.1116594102843436e-18, 2.1943348004293224e-04, -1.0931001556628331e-04, -6.9776478809210180e-05, -2.7934760190890894e-04,
        3.4451757069420533e-05, 1.3915602444573510e-04, 8.8828250006272802e-05, -1.3462461610690233e-05, 2.8472996610841432e-06, 2.7189423663698649e-18,
        3.7604074979956313e-05, -1.8072883786921466e-05, -3.0294674009569994e-05, -4.7871492379370643e-05, 1.5390330853863688e-05, 2.3007504344676501e-05,
        3.8566332418963669e-05, -2.8178579922984523e-06, 6.6491000063332310e-06, 9.1249273016952658e-19,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_DxNx_DyDy(void)
{
  static const double sol[640] = {
        2.9366410440283302e-04, 1.4679523967220964e-04, 1.5833417110396145e-04, -3.7384615745763466e-04, 7.8278512648638829e-05, -1.8687621490641275e-04,
        -2.0156573641115063e-04, -1.7623476288484310e-05, -8.6854584350381429e-06, 2.9951311795116681e-17, 7.1871323474270667e-04, 3.5755576442359968e-04,
        8.7710214994801027e-05, -9.1495071101377866e-04, 4.3774882151080975e-05, -4.5518279763380145e-04, -1.1165861388558839e-04, -4.4456825652005657e-05,
        -8.1880463207962555e-06, -3.9361218137492589e-17, 9.1041892369561194e-04, 4.5307301340353107e-04, 2.4407192087611417e-05, -1.1589997251881807e-03,
        1.2201146101977434e-05, -5.7678007822299421e-04, -3.1071332313249856e-05, -5.6202873171480712e-05, -7.0756666651934822e-06, -1.8044717852722420e-16,
        9.0124227629218500e-04, 4.4848460248163091e-04, -2.7574015762270924e-05, -1.1473174857904726e-03, -1.3619750416150265e-05, -5.7093884748939789e-04,
        3.5102825588767622e-05, -5.5653119508558761e-05, -5.4247562557377872e-06, -1.0868209604022288e-16, 7.3835886561839309e-04, 3.6743879842905105e-04,
        -6.3621467426616747e-05, -9.3996038534400355e-04, -3.1529332913601345e-05, -4.6776429544575414e-04, 8.0992674190034908e-05, -4.5587253648843827e-05,
        -3.2207980489953576e-06, -1.6498451698665985e-17, 4.8485609369309200e-04, 2.4126875678013410e-04, -7.9184942178653760e-05, -6.1724121140280101e-04,
        -3.9263369648527702e-05, -3.0714478305167691e-04, 1.0080560040600143e-04, -2.9948261566577983e-05, -4.6828117365321368e-07, -2.7907473466889570e-16,
        2.1943348004327650e-04, 1.0931001556619420e-04, -6.9776478809175255e-05, -2.7934760190897551e-04, -3.4451757069465494e-05, -1.3915602444572510e-04,
        8.8828250006243000e-05, -1.3462461610682388e-05, 2.8472996610262629e-06, -2.8517753439349179e-17, 3.7604074980109599e-05, 1.8072883786862083e-05,
        -3.0294674009660176e-05, -4.7871492379410237e-05, -1.5390330853845416e-05, -2.3007504344661451e-05, 3.8566332418977053e-05, -2.8178579922810013e-06,
        6.6491000063179548e-06, 1.8201534010187335e-16, 5.8489869202416793e-04, 2.8617271619977579e-05, 3.1367048785938028e-04, -7.4459944282304351e-04,
        1.5601216385000831e-05, -3.6430932046655128e-05, -3.9931508426120242e-04, -1.1993596362356830e-05, -1.8606531928369523e-05, 1.1609094708570108e-16,
        1.4277325347643946e-03, 7.0228009414118281e-05, 1.7441241606394765e-04, -1.8175606551731423e-03, 8.5314378018203992e-06, -8.9403066536814779e-05,
        -2.2203398570290480e-04, -3.0180009723319791e-05, -1.7465893456551608e-05, 7.9952680116957346e-18, 1.8087408305903605e-03, 8.8933756155068312e-05,
        4.8609262874017952e-05, -2.3025994638615327e-03, 2.3689535246137682e-06, -1.1321623074934033e-04, -6.1881537000456281e-05, -3.8108004663832317e-05,
        -1.5106133674201696e-05, -3.4437911683227506e-17, 1.7904898493185757e-03, 8.8038739116105775e-05, -5.4594815714643902e-05, -2.2793652342909749e-03,
        -2.7188239063507206e-06, -1.1207683823983096e-04, 6.9501385310806059e-05, -3.7747530559384820e-05, -1.1580436012161972e-05, -5.3751195984836405e-17,
        1.4669018971553392e-03, 7.2126947445379482e-05, -1.2615623684791947e-04, -1.8674248210698470e-03, -6.2473383435676977e-06, -9.1820490646863468e-05,
        1.6060193832966014e-04, -3.0915439591938444e-05, -6.8766677416264973e-06, -7.2829702745483415e-17, 9.6325796444333490e-04, 4.7364249466599229e-05,
        -1.5703540377681876e-04, -1.2262659386990920e-03, -7.7720809383183452e-06, -6.0296585106921430e-05, 1.9991235362654895e-04, -2.0321917006824303e-05,
        -9.9955532670582353e-07, -5.6107483444491660e-17, 4.3610842151573615e-04, 2.1413925090684689e-05, -1.3819686293538045e-04, -5.5518347382010706e-04,
        -6.8888362637714928e-06, -2.7260783634201552e-05, 1.7593013720956378e-04, -9.1041122183843021e-06, 6.0629566842222413e-06, -1.1946646468031830e-17,
        7.3270287126065377e-05, 3.7405038611191754e-06, -6.0708847522435975e-05, -9.3276007817201794e-05, -2.8745282594642544e-06, -4.7618111116780999e-06,
        7.7284792488002224e-05, -1.8716955469385451e-06, 1.4257402535388324e-05, 4.2151661090997407e-17, 5.4219548657403398e-04, -4.6964511263722479e-05,
        2.9091105641265626e-04, -6.9023655328565383e-04, -2.5099762187774952e-05, 5.9787702377004467e-05, -3.7034141718838637e-04, -7.1078530314491457e-06,
        -1.7138457671755167e-05, 7.2291469764503234e-17, 1.3236741350854942e-03, -1.1444951102741537e-04, 1.6169038423840810e-04, -1.6850901479239229e-03,
        -1.4004905074144496e-05, 1.4569880785234056e-04, -2.0583833004829050e-04, -1.7897506519463947e-05, -1.6042765456801554e-05, 8.1778256449745880e-17,
        1.6769291755919702e-03, -1.4502400361685578e-04, 4.5047880146293769e-05, -2.1347979518943441e-03, -3.8980585320151169e-06, 1.8462136052193719e-04,
        -5.7347754260101918e-05, -2.2607359364033162e-05, -1.3884451253030498e-05, 3.1613783950993451e-17, 1.6599996084288107e-03, -1.4355365040746436e-04,
        -5.0637183330118932e-05, -2.1132459353676147e-03, 4.3684028831022470e-06, 1.8274954204226714e-04, 6.4463160899275100e-05, -2.2392186132939819e-05,
        -1.0642774281307588e-05, -2.0569040606703768e-18, 1.3599971624205561e-03, -1.1761250922653675e-04, -1.1698967736129657e-04, -1.7313308153830715e-03,
        1.0103627331037176e-05, 1.4972543114426227e-04, 1.4893254125387590e-04, -1.8339153449488886e-05, -6.3208097101509027e-06, -1.6520376512049561e-17,
        8.9304807164383899e-04, -7.7225405099968049e-05, -1.4562879398960886e-04, -1.1368859353378985e-03, 1.2582857886104629e-05, 9.8311031283340679e-05,
        1.8539128287041925e-04, -1.2053450101604032e-05, -9.1822109816026481e-07, -1.1776408068596795e-17, 4.0432786569251865e-04, -3.4987563703472656e-05,
        -1.2816820392759742e-04, -5.1472555438692672e-04, 1.1047540869096497e-05, 4.4540568810563562e-05, 1.6316325294177080e-04, -5.4055108406492386e-06,
        5.5618567638242114e-06, 1.2197631094270155e-17, 6.8002926475168316e-05, -5.8069646962690798e-06, -5.6225462188887827e-05, -8.6570446906960555e-05,
        4.9002968245240032e-06, 7.3924984553178245e-06, 7.1577263531512648e-05, -1.1167339314809732e-06, 1.3140222484464619e-05, 5.8504681465728778e-17,
        3.0109676194992063e-04, -8.6992163559340913e-05, 1.6152938865410021e-04, -3.8330822797374379e-04, -4.6572459248559639e-05, 1.1074450567149812e-04,
        -2.0563337622629170e-04, -3.0476196624545761e-06, -9.5344390692410756e-06, 2.5658063100360113e-16, 7.3504616021197644e-04, -2.1219385883250536e-04,
        8.9786913323481145e-05, -9.3574317878682058e-04, -2.5939534318848258e-05, 2.7013127437548311e-04, -1.1430233396836765e-04, -7.6719397877277883e-06,
        -8.9337392931367838e-06, 1.3999285833550340e-16, 9.3120558827594559e-04, -2.6884245930155841e-04, 2.5018632896869424e-05, -1.1854619810896538e-03,
        -7.2279537638851598e-06, 3.4224720987189987e-04, -3.1849720933216418e-05, -9.6892163186646901e-06, -7.7309086961071439e-06, 8.5740907185919959e-17,
        9.2180607376459356e-04, -2.6612401553146941e-04, -2.8115105526496142e-05, -1.1734960229443305e-03, 8.1052313173897706e-06, 3.3878652215939243e-04,
        3.5791654512842690e-05, -9.5969568198558014e-06, -5.9258354008375010e-06, 4.8242096485927886e-17, 7.5521454677799423e-04, -2.1803021136888130e-04,
        -6.4959947881067893e-05, -9.6141834203184470e-04, 1.8739033755341866e-05, 2.7756118472744694e-04, 8.2696613375294535e-05, -7.8598788404704144e-06,
        -3.5194659963044798e-06, 2.6526019995138209e-17, 4.9591678847472865e-04, -1.4316667347180327e-04, -8.0861991798912851e-05, -6.3132191851341746e-04,
        2.3329868400285224e-05, 1.8225690491622097e-04, 1.0294055168875558e-04, -5.1659110589356979e-06, -5.1117197886743053e-07, 1.7518231999671528e-17,
        2.2452284971683406e-04, -6.4833668634043566e-05, -7.1167813182526515e-05, -2.8582657318240094e-04, 2.0511575929736604e-05, 8.2535854840080106e-05,
        9.0599474345229793e-05, -2.3154572454563152e-06, 3.0981066058906943e-06, 2.0293341889326816e-17, 3.7749886887509071e-05, -1.0834213338340566e-05,
        -3.1230076485337226e-05, -4.8057116773339659e-05, 9.0542379286492321e-06, 1.3792387169257779e-05, 3.9757137205697967e-05, -4.7736440054448337e-07,
        7.3084505808889546e-06, 5.1380996417200149e-17, -2.6369173463370998e-05, -9.7874810518086666e-05, -1.4156115092791256e-05, 3.3569013122789893e-05,
        -5.2437299075761091e-05, 1.2459855077792001e-04, 1.8021301046546994e-05, 2.0225634435349193e-07, 8.2737327434274039e-07, 2.9408986913590118e-16,
        -6.4401076442212866e-05, -2.3880385491961627e-04, -7.8652359084084544e-06, 8.1985147666903182e-05, -2.9185828660976991e-05, 3.0400686433693043e-04,
        1.0012760081321650e-05, 5.1040642098054282e-07, 7.7809221148924742e-07, 1.8801925825746067e-16, -8.1583353672089317e-05, -3.0254949112893596e-04,
        -2.1908826904423915e-06, 1.0385887422223798e-04, -8.1313422416431392e-06, 3.8515761035684808e-04, 2.7890813449623851e-06, 6.4546625765844063e-07,
        6.7280421817217426e-07, 1.2336526352308237e-16, -8.0757424528202224e-05, -2.9949136733266730e-04, 2.4649783637492719e-06, 1.0280743336817254e-04,
        9.1266848261748081e-06, 3.8126449637688130e-04, -3.1380161064321046e-06, 6.3927833037451373e-07, 5.1575107170323075e-07, 8.6261223990707393e-17,
        -6.6158684055452865e-05, -2.4536651900581988e-04, 5.6930861126347582e-06, 8.4222652499092417e-05, 2.1095228673489606e-05, 3.1236139835907060e-04,
        -7.2475264609154193e-06, 5.2383308233506750e-07, 3.0620542597858777e-07, 5.8835010400590615e-17, -4.3439872336614058e-05, -1.6111739776339496e-04,
        7.0856641229928402e-06, 5.5300695965266336e-05, 2.6262233432287653e-05, 2.0510889533447631e-04, -9.0203339996267705e-06, 3.4429239733317430e-07,
        4.4421321296056952e-08, 3.7118787564661390e-17, -1.9661496835279546e-05, -7.2956872225940153e-05, 6.2364545034698510e-06, 2.5029872332166819e-05,
        2.3096613004453422e-05, 9.2877018106378645e-05, -7.9392561682623682e-06, 1.5500833779623422e-07, -2.7032276650592434e-07, 1.8444121589638132e-17,
        -3.3149552044149548e-06, -1.2216240033765435e-05, 2.7320505890445826e-06, 4.2200706411741357e-06, 1.0173710041980746e-05, 1.5551762461853329e-05,
        -3.4780097375158219e-06, 3.2586873483438240e-08, -6.3374425049766073e-07, -1.3804174291048425e-17, -3.5046776690512922e-04, -8.6097341540895459e-05,
        -1.8804143086066668e-04, 4.4615949312742958e-04, -4.6138941825076984e-05, 1.0960536143114947e-04, 2.3938426697787092e-04, 2.6393373958459857e-06,
        1.1077684119986845e-05, 2.8754810140367706e-16, -8.5563118982646936e-04, -2.1009107982628570e-04, -1.0451452877109037e-04, 1.0892527473456870e-03,
        -2.5673892807405829e-05, 2.6745435254660658e-04, 1.3305117783823987e-04, 6.6467490018443551e-06, 1.0380762902200129e-05, 1.8034621692091753e-16,
        -1.0839644414273974e-03, -2.6616909306344407e-04, -2.9119379734335553e-05, 1.3799301146204357e-03, -7.1531646979175727e-06, 3.3884390766164969e-04,
        3.7070135770883186e-05, 8.3963955227189068e-06, 8.9821168968238431e-06, 1.3648290180391621e-16, -1.0730176634536303e-03, -2.6347918201769844e-04,
        3.2732351109902035e-05, 1.3659944281655622e-03, 8.0305360855978434e-06, 3.3541954324910601e-04, -4.1669592924459006e-05, 8.3164350440210235e-06,
        6.8851645536143642e-06, 1.0951181679474486e-16, -8.7908999816589000e-04, -2.1586222088963193e-04, 7.5622040169423250e-05, 1.1191167491928352e-03,
        1.8560279471198105e-05, 2.7480124608354631e-04, -9.6269883559473628e-05, 6.8115232787058816e-06, 4.0890266760255126e-06, 8.2352293683465305e-17,
        -5.7725205425149739e-04, -1.4174402118255575e-04, 9.4132221087056837e-05, 7.3486496691653816e-04, 2.3105816616275750e-05, 1.8044581161690863e-04,
        -1.1983408465234938e-04, 4.4769112358827155e-06, 5.9407987857321188e-07, 5.1878348453767501e-17, -2.6133624469446931e-04, -6.4181688197269951e-05,
        8.2846618815280525e-05, 3.3269149827545961e-04, 2.0323581324048325e-05, 8.1705857651539975e-05, -1.0546706130612097e-04, 2.0081691410597411e-06,
        -3.6005452741329560e-06, 1.9482062592479926e-17, -4.3958996718269349e-05, -1.0755488950552287e-05, 3.6342108184730914e-05, 5.5961562078827826e-05,
        8.9457858332771938e-06, 1.3692167872827874e-05, -4.6264958144590181e-05, 4.1533005657324436e-07, -8.4914149761139613e-06, -3.6675292072574968e-17,
        -6.0391518322665716e-04, -5.8133025437974522e-05, -3.2401990289680299e-04, 7.6880819717064326e-04, -3.1157234702201196e-05, 7.4005667889288912e-05,
        4.1249030379167482e-04, 4.2641058538617344e-06, 1.9094342915962170e-05, 1.4047071619921544e-16, -1.4743782859743583e-03, -1.4186131475806473e-04,
        -1.8009465973341197e-04, 1.8769425632444658e-03, -1.7335155668375508e-05, 1.8059513103273111e-04, 2.2926770929981871e-04, 1.0737709998697717e-05,
        1.7891815716318834e-05, 1.2513825497453909e-16, -1.8678324615433617e-03, -1.7972636255982402e-04, -5.0177806705344206e-05, 2.3778254749348320e-03,
        -4.8298278531647810e-06, 2.2879885225845213e-04, 6.3878356071490829e-05, 1.3563696440515171e-05, 1.5481586715522774e-05, 1.2872813679253942e-16,
        -1.8489709303059541e-03, -1.7791020436131418e-04, 5.6401394648706507e-05, 2.3538139908237531e-03, 5.4230258469302585e-06, 2.2648681018836461e-04,
        -7.1801232593804391e-05, 1.3434536925261634e-05, 1.1867224005867065e-05, 1.1871040135927971e-16, -1.5148068254141093e-03, -1.4575747047953298e-04,
        1.3030661548298330e-04, 1.9284097119174682e-03, 1.2533216551742667e-05, 1.8555509319178608e-04, -1.6588553643186879e-04, 1.1003322845519037e-05,
        7.0478959388720587e-06, 9.6340345162047416e-17, -9.9469578895184643e-04, -9.5710452873525712e-05, 1.6220276314313697e-04, 1.2662875474526423e-03,
        1.5602511542290004e-05, 1.2184323687798953e-04, -2.0649060889962924e-04, 7.2319800498276379e-06, 1.0239596690627012e-06, 6.3275699227690185e-17,
        -4.5032633391016435e-04, -4.3336922314424310e-05, 1.4275601491742489e-04, 5.7328344530474149e-04, 1.3724632897095396e-05, 5.5169636466947329e-05,
        -1.8173412014178784e-04, 3.2435488692378338e-06, -6.2054078588167645e-06, 2.8221989307501208e-17, -7.5742210399375013e-05, -7.2652034204510260e-06,
        6.2625581143304829e-05, 9.6422865071426825e-05, 6.0388313855867765e-06, 9.2488947105098056e-06, -7.9724871095863313e-05, 6.7044703543938449e-07,
        -1.4636620319136193e-05, 5.1691853407524822e-18, -7.4185325279727617e-04, -2.0456784067398124e-05, -3.9802597843889169e-04, 9.4440888007032834e-04,
        -1.0964684083796703e-05, 2.6042304806468720e-05, 5.0670299970890598e-04, 5.0764812236771154e-06, 2.3457263346578311e-05, -1.0073674248609138e-16,
        -1.8111297630613963e-03, -4.9921555573597583e-05, -2.2122899042755054e-04, 2.3056406705031186e-03, -6.1001783803646252e-06, 6.3552138125520804e-05,
        2.8163335848550429e-04, 1.2783176280236964e-05, 2.1979419790258758e-05, 6.6977617418593946e-17, -2.2944506126897760e-03, -6.3246271865093260e-05,
        -6.1638761923269263e-05, 2.9209274547706120e-03, -1.6996066367323673e-06, 8.0515035224969721e-05, 7.8468610735829105e-05, 1.6147340377878379e-05,
        1.9018655925890873e-05, 1.1423587423447701e-16, -2.2712814698432077e-03, -6.2607182589644026e-05, 6.9283273198982406e-05, 2.8914322086887149e-03,
        1.9084463650277778e-06, 7.9701449000711097e-05, -8.8200379526163659e-05, 1.5993584902367006e-05, 1.4578513220971289e-05, 1.1911881276741378e-16,
        -1.8607939734753178e-03, -5.1292513034896659e-05, 1.6006863502988669e-04, 2.3688651979412380e-03, 4.4105638650393805e-06, 6.5297421839984935e-05,
        -2.0377377840281700e-04, 1.3099218422382080e-05, 8.6581482607716177e-06, 1.0056404625464723e-16, -1.2218883130103921e-03, -3.3680823229325111e-05,
        1.9925006776766002e-04, 1.5555127229134004e-03, 5.4906542718873293e-06, 4.2877035890732338e-05, -2.5365331033436308e-04, 8.6095141065870595e-06,
        1.2579210421268802e-06, 6.9454607065394521e-17, -5.5318357068327587e-04, -1.5250285570936698e-05, 1.7536157032519589e-04, 7.0422482410400209e-04,
        4.8299508614765804e-06, 1.9414223854234138e-05, -2.2324229706303925e-04, 3.8612474615724914e-06, -7.6230039760333550e-06, 3.8651870276418839e-17,
        -9.3040091887944807e-05, -2.5570471988800072e-06, 7.6930184467910291e-05, 1.1844376047442972e-04, 2.1248537374226546e-06, 3.2552234182123116e-06,
        -9.7935203603902093e-05, 7.9801219628591905e-07, -1.7981070634122831e-05, 1.2404970576903513e-17,
  };
  return sol;
}

static const double *
get_sol_3x_bias_p2_NxDx_PyPy(void)
{
  static const double sol[640] = {
        2.6007376622747585e-02, -7.1800944496270788e-04, 1.1848784811247477e-02, -3.3108431266206981e-02, -3.2577961081452374e-04, 9.1405475845812636e-04,
        -1.5083977257749990e-02, -1.7357322650915152e-04, -2.4528231596292736e-03, -2.4878168330471732e-16, 2.6007376622748265e-02, -7.1800944496270181e-04,
        -1.1848784811247111e-02, -3.3108431266206918e-02, 3.2577961081452483e-04, 9.1405475845813167e-04, 1.5083977257749992e-02, -1.7357322650915244e-04,
        -2.4528231596292914e-03, -2.9270179092558335e-16, -2.6007376622746436e-02, 7.1800944496270365e-04, -1.1848784811247128e-02, 3.3108431266207619e-02,
        3.2577961081452374e-04, -9.1405475845813221e-04, 1.5083977257749993e-02, 1.7357322650914987e-04, 2.4528231596292996e-03, 4.7189503890956915e-16,
        -2.6007376622746225e-02, 7.1800944496270365e-04, 1.1848784811247218e-02, 3.3108431266207473e-02, -3.2577961081452266e-04, -9.1405475845812983e-04,
        -1.5083977257750094e-02, 1.7357322650915022e-04, 2.4528231596292701e-03, 3.0312274075053911e-16, 2.6007376622748459e-02, -7.1800944496270528e-04,
        1.1848784811247060e-02, -3.3108431266207390e-02, -3.2577961081452564e-04, 9.1405475845812646e-04, -1.5083977257750094e-02, -1.7357322650914743e-04,
        -2.4528231596292814e-03, -8.5723427941119742e-16, 2.6007376622748039e-02, -7.1800944496270603e-04, -1.1848784811247305e-02, -3.3108431266207404e-02,
        3.2577961081452428e-04, 9.1405475845812983e-04, 1.5083977257750113e-02, -1.7357322650914767e-04, -2.4528231596292857e-03, -8.5960190237677279e-16,
        -2.6007376622747401e-02, 7.1800944496269888e-04, -1.1848784811247390e-02, 3.3108431266207647e-02, 3.2577961081452456e-04, -9.1405475845813178e-04,
        1.5083977257750217e-02, 1.7357322650914784e-04, 2.4528231596292831e-03, 4.3789679653681983e-16, -2.6007376622747918e-02, 7.1800944496269563e-04,
        1.1848784811247176e-02, 3.3108431266207820e-02, -3.2577961081452510e-04, -9.1405475845813199e-04, -1.5083977257750136e-02, 1.7357322650914862e-04,
        2.4528231596293512e-03, 6.5295682528229116e-16, 2.1167581926037952e-02, -2.0404149472422264e-03, 9.6451162689074419e-03, -2.6947178926806876e-02,
        -9.2582118178574223e-04, 2.5975298860486684e-03, -1.2278619011668608e-02, -1.4581758526841442e-04, -1.9953554848409750e-03, -1.9846597865511841e-16,
        2.1167581926038660e-02, -2.0404149472422147e-03, -9.6451162689070637e-03, -2.6947178926806793e-02, 9.2582118178574689e-04, 2.5975298860486740e-03,
        1.2278619011668618e-02, -1.4581758526841315e-04, -1.9953554848410001e-03, -2.5019234191788001e-16, -2.1167581926036835e-02, 2.0404149472422121e-03,
        -9.6451162689070880e-03, 2.6947178926807511e-02, 9.2582118178574755e-04, -2.5975298860486731e-03, 1.2278619011668634e-02, 1.4581758526841532e-04,
        1.9953554848409988e-03, 3.9159096504921086e-16, -2.1167581926036603e-02, 2.0404149472422251e-03, 9.6451162689071834e-03, 2.6947178926807359e-02,
        -9.2582118178574538e-04, -2.5975298860486736e-03, -1.2278619011668735e-02, 1.4581758526841545e-04, 1.9953554848409628e-03, 2.0929696023269232e-16,
        2.1167581926038861e-02, -2.0404149472422139e-03, 9.6451162689070169e-03, -2.6947178926807269e-02, -9.2582118178574809e-04, 2.5975298860486775e-03,
        -1.2278619011668728e-02, -1.4581758526841361e-04, -1.9953554848409823e-03, -7.9443263974794385e-16, 2.1167581926038424e-02, -2.0404149472422199e-03,
        -9.6451162689072650e-03, -2.6947178926807272e-02, 9.2582118178574646e-04, 2.5975298860486822e-03, 1.2278619011668740e-02, -1.4581758526841242e-04,
        -1.9953554848409828e-03, -7.8880633582714310e-16, -2.1167581926037807e-02, 2.0404149472422139e-03, -9.6451162689073638e-03, 2.6947178926807511e-02,
        9.2582118178574147e-04, -2.5975298860486788e-03, 1.2278619011668837e-02, 1.4581758526841556e-04, 1.9953554848409784e-03, 3.3757240821831445e-16,
        -2.1167581926038358e-02, 2.0404149472422026e-03, 9.6451162689071418e-03, 2.6947178926807695e-02, -9.2582118178574679e-04, -2.5975298860486827e-03,
        -1.2278619011668761e-02, 1.4581758526841396e-04, 1.9953554848410500e-03, 5.6669789601608066e-16, 1.2275170226804905e-02, -3.0219442678672110e-03,
        5.5947305550585781e-03, -1.5626782955866884e-02, -1.3712873229270483e-03, 3.8470559923940166e-03, -7.1223158999083576e-03, -9.0296763116215595e-05,
        -1.1559633592517281e-03, -6.9718374565947281e-17, 1.2275170226805658e-02, -3.0219442678672001e-03, -5.5947305550581826e-03, -1.5626782955866773e-02,
        1.3712873229270539e-03, 3.8470559923940275e-03, 7.1223158999083758e-03, -9.0296763116215229e-05, -1.1559633592517569e-03, -1.4733928253897693e-16,
        -1.2275170226803833e-02, 3.0219442678672010e-03, -5.5947305550581930e-03, 1.5626782955867457e-02, 1.3712873229270600e-03, -3.8470559923940400e-03,
        7.1223158999083949e-03, 9.0296763116214091e-05, 1.1559633592517521e-03, 2.5233815904687563e-16, -1.2275170226803531e-02, 3.0219442678672201e-03,
        5.5947305550582936e-03, 1.5626782955867318e-02, -1.3712873229270563e-03, -3.8470559923940353e-03, -7.1223158999084843e-03, 9.0296763116214240e-05,
        1.1559633592516992e-03, 5.2185638882926514e-17, 1.2275170226805847e-02, -3.0219442678672036e-03, 5.5947305550581201e-03, -1.5626782955867217e-02,
        -1.3712873229270600e-03, 3.8470559923940396e-03, -7.1223158999084608e-03, -9.0296763116213142e-05, -1.1559633592517328e-03, -6.6939086156592258e-16,
        1.2275170226805410e-02, -3.0219442678672010e-03, -5.5947305550583725e-03, -1.5626782955867214e-02, 1.3712873229270585e-03, 3.8470559923940357e-03,
        7.1223158999084730e-03, -9.0296763116213142e-05, -1.1559633592517333e-03, -6.4814241569201302e-16, -1.2275170226804797e-02, 3.0219442678672027e-03,
        -5.5947305550584931e-03, 1.5626782955867457e-02, 1.3712873229270509e-03, -3.8470559923940418e-03, 7.1223158999085727e-03, 9.0296763116214592e-05,
        1.1559633592517172e-03, 1.6762567593882249e-16, -1.2275170226805394e-02, 3.0219442678671832e-03, 5.5947305550582476e-03, 1.5626782955867613e-02,
        -1.3712873229270585e-03, -3.8470559923940518e-03, -7.1223158999085137e-03, 9.0296763116213643e-05, 1.1559633592517944e-03, 4.2371697733289753e-16,
        9.0474567801809424e-04, -3.4352259624794690e-03, 4.1169560207639477e-04, -1.1517774563953767e-03, -1.5590467568280778e-03, 4.3731801293314995e-03,
        -5.2410497766300346e-04, -6.9908714440816747e-06, -8.5716525570030737e-05, 1.0087562543266326e-16, 9.0474567801889904e-04, -3.4352259624794495e-03,
        -4.1169560207597003e-04, -1.1517774563952127e-03, 1.5590467568280886e-03, 4.3731801293315194e-03, 5.2410497766303328e-04, -6.9908714440807184e-06,
        -8.5716525570061461e-05, -5.6432601174820048e-18, -9.0474567801707064e-04, 3.4352259624794521e-03, -4.1169560207595583e-04, 1.1517774563958598e-03,
        1.5590467568280977e-03, -4.3731801293315281e-03, 5.2410497766305550e-04, 6.9908714440814486e-06, 8.5716525570062870e-05, 8.6255396812375980e-17,
        -9.0474567801669475e-04, 3.4352259624794807e-03, 4.1169560207607108e-04, 1.1517774563957334e-03, -1.5590467568280932e-03, -4.3731801293315237e-03,
        -5.2410497766312977e-04, 6.9908714440819949e-06, 8.5716525569984293e-05, -1.3917550946240810e-16, 9.0474567801907370e-04, -3.4352259624794599e-03,
        4.1169560207588384e-04, -1.1517774563956219e-03, -1.5590467568280992e-03, 4.3731801293315272e-03, -5.2410497766310223e-04, -6.9908714440815824e-06,
        -8.5716525570027485e-05, -5.1400324451132017e-16, 9.0474567801863742e-04, -3.4352259624794608e-03, -4.1169560207614204e-04, -1.1517774563956321e-03,
        1.5590467568280975e-03, 4.3731801293315263e-03, 5.2410497766311502e-04, -6.9908714440814571e-06, -8.5716525570032567e-05, -4.7577337964055999e-16,
        -9.0474567801801639e-04, 3.4352259624794634e-03, -4.1169560207629687e-04, 1.1517774563958468e-03, 1.5590467568280832e-03, -4.3731801293315333e-03,
        5.2410497766318864e-04, 6.9908714440820102e-06, 8.5716525570000217e-05, -2.7030317045506741e-17, -9.0474567801871299e-04, 3.4352259624794265e-03,
        4.1169560207601508e-04, 1.1517774563959680e-03, -1.5590467568280971e-03, -4.3731801293315437e-03, -5.2410497766315698e-04, 6.9908714440803449e-06,
        8.5716525570093269e-05, 2.6678103786357335e-16, -1.0581162695290412e-02, -3.0527070726671960e-03, -4.8324237754404486e-03, 1.3470243573400322e-02,
        -1.3858207373338267e-03, 3.8862182740439268e-03, 6.1518688616377165e-03, 1.0413219047989221e-04, 9.8886045271099243e-04, 2.6471475052799496e-16,
        -1.0581162695289520e-02, -3.0527070726671652e-03, 4.8324237754409204e-03, 1.3470243573400565e-02, 1.3858207373338431e-03, 3.8862182740439549e-03,
        -6.1518688616376757e-03, 1.0413219047989387e-04, 9.8886045271095969e-04, 1.3007774156867251e-16, 1.0581162695291366e-02, 3.0527070726671726e-03,
        4.8324237754409638e-03, -1.3470243573399944e-02, 1.3858207373338529e-03, -3.8862182740439606e-03, -6.1518688616376228e-03, -1.0413219047989234e-04,
        -9.8886045271095080e-04, -7.2029495873022837e-17, 1.0581162695291859e-02, 3.0527070726672134e-03, -4.8324237754408294e-03, -1.3470243573400048e-02,
        -1.3858207373338460e-03, -3.8862182740439519e-03, 6.1518688616375673e-03, -1.0413219047989091e-04, -9.8886045271106789e-04, -3.3870964976905638e-16,
        -1.0581162695289399e-02, -3.0527070726671869e-03, -4.8324237754410488e-03, 1.3470243573400178e-02, -1.3858207373338583e-03, 3.8862182740439588e-03,
        6.1518688616376107e-03, 1.0413219047989210e-04, 9.8886045271101564e-04, -3.6264231760319537e-16, -1.0581162695289836e-02, -3.0527070726671861e-03,
        4.8324237754407869e-03, 1.3470243573400162e-02, 1.3858207373338581e-03, 3.8862182740439571e-03, -6.1518688616375950e-03, 1.0413219047989166e-04,
        9.8886045271100783e-04, -3.0792945033231846e-16, 1.0581162695290471e-02, 3.0527070726671943e-03, 4.8324237754405726e-03, -1.3470243573399975e-02,
        1.3858207373338386e-03, -3.8862182740439627e-03, -6.1518688616375317e-03, -1.0413219047989114e-04, -9.8886045271105163e-04, -2.1358911569739497e-16,
        1.0581162695289621e-02, 3.0527070726671401e-03, -4.8324237754409187e-03, -1.3470243573399878e-02, -1.3858207373338620e-03, -3.8862182740439693e-03,
        6.1518688616375326e-03, -1.0413219047989298e-04, -9.8886045271093757e-04, 1.3214781814909572e-16, -1.9031365111332042e-02, -1.6465866379909038e-03,
        -8.7078302006787928e-03, 2.4227689429447936e-02, -7.4805055566576264e-04, 2.0961706872077266e-03, 1.1085416336256877e-02, 2.4312033959708541e-04,
        1.7660289489446149e-03, 3.8057582919529090e-16, -1.9031365111331047e-02, -1.6465866379908704e-03, 8.7078302006793409e-03, 2.4227689429448307e-02,
        7.4805055566578898e-04, 2.0961706872077699e-03, -1.1085416336256830e-02, 2.4312033959708639e-04, 1.7660289489445911e-03, 2.3881830595770688e-16,
        1.9031365111332910e-02, 1.6465866379908717e-03, 8.7078302006794259e-03, -2.4227689429447686e-02, 7.4805055566580069e-04, -2.0961706872077682e-03,
        -1.1085416336256692e-02, -2.4312033959708937e-04, -1.7660289489445677e-03, -1.8584489782554045e-16, 1.9031365111333579e-02, 1.6465866379909394e-03,
        -8.7078302006792629e-03, -2.4227689429447780e-02, -7.4805055566579137e-04, -2.0961706872077665e-03, 1.1085416336256659e-02, -2.4312033959708498e-04,
        -1.7660289489447455e-03, -5.1718202387898003e-16, -1.9031365111331008e-02, -1.6465866379908984e-03, -8.7078302006795335e-03, 2.4227689429447918e-02,
        -7.4805055566581186e-04, 2.0961706872077656e-03, 1.1085416336256747e-02, 2.4312033959708631e-04, 1.7660289489446780e-03, -2.4150304858849214e-16,
        -1.9031365111331439e-02, -1.6465866379908991e-03, 8.7078302006792785e-03, 2.4227689429447880e-02, 7.4805055566581598e-04, 2.0961706872077573e-03,
        -1.1085416336256737e-02, 2.4312033959708406e-04, 1.7660289489446750e-03, -1.8347290836793516e-16, 1.9031365111332136e-02, 1.6465866379909175e-03,
        8.7078302006789923e-03, -2.4227689429447741e-02, 7.4805055566579007e-04, -2.0961706872077812e-03, -1.1085416336256640e-02, -2.4312033959708709e-04,
        -1.7660289489446971e-03, -3.4726000200140564e-16, 1.9031365111331040e-02, 1.6465866379908351e-03, -8.7078302006794485e-03, -2.4227689429447651e-02,
        -7.4805055566582986e-04, -2.0961706872077816e-03, 1.1085416336256615e-02, -2.4312033959708875e-04, -1.7660289489445636e-03, 6.9267354347634622e-17,
        -2.0505104876293125e-02, 1.0112205753418885e-03, -9.4250538894833577e-03, 2.6103819129888628e-02, 4.5837908688127320e-04, -1.2873242618557583e-03,
        1.1998470795679696e-02, 4.1004436283959206e-04, 1.8695463276148204e-03, 4.3637760805796555e-16, -2.0505104876292011e-02, 1.0112205753419238e-03,
        9.4250538894840395e-03, 2.6103819129889190e-02, -4.5837908688122088e-04, -1.2873242618556893e-03, -1.1998470795679674e-02, 4.1004436283959467e-04,
        1.8695463276148486e-03, 2.9371309167686647e-16, 2.0505104876293870e-02, -1.0112205753419258e-03, 9.4250538894841800e-03, -2.6103819129888555e-02,
        -4.5837908688119828e-04, 1.2873242618556939e-03, -1.1998470795679395e-02, -4.1004436283959342e-04, -1.8695463276147601e-03, -2.1445411480602755e-16,
        2.0505104876294828e-02, -1.0112205753418297e-03, -9.4250538894839822e-03, -2.6103819129888618e-02, 4.5837908688120944e-04, 1.2873242618557062e-03,
        1.1998470795679393e-02, -4.1004436283959196e-04, -1.8695463276150368e-03, -6.4377669631513899e-16, -2.0505104876292091e-02, 1.0112205753418824e-03,
        -9.4250538894843239e-03, 2.6103819129888756e-02, 4.5837908688119063e-04, -1.2873242618557058e-03, 1.1998470795679535e-02, 4.1004436283959304e-04,
        1.8695463276149720e-03, -1.9183169463589791e-16, -2.0505104876292531e-02, 1.0112205753418854e-03, 9.4250538894840898e-03, 2.6103819129888715e-02,
        -4.5837908688118191e-04, -1.2873242618557073e-03, -1.1998470795679519e-02, 4.1004436283959380e-04, 1.8695463276149843e-03, -1.4415816212542290e-16,
        2.0505104876293291e-02, -1.0112205753418581e-03, 9.4250538894836821e-03, -2.6103819129888638e-02, -4.5837908688122555e-04, 1.2873242618556954e-03,
        -1.1998470795679428e-02, -4.1004436283959196e-04, -1.8695463276149401e-03, -4.2175518614235540e-16, 2.0505104876291865e-02, -1.0112205753419709e-03,
        -9.4250538894843309e-03, -2.6103819129888559e-02, 4.5837908688115355e-04, 1.2873242618556950e-03, 1.1998470795679389e-02, -4.1004436283959678e-04,
        -1.8695463276148037e-03, 1.0082365576140328e-16, -1.0271215618252966e-02, 5.1490211242928427e-03, -4.8319444906675807e-03, 1.3075668540126886e-02,
        2.3387746451529658e-03, -6.5549099570769986e-03, 6.1512587129415961e-03, 6.0501265406261173e-04, 8.5061915836237459e-04, 4.1782997179837206e-16,
        -1.0271215618251738e-02, 5.1490211242928722e-03, 4.8319444906685131e-03, 1.3075668540127792e-02, -2.3387746451528726e-03, -6.5549099570768746e-03,
        -6.1512587129416560e-03, 6.0501265406261542e-04, 8.5061915836254654e-04, 2.6693706628250703e-16, 1.0271215618253567e-02, -5.1490211242928869e-03,
        4.8319444906687256e-03, -1.3075668540127089e-02, -2.3387746451528505e-03, 6.5549099570769084e-03, -6.1512587129411147e-03, -6.0501265406261390e-04,
        -8.5061915836233762e-04, -1.7021528662009007e-16, 1.0271215618254951e-02, -5.1490211242927377e-03, -4.8319444906684767e-03, -1.3075668540127100e-02,
        2.3387746451528696e-03, 6.5549099570769301e-03, 6.1512587129411772e-03, -6.0501265406261130e-04, -8.5061915836276327e-04, -7.1891925493360827e-16,
        -1.0271215618251993e-02, 5.1490211242928132e-03, -4.8319444906688921e-03, 1.3075668540127258e-02, 2.3387746451528453e-03, -6.5549099570769154e-03,
        6.1512587129413984e-03, 6.0501265406261216e-04, 8.5061915836275167e-04, -2.2397369512977885e-16, -1.0271215618252423e-02, 5.1490211242928097e-03,
        4.8319444906687308e-03, 1.3075668540127205e-02, -2.3387746451528127e-03, -6.5549099570769162e-03, -6.1512587129413689e-03, 6.0501265406261205e-04,
        8.5061915836281857e-04, -1.9582901219255179e-16, 1.0271215618253280e-02, -5.1490211242927802e-03, 4.8319444906680890e-03, -1.3075668540127176e-02,
        -2.3387746451529043e-03, 6.5549099570769015e-03, -6.1512587129412370e-03, -6.0501265406260999e-04, -8.5061915836256996e-04, -4.2037013883311255e-16,
        1.0271215618251405e-02, -5.1490211242929294e-03, -4.8319444906691081e-03, -1.3075668540127058e-02, 2.3387746451527603e-03, 6.5549099570769223e-03,
        6.1512587129412049e-03, -6.0501265406261661e-04, -8.5061915836252063e-04, 2.5108205205904367e-16,
  };
  return sol;
}

static const double *
get_sol_3x_biasB_p1_DxDx_DyDy(void)
{
  static const double sol[512] = {
        -4.4041922402994508e-03, -2.5427615754988948e-03, -1.0665671304152860e-04,  3.7718620629148253e-05,
        -6.1578281985007731e-05,  2.1776855773916784e-05, -1.6454514526474550e-05, -9.5000183910323872e-06,
        -2.6096743371754640e-03, -1.5066961810654262e-03, -6.9169901643607081e-05,  6.1409551253443335e-04,
        -3.9935261333768437e-05,  3.5454820946981773e-04,  1.6276715296146205e-05,  9.3973659577225150e-06,
        -4.8238508035630792e-04, -2.7850515599675828e-04, -1.2785705897272424e-05,  6.1409551253432623e-04,
        -7.3818307414892410e-06,  3.5454820946985161e-04,  1.6276715296142577e-05,  9.3973659577230299e-06,
         1.6449041764628515e-03,  9.4968586907191971e-04,  4.3598489849228317e-05,  6.1409551253429771e-04,
         2.5171599850703604e-05,  3.5454820946987070e-04,  1.6276715296122763e-05,  9.3973659577310852e-06,
         3.7721934332819408e-03,  2.1778768941405919e-03,  9.9982685595574256e-05,  6.1409551253416305e-04,
         5.7725030442970594e-05,  3.5454820946990772e-04,  1.6276715296088936e-05,  9.3973659577540534e-06,
         5.8994826901014163e-03,  3.4060679192091698e-03,  1.5636688134197224e-04,  6.1409551253423353e-04,
         9.0278461035230045e-05,  3.5454820946988067e-04,  1.6276715296172321e-05,  9.3973659576982136e-06,
         8.0267719469204173e-03,  4.6342589442779272e-03,  2.1275107708848732e-04,  6.1409551253416761e-04,
         1.2283189162741587e-04,  3.5454820946994730e-04,  1.6276715296031209e-05,  9.3973659577826543e-06,
         1.0154061203738809e-02,  5.8624499693468073e-03,  2.6913527283490597e-04,  6.1409551253420447e-04,
         1.5538532221965102e-04,  3.5454820946989601e-04,  1.6276715296035092e-05,  9.3973659577696252e-06,
                                                         
        -1.1018979264020300e-02, -1.2762874932783809e-03, -2.6494102141839070e-04, -2.8994604615678344e-04,
        -2.9807206064867279e-05, -2.1095413934670354e-04, -5.7460710141012667e-05, -1.4174919685463714e-05,
        -6.6613512061697421e-03, -8.3254054991779977e-04, -1.7752438980077028e-04,  1.5675158485682092e-03,
        -2.2623231571673787e-05,  1.9590927819014637e-04,  4.1774151505839880e-05,  5.3235857015895280e-06,
        -1.2313170233906152e-03, -1.5389090291880844e-04, -3.2814484098407900e-05,  1.5675158485681827e-03,
        -4.1817897444755747e-06,  1.9590927819015892e-04,  4.1774151505836052e-05,  5.3235857015888046e-06,
         4.1987171593885252e-03,  5.2475874408017883e-04,  1.1189542160396868e-04,  1.5675158485681906e-03,
         1.4259652082721208e-05,  1.9590927819016068e-04,  4.1774151505830529e-05,  5.3235857015890257e-06,
         9.6287513421676133e-03,  1.2034083910791829e-03,  2.5660532730632313e-04,  1.5675158485681483e-03,
         3.2701093909920751e-05,  1.9590927819017681e-04,  4.1774151505832115e-05,  5.3235857015865329e-06,
         1.5058785524946838e-02,  1.8820580380781381e-03,  4.0131523300871499e-04,  1.5675158485681689e-03,
         5.1142535737119170e-05,  1.9590927819017611e-04,  4.1774151505820053e-05,  5.3235857015869310e-06,
         2.0488819707725963e-02,  2.5607076850771271e-03,  5.4602513871109741e-04,  1.5675158485682298e-03,
         6.9583977564318124e-05,  1.9590927819018239e-04,  4.1774151505813304e-05,  5.3235857015801149e-06,
         2.5918853890504821e-02,  3.2393573320761929e-03,  6.9073504441345018e-04,  1.5675158485681817e-03,
         8.8025419391506493e-05,  1.9590927819018424e-04,  4.1774151505804312e-05,  5.3235857015857799e-06,
                                                         
        -1.3366720692763396e-02, -7.9181652594083197e-05, -3.2014002221508991e-04, -1.5665805187687065e-03,
        -2.0619519041054335e-06, -5.2611111707254622e-04, -1.0875897838019112e-04, -1.5442149291387570e-05,
        -8.4983155591709623e-03, -2.2803131377920161e-04, -2.2803131377920622e-04,  1.9997811123959363e-03,
        -6.5369545832163292e-06,  5.3659188242116111e-05,  5.3659188242125605e-05,  1.5382434574337128e-06,
        -1.5708705777981867e-03, -4.2150433122697153e-05, -4.2150433122704234e-05,  1.9997811123959380e-03,
        -1.2083229378424886e-06,  5.3659188242120123e-05,  5.3659188242120949e-05,  1.5382434574338455e-06,
         5.3565744035745932e-03,  1.4373044753380644e-04,  1.4373044753380633e-04,  1.9997811123959493e-03,
         4.1203087075293580e-06,  5.3659188242120502e-05,  5.3659188242116450e-05,  1.5382434574343029e-06,
         1.2284019384947355e-02,  3.2961132819031189e-04,  3.2961132819030463e-04,  1.9997811123959436e-03,
         9.4489403529038885e-06,  5.3659188242125564e-05,  5.3659188242113645e-05,  1.5382434574343915e-06,
         1.9211464366320145e-02,  5.1549220884680247e-04,  5.1549220884683142e-04,  1.9997811123959662e-03,
         1.4777571998274944e-05,  5.3659188242126635e-05,  5.3659188242104951e-05,  1.5382434574354242e-06,
         2.6138909347692911e-02,  7.0137308950330574e-04,  7.0137308950334932e-04,  1.9997811123960299e-03,
         2.0106203643644923e-05,  5.3659188242122142e-05,  5.3659188242088599e-05,  1.5382434574373129e-06,
         3.3066354329065621e-02,  8.8725397015985226e-04,  8.8725397015982602e-04,  1.9997811123959944e-03,
         2.5434835289019607e-05,  5.3659188242128064e-05,  5.3659188242090029e-05,  1.5382434574371947e-06,
                                                         
        -1.1772420545074049e-02,  9.9965127203158893e-04, -2.9771347337083224e-04, -4.1374953013227720e-03,
         1.5009925916331957e-05, -9.5820722469863165e-04, -1.4619029680330181e-04, -6.1688324763178938e-06,
        -8.4328227495495626e-03,  2.6584360504410145e-04, -2.2754565363904947e-04,  1.9843696720032210e-03,
         6.8173505958705989e-06, -6.2556987501449985e-05,  5.3544905127008731e-05, -1.6042248447103262e-06,
        -1.5587645637327799e-03,  4.9139843602173047e-05, -4.2060661306196192e-05,  1.9843696720032301e-03,
         1.2601527202690473e-06, -6.2556987501449782e-05,  5.3544905127005309e-05, -1.6042248447096918e-06,
         5.3152936220840080e-03, -1.6756391783975450e-04,  1.4342433102665958e-04,  1.9843696720032410e-03,
        -4.2970451553338836e-06, -6.2556987501450744e-05,  5.3544905127002904e-05, -1.6042248447089394e-06,
         1.2189351807900777e-02, -3.8426767928168319e-04,  3.2890932335951224e-04,  1.9843696720032457e-03,
        -9.8542430309342842e-06, -6.2556987501449619e-05,  5.3544905127001169e-05, -1.6042248447085263e-06,
         1.9063409993717553e-02, -6.0097144072361120e-04,  5.1439431569237854e-04,  1.9843696720032670e-03,
        -1.5411440906539592e-05, -6.2556987501451611e-05,  5.3544905126996466e-05, -1.6042248447072574e-06,
         2.5937468179534327e-02, -8.1767520216553769e-04,  6.9987930802523882e-04,  1.9843696720033108e-03,
        -2.0968638782141884e-05, -6.2556987501459797e-05,  5.3544905126988287e-05, -1.6042248447047786e-06,
         3.2811526365351132e-02, -1.0343789636074571e-03,  8.8536430035807004e-04,  1.9843696720032926e-03,
        -2.6525836657742338e-05, -6.2556987501455094e-05,  5.3544905126988287e-05, -1.6042248447056049e-06,
                                                         
        -9.0272890314765722e-03,  5.8525114630483929e-04, -2.2574513915869974e-04, -4.4933715731932363e-03,
         2.6541011214171422e-05,  7.5274196333601873e-04, -1.5552020407559827e-04,  7.8220800114296213e-07,
        -6.9204828041277604e-03,  6.0730626955139981e-04, -1.8789706627799494e-04,  1.6284934001327706e-03,
         1.6073771990022628e-05, -1.4290827386115464e-04,  4.4214997854705820e-05, -3.7823996304681823e-06,
        -1.2792161864866225e-03,  1.1225748725241897e-04, -3.4731820796203245e-05,  1.6284934001327780e-03,
         2.9711553210232943e-06, -1.4290827386115610e-04,  4.4214997854704085e-05, -3.7823996304676809e-06,
         4.3620504311545211e-03, -3.8279129504656216e-04,  1.1843342468558873e-04,  1.6284934001327851e-03,
        -1.0131461347976333e-05, -1.4290827386115716e-04,  4.4214997854703523e-05, -3.7823996304675407e-06,
         1.0003317048795649e-02, -8.7784007734554280e-04,  2.7159867016738124e-04,  1.6284934001327919e-03,
        -2.3234078016975833e-05, -1.4290827386115735e-04,  4.4214997854702446e-05, -3.7823996304672675e-06,
         1.5644583666436781e-02, -1.3728888596445214e-03,  4.2476391564917615e-04,  1.6284934001328049e-03,
        -3.6336694685976650e-05, -1.4290827386115966e-04,  4.4214997854700812e-05, -3.7823996304667661e-06,
         2.1285850284077916e-02, -1.8679376419435008e-03,  5.7792916113097177e-04,  1.6284934001328227e-03,
        -4.9439311354977518e-05, -1.4290827386116670e-04,  4.4214997854699179e-05, -3.7823996304657628e-06,
         2.6927116901719089e-02, -2.3629864242424841e-03,  7.3109440661274961e-04,  1.6284934001328183e-03,
        -6.2541928023973379e-05, -1.4290827386116294e-04,  4.4214997854697953e-05, -3.7823996304655858e-06,
                                                         
        -6.2247900453550730e-03,  1.0327723977360613e-03, -1.3371278969370395e-04, -2.4924537019540933e-03,
         2.6593890523597406e-05,  4.0248850825023367e-04, -1.3292722975764579e-04,  1.2261851803121024e-05,
        -4.5637037010479806e-03,  7.5338078003219552e-04, -1.2532535988959181e-04,  1.0739079292105612e-03,
         2.0052019536976041e-05, -1.7728179706443800e-04,  2.9490936864655919e-05, -4.7185409457015368e-06,
        -8.4357750896041491e-04,  1.3925862048674964e-04, -2.3165757864813873e-05,  1.0739079292105657e-03,
         3.7065142258784086e-06, -1.7728179706443867e-04,  2.9490936864655435e-05, -4.7185409457013454e-06,
         2.8765486831271551e-03, -4.7486353905869716e-04,  7.8993844159963767e-05,  1.0739079292105696e-03,
        -1.2638991085219166e-05, -1.7728179706443922e-04,  2.9490936864655282e-05, -4.7185409457013344e-06,
         6.5966748752147126e-03, -1.0889856986041411e-03,  1.8115344618474199e-04,  1.0739079292105755e-03,
        -2.8984496396316930e-05, -1.7728179706443943e-04,  2.9490936864654591e-05, -4.7185409457010650e-06,
         1.0316801067302282e-02, -1.7031078581495856e-03,  2.8331304820951964e-04,  1.0739079292105826e-03,
        -4.5330001707415101e-05, -1.7728179706444065e-04,  2.9490936864654490e-05, -4.7185409457012421e-06,
         1.4036927259389841e-02, -2.3172300176950325e-03,  3.8547265023429852e-04,  1.0739079292105855e-03,
        -6.1675507018512700e-05, -1.7728179706444239e-04,  2.9490936864654900e-05, -4.7185409457009464e-06,
         1.7757053451477430e-02, -2.9313521772404848e-03,  4.8763225225907262e-04,  1.0739079292105883e-03,
        -7.8021012329608367e-05, -1.7728179706444209e-04,  2.9490936864654080e-05, -4.7185409457005339e-06,
                                                         
        -2.7685660836065627e-03,  9.6267943695907572e-04, -5.1302926806916285e-05, -1.3595218641067415e-03,
         2.0985465997969670e-05,  2.5160999330443233e-04, -8.5156552130073984e-05,  1.5318561784528383e-05,
        -2.1160573442293054e-03,  6.5976850295807366e-04, -5.9827506156292608e-05,  4.9794003061817010e-04,
         1.7763183947286998e-05, -1.5525342423246231e-04,  1.4078309516760695e-05, -4.1799436025251460e-06,
        -3.9114247992314439e-04,  1.2195486531873465e-04, -1.1058811420875967e-05,  4.9794003061817314e-04,
         3.2834345625935113e-06, -1.5525342423246226e-04,  1.4078309516760042e-05, -4.1799436025254264e-06,
         1.3337723843830172e-03, -4.1585877232060533e-04,  3.7709883314540850e-05,  4.9794003061817682e-04,
        -1.1196314822099668e-05, -1.5525342423246185e-04,  1.4078309516759505e-05, -4.1799436025257211e-06,
         3.0586872486891748e-03, -9.5367240995994330e-04,  8.6478578049958005e-05,  4.9794003061818224e-04,
        -2.5676064206792951e-05, -1.5525342423246198e-04,  1.4078309516758496e-05, -4.1799436025261192e-06,
         4.7836021129953393e-03, -1.4914860475992832e-03,  1.3524727278537308e-04,  4.9794003061818745e-04,
        -4.0155813591486747e-05, -1.5525342423246160e-04,  1.4078309516757525e-05, -4.1799436025265182e-06,
         6.5085169773014956e-03, -2.0292996852386224e-03,  1.8401596752079196e-04,  4.9794003061818853e-04,
        -5.4635562976179490e-05, -1.5525342423246120e-04,  1.4078309516757781e-05, -4.1799436025266952e-06,
         8.2334318416076614e-03, -2.5671133228779642e-03,  2.3278466225620765e-04,  4.9794003061819168e-04,
        -6.9115312360872741e-05, -1.5525342423246101e-04,  1.4078309516756759e-05, -4.1799436025271077e-06,
                                                         
        -5.5057819369581863e-04,  3.1787646834022232e-04, -7.4775167389616631e-06, -4.6186028600569848e-04,
         4.3171463021096952e-06,  2.6665516045338905e-04, -2.9312012410193478e-05,  1.6923298255515372e-05,
        -4.8665238793613104e-04,  2.8096888717670101e-04, -1.4530384527701296e-05,  1.1451660589924915e-04,
         8.3891214184953722e-06, -6.6116193242614837e-05,  3.4192174122083424e-06, -1.9740860933551212e-06,
        -8.9955228480437954e-05,  5.1935675378195315e-05, -2.6858679675683887e-06,  1.1451660589925117e-04,
         1.5506865940831339e-06, -6.6116193242615583e-05,  3.4192174122070139e-06, -1.9740860933552669e-06,
         3.0674193097525530e-04, -1.7709753642030946e-04,  9.1586485925635640e-06,  1.1451660589925462e-04,
        -5.2877482303300981e-06, -6.6116193242616098e-05,  3.4192174122056849e-06, -1.9740860933554109e-06,
         7.0343909043094825e-04, -4.0613074821881383e-04,  2.1003165152696916e-05,  1.1451660589925732e-04,
        -1.2126183054742609e-05, -6.6116193242617480e-05,  3.4192174122050653e-06, -1.9740860933548040e-06,
         1.1001362498866394e-03, -6.3516396001732215e-04,  3.2847681712829456e-05,  1.1451660589926117e-04,
        -1.8964617879153966e-05, -6.6116193242618767e-05,  3.4192174122030078e-06, -1.9740860933549811e-06,
         1.4968334093423308e-03, -8.6419717181582527e-04,  4.4692198272961002e-05,  1.1451660589926484e-04,
        -2.5803052703569116e-05, -6.6116193242617331e-05,  3.4192174122015514e-06, -1.9740860933558666e-06,
         1.8935305687980254e-03, -1.0932303836143303e-03,  5.6536714833096980e-05,  1.1451660589926414e-04,
        -3.2641487527979211e-05, -6.6116193242619960e-05,  3.4192174122038129e-06, -1.9740860933535500e-06,
  };
  return sol;
}

static const double *
get_sol_3x_biasB_p1_DxDx_PyPy(void)
{
  static const double sol[512] = {
        -4.4584195633810733e-02, -2.5740697350783534e-02, -2.5740697350783579e-02, -1.8124746495123774e-03,
        -1.4861398544603536e-02, -1.0464327267955061e-03, -1.0464327267952910e-03, -6.0415821650423960e-04,
        -2.7172476676301741e-02, -1.5688036723611624e-02, -1.5688036723611940e-02,  6.3940913062058213e-03,
        -9.0574922254338681e-03,  3.6916303368608886e-03,  3.6916303368611492e-03,  2.1313637687352016e-03,
        -5.0226946551364480e-03, -2.8998541112002685e-03, -2.8998541112003865e-03,  6.3940913062057042e-03,
        -1.6742315517121146e-03,  3.6916303368608925e-03,  3.6916303368609532e-03,  2.1313637687352198e-03,
         1.7127087366028894e-02,  9.8883285012110831e-03,  9.8883285012111924e-03,  6.3940913062055550e-03,
         5.7090291220096342e-03,  3.6916303368608968e-03,  3.6916303368608166e-03,  2.1313637687352294e-03,
         3.9276869387194165e-02,  2.2676511113622431e-02,  2.2676511113622465e-02,  6.3940913062055810e-03,
         1.3092289795731417e-02,  3.6916303368608795e-03,  3.6916303368607927e-03,  2.1313637687352285e-03,
         6.1426651408359283e-02,  3.5464693726033808e-02,  3.5464693726033690e-02,  6.3940913062052454e-03,
         2.0475550469453228e-02,  3.6916303368609051e-03,  3.6916303368603937e-03,  2.1313637687352597e-03,
         8.3576433429524880e-02,  4.8252876338445126e-02,  4.8252876338445223e-02,  6.3940913062050780e-03,
         2.7858811143174983e-02,  3.6916303368609116e-03,  3.6916303368605278e-03,  2.1313637687352407e-03,
         1.0572621545068980e-01,  6.1041058950856514e-02,  6.1041058950856819e-02,  6.3940913062048828e-03,
         3.5242071816896710e-02,  3.6916303368609302e-03,  3.6916303368602554e-03,  2.1313637687352754e-03,
                                                          
        -1.1102433537931264e-01, -1.2618535216278336e-02, -6.4099929917845411e-02, -8.9682417343014686e-03,
        -7.2853147038970641e-03, -3.0849513258657400e-03, -5.1778167794565748e-03, -1.7810974784254808e-03,
        -6.9196949126305612e-02, -8.5748037579501200e-03, -3.9950877205173575e-02,  1.6283079974463717e-02,
        -4.9506652579006920e-03,  2.0177799327709633e-03,  9.4010406064929421e-03,  1.1649657873507272e-03,
        -1.2790705487349173e-02, -1.5850090338456765e-03, -7.3847172562463191e-03,  1.6283079974463602e-02,
        -9.1510539235878645e-04,  2.0177799327709620e-03,  9.4010406064927929e-03,  1.1649657873507361e-03,
         4.3615538151607305e-02,  5.4047856902587656e-03,  2.5181442692680946e-02,  1.6283079974463464e-02,
         3.1204544731831146e-03,  2.0177799327709646e-03,  9.4010406064926819e-03,  1.1649657873507419e-03,
         1.0002178179056370e-01,  1.2394580414363204e-02,  5.7747602641608005e-02,  1.6283079974463429e-02,
         7.1560143387250389e-03,  2.0177799327709498e-03,  9.4010406064926472e-03,  1.1649657873507367e-03,
         1.5642802542952000e-01,  1.9384375138467660e-02,  9.0313762590535099e-02,  1.6283079974463169e-02,
         1.1191574204266976e-02,  2.0177799327709732e-03,  9.4010406064923436e-03,  1.1649657873507619e-03,
         2.1283426906847663e-01,  2.6374169862572072e-02,  1.2287992253946231e-01,  1.6283079974463037e-02,
         1.5227134069808879e-02,  2.0177799327709767e-03,  9.4010406064924355e-03,  1.1649657873507486e-03,
         2.6924051270743282e-01,  3.3363964586676544e-02,  1.5544608248838956e-01,  1.6283079974462881e-02,
         1.9262693935350784e-02,  2.0177799327709745e-03,  9.4010406064922517e-03,  1.1649657873507658e-03,
                                                          
        -1.3411820617565154e-01, -7.1471730461837319e-04, -7.7433182438742088e-02, -2.3541399004220864e-02,
        -4.1264222821590300e-04, -5.3288649468649753e-03, -1.3591633052187288e-02, -3.0766216115476501e-03,
        -8.7959388364324778e-02, -2.2576955867742456e-03, -5.0783376549897868e-02,  2.0698163334151001e-02,
        -1.3034811547723213e-03,  5.3126963343912994e-04,  1.1950090172702994e-02,  3.0672866587834431e-04,
        -1.6258847328107777e-02, -4.1732359150410962e-04, -9.3870498815960927e-03,  2.0698163334150890e-02,
        -2.4094188789407537e-04,  5.3126963343912972e-04,  1.1950090172702871e-02,  3.0672866587835244e-04,
         5.5441693708109245e-02,  1.4230484037660233e-03,  3.2009276786705683e-02,  2.0698163334150758e-02,
         8.2159737898417041e-04,  5.3126963343913048e-04,  1.1950090172702781e-02,  3.0672866587835640e-04,
         1.2714223474432620e-01,  3.2634203990361591e-03,  7.3405603455007323e-02,  2.0698163334150685e-02,
         1.8841366458624319e-03,  5.3126963343912452e-04,  1.1950090172702732e-02,  3.0672866587835499e-04,
         1.9884277578054310e-01,  5.1037923943063040e-03,  1.1480193012330903e-01,  2.0698163334150494e-02,
         2.9466759127406974e-03,  5.3126963343914316e-04,  1.1950090172702516e-02,  3.0672866587837911e-04,
         2.7054331681676019e-01,  6.9441643895764137e-03,  1.5619825679161076e-01,  2.0698163334150376e-02,
         4.0092151796189430e-03,  5.3126963343914295e-04,  1.1950090172702569e-02,  3.0672866587836512e-04,
         3.4224385785297706e-01,  8.7845363848465582e-03,  1.9759458345991249e-01,  2.0698163334150217e-02,
         5.0717544464971807e-03,  5.3126963343913763e-04,  1.1950090172702404e-02,  3.0672866587837380e-04,
                                                          
        -1.1928578703942408e-01,  9.2782184856526548e-03, -6.8869681257707843e-02, -4.6181363305459550e-02,
         5.3567819402917122e-03, -7.7423245368986879e-03, -2.6662822535950936e-02, -4.4700331555318972e-03,
        -8.6910597817874607e-02,  2.8632150910907813e-03, -5.0177857045581288e-02,  2.0451367188368066e-02,
         1.6530780035890727e-03, -6.7375745464186710e-04,  1.1807602351500225e-02, -3.8899404780600724e-04,
        -1.6064983708874006e-02,  5.2925080425480250e-04, -9.2751226688453908e-03,  2.0451367188367962e-02,
         3.0556309430533808e-04, -6.7375745464186331e-04,  1.1807602351500130e-02, -3.8899404780599808e-04,
         5.4780630400126623e-02, -1.8047134825811763e-03,  3.1627611707890503e-02,  2.0451367188367844e-02,
        -1.0419518149783942e-03, -6.7375745464185745e-04,  1.1807602351500053e-02, -3.8899404780599412e-04,
         1.2562624450912718e-01, -4.1386777694171518e-03,  7.2530346084626315e-02,  2.0451367188367760e-02,
        -2.3894667242621149e-03, -6.7375745464186157e-04,  1.1807602351499998e-02, -3.8899404780599623e-04,
         1.9647185861812774e-01, -6.4726420562531166e-03,  1.1343308046136219e-01,  2.0451367188367618e-02,
        -3.7369816335458337e-03, -6.7375745464185604e-04,  1.1807602351499856e-02, -3.8899404780597677e-04,
         2.6731747272712841e-01, -8.8066063430890917e-03,  1.5433581483809802e-01,  2.0451367188367504e-02,
        -5.0844965428295582e-03, -6.7375745464184748e-04,  1.1807602351499857e-02, -3.8899404780599114e-04,
         3.3816308683612889e-01, -1.1140570629925068e-02,  1.9523854921483388e-01,  2.0451367188367358e-02,
        -6.4320114521133013e-03, -6.7375745464184401e-04,  1.1807602351499738e-02, -3.8899404780597807e-04,
                                                          
        -9.0256922639355669e-02,  7.4816041899959870e-03, -5.2109858582059222e-02, -4.9944638029515029e-02,
         4.3195061930644124e-03,  5.5695968619307399e-03, -2.8835550210918881e-02,  3.2156082475134026e-03,
        -7.0918099467613696e-02,  6.3700581364470044e-03, -4.0944583818043473e-02,  1.6688092464312621e-02,
         3.6777547798312509e-03, -1.4989702203260580e-03,  9.6348746765322780e-03, -8.6543086021248668e-04,
        -1.3108851408420578e-02,  1.1774729751721309e-03, -7.5683988894184482e-03,  1.6688092464312531e-02,
         6.7981433917914131e-04, -1.4989702203260537e-03,  9.6348746765322121e-03, -8.6543086021248039e-04,
         4.4700396650772564e-02, -4.0151121861027438e-03,  2.5807786039206568e-02,  1.6688092464312434e-02,
        -2.3181261014729678e-03, -1.4989702203260457e-03,  9.6348746765321479e-03, -8.6543086021247551e-04,
         1.0250964470996564e-01, -9.2076973473776153e-03,  5.9183970967831553e-02,  1.6688092464312344e-02,
        -5.3160665421250692e-03, -1.4989702203260431e-03,  9.6348746765320993e-03, -8.6543086021247031e-04,
         1.6031889276915876e-01, -1.4400282508652482e-02,  9.2560155896456567e-02,  1.6688092464312219e-02,
        -8.3140069827771709e-03, -1.4989702203260405e-03,  9.6348746765320039e-03, -8.6543086021246185e-04,
         2.1812814082835194e-01, -1.9592867669927353e-02,  1.2593634082508157e-01,  1.6688092464312139e-02,
        -1.1311947423429281e-02, -1.4989702203260298e-03,  9.6348746765319727e-03, -8.6543086021246326e-04,
         2.7593738888754493e-01, -2.4785452831202218e-02,  1.5931252575370652e-01,  1.6688092464312014e-02,
        -1.4309887864081385e-02, -1.4989702203260164e-03,  9.6348746765319068e-03, -8.6543086021245220e-04,
                                                          
        -5.8558107995076741e-02,  1.0819714977870510e-02, -3.3808539414192727e-02, -3.3336939184891466e-02,
         6.2467653550285626e-03,  4.0188625366329098e-03, -1.9247090812355247e-02,  2.3202913673611098e-03,
        -4.6332036744757193e-02,  7.8247117949072180e-03, -2.6749813886689226e-02,  1.0902623153480511e-02,
         4.5175994611209259e-03, -1.8412720436711764e-03,  6.2946324125350284e-03, -1.0630589100648897e-03,
        -8.5642422695474067e-03,  1.4463583345179297e-03, -4.9445675797283829e-03,  1.0902623153480446e-02,
         8.3505537377858677e-04, -1.8412720436711710e-03,  6.2946324125349841e-03, -1.0630589100648835e-03,
         2.9203552205662404e-02, -4.9319951258713617e-03,  1.6860678727232462e-02,  1.0902623153480373e-02,
        -2.8474887135637528e-03, -1.8412720436711636e-03,  6.2946324125349408e-03, -1.0630589100648780e-03,
         6.6971346680872179e-02, -1.1310348586260647e-02,  3.8665925034193284e-02,  1.0902623153480289e-02,
        -6.5300328009060870e-03, -1.8412720436711597e-03,  6.2946324125349000e-03, -1.0630589100648752e-03,
         1.0473914115608197e-01, -1.7688702046649935e-02,  6.0471171341154130e-02,  1.0902623153480199e-02,
        -1.0212576888248429e-02, -1.8412720436711517e-03,  6.2946324125348506e-03, -1.0630589100648624e-03,
         1.4250693563129183e-01, -2.4067055507039227e-02,  8.2276417648114969e-02,  1.0902623153480149e-02,
        -1.3895120975590763e-02, -1.8412720436711443e-03,  6.2946324125348081e-03, -1.0630589100648648e-03,
         1.8027473010650152e-01, -3.0445408967428502e-02,  1.0408166395507577e-01,  1.0902623153480041e-02,
        -1.7577665062933099e-02, -1.8412720436711389e-03,  6.2946324125347829e-03, -1.0630589100648566e-03,
                                                          
        -2.3652444410766166e-02,  9.3330792887737161e-03, -1.3655745147548512e-02, -2.0295051641274243e-02,
         5.3884558397416238e-03,  3.5108747474153165e-03, -1.1717353528307028e-02,  2.0270044805112890e-03,
        -2.1062278650861084e-02,  6.7647898429604248e-03, -1.2160312248821570e-02,  4.9562700674910870e-03,
         3.9056532368444519e-03, -1.5918565112213788e-03,  2.8615038576424623e-03, -9.1905878526492464e-04,
        -3.8932555050065249e-03,  1.2504371313709281e-03, -2.2477721138395197e-03,  4.9562700674910401e-03,
         7.2194021440171005e-04, -1.5918565112213727e-03,  2.8615038576424393e-03, -9.1905878526491835e-04,
         1.3275767640848053e-02, -4.2639155802185694e-03,  7.6647680211425293e-03,  4.9562700674909976e-03,
        -2.4617728080410309e-03, -1.5918565112213643e-03,  2.8615038576424133e-03, -9.1905878526491358e-04,
         3.0444790786702609e-02, -9.7782682918080652e-03,  1.7577308156124568e-02,  4.9562700674909361e-03,
        -5.6454858304837699e-03, -1.5918565112213565e-03,  2.8615038576423899e-03, -9.1905878526490664e-04,
         4.7613813932557182e-02, -1.5292621003397559e-02,  2.7489848291106629e-02,  4.9562700674908727e-03,
        -8.8291988529265123e-03, -1.5918565112213502e-03,  2.8615038576423760e-03, -9.1905878526489938e-04,
         6.4782837078411787e-02, -2.0806973714987066e-02,  3.7402388426088666e-02,  4.9562700674908502e-03,
        -1.2012911875369252e-02, -1.5918565112213396e-03,  2.8615038576423270e-03, -9.1905878526489797e-04,
         8.1951860224266287e-02, -2.6321326426576551e-02,  4.7314928561070724e-02,  4.9562700674907565e-03,
        -1.5196624897811980e-02, -1.5918565112213333e-03,  2.8615038576423326e-03, -9.1905878526488875e-04,
                                                          
        -3.7435384457706450e-03,  2.1613329293873962e-03, -2.1613329293874040e-03, -7.1070190998701816e-03,
         1.2478461485902145e-03,  4.1032393904458621e-03, -4.1032393904458613e-03,  2.3690063666233913e-03,
        -4.6726594701638740e-03,  2.6977612029305646e-03, -2.6977612029305677e-03,  1.0995468558481671e-03,
         1.5575531567212943e-03, -6.3482367321053882e-04,  6.3482367321053459e-04, -3.6651561861605436e-04,
        -8.6371743090069562e-04,  4.9866749123429245e-04, -4.9866749123429516e-04,  1.0995468558481431e-03,
         2.8790581030023514e-04, -6.3482367321053210e-04,  6.3482367321053307e-04, -3.6651561861604824e-04,
         2.9452246083624873e-03, -1.7004262204619832e-03,  1.7004262204619765e-03,  1.0995468558481299e-03,
        -9.8174153612082390e-04, -6.3482367321052331e-04,  6.3482367321052440e-04, -3.6651561861604314e-04,
         6.7541666476256586e-03, -3.8995199321582566e-03,  3.8995199321582427e-03,  1.0995468558480971e-03,
        -2.2513888825418833e-03, -6.3482367321051475e-04,  6.3482367321052418e-04, -3.6651561861603663e-04,
         1.0563108686888851e-02, -6.0986136438545271e-03,  6.0986136438545184e-03,  1.0995468558480518e-03,
        -3.5210362289629427e-03, -6.3482367321051041e-04,  6.3482367321053936e-04, -3.6651561861602774e-04,
         1.4372050726152056e-02, -8.2977073555508037e-03,  8.2977073555507794e-03,  1.0995468558480691e-03,
        -4.7906835753840022e-03, -6.3482367321049816e-04,  6.3482367321049653e-04, -3.6651561861602606e-04,
         1.8180992765415183e-02, -1.0496801067247074e-02,  1.0496801067247093e-02,  1.0995468558479870e-03,
        -6.0603309218050499e-03, -6.3482367321049653e-04,  6.3482367321053340e-04, -3.6651561861601522e-04,
  };
  return sol;
}

