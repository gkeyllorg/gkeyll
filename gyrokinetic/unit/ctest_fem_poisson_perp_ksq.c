/* Test the perpendicular FEM Helmholtz solver with a nonzero kSq,
 * essentially solving the Helmholtz equation:
 *   - nabla . (epsilon * nabla phi) + kSq * phi = rho
 * using the gkyl_fem_poisson_perp struct. Also test the LHS update routine,
 * 
 * This test contains additional optional environment variables:
 *   TEST_NX, TEST_NY, TEST_NZ: to set the grid resolution (default: 32x16x2 for 3x tests, 32x48 for 2x tests)
 *   TEST_OUTPUT: when set, enables writing the solution in the ctest_fem_poisson_perp_ksq.txt file (default: no output file written)
 *   TEST_VERBOSE: when set, enables writing the DG L2 error in standard output (default: no extra output)
 * 
 *  Example:
 *   TEST_NX=32 TEST_NY=32 TEST_NZ=16 TEST_OUTPUT=1 TEST_VERBOSE=1 ./ctest_fem_poisson_perp_kSq test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_ho
*/

#include <acutest.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>

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

static double field_L2norm(struct gkyl_rect_grid grid, struct gkyl_range range,
  struct gkyl_basis basis, struct gkyl_array *field)
{
  // Compute the L2 norm of a single field.
  struct gkyl_array *l2_cell = gkyl_array_new(GKYL_DOUBLE, 1, field->size);
  gkyl_dg_calc_l2_range(&basis, 0, l2_cell, 0, field, range);
  gkyl_array_scale_range(l2_cell, grid.cellVolume, &range);
  double l2[1];
  gkyl_array_reduce_range(l2, l2_cell, GKYL_SUM, &range);
  gkyl_array_release(l2_cell);
  return sqrt(l2[0]);
}

// Return true when the TEST_OUTPUT env var is set (enables results file writing).
static bool helmholtz_write_output(void) { return getenv("TEST_OUTPUT") != NULL; }
// Return true when the HELMHOLTZ_VERBOSE env var is set (enables extra printf output).
static bool helmholtz_verbose(void) { return getenv("HELMHOLTZ_VERBOSE") != NULL; }

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // allocate array (filled with zeros)
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// We want to test the solution to the Helmholtz equation: -d^2(phi)/dx^2 - kSq*phi = rho

static double ksquare() { return 20.0; }  // To get kSq everywhere in one place for easy editing.
static double x_dirichletbc_lo() { return 0.0; } // to set BC values in one place for easy editing.
static double x_dirichletbc_up() { return 0.5; } // Only effective for 2x problem.
static double ksq_factor() { return 3.0; } // To set the factor we introduce on kSq for testing, in one place for easy editing.
static double eps_factor() { return 4.0; } // To set the factor we introduce on epsilon for testing, in one place for easy editing.

static double p1_func(double x)
{
  // This is a linear polynomial that is compatible with the BCs above. (assumes domain is [0,1])
  double a = -ksquare()*x_dirichletbc_lo();
  double b = -ksquare()*x_dirichletbc_up();
  return (b - a)*x + a;
}

static void evalFunc_ksquare(double t, const double *xn, double* restrict fout, void *ctx) 
{ 
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = ksquare();
  fout[0] *= cos(2*M_PI*x); 
  fout[0] *= exp(4*M_PI*y); 
  fout[0] *= sin(z)+2;
}

// RHS: rho(x,z) = (1.+kz*z+0.5*pow(z,2)) * sum_{m=1}^2 b_m sin(2*pi*m*x)
void evalFunc_rhs_dirichletx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  // These values have to match those in the test below.
  double bmn[] = {1.5, 1.5, -0.5};
  fout[0] = 0.;
  for (int m=1; m<3; m++) {
    double b = bmn[m-1];
    double t2 = b*sin(2*M_PI*m*x);
    fout[0] += t2;
  }
  double kz = 1.;
  fout[0] *= (1.+kz*z+0.5*pow(z,2));
  
  // fout[0] *= 0.0; // To test exact integration (linear polynomial)
  fout[0] += p1_func(x);
}
// Solution: phi(x,z) = (1.+kz*z+0.5*pow(z,2)) * sum_{m=1}^2 b_m*sin(2*pi*m*x)/( (2*pi*m)^2 - kSq )
void evalFunc_sol_dirichletx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  // These values have to match those in the test below.
  double bmn[] = {1.5, 1.5, -0.5};
  fout[0] = 0.;
  for (int m=1; m<3; m++) {
    double b = bmn[m-1];
    double kx = 2*M_PI*m;
    double t2 = b*sin(kx*x)/(kx*kx - ksquare());
    fout[0] += t2;
  }
  double kz = 1.;
  fout[0] *= (1.+kz*z+0.5*pow(z,2));
  
  // fout[0] *= 0.0; // To test exact integration (linear polynomial)
  fout[0] += -p1_func(x)/ksquare();
}
// Periodic case (This is not suitable for convergence test as it leaves a constant offset in the solution)
void evalFunc_rhs_periodicx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  // These values have to match those in the test below.
  double gxx = 1.0;
  double amn[] = {-0.6, -2., 1.5};
  double bmn[] = { 1.,  0.7, -0.3};
  fout[0] = 0.;
  for (int m=1; m<3; m++) {
    double a = amn[m-1];
    double b = bmn[m-1];
    double kx = 2*M_PI*m;
    double t1 = a*cos(kx*x);
    double t2 = b*sin(kx*x);
    fout[0] += t1+t2;
  }
  double kz = 0.;
  fout[0] *= (1.+kz*z);

}
void evalFunc_sol_periodicx_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  // These values have to match those in the test below.
  double amn[] = {-0.6, -2., 1.5};
  double bmn[] = { 1.,  0.7, -0.3};
  fout[0] = 0.;
  for (int m=1; m<3; m++) {
    double a = amn[m-1];
    double b = bmn[m-1];
    double kx = 2*M_PI*m;
    double t1 = a*cos(kx*x)/(kx*kx - ksquare());
    double t2 = b*sin(kx*x)/(kx*kx - ksquare());
    fout[0] += t1+t2;
  }
  double kz = 0.;
  fout[0] *= (1.+kz*z);
}

double trig_func(double x, double y, bool laplacian)
{
  double amn[] = {0.0, -1., 0., -0.0};
  double out = 0.;
  for (int m=0; m<4; m++) {
    double a = amn[m];
    double kx = M_PI*m;
    double ky = 2*M_PI*m*0;
    double t2 = a*cos(ky*y);
    if (m > 0)
      t2 *= sin(kx*x);

    if (laplacian) {
      t2 *= -(kx*kx + ky*ky);
    }
    out += t2;
  }
  return out;
}

double z_envelope(double z)
{
  double kz = 1.;
  return 1.0;
}

// 3x Dirichlet-Dirichlet cases — manufactured solution approach.
// Choose phi = sin(kx*x)*sin(ky*y)*z_envelope(z), which vanishes at x=0,1 and y=0,1.
// Derive rho from: -nabla_perp^2 phi + kSq(x,y)*phi = rho
//   => rho = (kx^2 + ky^2)*phi + kSq(x,y)*phi
// This is consistent with any kSq(x,y), including the spatially varying evalFunc_ksquare.
void evalFunc_sol_dirichletx_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double kx = M_PI, ky = M_PI;
  fout[0] = sin(kx*x) * sin(ky*y) * z_envelope(z);
}
void evalFunc_rhs_dirichletx_dirichlety_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double kx = M_PI, ky = M_PI;
  double kSq_val[1];
  evalFunc_ksquare(t, xn, kSq_val, ctx);
  double phi_val = sin(kx*x) * sin(ky*y) * z_envelope(z);
  fout[0] = (kx*kx + ky*ky - kSq_val[0]) * phi_val;
}

void evalFunc_sol_dirichletx_periodicy_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = trig_func(x,y,false);
  fout[0] *= z_envelope(z);
}

// Periodic case (This is not suitable for convergence test as it leaves a constant offset in the solution)
void evalFunc_rhs_dirichletx_periodicy_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double kSq[1];
  evalFunc_ksquare(t, xn, kSq, ctx);
  double fsol[1];
  evalFunc_sol_dirichletx_periodicy_3x(t, xn, fsol, ctx);
  fout[0] = -trig_func(x,y,true)*z_envelope(z) - kSq[0]*fsol[0];
}

// Read 2x grid resolution from environment variables (or defaults).
static void get_2x_cells(int *cells) {
  int nx = 32, nz = 48;
  char *env;
  if ((env = getenv("TEST_NX"))) nx = atoi(env);
  if ((env = getenv("TEST_NZ"))) nz = atoi(env);
  cells[0] = nx; cells[1] = nz;
}

// Read 3x grid resolution from environment variables (or defaults).
static void get_3x_cells(int *cells) {
  int nx = 8, ny = 8, nz = 4;
  char *env;
  if ((env = getenv("TEST_NX"))) nx = atoi(env);
  if ((env = getenv("TEST_NY"))) ny = atoi(env);
  if ((env = getenv("TEST_NZ"))) nz = atoi(env);
  cells[0] = nx; cells[1] = ny; cells[2] = nz;
}

// Introduce eps and ksq factors.
static void introduce_factors(struct gkyl_array *eps, struct gkyl_array *kSqFld) {
  gkyl_array_scale(eps, eps_factor());
  gkyl_array_scale(kSqFld, ksq_factor());
}
// Remove eps and ksq factors.
static void remove_factors(struct gkyl_array *eps, struct gkyl_array *kSqFld) {
  gkyl_array_scale(eps, 1.0/eps_factor());
  gkyl_array_scale(kSqFld, 1.0/ksq_factor());
}

static double
solve_fem_helmholtz_perp_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu,
  double *sol_L2_out)
{
  double epsilon_0 = 1.0;
  double kSq = ksquare();  // Helmholtz wave number squared.

  double lower[] = {0.0, -M_PI}, upper[] = {1.0, M_PI};
  int dim = sizeof(lower)/sizeof(lower[0]);
  int dim_perp = dim - 1;

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
  gkyl_proj_on_basis *projob = NULL, *projob_sol = NULL;
  if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order + 1, 1, evalFunc_rhs_dirichletx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order + 1, 1, evalFunc_sol_dirichletx_2x, NULL);
  } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order + 1, 1, evalFunc_rhs_periodicx_2x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order + 1, 1, evalFunc_sol_periodicx_2x, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor (1 component for 1D perp).
  int epsnum = dim_perp + (int)ceil((pow(3., dim_perp - 1) - dim_perp) / 2);
  struct gkyl_array *eps = mkarr(use_gpu, epsnum * basis.num_basis, localRange_ext.volume);
  // Create DG field for kSq.
  struct gkyl_array *kSqFld = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Analytic solution.
  struct gkyl_array *phisol_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, rho->ncomp, rho->size);
    phi_ho = mkarr(false, phi->ncomp, phi->size);
  } else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.), dim);
  gkyl_array_shiftc(eps, epsilon_0 * dg0norm, 0 * basis.num_basis);

  // Project kSq onto the basis.
  gkyl_array_shiftc(kSqFld, kSq * dg0norm, 0);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol_ho);

  // Here we test if we can change kSq and epsilon after creating the solver.
  introduce_factors(eps, kSqFld);

  // FEM Helmholtz solver.
  struct gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&localRange, &grid, basis,
    &bcs, NULL, eps, kSqFld, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_perp_set_rhs(poisson, rho);

  // Solve the problem once, to allow kSq and epsilon modification.
  gkyl_fem_poisson_perp_solve(poisson, phi);

  // Remove the factor and correct the solver's internal kSq and epsilon field accordingly.
  remove_factors(eps, kSqFld);
  gkyl_fem_poisson_perp_update_lhs(poisson, eps, kSqFld);

  // Set the RHS source.
  gkyl_fem_poisson_perp_set_rhs(poisson, rho);

  // Solve the problem.
  gkyl_fem_poisson_perp_solve(poisson, phi);
  gkyl_array_copy(phi_ho, phi);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) cudaDeviceSynchronize();
#endif

  // Write data to text file for visualization (Python notebook can read this)
  FILE *fp = helmholtz_write_output() ? fopen("ctest_fem_poisson_perp_ksq_results.txt", "w") : NULL;
  if (fp) {
    fprintf(fp, "# x z phi_numerical phi_analytical\n");
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &localRange);
    while (gkyl_range_iter_next(&iter)) {
      double xc[2];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);
      long loc = gkyl_range_idx(&localRange, iter.idx);
      const double *phi_p = gkyl_array_cfetch(phi_ho, loc);
      const double *phisol_p = gkyl_array_cfetch(phisol_ho, loc);
      // Write cell center coordinates and the 0th basis coefficient (cell average)
      fprintf(fp, "%.16e %.16e %.16e %.16e\n", xc[0], xc[1], phi_p[0]/dg0norm, phisol_p[0]/dg0norm);
    }
    fclose(fp);
  }

  double err_L2 = error_L2norm(grid, localRange, basis, phi_ho, phisol_ho);
  if (sol_L2_out)
    *sol_L2_out = field_L2norm(grid, localRange, basis, phisol_ho);

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(kSqFld);
  gkyl_array_release(phi);
  gkyl_array_release(phisol_ho);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  return err_L2;
}

void
test_fem_helmholtz_perp_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  if (helmholtz_write_output()) {
    // Single-resolution mode (called from Python): check relative L2 error < 1%.
    double sol_L2;
    double err_L2 = solve_fem_helmholtz_perp_2x(poly_order, cells, bcs, use_gpu, &sol_L2);
    double rel_err = (sol_L2 > 0.0) ? err_L2 / sol_L2 : err_L2;
    TEST_CHECK(rel_err < 0.01);
    TEST_MSG("Relative L2 error = %.6e (threshold 1%%)", rel_err);
    if (helmholtz_verbose())
      printf("L2 error = %.6e, relative L2 error = %.6e\n", err_L2, rel_err);
  } else {
    // Default convergence check: run at baseline and 2x resolution, verify error ratio >= 3.5
    // (consistent with 2nd-order convergence: factor of 4 = 2^2).
    int cells2[2] = {2*cells[0], 2*cells[1]};
    double err_coarse = solve_fem_helmholtz_perp_2x(poly_order, cells,  bcs, use_gpu, NULL);
    double err_fine   = solve_fem_helmholtz_perp_2x(poly_order, cells2, bcs, use_gpu, NULL);
    double ratio = err_coarse / err_fine;
    TEST_CHECK(ratio >= 3.5);
    TEST_MSG("L2 error ratio (2x refinement) = %.4f (expected >= 3.5 for 2nd order)", ratio);
    if (helmholtz_verbose())
      printf("L2 error: coarse = %.6e, fine = %.6e, ratio = %.4f\n", err_coarse, err_fine, ratio);
  }
}

static double
solve_fem_helmholtz_perp_3x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu,
  double *sol_L2_out)
{
  double epsilon_0 = 1.0;
  double kSq = ksquare();  // Helmholtz wave number squared.

  double lower[] = {0.0, 0.0, -M_PI}, upper[] = {1.0, 1.0, M_PI};
  int dim = sizeof(lower)/sizeof(lower[0]);
  int dim_perp = dim - 1;

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
  gkyl_proj_on_basis *projob = NULL, *projob_sol = NULL, *projob_kSq = NULL;
  if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
      (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order + 1, 1, evalFunc_rhs_dirichletx_dirichlety_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2 * (poly_order + 1), 1, evalFunc_sol_dirichletx_dirichlety_3x, NULL);
  } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order + 1, 1, evalFunc_rhs_dirichletx_periodicy_3x, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2 * (poly_order + 1), 1, evalFunc_sol_dirichletx_periodicy_3x, NULL);
  }
  projob_kSq = gkyl_proj_on_basis_new(&grid, &basis, 
    2 * (poly_order + 1), 1, evalFunc_ksquare, NULL);

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim_perp + (int)ceil((pow(3., dim_perp - 1) - dim_perp) / 2);
  struct gkyl_array *eps = mkarr(use_gpu, epsnum * basis.num_basis, localRange_ext.volume);
  // Create DG field for kSq.
  struct gkyl_array *kSqFld = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Analytic solution.
  struct gkyl_array *phisol_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho, *kSq_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, rho->ncomp, rho->size);
    phi_ho = mkarr(false, phi->ncomp, phi->size);
    kSq_ho = mkarr(false, kSqFld->ncomp, kSqFld->size);
  } else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
    kSq_ho = gkyl_array_acquire(kSqFld);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.), dim);
  gkyl_array_shiftc(eps, epsilon_0 * dg0norm, 0 * basis.num_basis);
  gkyl_array_shiftc(eps, 0.0 * dg0norm, 1 * basis.num_basis);           // gxy = 0
  gkyl_array_shiftc(eps, epsilon_0 * dg0norm, 2 * basis.num_basis);

  // Project kSq onto the basis.
  gkyl_proj_on_basis_advance(projob_kSq, 0.0, &localRange, kSq_ho);
  gkyl_array_copy(kSqFld, kSq_ho);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol_ho);

  // Here we test if we can change kSq and epsilon after creating the solver.
  introduce_factors(eps, kSqFld);

  // FEM Helmholtz solver.
  struct gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&localRange, &grid, basis,
    &bcs, NULL, eps, kSqFld, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_perp_set_rhs(poisson, rho);

  // Solve the problem once, to allow kSq and epsilon modification.
  gkyl_fem_poisson_perp_solve(poisson, phi);

  // Remove the factor and correct the solver's internal kSq and epsilon field accordingly.
  remove_factors(eps, kSqFld);
  gkyl_fem_poisson_perp_update_lhs(poisson, eps, kSqFld);

  // Set the RHS source.
  gkyl_fem_poisson_perp_set_rhs(poisson, rho);

  // Solve the problem.
  gkyl_fem_poisson_perp_solve(poisson, phi);
  gkyl_array_copy(phi_ho, phi);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) cudaDeviceSynchronize();
#endif

  // Write data to text file for visualization (Python notebook can read this)
  FILE *fp = helmholtz_write_output() ? fopen("ctest_fem_poisson_perp_ksq_results.txt", "w") : NULL;
  if (fp) {
    fprintf(fp, "# x y z phi_numerical phi_analytical\n");
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &localRange);
    while (gkyl_range_iter_next(&iter)) {
      double xc[3];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);
      long loc = gkyl_range_idx(&localRange, iter.idx);
      const double *phi_p = gkyl_array_cfetch(phi_ho, loc);
      const double *phisol_p = gkyl_array_cfetch(phisol_ho, loc);
      // Write cell center coordinates and the 0th basis coefficient (cell average)
      fprintf(fp, "%.16e %.16e %.16e %.16e %.16e\n", xc[0], xc[1], xc[2], phi_p[0]/dg0norm, phisol_p[0]/dg0norm);
    }
    fclose(fp);
  }

  double err_L2 = error_L2norm(grid, localRange, basis, phi_ho, phisol_ho);
  if (sol_L2_out)
    *sol_L2_out = field_L2norm(grid, localRange, basis, phisol_ho);

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_proj_on_basis_release(projob_kSq);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(kSqFld);
  gkyl_array_release(phi);
  gkyl_array_release(phisol_ho);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(kSq_ho);
  return err_L2;
}

void
test_fem_helmholtz_perp_3x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  // Validate poly_order and BC combination before running.
  if (poly_order != 1) {
    TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
    TEST_MSG("This poly_order is not available");
    return;
  }
  bool dd = (bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET);
  bool ydir_ok = (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET) ||
                 (bcs.lo_type[1] == GKYL_POISSON_PERIODIC  && bcs.up_type[1] == GKYL_POISSON_PERIODIC);
  if (!dd || !ydir_ok) {
    TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
    TEST_MSG("This BC combination is not available");
    return;
  }

  if (helmholtz_write_output()) {
    // Single-resolution mode (called from Python): check relative L2 error < 1%.
    double sol_L2;
    double err_L2 = solve_fem_helmholtz_perp_3x(poly_order, cells, bcs, use_gpu, &sol_L2);
    double rel_err = (sol_L2 > 0.0) ? err_L2 / sol_L2 : err_L2;
    TEST_CHECK(rel_err < 0.01);
    TEST_MSG("Relative L2 error = %.6e (threshold 1%%)", rel_err);
    if (helmholtz_verbose())
      printf("L2 error = %.6e, relative L2 error = %.6e\n", err_L2, rel_err);
  } else {
    // Default convergence check: run at baseline and 2x resolution, verify error ratio >= 3.5
    int cells2[3] = {2*cells[0], 2*cells[1], 2*cells[2]};
    double err_coarse = solve_fem_helmholtz_perp_3x(poly_order, cells,  bcs, use_gpu, NULL);
    double err_fine   = solve_fem_helmholtz_perp_3x(poly_order, cells2, bcs, use_gpu, NULL);
    double ratio = err_coarse / err_fine;
    TEST_CHECK(ratio >= 3.5);
    TEST_MSG("L2 error ratio (2x refinement) = %.4f (expected >= 3.5 for 2nd order)", ratio);
    if (helmholtz_verbose())
      printf("L2 error: coarse = %.6e, fine = %.6e, ratio = %.4f\n", err_coarse, err_fine, ratio);
  }
}

// 2x test wrappers
void test_fem_poisson_perp_ksq_2x_p1_dirichletx_ho() {
  int cells[2]; get_2x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET; bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = x_dirichletbc_lo(); bc_tv.up_value[0].v[0] = x_dirichletbc_up();
  test_fem_helmholtz_perp_2x(1, cells, bc_tv, false);
}

void test_fem_poisson_perp_ksq_2x_p1_periodicx_ho() {
  int cells[2]; get_2x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC; bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.; bc_tv.up_value[0].v[0] = 0.;
  test_fem_helmholtz_perp_2x(1, cells, bc_tv, false);
}

// 3x test wrappers
void test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_ho() {
  int cells[3]; get_3x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET; bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET; bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.; bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.; bc_tv.up_value[1].v[0] = 0.;
  test_fem_helmholtz_perp_3x(1, cells, bc_tv, false);
}

void test_fem_poisson_perp_ksq_3x_p1_dirichletx_periodicy_ho() {
  int cells[3]; get_3x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET; bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.; bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.; bc_tv.up_value[1].v[0] = 0.;
  test_fem_helmholtz_perp_3x(1, cells, bc_tv, false);
}

#ifdef GKYL_HAVE_CUDA
void test_fem_poisson_perp_ksq_2x_p1_dirichletx_dev() {
  int cells[2]; get_2x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = x_dirichletbc_lo();
  bc_tv.up_value[0].v[0] = x_dirichletbc_up();
  test_fem_helmholtz_perp_2x(1, cells, bc_tv, true);
}

void test_fem_poisson_perp_ksq_2x_p1_periodicx_dev() {
  int cells[2]; get_2x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_fem_helmholtz_perp_2x(1, cells, bc_tv, true);
}

void test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_dev() {
  int cells[3]; get_3x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_helmholtz_perp_3x(1, cells, bc_tv, true);
}

void test_fem_poisson_perp_ksq_3x_p1_dirichletx_periodicy_dev() {
  int cells[3]; get_3x_cells(cells);
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_helmholtz_perp_3x(1, cells, bc_tv, true);
}
#endif

TEST_LIST = {
  // 2x tests
  { "test_fem_poisson_perp_ksq_2x_p1_dirichletx_ho", test_fem_poisson_perp_ksq_2x_p1_dirichletx_ho },
  { "test_fem_poisson_perp_ksq_2x_p1_periodicx_ho", test_fem_poisson_perp_ksq_2x_p1_periodicx_ho },

  // 3x tests
  { "test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_ho", test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_ho },
  { "test_fem_poisson_perp_ksq_3x_p1_dirichletx_periodicy_ho", test_fem_poisson_perp_ksq_3x_p1_dirichletx_periodicy_ho },

#ifdef GKYL_HAVE_CUDA
  { "test_fem_poisson_perp_ksq_2x_p1_dirichletx_dev", test_fem_poisson_perp_ksq_2x_p1_dirichletx_dev },
  { "test_fem_poisson_perp_ksq_2x_p1_periodicx_dev", test_fem_poisson_perp_ksq_2x_p1_periodicx_dev },
  { "test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_dev", test_fem_poisson_perp_ksq_3x_p1_dirichletx_dirichlety_dev },
  { "test_fem_poisson_perp_ksq_3x_p1_dirichletx_periodicy_dev", test_fem_poisson_perp_ksq_3x_p1_dirichletx_periodicy_dev },
#endif
  { NULL, NULL },
};
