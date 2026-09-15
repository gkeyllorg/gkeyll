// Test the FEM Poisson solver with spatially varying permittivity.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_reduce.h>
#include <gkyl_fem_poisson.h>

void evalFunc2x_periodicx_periodicy_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
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
}
void evalFunc2x_periodicx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
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
}

double poly_test_func_1x(double x, double a, double *c)
{
  // Function that can be used to produce homogeneous Dirichlet or Neumann
  // boundary values depending on the choice of a and c. It assumes x \in [0,1].
  return pow(x,2)/2.-a*pow(x,4)/12.+c[0]*x+c[1];
}

void evalFunc2x_dirichletx_dirichlety_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double a = M_PI/3.;
  double b = 4.*M_PI/5.;
  // RHSsource:
  //   gxx*(2*(x-a)*(y-M_PI)*(y+M_PI)*(y-b)+2*(x+M_PI)*(y-M_PI)*(y+M_PI)*(y-b)+2*(x-M_PI)*(y-M_PI)*(y+M_PI)*(y-b))+2*gxy*((x+M_PI)*(x-a)*(y+M_PI)*(y-b)+(x-M_PI)*(x-a)*(y+M_PI)*(y-b)+(x-M_PI)*(x+M_PI)*(y+M_PI)*(y-b)+(x+M_PI)*(x-a)*(y-M_PI)*(y-b)+(x-M_PI)*(x-a)*(y-M_PI)*(y-b)+(x-M_PI)*(x+M_PI)*(y-M_PI)*(y-b)+(x+M_PI)*(x-a)*(y-M_PI)*(y+M_PI)+(x-M_PI)*(x-a)*(y-M_PI)*(y+M_PI)+(x-M_PI)*(x+M_PI)*(y-M_PI)*(y+M_PI))+gyy*(2*(x-M_PI)*(x+M_PI)*(x-a)*(y-b)+2*(x-M_PI)*(x+M_PI)*(x-a)*(y+M_PI)+2*(x-M_PI)*(x+M_PI)*(x-a)*(y-M_PI))
  fout[0] = ((x+M_PI)*(x-a)*(x-M_PI))*((y+M_PI)*(y-b)*(y-M_PI));
}
void evalFunc2x_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double a = M_PI/3.;
  double b = 4.*M_PI/5.;
  // Expected solution: ((x+M_PI)*(x-a)*(x-M_PI))*((x+M_PI)*(x-b.)*(x-M_PI));
  fout[0] = -( gxx*(2*(x-a)*(y-M_PI)*(y+M_PI)*(y-b)+2*(x+M_PI)*(y-M_PI)*(y+M_PI)*(y-b)+2*(x-M_PI)*(y-M_PI)*(y+M_PI)*(y-b))+2*gxy*((x+M_PI)*(x-a)*(y+M_PI)*(y-b)+(x-M_PI)*(x-a)*(y+M_PI)*(y-b)+(x-M_PI)*(x+M_PI)*(y+M_PI)*(y-b)+(x+M_PI)*(x-a)*(y-M_PI)*(y-b)+(x-M_PI)*(x-a)*(y-M_PI)*(y-b)+(x-M_PI)*(x+M_PI)*(y-M_PI)*(y-b)+(x+M_PI)*(x-a)*(y-M_PI)*(y+M_PI)+(x-M_PI)*(x-a)*(y-M_PI)*(y+M_PI)+(x-M_PI)*(x+M_PI)*(y-M_PI)*(y+M_PI))+gyy*(2*(x-M_PI)*(x+M_PI)*(x-a)*(y-b)+2*(x-M_PI)*(x+M_PI)*(x-a)*(y+M_PI)+2*(x-M_PI)*(x+M_PI)*(x-a)*(y-M_PI)) );
}

void evalFunc2x_dirichletx_periodicy_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double m = 5.;
  // RHSsource: (2*(2*M_PI - 45*x + 5*(M_PI - 9*x)*sin(m*y))*gxx)/15.
  //           +(m*((3*pow(M_PI,2) + 2*M_PI*x - 9*pow(x,2))*cos(m*y)*2.*gxy
  //                +m*(M_PI - 3*x)*(M_PI - x)*(M_PI + x)*sin(m*y)*gyy))/3.;
  fout[0] = ((x-4.*M_PI/5.) + (x-M_PI)*sin(m*y))*(x+M_PI)*(x-M_PI/3.);
}
void evalFunc2x_dirichletx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double m = 5.;
  // Expected solution: ((x-4.*M_PI/5.) + (x-M_PI)*sin(m*y))*(x+M_PI)*(x-M_PI/3.);
  fout[0] = (2*(2*M_PI - 45*x + 5*(M_PI - 9*x)*sin(m*y))*gxx)/15.
           +(m*((3*pow(M_PI,2) + 2*M_PI*x - 9*pow(x,2))*cos(m*y)*2.*gxy
                +m*(M_PI - 3*x)*(M_PI - x)*(M_PI + x)*sin(m*y)*gyy))/3.;
}

void evalFunc2x_periodicx_dirichlety_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double m = 5.;
  // RHSsource: (2*(2*M_PI - 45*y + 5*(M_PI - 9*y)*sin(m*x))*gxx)/15.
  //           +(m*((3*pow(M_PI,2) + 2*M_PI*y - 9*pow(y,2))*cos(m*x)*2.*gxy
  //                +m*(M_PI - 3*y)*(M_PI - y)*(M_PI + y)*sin(m*x)*gyy))/3.;
  fout[0] = ((y-4.*M_PI/5.) + (y-M_PI)*sin(m*x))*(y+M_PI)*(y-M_PI/3.);
}
void evalFunc2x_periodicx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double m = 5.;
  // Expected solution: ((y-4.*M_PI/5.) + (y-M_PI)*sin(m*x))*(y+M_PI)*(y-M_PI/3.);
  fout[0] = (5*pow(m,2)*(M_PI - 3*y)*(M_PI - y)*(M_PI + y)*sin(m*x)*gxx
    +5*m*(3*pow(M_PI,2) + 2*M_PI*y - 9*pow(y,2))*cos(m*x)*(gxy + gxy)
    +2*(2*M_PI - 45*y + 5*(M_PI - 9*y)*sin(m*x))*gyy)/15.;
}

void evalFunc2x_dirichletx_neumanny_dirichlety_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 0.;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double a = 2.;
  double c[] = {a/12.-0.5, 0.};
  double b = 5.;
  double d[] = {0., b/12.-0.5};
  // RHSsource:  (-3*((6 - 6*a*pow(x,2))*pow(y,2) + b*(-1 + a*pow(x,2))*pow(y,4)
  //     +4*(9*x - 5*a*pow(x,3) + 6*c[0])*d[0] - 12*(-1 + a*pow(x,2))*d[1])*gxx
  //     -4*y*(-3 + b*pow(y,2))*(-3*x + a*pow(x,3) - 3*c[0])*(gxy + gxy)
  //     -3*(-1 + b*pow(y,2))*(-6*pow(x,2) + a*pow(x,4) - 12*x*c[0] - 12*c[1])*gyy)/36.
  fout[0] = poly_test_func_1x(x, a, c)*poly_test_func_1x(y, b, d);
}
void evalFunc2x_dirichletx_neumanny_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 0.;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double a = 2.;
  double c[] = {a/12.-0.5, 0.};
  double b = 5.;
  double d[] = {0., b/12.-0.5};
  // Expected solution: (x^2/2 - (a*x^4)/12 + c*x + q)*(y^2/2 - (b*y^4)/12 + f*x + h)
  fout[0] = (-3*((6 - 6*a*pow(x,2))*pow(y,2) + b*(-1 + a*pow(x,2))*pow(y,4)
    +4*(9*x - 5*a*pow(x,3) + 6*c[0])*d[0] - 12*(-1 + a*pow(x,2))*d[1])*gxx
    -4*y*(-3 + b*pow(y,2))*(-3*x + a*pow(x,3) - 3*c[0])*(gxy + gxy)
    -3*(-1 + b*pow(y,2))*(-6*pow(x,2) + a*pow(x,4) - 12*x*c[0] - 12*c[1])*gyy)/36.;
}

void evalFunc2x_neumannx_dirichletx_dirichlety_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 0.;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double a = 5.;
  double c[] = {0., a/12.-0.5};
  double b = 2.;
  double d[] = {b/12.-0.5, 0.};
  // RHSsource:  
  //   gxx*(1.-a*pow(x,2))*(-((b*pow(y,4))/12)+pow(y,2)/2+d[0]*y+d[1])+2*gxy*(-((a*pow(x,3))/3)+x+c[0])*(-((b*pow(y,3))/3)+y+d[0])+gyy*(-((a*pow(x,4))/12)+pow(x,2)/2+c[0]*x+c[1])*(1.-b*pow(y,2));
  fout[0] = poly_test_func_1x(x, a, c)*poly_test_func_1x(y, b, d);
}
void evalFunc2x_neumannx_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 0.;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double a = 5.;
  double c[] = {0., a/12.-0.5};
  double b = 2.;
  double d[] = {b/12.-0.5, 0.};
  // Expected solution: (x^2/2 - (a*x^4)/12 + c*x + q)*(y^2/2 - (b*y^4)/12 + f*x + h)
  fout[0] = -( gxx*(1.-a*pow(x,2))*(-((b*pow(y,4))/12)+pow(y,2)/2+d[0]*y+d[1])+2*gxy*(-((a*pow(x,3))/3)+x+c[0])*(-((b*pow(y,3))/3)+y+d[0])+gyy*(-((a*pow(x,4))/12)+pow(x,2)/2+c[0]*x+c[1])*(1.-b*pow(y,2)) );
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

// Apply periodic BCs in one direction.
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
}

void
test_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  // Determinant of g tensor has to be >0. Diagonal entries have to be  >0.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  if (
      ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
       (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) ||
      ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
       (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
     ){
    lower[0] = 0.;  lower[1] = 0.;
    upper[0] = 1.;  upper[1] = 1.;
    gxy = 0.;
  } else if (
      ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
       (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
    )
    gxy = 2.;
  {
  }
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // Local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // Skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_periodicy, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc2x_periodicx_periodicy_sol, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_dirichlety, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc2x_dirichletx_dirichlety_sol, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc2x_dirichletx_periodicy_sol, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_dirichlety, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc2x_periodicx_dirichlety_sol, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_neumanny_dirichlety, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc2x_dirichletx_neumanny_dirichlety_sol, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_neumannx_dirichletx_dirichlety, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc2x_neumannx_dirichletx_dirichlety_sol, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim+ceil((pow(3.,dim-1)-dim)/2);
  struct gkyl_array *eps = use_gpu? mkarr_cu(epsnum*basis.num_basis, localRange_ext.volume)
                                  : mkarr(epsnum*basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_cu, *phi_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  }

  // Project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  if (use_gpu) gkyl_array_copy(rho_cu, rho);

  // Project the expected solution.
  struct gkyl_array *phisol = mkarr(basis.num_basis, localRange_ext.volume);
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol);

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, gxx*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps, gxy*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, gyy*dg0norm, 2*basis.num_basis);

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bcs, NULL, eps, NULL, false, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_poisson_set_rhs(poisson, rho_cu, NULL);
  else
    gkyl_fem_poisson_set_rhs(poisson, rho, NULL);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_poisson_solve(poisson, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_poisson_solve(poisson, phi);
  }

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
    struct gkyl_array *sol_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, localRange_ext.volume);
    double* sol_avg = (double*) gkyl_malloc(sizeof(double));
    gkyl_array_clear(sol_cellavg, 0.0);
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim)/localRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_dg_calc_average_range(&basis, 0, sol_cellavg, 0, phi, localRange);
    gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &localRange);
    gkyl_array_shiftc(phi, mavgfac*sol_avg[0], 0);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8:
      const double sol[256] = {
        -1.0487258355542172e+01, -3.4672194615719323e+00, 1.2487089948838117e+00, -7.3829839001277786e-01, -1.5453564258762444e+00, -3.3771781285781359e+00,
        3.9139004912760469e+00, 7.9028377785493853e-01, 5.2215272976960074e+00, 4.6956135527195275e-01, -7.0383518967097811e-03, 1.4306322986482731e+00,
        1.4932886788449311e+00, 3.3046036633858393e+00, -2.1454612183001176e+00, 2.0618014110525743e-01, -2.9431153863741066e-02, 2.4254609580711044e+00,
        1.2663185129853838e+00, -7.1375341867481512e-01, 5.9015118519390191e+00, 8.6215062070894122e-01, 2.1579130279631680e+00, -1.8882422542815416e-01,
        5.2951622117098962e+00, 5.7219714822887291e-01, -2.5079891559724854e+00, 2.1419510039319974e-02, -5.8494441049077164e+00, -7.8957615551664595e-01,
        -3.9263523009390977e+00, -8.0763969353204146e-01, -6.7270067596276100e+00, 5.6382017326939247e+00, -3.3771781285781368e+00, -1.9324587859588025e+00,
        -9.4839930397771006e+00, -1.2061958574561431e+00, 1.7854313575817238e+00, -2.0191560051551658e+00, -3.6051177810056378e+00, -5.5656272675016334e+00,
        1.6087388556022115e+00, -4.9776289294830867e-01, -3.0508358878078150e+00, -5.9281552051739945e+00, -1.2887240553575130e+00, 2.8845729054441199e-01,
        -3.7896827497783017e+00, -4.5964432291930963e+00, 8.6215062070894188e-01, 4.8040697727120762e-01, 5.1278376137143153e+00, -1.3088316504130526e+00,
        4.2863821616574915e+00, 1.4176964527283813e+00, 1.4121807290411542e+01, 4.5238687640008095e+00, 9.0628865226698430e-01, 1.9498147016359040e+00,
        7.4069913138705905e+00, 8.4431827130431927e+00, -4.7830894638817023e+00, 3.1300226188237157e-01, 1.0319933757546561e+01, 4.2038539637561216e+00,
        -3.0271020005688309e+00, 2.1345753327018371e+00, 3.6051177810055952e+00, 8.7631975132878779e+00, -8.4969881104589473e-01, 4.9776289294830606e-01,
        2.2148951700382993e+00, 8.9258133114739593e+00, 4.7053412370420863e-02, -4.0387661809108605e-01, -3.5530500339743374e+00, 5.6382017326939149e+00,
        -3.3771781285781279e+00, -1.4942268119085045e+00, -1.4676089183609335e+01, -1.6888264558869253e+00, -3.0447115186703853e+00, -2.7360348851286194e+00,
        -1.4121807290411509e+01, -9.8049560167887009e+00, 3.3647263189150864e+00, -1.9498147016359004e+00, 2.1412602560244745e+00, -1.1440840819343155e+01,
        6.0247601068687962e+00, 1.0053361705178678e+00, 1.4069739543380248e+01, -4.5964432291930901e+00, 8.6215062070893500e-01, 2.9462786205960989e+00,
        -1.4932886788448974e+00, -1.1024221117403657e+01, 3.1297644834797150e-01, -2.0618014110525401e-01, 2.0475705720427344e+00, -9.6624478136581367e+00,
        1.7313395933145797e+00, 9.9240032459797578e-01, 1.1825413219478111e+01, -3.3771781285781262e+00, 3.9139004912760385e+00, 2.6364018200123627e+00,
        1.7756356225280861e+01, 6.6647897073581737e+00, -4.8966895032749230e-01, 3.1613310132590238e+00, 5.8494441049076764e+00, 1.3539248625272847e+01,
        -6.3847899675871842e+00, 8.0763969353203824e-01, -1.2564260081448641e+01, 1.0704206317158961e+01, -4.2463671011837762e+00, -2.4444521332855680e+00,
        -1.6181568645540896e+01, 8.6215062070893345e-01, 2.1579130279631737e+00, -3.2378613724391472e+00, -7.2396667158749599e+00, -7.7065482108590002e+00,
        3.0046964581966891e+00, -1.7092792045714316e+00, -7.2396667158749564e+00, 7.7065482108589975e+00, -3.0046964581966829e+00, -1.7092792045714302e+00,
        -1.6181568645540882e+01, -8.6215062070893089e-01, -2.1579130279631746e+00, -3.2378613724391458e+00, -1.2564260081448637e+01, -1.0704206317158961e+01,
        4.2463671011837727e+00, -2.4444521332855706e+00, 5.8494441049076782e+00, -1.3539248625272851e+01, 6.3847899675871851e+00, 8.0763969353203946e-01,
        1.7756356225280864e+01, -6.6647897073581719e+00, 4.8966895032749130e-01, 3.1613310132590242e+00, 1.1825413219478115e+01, 3.3771781285781279e+00,
        -3.9139004912760367e+00, 2.6364018200123640e+00, 2.0475705720427291e+00, 9.6624478136581367e+00, -1.7313395933145861e+00, 9.9240032459797478e-01,
        -1.4932886788449089e+00, 1.1024221117403654e+01, -3.1297644834796867e-01, -2.0618014110525504e-01, 1.4069739543380248e+01, 4.5964432291930901e+00,
        -8.6215062070892989e-01, 2.9462786205960971e+00, 2.1412602560244904e+00, 1.1440840819343155e+01, -6.0247601068687908e+00, 1.0053361705178696e+00,
        -1.4121807290411496e+01, 9.8049560167887044e+00, -3.3647263189150931e+00, -1.9498147016358995e+00, -1.4676089183609333e+01, 1.6888264558869279e+00,
        3.0447115186703857e+00, -2.7360348851286211e+00, -3.5530500339743369e+00, -5.6382017326939131e+00, 3.3771781285781262e+00, -1.4942268119085040e+00,
        2.2148951700383011e+00, -8.9258133114739593e+00, -4.7053412370418296e-02, -4.0387661809108694e-01, 3.6051177810055952e+00, -8.7631975132878761e+00,
        8.4969881104589007e-01, 4.9776289294830811e-01, 1.0319933757546551e+01, -4.2038539637561172e+00, 3.0271020005688305e+00, 2.1345753327018371e+00,
        7.4069913138705923e+00, -8.4431827130431909e+00, 4.7830894638817041e+00, 3.1300226188237157e-01, 1.4121807290411557e+01, -4.5238687640008104e+00,
        -9.0628865226697919e-01, 1.9498147016359024e+00, 5.1278376137143313e+00, 1.3088316504130499e+00, -4.2863821616574942e+00, 1.4176964527283822e+00,
        -3.7896827497782946e+00, 4.5964432291930963e+00, -8.6215062070894288e-01, 4.8040697727120818e-01, -3.0508358878078097e+00, 5.9281552051739954e+00,
        1.2887240553575130e+00, 2.8845729054441255e-01, -3.6051177810056338e+00, 5.5656272675016352e+00, -1.6087388556022122e+00, -4.9776289294830967e-01,
        -9.4839930397771006e+00, 1.2061958574561427e+00, -1.7854313575817249e+00, -2.0191560051551658e+00, -6.7270067596276117e+00, -5.6382017326939247e+00,
        3.3771781285781368e+00, -1.9324587859588018e+00, -5.8494441049077155e+00, 7.8957615551664562e-01, 3.9263523009390977e+00, -8.0763969353204135e-01,
        5.2951622117098980e+00, -5.7219714822887602e-01, 2.5079891559724858e+00, 2.1419510039318492e-02, 5.9015118519390226e+00, -8.6215062070894533e-01,
        -2.1579130279631680e+00, -1.8882422542815386e-01, -2.9431153863737514e-02, -2.4254609580711071e+00, -1.2663185129853838e+00, -7.1375341867481512e-01,
        1.4932886788449364e+00, -3.3046036633858407e+00, 2.1454612183001189e+00, 2.0618014110525712e-01, 5.2215272976960128e+00, -4.6956135527195375e-01,
        7.0383518967087558e-03, 1.4306322986482738e+00, -1.5453564258762436e+00, 3.3771781285781368e+00, -3.9139004912760487e+00, 7.9028377785493842e-01,
        -1.0487258355542172e+01, 3.4672194615719323e+00, -1.2487089948838108e+00, -7.3829839001277808e-01,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
              (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[256] = {
        1.6461612387768028e+02, 9.5041163433732649e+01, 9.5041163433727860e+01, 5.4872041292561697e+01, 3.9412893541820802e+02, 2.2755044695912608e+02,
        3.7468120091663707e+01, 2.1632229220948631e+01, 4.5713749046729765e+02, 2.6392845317796167e+02, -1.0901138728224475e+00, -6.2937753792202944e-01,
        4.0993267971799736e+02, 2.3667474298481230e+02, -2.6163596320331422e+01, -1.5105559378508811e+01, 2.9772447182984757e+02, 1.7189130395530435e+02,
        -3.8619842709177362e+01, -2.2297176584207307e+01, 1.6390341560120643e+02, 9.4629681118453576e+01, -3.8641780127668660e+01, -2.2309842158675828e+01,
        5.1409336898312709e+01, 2.9681194497099966e+01, -2.6306706493685997e+01, -1.5188184075621486e+01, 2.9223923357232273e+00, 1.6872440017075048e+00,
        -1.6872440017075048e+00, -9.7413077857440933e-01, 3.7701335578602766e+02, 2.7586435583683095e+01, 2.1766876245114722e+02, 1.5927036010222864e+01,
        9.0866289441247636e+02, 6.9515872773418664e+01, 8.9279241789373359e+01, 8.2809358382405023e+00, 1.0595914450098408e+03, 8.3898499651530244e+01,
        -2.1406024735849591e+00, 2.2877661491188903e-02, 9.5034234310563488e+02, 7.5330921644715289e+01, -6.0934395919533266e+01, -4.9693711300289287e+00,
        6.8915221373938357e+02, 5.4099608204453631e+01, -8.9863795579745130e+01, -7.2885333999554902e+00, 3.7806552018444836e+02, 2.9016867612897943e+01,
        -8.9742190685507339e+01, -7.1929936325927386e+00, 1.1722533389750654e+02, 8.3176890873682030e+00, -6.0853961082730649e+01, -4.7576826611259841e+00,
        5.9115907301989470e+00, 3.8570496004249159e-02, -3.4130584994192588e+00, -2.2268686250830640e-02, 3.9532494908836981e+02, -1.7014232261284924e+01,
        2.2824096577354720e+02, -9.8231715761077769e+00, 9.5975101530934739e+02, -4.0020132421218086e+01, 9.7630575496775776e+01, -3.4592910741793057e+00,
        1.1268165767066630e+03, -4.5086051770062070e+01, -1.1752286517197190e+00, 5.3448117443131593e-01, 1.0117621212132726e+03, -3.9870196220703590e+01,
        -6.5251492198855004e+01, 2.4768944310450070e+00, 7.3179804407702431e+02, -2.9478026579201668e+01, -9.6385843099183120e+01, 3.5230275089404217e+00,
        3.9862330008493541e+02, -1.7147827853410664e+01, -9.5972685031834018e+01, 3.5958160445566052e+00, 1.2055374795950144e+02, -6.3960283327009160e+00,
        -6.4570845741224005e+01, 2.6117383029980612e+00, 4.3568812240045069e+00, -9.3618244791727923e-01, -2.5154465475062286e+00, 5.4050518831564409e-01,
        2.9248441181369833e+02, -4.2360779617852387e+01, 1.6886595389440930e+02, -2.4457007515449735e+01, 7.1237128462567580e+02, -1.0280462168105679e+02,
        7.3555845152776087e+01, -1.0440260970596553e+01, 8.3987621279525979e+02, -1.2057904457555561e+02, 5.9159448937202881e-02, 1.7819312776739105e-01,
        7.5508761769510738e+02, -1.0832089747962264e+02, -4.9011877654220484e+01, 6.8990513978354988e+00, 5.4401332684432566e+02, -7.8939530464794643e+01,
        -7.2851920987487702e+01, 1.0064288757334635e+01, 2.9258586917084318e+02, -4.4072911429016834e+01, -7.2309789381960641e+01, 1.0065963128703812e+01,
        8.4183252840990875e+01, -1.4602486814531128e+01, -4.8011517255901900e+01, 6.9487944556020906e+00, 5.1243280264945124e-01, -1.2834108830376940e+00,
        -2.9585321655125507e-01, 7.4097761880270807e-01, 1.3335338589153289e+02, -4.9513561034731154e+01, 7.6991613241823984e+01, -2.8586667791939774e+01,
        3.2521042488389821e+02, -1.2072280491054548e+02, 3.3777099866342660e+01, -1.2526008328551340e+01, 3.8380134263736056e+02, -1.4273590448073116e+02,
        5.0382270686085653e-02, -1.8326063399341785e-01, 3.4393929741895795e+02, -1.2905569590867466e+02, -2.3064744807979871e+01, 8.0815327356404634e+00,
        2.4445077485208702e+02, -9.4012989567055612e+01, -3.4374980477281234e+01, 1.2150383203826303e+01, 1.2638602665250625e+02, -5.1882612388218419e+01,
        -3.3789733677552128e+01, 1.2173601401435823e+01, 3.0388116477691362e+01, -1.6456149645912479e+01, -2.1634685603518761e+01, 8.2798763992707602e+00,
        -3.5421290966910353e+00, -1.0574915209925362e+00, 2.0450491874789756e+00, 6.1054301431078772e-01, -1.9656955962165007e+01, -3.8827001023298003e+01,
        -1.1348948816205576e+01, -2.2416779492626443e+01, -4.8069780771266380e+01, -9.4790622307622698e+01, -5.0552032357667791e+00, -9.8938323207046963e+00,
        -5.8092011807671106e+01, -1.1239134266121465e+02, -7.3113455098211444e-01, -2.6794831337293828e-01, -5.6147582965799181e+01, -1.0193457218070250e+02,
        1.8537510662569450e+00, 6.3051675651507413e+00, -4.7047247424851960e+01, -7.4283472062778671e+01, 3.4003301080248480e+00, 9.6592025313218883e+00,
        -3.3773606129953968e+01, -4.0585594712046607e+01, 4.2632102667110452e+00, 9.7962760282422643e+00, -1.9158284726430185e+01, -1.2149478426664299e+01,
        4.1749494799063092e+00, 6.6213233638307170e+00, -5.9635300540995999e+00, -3.4050497358331655e-01, 3.4430456820548279e+00, 1.9659063815873426e-01,
        -1.0500433735125506e+02, -1.0448332596322809e+01, -6.0624282435825428e+01, -6.0323476370696971e+00, -2.5645491158351632e+02, -2.5520589061945198e+01,
        -2.6815747366094026e+01, -2.6696236906524677e+00, -3.0514859202734101e+02, -3.0246840433680184e+01, -1.2975621459817959e+00, -5.9078811076562487e-02,
        -2.7996069295136755e+02, -2.7283987117632709e+01, 1.5839802457816219e+01, 1.7696829706659454e+00, -2.0950900058815870e+02, -1.9513864859075049e+01,
        2.4835501092947101e+01, 2.7163992069486280e+00, -1.2151783195207202e+02, -1.0073557686188563e+01, 2.5966223805407331e+01, 2.7339646805501356e+00,
        -4.3678446429614837e+01, -2.0072435332425118e+00, 1.8974366379538381e+01, 1.9231239670213580e+00, -5.4069399094138149e+00, 6.6185244344594785e-01,
        3.1216982121921970e+00, -3.8212068638732960e-01, -6.1550690131232102e+01, 3.5536307516074103e+01, -3.5536307516074103e+01, 2.0516896710410705e+01,
        -1.5032893423894598e+02, 8.6792450649845023e+01, -1.5719835617696818e+01, 9.0758513254939270e+00, -1.7876882821345188e+02, 1.0321223109175043e+02,
        -6.9994482420858450e-01, 4.0411333267471161e-01, -1.6360897243608105e+02, 9.4459684277809473e+01, 9.4524916381495334e+00, -5.4573985917983201e+00,
        -1.2165400298805480e+02, 7.0236971373148990e+01, 1.4770221266510960e+01, -8.5275912242104397e+00, -6.9482872838763285e+01, 4.0115955337528526e+01,
        1.5350794769109491e+01, -8.8627854922200644e+00, -2.3577547106177466e+01, 1.3612503168582643e+01, 1.1152657399836393e+01, -6.4389897519752148e+00,
        -2.1302889251259129e+00, 1.2299228843731245e+00, 1.2299228843731245e+00, -7.1009630837530446e-01,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
              (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8:
      const double sol[256] = {
        22.734370830726306 ,  13.125695118976289 ,  -3.945074572528958 , -2.2776898664234544,
         4.1430157490856105,   2.391971257991416 ,  -6.788649288456215 , -3.919428494123271 ,
        15.846526665407206 ,   9.148996435993684 ,  13.545674466458212 ,  7.8205987995639585,
        17.886617881804142 ,  10.32684364895111  , -12.367827253500051 , -7.1405683940987394,
         3.297982298775476 ,   1.9040909679812303,   3.945074572529099 ,  2.277689866423136 ,
        21.88933738041535  ,  12.63781482896598  ,   6.788649288455599 ,  3.9194284941235193,
        10.185826464094799 ,   5.880789650963231 , -13.545674466456994 , -7.820598799564486 ,
         8.145735247698687 ,   4.702942438005926 ,  12.367827253499309 ,  7.1405683940993345,
        52.060758559828486 ,   3.8059027307801987,   1.6206081646797366,  5.491038292974948 ,
        14.900603982252303 ,   3.8189252042585258, -23.075033413171454 , -5.48351976440672  ,
        28.648555703845823 ,  -1.7577415244777435,  31.012417040439452 ,  2.263829727390006 ,
        46.36617050189349  ,   6.115833722657592 , -20.78314736738853  ,  2.2819810610283544,
         7.561727639998614 ,   0.5575835525299617,  -1.620608164679934 , -5.4910382929748245,
        44.721882217574105 ,   0.5445610790518283,  23.075033413171248 ,  5.48351976440671  ,
        30.973930495980493 ,   6.121227807787921 , -31.0124170404393   , -2.2638297273900965,
        13.256315697933521 ,  -1.7523474393476068,  20.78314736738878  , -2.2819810610283753,
        52.42701650766966  ,  -3.594443606001256 ,  15.739268866564824 ,  2.6603742638553998,
        23.899287301995784 ,   1.3764670327475685, -32.20976100282939  ,  0.2095823315438092,
        19.74661058849544  ,  -3.381798884181076 ,  29.81221198443248  , -2.9567684395583833,
        54.14711152255226  ,  -1.6235053298316517,  -9.951073509896792 ,  3.971919696476408 ,
         9.650133262681052 ,   0.648157995902873 , -15.739268866564776 , -2.6603742638553802,
        38.177862468354775 ,  -4.322752642845835 ,  32.209761002829254 , -0.2095823315437611,
        42.3305391818548   ,   0.4355132740828546, -29.81221198443252  ,  2.956768439558362 ,
         7.930038247798141 ,  -1.322780280266687 ,   9.951073509896931 , -3.9719196964764536,
        35.07594487422056  ,  -6.423202272299062 ,  20.102030185321045 , -0.1415328420614093,
        23.825570505411957 ,  -1.4190274450990383, -26.597436855636765 ,  3.030694525617292 ,
         8.089869666290408 ,  -3.3482236251277255,  17.512425740282517 , -4.144516459476453 ,
        41.59388557523048  ,  -5.624103050052774 ,   1.8311268636775433,  2.830536860852845 ,
         9.947752721965076 ,  -0.4763273209691713, -20.102030185320988 ,  0.1415328420613958,
        21.198127090773767 ,  -5.480502148169174 ,  26.59743685563675  , -3.030694525617267 ,
        36.93382792989518  ,  -3.5513059681404244, -17.51242574028258  ,  4.144516459476464 ,
         3.429812020955032 ,  -1.2754265432153953,  -1.8311268636775244, -2.8305368608528676,
        12.487424458224453 ,  -6.618286070471654 ,  16.734971692483434 , -1.8024392851556001,
        16.8072546664208   ,  -2.6329990938003207, -14.240916558927797 ,  4.103345793934705 ,
        -1.9615785556745362,  -2.45498271156932  ,   3.404725645775862 , -4.000567987733473 ,
        20.261729728735126 ,  -6.692022870316367 ,   9.425907374512102 ,  1.5543117115136114,
         7.602058974220037 ,  -0.8779595957278498, -16.734971692483423 ,  1.8024392851555864,
         3.2822287660237572,  -4.863246572399208 ,  14.240916558927825 , -4.103345793934706 ,
        22.0510619881191   ,  -5.041262954630193 ,  -3.4047256457758874,  4.000567987733484 ,
        -0.1722462962906359,  -0.8042227958831216,  -9.425907374512114 , -1.5543117115136071,
        -6.229569715633301 ,  -4.18797555422578  ,   8.147483378324223 , -3.155549404686988 ,
         7.600917434960932 ,  -2.6822821850335097,  -0.1624478989049808,  4.02486187670718  ,
        -6.39439128250139  ,  -0.1043029095313632,  -7.917747356513808 , -2.536464848030739 ,
        -0.4325230352620455,  -5.25580953365556  ,  11.359833593930528 , -0.4377588881395075,
         5.131430749872889 ,  -0.5484582746664621,  -8.147483378324235 ,  3.1555494046869876,
        -8.699056400721336 ,  -2.0541516438587415,   0.1624478989049978, -4.024861876707185 ,
         5.296252316741011 ,  -4.632130919360892 ,   7.9177473565138055,  2.5364648480307412,
        -0.6656159304983422,   0.5193757047633152, -11.35983359393053  ,  0.4377588881395105,
        -10.810726987224193,   1.5430431702727658,  -2.4776533786098813, -2.978876162105506 ,
         0.4878062472448956,  -1.4244744739700843,   9.000864582978544 ,  1.2655790513769316,
        -1.6783429198783772,   2.8271147012030537, -10.251491387722048 ,  1.189077103392979 ,
        -9.91347862407862  ,  -0.2180227277225759,   5.496933572088938 , -2.947188017702597 ,
         3.898911143907491 ,  -0.1631372516191521,   2.477653378609876 ,  2.9788761621055095,
        -7.399622090561606 ,   2.804380392623697 ,  -9.000864582978544 , -1.2655790513769358,
        -5.233472923438328 ,  -1.4472087825494466,  10.251491387722051 , -1.1890771033929781,
         3.001662780761924 ,   1.5979286463761837,  -5.4969335720889365,  2.9471880177025964,
         4.199291539060148 ,   7.122995066418769 ,  -3.818609120416201 ,  2.2046750036022527,
         7.27861249020183  ,   5.345148286355863 ,   5.596455900479105 , -3.2311153206494816,
         9.87752213879506  ,   3.844667101141045 ,  -4.095974715264289 ,  2.364812104451738 ,
         3.122787915230512 ,   7.744514723420415 ,   0.1961270929849189, -0.1132340299302217,
        10.076515015827935 ,   3.729778510019296 ,   3.8186091204161996, -2.2046750036022535,
         6.997194064686247 ,   5.5076252900822045,  -5.59645590047911  ,  3.231115320649484 ,
         4.398284416093016 ,   7.008106475297023 ,   4.095974715264291 , -2.3648121044517394,
        11.15301863965757  ,   3.108258853017651 ,  -0.1961270929849173,  0.1132340299302219,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[256] = {
        22.546296767842215 ,  -2.171251024451982 ,  13.017110508143     ,  -1.253572363445537 ,
        50.92375002609591  ,   3.97827034412882  ,   3.3666197694257964 ,   4.80400018098302  ,
        51.06156064013447  ,  14.536294440523537 ,  -3.2870547743134426 ,   1.2916778731810277,
        35.161182830840765 ,  15.854518698111484 ,  -5.893032634099033  ,  -0.5306007432103331,
        13.982658785388029 ,  12.365469949376859 ,  -6.334393924582068  ,  -1.4838024910873295,
        -4.468909270163868 ,   5.717608387446057 ,  -4.318623859261521  ,  -2.3543421712287853,
        -9.9733872194679   ,  -2.794128262242933 ,   1.1406120334824776 ,  -2.559911274740357 ,
         4.269445835204085 ,  -3.6140123264807964,   7.082491464655293  ,   2.0865509895483414,
         5.175875409986401 ,  -7.857566422443296 ,   2.9882930612475262 ,  -4.5365680891731825,
        16.900646445855752 , -23.621518354846998 ,   3.781006652498372  ,  -4.564753802492417 ,
        24.254534337939802 , -30.01333829226902  ,   0.4647625015866842 ,   0.8744348410101994,
        21.61051037682827  , -23.678003087136048 ,  -1.9912904472782482 ,   2.783272645079833 ,
        13.53329016180549  , -12.624913045167608 ,  -2.672094818169042  ,   3.5982318660278167,
         5.1232524875614045,  -0.1795712145129375,  -2.183442696950677  ,   3.5870895900575697,
        -0.2648319678107458,   8.399365250229547 ,  -0.9273693137748725 ,   1.3659616872222011,
         7.332799079755453 ,   5.382640146845442 ,   5.3138636442906755 ,  -3.1076687377321153,
        14.573916587816605 ,  13.28352802631875  ,   8.414254665123067  ,   7.6692484817833   ,
        26.95707706022673  ,  29.427601277140873 ,  -1.2648336329587488 ,   1.6515365553958905,
        20.609671279804072 ,  27.908975624495028 ,  -2.3998428030240584 ,  -2.5283154847492426,
        11.13719978498403  ,  17.631234397628216 ,  -3.06909116440135   ,  -3.405541179243627 ,
         1.173272122316551 ,   5.488853302880039 ,  -2.6835851538256796 ,  -3.604865814412256 ,
        -4.651104356342016 ,  -5.463656340470164 ,  -0.6791201739895689 ,  -2.7185685764779564,
        -1.4512915163426436,  -9.08436798995704  ,   2.526532978519627  ,   0.6281497309886949,
         9.730736432805806 ,  -3.99819037056196  ,   3.9294138680071504 ,   2.308356286715305 ,
        18.653500652244155 , -10.9281790685398   ,  10.769603622902233  ,  -6.309387126973994 ,
        46.75824007655975  , -17.99539447939334  ,   5.456675249230666  ,   2.2291284070966295,
        52.57131235180712  ,  -9.455913547831864 ,  -2.1004964062978475 ,   2.7011432074800648,
        39.499370120928276 ,  -1.2563277193687856,  -5.446592959531293  ,   2.0328898778265354,
        19.102345888520077 ,   4.862502262358224 ,  -6.3296344717169015 ,   1.4998182592491074,
        -0.4202381019451584,   7.906348111351493 ,  -4.941734650455281  ,   0.2575469610388816,
        -9.481939583259177 ,   4.44787116675563  ,  -0.2900411390973994 ,  -2.254299555987418 ,
         3.2761876618493564,   0.2716549001526837,   7.655949338416311  ,  -0.1568400297298874,
         3.4860563616592164,   2.1712510244519323,   2.0126755788141875 ,   1.2535723634455773,
         8.698736173730346 ,  -3.9782703441288323,   0.9968665138844068 ,  -4.804000180983039 ,
        11.015589130215647 , -14.536294440523559 ,   0.340769164215163  ,  -1.2916778731810132,
         9.86251476534455  , -15.85451869811148  ,  -1.0064969591691462 ,   0.5306007432103317,
         6.106824647056244 , -12.365469949376868 ,  -1.1618517416174248 ,   1.483802491087325 ,
         3.3707703044033037,  -5.717608387446074 ,  -0.4178099696306948 ,   2.3543421712287853,
         3.061571376151134 ,   2.794128262242935 ,   0.2392938851711605 ,   2.559911274740368 ,
        10.006360719683936 ,   3.614012326480878 ,   3.7702821117827483 ,  -2.0865509895483068,
        20.85647771951526  ,   7.8575664224434805,  12.041493025709583  ,   4.536568089173095 ,
        42.721839753970514 ,  23.62151835484702  ,   0.5824796308117849 ,   4.564753802492409 ,
        37.822615432410274 ,  30.01333829226901  ,  -3.4110481116849485 ,  -0.8744348410102045,
        23.41318721935705  ,  23.678003087136048 ,  -4.908239145989923  ,  -2.7832726450798257,
         6.556193270638787 ,  12.62491304516762  ,  -4.824150848030459  ,  -3.5982318660278167,
        -6.221391453322001 ,   0.1795712145129341,  -2.5529911319415532 ,  -3.5870895900575777,
        -6.6469838755060575,  -8.399365250229572 ,   2.307275232428523  ,  -1.3659616872222047,
         6.943007475132667 ,  -5.3826401468454685,   5.538909932147433  ,   3.1076687377321184,
        11.45843654168509  , -13.283528026318917 ,   6.615531421834127  ,  -7.669248481783164 ,
        32.66540913959966  , -29.427601277140827 ,   5.628319916268867  ,  -1.651536555395904 ,
        41.467478490546014 , -27.90897562449502  ,  -0.5464428070742327 ,   2.5283154847492315,
        33.886497811201274 , -17.631234397628223 ,  -3.830438428866806  ,   3.4055411792436296,
        18.916211310127746 ,  -5.488853302880039 ,  -4.812660512373815  ,   3.604865814412259 ,
         3.5529653905814294,   5.463656340470172 ,  -4.057313654902672  ,   2.718568576477959 ,
        -5.460524326974218 ,   9.08436798995703  ,  -1.1466270598660062 ,  -0.6281497309887072,
         4.54507012208222  ,   3.9981903705619337,   6.923359708430967  ,  -2.3083562867153025,
         7.3788524772573085,  10.92817906853983  ,   4.260182464055043  ,   6.3093871269739035,
        12.864246123266621 ,  17.995394479393287 ,  -1.093188965920503  ,  -2.2291284070965895,
         9.505837418543004 ,   9.455913547831882 ,  -0.8457892038004596 ,  -2.7011432074800634,
         5.524327475257027 ,   1.2563277193687925,  -1.452936633736871  ,  -2.032889877826542 ,
         0.9871375439242183,  -4.862502262358228 ,  -1.166611194482586  ,  -1.4998182592491058,
        -0.6779008638153945,  -7.906348111351482 ,   0.2053008215630547 ,  -0.2575469610388753,
         2.570123739942353 ,  -4.447871166755599 ,   1.6699470577510085 ,   2.2542995559874233,
        10.999618893038576 ,  -0.2716549001527133,   3.1968242380217413 ,   0.156840029729847 ,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[256] = {
        2.7407288610112436e-03, 1.5823605456806854e-03, -7.6809941479703811e-05, -4.4346240389711878e-05, 2.2293997703677062e-03, 1.2871445575538851e-03,
        -2.1840604664782911e-04, -1.2609678982437782e-04, 1.2900949949637090e-03, 7.4483669262332460e-04, -3.2390181828276970e-04, -1.8700480197693243e-04,
        8.9754280741095256e-05, 5.1819658147058668e-05, -3.6911521619292412e-04, -2.1310876943093592e-04, -1.1204747691240249e-03, -6.4690640957332116e-04,
        -3.2961085152780805e-04, -1.9030091385723410e-04, -2.0037116078668389e-03, -1.1568434361802029e-03, -1.8032617507854991e-04, -1.0411136572400634e-04,
        -2.1335600473097673e-03, -1.2318114676469632e-03, 1.0535814361232252e-04, 6.0828552576055904e-05, -9.7553719479136452e-04, -5.6322666201722511e-04,
        5.6322666201691568e-04, 3.2517906493036695e-04, 6.9996983928866515e-03, 8.7655666001799676e-04, -1.9623045690108916e-04, -2.4601226335588913e-05,
        5.6934333929063554e-03, 7.1281618690071866e-04, -5.5794199247052195e-04, -6.9934379895908092e-05, 3.2940753226675251e-03, 4.1216188902718821e-04,
        -8.2732803526416768e-04, -1.0364845991439264e-04, 2.2858332724289904e-04, 2.8333329222096790e-05, -9.4253459349288350e-04, -1.1795506241174320e-04,
        -2.8606716115302335e-03, -3.5779670583284043e-04, -8.4104757700316577e-04, -1.0497721726942451e-04, -5.1120999777412111e-03, -6.3778542591266818e-04,
        -4.5881519628990138e-04, -5.6674345638722714e-05, -5.4353257005495745e-03, -6.7446382105210059e-04, 2.7220073821748091e-04, 3.5498064331527943e-05,
        -2.4819300960496754e-03, -3.0648968502960212e-04, 1.4329430090639938e-03, 1.7695190215588161e-04, 8.9363807954811747e-03, 2.4158744645476202e-04,
        -2.5062337671642466e-04, -6.8025405618047183e-06, 7.2681427235590113e-03, 1.9634266904700807e-04, -7.1253432318021828e-04, -1.9319543853986581e-05,
        4.2043683326866751e-03, 1.1339602534880278e-04, -1.0563366461262076e-03, -2.8569723213548219e-05, 2.9125138296799517e-04, 7.8480896203779337e-06,
        -1.2029024781643599e-03, -3.2368405891667179e-05, -3.6495085404055134e-03, -9.7638507400020033e-05, -1.0722963244070246e-03, -2.8534309293957133e-05,
        -6.5162426599213348e-03, -1.7289672962466951e-04, -5.8281339119051971e-04, -1.4916045566177084e-05, -6.9173360428791955e-03, -1.8117524903371512e-04,
        3.5124201856962114e-04, 1.0136440156872902e-05, -3.1544835104815490e-03, -8.1809209837042127e-05, 1.8212419039310252e-03, 4.7232569321633170e-05,
        8.8731940145174073e-03, -2.7806835145341907e-04, -2.4893599956974554e-04, 7.7767482116643509e-06, 7.2163153814418703e-03, -2.2626519896972250e-04,
        -7.0766332525095505e-04, 2.2131815819668658e-05, 4.1739282906140922e-03, -1.3097059183355206e-04, -1.0488596806175051e-03, 3.2886551262719056e-05,
        2.8950500448299775e-04, -8.8563617087981024e-06, -1.1938131492767329e-03, 3.7616132371691645e-05, -3.6196437585844237e-03, 1.1488094722371538e-04,
        -1.0631349413825358e-03, 3.3823636249294967e-05, -6.4588912082836443e-03, 2.0600860566610596e-04, -5.7610533799726759e-04, 1.8788941883039502e-05,
        -6.8499677692227890e-03, 2.2007033996408154e-04, 3.5031718026529878e-04, -1.0670395800830287e-05, -3.1216003071193998e-03, 1.0079433614999606e-04,
        1.8022567776179350e-03, -5.8193637109094814e-05, 7.2872894092721203e-03, -6.3755409929401423e-04, -2.0452321274956720e-04, 1.7864986214426781e-05,
        5.9261546032858530e-03, -5.1860947359655071e-04, -5.8132833388997338e-04, 5.0807725450663175e-05, 3.4273853696636623e-03, -3.0004616453878808e-04,
        -8.6133675578455374e-04, 7.5379859868811500e-05, 2.3849770473808284e-04, -2.0592716529564491e-05, -9.7976839597570709e-04, 8.5962663565282760e-05,
        -2.9679186382970046e-03, 2.6139272641186739e-04, -8.7145694280961646e-04, 7.6841707824508809e-05, -5.2926770021366433e-03, 4.6730548018563041e-04,
        -4.7074292435403327e-04, 4.2042075996350450e-05, -5.6080254728605207e-03, 4.9696537915865666e-04, 2.8867639989306633e-04, -2.4917925340136650e-05,
        -2.5540116406498441e-03, 2.2690313322518268e-04, 1.4745593082426082e-03, -1.3100258504750940e-04, 4.8119947747803796e-03, -7.9155792425342957e-04,
        -1.3514680468037897e-04, 2.2189501659728380e-05, 3.9127204091128362e-03, -6.4384730040483818e-04, -3.8404949241310673e-04, 6.3091266781423283e-05,
        2.2624464264025014e-03, -3.7253164794252147e-04, -5.6873663574134100e-04, 9.3552898203054172e-05, 1.5799053387988041e-04, -2.5888120237113899e-05,
        -6.4627154030428863e-04, 1.0658183583050446e-04, -1.9549245047614456e-03, 3.2345970927239950e-04, -5.7362052603008966e-04, 9.5114227577624630e-05,
        -3.4818563975266038e-03, 5.7817228334009743e-04, -3.0795401329209934e-04, 5.1944145626343161e-05, -3.6823612630468537e-03, 6.1481737074609174e-04,
        1.9219247521015334e-04, -3.0787094548014825e-05, -1.6747370655752166e-03, 2.8074627938574746e-04, 9.6690989563170287e-04, -1.6208893997735212e-04,
        2.2385331403663865e-03, -6.9423084312466362e-04, -6.2981607343844866e-05, 1.9475094448642045e-05, 1.8195934251495690e-03, -5.6462012723442400e-04,
        -1.7889335001080992e-04, 5.5355687260447954e-05, 1.0513854240675970e-03, -3.2667474776059384e-04, -2.6463174620749096e-04, 8.2022141631195658e-05,
        7.3345062773665792e-05, -2.2981965291741798e-05, -3.0004011966387918e-04, 9.3314968078142519e-05, -9.0535782548232732e-04, 2.8250789554186103e-04,
        -2.6501425632738276e-04, 8.3059685308839399e-05, -1.6067313238332303e-03, 5.0443168272098875e-04, -1.3992392174798064e-04, 4.5068072958611276e-05,
        -1.6884773405719576e-03, 5.3635204862765651e-04, 9.2727836978696821e-05, -2.6638841109765338e-05, -7.6393400782444570e-04, 2.4510611118539448e-04,
        4.4105750506055297e-04, -1.4151207927290955e-04, 5.1804502394653374e-04, -2.9909343402789056e-04, -1.4624877138291715e-05, 8.4436767526642505e-06,
        4.2082133890177203e-04, -2.4296131329567598e-04, -4.1507243593908907e-05, 2.3964218262264415e-05, 2.4278408169825077e-04, -1.4017145492344109e-04,
        -6.1282614778326744e-05, 3.5381534138908897e-05, 1.6769565615290283e-05, -9.6819132221832068e-06, -6.9206926922933511e-05, 3.9956637888743246e-05,
        -2.0801989843223189e-04, 1.2010034435664730e-04, -6.0575330655895149e-05, 3.4973183460431897e-05, -3.6651501020650706e-04, 2.1160753980476531e-04,
        -3.0931864792222451e-05, 1.7858520464326741e-05, -3.7974417080260177e-04, 2.1924539923607339e-04, 2.3294005360914364e-05, -1.3448800265628496e-05,
        -1.6969888500285816e-04, 9.7975696937579516e-05, 9.7975696937579516e-05, -5.6566295000952736e-05,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[256] = {
        2.7455274977581101e-03, -7.7016366281726489e-05, 1.5851310398980775e-03, -4.4465419804975561e-05, 7.0070993653027344e-03, -1.9661720882614239e-04,
        8.7528862499974770e-04, -2.4586158833349093e-05, 8.9378438687396887e-03, -2.5085151156538848e-04, 2.3942723379596830e-04, -6.7260304524665174e-06,
        8.8676660657544103e-03, -2.4887487409620801e-04, -2.7994440724065528e-04, 7.8672426273880338e-06, 7.2760587231604258e-03, -2.0416960868078495e-04,
        -6.3897052045016646e-04, 1.7943354394400559e-05, 4.7967146726859197e-03, -1.3458026370479623e-04, -7.9247943450498712e-04, 2.2234072660216491e-05,
        2.2223314203383862e-03, -6.2382824018365708e-05, -6.9384142923512612e-04, 1.9449138577545980e-05, 5.1028140625373726e-04, -1.4347963919079579e-05,
        -2.9461110726274856e-04, 8.2838008312629473e-06, 2.2328620509273703e-03, -2.1897116745022105e-04, 1.2891435061660986e-03, -1.2642306247195596e-04,
        5.6983929678070212e-03, -5.5896478205819263e-04, 7.1168170157933893e-04, -6.9872342444710309e-05, 7.2683026867494246e-03, -7.1305853927985659e-04,
        1.9470609725546733e-04, -1.9093729767658177e-05, 7.2114189314569262e-03, -7.0735985473499086e-04, -2.2754794868610809e-04, 2.2383866823663154e-05,
        5.9174542262400379e-03, -5.8022106333843716e-04, -5.1952292219273723e-04, 5.1019748613580591e-05, 3.9013337648724092e-03, -3.8236814444858586e-04,
        -6.4448476889657604e-04, 6.3210687367428438e-05, 1.8074675127103821e-03, -1.7713896472772248e-04, -5.6440947543624320e-04, 5.5278434789978199e-05,
        4.1494081249057832e-04, -4.0696953554778999e-05, -2.3956618978940173e-04, 2.3496397089494604e-05, 1.2912653889166565e-03, -3.2465991882971775e-04,
        7.4551241988710028e-04, -1.8744249153107242e-04, 3.2951101213161823e-03, -8.2857121630210346e-04, 4.1140787577797973e-04, -1.0349083171230310e-04,
        4.2030153009268521e-03, -1.0566859580683698e-03, 1.1277142406889629e-04, -2.8211275852572628e-05, 4.1710745676331838e-03, -1.0479837821478085e-03,
        -1.3121241503410969e-04, 3.3235479462844166e-05, 3.4239828477106065e-03, -8.5938530827215753e-04, -3.0012119057253115e-04, 7.5651900198015479e-05,
        2.2586220374435943e-03, -5.6605191358339549e-04, -3.7270018693816884e-04, 9.3704214187848050e-05, 1.0470395366128852e-03, -2.6189433197148961e-04,
        -3.2680729606188115e-04, 8.1901247431851484e-05, 2.4049634777412908e-04, -6.0018605107903177e-05, -1.3885063112723396e-04, 3.4651757815691510e-05,
        8.8453402363655091e-05, -3.6978390539116611e-04, 5.1068595666166369e-05, -2.1349483731995697e-04, 2.2612198365319175e-04, -9.4330991111739261e-04,
        2.8414396800301509e-05, -1.1763055647333072e-04, 2.9025319251122113e-04, -1.2023482985003631e-03, 8.6117738973378721e-06, -3.1925326212673881e-05,
        2.9141683433858939e-04, -1.1919376545334395e-03, -7.9399449750664991e-06, 3.7935914309414511e-05, 2.4327891032320347e-04, -9.7699496619096243e-04,
        -1.9852498413446320e-05, 8.6161304665456563e-05, 1.6455178572971277e-04, -6.4296010994560149e-04, -2.5600628163130763e-05, 1.0669380950652052e-04,
        7.9413826032911092e-05, -2.9676463250666145e-04, -2.3553795786073190e-05, 9.3182242585082505e-05, 1.9308727509790826e-05, -6.7684126993125255e-05,
        -1.1147899025919274e-05, 3.9077448939253196e-05, -1.1232650267591403e-03, -3.2980205584491026e-04, -6.4851736557052430e-04, -1.9041130572132949e-04,
        -2.8634919921787291e-03, -8.4047954952119720e-04, -3.5620314136552121e-04, -1.0442848272175826e-04, -3.6456758182267508e-03, -1.0700613753604627e-03,
        -9.5390901159075347e-05, -2.8120646227587442e-05, -3.6090852016126456e-03, -1.0600182458976896e-03, 1.1651650351104316e-04, 3.3919049726425768e-05,
        -2.9528322367580384e-03, -8.6828066493635644e-04, 2.6237132240425440e-04, 7.6780694255367584e-05, -1.9373743803905745e-03, -5.7058752788066453e-04,
        3.2390354432022150e-04, 9.5092518559640354e-05, -8.8863742870235580e-04, -2.6214002000409995e-04, 2.8158501704607912e-04, 8.2989733143763419e-05,
        -2.0045893626442158e-04, -5.9198792846091458e-05, 1.1573502081347314e-04, 3.4178438985585727e-05, -2.0056573014876579e-03, -1.7964736150044591e-04,
        -1.1579667829159270e-03, -1.0371945252153687e-04, -5.1094538519957291e-03, -4.5622713483378034e-04, -6.3401099102973989e-04, -5.5963954064869424e-05,
        -6.5015813752730047e-03, -5.7879646678034042e-04, -1.6973420928068230e-04, -1.4801482728869664e-05, -6.4365357859305066e-03, -5.7241111007859250e-04,
        2.0728829779717117e-04, 1.8488070139495587e-05, -5.2679135318597775e-03, -4.6833214398650004e-04, 4.6741607517155210e-04, 4.1601948950751952e-05,
        -3.4570580235829475e-03, -3.0680223259952336e-04, 5.7808183999559222e-04, 5.1657388870695265e-05, -1.5835714386286208e-03, -1.3908031769585528e-04,
        5.0357614415094112e-04, 4.5176903847937117e-05, -3.5567598573975194e-04, -3.0415812451309997e-05, 2.0534962611116486e-04, 1.7560577506358637e-05,
        -2.1313761952351938e-03, 1.0706352435308422e-04, -1.2305506200633862e-03, 6.1813154605635436e-05, -5.4220535781736440e-03, 2.7574759877635788e-04,
        -6.6932285279234149e-04, 3.5576641170644565e-05, -6.8908228269317136e-03, 3.5406780988542355e-04, -1.7867146835557504e-04, 9.6415537961626971e-06,
        -6.8192916308551501e-03, 3.5142691997744688e-04, 2.1997002366583540e-04, -1.1166272295433120e-05, -5.5801083789985334e-03, 2.8808636495132524e-04,
        4.9547276036885821e-04, -2.5403414199518600e-05, -3.6595058178618030e-03, 1.8991894407578038e-04, 6.1338764530971795e-04, -3.1273572668628301e-05,
        -1.6691578563090108e-03, 8.9666976409106446e-05, 5.3574028607380976e-04, -2.6606927850515667e-05, -3.7061423058384901e-04, 2.1791212769350009e-05,
        2.1397422579308521e-04, -1.2581162558349135e-05, -9.7296836570905717e-04, 5.6174354785520095e-04, -5.6174354785511985e-04, 3.2432278856973267e-04,
        -2.4722223635138896e-03, 1.4273382470714124e-03, -3.0385115136112275e-04, 1.7542854403186008e-04, -3.1387796954427919e-03, 1.8121753020907926e-03,
        -8.0985903658347883e-05, 4.6757233277722083e-05, -3.1053011751533922e-03, 1.7928464693896641e-03, 1.0031473635950300e-04, -5.7916740040839907e-05,
        -2.5405640789675039e-03, 1.4667953548853823e-03, 2.2573637814478725e-04, -1.3032895868778465e-04, -1.6652782787013596e-03, 9.6144886248386639e-04,
        2.7961011425672833e-04, -1.6143297473426388e-04, -7.5692504870367920e-04, 4.3701088062543986e-04, 2.4482786760169858e-04, -1.4135143526496275e-04,
        -1.6643537145439555e-04, 9.6091506511870639e-05, 9.6091506511870639e-05, -5.5478457151465202e-05,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-7) );
      TEST_MSG("This BC combination is not available");
    }
  } else if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8:
      const double sol[512] = {
        -2.1055183187473453e+01, -5.6960053628737226e+00, 3.8393108916532830e+00, -1.6199865671305258e+00, -5.0628141590584352e-01, 1.6107399622456251e+00,
        -4.8925064333966484e-01, 2.8533276751408687e-01, 5.1858175004515497e+00, -6.5507355610382758e+00, 8.6071381884585936e+00, 1.3050637591200724e+00,
        -1.9423809627420625e+00, -4.8361430357980550e-01, -3.3988181660932626e-01, 4.2364157446457607e-01, 1.9307135289495562e+01, 1.4237688169116351e+00,
        -1.7079115122656434e+00, 2.5311481566680722e+00, -2.0458214044106584e+00, -1.4547416014946468e+00, 2.8016044976686794e-01, -1.7114822825457343e-01,
        1.9294569151566443e-01, 6.4965350024656203e+00, -6.8633995246341426e+00, 1.8492573241138191e-02, -4.8046522865811059e-01, 4.5407487347050113e-01,
        6.2359835968150901e-01, -4.6481501607667941e-01, -9.9523675366444131e+00, 2.7718165640288226e+00, 2.0425310004928692e+00, -1.5705184304155917e+00,
        1.1715022813050320e+00, 1.2569681780583184e+00, 3.3016552688822687e-01, -1.2643505423602566e-03, 7.6579305454496573e+00, -5.0888069607267306e-01,
        5.5928842792949549e+00, 1.6486681682650925e-02, 1.6668374583201608e+00, -7.0422930106180925e-01, -4.4183629099451172e-02, 2.6216025043043689e-01,
        1.1700415434622263e+01, 1.5004199819332511e+00, -4.1739303798805034e+00, 6.5935684087803981e-01, 1.3806005390114717e+00, -1.4129665388092891e+00,
        -1.2107533331542997e-01, -1.1292018871715381e-01, -1.3036693737416885e+01, 5.6308125464531433e-01, -7.3366229431194112e+00, -1.3400430140438562e+00,
        7.5600873308001226e-01, 7.3376873117110675e-01, -2.3953291397273174e-01, -2.2098680881833371e-01, -2.3622486592410571e+01, 8.7552583778672481e+00,
        -6.5794472800424346e+00, -3.2001707637442136e+00, 3.0115387663363169e+00, 1.5213199948472833e+00, 4.3648184555345537e-01, -3.3695940976244687e-01,
        -2.3454873379319217e+01, -3.4914825460446455e+00, 5.8670018000693345e+00, -2.9552064655813575e+00, 3.0874566856601326e+00, 8.9450319136408574e-01,
        -3.9265061439553373e-01, 3.7201493221621273e-01, -4.6581210525804884e-01, -8.7010578192946273e+00, 4.8357086162180192e+00, -4.4645021121903543e-02,
        1.1599456712750884e+00, -1.0962337178653425e+00, -7.2019838862564134e-01, 3.7813285137455588e-01, 8.7046619753618626e-01, -6.7459576998612683e+00,
        -3.1873673791330877e+00, 6.3604115204390721e-01, -9.7666360081639436e-01, -4.1703591022715947e-01, -5.1337354976942640e-01, -3.8121029385150496e-02,
        -7.3850641317072405e+00, -5.8310695790223193e+00, 6.9760538789625404e-01, 9.6657661981126974e-03, -2.3463179009371360e+00, 1.3463881454566500e+00,
        -2.7739672910201180e-01, 5.2890992790714852e-02, 1.0611125333417913e+01, -3.5681337110662690e+00, 8.3330206676841900e+00, 1.6336560247786223e+00,
        -2.8119131812382241e+00, 2.9334041327753013e-01, 8.5851686867588969e-03, 3.1378689267879872e-01, 3.1473362829375780e+01, 5.7768690204496957e+00,
        1.0461332759281459e+00, 3.2351500186680004e+00, -1.8251665366742689e+00, -1.7714744224385899e+00, 5.6111327217419982e-01, -9.4064434402826969e-02,
        1.1973281848365051e+01, 1.3805573956972166e+01, -1.1012655088620420e+01, 6.8550928875883277e-01, 7.0112009639448569e-01, -7.7080769441445973e-01,
        8.9743899547819916e-01, -6.4768079550985869e-01, 1.9965893086268906e+01, 7.7637190920066006e+00, -7.9984011794568355e+00, 2.5942857510310846e+00,
        -3.6862323043438461e+00, 9.5566989841210148e-02, 6.0174080796833174e-01, -4.8619947147572762e-01, 4.6581210525840011e-01, 1.1818723885415153e+01,
        -1.9681206332353183e+00, 4.4645021121926781e-02, -1.1599456712750962e+00, 1.0962337178653341e+00, 8.5681145968406458e-01, -2.5554575843695271e-01,
        2.6185140955145840e+00, 1.1018194245823352e+01, 1.0559679997454556e+00, -2.7512043749361204e-01, 1.5754392195000968e+00, -5.7303427097814530e-01,
        7.2246374334222785e-01, -7.6063509874369384e-02, -5.8445752972248339e+00, 8.7552583778673423e+00, -6.5794472800425243e+00, -1.3682013534830941e+00,
        3.5827918626752542e+00, -1.0666942877560484e+00, 4.3648184555345626e-01, -3.3695940976244909e-01, -3.2809641132170270e+01, -7.0410283489569281e-01,
        -6.2016212882966819e+00, -3.9158361918338302e+00, 3.9617758087657560e+00, 1.0922766148004119e+00, -2.1767536225955730e-01, -1.9960235341928123e-01,
        -3.1473362829376207e+01, -1.4742912684260107e+01, 7.8499625253815193e+00, -3.2351500186680293e+00, 1.8251665366742769e+00, 1.7714744224385992e+00,
        -1.0158965761355097e+00, 5.3961417540868561e-01, 1.0225233950386722e+01, -1.8077810502934248e+01, 1.3144054468008058e+01, 1.5966708782963670e+00,
        -1.8509827239220087e+00, -6.1480933366348278e-01, -1.1065291890510007e+00, 7.6186533476937779e-01, 3.6852126021342571e+01, -5.8310695790223797e+00,
        6.9760538789632731e-01, 4.5587063510291861e+00, -4.2480127280744373e+00, -1.8010138525478823e+00, -2.7739672910201296e-01, 5.2890992790716518e-02,
        -1.9294569151607144e-01, -1.4023246702244311e+01, -5.9570275180648502e-02, -1.8492573241165825e-02, 4.8046522865812202e-01, -4.5407487347048825e-01,
        -9.5341148862819858e-01, 1.6886359373483151e-01, 1.5292239943095500e+00, -1.3085907974956481e+01, 3.1031222880583682e+00, 6.9917874640692146e-01,
        -2.6170745007495793e+00, 1.1332726811086977e+00, -8.3495390797089808e-01, 2.7693021383601302e-01, 2.4281244389184572e+01, -6.5507355610383557e+00,
        8.6071381884586700e+00, 3.2633083581072779e+00, -4.6519496662695241e+00, 2.8988596488549678e-02, -3.3988181660932693e-01, 4.2364157446457562e-01,
        4.1891542471279102e+01, 8.8136714289942564e+00, -9.7172290867058742e-01, 4.8503134702055455e+00, -4.1566144892544088e+00, -1.9322088826315809e+00,
        6.2586371439809296e-01, -1.6274567457648992e-01, 1.3036693737417270e+01, 2.1082862959355285e+01, -1.4140452192572930e+01, 1.3400430140438870e+00,
        -7.5600873308002603e-01, -7.3376873117112851e-01, 1.3374769343369795e+00, -8.5466541862984968e-01, -3.2536774718427473e+01, 1.6010096773801383e+01,
        -8.9849641802045159e+00, -3.8896837439530452e+00, 3.2822953661487722e+00, 1.7344354591952476e+00, 9.9403902442233716e-01, -5.6099863080773993e-01,
        -3.7124992435085780e+01, -5.0888069607260489e-01, 5.5928842792949158e+00, -4.5848587989099991e+00, 4.9274931706914282e+00, 1.1588550081530569e+00,
        -4.4183629099453725e-02, 2.6216025043044233e-01, -1.0883991747161211e+01, -1.1737860227839144e+01, 6.8535648008167245e+00, -1.6598084726594220e+00,
        3.4913936238552177e+00, -9.3549925767235820e-01, -7.8494883084953226e-01, 4.4681409154821727e-01, -1.0883991747161195e+01, 1.1737860227839370e+01,
        -6.8535648008169412e+00, -1.6598084726594260e+00, 3.4913936238552186e+00, -9.3549925767235553e-01, 7.8494883084954370e-01, -4.4681409154822904e-01,
        -3.7124992435086313e+01, 5.0888069607270747e-01, -5.5928842792949904e+00, -4.5848587989100480e+00, 4.9274931706914433e+00, 1.1588550081530800e+00,
        4.4183629099451505e-02, -2.6216025043044144e-01, -3.2536774718428056e+01, -1.6010096773801425e+01, 8.9849641802045515e+00, -3.8896837439530723e+00,
        3.2822953661487810e+00, 1.7344354591952582e+00, -9.9403902442234049e-01, 5.6099863080774659e-01, 1.3036693737416901e+01, -2.1082862959355367e+01,
        1.4140452192573010e+01, 1.3400430140438768e+00, -7.5600873308002070e-01, -7.3376873117111607e-01, -1.3374769343369779e+00, 8.5466541862984524e-01,
        4.1891542471279116e+01, -8.8136714289943683e+00, 9.7172290867070843e-01, 4.8503134702055508e+00, -4.1566144892544088e+00, -1.9322088826315902e+00,
        -6.2586371439809696e-01, 1.6274567457649522e-01, 2.4281244389184916e+01, 6.5507355610382714e+00, -8.6071381884585936e+00, 3.2633083581072806e+00,
        -4.6519496662695277e+00, 2.8988596488547458e-02, 3.3988181660932781e-01, -4.2364157446457673e-01, 1.5292239943101598e+00, 1.3085907974956429e+01,
        -3.1031222880583029e+00, 6.9917874640694755e-01, -2.6170745007495921e+00, 1.1332726811086875e+00, 8.3495390797089153e-01, -2.7693021383600680e-01,
        -1.9294569151546725e-01, 1.4023246702244418e+01, 5.9570275180559962e-02, -1.8492573241108392e-02, 4.8046522865810043e-01, -4.5407487347051267e-01,
        9.5341148862820035e-01, -1.6886359373483195e-01, 3.6852126021343210e+01, 5.8310695790224907e+00, -6.9760538789643034e-01, 4.5587063510292465e+00,
        -4.2480127280744560e+00, -1.8010138525479076e+00, 2.7739672910201563e-01, -5.2890992790721292e-02, 1.0225233950386755e+01, 1.8077810502934451e+01,
        -1.3144054468008258e+01, 1.5966708782963559e+00, -1.8509827239220089e+00, -6.1480933366347124e-01, 1.1065291890510089e+00, -7.6186533476938623e-01,
        -3.1473362829376690e+01, 1.4742912684260212e+01, -7.8499625253816108e+00, -3.2351500186680662e+00, 1.8251665366742946e+00, 1.7714744224386161e+00,
        1.0158965761355114e+00, -5.3961417540868917e-01, -3.2809641132170846e+01, 7.0410283489566261e-01, 6.2016212882967103e+00, -3.9158361918338600e+00,
        3.9617758087657693e+00, 1.0922766148004239e+00, 2.1767536225955242e-01, 1.9960235341928756e-01, -5.8445752972251830e+00, -8.7552583778674418e+00,
        6.5794472800426256e+00, -1.3682013534831063e+00, 3.5827918626752551e+00, -1.0666942877560404e+00, -4.3648184555345848e-01, 3.3695940976245353e-01,
        2.6185140955146124e+00, -1.1018194245823448e+01, -1.0559679997453570e+00, -2.7512043749360404e-01, 1.5754392195000948e+00, -5.7303427097815152e-01,
        -7.2246374334222718e-01, 7.6063509874368496e-02, 4.6581210525872974e-01, -1.1818723885415245e+01, 1.9681206332354051e+00, 4.4645021121922049e-02,
        -1.1599456712750968e+00, 1.0962337178653363e+00, -8.5681145968406436e-01, 2.5554575843695337e-01, 1.9965893086269553e+01, -7.7637190920066717e+00,
        7.9984011794569154e+00, 2.5942857510311130e+00, -3.6862323043438567e+00, 9.5566989841200378e-02, -6.0174080796833818e-01, 4.8619947147573656e-01,
        1.1973281848365744e+01, -1.3805573956972246e+01, 1.1012655088620505e+01, 6.8550928875886941e-01, 7.0112009639446882e-01, -7.7080769441447927e-01,
        -8.9743899547820727e-01, 6.4768079550986690e-01, 3.1473362829376391e+01, -5.7768690204495927e+00, -1.0461332759282804e+00, 3.2351500186680466e+00,
        -1.8251665366742902e+00, -1.7714744224386076e+00, -5.6111327217419404e-01, 9.4064434402818808e-02, 1.0611125333417949e+01, 3.5681337110664368e+00,
        -8.3330206676843588e+00, 1.6336560247786172e+00, -2.8119131812382241e+00, 2.9334041327753324e-01, -8.5851686867530752e-03, -3.1378689267880400e-01,
        -7.3850641317077024e+00, 5.8310695790224178e+00, -6.9760538789635351e-01, 9.6657661980786500e-03, -2.3463179009371222e+00, 1.3463881454566655e+00,
        2.7739672910201407e-01, -5.2890992790719071e-02, 8.7046619753558563e-01, 6.7459576998612381e+00, 3.1873673791331085e+00, 6.3604115204387313e-01,
        -9.7666360081638093e-01, -4.1703591022714370e-01, 5.1337354976942440e-01, 3.8121029385151051e-02, -4.6581210525837857e-01, 8.7010578192945207e+00,
        -4.8357086162179010e+00, -4.4645021121907096e-02, 1.1599456712750904e+00, -1.0962337178653403e+00, 7.2019838862563668e-01, -3.7813285137454988e-01,
        -2.3454873379319199e+01, 3.4914825460445620e+00, -5.8670018000692536e+00, -2.9552064655813530e+00, 3.0874566856601318e+00, 8.9450319136408307e-01,
        3.9265061439553617e-01, -3.7201493221621629e-01, -2.3622486592410262e+01, -8.7552583778673565e+00, 6.5794472800425332e+00, -3.2001707637442247e+00,
        3.0115387663363178e+00, 1.5213199948472917e+00, -4.3648184555345715e-01, 3.3695940976244865e-01, -1.3036693737416577e+01, -5.6308125464543024e-01,
        7.3366229431195347e+00, -1.3400430140438684e+00, 7.5600873308001537e-01, 7.3376873117111230e-01, 2.3953291397272691e-01, 2.2098680881833999e-01,
        1.1700415434622929e+01, -1.5004199819333131e+00, 4.1739303798805487e+00, 6.5935684087808211e-01, 1.3806005390114564e+00, -1.4129665388093109e+00,
        1.2107533331542447e-01, 1.1292018871715899e-01, 7.6579305454501920e+00, 5.0888069607277087e-01, -5.5928842792950695e+00, 1.6486681682687045e-02,
        1.6668374583201444e+00, -7.0422930106182391e-01, 4.4183629099455723e-02, -2.6216025043044228e-01, -9.9523675366444024e+00, -2.7718165640286627e+00,
        -2.0425310004930370e+00, -1.5705184304155926e+00, 1.1715022813050326e+00, 1.2569681780583182e+00, -3.3016552688822121e-01, 1.2643505423549639e-03,
        1.9294569151517549e-01, -6.4965350024655244e+00, 6.8633995246340387e+00, 1.8492573241105727e-02, -4.8046522865809577e-01, 4.5407487347051401e-01,
        -6.2359835968150668e-01, 4.6481501607667730e-01, 1.9307135289494944e+01, -1.4237688169116791e+00, 1.7079115122656794e+00, 2.5311481566680309e+00,
        -2.0458214044106442e+00, -1.4547416014946284e+00, -2.8016044976687038e-01, 1.7114822825457676e-01, 5.1858175004512272e+00, 6.5507355610381603e+00,
        -8.6071381884584799e+00, 1.3050637591200762e+00, -1.9423809627420634e+00, -4.8361430357980417e-01, 3.3988181660931938e-01, -4.2364157446456985e-01,
        -2.1055183187473450e+01, 5.6960053628736542e+00, -3.8393108916532168e+00, -1.6199865671305218e+00, -5.0628141590584441e-01, 1.6107399622456211e+00,
        4.8925064333967161e-01, -2.8533276751409536e-01,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
              (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[512] = {
        1.9225897306558122e+02, 9.8696543926233048e+01, 1.0099327346116128e+02, 5.1204645467852352e+01, -9.5308124208910101e+00, -7.7517733729732239e+00,
        -5.5026171167895948e+00, -4.4754884435817042e+00, 4.4002429450586322e+02, 2.2465172671550962e+02, 4.4092504640193397e+01, 2.2692485783685946e+01,
        -2.2770367600114788e+01, -6.1728331428662724e+00, -2.1412436298861719e+00, -3.5638868766862246e+00, 5.1479752556826054e+02, 2.6302301520217674e+02,
        1.3096196553750585e+00, 7.4970628292885422e-01, -2.6487700854411869e+01, -4.4441043975949936e+00, -4.9597251497860732e-03, -2.5658048702558278e+00,
        4.6564915159799369e+02, 2.3793642874117083e+02, -2.7477709662214245e+01, -1.3958816215443175e+01, -2.3939866133519576e+01, -2.7339987210843586e+00,
        1.4759527871076226e+00, -1.5784748975789984e+00, 3.4090726873506662e+02, 1.7422577999436439e+02, -4.2341277783958972e+01, -2.1553927564394165e+01,
        -1.7503656503446397e+01, -1.0292879021388917e+00, 2.2399945753759849e+00, -5.9425964737300341e-01, 1.8885763536205437e+02, 9.6523834483826846e+01,
        -4.3237014899218686e+01, -2.2032750785757997e+01, -9.6926614542966671e+00, 6.8071486129409275e-01, 2.2696855188894580e+00, 3.9301090840859615e-01,
        5.7656197772677913e+01, 2.9560696015497196e+01, -3.0323877274044818e+01, -1.5365023602984083e+01, -2.8870188347116104e+00, 2.3757594365606702e+00,
        1.6595540795363681e+00, 1.3716453502277908e+00, -3.9614453538742946e+00, -2.2952663188533133e+00, -3.0155181372963034e+00, -1.7363193658871419e+00,
        -6.2934255232094897e-03, 4.1074225252697261e+00, 3.6335109199499257e-03, 2.3714215006400075e+00, 4.1237707730420999e+02, 3.1423441173067523e+01,
        2.1530543759641580e+02, 1.6545620278432789e+01, -7.1801141924131535e+00, -1.7645760611223963e+01, -4.1454408618002923e+00, -1.2368077517810760e+00,
        9.4035261607490986e+02, 7.1616355146582080e+01, 9.4415300291442733e+01, 7.1315294353256338e+00, -1.7035725806528042e+01, -1.3854947677946351e+01,
        -1.5446991566373500e+00, -8.7138401808781685e-01, 1.1004324789476625e+03, 8.3789825667683985e+01, 3.0122036698565835e+00, 2.0963173523695869e-01,
        -1.9751542617693776e+01, -9.9778186594763714e+00, -2.3278410358856350e-02, -6.2908654845986633e-01, 9.9552492353974560e+02, 7.5822866592913570e+01,
        -5.8587959426371349e+01, -4.4761251397385058e+00, -1.7870610162569566e+01, -6.1105289099653648e+00, 1.1092352696522505e+00, -3.7096571589837640e-01,
        7.2889812853782894e+02, 5.5523529929299663e+01, -9.0358376374595977e+01, -6.9035439316519689e+00, -1.3055393709560402e+01, -2.2447368676808148e+00,
        1.6708312456989478e+00, -1.0748014006887091e-01, 4.0380469128737400e+02, 3.0767135046854325e+01, -9.2341034714787938e+01, -7.0557640214800372e+00,
        -7.2207461538046589e+00, 1.6231127763177817e+00, 1.6978040912432262e+00, 1.5108278151407334e-01, 1.2348873758767802e+02, 9.3896048365970657e+00,
        -6.4521415290046576e+01, -4.9393085662844278e+00, -2.1574523602860203e+00, 5.4790870035551107e+00, 1.2254899434309243e+00, 4.2006165596002820e-01,
        -8.9597955316263551e+00, -6.0489782500546929e-01, -6.9241557519082599e+00, -5.1203978984010667e-01, -1.7420757049447295e-02, 9.3703703573804944e+00,
        1.0057878772022792e-02, 6.6714284696010184e-01, 4.3312460696061254e+02, -1.6329942592688518e+01, 2.2626137357931032e+02, -8.4218154398872986e+00,
        -4.7673178002322407e+00, -1.8437906425823950e+01, -2.7524122152727362e+00, 7.7946215248432915e-01, 9.8742492719370910e+02, -3.6997350288963204e+01,
        9.9038382508956204e+01, -3.7626547710595455e+00, -1.1271337455292487e+01, -1.4352438053006937e+01, -1.0026852833909299e+00, 5.8415781612733753e-01,
        1.1554243309778467e+03, -4.3372668467236323e+01, 3.1476816693852494e+00, -1.2343253176582421e-01, -1.3037651139637806e+01, -1.0331070505720698e+01,
        -1.7096397739202965e-02, 4.2513649993896863e-01, 1.0453143331037525e+03, -3.9230622574556143e+01, -6.1516018215434102e+01, 2.3034614991152207e+00,
        -1.1792880694817788e+01, -6.3003015336895656e+00, 7.3576494913538548e-01, 2.6140044050642253e-01, 7.6537795102827488e+02, -2.8719053066087557e+01,
        -9.4892763923262279e+01, 3.5533211931059641e+00, -8.6070077897084261e+00, -2.2626672792163265e+00, 1.1035996302334381e+00, 9.7128012142172301e-02,
        4.2402771143586085e+02, -1.5907654183146764e+01, -9.6978027225107240e+01, 3.6333754723337215e+00, -4.7546513225920846e+00, 1.7713249454239182e+00,
        1.1205594130705743e+00, -6.5512445783472192e-02, 1.2963668872001253e+02, -4.8781509178848657e+00, -6.7776483880235958e+01, 2.5224670443672474e+00,
        -1.4123425122644739e+00, 5.8086971523804172e+00, 8.0912347828693332e-01, -2.2976114780810064e-01, -9.3963126923480988e+00, 3.6832900369737825e-01,
        -7.2841445136128051e+00, 2.9527753378829469e-01, -5.4497691372753345e-03, 9.8444328621049895e+00, 3.1464256784867707e-03, -3.9344273224470699e-01,
        3.2235376963293595e+02, -4.4494858976233203e+01, 1.6844698012137096e+02, -2.3150943302758790e+01, -2.3437820649419501e+00, -1.3682518594711942e+01,
        -1.3531832061170055e+00, 1.9660622919092583e+00, 7.3508173189064439e+02, -1.0125236951682614e+02, 7.3719237809870862e+01, -1.0172268514182901e+01,
        -5.5078022689936388e+00, -1.0678928747005607e+01, -4.7356471041367887e-01, 1.5367437705631224e+00, 8.6007117816312916e+02, -1.1849395147585753e+02,
        2.3274595486845953e+00, -3.3178580065447261e-01, -6.3330503964170308e+00, -7.6703899333200800e+00, -2.8925181024576953e-03, 1.1110081447641107e+00,
        7.7808245883342590e+02, -1.0720609301985657e+02, -4.5784822953317146e+01, 6.2953111048651751e+00, -5.7125455242868899e+00, -4.6658489186918173e+00,
        3.6114117306029364e-01, 6.8225121674022804e-01, 5.6972002350277933e+02, -7.8501356723755492e+01, -7.0631276227126421e+01, 9.7212549390903220e+00,
        -4.1586948308806937e+00, -1.6588577186035227e+00, 5.3597494305825177e-01, 2.5148160021688282e-01, 3.1562735900273304e+02, -4.3491146307054201e+01,
        -7.2186283511650160e+01, 9.9369053925263824e+00, -2.2866529741594457e+00, 1.3488302797737264e+00, 5.4484892685402475e-01, -1.7841496316083455e-01,
        9.6480436626325343e+01, -1.3310160655960640e+01, -5.0457360744270169e+01, 6.9314924915758214e+00, -6.7302136437066695e-01, 4.3549550342039485e+00,
        3.8678171743042289e-01, -6.0955725545338146e-01, -7.0533281398645284e+00, 9.8943106376197598e-01, -5.4339340435613721e+00, 7.7003368953857587e-01,
        -1.5478891712049066e-03, 7.3634514486903964e+00, 8.9367422968142883e-04, -1.0389525546446661e+00, 1.4781932454085802e+02, -5.3152436916607961e+01,
        7.7285775564560254e+01, -2.7679546893108320e+01, 7.3124504243893246e-02, -6.2415072802381246e+00, 4.2218452210944883e-02, 2.3300075935452469e+00,
        3.3719888074644251e+02, -1.2102197156840140e+02, 3.3831173009631321e+01, -1.2162561636437315e+01, 2.5785391651682188e-01, -4.8638012632779883e+00,
        6.4435123692078530e-02, 1.8206216475390413e+00, 3.9454240506876926e+02, -1.4161675895271290e+02, 1.0611236475995243e+00, -3.9008574745726804e-01,
        3.7685609578913409e-01, -3.4811481812282281e+00, 4.2708165449212854e-03, 1.3076517085065784e+00, 3.5690835858470405e+02, -1.2811427338762630e+02,
        -2.1009632112402379e+01, 7.5272250164294334e+00, 3.6386615548210294e-01, -2.1027391530487489e+00, -1.1770562077925548e-02, 7.9756089641636574e-01,
        2.6131316758286533e+02, -9.3814028899129042e+01, -3.2404608609185154e+01, 1.1617210053621342e+01, 2.9012031935336247e-01, -7.2576302561570505e-01,
        -3.0806616262630482e-02, 2.8724087195906034e-01, 1.4474419004910749e+02, -5.1981366770667925e+01, -3.3116217755997077e+01, 1.1875943595976237e+01,
        1.8193278990071829e-01, 6.5325467259024350e-01, -3.1655482989800067e-02, -2.2317580088828998e-01, 4.4218492015554460e+01, -1.5909172559360719e+01,
        -2.3143526912924656e+01, 8.2933039302731828e+00, 6.6032606786869341e-02, 2.0311026973954731e+00, -3.5259518930106808e-02, -7.3211951675993081e-01,
        -3.2939264042909402e+00, 1.1862613803985360e+00, -2.5040868312824576e+00, 9.1851168070307543e-01, 2.4806642747430114e-03, 3.4127459895780490e+00,
        -1.4322121867907276e-03, -1.2419883056627714e+00, -2.2905469721302193e+01, -4.2294386366573058e+01, -1.1923860184495274e+01, -2.2023644682704976e+01,
        2.4907821823826710e+00, 1.0074550872667674e+00, 1.4380537634913269e+00, 1.8551827806791803e+00, -5.2132073344659631e+01, -9.6314262014631410e+01,
        -5.2193085785319635e+00, -9.6894239636891779e+00, 6.0240270961214888e+00, 7.9894815383500928e-01, 6.0186613856865634e-01, 1.4487682527845545e+00,
        -6.0985102507581033e+01, -1.1271850163655034e+02, -1.6632645823316813e-01, -3.0828373173151635e-01, 7.0877048850350741e+00, 5.8644819562257022e-01,
        1.2248519191625565e-02, 1.0407761546229803e+00, -5.5192999986062240e+01, -1.0196748230317775e+02, 3.2419906133774545e+00, 5.9915413028102469e+00,
        6.4406348210024218e+00, 3.7854041174995884e-01, -3.8583459484542287e-01, 6.3500652825488857e-01, -4.0441816349118305e+01, -7.4663730728369558e+01,
        5.0019542301563984e+00, 9.2473524916175407e+00, 4.7367401948801247e+00, 1.6734275680394184e-01, -5.9790942621704324e-01, 2.2839399193574725e-01,
        -2.2464622052456324e+01, -4.1366425667536916e+01, 5.1077972430478011e+00, 9.4524111702174594e+00, 2.6530977454391556e+00, -4.1323123604449410e-02,
        -6.0508210286263142e-01, -1.7783887671785187e-01, -6.9566606605275565e+00, -1.2681611910133302e+01, 3.5810948250901689e+00, 6.5860437146631012e+00,
        8.0593021424560563e-01, -2.4630744118503087e-01, -4.6138056851030101e-01, -5.8274383980470368e-01, 4.1918422440863951e-01, 9.5868807579065729e-01,
        3.7665830958924651e-01, 7.4400369881679473e-01, 3.3978139803747258e-03, -4.7922315622460576e-01, -1.9617288162267353e-03, -1.0050411283441087e+00,
        -1.2240397998182672e+02, -1.2045347867162372e+01, -6.3842602653506603e+01, -6.1585453584632672e+00, 4.8964915579624364e+00, 5.2884566200641938e+00,
        2.8269907190748351e+00, 6.1645460668261620e-01, -2.7882058521764810e+02, -2.7131398799817458e+01, -2.7986548205290813e+01, -2.7500159726458508e+00,
        1.1781613850189567e+01, 4.1095346489752336e+00, 1.1481364897459991e+00, 4.6259975136019243e-01, -3.2629540518948824e+02, -3.1790829279440789e+01,
        -8.9719028933284273e-01, -1.0609802238104514e-01, 1.3801632868750907e+01, 2.9676719511176075e+00, 1.8122034388553452e-02, 3.3402402161281725e-01,
        -2.9520438209736915e+02, -2.8757935761158393e+01, 1.7352737906438872e+01, 1.6728920416476478e+00, 1.2517514558079052e+01, 1.8097561010737333e+00,
        -7.5950808672627546e-01, 1.9130623524461565e-01, -2.1618383094082404e+02, -2.1059077935631542e+01, 2.6778814297943534e+01, 2.5935773873711301e+00,
        9.1843896592678664e+00, 6.5484464593809988e-01, -1.1648724708450036e+00, 5.3065354986315823e-02, -1.1983010684958269e+02, -1.1662017462424414e+01,
        2.7358298074208111e+01, 2.6499611199416573e+00, 5.1206108556290539e+00, -5.0142381889631171e-01, -1.1813513153629656e+00, -8.7800383563236981e-02,
        -3.6762367412702822e+01, -3.5725733473969647e+00, 1.9101945607910103e+01, 1.8305856840783836e+00, 1.5450097355863150e+00, -1.6624962342005629e+00,
        -8.8302295380844087e-01, -2.3489314106614415e-01, 2.5362522967644598e+00, 2.6926500742843040e-01, 2.1345452616294498e+00, 2.6764312045172245e-01,
        7.7845576702786929e-03, -2.7876581587344296e+00, -4.4944164664543979e-03, -3.2773444176169236e-01, -8.2574478177290928e+01, 3.8196987619970166e+01,
        -4.3571495931617385e+01, 1.9684229930639169e+01, 7.3411698967120085e+00, 3.1780936596991651e+00, 4.2384264160337555e+00, -1.8348732299371819e+00,
        -1.8901691596793813e+02, 8.6449270326545843e+01, -1.8816079997681733e+01, 8.7129096390028806e+00, 1.7567617747035836e+01, 2.4553904609498884e+00,
        1.6658160032043772e+00, -1.4176203435950809e+00, -2.2110473228560673e+02, 1.0117986894194146e+02, -5.9065062909415778e-01, 3.0034781582304987e-01,
        2.0507452830784278e+01, 1.7731092637493990e+00, 3.1498573770896779e-02, -1.0237051107287320e+00, -2.0009426246694153e+02, 9.1517111712290941e+01,
        1.1814056231487582e+01, -5.3554896096866162e+00, 1.8596024645811884e+01, 1.0705541101609184e+00, -1.1350621509013501e+00, -6.1808470368363944e-01,
        -1.4655150069660706e+02, 6.7007937429309848e+01, 1.8198821548326929e+01, -8.2752445687609306e+00, 1.3635698406959255e+01, 3.7337826844797695e-01,
        -1.7287835383685790e+00, -2.1557004379807332e-01, -8.1234074284841000e+01, 3.7117640860149550e+01, 1.8608351405056347e+01, -8.4600819437576860e+00,
        7.5777814576609375e+00, -3.2674927207593346e-01, -1.7687564430372047e+00, 1.8864878019055645e-01, -2.4889329326410696e+01, 1.1401368538642078e+01,
        1.3008378540357970e+01, -5.8595606956889412e+00, 2.2993843628198807e+00, -1.0346715444382852e+00, -1.2787275405590182e+00, 5.9736789470428842e-01,
        1.4154752488564883e+00, -8.7181038721819559e-01, 1.3486188331646116e+00, -7.4711056757121341e-01, 4.2281646767034310e-02, -1.6776554316279519e+00,
        -2.4411320142727996e-02, 9.6859481505783585e-01,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8:
      const double sol[512] = {
        7.5711738636379700e+00, 2.5568324740391297e+00, 5.1167331614558211e+00, 2.4727832928346594e+00, -1.4054179675134042e+00, 6.1419700713567691e+00,
        -3.7286293196278525e-01, 3.5460680740521866e+00, 2.0758894207673666e+01, 1.1267795520680716e+01, -1.0031441046941417e+01, -4.6769211349736635e+00,
        -5.5566283720313869e-01, -3.5626741185393747e+00, 8.6346928519268240e-01, -2.0569108613736473e+00, 1.3594715095713351e+01, 7.1655480313952662e+00,
        9.0698668172746899e+00, 4.1413820063944256e+00, -5.2933177934331999e-01, -1.1035880146021995e+00, -8.4826704184926838e-01, -6.3715683730558981e-01,
        1.0538674015081105e+01, 4.2560390203120635e+00, -2.7952876149664085e+00, -1.1798774654372435e+00, -1.4163246487905283e+00, 5.1233832560628052e+00,
        3.3616146990464424e-01, 2.9579867020494177e+00, 2.2024747870440134e+01, 1.2472953612918014e+01, -5.1167331614558895e+00, -2.4727832928346856e+00,
        -1.8825836094986273e-01, -6.1419700713565017e+00, 3.7286293196286846e-01, -3.5460680740524686e+00, 8.8370275264048548e+00, 3.7619905662762716e+00,
        1.0031441046941339e+01, 4.6769211349739148e+00, -1.0380134912601533e+00, 3.5626741185393436e+00, -8.6346928519278077e-01, 2.0569108613736105e+00,
        1.6001206638363609e+01, 7.8642380555622040e+00, -9.0698668172749777e+00, -4.1413820063947195e+00, -1.0643445491199175e+00, 1.1035880146025827e+00,
        8.4826704184939872e-01, 6.3715683730530448e-01, 1.9057247718996599e+01, 1.0773747066644622e+01, 2.7952876149668433e+00, 1.1798774654373116e+00,
        -1.7735167967268334e-01, -5.1233832560626400e+00, -3.3616146990475915e-01, -2.9579867020495270e+00, 1.0547053185519411e+01, 5.6121235077038900e-01,
        8.0357928214499879e+00, -4.5108902511093207e-01, -3.2104452298288144e-01, 1.3775209244937894e+01, -1.1230845198426297e-01, 8.6098461760367073e-01,
        4.3173120529816174e+01, 1.3313117726991164e+00, -1.5775484413258026e+01, 1.8666732294070940e-02, -8.2038303621019348e-01, -6.8107509670996524e+00,
        -1.7598477304427174e-01, 1.8163281850877949e-01, 3.4608182327383652e+01, 5.0049655083846334e+00, 1.4274111188784909e+01, 4.2469027913556817e-01,
        -4.9960076306725493e-01, -4.1433528573200746e+00, 3.6118850479460851e-01, -1.1178522129107362e+00, 1.4094766749853889e+01, -9.6046485001291515e-01,
        -4.4111572207434300e+00, -6.1926948485557964e-01, -4.5391689108755867e-01, 1.2670336771618857e+01, -3.3481290900950705e-01, 1.3992489417184601e+00,
        5.1670055472623496e+01, 3.8022739325395261e+00, -8.0357928214498031e+00, 4.5108902511095250e-01, -8.3930591556586931e-01, -1.3775209244937987e+01,
        1.1230845198422829e-01, -8.6098461760359934e-01, 1.9043988128327072e+01, 3.0321745106108899e+00, 1.5775484413258118e+01, -1.8666732294083995e-02,
        -3.3996740233860129e-01, 6.8107509670996196e+00, 1.7598477304428253e-01, -1.8163281850874280e-01, 2.7608926330759374e+01, -6.4147922507440869e-01,
        -1.4274111188785204e+01, -4.2469027913545798e-01, -6.6074967548153596e-01, 4.1433528573199894e+00, -3.6118850479461617e-01, 1.1178522129107504e+00,
        4.8122341908288419e+01, 5.3239511333231038e+00, 4.4111572207434486e+00, 6.1926948485546185e-01, -7.0643354746118825e-01, -1.2670336771618837e+01,
        3.3481290900953997e-01, -1.3992489417184355e+00, 1.2834086265285496e+01, 6.6660597862113313e-01, 5.4592934411306917e+00, -5.3613267888178739e-01,
        -3.9277282123210222e-01, 1.3651686305630548e+01, 2.7523845597608154e-01, -9.3230061986385970e-01, 3.8079683584953059e+01, -3.4752185565637603e+00,
        -1.5203332332988003e+01, 3.2483479157023204e-01, -2.0319184059141193e-01, -5.1861881128342349e+00, -1.6578382576994696e-01, 7.5630898271682545e-01,
        4.2060635144023266e+01, -7.8152890802934238e-01, 1.6041465337446521e+01, 7.6746911112531901e-02, -5.6097905910266699e-01, -6.3173087394422609e+00,
        -4.0784721150132874e-02, -1.3728180083887725e-01, 1.1185122138368454e+01, -4.4915680662569779e-01, -7.4827255075668804e+00, -4.3337131413575919e-01,
        -2.4457250288101104e-01, 1.4120211809851538e+01, 2.2346213155806885e-01, -5.6216319810350202e-01, 5.0868739817121380e+01, -3.6128915887194353e+00,
        -5.4592934411306455e+00, 5.3613267888173377e-01, -3.3425172740217313e-01, -1.3651686305630527e+01, -2.7523845597607949e-01, 9.3230061986385515e-01,
        2.5623142497454019e+01, 5.2893294646533351e-01, 1.5203332332988076e+01, -3.2483479157024980e-01, -5.2383270804286930e-01, 5.1861881128342597e+00,
        1.6578382576994177e-01, -7.5630898271682989e-01, 2.1642190938384100e+01, -2.1647567020690461e+00, -1.6041465337446517e+01, -7.6746911112476141e-02,
        -1.6604548953164397e-01, 6.3173087394422121e+00, 4.0784721150120259e-02, 1.3728180083888500e-01, 5.2517703944038566e+01, -2.4971288034724965e+00,
        7.4827255075667578e+00, 4.3337131413577545e-01, -4.8245204575329603e-01, -1.4120211809851520e+01, -2.2346213155805197e-01, 5.6216319810347881e-01,
        1.1389883870962576e+01, -1.3868221174744275e+00, 4.4760604474563914e+00, -4.1951753830011562e-01, -3.0478289812702114e-01, 9.5116421105830735e+00,
        -2.5289856897649846e-02, -1.4579550106037478e+00, 2.4658979651651553e+01, -4.0728807035622605e+00, -1.2842002664550789e+01, 1.4765469333527865e+00,
        -4.8003126436937500e-02, -2.9027088989601997e+00, 1.7354172720537125e-01, 5.6205835610226373e-01, 3.1718420671373405e+01, -4.6314885092856359e+00,
        1.3685273888782918e+01, -1.6686351603278100e+00, -1.2870540135553579e-01, -5.4065918180524992e+00, -2.2013520735383832e-01, 6.6308446055880155e-01,
        8.4657676578212353e+00, -1.1554391882963218e+00, -6.5118972737563556e+00, 8.8325954103543369e-01, -2.7135492134138733e-01, 1.0548784374265418e+01,
        1.3777646859023468e-01, -1.4998013932233591e+00, 3.4290543891020292e+01, -5.5127074757936700e+00, -4.4760604474564500e+00, 4.1951753830011534e-01,
        1.1084239407202070e-02, -9.5116421105830771e+00, 2.5289856897656619e-02, 1.4579550106037384e+00, 2.1021448110331136e+01, -2.8266488897059010e+00,
        1.2842002664550764e+01, -1.4765469333528143e+00, -2.4569553228286417e-01, 2.9027088989602130e+00, -1.7354172720536756e-01, -5.6205835610226673e-01,
        1.3962007090609443e+01, -2.2680410839826166e+00, -1.3685273888782840e+01, 1.6686351603278042e+00, -1.6499325736427489e-01, 5.4065918180524877e+00,
        2.2013520735383085e-01, -6.6308446055878967e-01, 3.7214660104161752e+01, -5.7440904049718569e+00, 6.5118972737563618e+00, -8.8325954103540005e-01,
        -2.2343737378438464e-02, -1.0548784374265438e+01, -1.3777646859023721e-01, 1.4998013932233598e+00, 5.2864798717669776e+00, -1.5541361622449834e+00,
        1.9459080267195821e+00, -1.0078353910637479e+00, 1.4668586129324565e-01, 4.0836858666066682e+00, 6.0587942565820945e-04, -1.6758769880055386e+00,
        9.0769379769636984e+00, -4.8410921221688952e+00, -5.8399824469108212e+00, 2.2436644081508579e+00, 1.5764095609959970e-02, -6.2779560686890912e-01,
        -7.6193596085684684e-02, 7.5136344546970124e-01, 1.5638698717011067e+01, -4.3964175120957218e+00, 6.3130743537226452e+00, -2.1651852443570010e+00,
        6.9378934107321086e-02, -3.1958488049744180e+00, 1.0714813752469497e-01, 6.1328861315091376e-01, 2.5685095801924365e+00, -1.7383264165802581e+00,
        -3.0880529243934101e+00, 8.1836992946892340e-01, 1.2447786804319297e-01, 5.1474083301575693e+00, -7.5336753184758354e-02, -1.6186845198367033e+00,
        1.4490787580215658e+01, -5.9421095039544918e+00, -1.9459080267195892e+00, 1.0078353910637696e+00, -7.0586300985872683e-03, -4.0836858666066851e+00,
        -6.0587942565855347e-04, 1.6758769880055411e+00, 1.0700329475018780e+01, -2.6551535440305218e+00, 5.8399824469107600e+00, -2.2436644081508517e+00,
        1.2386313558470298e-01, 6.2779560686891089e-01, 7.6193596085687210e-02, -7.5136344546970324e-01, 4.1385687349712983e+00, -3.0998281541037285e+00,
        -6.3130743537226417e+00, 2.1651852443569783e+00, 7.0248297087352390e-02, 3.1958488049744251e+00, -1.0714813752469037e-01, -6.1328861315091387e-01,
        1.7208757871790073e+01, -5.7579192496192562e+00, 3.0880529243934753e+00, -8.1836992946892850e-01, 1.5149363151478083e-02, -5.1474083301575719e+00,
        7.5336753184751998e-02, 1.6186845198367106e+00, 2.0851644586908202e+00, -1.5289639520143722e-01, -1.0190853525359438e+00, -4.4200250843767624e-01,
        2.5609590421969947e-01, -1.3111543673389618e+00, 2.0355162130222118e-01, -1.4388354732979942e+00, -5.2199298577201665e+00, -2.9026745735353474e+00,
        3.8255692666894958e-01, 1.3001259421911353e+00, 4.1122187447099579e-01, 1.4627303288154059e+00, -1.1398960061933343e-01, 4.5560226624553257e-01,
        1.2352101112731999e+00, -3.8277001636926316e+00, 4.7806815846104928e-01, -1.3966532318021432e+00, 1.4044070210636231e-01, -7.5745870176624008e-01,
        -4.2345982136870186e-02, 7.9451656932564751e-01, -5.8864206348251669e-01, 2.3026174978385772e-01, -1.0586474001031971e+00, 6.7504000015566690e-01,
        3.6825713824841150e-01, -3.9152195984007238e-01, 1.7387586286930351e-01, -1.5792183741160062e+00, -4.4644655511723803e+00, -4.5835374336907790e+00,
        1.0190853525359791e+00, 4.4200250843767314e-01, 3.1685721688944901e-01, 1.3111543673389592e+00, -2.0355162130222612e-01, 1.4388354732980002e+00,
        2.8406287652386411e+00, -1.8337592553568307e+00, -3.8255692666896990e-01, -1.3001259421911211e+00, 1.6173124663814459e-01, -1.4627303288154128e+00,
        1.1398960061933365e-01, -4.5560226624553479e-01, -3.6145112037548346e+00, -9.0873366519952370e-01, -4.7806815846107620e-01, 1.3966532318021447e+00,
        4.3251241900278042e-01, 7.5745870176624530e-01, 4.2345982136871282e-02, -7.9451656932564774e-01, -1.7906590289991435e+00, -4.9666955786760303e+00,
        1.0586474001032093e+00, -6.7504000015567955e-01, 2.0469598286073923e-01, 3.9152195984007748e-01, -1.7387586286929990e-01, 1.5792183741160046e+00,
        1.9980167436569440e+00, 2.6279956734987525e-01, -9.1824150266946936e-01, 7.6348475138740574e-02, 3.8020017828485891e-01, -4.0973051386771191e+00,
        -1.2478151566869267e-01, -1.6974942453697403e-01, -7.0245946979708958e+00, 2.0099206798419127e+00, 3.1949717001217537e+00, 7.4112012303711705e-01,
        5.2676833384458188e-01, 1.7012769926403721e+00, 2.0940267973572577e-01, -3.1787728567190099e-01, -7.7041224392599732e+00, -7.4960990094331026e-01,
        -3.6001308070409044e+00, -1.1244506044854634e+00, 5.9266262324759322e-01, 1.6913361423317976e+00, -1.7135859401083814e-01, 6.1929579310452265e-01,
        2.2794863501078102e+00, 1.4058345596943498e+00, 1.8963821135127572e+00, 8.4909317204483359e-01, 3.5290586993119349e-01, -4.0931875036577843e+00,
        3.2934967943583970e-02, -5.5793922405711172e-01, -1.1159940859954022e+01, 1.1171063513037944e+00, 9.1824150266947813e-01, -7.6348475138740976e-02,
        6.2607883273876874e-01, 4.0973051386771173e+00, 1.2478151566869694e-01, 1.6974942453696854e-01, -2.1373294183261322e+00, -6.3001476118826361e-01,
        -3.1949717001217364e+00, -7.4112012303711539e-01, 4.7951067717904555e-01, -1.7012769926403750e+00, -2.0940267973573007e-01, 3.1787728567190571e-01,
        -1.4578016770370641e+00, 2.1295158195969872e+00, 3.6001308070408862e+00, 1.1244506044854694e+00, 4.1361638777603105e-01, -1.6913361423317970e+00,
        1.7135859401084064e-01, -6.1929579310452432e-01, -1.1441410466404900e+01, -2.5928641040681811e-02, -1.8963821135127661e+00, -8.4909317204484103e-01,
        6.5337314109243083e-01, 4.0931875036577887e+00, -3.2934967943586656e-02, 5.5793922405711294e-01, 7.5087282814798053e+00, 3.8201328875905691e+00,
        -7.4192744296018398e-01, 3.0744537382085513e-01, 1.0783609380268158e+00, -2.1956598832653493e+00, 9.3653876830995159e-02, 1.2676648246521331e+00,
        5.9367957061953227e+00, 5.4565665660586538e+00, 3.4117284350768085e+00, -1.4280377956224548e+00, 5.1377448012734850e-01, 5.7534869164227787e-01,
        -4.1961802028007000e-01, -3.3217772199756690e-01, 2.9704312624447171e+00, 6.9899600377174531e+00, -4.0829851810593842e+00, 1.7121050443297883e+00,
        6.5261151279673846e-01, 1.3819939604512403e+00, 4.9977561846523150e-01, -7.9789458508496924e-01, 8.7374366650226367e+00, 3.1849805151751132e+00,
        2.3624845829457257e+00, -9.9324437827614753e-01, 1.0208527561355243e+00, -2.5297832936301354e+00, -2.8717143749686075e-01, 1.4605710655687689e+00,
        3.5480238541686964e+00, 7.0326406888475326e+00, 7.4192744296022595e-01, -3.0744537382083909e-01, 3.6124396291129257e-01, 2.1956598832653560e+00,
        -9.3653876830993660e-02, -1.2676648246521234e+00, 5.1199564294532207e+00, 5.3962070103794693e+00, -3.4117284350768302e+00, 1.4280377956224359e+00,
        9.2583042081076283e-01, -5.7534869164227420e-01, 4.1961802028007011e-01, 3.3217772199756546e-01, 8.0863208732038245e+00, 3.8628135387206504e+00,
        4.0829851810593940e+00, -1.7121050443297823e+00, 7.8699338814137187e-01, -1.3819939604512455e+00, -4.9977561846523177e-01, 7.9789458508496636e-01,
        2.3193154706258663e+00, 7.6677930612630050e+00, -2.3624845829457559e+00, 9.9324437827614442e-01, 4.1875214480258355e-01, 2.5297832936301310e+00,
        2.8717143749685969e-01, -1.4605710655687754e+00,        
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[512] = {
        7.2565136902344483e+00, 4.8231040378017838e+00, 2.4660327201331449e+00, 2.2399474225899145e+00, 6.1343565305018393e+00, -1.3350308469482448e+00,
        3.5416723941906159e+00, -4.2190188540145057e-01, 1.1949813331139318e+01, 6.9702547633583389e+00, 1.5068482148949991e+00, -2.4695883553587722e-01,
        1.3436531385394451e+01, -3.5655791455018548e-01, 6.7424022395134853e-01, 1.6162543991040532e-01, 1.4716264229951848e+01, 6.6338750059206353e+00,
        -1.8320658369255069e-01, -7.0972051070675959e-02, 1.3403839983501085e+01, -5.6846396898344476e-01, -6.9311461363466975e-01, 6.5790834873097587e-02,
        1.0549934521247797e+01, 5.3526143020031389e+00, -1.6852857041035738e+00, -7.7274514171037867e-01, 9.7561428635705152e+00, -1.5255258254589169e-01,
        -1.4128843004794656e+00, -1.4752486265131429e-02, 4.5889680258859959e+00, 2.1817893301442188e+00, -1.4933184581198531e+00, -9.7285807268042979e-01,
        4.3495411425699286e+00, 5.1136479604261972e-02, -1.7086186585413166e+00, 5.1145114416602425e-02, 7.6872761317398586e-01, -8.8753538357070305e-01,
        -3.3212657227395292e-01, -6.9330142199440248e-01, -1.0482760888626736e+00, 3.4561628982098297e-01, -1.4078125730626960e+00, 1.3318726039373399e-01,
        1.5370523987374454e+00, -1.6273654716803474e+00, 8.0187784881376467e-01, 1.2612521825998468e-01, -4.0200641079134964e+00, 3.6587877413708753e-01,
        -3.0795003971080093e-01, 2.4716535870589681e-02, 7.7578733640628128e+00, -8.4543282675032450e-01, 3.6929071177981427e+00, 3.9076288214178517e-01,
        -2.2767246114428010e+00, 1.0654885273666856e+00, 1.3144675672870636e+00, 7.5405427234410613e-02, 2.0268549938791278e+01, -9.4969060696722671e+00,
        1.0888323728929288e+01, -4.4131017188939587e+00, -3.3051144715720029e+00, -6.3031179964622264e-01, -1.9082087298650074e+00, 8.2877161706433056e-01,
        4.0767021158311643e+01, -1.5658193196918768e+01, 9.1291966939315738e-01, -6.4613295094788326e-01, -6.1805639124758018e+00, -6.5630298673172449e-01,
        2.4806722111804072e-01, -3.3468333802268019e-01, 3.7758585266463847e+01, -1.7342738885110744e+01, -2.0416669776853111e+00, 3.0635683668654923e-01,
        -5.1960495480828603e+00, -1.8521346012573825e-01, 3.2034241218527731e-01, 1.5547894958296166e-01, 2.6668131237388220e+01, -1.3172813714297360e+01,
        -4.3163372608041159e+00, 1.9211520705163905e+00, -3.5096342603975730e+00, -1.5030054421297667e-01, 6.5330990812532364e-01, 1.6052701202864460e-02,
        1.0437380332720705e+01, -5.9052182837064375e+00, -4.7197787440845147e+00, 2.2311481016558434e+00, -1.1501040187224223e+00, 1.0896656526118718e-01,
        7.0896551206688718e-01, -1.7756898895317105e-02, -3.7255517594536443e+00, 1.0076669738259358e+00, -3.1947882140381423e+00, 1.5860998166900273e+00,
        1.0546802782338107e+00, 3.1222514477371105e-01, 5.6396729528585388e-01, -1.5246564697532630e-01, -7.1988691712637394e+00, 3.9804950282612692e+00,
        1.5774257623995511e+00, 4.7918342366710176e-01, 1.7094851042336434e+00, 6.1273634482641992e-01, -1.8591555272818444e-01, 1.1780674901839330e-01,
        5.5715694778914706e+00, 3.4180165071055630e+00, 5.6676306193408292e+00, -1.4647055793739878e+00, 6.9373496049558103e-01, 5.1361910218894236e-01,
        -4.0052806618838704e-01, -3.9402738841830454e-01, 1.4602826784400772e+01, 8.6075493265117853e+00, 7.7930019397339709e+00, 4.0011208804020413e+00,
        -1.4602188196092223e+00, -4.9414932480429263e-01, -8.4305772857675976e-01, -7.5015817556081665e-01, 3.6608160545361521e+01, 1.5173774417982225e+01,
        4.6510254311452552e+00, 1.1607288178625641e+00, -4.6958940772574875e+00, -6.9613160278726083e-01, -1.0250602524366783e+00, 3.1168827582156533e-01,
        4.0632558776218822e+01, 1.7892461534098196e+01, -1.9590644309267153e+00, -3.6228194229717742e-01, -6.0555162418397277e+00, -4.1071318083536745e-01,
        2.4008202971894158e-01, -2.8567127403685122e-01, 2.9717000599675977e+01, 1.3276557507370477e+01, -3.9887253575562287e+00, -1.9441741717949959e+00,
        -4.7927704935470272e+00, -1.3626532132207336e-01, 4.8896456797590365e-01, -7.9494614886657690e-03, 1.4412254532553726e+01, 6.1694504554472314e+00,
        -4.6287982408052128e+00, -2.1824618323441993e+00, -2.7230484411529030e+00, 3.3120119168853060e-02, 7.0599001612154855e-01, -2.6033067173168618e-02,
        4.3827718290769346e-01, -5.3752091716913175e-01, -3.2353589144670192e+00, -1.5497824500461481e+00, -4.4326706458301257e-01, 1.9092283427357970e-01,
        6.1024237533457426e-01, 8.2431725354751781e-02, -6.9966967706349799e+00, -4.0019045822453716e+00, -6.7704800450272062e-01, -8.0379291487445303e-01,
        1.6024870888315665e+00, 4.8540690033464212e-01, 5.7087433583512792e-01, -1.9132043807148655e-01, 3.2377940652632744e+00, -3.9883724740136244e+00,
        6.8186961608289733e+00, 1.6806436130924369e+00, 1.2956352216174452e+00, 6.6570366104462586e-01, -7.4803534397254068e-01, 4.8183344941321854e-01,
        9.6033330614339949e+00, -2.6760069266761515e+00, 3.7481569851264216e+00, -1.2453376948648369e+00, 5.3701757302955393e+00, -1.3914311907140546e+00,
        3.1004724034819504e+00, 2.3211224873884562e-01, 1.3672469801041744e+01, -5.8007643773818742e+00, -4.1525889207611438e-02, -9.9538548551070161e-01,
        1.2821561004000882e+01, -3.4006036160939240e-01, 1.2015868901274416e+00, -1.0611044887684108e-01, 1.3525825424309970e+01, -7.9610228806500452e+00,
        -2.1742167884647040e-01, 2.0598719951303407e-01, 1.3759842744463120e+01, -4.7505892635419916e-01, -6.5986967469585844e-01, 2.4852124054034622e-01,
        9.2870514814841503e+00, -5.6030739742524212e+00, -1.8209869976236110e+00, 8.2832541085158895e-01, 1.0287635293713352e+01, -1.5836616221819666e-01,
        -1.3448102316767432e+00, -4.8104649520278034e-03, 2.9425212235886828e+00, -2.8197022227758892e+00, -1.5310038164894428e+00, 8.5531902100701973e-01,
        5.0010760551997526e+00, 8.2553106233491144e-02, -1.7073861677659867e+00, 5.4573215561784238e-02, -9.5598680615284481e-01, -2.4749760270591306e-01,
        -3.1532163792134243e-01, 6.0562354289307252e-01, -4.2780598374722406e-01, 3.9586135197729472e-01, -1.4269803388187741e+00, 3.5889583008811080e-02,
        1.4533098484599591e+00, 1.6790527072733459e+00, 1.7357114590788620e+00, 6.5755141788770677e-01, -3.9757440787869918e+00, 4.1862035693500366e-01,
        -6.2142267542056351e-01, 1.5276120926146078e-01, 8.7245547915061294e+00, 2.2223939374401298e+00, 3.2161201593333195e+00, -9.1208341177703689e-01,
        -2.5260398627953622e+00, 1.0024930404611774e+00, 1.4584097947686661e+00, -2.8738801054686564e-01, 2.2339408043844461e+01, -4.8231040378014391e+00,
        1.2563753366823573e+01, -2.2399474225897751e+00, -6.1343565305019432e+00, -2.5864548151508043e-01, -3.5416723941905062e+00, 4.2190188540134388e-01,
        5.0267295327003161e+01, -6.9702547633581187e+00, 2.8566380684147634e+00, 2.4695883553579936e-01, -1.3436531385394398e+01, -8.0379252399852275e-01,
        -6.7424022395136907e-01, -1.6162543991040920e-01, 4.8986561852454606e+01, -6.6338750059207348e+00, -2.7630790264056375e+00, 7.0972051070600173e-02,
        -1.3403839983501079e+01, -1.5856057965081735e-01, 6.9311461363466331e-01, -6.5790834873076742e-02, 3.5130493240734950e+01, -5.3526143020032064e+00,
        -5.2142438891645240e+00, 7.7274514171045117e-01, -9.7561428635705383e+00, -1.4114607617396283e-01, 1.4128843004794569e+00, 1.4752486265134764e-02,
        1.5188299426096288e+01, -2.1817893301441087e+00, -6.0029272080797043e+00, 9.7285807268043278e-01, -4.3495411425699322e+00, 8.8490751590401054e-02,
        1.7086186585413363e+00, -5.1145114416619598e-02, -3.1480287056559360e+00, 8.8753538357069084e-01, -4.4043072566181865e+00, 6.9330142199436606e-01,
        1.0482760888626790e+00, 2.2733683128818291e-01, 1.4078125730626825e+00, -1.3318726039372250e-01, -1.0698976515034619e+01, 1.6273654716803183e+00,
        5.7802806983996791e-01, -1.2612521825998788e-01, 4.0200641079135284e+00, 6.4040023688653513e-01, 3.0795003971083024e-01, -2.4716535870601318e-02,
        3.2988787715858732e+00, 8.4543282675016818e-01, 7.1598664586400513e+00, -3.9076288214185306e-01, 2.2767246114427619e+00, 3.7411637357141198e-01,
        -1.3144675672871342e+00, -7.5405427234420203e-02, 9.3273717952862238e+00, 9.4969060696714962e+00, 4.1414623580284928e+00, 4.4131017188939197e+00,
        3.3051144715721992e+00, -9.6336452881698631e-01, 1.9082087298647199e+00, -8.2877161706415614e-01, 2.1450087499831831e+01, 1.5658193196918884e+01,
        3.4505666139165925e+00, 6.4613295094810441e-01, 6.1805639124756686e+00, -5.0404745181712973e-01, -2.4806722111794385e-01, 3.3468333802259947e-01,
        2.5944240815942624e+01, 1.7342738885110894e+01, -9.0461863241320528e-01, -3.0635683668665153e-01, 5.1960495480828959e+00, -5.4181108850849746e-01,
        -3.2034241218527509e-01, -1.5547894958296579e-01, 1.9012296524594099e+01, 1.3172813714297234e+01, -2.5831923324638790e+00, -1.9211520705164102e+00,
        3.5096342603975881e+00, -1.4339811450682732e-01, -6.5330990812533796e-01, -1.6052701202840355e-02, 9.3398871192617534e+00, 5.9052182837064073e+00,
        -2.7764669221148899e+00, -2.2311481016557915e+00, 1.1501040187224039e+00, 3.0660665933457315e-02, -7.0896551206689318e-01, 1.7756898895323451e-02,
        1.3462506669718795e+00, -1.0076669738258186e+00, -1.5416456148541455e+00, -1.5860998166900340e+00, -1.0546802782338074e+00, 2.6072797633543321e-01,
        -5.6396729528583400e-01, 1.5246564697530224e-01, -1.9630549450335482e+00, -3.9804950282613318e+00, -1.9751984374577769e-01, -4.7918342366713162e-01,
        -1.7094851042336310e+00, 3.9354266619723177e-01, 1.8591555272816998e-01, -1.1780674901836495e-01, 5.4851826577571474e+00, -3.4180165071054893e+00,
        5.1851429570973018e+00, 1.4647055793740646e+00, -6.9373496049565408e-01, 9.2598579874914499e-01, 4.0052806618835141e-01, 3.9402738841830831e-01,
        1.4993094949676593e+01, -8.6075493265116485e+00, 7.2367841472226750e+00, -4.0011208804019649e+00, 1.4602188196089905e+00, -1.0995270036588776e+00,
        8.4305772857700834e-01, 7.5015817556066489e-01, 2.5608948112781171e+01, -1.5173774417982559e+01, -2.8753914783475915e-01, -1.1607288178625927e+00,
        4.6958940772575373e+00, -4.6421883576156686e-01, 1.0250602524365906e+00, -3.1168827582146935e-01, 2.3070267306188306e+01, -1.7892461534098082e+01,
        -9.8722117917186647e-01, 3.6228194229727734e-01, 6.0555162418396753e+00, -3.1631136779895441e-01, -2.4008202971891188e-01, 2.8567127403680437e-01,
        1.5963427162306223e+01, -1.3276557507370432e+01, -2.9108042357119785e+00, 1.9441741717949204e+00, 4.7927704935470388e+00, -1.5743333739768528e-01,
        -4.8896456797589682e-01, 7.9494614886675818e-03, 5.3650129194284357e+00, -6.1694504554473397e+00, -2.8674474253941060e+00, 2.1824618323441931e+00,
        2.7230484411529106e+00, 1.0650711202581434e-01, -7.0599001612155710e-01, 2.6033067173175637e-02, -2.8175782753892253e+00, 5.3752091716910833e-01,
        -1.5010749144251061e+00, 1.5497824500462098e+00, 4.4326706458298470e-01, 3.8203028683554430e-01, -6.1024237533458470e-01, -8.2431725354739485e-02,
        -2.1652273456622262e+00, 4.0019045822455137e+00, 2.0569539231562413e+00, 8.0379291487441784e-01, -1.6024870888315290e+00, 5.2087211068898720e-01,
        -5.7087433583507952e-01, 1.9132043807144530e-01, 7.8189580703851567e+00, 3.9883724740134583e+00, 4.0340774156093291e+00, -1.6806436130925051e+00,
        -1.2956352216175056e+00, 7.7390123989349768e-01, 7.4803534397243487e-01, -4.8183344941320183e-01, 1.9992588672643731e+01, 2.6760069266764401e+00,
        1.1281629101830370e+01, 1.2453376948646588e+00, -5.3701757302956015e+00, -2.0224513774923147e-01, -3.1004724034818452e+00, -2.3211224873876074e-01,
        4.8544638857100331e+01, 5.8007643773818707e+00, 4.4050121725177291e+00, 9.9538548551058659e-01, -1.2821561004000813e+01, -8.2029007693928679e-01,
        -1.2015868901274718e+00, 1.0611044887682913e-01, 5.0177000658096993e+01, 7.9610228806498817e+00, -2.7288639312516918e+00, -2.0598719951295591e-01,
        -1.3759842744463137e+01, -2.5196562228015207e-01, 6.5986967469583846e-01, -2.4852124054031519e-01, 3.6393376280498451e+01, 5.6030739742525686e+00,
        -5.0785425956447172e+00, -8.2832541085156619e-01, -1.0287635293713372e+01, -1.3533249650161161e-01, 1.3448102316767636e+00, 4.8104649519981908e-03,
        1.6834746228393318e+01, 2.8197022227759163e+00, -5.9652418497099902e+00, -8.5531902100706858e-01, -5.0010760551997357e+00, 5.7074124961195356e-02,
        1.7073861677659865e+00, -5.4573215561780575e-02, -1.4233142863289567e+00, 2.4749760270583140e-01, -4.4211121909707272e+00, -6.0562354289309139e-01,
        4.2780598374723527e-01, 1.7709176913185090e-01, 1.4269803388187716e+00, -3.5889583008810913e-02, -1.0615233964756841e+01, -1.6790527072733965e+00,
        -3.5580554042511986e-01, -6.5755141788763871e-01, 3.9757440787869553e+00, 5.8765865408859663e-01, 6.2142267542053842e-01, -1.5276120926143633e-01,
        2.3321973441423172e+00, -2.2223939374398811e+00, 7.6366534171046352e+00, 9.1208341177709618e-01, 2.5260398627953577e+00, 4.3711186047695644e-01,
        -1.4584097947686228e+00, 2.8738801054685470e-01,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[512] = {
        2.9454144508192573e-03, 1.5802962885446991e-03, -8.1220932092995003e-05, -4.3533810263315391e-05, -9.3137145266585083e-05, -2.0139117638263735e-05,
        2.6019606427078908e-06, -1.1627324989882401e-05, 2.3977549337022169e-03, 1.2864677832761426e-03, -2.3080949477184936e-04, -1.2370583492221198e-04,
        -7.5814945013522543e-05, -1.6916065230087864e-05, 7.3990163363543645e-06, -9.7664948141897960e-06, 1.3914848579091601e-03, 7.4658885947776258e-04,
        -3.4183873168162134e-04, -1.8318865615912528e-04, -4.3985702544433978e-05, -1.0469815448560873e-05, 1.0977605371275791e-05, -6.0447507678323708e-06,
        1.0470809822782509e-04, 5.6205957086642959e-05, -3.8859985852632078e-04, -2.0819742850286847e-04, -3.2899379366231808e-06, -8.0105320620244238e-07,
        1.2518105279921458e-05, -4.6248828444866835e-07, -1.1954480557012073e-03, -6.4132020832850215e-04, -3.4539864687539403e-04, -1.8490834090177961e-04,
        3.7856125934796355e-05, 1.2093657591401790e-05, 1.1237585772336225e-05, 6.9822764657059424e-06, -2.1526562979653596e-03, -1.1550179533269416e-03,
        -1.8645540951630360e-04, -9.9671858995086963e-05, 6.8024103326402879e-05, 2.8197816652966905e-05, 6.1799040956142929e-06, 1.6280017034490287e-05,
        -2.3223780125947308e-03, -1.2453849494040894e-03, 1.1351008767306690e-04, 6.1957353160490669e-05, 7.3927987392683565e-05, 4.7596471550266739e-05,
        -2.7712950406826957e-06, 2.7479835663115355e-05, -1.1659088106909560e-03, -6.2851585708619104e-04, 5.8293435125386509e-04, 3.1079483366484276e-04,
        3.4563981789703322e-05, 6.9871264386467942e-05, -1.9955524190579951e-05, 4.0340193302214889e-05, 7.1907956051621390e-03, 8.7835443511798811e-04,
        -1.9820546892075793e-04, -2.4219105429680350e-05, -8.7266654190746009e-05, -5.0742376387900417e-05, 2.4378531481555240e-06, -6.0414746873002307e-06,
        5.8537814702808158e-03, 7.1503760511787125e-04, -5.6323747400306148e-04, -6.8824462963351944e-05, -7.1037591261220273e-05, -4.2622645299582466e-05,
        6.9320007029016667e-06, -5.0752061088776479e-06, 3.3971515566350342e-03, 4.1495965265860577e-04, -8.3413313207031533e-04, -1.0193078316635971e-04,
        -4.1215509592931909e-05, -2.6382965965691905e-05, 1.0285786176079738e-05, -3.1427109668882989e-06, 2.5570802352229950e-04, 3.1230035269336657e-05,
        -9.4813318581610641e-04, -1.1587304851261703e-04, -3.0915289970245487e-06, -2.0242619983592091e-06, 1.1725104283547521e-05, -2.4373164097814988e-07,
        -2.9183643369499615e-03, -3.5648089485149369e-04, -8.4248417026401253e-04, -1.0298355368059811e-04, 3.5474299728766200e-05, 3.0457810799828251e-05,
        1.0540887312808836e-05, 3.6202723326186516e-06, -5.2554104848072606e-03, -6.4196098795661009e-04, -4.5441272355273087e-04, -5.5593505979902678e-05,
        6.3684056561191883e-05, 7.1044607132342305e-05, 5.7460233881661060e-06, 8.4575889826887134e-06, -5.6685699947737729e-03, -6.9244248203468358e-04,
        2.7884442364894062e-04, 3.3885625879175091e-05, 6.9355861018852180e-05, 1.1982053676960119e-04, -2.4714055577450204e-06, 1.4218747833236671e-05,
        -2.8533343234986215e-03, -3.4833572521995930e-04, 1.4196091704694591e-03, 1.7376995150007819e-04, 3.2537630513163073e-05, 1.7642536690424146e-04,
        -1.8785609735800377e-05, 2.1178846469680725e-05, 9.1097460852529299e-03, 2.4472708031813889e-04, -2.5110830227391921e-04, -6.7482221991316338e-06,
        -7.5512196835018694e-05, -6.4112513957384137e-05, 2.1095306275704184e-06, -1.6777778375435766e-06, 7.4159318202933104e-03, 1.9922344402943450e-04,
        -7.1357220381875082e-04, -1.9176646944203308e-05, -6.1468902515097800e-05, -5.3854053317244781e-05, 5.9983691283449032e-06, -1.4092503334982851e-06,
        4.3037179995065785e-03, 1.1561467930872094e-04, -1.0567779427295731e-03, -2.8401605422880185e-05, -3.5663237742205277e-05, -3.3337047466052149e-05,
        8.9005383749018049e-06, -8.7222985931139863e-07, 3.2393510846188414e-04, 8.7001219740704852e-06, -1.2012200053436138e-03, -3.2285441357682172e-05,
        -2.6738421603053111e-06, -2.5617211843318249e-06, 1.0145898044711315e-05, -6.6570564721603515e-08, -3.6971709813012699e-03, -9.9330926542001054e-05,
        -1.0673912986844585e-03, -2.8700079352181221e-05, 3.0696900696076708e-05, 3.8472831201701916e-05, 9.1207093264786159e-06, 1.0072018539634889e-06,
        -6.6579376458280673e-03, -1.7885207668121318e-04, -5.7579088410443090e-04, -1.5477492943621758e-05, 5.5114205089303011e-05, 8.9764490935925406e-05,
        4.9766279378358783e-06, 2.3503409705086686e-06, -7.1812728008930288e-03, -1.9299951105228473e-04, 3.5310329741445771e-04, 9.3927259295886449e-06,
        5.9996703214613967e-05, 1.5131402435257719e-04, -2.1577163308699145e-06, 3.9640257005147641e-06, -3.6146488029613457e-03, -9.6899989605516777e-05,
        1.7988304609254128e-03, 4.8459029277639108e-05, 2.8129714448702223e-05, 2.2315200407948573e-04, -1.6240698211392054e-05, 5.7987900817725341e-06,
        9.0177565111321835e-03, -2.7507798332316246e-04, -2.4857112966730646e-04, 7.5772432740715380e-06, -5.7882917020132727e-05, -6.3496060694986082e-05,
        1.6170295071426932e-06, 2.0336872945321556e-06, 7.3410454860559750e-03, -2.2393242527812523e-04, -7.0636214180566301e-04, 2.1531479484430446e-05,
        -4.7118210302593887e-05, -5.3336225987068894e-05, 4.5979768139757526e-06, 1.7082180819688677e-06, 4.2602569175799734e-03, -1.2995818786303053e-04,
        -1.0460993447132367e-03, 3.1884219912884969e-05, -2.7337284980894984e-05, -3.3016523735266036e-05, 6.8225457453270383e-06, 1.0572843215624519e-06,
        3.2065845040763071e-04, -9.7859520015234814e-06, -1.1890799666174055e-03, 3.6236689242385506e-05, -2.0495564532662268e-06, -2.5369781191878596e-06,
        7.7773311272936553e-06, 8.0855980043051561e-08, -3.6598484720956432e-03, 1.1162641418891941e-04, -1.0566012181598438e-03, 3.2180193597608133e-05,
        2.3529810939151314e-05, 3.8102668699248937e-05, 6.9909235224191676e-06, -1.2209152743986281e-06, -6.5907198832517557e-03, 2.0104898149736030e-04,
        -5.6995716625598430e-04, 1.7346634574953829e-05, 4.2247156344146251e-05, 8.8901503535101169e-05, 3.8155408856701057e-06, -2.8485869786818794e-06,
        -7.1088121777639970e-03, 2.1675429654359875e-04, 3.4954002903375870e-04, -1.0799200306602247e-05, 4.5991705133030206e-05, 1.4986504102866444e-04,
        -1.6536246344143250e-06, -4.8005966126270943e-06, -3.5781972670616492e-03, 1.0946869793234542e-04, 1.7806298296942254e-03, -5.4073203785811124e-05,
        2.1563771624132259e-05, 2.2094852997700201e-04, -1.2449849352694297e-05, -7.0709664479939801e-06, 7.3879809458021967e-03, -6.3552715512547452e-04,
        -2.0364783787258774e-04, 1.7511480650208876e-05, -3.4376833827426706e-05, -5.2008579772633332e-05, 9.6036464919927713e-07, 4.5986129082988774e-06,
        6.0142998836003221e-03, -5.1736180879796451e-04, -5.7870422787173416e-04, 4.9761310767244505e-05, -2.7983583757900193e-05, -4.3687037174877011e-05,
        2.7307799994384960e-06, 3.8627436762118259e-06, 3.4902979044166751e-03, -3.0024557864690533e-04, -8.5704280104083373e-04, 7.3690802618229373e-05,
        -1.6235550950274673e-05, -2.7043903484178800e-05, 4.0519499044925716e-06, 2.3910095881701911e-06, 2.6269656745679838e-04, -2.2603494454847577e-05,
        -9.7418454781392769e-04, 8.3755890534553272e-05, -1.2169717205572190e-06, -2.0793097886905229e-06, 4.6190308566307569e-06, 1.8337895376917601e-07,
        -2.9984238446050604e-03, 2.5791176006065172e-04, -8.6565107141790289e-04, 7.4399789671366034e-05, 1.3974764058163703e-05, 3.1207532438979965e-05,
        4.1519218846715914e-06, -2.7599935015668652e-06, -5.3996095123790885e-03, 4.6449045764734116e-04, -4.6696068084749697e-04, 4.0118044862483829e-05,
        2.5091037559222535e-05, 7.2813783690334763e-05, 2.2660616135503387e-06, -6.4396624043412238e-06, -5.8240722780158504e-03, 5.0087577604735996e-04,
        2.8636488989313314e-04, -2.4809915859630460e-05, 2.7312419260886122e-05, 1.2275040511432561e-04, -9.8354629012194350e-07, -1.0854045731495190e-05,
        -2.9315113588201957e-03, 2.5258732845322358e-04, 1.4588473580084335e-03, -1.2517916774314479e-04, 1.2804433557384362e-05, 1.8099343058216964e-04,
        -7.3926431611940980e-06, -1.5997120943110683e-05, 4.8513053211664955e-03, -7.9109106616121524e-04, -1.3372380484153113e-04, 2.1799490283123503e-05,
        -4.9947296048720255e-06, -3.4177209945446541e-05, 1.3953357513057844e-07, 5.6963332614470956e-06, 3.9492789724229925e-03, -6.4400082304249954e-04,
        -3.8000161541640629e-04, 6.1946532533867433e-05, -4.0657836928335455e-06, -2.8708762724645629e-05, 3.9679359724754937e-07, 4.7849671096256083e-06,
        2.2918948796507767e-03, -3.7373806306923346e-04, -5.6277007512186953e-04, 9.1736560736739719e-05, -2.3589093832583682e-06, -1.7771874817657828e-05,
        5.8867074485852490e-07, 2.9621986583796967e-06, 1.7249106136722970e-04, -2.8133448726912747e-05, -6.3968902507017535e-04, 1.0426874001204558e-04,
        -1.7651020443375028e-07, -1.3664868346019459e-06, 6.7133800851515766e-07, 2.2816957065840111e-07, -1.9689279834608363e-03, 3.2104735935276969e-04,
        -5.6842056737624140e-04, 9.2624173569731219e-05, 2.0300239681481928e-06, 2.0507557287362720e-05, 6.0260509000118129e-07, -3.4176400325416493e-06,
        -3.5456687095027318e-03, 5.7819743192140838e-04, -3.0661400038094287e-04, 4.9959252250043811e-05, 3.6460715791936043e-06, 4.7848584895897408e-05,
        3.3042043325918573e-07, -7.9740018399998324e-06, -3.8244430614088800e-03, 6.2347426592679919e-04, 1.8803933619758128e-04, -3.0874605919218922e-05,
        3.9688452100248159e-06, 8.0672725212157189e-05, -1.4406699061151716e-07, -1.3439514086896728e-05, -1.9250048499072851e-03, 3.1438982954188488e-04,
        9.5792792719695454e-04, -1.5586904563329022e-04, 1.8596569313441248e-06, 1.1888055850170400e-04, -1.0736734298563878e-06, -1.9863762472686977e-05,
        2.1963203666646599e-03, -6.9624461964085159e-04, -6.0547835005102099e-05, 1.9177068354589424e-05, 3.0265369629908206e-05, -1.5338588673617838e-05,
        -8.4544916186461430e-07, 5.1801498010046700e-06, 1.7879462173578853e-03, -5.6678994487324719e-04, -1.7205947860693760e-04, 5.4493073013620596e-05,
        2.4637112795788217e-05, -1.2884967417993727e-05, -2.4040264363828245e-06, 4.3509053702713410e-06, 1.0375990840386848e-03, -3.2893020251966809e-04,
        -2.5481861647818450e-04, 8.0694150272559692e-05, 1.4294725086052851e-05, -7.9778786783124563e-06, -3.5671538918963875e-06, 2.6923676491135237e-06,
        7.8077445966144823e-05, -2.4762485623037474e-05, -2.8965803293401466e-04, 9.1705731007558593e-05, 1.0735276393480062e-06, -6.1655607033487706e-07,
        -4.0661080129677908e-06, 2.0480315796477439e-07, -8.9140534094971091e-04, 2.8255712762494871e-04, -2.5740565873720605e-04, 8.1441659842484913e-05,
        -1.2302726069713017e-05, 9.1956134101472709e-06, -3.6566756667073014e-06, -3.1133138100266949e-06, -1.6052853325798863e-03, 5.0887195474260045e-04,
        -1.3890215883913964e-04, 4.3876277596630102e-05, -2.2079428918499198e-05, 2.1473735733559101e-05, -1.9879063548263463e-06, -7.2535244237120365e-06,
        -1.7314656488376526e-03, 5.4871225335417057e-04, 8.4998876851298603e-05, -2.7345163128415205e-05, -2.4067335201346147e-05, 3.6149985959550321e-05,
        8.4018812730091390e-07, -1.2265701405655386e-05, -8.7137255259335414e-04, 2.7692821497714472e-04, 4.3399692520102665e-04, -1.3680952581435494e-04,
        -1.1306043338472012e-05, 5.3517079543489643e-05, 6.5275471649367974e-06, -1.7873859699008542e-05, 3.6938986961213048e-04, -3.0543697973039519e-04,
        -1.0151390511348573e-05, 8.4358553295444398e-06, 7.1394295431845384e-05, -3.1831530137326966e-06, -1.9945455905808060e-06, 1.8377942493472945e-06,
        3.0071300409620594e-04, -2.4864383043358133e-04, -2.8842534029753127e-05, 2.3974549030357833e-05, 5.8115739105822045e-05, -2.6744891288801269e-06,
        -5.6718324786983921e-06, 1.5441170185031866e-06, 1.7452625589438092e-04, -1.4429115934672037e-04, -4.2698208395339102e-05, 3.5515311208936911e-05,
        3.3716937023161008e-05, -1.6572805586334360e-06, -8.4148224716303384e-06, 9.5683137672264763e-07, 1.3158335599317376e-05, -1.0850733523082882e-05,
        -4.8497072164827641e-05, 4.0388867576268002e-05, 2.5203554716036562e-06, -1.3091329757515699e-07, -9.5965322849574142e-06, 7.5582827618046964e-08,
        -1.4987078858603924e-04, 1.2398936849488513e-04, -4.3003875953551950e-05, 3.5946157253459655e-05, -2.9017497633503618e-05, 1.9015978556327356e-06,
        -8.6118556949391419e-06, -1.0978880338456148e-06, -2.7001419587187905e-04, 2.2322219215679677e-04, -2.3064354890750821e-05, 1.9442648634427122e-05,
        -5.2153147173544714e-05, 4.4551314488746577e-06, -4.7455177948804227e-06, -2.5721713412822791e-06, -2.9075952304168459e-04, 2.4101035061782007e-04,
        1.4956729226399596e-05, -1.1406819755368516e-05, -5.6654203115960673e-05, 7.4525839672430929e-06, 2.1468319348890182e-06, -4.3027513596461785e-06,
        -1.4931594530876401e-04, 1.2037750100710676e-04, 7.1646057570357943e-05, -6.1092871465861189e-05, -2.6467890564710737e-05, 1.1279323208724605e-05,
        1.5281243742417297e-05, -6.5121202908339018e-06,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[512] = {
        2.9454012202039959e-03, -8.1220101364014335e-05, 1.5803009873427008e-03, -4.3533378910194729e-05, -2.0166877752327100e-05, -9.3127588681861211e-05,
        -1.1643352297484564e-05, 2.6019232525726185e-06, 7.1908056327863751e-03, -1.9820583297107370e-04, 8.7835370939146529e-04, -2.4220287994320529e-05,
        -5.0761109464237829e-05, -8.7264421478999864e-05, -6.0202356170370377e-06, 2.4377681936242949e-06, 9.1097483696012064e-03, -2.5111007451780168e-04,
        2.4472587255051442e-04, -6.7478525891956605e-06, -6.4108763129443647e-05, -7.5507998903114623e-05, -1.6860358196197614e-06, 2.1094457243410285e-06,
        9.0177576000672853e-03, -2.4857247581879416e-04, -2.7507878151668448e-04, 7.5771392648586923e-06, -6.3502919856844940e-05, -5.7879738294026940e-05,
        2.0358195961416349e-06, 1.6169597830366666e-06, 7.3879826190809909e-03, -2.0364953845098299e-04, -6.3552670738178235e-04, 1.7511406628100266e-05,
        -5.2007975029575144e-05, -3.4374187885035892e-05, 4.6007898942026179e-06, 9.6031552715662194e-07, 4.8513141176030304e-03, -1.3372486950495100e-04,
        -7.9108960093056456e-04, 2.1799972909633512e-05, -3.4186754722054269e-05, -4.9937874519117131e-06, 5.6882964476321752e-06, 1.3951656687378755e-07,
        2.1963303396407945e-03, -6.0549680511663051e-05, -6.9624475060609892e-04, 1.9176179028479887e-05, -1.5325073298283675e-05, 3.0266819135900875e-05,
        5.2015003997508587e-06, -8.4543199460802534e-07, 3.6941480543372203e-04, -1.0153105900652723e-05, -3.0544072884953544e-04, 8.4368216630257555e-06,
        -3.1579051656142495e-06, 7.1386047849479530e-05, 1.8232173969212332e-06, -1.9945269638768098e-06, 2.3977384860881484e-03, -2.3080698796021753e-04,
        1.2864708014675135e-03, -1.2370423731258926e-04, -1.6939801410706055e-05, -7.5805251541327030e-05, -9.7801989038748765e-06, 7.3991327585082440e-06,
        5.8537837476602588e-03, -5.6323841000574222e-04, 7.1503602061339885e-04, -6.8827941909748464e-05, -4.2638919377799938e-05, -7.1035161331602551e-05,
        -5.0571937723627049e-06, 6.9321995212242636e-06, 7.4159243241399952e-03, -7.1357723361086141e-04, 1.9922181550321281e-04, -1.9175593602971073e-05,
        -5.3851197350447070e-05, -6.1464590531628527e-05, -1.4162179333740790e-06, 5.9985198792762209e-06, 7.3410359827878148e-03, -7.0636585601902880e-04,
        -2.2393338554223236e-04, 2.1531100167414857e-05, -5.3342267913384456e-05, -4.7115005978443583e-05, 1.7100484808606391e-06, 4.5980613171196657e-06,
        6.0142908404055067e-03, -5.7870901166612787e-04, -5.1736136667162170e-04, 4.9760990572992248e-05, -4.3686662652689917e-05, -2.7980986532591139e-05,
        3.8646178155903161e-06, 2.7308009946611132e-06, 3.9492776173771678e-03, -3.8000453723856020e-04, -6.4399912438797466e-04, 6.1947829473316715e-05,
        -2.8716855538774154e-05, -4.0649664705468967e-06, 4.7782043513454510e-06, 3.9673847674617993e-07, 1.7879480608208559e-03, -1.7206478508547871e-04,
        -5.6678921739335923e-04, 5.4490311572900552e-05, -1.2873356317686101e-05, 2.4638378877090847e-05, 4.3690441888556233e-06, -2.4041495039335977e-06,
        3.0073461796833155e-04, -2.8847556652659043e-05, -2.4864554358644894e-04, 2.3977541039039716e-05, -2.6529749016121289e-06, 5.8107400092355125e-05,
        1.5316957729353425e-06, -5.6719038931731779e-06, 1.3914632549221905e-03, -3.4183397060441980e-04, 7.4658797035895882e-04, -1.8318659507990241e-04,
        -1.0485616393844903e-05, -4.3976730104647970e-05, -6.0538734470919341e-06, 1.0977072660866134e-05, 3.3971404391524500e-03, -8.3413379085504555e-04,
        4.1495638642054895e-04, -1.0193656011686751e-04, -2.6394484996517297e-05, -4.1214445138846047e-05, -3.1311161231633720e-06, 1.0284799000091009e-05,
        4.3036929380287376e-03, -1.0567854374164908e-03, 1.1561251928664079e-04, -2.8399600782053514e-05, -3.3336037467056137e-05, -3.5660140506894260e-05,
        -8.7659106429607982e-07, 8.8996862887946879e-06, 4.2602284462730728e-03, -1.0461044936674639e-03, -1.2995915103088186e-04, 3.1883814142634087e-05,
        -3.3020948276874124e-05, -2.7335082036188608e-05, 1.0585078930664303e-06, 6.8218830954918780e-06, 3.4902698295759737e-03, -8.5704963373963946e-04,
        -3.0024483422289489e-04, 7.3690512324013055e-05, -2.7044015848146516e-05, -1.6233694748892458e-05, 2.3922756535875034e-06, 4.0515010789063269e-06,
        2.2918752967835201e-03, -5.6277390699299343e-04, -3.7373562454794834e-04, 9.1739015577030548e-05, -1.7777047290544649e-05, -2.3583854090292803e-06,
        2.9580111377158402e-06, 5.8855655851499056e-07, 1.0375854517644043e-03, -2.5482699179924987e-04, -3.2892780489151413e-04, 8.0689634050916723e-05,
        -7.9699444943101593e-06, 1.4296333939449161e-05, 2.7041223016604654e-06, -3.5668329247850783e-06, 1.7454101159724596e-04, -4.2706656492167086e-05,
        -1.4428946137282037e-04, 3.5519779882484225e-05, -1.6431336394700643e-06, 3.3709022826938628e-05, 9.4866364849333451e-07, -8.4145057888050827e-06,
        1.0467737661037531e-04, -3.8859382202929388e-04, 5.6203300802310276e-05, -2.0819123633224822e-04, -8.0441241093018447e-07, -3.2782563589691267e-06,
        -4.6442772341934260e-07, 1.2520202111808520e-05, 2.5567864186943284e-04, -9.4813297672244443e-04, 3.1224884744875965e-05, -1.1588107807708817e-04,
        -2.0279613731276092e-06, -3.0864937978382444e-06, -2.4198859927203242e-07, 1.1728383970288819e-05, 3.2388834732881148e-04, -1.2012284440763674e-03,
        8.6965617261401824e-06, -3.2283330170791101e-05, -2.5633487783052417e-06, -2.6677819565688528e-06, -6.7117463228011717e-08, 1.0148460801436374e-05,
        3.2060674566349419e-04, -1.1890850795523928e-03, -9.7875823604988555e-06, 3.6235571281476438e-05, -2.5389611769987272e-06, -2.0453060262715725e-06,
        8.1197651407210826e-08, 7.7791758915792484e-06, 2.6264541286741977e-04, -9.7419182959743796e-04, -2.2603361393226133e-05, 8.3754690861529983e-05,
        -2.0800336062595361e-06, -1.2141271317931038e-06, 1.8376430509761054e-07, 4.6200503279378238e-06, 1.7244890955450532e-04, -6.3969254246846118e-04,
        -2.8131778414908485e-05, 1.0427071427267564e-04, -1.3672588878457042e-06, -1.7629489258284240e-07, 2.2775637045017820e-07, 6.7127398855154552e-07,
        7.8045369342150914e-05, -2.8966811930850869e-04, -2.4759193610949960e-05, 9.1698100347770243e-05, -6.1467369026941733e-07, 1.0730809895346153e-06,
        2.0674889595862333e-07, -4.0676157254105289e-06, 1.3166553788989042e-05, -4.8507271627262026e-05, -1.0842831493873463e-05, 4.0396567815994165e-05,
        -1.2828705092508016e-07, 2.5105593003948413e-06, 7.4066561197068482e-08, -9.5979355265476868e-06, -1.1954709182715321e-03, -3.4538151677364744e-04,
        -6.4133408421936382e-04, -1.8491025515425375e-04, 1.2104728484199687e-05, 3.7855602169716078e-05, 6.9886682484584207e-06, 1.1228442182535923e-05,
        -2.9184032434070024e-03, -8.4247985755523552e-04, -3.5648875756049873e-04, -1.0299281551868266e-04, 3.0461873620089127e-05, 3.5464112360899039e-05,
        3.6098344373003926e-06, 1.0528818872881475e-05, -3.6972393027735668e-03, -1.0674005483978549e-03, -9.9333719803364664e-05, -2.8694833445077218e-05,
        3.8467144273081508e-05, 3.0691614043789165e-05, 1.0120107291407189e-06, 9.1115954593737303e-06, -3.6599210540769495e-03, -1.0566042079448816e-03,
        1.1162724518073069e-04, 3.2181124966835425e-05, 3.8103745885793435e-05, 2.3524909676956877e-05, -1.2218188858647094e-06, 6.9837950279166614e-06,
        -2.9984939624715938e-03, -8.6565634565169132e-04, 2.5791494506950463e-04, 7.4401149941114939e-05, 3.1205563407591607e-05, 1.3971871591080566e-05,
        -2.7608486248442764e-06, 4.1475901226266124e-06, -1.9689875593994753e-03, -5.6841961601480321e-04, 3.2105203678306428e-04, 9.2630791582974501e-05,
        2.0512589896522284e-05, 2.0285072267632719e-06, -3.4127425102091094e-06, 6.0166910856288977e-07, -8.9146158245425298e-04, -2.5741717340187448e-04,
        2.8256337544435052e-04, 8.1435397979603627e-05, 9.1915705603400686e-06, -1.2297271353333881e-05, -3.1234510510367390e-06, -3.6517607989047280e-06,
        -1.4988522929525297e-04, -4.3024092511512739e-05, 1.2399895910463496e-04, 3.5949439649306087e-05, 1.8907973218337211e-06, -2.9018468406631176e-05,
        -1.0916523435073203e-06, -8.6053571073909431e-06, -2.1527369007017694e-03, -1.8646896929811231e-04, -1.1550064332283025e-03, -9.9633974313152713e-05,
        2.8237475290976773e-05, 6.8069073396300730e-05, 1.6302913960317764e-05, 6.2153135632857116e-06, -5.2554636227202340e-03, -4.5439874512202377e-04,
        -6.4197688667954006e-04, -5.5616195618962825e-05, 7.1069202708901050e-05, 6.3713352320859192e-05, 8.4259953942776365e-06, 5.7808874224038826e-06,
        -6.6580241381562722e-03, -5.7579789951344741e-04, -1.7885832254778354e-04, -1.5476058264713389e-05, 8.9754331351462769e-05, 5.5141264271663295e-05,
        2.3618686573480313e-06, 5.0044166812814061e-06, -6.5908068798421638e-03, -5.6995646609145616e-04, 2.0104511865453936e-04, 1.7341135947997511e-05,
        8.8907253598389082e-05, 4.2266610866794477e-05, -2.8509292261096717e-06, 3.8367311991075141e-06, -5.3996900516909870e-03, -4.6696178785538418e-04,
        4.6448717187077420e-04, 4.0111112600911183e-05, 7.2811532382451130e-05, 2.5102067002363931e-05, -6.4419397507130650e-06, 2.2784311942111144e-06,
        -3.5457360568204046e-03, -3.0660800028773087e-04, 5.7819348891231510e-04, 4.9956307586287455e-05, 4.7862477545034580e-05, 3.6456019898696849e-06,
        -7.9624037756962899e-06, 3.3196098822177080e-07, -1.6053531486772204e-03, -1.3891805426815166e-04, 5.0886290890303216e-04, 4.3846488387641553e-05,
        2.1453637431742436e-05, -2.2089749982095302e-05, -7.2847471726987051e-06, -2.0019293734442781e-06, -2.6998436645517406e-04, -2.3052517026955414e-05,
        2.2325295411531143e-04, 1.9471183673631829e-05, 4.4180426042971720e-06, -5.2190315405154041e-05, -2.5507580867322769e-06, -4.7729149948270460e-06,
        -2.3220965407560255e-03, 1.1372626845228194e-04, -1.2453897805284240e-03, 6.1906319856063524e-05, 4.7631155065409337e-05, 7.3798367185252568e-05,
        2.7499860198971598e-05, -2.9075042519679062e-06, -5.6686269193926926e-03, 2.7880432247664029e-04, -6.9250187111350413e-04, 3.3876803090129302e-05,
        1.1982659122391187e-04, 6.9327831552477994e-05, 1.4182194301402512e-05, -2.5393663266690694e-06, -7.1813597590535769e-03, 3.5309871945119663e-04,
        -1.9298642871185722e-04, 9.4161352330371667e-06, 1.5129509290712090e-04, 5.9946236111054975e-05, 3.9861536163922774e-06, -2.2302648963597872e-06,
        -7.1088929462936505e-03, 3.4954800622020494e-04, 2.1676647742164119e-04, -1.0788651660103704e-05, 1.4987472726567217e-04, 4.5958038747264697e-05,
        -4.8062021018302064e-06, -1.7054843186236719e-06, -5.8241256952107617e-03, 2.8637215473879330e-04, 5.0089152032476694e-04, -2.4796203091292665e-05,
        1.2274079426761561e-04, 2.7288151718768185e-05, -1.0859581418771334e-05, -1.0162945947238439e-06, -3.8244790783418998e-03, 1.8805805739426196e-04,
        6.2349876160806275e-04, -3.0841720811719915e-05, 8.0691755516135795e-05, 3.9679657721800867e-06, -1.3417442423560333e-05, -1.4584417172785917e-07,
        -1.7315041659079929e-03, 8.5051747437732335e-05, 5.4867496361982941e-04, -2.7337450004950627e-05, 3.6145425568310897e-05, -2.4077006774690238e-05,
        -1.2301393163226440e-05, 8.5458612929063640e-07, -2.9107657039730310e-04, 1.4749655139206393e-05, 2.4104635737920416e-04, -1.1435232611070181e-05,
        7.4193938028526494e-06, -5.6540305945549369e-05, -4.2835890093106258e-06, 2.2614467853574804e-06, -1.1659391664103568e-03, 5.8272261883050261e-04,
        -6.2876933898087868e-04, 3.1080881916009130e-04, 7.0048847106674077e-05, 3.4381211048814196e-05, 4.0442720733512189e-05, -1.9850001454094279e-05,
        -2.8537095930345883e-03, 1.4196417456443398e-03, -3.4813949878996714e-04, 1.7381510224882214e-04, 1.7656795992097995e-04, 3.2464760028068929e-05,
        2.1056117723667480e-05, -1.8743537941105015e-05, -3.6146775624239977e-03, 1.7988700727608567e-03, -9.6915771094449986e-05, 4.8429264571519661e-05,
        2.2313418250579205e-04, 2.8041651998169706e-05, 5.8289034208178767e-06, -1.6189855329597622e-05, -3.5782740283023416e-03, 1.7806429045283578e-03,
        1.0949074126933192e-04, -5.4078376098771870e-05, 2.2097273092492643e-04, 2.1502026627993337e-05, -7.0768180728710684e-06, -1.2414200861742585e-05,
        -2.9315513747036792e-03, 1.4588653495762290e-03, 2.5261393415669178e-04, -1.2518699583789750e-04, 1.8099739002053002e-04, 1.2763938922627853e-05,
        -1.6002955759229244e-05, -7.3692635728953023e-06, -1.9250880771233487e-03, 9.5791667978322391e-04, 3.1438799948895304e-04, -1.5590680667852328e-04,
        1.1892649105348721e-04, 1.8576781283812511e-06, -1.9833694468234978e-05, -1.0725309674870637e-06, -8.7109635416753701e-04, 4.3398455254816179e-04,
        2.7714996862665574e-04, -1.3677957882656216e-04, 5.3403143668107477e-05, -1.1298410089657595e-05, -1.7996227782919417e-05, 6.5231401066786274e-06,
        -1.4932655918166641e-04, 7.1862542461001002e-05, 1.2018196532368101e-04, -6.1101428538653195e-05, 1.1116381401754253e-05, -2.6311682607348453e-05,
        -6.4180457947173718e-06, 1.5191057036184607e-05,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-7) );
      TEST_MSG("This BC combination is not available");
    }
  } else {
    TEST_CHECK( gkyl_compare(1., 2., 1e-7) );
    TEST_MSG("This poly_order is not available");
  }

//  gkyl_grid_sub_array_write(&grid, &localRange, rho,    "ctest_fem_poisson_vareps_2x_rho_1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, eps,    "ctest_fem_poisson_vareps_2x_eps_1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, phisol, "ctest_fem_poisson_vareps_2x_phisol_256x256_p1.gkyl");
//  gkyl_grid_sub_array_write(&grid, &localRange, phi,       "ctest_fem_poisson_vareps_2x_phi_256x256_p1.gkyl");

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  gkyl_array_release(phisol);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }
}

void test_poisson_vareps_2x_p1_periodicx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(1, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p1_dirichletx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p1_dirichletx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(1, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p1_periodicx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(1, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p1_dirichletx_neumanny_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p1_neumannx_dirichletx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p2_periodicx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(2, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p2_dirichletx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p2_dirichletx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(2, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p2_periodicx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(2, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p2_dirichletx_neumanny_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, cells, bc_tv, false);
}

void test_poisson_vareps_2x_p2_neumannx_dirichletx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, cells, bc_tv, false);
}

#ifdef GKYL_HAVE_CUDA
void test_poisson_vareps_2x_p1_periodicx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(1, cells, bc_tv, true);
}

void test_poisson_vareps_2x_p1_dirichletx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(1, cells, bc_tv, true);
}

void test_poisson_vareps_2x_p1_periodicx_dirichlety_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(1, cells, bc_tv, true);
}

void test_poisson_vareps_2x_p2_periodicx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(2, cells, bc_tv, true);
}

void test_poisson_vareps_2x_p2_dirichletx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(2, cells, bc_tv, true);
}

void test_poisson_vareps_2x_p2_periodicx_dirichlety_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(2, cells, bc_tv, true);
}
#endif


TEST_LIST = {
  // 2x tests
  { "test_poisson_vareps_2x_p1_periodicx_periodicy_ho", test_poisson_vareps_2x_p1_periodicx_periodicy_ho },
  { "test_poisson_vareps_2x_p1_dirichletx_dirichlety_ho", test_poisson_vareps_2x_p1_dirichletx_dirichlety_ho },
  { "test_poisson_vareps_2x_p1_dirichletx_periodicy_ho", test_poisson_vareps_2x_p1_dirichletx_periodicy_ho },
  { "test_poisson_vareps_2x_p1_periodicx_dirichlety_ho", test_poisson_vareps_2x_p1_periodicx_dirichlety_ho },
  { "test_poisson_vareps_2x_p1_dirichletx_neumanny_dirichlety_ho", test_poisson_vareps_2x_p1_dirichletx_neumanny_dirichlety_ho },
  { "test_poisson_vareps_2x_p1_neumannx_dirichletx_dirichlety_ho", test_poisson_vareps_2x_p1_neumannx_dirichletx_dirichlety_ho },
  { "test_poisson_vareps_2x_p2_periodicx_periodicy_ho", test_poisson_vareps_2x_p2_periodicx_periodicy_ho },
  { "test_poisson_vareps_2x_p2_dirichletx_dirichlety_ho", test_poisson_vareps_2x_p2_dirichletx_dirichlety_ho },
  { "test_poisson_vareps_2x_p2_dirichletx_periodicy_ho", test_poisson_vareps_2x_p2_dirichletx_periodicy_ho },
  { "test_poisson_vareps_2x_p2_periodicx_dirichlety_ho", test_poisson_vareps_2x_p2_periodicx_dirichlety_ho },
  { "test_poisson_vareps_2x_p2_dirichletx_neumanny_dirichlety_ho", test_poisson_vareps_2x_p2_dirichletx_neumanny_dirichlety_ho },
  { "test_poisson_vareps_2x_p2_neumannx_dirichletx_dirichlety_ho", test_poisson_vareps_2x_p2_neumannx_dirichletx_dirichlety_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_poisson_vareps_2x_p1_periodicx_periodicy_dev", test_poisson_vareps_2x_p1_periodicx_periodicy_dev },
  { "test_poisson_vareps_2x_p1_dirichletx_periodicy_dev", test_poisson_vareps_2x_p1_dirichletx_periodicy_dev },
  { "test_poisson_vareps_2x_p1_periodicx_dirichlety_dev", test_poisson_vareps_2x_p1_periodicx_dirichlety_dev },
  { "test_poisson_vareps_2x_p2_periodicx_periodicy_dev", test_poisson_vareps_2x_p2_periodicx_periodicy_dev },
  { "test_poisson_vareps_2x_p2_dirichletx_periodicy_dev", test_poisson_vareps_2x_p2_dirichletx_periodicy_dev },
  { "test_poisson_vareps_2x_p2_periodicx_dirichlety_dev", test_poisson_vareps_2x_p2_periodicx_dirichlety_dev },
#endif
  { NULL, NULL },
};