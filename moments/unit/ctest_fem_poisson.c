// Test the FEM Poisson solver.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_array_reduce.h>

void evalFunc1x_periodicx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sin((2.*M_PI/(2.*M_PI))*x);
}
void evalFunc1x_dirichletx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 2.0;
  double c0 = a/12. - 1./2.;
  double c1 = 0.;
  fout[0] = -(1.-a*pow(x,2));
}
void evalFunc1x_neumannx_dirichletx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 5.0;
  double c0 = 0.;
  double c1 = a/12. - 1./2.;
  fout[0] = -(1.-a*pow(x,2));
}
void evalFunc1x_dirichletx_neumannx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 5.0;
  double c0 = 0.;
  double c1 = a/12. - 1./2.;
  fout[0] = -(1.-a*pow(x-1.,2));
}

void evalFunc2x_periodicx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = sin((2.*M_PI/(2.*M_PI))*x)*cos((2.*2.*M_PI/(2.*M_PI))*y);
}
void evalFunc2x_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 2.*M_PI;
  double kx = 2.*M_PI/Lx, ky = 2.*M_PI/Ly;
  double sig = 0.3*sqrt(Lx*Lx+Ly*Ly);
  fout[0] = exp(-(pow(kx*x,2)+pow(ky*y,2))/(2.*(sig*sig)));
}
void evalFunc2x_dirichletx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = (y/(2.*M_PI)+0.5)*sin(5.0*M_PI*(y/(2.*M_PI)+0.5))+exp(-(pow(x/(2.*M_PI)+0.5-0.5,2)+pow(y/(2.*M_PI)+0.5-0.5,2))/0.02);
}
void evalFunc2x_periodicx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = (x/(2.*M_PI)+0.5)*sin(5.0*M_PI*(x/(2.*M_PI)+0.5))+exp(-(pow(x/(2.*M_PI)+0.5-0.5,2)+pow(y/(2.*M_PI)+0.5-0.5,2))/0.02);
}
void evalFunc2x_dirichletx_neumanny_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double a = 2.0, b = 5.0;
  double c1 = 0., c0 = a/12. - 1./2.;
  double d1 = b/12. - 1./2., d0 = 0.;
  fout[0] = -( (1.-a*pow(x,2))*(-b*pow(y,4)/12.+pow(y,2)/2.+d0*y+d1)
              +(1.-b*pow(y,2))*(-a*pow(x,4)/12.+pow(x,2)/2.+c0*x+c1) );
}
void evalFunc2x_dirichletx_dirichlety_neumanny(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double a = 2.0, b = 5.0;
  double c1 = 0., c0 = a/12. - 1./2.;
  double d1 = b/12. - 1./2., d0 = 0.;
  fout[0] = -( (1.-a*pow(x,2))*(-b*pow(y-1.,4)/12.+pow(y-1.,2)/2.+d0*y+d1)
              +(1.-b*pow(y-1.,2))*(-a*pow(x,4)/12.+pow(x,2)/2.+c0*x+c1) );
}
void evalFunc2x_neumannx_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double a = 2.0, b = 5.0;
  double c1 = 0., c0 = a/12. - 1./2.;
  double d1 = b/12. - 1./2., d0 = 0.;
  fout[0] = -( (1.-a*pow(y,2))*(-b*pow(x,4)/12.+pow(x,2)/2.+d0*x+d1)
              +(1.-b*pow(x,2))*(-a*pow(y,4)/12.+pow(y,2)/2.+c0*y+c1) );
}
void evalFunc2x_dirichletx_neumannx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double a = 2.0, b = 5.0;
  double c1 = 0., c0 = a/12. - 1./2.;
  double d1 = b/12. - 1./2., d0 = 0.;
  fout[0] = -( (1.-a*pow(y,2))*(-b*pow(x-1.,4)/12.+pow(x-1.,2)/2.+d0*x+d1)
              +(1.-b*pow(x-1.,2))*(-a*pow(y,4)/12.+pow(y,2)/2.+c0*y+c1) );
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
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
  gkyl_array_copy_to_buffer(buff->data, fld,   &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld,   &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
}

void
test_1x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {0.0}, upper[] = {1.0};
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    lower[0] = -M_PI;
    upper[0] =  M_PI;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob;
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_periodicx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_dirichletx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_neumannx_dirichletx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_dirichletx_neumannx, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
    phi_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  gkyl_array_clear(epsilon, 0.);
  gkyl_array_shiftc(epsilon, epsilon_0*pow(sqrt(2.),dim), 0);

  // Project the right-side source on the basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  struct gkyl_array *perbuff = mkarr(use_gpu, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);

//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, rho_ho, "ctest_fem_poisson_1x_rho_1.gkyl");

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bcs, NULL, epsilon, NULL, true, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho, NULL);
  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    cudaDeviceSynchronize();
  }
#endif

  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);

  gkyl_array_copy(phi_ho, phi);
//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, phi_ho, "ctest_fem_poisson_1x_phi_1.gkyl");

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    struct gkyl_array *sol_cellavg = mkarr(false, 1, localRange_ext.volume);
    double* sol_avg = (double*) gkyl_malloc(sizeof(double)); 
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim)/localRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_dg_calc_average_range(&basis, 0, sol_cellavg, 0, phi_ho, localRange);
    gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &localRange);
    gkyl_array_shiftc(phi_ho, mavgfac*sol_avg[0], 0);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
	-2.7060255565796842e-01, -1.5623245835252897e-01,
	-7.7061088089325591e-01, -1.3244748281911814e-01,
	-1.1533006851175758e+00, -8.8498578665918093e-02,
	-1.3604109147294894e+00, -3.1076568152445378e-02,
	-1.3604109147294894e+00, 3.1076568152445333e-02,
	-1.1533006851175760e+00, 8.8498578665918051e-02,
	-7.7061088089325636e-01, 1.3244748281911811e-01,
	-2.7060255565796909e-01, 1.5623245835252891e-01,
	2.7060255565796831e-01, 1.5623245835252891e-01,
	7.7061088089325547e-01, 1.3244748281911814e-01,
	1.1533006851175758e+00, 8.8498578665918204e-02,
	1.3604109147294898e+00, 3.1076568152445395e-02,
	1.3604109147294898e+00, -3.1076568152445395e-02,
	1.1533006851175758e+00, -8.8498578665918023e-02,
	7.7061088089325658e-01, -1.3244748281911806e-01,
	2.7060255565796898e-01, -1.5623245835252900e-01,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
        -0.0133521216082515, -0.0077088510047765,
        -0.0373194046782545, -0.0061286663274337,
        -0.0558775050145713, -0.0045858578973869,
        -0.0691990561087026, -0.0031053432128332,
        -0.0775430081978993, -0.0017120397719699,
        -0.0812546282651626, -0.0004308650729942,
        -0.0807655000392442,  0.0007132633858965,
        -0.0765935239946458,  0.0016954281065051,
        -0.0693429173516197,  0.0024907115906342,
        -0.0597042140761683,  0.0030741963400866,
        -0.0484542648800444,  0.0034209648566651,
        -0.0364562372207511,  0.0035060996421724,
        -0.024659615301542 ,  0.0033046831984112,
        -0.0141002000714206,  0.0027917980271844,
        -0.0059001092251411,  0.0019425266302945,
        -0.0012677772032077,  0.0007319515095444,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
        -0.1164745579295943,  0.0007947643695833,
        -0.1110222168230329,  0.0023531462360035,
        -0.1003333264742861,  0.0038180874841839,
        -0.0848394706121053,  0.0051272943686315,
        -0.0651880248296178,  0.006218473143853 ,
        -0.0422421565843268,  0.0070293300643554,
        -0.0170808251981111,  0.0074975713846456,
         0.0090012181427745,  0.0075609033592303,
         0.0344934303876997,  0.0071570322426164,
         0.0576694766216582,  0.0062236642893108,
         0.0765872300652681,  0.0046985057538203,
         0.0890887720747714,  0.0025192628906519,
         0.0928003921420348, -0.0003763580456877,
         0.0851325878945489, -0.0040506508006916,
         0.0632800650954285, -0.008565909119853 ,
         0.0242217376434128, -0.0139844267486649,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) {
      // Solution with N=16 (checked visually against g2):
      const double sol[32] = {
         0.0242217376434127,  0.0139844267486649,
         0.0632800650954284,  0.008565909119853 ,
         0.0851325878945488,  0.0040506508006916,
         0.0928003921420348,  0.0003763580456877,
         0.0890887720747714, -0.0025192628906519,
         0.076587230065268 , -0.0046985057538203,
         0.0576694766216582, -0.0062236642893108,
         0.0344934303876997, -0.0071570322426164,
         0.0090012181427746, -0.0075609033592302,
        -0.0170808251981111, -0.0074975713846456,
        -0.0422421565843267, -0.0070293300643554,
        -0.0651880248296177, -0.006218473143853 ,
        -0.0848394706121052, -0.0051272943686315,
        -0.100333326474286 , -0.0038180874841839,
        -0.1110222168230328, -0.0023531462360035,
        -0.1164745579295941, -0.0007947643695833,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
  } if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
	-5.2739372946099083e-01, -2.8867496150229432e-01, 1.2250982323291082e-02, 
	-1.2732410943752506e+00, -1.1957308417178132e-01, 2.9576487677282437e-02, 
	-1.2732410943752503e+00, 1.1957308417178146e-01, 2.9576487677282326e-02, 
	-5.2739372946099017e-01, 2.8867496150229438e-01, 1.2250982323291082e-02, 
	5.2739372946099117e-01, 2.8867496150229432e-01, -1.2250982323291082e-02, 
	1.2732410943752508e+00, 1.1957308417178118e-01, -2.9576487677282381e-02, 
	1.2732410943752503e+00, -1.1957308417178146e-01, -2.9576487677282326e-02, 
	5.2739372946099039e-01, -2.8867496150229427e-01, -1.2250982323291082e-02,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
        -2.5791443630192343e-02, -1.3837517332209948e-02,  8.1578940289174257e-04,
        -6.2965188799992161e-02, -7.6912011102198800e-03,  7.6432003993326847e-04,
        -7.9768181972718816e-02, -2.1429048449640056e-03,  6.6138131401632689e-04,
        -7.8962559012382361e-02,  2.4086914924015437e-03,  5.0697322514091727e-04,
        -6.4691523714997728e-02,  5.5649079307207141e-03,  3.0109577330704072e-04,
        -4.2479347808584755e-02,  6.9270644988373735e-03,  4.3748958514695155e-05,
        -1.9231370953168247e-02,  6.0964812255954216e-03, -2.6506721923612130e-04,
        -3.2340007407779619e-03,  2.6744781398387418e-03, -6.2535275994540806e-04,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
        -1.1419705462765599e-01,  3.1479106055867533e-03,  8.0420879622607536e-04,
        -9.2963135173079955e-02,  8.9453818528151628e-03,  6.7553538882990415e-04,
        -5.3947966093940251e-02,  1.3247803208208148e-02,  4.1818857403756060e-04,
        -4.0568870502616112e-03,  1.5058474743875470e-02,  3.2168351849042897e-05,
         4.6352092467918868e-02,  1.3380696531926833e-02, -4.8252527773564799e-04,
         8.3468293140551653e-02,  7.2177686444719596e-03, -1.1258923147165170e-03,
         9.0028365817574674e-02, -4.4270088463794948e-03, -1.8979327590935657e-03,
         4.5316291518913317e-02, -2.2550335868517860e-02, -2.7986466108667861e-03,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8 (checked visually against g2):
      const double sol[24] =  {
         4.5316291518912914e-02,  2.2550335868517835e-02, -2.7986466108667865e-03,
         9.0028365817574202e-02,  4.4270088463794748e-03, -1.8979327590935605e-03,
         8.3468293140551125e-02, -7.2177686444719657e-03, -1.1258923147165190e-03,
         4.6352092467918285e-02, -1.3380696531926861e-02, -4.8252527773564712e-04,
        -4.0568870502622808e-03, -1.5058474743875493e-02,  3.2168351849043662e-05,
        -5.3947966093940987e-02, -1.3247803208208161e-02,  4.1818857403755978e-04,
        -9.2963135173080691e-02, -8.9453818528151541e-03,  6.7553538882990382e-04,
        -1.1419705462765670e-01, -3.1479106055867477e-03,  8.0420879622607980e-04,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(perbuff);
}

void
test_1x_bias(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob;
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_periodicx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_dirichletx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_neumannx_dirichletx, NULL);
  } else if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc1x_dirichletx_neumannx, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
    phi_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  gkyl_array_clear(epsilon, 0.);
  gkyl_array_shiftc(epsilon, epsilon_0*pow(sqrt(2.),dim), 0);

  // Project the right-side source on the basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, rho_ho, "ctest_fem_poisson_1x_rho_1.gkyl");

  // Specify the bias:
  struct gkyl_poisson_bias_plane bias = {
    .dir = 0,
    .loc = 0.5, // Location of the plane in the 'dir' dimension.
    .val = 0., // Biasing value.
  };
  struct gkyl_poisson_bias_plane_list bpl = {
    .num_bias_plane = 1, // Number of bias planes.
    .bp = &bias,
  };

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bcs, &bpl, epsilon, NULL, true, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho, NULL);
  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    cudaDeviceSynchronize();
  }
#endif

  gkyl_array_copy(phi_ho, phi);
//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, phi_ho, "ctest_fem_poisson_1x_phi_1.gkyl");

  if (poly_order == 1) {
    if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=16:
      const double sol[32] = {
        -0.0087485618349013, -0.0050509845304024,
        -0.0235087253582045, -0.0034707998530596,
        -0.0328597061478215, -0.0019279914230128,
        -0.036974137695253 , -0.0004474767384591,
        -0.0361109702377499,  0.0009458267024042,
        -0.0306154707583134,  0.0022270014013799,
        -0.0209192229856952,  0.0033711298602706,
        -0.007540127394397 ,  0.0043532945808792,
        -0.0002895207513708, -0.0001671548837399,
         0.0001420629773807,  0.0004163298657125,
         0.0021848926268047,  0.000763098382291 ,
         0.0049758007393981,  0.0008482331677983,
         0.0075653031119075,  0.0006468167240371,
         0.008917598795329 ,  0.0001339315528103,
         0.0079105700949087, -0.0007153398440796,
         0.0033357825701422, -0.0019259149648297,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
    else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=16:
      const double sol[32] = {
        -0.1385716448416738,  0.0007947643695833,
        -0.1331193037351125,  0.0023531462360035,
        -0.1224304133863657,  0.0038180874841839,
        -0.1069365575241848,  0.0051272943686315,
        -0.0872851117416974,  0.006218473143853 ,
        -0.0643392434964064,  0.0070293300643554,
        -0.0391779121101907,  0.0074975713846456,
        -0.013095868769305 ,  0.0075609033592303,
         0.0137774114076251,  0.0079543921849286,
         0.0397155935055936,  0.007021024231623 ,
         0.0613954828132133,  0.0054958656961326,
         0.0766591606867267,  0.0033166228329641,
         0.083132916618    ,  0.0004210018966245,
         0.078227248234524 , -0.0032532908583794,
         0.0591368612994136, -0.0077685491775407,
         0.0228406697114078, -0.0131870668063526,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
  }
  else if (poly_order == 2) {
    if (bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8:
      const double sol[24] =  {
        -1.6584324083492974e-02, -8.5217843834619264e-03,  8.1578940289174376e-04,
        -3.5343830159893758e-02, -2.3754681614718402e-03,  7.6432003993327194e-04,
        -3.3732584239221292e-02,  3.1728281037840698e-03,  6.6138131401633307e-04,
        -1.4512722185485518e-02,  7.7244244411496937e-03,  5.0697322514092562e-04,
        -2.4168688810088715e-04,  2.4917498197256333e-04,  3.0109577330704798e-04,
         3.5562499249127679e-03,  1.6113315500892961e-03,  4.3748958514699756e-05,
         8.3899876869301757e-03,  7.8074827684739695e-04, -2.6506721923611772e-04,
         5.9731188059214974e-03, -2.6412548089092564e-03, -6.2535275994540730e-04,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
    else if (bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) {
      // Solution with N=8:
      const double sol[24] =  {
        -1.3629414153973707e-01,  3.1479106055867416e-03,  8.0420879622607210e-04,
        -1.1506022208516119e-01,  8.9453818528150968e-03,  6.7553538882990502e-04,
        -7.6045053006021726e-02,  1.3247803208208071e-02,  4.1818857403755778e-04,
        -2.6153973962343396e-02,  1.5058474743875368e-02,  3.2168351849040308e-05,
         2.7017141419847148e-02,  1.4975416416551469e-02, -4.8252527773565146e-04,
         6.9657613820500383e-02,  8.8124885290965828e-03, -1.1258923147165133e-03,
         8.1741958225543906e-02, -2.8322889617548573e-03, -1.8979327590935668e-03,
         4.2554155654903059e-02, -2.0955615983893214e-02, -2.7986466108667878e-03,
      };

      for (int k=0; k<cells[0]; k++) {
        long linidx;
        const double *phi_p;
        int idx0[] = {k+1};
        linidx = gkyl_range_idx(&localRange, idx0);
        phi_p  = gkyl_array_cfetch(phi_ho, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[k*basis.num_basis+m], phi_p[m], 1e-12) );
          TEST_MSG("Expected: %.13e in cell (%d)", sol[k*basis.num_basis+m], k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
	}
      }
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
}

void
test_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  if ( ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_NEUMANN))
      || ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
     ) {
    lower[0] = 0.; lower[1] = 0.;
    upper[0] = 1.; upper[1] = 1.;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob;
  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_periodicy, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_PERIODIC && bcs.up_type[0]==GKYL_POISSON_PERIODIC) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_periodicx_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_neumanny_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_NEUMANN)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_dirichlety_neumanny, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_neumannx_dirichletx_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) &&
             (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_neumannx_dirichlety, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
    phi_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  gkyl_array_clear(epsilon, 0.);
  gkyl_array_shiftc(epsilon, epsilon_0*pow(sqrt(2.),dim), 0);

  // Only used in the fully periodic case.
  struct gkyl_array *rho_cellavg = mkarr(use_gpu, 1, localRange_ext.volume);

  // Project the right-side source on the basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  struct gkyl_array *perbuff = mkarr(use_gpu, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);

//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, rho_ho, "ctest_fem_poisson_2x_rho_1.gkyl");

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bcs, NULL, epsilon, NULL, true, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho, NULL);
  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    cudaDeviceSynchronize();
  }
#endif

  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);

//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, phi_ho, "ctest_fem_poisson_2x_phi_1.gkyl");

  gkyl_array_copy(phi_ho, phi);

  if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
    struct gkyl_array *sol_cellavg = mkarr(false, 1, localRange_ext.volume);
    double* sol_avg = (double*) gkyl_malloc(sizeof(double)); 
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    double mavgfac = -pow(sqrt(2.),dim)/localRange.volume;
    // Subtract the volume averaged sol from the sol.
    gkyl_dg_calc_average_range(&basis, 0, sol_cellavg, 0, phi_ho, localRange);
    gkyl_array_reduce_range(sol_avg, sol_cellavg, GKYL_SUM, &localRange);
    gkyl_array_shiftc(phi_ho, mavgfac*sol_avg[0], 0);

    gkyl_free(sol_avg);
    gkyl_array_release(sol_cellavg);
  }

  if (poly_order == 1) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
	-7.6861958028898086e-02, -4.4376272158426031e-02, 4.4376272158426031e-02, 2.5620652676299389e-02,
	7.6861958028898114e-02, 4.4376272158426100e-02, 4.4376272158426051e-02, 2.5620652676299389e-02,
	7.6861958028898059e-02, 4.4376272158426072e-02, -4.4376272158426079e-02, -2.5620652676299396e-02,
	-7.6861958028898225e-02, -4.4376272158426079e-02, -4.4376272158426065e-02, -2.5620652676299392e-02,
	-7.6861958028898170e-02, -4.4376272158426079e-02, 4.4376272158426100e-02, 2.5620652676299396e-02,
	7.6861958028898197e-02, 4.4376272158426065e-02, 4.4376272158426086e-02, 2.5620652676299392e-02,
	7.6861958028898253e-02, 4.4376272158426079e-02, -4.4376272158426051e-02, -2.5620652676299382e-02,
	-7.6861958028898031e-02, -4.4376272158426031e-02, -4.4376272158426072e-02, -2.5620652676299389e-02,
	-1.8556118150391737e-01, -1.8381253775579646e-02, 1.0713379809243175e-01, 1.0612421815373736e-02,
	1.8556118150391748e-01, 1.8381253775579632e-02, 1.0713379809243175e-01, 1.0612421815373724e-02,
	1.8556118150391740e-01, 1.8381253775579632e-02, -1.0713379809243180e-01, -1.0612421815373727e-02,
	-1.8556118150391759e-01, -1.8381253775579646e-02, -1.0713379809243177e-01, -1.0612421815373736e-02,
	-1.8556118150391751e-01, -1.8381253775579632e-02, 1.0713379809243183e-01, 1.0612421815373741e-02,
	1.8556118150391762e-01, 1.8381253775579681e-02, 1.0713379809243183e-01, 1.0612421815373745e-02,
	1.8556118150391770e-01, 1.8381253775579684e-02, -1.0713379809243177e-01, -1.0612421815373741e-02,
	-1.8556118150391729e-01, -1.8381253775579632e-02, -1.0713379809243180e-01, -1.0612421815373745e-02,
	-1.8556118150391737e-01, 1.8381253775579646e-02, 1.0713379809243175e-01, -1.0612421815373736e-02,
	1.8556118150391746e-01, -1.8381253775579646e-02, 1.0713379809243173e-01, -1.0612421815373734e-02,
	1.8556118150391737e-01, -1.8381253775579646e-02, -1.0713379809243177e-01, 1.0612421815373736e-02,
	-1.8556118150391757e-01, 1.8381253775579663e-02, -1.0713379809243176e-01, 1.0612421815373745e-02,
	-1.8556118150391748e-01, 1.8381253775579656e-02, 1.0713379809243181e-01, -1.0612421815373750e-02,
	1.8556118150391762e-01, -1.8381253775579656e-02, 1.0713379809243183e-01, -1.0612421815373741e-02,
	1.8556118150391770e-01, -1.8381253775579670e-02, -1.0713379809243177e-01, 1.0612421815373731e-02,
	-1.8556118150391729e-01, 1.8381253775579632e-02, -1.0713379809243180e-01, 1.0612421815373745e-02,
	-7.6861958028898114e-02, 4.4376272158426017e-02, 4.4376272158426044e-02, -2.5620652676299378e-02,
	7.6861958028898086e-02, -4.4376272158426079e-02, 4.4376272158426031e-02, -2.5620652676299378e-02,
	7.6861958028898045e-02, -4.4376272158426058e-02, -4.4376272158426065e-02, 2.5620652676299389e-02,
	-7.6861958028898197e-02, 4.4376272158426065e-02, -4.4376272158426051e-02, 2.5620652676299378e-02,
	-7.6861958028898142e-02, 4.4376272158426051e-02, 4.4376272158426079e-02, -2.5620652676299389e-02,
	7.6861958028898197e-02, -4.4376272158426107e-02, 4.4376272158426100e-02, -2.5620652676299396e-02,
	7.6861958028898253e-02, -4.4376272158426114e-02, -4.4376272158426065e-02, 2.5620652676299396e-02,
	-7.6861958028898059e-02, 4.4376272158426017e-02, -4.4376272158426072e-02, 2.5620652676299378e-02,
	7.6861958028898031e-02, 4.4376272158426031e-02, -4.4376272158426051e-02, -2.5620652676299375e-02,
	-7.6861958028898170e-02, -4.4376272158426051e-02, -4.4376272158426051e-02, -2.5620652676299375e-02,
	-7.6861958028898170e-02, -4.4376272158426044e-02, 4.4376272158426051e-02, 2.5620652676299378e-02,
	7.6861958028898072e-02, 4.4376272158426065e-02, 4.4376272158426065e-02, 2.5620652676299389e-02,
	7.6861958028898114e-02, 4.4376272158426065e-02, -4.4376272158426051e-02, -2.5620652676299382e-02,
	-7.6861958028898114e-02, -4.4376272158426058e-02, -4.4376272158426051e-02, -2.5620652676299392e-02,
	-7.6861958028898114e-02, -4.4376272158426079e-02, 4.4376272158426051e-02, 2.5620652676299378e-02,
	7.6861958028898072e-02, 4.4376272158426017e-02, 4.4376272158426031e-02, 2.5620652676299378e-02,
	1.8556118150391737e-01, 1.8381253775579677e-02, -1.0713379809243177e-01, -1.0612421815373757e-02,
	-1.8556118150391754e-01, -1.8381253775579670e-02, -1.0713379809243175e-01, -1.0612421815373741e-02,
	-1.8556118150391748e-01, -1.8381253775579663e-02, 1.0713379809243177e-01, 1.0612421815373750e-02,
	1.8556118150391743e-01, 1.8381253775579656e-02, 1.0713379809243177e-01, 1.0612421815373731e-02,
	1.8556118150391748e-01, 1.8381253775579656e-02, -1.0713379809243175e-01, -1.0612421815373731e-02,
	-1.8556118150391740e-01, -1.8381253775579625e-02, -1.0713379809243176e-01, -1.0612421815373736e-02,
	-1.8556118150391743e-01, -1.8381253775579632e-02, 1.0713379809243175e-01, 1.0612421815373731e-02,
	1.8556118150391743e-01, 1.8381253775579687e-02, 1.0713379809243175e-01, 1.0612421815373750e-02,
	1.8556118150391743e-01, -1.8381253775579642e-02, -1.0713379809243181e-01, 1.0612421815373738e-02,
	-1.8556118150391765e-01, 1.8381253775579615e-02, -1.0713379809243181e-01, 1.0612421815373708e-02,
	-1.8556118150391762e-01, 1.8381253775579601e-02, 1.0713379809243183e-01, -1.0612421815373717e-02,
	1.8556118150391746e-01, -1.8381253775579642e-02, 1.0713379809243180e-01, -1.0612421815373722e-02,
	1.8556118150391751e-01, -1.8381253775579632e-02, -1.0713379809243176e-01, 1.0612421815373727e-02,
	-1.8556118150391737e-01, 1.8381253775579656e-02, -1.0713379809243175e-01, 1.0612421815373750e-02,
	-1.8556118150391740e-01, 1.8381253775579663e-02, 1.0713379809243173e-01, -1.0612421815373745e-02,
	1.8556118150391748e-01, -1.8381253775579642e-02, 1.0713379809243177e-01, -1.0612421815373736e-02,
	7.6861958028898114e-02, -4.4376272158426051e-02, -4.4376272158426100e-02, 2.5620652676299389e-02,
	-7.6861958028898281e-02, 4.4376272158426114e-02, -4.4376272158426114e-02, 2.5620652676299406e-02,
	-7.6861958028898308e-02, 4.4376272158426107e-02, 4.4376272158426107e-02, -2.5620652676299410e-02,
	7.6861958028898086e-02, -4.4376272158426079e-02, 4.4376272158426086e-02, -2.5620652676299396e-02,
	7.6861958028898156e-02, -4.4376272158426079e-02, -4.4376272158426051e-02, 2.5620652676299396e-02,
	-7.6861958028898059e-02, 4.4376272158426044e-02, -4.4376272158426031e-02, 2.5620652676299378e-02,
	-7.6861958028898059e-02, 4.4376272158426058e-02, 4.4376272158426044e-02, -2.5620652676299368e-02,
	7.6861958028898170e-02, -4.4376272158426051e-02, 4.4376272158426058e-02, -2.5620652676299389e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        2.3681559834188498e-01,  1.3672554945099072e-01,  1.3672554945098905e-01,  7.8938532780630122e-02,
        6.3388180537571182e-01,  3.6597183096807612e-01,  9.2520732066094527e-02,  5.3416869563979862e-02,
        8.9067264750613684e-01,  5.1423009279750498e-01,  5.5737529763337020e-02,  3.2180077812826891e-02,
        1.0196938687973565e+00,  5.8872052964116539e-01,  1.8752907080322980e-02,  1.0826995950913283e-02,
        1.0196938687973551e+00,  5.8872052964116617e-01, -1.8752907080323782e-02, -1.0826995950912827e-02,
        8.9067264750613595e-01,  5.1423009279750564e-01, -5.5737529763335895e-02, -3.2180077812827390e-02,
        6.3388180537571470e-01,  3.6597183096807495e-01, -9.2520732066093542e-02, -5.3416869563980424e-02,
        2.3681559834188717e-01,  1.3672554945099089e-01, -1.3672554945099044e-01, -7.8938532780628831e-02,
                                                        
        6.3388180537571459e-01,  9.2520732066094485e-02,  3.6597183096807473e-01,  5.3416869563980014e-02,
        1.7090814385197843e+00,  2.5479496666024115e-01,  2.5479496666024093e-01,  4.0272203461491046e-02,
        2.4167368175137374e+00,  3.6684346655702615e-01,  1.5377005689540441e-01,  2.4419028116211553e-02,
        2.7726579759463506e+00,  4.2335376950105846e-01,  5.1721119735948180e-02,  8.2072105005127904e-03,
        2.7726579759463510e+00,  4.2335376950105857e-01, -5.1721119735948097e-02, -8.2072105005127401e-03,
        2.4167368175137383e+00,  3.6684346655702638e-01, -1.5377005689540432e-01, -2.4419028116211581e-02,
        1.7090814385197852e+00,  2.5479496666024121e-01, -2.5479496666024098e-01, -4.0272203461491102e-02,
        6.3388180537571515e-01,  9.2520732066093431e-02, -3.6597183096807484e-01, -5.3416869563980611e-02,
                                                        
        8.9067264750613706e-01,  5.5737529763335694e-02,  5.1423009279750498e-01,  3.2180077812827605e-02,
        2.4167368175137374e+00,  1.5377005689540435e-01,  3.6684346655702599e-01,  2.4419028116211518e-02,
        3.4373208199425331e+00,  2.2239098197586474e-01,  2.2239098197586490e-01,  1.5199281451033411e-02,
        3.9521264825956046e+00,  2.5761269031357442e-01,  7.4832205937229848e-02,  5.1359813390617843e-03,
        3.9521264825956055e+00,  2.5761269031357470e-01, -7.4832205937229418e-02, -5.1359813390615995e-03,
        3.4373208199425345e+00,  2.2239098197586496e-01, -2.2239098197586502e-01, -1.5199281451033622e-02,
        2.4167368175137383e+00,  1.5377005689540424e-01, -3.6684346655702615e-01, -2.4419028116211505e-02,
        8.9067264750613750e-01,  5.5737529763336638e-02, -5.1423009279750520e-01, -3.2180077812827015e-02,
                                                        
        1.0196938687973565e+00,  1.8752907080324045e-02,  5.8872052964116561e-01,  1.0826995950912709e-02,
        2.7726579759463501e+00,  5.1721119735948284e-02,  4.2335376950105824e-01,  8.2072105005128095e-03,
        3.9521264825956042e+00,  7.4832205937229584e-02,  2.5761269031357442e-01,  5.1359813390616359e-03,
        4.5485900655688685e+00,  8.6755719877848164e-02,  8.6755719877848719e-02,  1.7480626442406752e-03,
        4.5485900655688685e+00,  8.6755719877847776e-02, -8.6755719877848275e-02, -1.7480626442409205e-03,
        3.9521264825956046e+00,  7.4832205937228988e-02, -2.5761269031357459e-01, -5.1359813390615587e-03,
        2.7726579759463506e+00,  5.1721119735947924e-02, -4.2335376950105819e-01, -8.2072105005127505e-03,
        1.0196938687973565e+00,  1.8752907080322814e-02, -5.8872052964116572e-01, -1.0826995950913227e-02,
                                                        
        1.0196938687973578e+00, -1.8752907080323244e-02,  5.8872052964116495e-01, -1.0826995950913076e-02,
        2.7726579759463510e+00, -5.1721119735948159e-02,  4.2335376950105824e-01, -8.2072105005127211e-03,
        3.9521264825956046e+00, -7.4832205937229349e-02,  2.5761269031357431e-01, -5.1359813390616368e-03,
        4.5485900655688685e+00, -8.6755719877848525e-02,  8.6755719877847970e-02, -1.7480626442410686e-03,
        4.5485900655688667e+00, -8.6755719877848803e-02, -8.6755719877848247e-02,  1.7480626442408095e-03,
        3.9521264825956033e+00, -7.4832205937229904e-02, -2.5761269031357409e-01,  5.1359813390617765e-03,
        2.7726579759463501e+00, -5.1721119735948243e-02, -4.2335376950105824e-01,  8.2072105005128078e-03,
        1.0196938687973569e+00, -1.8752907080322453e-02, -5.8872052964116517e-01,  1.0826995950913517e-02,
                                                        
        8.9067264750613750e-01, -5.5737529763336964e-02,  5.1423009279750509e-01, -3.2180077812826779e-02,
        2.4167368175137387e+00, -1.5377005689540407e-01,  3.6684346655702621e-01, -2.4419028116211449e-02,
        3.4373208199425340e+00, -2.2239098197586465e-01,  2.2239098197586460e-01, -1.5199281451033603e-02,
        3.9521264825956037e+00, -2.5761269031357431e-01,  7.4832205937229085e-02, -5.1359813390615770e-03,
        3.9521264825956024e+00, -2.5761269031357420e-01, -7.4832205937229654e-02,  5.1359813390616099e-03,
        3.4373208199425314e+00, -2.2239098197586477e-01, -2.2239098197586432e-01,  1.5199281451033496e-02,
        2.4167368175137374e+00, -1.5377005689540416e-01, -3.6684346655702588e-01,  2.4419028116211574e-02,
        8.9067264750613806e-01, -5.5737529763336978e-02, -5.1423009279750453e-01,  3.2180077812826766e-02,
                                                        
        6.3388180537571459e-01, -9.2520732066093431e-02,  3.6597183096807528e-01, -5.3416869563980569e-02,
        1.7090814385197859e+00, -2.5479496666024121e-01,  2.5479496666024121e-01, -4.0272203461491109e-02,
        2.4167368175137387e+00, -3.6684346655702615e-01,  1.5377005689540399e-01, -2.4419028116211421e-02,
        2.7726579759463506e+00, -4.2335376950105813e-01,  5.1721119735947743e-02, -8.2072105005127610e-03,
        2.7726579759463492e+00, -4.2335376950105796e-01, -5.1721119735948368e-02,  8.2072105005128824e-03,
        2.4167368175137369e+00, -3.6684346655702577e-01, -1.5377005689540402e-01,  2.4419028116211414e-02,
        1.7090814385197848e+00, -2.5479496666024104e-01, -2.5479496666024082e-01,  4.0272203461491039e-02,
        6.3388180537571470e-01, -9.2520732066093667e-02, -3.6597183096807501e-01,  5.3416869563980458e-02,
                                                        
        2.3681559834188717e-01, -1.3672554945099052e-01,  1.3672554945099086e-01, -7.8938532780628887e-02,
        6.3388180537571537e-01, -3.6597183096807490e-01,  9.2520732066093569e-02, -5.3416869563980528e-02,
        8.9067264750613795e-01, -5.1423009279750498e-01,  5.5737529763336596e-02, -3.2180077812827008e-02,
        1.0196938687973571e+00, -5.8872052964116517e-01,  1.8752907080323126e-02, -1.0826995950912990e-02,
        1.0196938687973565e+00, -5.8872052964116517e-01, -1.8752907080323424e-02,  1.0826995950912984e-02,
        8.9067264750613695e-01, -5.1423009279750498e-01, -5.5737529763336513e-02,  3.2180077812827092e-02,
        6.3388180537571459e-01, -3.6597183096807484e-01, -9.2520732066093569e-02,  5.3416869563980382e-02,
        2.3681559834188676e-01, -1.3672554945099064e-01, -1.3672554945099064e-01,  7.8938532780628942e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        3.3943303570130012e-01,  1.9597175453399679e-01, -2.5382106732311723e-02, -1.4654366154499068e-02,
        2.5852815051768180e-01,  1.4926129729448101e-01, -2.1328350507206096e-02, -1.2313928906705882e-02,
        2.3704347887786692e-01,  1.3685711633978268e-01,  8.9241695525082501e-03,  5.1523716934339632e-03,
        2.9838358920906116e-01,  1.7227184555161776e-01,  2.6490559659328117e-02,  1.5294331750296710e-02,
        2.7855144864018777e-01,  1.6082175385557235e-01, -3.7940651355373674e-02, -2.1905045273254424e-02,
        2.0381615307077769e-01,  1.1767331084060774e-01, -5.2077916595914715e-03, -3.0067199165486271e-03,
        3.1270423370037509e-01,  1.8053987350364709e-01,  6.8074354322631203e-02,  3.9302746793080885e-02,
        4.0700430412220923e-01,  2.3498404454629240e-01, -1.3630183279984597e-02, -7.8693899858035734e-03,
                                                        
        9.2471721182926658e-01,  1.4194222210591301e-01, -5.7372732217365302e-02, -3.8154300808411140e-03,
        7.5028336345656188e-01,  1.3465370727118317e-01, -4.3336697096398154e-02, -3.9259592098273228e-04,
        7.2217959325319470e-01,  1.4323634988842188e-01,  2.7110977804240665e-02,  5.3477869464038522e-03,
        8.6863423587078081e-01,  1.5696251880408427e-01,  5.7444649535093165e-02,  2.5770203719961103e-03,
        8.2578551259476063e-01,  1.5512398027833324e-01, -8.2183371452935167e-02, -3.6385010847539575e-03,
        6.5796549632452850e-01,  1.4452993473923709e-01, -1.4707560116090753e-02, -2.4779739590503548e-03,
        8.7197741194623291e-01,  1.4235664650713659e-01,  1.3826739721006934e-01,  1.2232254132204976e-03,
        1.0677766339715837e+00,  1.4651303796530019e-01, -2.5222663666613812e-02,  1.1764683140077245e-03,
                                                        
        1.3236042788217255e+00,  8.8355333398443453e-02, -6.3319731576856511e-02,  3.8192839976832317e-04,
        1.1500636120564152e+00,  9.6159526874637705e-02, -3.6874019095442509e-02,  4.1238248045203435e-03,
        1.1823773590563609e+00,  1.2245895407850058e-01,  5.5530369624386594e-02,  1.1060156571162853e-02,
        1.4109181535507069e+00,  1.5612524704562614e-01,  7.6417719597735217e-02,  8.3770867360238831e-03,
        1.3654745745531423e+00,  1.5646564492203485e-01, -1.0265458216491942e-01, -8.1805579304477273e-03,
        1.1227306086903228e+00,  1.2380232799510095e-01, -3.7493711870136587e-02, -1.0677616889943637e-02,
        1.2888335450678960e+00,  9.8315354164002000e-02,  1.3339328690091568e-01, -4.0372943122702881e-03,
        1.4765773832016793e+00,  8.9508184647616132e-02, -2.4999331415682426e-02, -1.0475273788137696e-03,
                                                        
        1.5295214682804199e+00,  3.0531011366304865e-02, -6.1912475511358124e-02,  4.3055126846589144e-04,
        1.3790145293808918e+00,  3.6025346873860498e-02, -2.4982746177190984e-02,  2.7416048158394697e-03,
        1.4901027388499875e+00,  5.5206376781829561e-02,  8.9119553817967004e-02,  8.3325679648607077e-03,
        1.8232429688045224e+00,  8.1930596034726758e-02,  1.0321904762417754e-01,  7.0966672146820734e-03,
        1.7786475086544757e+00,  8.2079859803237842e-02, -1.2896624854644198e-01, -7.0104897377519498e-03,
        1.4335217115489933e+00,  5.5632998902093111e-02, -7.0292223316692631e-02, -8.2586125227445510e-03,
        1.5224137527023140e+00,  3.6542241591097972e-02,  1.2161406720537776e-01, -2.7634413497923073e-03,
        1.6849062587527901e+00,  3.0770547731789839e-02, -2.7798975095838528e-02, -5.6884765355937322e-04,
                                                        
        1.5295214682804197e+00, -3.0531011366304792e-02, -6.1912475511358096e-02, -4.3055126846587653e-04,
        1.3790145293808918e+00, -3.6025346873860394e-02, -2.4982746177190963e-02, -2.7416048158394866e-03,
        1.4901027388499877e+00, -5.5206376781829401e-02,  8.9119553817967018e-02, -8.3325679648606973e-03,
        1.8232429688045226e+00, -8.1930596034726688e-02,  1.0321904762417747e-01, -7.0966672146820916e-03,
        1.7786475086544760e+00, -8.2079859803237690e-02, -1.2896624854644190e-01,  7.0104897377519706e-03,
        1.4335217115489931e+00, -5.5632998902093125e-02, -7.0292223316692812e-02,  8.2586125227444539e-03,
        1.5224137527023138e+00, -3.6542241591098125e-02,  1.2161406720537783e-01,  2.7634413497923224e-03,
        1.6849062587527901e+00, -3.0770547731789867e-02, -2.7798975095838487e-02,  5.6884765355936563e-04,
                                                        
        1.3236042788217257e+00, -8.8355333398443481e-02, -6.3319731576856414e-02, -3.8192839976835028e-04,
        1.1500636120564156e+00, -9.6159526874637732e-02, -3.6874019095442509e-02, -4.1238248045203391e-03,
        1.1823773590563613e+00, -1.2245895407850059e-01,  5.5530369624386615e-02, -1.1060156571162844e-02,
        1.4109181535507074e+00, -1.5612524704562611e-01,  7.6417719597735093e-02, -8.3770867360238987e-03,
        1.3654745745531423e+00, -1.5646564492203499e-01, -1.0265458216491942e-01,  8.1805579304476995e-03,
        1.1227306086903228e+00, -1.2380232799510095e-01, -3.7493711870136739e-02,  1.0677616889943734e-02,
        1.2888335450678954e+00, -9.8315354164002056e-02,  1.3339328690091565e-01,  4.0372943122701779e-03,
        1.4765773832016789e+00, -8.9508184647616243e-02, -2.4999331415682322e-02,  1.0475273788138081e-03,
                                                        
        9.2471721182926658e-01, -1.4194222210591306e-01, -5.7372732217365198e-02,  3.8154300808411362e-03,
        7.5028336345656221e-01, -1.3465370727118314e-01, -4.3336697096398057e-02,  3.9259592098277505e-04,
        7.2217959325319514e-01, -1.4323634988842177e-01,  2.7110977804240696e-02, -5.3477869464038444e-03,
        8.6863423587078115e-01, -1.5696251880408424e-01,  5.7444649535093005e-02, -2.5770203719961463e-03,
        8.2578551259476063e-01, -1.5512398027833330e-01, -8.2183371452935264e-02,  3.6385010847539753e-03,
        6.5796549632452839e-01, -1.4452993473923711e-01, -1.4707560116090757e-02,  2.4779739590503531e-03,
        8.7197741194623257e-01, -1.4235664650713659e-01,  1.3826739721006917e-01, -1.2232254132204556e-03,
        1.0677766339715833e+00, -1.4651303796530019e-01, -2.5222663666613656e-02, -1.1764683140077696e-03,
                                                        
        3.3943303570129962e-01, -1.9597175453399704e-01, -2.5382106732311140e-02,  1.4654366154499322e-02,
        2.5852815051768296e-01, -1.4926129729448057e-01, -2.1328350507205690e-02,  1.2313928906706016e-02,
        2.3704347887786834e-01, -1.3685711633978229e-01,  8.9241695525079847e-03, -5.1523716934341332e-03,
        2.9838358920906127e-01, -1.7227184555161790e-01,  2.6490559659327607e-02, -1.5294331750296876e-02,
        2.7855144864018744e-01, -1.6082175385557249e-01, -3.7940651355373403e-02,  2.1905045273254604e-02,
        2.0381615307077783e-01, -1.1767331084060755e-01, -5.2077916595914967e-03,  3.0067199165486245e-03,
        3.1270423370037476e-01, -1.8053987350364703e-01,  6.8074354322630953e-02, -3.9302746793080955e-02,
        4.0700430412220806e-01, -2.3498404454629279e-01, -1.3630183279984807e-02,  7.8693899858033965e-03,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against dirichletx-periodicy,
      // rotated by pgkyl):
      const double sol[256] = {
         3.3943303570129929e-01, -2.5382106732310627e-02,  1.9597175453399684e-01, -1.4654366154499540e-02,
         9.2471721182926570e-01, -5.7372732217365149e-02,  1.4194222210591295e-01, -3.8154300808411951e-03,
         1.3236042788217244e+00, -6.3319731576856414e-02,  8.8355333398443286e-02,  3.8192839976837950e-04,
         1.5295214682804186e+00, -6.1912475511358069e-02,  3.0531011366304959e-02,  4.3055126846582503e-04,
         1.5295214682804188e+00, -6.1912475511358138e-02, -3.0531011366304928e-02, -4.3055126846586200e-04,
         1.3236042788217248e+00, -6.3319731576856483e-02, -8.8355333398443120e-02, -3.8192839976834253e-04,
         9.2471721182926625e-01, -5.7372732217365198e-02, -1.4194222210591306e-01,  3.8154300808411765e-03,
         3.3943303570129940e-01, -2.5382106732310755e-02, -1.9597175453399698e-01,  1.4654366154499521e-02,
                                                         
         2.5852815051768291e-01, -2.1328350507206047e-02,  1.4926129729448026e-01, -1.2313928906705863e-02,
         7.5028336345656133e-01, -4.3336697096398091e-02,  1.3465370727118300e-01, -3.9259592098274870e-04,
         1.1500636120564143e+00, -3.6874019095442460e-02,  9.6159526874637649e-02,  4.1238248045203573e-03,
         1.3790145293808904e+00, -2.4982746177191012e-02,  3.6025346873860394e-02,  2.7416048158394237e-03,
         1.3790145293808904e+00, -2.4982746177190981e-02, -3.6025346873860421e-02, -2.7416048158394055e-03,
         1.1500636120564147e+00, -3.6874019095442495e-02, -9.6159526874637552e-02, -4.1238248045204129e-03,
         7.5028336345656177e-01, -4.3336697096398126e-02, -1.3465370727118303e-01,  3.9259592098280421e-04,
         2.5852815051768296e-01, -2.1328350507205975e-02, -1.4926129729448043e-01,  1.2313928906705860e-02,
                                                         
         2.3704347887786717e-01,  8.9241695525077245e-03,  1.3685711633978231e-01,  5.1523716934341627e-03,
         7.2217959325319403e-01,  2.7110977804240554e-02,  1.4323634988842177e-01,  5.3477869464038947e-03,
         1.1823773590563600e+00,  5.5530369624386580e-02,  1.2245895407850053e-01,  1.1060156571162867e-02,
         1.4901027388499859e+00,  8.9119553817966893e-02,  5.5206376781829269e-02,  8.3325679648606505e-03,
         1.4901027388499863e+00,  8.9119553817966962e-02, -5.5206376781829172e-02, -8.3325679648606332e-03,
         1.1823773590563604e+00,  5.5530369624386740e-02, -1.2245895407850046e-01, -1.1060156571162811e-02,
         7.2217959325319470e-01,  2.7110977804240727e-02, -1.4323634988842171e-01, -5.3477869464039407e-03,
         2.3704347887786764e-01,  8.9241695525078789e-03, -1.3685711633978243e-01, -5.1523716934341384e-03,
                                                         
         2.9838358920905961e-01,  2.6490559659327586e-02,  1.7227184555161815e-01,  1.5294331750296947e-02,
         8.6863423587077970e-01,  5.7444649535092998e-02,  1.5696251880408424e-01,  2.5770203719960869e-03,
         1.4109181535507060e+00,  7.6417719597735009e-02,  1.5612524704562616e-01,  8.3770867360239074e-03,
         1.8232429688045209e+00,  1.0321904762417757e-01,  8.1930596034726591e-02,  7.0966672146821827e-03,
         1.8232429688045211e+00,  1.0321904762417759e-01, -8.1930596034726424e-02, -7.0966672146822009e-03,
         1.4109181535507065e+00,  7.6417719597734982e-02, -1.5612524704562597e-01, -8.3770867360238883e-03,
         8.6863423587078070e-01,  5.7444649535092970e-02, -1.5696251880408427e-01, -2.5770203719961051e-03,
         2.9838358920906094e-01,  2.6490559659327881e-02, -1.7227184555161787e-01, -1.5294331750296758e-02,
                                                         
         2.7855144864018666e-01, -3.7940651355372897e-02,  1.6082175385557240e-01, -2.1905045273254868e-02,
         8.2578551259475919e-01, -8.2183371452935167e-02,  1.5512398027833305e-01, -3.6385010847539805e-03,
         1.3654745745531405e+00, -1.0265458216491943e-01,  1.5646564492203499e-01, -8.1805579304476978e-03,
         1.7786475086544749e+00, -1.2896624854644190e-01,  8.2079859803237828e-02, -7.0104897377519255e-03,
         1.7786475086544753e+00, -1.2896624854644170e-01, -8.2079859803237634e-02,  7.0104897377520183e-03,
         1.3654745745531418e+00, -1.0265458216491913e-01, -1.5646564492203482e-01,  8.1805579304476423e-03,
         8.2578551259476041e-01, -8.2183371452935056e-02, -1.5512398027833324e-01,  3.6385010847539249e-03,
         2.7855144864018744e-01, -3.7940651355373500e-02, -1.6082175385557240e-01,  2.1905045273254500e-02,
                                                         
         2.0381615307077769e-01, -5.2077916595916198e-03,  1.1767331084060717e-01, -3.0067199165485217e-03,
         6.5796549632452717e-01, -1.4707560116090692e-02,  1.4452993473923695e-01, -2.4779739590503344e-03,
         1.1227306086903213e+00, -3.7493711870136663e-02,  1.2380232799510092e-01, -1.0677616889943762e-02,
         1.4335217115489918e+00, -7.0292223316692895e-02,  5.5632998902093299e-02, -8.2586125227444972e-03,
         1.4335217115489924e+00, -7.0292223316692923e-02, -5.5632998902092813e-02,  8.2586125227444782e-03,
         1.1227306086903228e+00, -3.7493711870136760e-02, -1.2380232799510088e-01,  1.0677616889943752e-02,
         6.5796549632452839e-01, -1.4707560116090772e-02, -1.4452993473923714e-01,  2.4779739590503622e-03,
         2.0381615307077774e-01, -5.2077916595914247e-03, -1.1767331084060761e-01,  3.0067199165486531e-03,
                                                         
         3.1270423370037409e-01,  6.8074354322630801e-02,  1.8053987350364703e-01,  3.9302746793081066e-02,
         8.7197741194623157e-01,  1.3826739721006912e-01,  1.4235664650713634e-01,  1.2232254132204063e-03,
         1.2888335450678943e+00,  1.3339328690091576e-01,  9.8315354164002167e-02, -4.0372943122700409e-03,
         1.5224137527023127e+00,  1.2161406720537798e-01,  3.6542241591098042e-02, -2.7634413497924326e-03,
         1.5224137527023132e+00,  1.2161406720537785e-01, -3.6542241591097750e-02,  2.7634413497923402e-03,
         1.2888335450678956e+00,  1.3339328690091576e-01, -9.8315354164001792e-02,  4.0372943122702351e-03,
         8.7197741194623291e-01,  1.3826739721006936e-01, -1.4235664650713664e-01, -1.2232254132204989e-03,
         3.1270423370037487e-01,  6.8074354322631009e-02, -1.8053987350364709e-01, -3.9302746793080982e-02,
                                                         
         4.0700430412220695e-01, -1.3630183279984911e-02,  2.3498404454629299e-01, -7.8693899858033878e-03,
         1.0677766339715822e+00, -2.5222663666613573e-02,  1.4651303796530005e-01,  1.1764683140078709e-03,
         1.4765773832016775e+00, -2.4999331415682360e-02,  8.9508184647616285e-02, -1.0475273788140100e-03,
         1.6849062587527890e+00, -2.7798975095838598e-02,  3.0770547731789843e-02, -5.6884765355922739e-04,
         1.6849062587527892e+00, -2.7798975095838566e-02, -3.0770547731789812e-02,  5.6884765355924582e-04,
         1.4765773832016786e+00, -2.4999331415682616e-02, -8.9508184647615868e-02,  1.0475273788138251e-03,
         1.0677766339715835e+00, -2.5222663666613909e-02, -1.4651303796530027e-01, -1.1764683140077321e-03,
         4.0700430412220773e-01, -1.3630183279985131e-02, -2.3498404454629301e-01,  7.8693899858033340e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_NEUMANN && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
         2.7426849265785821e-03,  1.5834898806625446e-03, -7.6901868952472112e-05, -4.4399314740841794e-05,
         2.2307627582822381e-03,  1.2879314789925147e-03, -2.1865653271755913e-04, -1.2624140802460550e-04,
         1.2904511538739966e-03,  7.4504232106521533e-04, -3.2423262520966782e-04, -1.8719579344485024e-04,
         8.9060811742788660e-05,  5.1419283634405537e-05, -3.6939041222160184e-04, -2.1326765393195373e-04,
        -1.1216934698784315e-03, -6.4761002678181273e-04, -3.2963889819490169e-04, -1.9031710660828843e-04,
        -2.0042679072861938e-03, -1.1571646157994585e-03, -1.7991569082235231e-04, -1.0387437252785848e-04,
        -2.1317958699093487e-03, -1.2307929193492778e-03,  1.0628738727266920e-04,  6.1365051653394852e-05,
        -9.7385035747464554e-04, -5.6225276603835304e-04,  5.6225276603842590e-04,  3.2461678582474123e-04,
                                                         
         7.0025973717426759e-03,  8.7597171627718244e-04, -1.9640674915166975e-04, -2.4596860011638872e-05,
         5.6952244942972364e-03,  7.1227643689290312e-04, -5.5840533357267067e-04, -6.9912653605945670e-05,
         3.2941110650197636e-03,  4.1177126799925512e-04, -8.2787815117550407e-04, -1.0358408654767605e-04,
         2.2715075512227500e-04,  2.8306982348121067e-05, -9.4283220933770928e-04, -1.1780912199761380e-04,
        -2.8619233323846529e-03, -3.5711215278797967e-04, -8.4064554663111309e-04, -1.0471271942404679e-04,
        -5.1101105661364738e-03, -6.3599447934837444e-04, -4.5734595796418119e-04, -5.6300066887826217e-05,
        -5.4273249985191275e-03, -6.7188171017583094e-04,  2.7420212003722227e-04,  3.5580564509123963e-05,
        -2.4761964975358519e-03, -3.0512718234197955e-04,  1.4296327144186487e-03,  1.7616526086226424e-04,
                                                         
         8.9364896785849718e-03,  2.4056152766196369e-04, -2.5073307216640474e-04, -6.7684572050010251e-06,
         7.2676304559713464e-03,  1.9555256835501660e-04, -7.1278324924916182e-04, -1.9217477566809790e-05,
         4.2031648967955498e-03,  1.1307120648435601e-04, -1.0564867662633216e-03, -2.8403158912342905e-05,
         2.9035326851860549e-04,  8.1830057747113072e-06, -1.2025760806106706e-03, -3.2154072002185963e-05,
        -3.6477085693248867e-03, -9.6561165284698010e-05, -1.0710649815536620e-03, -2.8320003355008237e-05,
        -6.5084005513759806e-03, -1.7130862003481248e-04, -5.8055630435214156e-04, -1.4835459766209707e-05,
        -6.9026112131815769e-03, -1.7987518299127813e-04,  3.5295867264126013e-04,  9.8895523372643026e-06,
        -3.1456344295972740e-03, -8.1372987939429096e-05,  1.8161328847003111e-03,  4.6980716491689845e-05,
                                                         
         8.8699248778118141e-03, -2.7899273330690011e-04, -2.4891211971957389e-04,  7.8197845903604755e-06,
         7.2133278602044600e-03, -2.2690418663872376e-04, -7.0752461433476357e-04,  2.2253551850207463e-05,
         4.1719433931248044e-03, -1.3109695003316613e-04, -1.0484195264428224e-03,  3.3060781994325324e-05,
         2.9000061244580100e-04, -8.3866118532763377e-06, -1.1928211829609248e-03,  3.7786064786184074e-05,
        -3.6144611318088106e-03,  1.1575658228446954e-04, -1.0614208561250714e-03,  3.3888041767303515e-05,
        -6.4470455480175864e-03,  2.0673194773991214e-04, -5.7397251907540382e-04,  1.8636609968020875e-05,
        -6.8330472413997891e-03,  2.2003796082539271e-04,  3.5111433749353732e-04, -1.0954379731277978e-05,
        -3.1124496847974503e-03,  1.0053220928259660e-04,  1.7969736633570872e-03, -5.8042298091663899e-05,
                                                         
         7.2812441275884928e-03, -6.3823252549091139e-04, -2.0435956575817752e-04,  1.7902644432336832e-05,
         5.9213079724601374e-03, -5.1904384334881132e-04, -5.8079993948588330e-04,  5.0910973286761385e-05,
         3.4251487633341262e-03, -3.0006513050584904e-04, -8.6035825184318300e-04,  7.5516445519920679e-05,
         2.4050487889585632e-04, -2.0189763235522183e-05, -9.7829675211034734e-04,  8.6069673113147962e-05,
        -2.9602655566508232e-03,  2.6194340913565666e-04, -8.6966892046672398e-04,  7.6819989902675276e-05,
        -5.2793941934909418e-03,  4.6741187611567658e-04, -4.6928062229828889e-04,  4.1807284818218978e-05,
        -5.5925404103385045e-03,  4.9616899201940950e-04,  2.8848556970563607e-04, -2.5204356210081428e-05,
        -2.5464343731789914e-03,  2.2625688324577048e-04,  1.4701845708285533e-03, -1.3062947244791865e-04,
                                                         
         4.8041572544600909e-03, -7.9191424751586025e-04, -1.3488355963031778e-04,  2.2209346407803159e-05,
         3.9067216146501721e-03, -6.4407813263852004e-04, -3.8325114859430275e-04,  6.3143874312577674e-05,
         2.2601283795129056e-03, -3.7255970170473366e-04, -5.6740989895801485e-04,  9.3617364876988726e-05,
         1.6096491402340063e-04, -2.5732656894923494e-05, -6.4454269291537862e-04,  1.0662332280619761e-04,
        -1.9460010749482062e-03,  3.2364246240482560e-04, -5.7191468799076492e-04,  9.5088496369670188e-05,
        -3.4686181223184866e-03,  5.7804017601787404e-04, -3.0716867418117987e-04,  5.1788092066048652e-05,
        -3.6697476510151273e-03,  6.1395592517180572e-04,  1.9104648664619163e-04, -3.1052124623879169e-05,
        -1.6694227147681959e-03,  2.8008603382014246e-04,  9.6384165376269546e-04, -1.6170774702231475e-04,
                                                         
         2.2306540994137579e-03, -6.9389849181049129e-04, -6.2708495935571388e-05,  1.9460946045137133e-05,
         1.8135960357380603e-03, -5.6438848379407899e-04, -1.7808008939529708e-04,  5.5311691945889000e-05,
         1.0490629032274506e-03, -3.2664927703493615e-04, -2.6332332046409134e-04,  8.1947103073430277e-05,
         7.6120474946445585e-05, -2.3252302845399339e-05, -2.9840525240995520e-04,  9.3219221646216737e-05,
        -8.9710993891016207e-04,  2.8193511733735156e-04, -2.6348958901368814e-04,  8.2980817549581983e-05,
        -1.5948876788955354e-03,  5.0375859988116372e-04, -1.3937257700139577e-04,  4.5089029809668238e-05,
        -1.6783833458873889e-03,  5.3575879244830630e-04,  9.1166331187481845e-05, -2.6613710016908576e-05,
        -7.6023931416550959e-04,  2.4483124726055788e-04,  4.3892437268199354e-04, -1.4135338651191499e-04,
                                                         
         5.1439332815128255e-04, -2.9698512647749389e-04, -1.4500574311018977e-05,  8.3719104818705462e-06,
         4.1802325329997968e-04, -2.4134583782026186e-04, -4.1138714346205186e-05,  2.3751447801901673e-05,
         2.4164487957365318e-04, -1.3951373627013949e-04, -6.0693387203917837e-05,  3.5041343440206787e-05,
         1.7923152512615924e-05, -1.0347936927887651e-05, -6.8472412138340007e-05,  3.9532565580135514e-05,
        -2.0439199562198769e-04,  1.1800577369255985e-04, -5.9881298482104929e-05,  3.4572483798067180e-05,
        -3.6117609457579966e-04,  2.0852511542819542e-04, -3.0638043253531359e-05,  1.7688882519869748e-05,
        -3.7521094838258683e-04,  2.1662814205158118e-04,  2.2535016630145747e-05, -1.3010597917607269e-05,
        -1.6808957731488241e-04,  9.7046562710717776e-05,  9.7046562710717776e-05, -5.6029859104960818e-05,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_NEUMANN)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        -9.7385035747442393e-04, -5.6225276603837191e-04, -5.6225276603831217e-04, -3.2461678582478790e-04,
        -2.1317958699089319e-03, -1.2307929193493136e-03, -1.0628738727267024e-04, -6.1365051653357894e-05,
        -2.0042679072855393e-03, -1.1571646157995700e-03,  1.7991569082249033e-04,  1.0387437252777771e-04,
        -1.1216934698778211e-03, -6.4761002678188971e-04,  3.2963889819473852e-04,  1.9031710660838912e-04,
         8.9060811743257198e-05,  5.1419283634326885e-05,  3.6939041222168315e-04,  2.1326765393185208e-04,
         1.2904511538740504e-03,  7.4504232106528234e-04,  3.2423262520934700e-04,  1.8719579344503596e-04,
         2.2307627582821917e-03,  1.2879314789926051e-03,  2.1865653271782218e-04,  1.2624140802443317e-04,
         2.7426849265788397e-03,  1.5834898806624822e-03,  7.6901868952384753e-05,  4.4399314740925765e-05,
                                                         
        -2.4761964975356984e-03, -3.0512718234200026e-04, -1.4296327144186216e-03, -1.7616526086226754e-04,
        -5.4273249985188473e-03, -6.7188171017587420e-04, -2.7420212003717614e-04, -3.5580564509133681e-05,
        -5.1101105661361026e-03, -6.3599447934842561e-04,  4.5734595796418823e-04,  5.6300066887831387e-05,
        -2.8619233323842678e-03, -3.5711215278803258e-04,  8.4064554663111364e-04,  1.0471271942404063e-04,
         2.2715075512256519e-04,  2.8306982348096778e-05,  9.4283220933765398e-04,  1.1780912199763648e-04,
         3.2941110650199410e-03,  4.1177126799925935e-04,  8.2787815117549421e-04,  1.0358408654766980e-04,
         5.6952244942973639e-03,  7.1227643689291298e-04,  5.5840533357265180e-04,  6.9912653605955252e-05,
         7.0025973717428164e-03,  8.7597171627717746e-04,  1.9640674915169629e-04,  2.4596860011620800e-05,
                                                         
        -3.1456344295972402e-03, -8.1372987939477303e-05, -1.8161328847002679e-03, -4.6980716491677343e-05,
        -6.9026112131814199e-03, -1.7987518299130643e-04, -3.5295867264123259e-04, -9.8895523372652428e-06,
        -6.5084005513757525e-03, -1.7130862003484438e-04,  5.8055630435215532e-04,  1.4835459766208479e-05,
        -3.6477085693246425e-03, -9.6561165284726890e-05,  1.0710649815536570e-03,  2.8320003355011235e-05,
         2.9035326851882027e-04,  8.1830057746921219e-06,  1.2025760806106591e-03,  3.2154072002188565e-05,
         4.2031648967957207e-03,  1.1307120648434805e-04,  1.0564867662633073e-03,  2.8403158912346775e-05,
         7.2676304559714851e-03,  1.9555256835501359e-04,  7.1278324924915759e-04,  1.9217477566808923e-05,
         8.9364896785850984e-03,  2.4056152766196043e-04,  2.5073307216640127e-04,  6.7684572050017476e-06,
                                                         
        -3.1124496847974941e-03,  1.0053220928259961e-04, -1.7969736633570369e-03,  5.8042298091655517e-05,
        -6.8330472413997067e-03,  2.2003796082537783e-04, -3.5111433749351466e-04,  1.0954379731276026e-05,
        -6.4470455480174459e-03,  2.0673194773989311e-04,  5.7397251907541488e-04, -1.8636609968021309e-05,
        -3.6144611318086493e-03,  1.1575658228445096e-04,  1.0614208561250729e-03, -3.3888041767302810e-05,
         2.9000061244595664e-04, -8.3866118532913302e-06,  1.1928211829609200e-03, -3.7786064786182719e-05,
         4.1719433931249415e-03, -1.3109695003317735e-04,  1.0484195264428165e-03, -3.3060781994324497e-05,
         7.2133278602045797e-03, -2.2690418663873213e-04,  7.0752461433475978e-04, -2.2253551850206742e-05,
         8.8699248778119251e-03, -2.7899273330690640e-04,  2.4891211971957215e-04, -7.8197845903600418e-06,
                                                         
        -2.5464343731789896e-03,  2.2625688324579403e-04, -1.4701845708285403e-03,  1.3062947244790564e-04,
        -5.5925404103384559e-03,  4.9616899201940538e-04, -2.8848556970562182e-04,  2.5204356210078464e-05,
        -5.2793941934908550e-03,  4.6741187611566514e-04,  4.6928062229829702e-04, -4.1807284818220205e-05,
        -2.9602655566507174e-03,  2.6194340913564322e-04,  8.6966892046672658e-04, -7.6819989902675290e-05,
         2.4050487889596436e-04, -2.0189763235534736e-05,  9.7829675211034582e-04, -8.6069673113147433e-05,
         3.4251487633342255e-03, -3.0006513050585956e-04,  8.6035825184317943e-04, -7.5516445519920001e-05,
         5.9213079724602276e-03, -5.1904384334881999e-04,  5.8079993948588113e-04, -5.0910973286761026e-05,
         7.2812441275885778e-03, -6.3823252549091963e-04,  2.0435956575817714e-04, -1.7902644432336977e-05,
                                                         
        -1.6694227147681640e-03,  2.8008603382013634e-04, -9.6384165376269849e-04,  1.6170774702231849e-04,
        -3.6697476510150909e-03,  6.1395592517180268e-04, -1.9104648664618608e-04,  3.1052124623877217e-05,
        -3.4686181223184324e-03,  5.7804017601786623e-04,  3.0716867418118448e-04, -5.1788092066049479e-05,
        -1.9460010749481403e-03,  3.2364246240481601e-04,  5.7191468799076719e-04, -9.5088496369670378e-05,
         1.6096491402347039e-04, -2.5732656894933102e-05,  6.4454269291537851e-04, -1.0662332280619741e-04,
         2.2601283795129724e-03, -3.7255970170474228e-04,  5.6740989895801334e-04, -9.3617364876988346e-05,
         3.9067216146502345e-03, -6.4407813263852741e-04,  3.8325114859430194e-04, -6.3143874312577308e-05,
         4.8041572544601507e-03, -7.9191424751586664e-04,  1.3488355963031703e-04, -2.2209346407802942e-05,
                                                         
        -7.6023931416549680e-04,  2.4483124726055295e-04, -4.3892437268199062e-04,  1.4135338651191469e-04,
        -1.6783833458873679e-03,  5.3575879244830056e-04, -9.1166331187479962e-05,  2.6613710016908431e-05,
        -1.5948876788955061e-03,  5.0375859988115732e-04,  1.3937257700139864e-04, -4.5089029809668454e-05,
        -8.9710993891012531e-04,  2.8193511733734435e-04,  2.6348958901368971e-04, -8.2980817549582240e-05,
         7.6120474946484617e-05, -2.3252302845407416e-05,  2.9840525240995498e-04, -9.3219221646216954e-05,
         1.0490629032274879e-03, -3.2664927703494450e-04,  2.6332332046409047e-04, -8.1947103073430277e-05,
         1.8135960357380959e-03, -5.6438848379408713e-04,  1.7808008939529699e-04, -5.5311691945888858e-05,
         2.2306540994137948e-03, -6.9389849181049855e-04,  6.2708495935571632e-05, -1.9460946045136808e-05,
                                                         
        -1.6808957731488268e-04,  9.7046562710715215e-05, -9.7046562710717126e-05,  5.6029859104959808e-05,
        -3.7521094838258699e-04,  2.1662814205157473e-04, -2.2535016630146340e-05,  1.3010597917606041e-05,
        -3.6117609457579625e-04,  2.0852511542818696e-04,  3.0638043253534035e-05, -1.7688882519869666e-05,
        -2.0439199562196969e-04,  1.1800577369255624e-04,  5.9881298482110635e-05, -3.4572483798064477e-05,
         1.7923152512632370e-05, -1.0347936927892654e-05,  6.8472412138333421e-05, -3.9532565580139031e-05,
         2.4164487957365534e-04, -1.3951373627015142e-04,  6.0693387203916191e-05, -3.5041343440207282e-05,
         4.1802325329998662e-04, -2.4134583782027048e-04,  4.1138714346209597e-05, -2.3751447801899271e-05,
         5.1439332815129729e-04, -2.9698512647749915e-04,  1.4500574311019040e-05, -8.3719104818709799e-06,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
         2.7426849265786471e-03, -7.6901868952414880e-05,  1.5834898806625407e-03, -4.4399314740884376e-05,
         7.0025973717427374e-03, -1.9640674915167664e-04,  8.7597171627718396e-04, -2.4596860011633376e-05,
         8.9364896785850360e-03, -2.5073307216640149e-04,  2.4056152766196320e-04, -6.7684572050005914e-06,
         8.8699248778118644e-03, -2.4891211971957063e-04, -2.7899273330690813e-04,  7.8197845903598978e-06,
         7.2812441275885067e-03, -2.0435956575817513e-04, -6.3823252549092418e-04,  1.7902644432337120e-05,
         4.8041572544600519e-03, -1.3488355963031391e-04, -7.9191424751587803e-04,  2.2209346407803664e-05,
         2.2306540994136361e-03, -6.2708495935559299e-05, -6.9389849181052132e-04,  1.9460946045141362e-05,
         5.1439332815105606e-04, -1.4500574311105517e-05, -2.9698512647752414e-04,  8.3719104818094158e-06,
                                                         
         2.2307627582822584e-03, -2.1865653271764226e-04,  1.2879314789925184e-03, -1.2624140802455858e-04,
         5.6952244942972980e-03, -5.5840533357266384e-04,  7.1227643689292339e-04, -6.9912653605940683e-05,
         7.2676304559714392e-03, -7.1278324924914881e-04,  1.9555256835501460e-04, -1.9217477566811237e-05,
         7.2133278602045312e-03, -7.0752461433475447e-04, -2.2690418663873414e-04,  2.2253551850206742e-05,
         5.9213079724601695e-03, -5.8079993948587528e-04, -5.1904384334882357e-04,  5.0910973286761385e-05,
         3.9067216146501582e-03, -3.8325114859429212e-04, -6.4407813263853446e-04,  6.3143874312579260e-05,
         1.8135960357379883e-03, -1.7808008939528016e-04, -5.6438848379409851e-04,  5.5311691945890992e-05,
         4.1802325329974165e-04, -4.1138714346125267e-05, -2.4134583782033838e-04,  2.3751447801936089e-05,
                                                         
         1.2904511538738648e-03, -3.2423262520967281e-04,  7.4504232106538328e-04, -1.8719579344480202e-04,
         3.2941110650199284e-03, -8.2787815117545095e-04,  4.1177126799925859e-04, -1.0358408654769062e-04,
         4.2031648967957016e-03, -1.0564867662633006e-03,  1.1307120648434487e-04, -2.8403158912346883e-05,
         4.1719433931249120e-03, -1.0484195264428109e-03, -1.3109695003318017e-04,  3.3060781994323881e-05,
         3.4251487633341878e-03, -8.6035825184317411e-04, -3.0006513050586156e-04,  7.5516445519920543e-05,
         2.2601283795129290e-03, -5.6740989895800401e-04, -3.7255970170474326e-04,  9.3617364876989878e-05,
         1.0490629032274485e-03, -2.6332332046406787e-04, -3.2664927703494130e-04,  8.1947103073436511e-05,
         2.4164487957369405e-04, -6.0693387203836719e-05, -1.3951373627010949e-04,  3.5041343440233865e-05,
                                                         
         8.9060811743351049e-05, -3.6939041222119586e-04,  5.1419283634325340e-05, -2.1326765393214525e-04,
         2.2715075512261935e-04, -9.4283220933765886e-04,  2.8306982348075402e-05, -1.1780912199762746e-04,
         2.9035326851882456e-04, -1.2025760806106524e-03,  8.1830057746845782e-06, -3.2154072002190882e-05,
         2.9000061244594167e-04, -1.1928211829609172e-03, -8.3866118532948674e-06,  3.7786064786182841e-05,
         2.4050487889594224e-04, -9.7829675211034257e-04, -2.0189763235535268e-05,  8.6069673113147745e-05,
         1.6096491402345218e-04, -6.4454269291537309e-04, -2.5732656894930253e-05,  1.0662332280619814e-04,
         7.6120474946493859e-05, -2.9840525240994956e-04, -2.3252302845394443e-05,  9.3219221646216344e-05,
         1.7923152512759621e-05, -6.8472412138361773e-05, -1.0347936927837471e-05,  3.9532565580120085e-05,
                                                         
        -1.1216934698773302e-03, -3.2963889819499645e-04, -6.4761002678212140e-04, -1.9031710660822880e-04,
        -2.8619233323842136e-03, -8.4064554663110865e-04, -3.5711215278805313e-04, -1.0471271942404913e-04,
        -3.6477085693246369e-03, -1.0710649815536626e-03, -9.6561165284734344e-05, -2.8320003355008816e-05,
        -3.6144611318086606e-03, -1.0614208561250736e-03,  1.1575658228444842e-04,  3.3888041767303257e-05,
        -2.9602655566507317e-03, -8.6966892046672539e-04,  2.6194340913564425e-04,  7.6819989902675873e-05,
        -1.9460010749481468e-03, -5.7191468799076578e-04,  3.2364246240481920e-04,  9.5088496369669835e-05,
        -8.9710993891011512e-04, -2.6348958901369470e-04,  2.8193511733735086e-04,  8.2980817549579137e-05,
        -2.0439199562196194e-04, -5.9881298482151265e-05,  1.1800577369254833e-04,  3.4572483798047007e-05,
                                                         
        -2.0042679072857297e-03, -1.7991569082262559e-04, -1.1571646157994585e-03, -1.0387437252773997e-04,
        -5.1101105661361095e-03, -4.5734595796422857e-04, -6.3599447934843125e-04, -5.6300066887814325e-05,
        -6.5084005513757707e-03, -5.8055630435216378e-04, -1.7130862003484528e-04, -1.4835459766207032e-05,
        -6.4470455480174641e-03, -5.7397251907541792e-04,  2.0673194773989412e-04,  1.8636609968022898e-05,
        -5.2793941934908672e-03, -4.6928062229829713e-04,  4.6741187611566747e-04,  4.1807284818220280e-05,
        -3.4686181223184381e-03, -3.0716867418118492e-04,  5.7804017601786786e-04,  5.1788092066049120e-05,
        -1.5948876788955102e-03, -1.3937257700140189e-04,  5.0375859988115623e-04,  4.5089029809667153e-05,
        -3.6117609457582775e-04, -3.0638043253516085e-05,  2.0852511542817227e-04,  1.7688882519883192e-05,
                                                         
        -2.1317958699091887e-03,  1.0628738727276690e-04, -1.2307929193492383e-03,  6.1365051653299117e-05,
        -5.4273249985189401e-03,  2.7420212003716692e-04, -6.7188171017585436e-04,  3.5580564509131336e-05,
        -6.9026112131814685e-03,  3.5295867264122359e-04, -1.7987518299130041e-04,  9.8895523372676992e-06,
        -6.8330472413997388e-03,  3.5111433749350989e-04,  2.2003796082538133e-04, -1.0954379731276097e-05,
        -5.5925404103384663e-03,  2.8848556970562306e-04,  4.9616899201941394e-04, -2.5204356210074923e-05,
        -3.6697476510150926e-03,  1.9104648664618846e-04,  6.1395592517179899e-04, -3.1052124623879962e-05,
        -1.6783833458873766e-03,  9.1166331187480531e-05,  5.3575879244830012e-04, -2.6613710016906733e-05,
        -3.7521094838258639e-04,  2.2535016630146950e-05,  2.1662814205158050e-04, -1.3010597917607730e-05,
                                                         
        -9.7385035747447619e-04,  5.6225276603833363e-04, -5.6225276603832789e-04,  3.2461678582482883e-04,
        -2.4761964975357174e-03,  1.4296327144186732e-03, -3.0512718234202525e-04,  1.7616526086224407e-04,
        -3.1456344295973885e-03,  1.8161328847002194e-03, -8.1372987939527000e-05,  4.6980716491642716e-05,
        -3.1124496847975457e-03,  1.7969736633570307e-03,  1.0053220928270552e-04, -5.8042298091596394e-05,
        -2.5464343731789350e-03,  1.4701845708285763e-03,  2.2625688324574906e-04, -1.3062947244794012e-04,
        -1.6694227147681826e-03,  9.6384165376268624e-04,  2.8008603382013915e-04, -1.6170774702231198e-04,
        -7.6023931416550460e-04,  4.3892437268199062e-04,  2.4483124726055647e-04, -1.4135338651191412e-04,
        -1.6808957731488113e-04,  9.7046562710717045e-05,  9.7046562710717045e-05, -5.6029859104960391e-05,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_NEUMANN) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
        -9.7385035747447923e-04, -5.6225276603807397e-04, -5.6225276603855037e-04, -3.2461678582481864e-04,
        -2.4761964975359234e-03, -1.4296327144185548e-03, -3.0512718234191976e-04, -1.7616526086233577e-04,
        -3.1456344295974384e-03, -1.8161328847002257e-03, -8.1372987939542246e-05, -4.6980716491623221e-05,
        -3.1124496847975899e-03, -1.7969736633570341e-03,  1.0053220928272368e-04,  5.8042298091578539e-05,
        -2.5464343731790439e-03, -1.4701845708285615e-03,  2.2625688324569369e-04,  1.3062947244796877e-04,
        -1.6694227147683789e-03, -9.6384165376262097e-04,  2.8008603382014403e-04,  1.6170774702231231e-04,
        -7.6023931416560110e-04, -4.3892437268195917e-04,  2.4483124726060916e-04,  1.4135338651189423e-04,
        -1.6808957731484405e-04, -9.7046562710695645e-05,  9.7046562710741547e-05,  5.6029859104974540e-05,
                                                         
        -2.1317958699087467e-03, -1.0628738727276966e-04, -1.2307929193494205e-03, -6.1365051653285957e-05,
        -5.4273249985188976e-03, -2.7420212003714215e-04, -6.7188171017590293e-04, -3.5580564509128585e-05,
        -6.9026112131815023e-03, -3.5295867264120820e-04, -1.7987518299129591e-04, -9.8895523372758680e-06,
        -6.8330472413997761e-03, -3.5111433749350236e-04,  2.2003796082537493e-04,  1.0954379731279712e-05,
        -5.5925404103385184e-03, -2.8848556970560593e-04,  4.9616899201941178e-04,  2.5204356210076946e-05,
        -3.6697476510151351e-03, -1.9104648664616480e-04,  6.1395592517180745e-04,  3.1052124623881629e-05,
        -1.6783833458873603e-03, -9.1166331187446853e-05,  5.3575879244832539e-04,  2.6613710016910782e-05,
        -3.7521094838252183e-04, -2.2535016630152493e-05,  2.1662814205158321e-04,  1.3010597917581004e-05,
                                                         
        -2.0042679072855732e-03,  1.7991569082246353e-04, -1.1571646157995127e-03,  1.0387437252780057e-04,
        -5.1101105661360731e-03,  4.5734595796419967e-04, -6.3599447934844697e-04,  5.6300066887830588e-05,
        -6.5084005513757716e-03,  5.8055630435216736e-04, -1.7130862003485040e-04,  1.4835459766209634e-05,
        -6.4470455480174728e-03,  5.7397251907542681e-04,  2.0673194773989449e-04, -1.8636609968022539e-05,
        -5.2793941934908715e-03,  4.6928062229830824e-04,  4.6741187611567023e-04, -4.1807284818219338e-05,
        -3.4686181223184172e-03,  3.0716867418119733e-04,  5.7804017601787903e-04, -5.1788092066049262e-05,
        -1.5948876788954261e-03,  1.3937257700140739e-04,  5.0375859988118171e-04, -4.5089029809671057e-05,
        -3.6117609457557334e-04,  3.0638043253631261e-05,  2.0852511542824513e-04, -1.7688882519815954e-05,
                                                         
        -1.1216934698776346e-03,  3.2963889819489237e-04, -6.4761002678198143e-04,  1.9031710660828038e-04,
        -2.8619233323842431e-03,  8.4064554663110019e-04, -3.5711215278803448e-04,  1.0471271942405270e-04,
        -3.6477085693246277e-03,  1.0710649815536646e-03, -9.6561165284730644e-05,  2.8320003355011381e-05,
        -3.6144611318086415e-03,  1.0614208561250812e-03,  1.1575658228445058e-04, -3.3888041767302519e-05,
        -2.9602655566507013e-03,  8.6966892046673385e-04,  2.6194340913564832e-04, -7.6819989902676090e-05,
        -1.9460010749480965e-03,  5.7191468799077055e-04,  3.2364246240482674e-04, -9.5088496369671895e-05,
        -8.9710993891002904e-04,  2.6348958901369036e-04,  2.8193511733736398e-04, -8.2980817549582389e-05,
        -2.0439199562173984e-04,  5.9881298482017428e-05,  1.1800577369261382e-04, -3.4572483798118524e-05,
                                                         
         8.9060811743205156e-05,  3.6939041222139156e-04,  5.1419283634381339e-05,  2.1326765393204521e-04,
         2.2715075512260449e-04,  9.4283220933767588e-04,  2.8306982348095057e-05,  1.1780912199762439e-04,
         2.9035326851885426e-04,  1.2025760806106624e-03,  8.1830057746907124e-06,  3.2154072002189812e-05,
         2.9000061244598775e-04,  1.1928211829609250e-03, -8.3866118532915809e-06, -3.7786064786182936e-05,
         2.4050487889599971e-04,  9.7829675211034972e-04, -2.0189763235532012e-05, -8.6069673113148070e-05,
         1.6096491402351967e-04,  6.4454269291537830e-04, -2.5732656894927719e-05, -1.0662332280619897e-04,
         7.6120474946557394e-05,  2.9840525240994078e-04, -2.3252302845399292e-05, -9.3219221646223432e-05,
         1.7923152512686993e-05,  6.8472412138325452e-05, -1.0347936927911219e-05, -3.9532565580128941e-05,
                                                         
         1.2904511538740524e-03,  3.2423262520966972e-04,  7.4504232106531823e-04,  1.8719579344483216e-04,
         3.2941110650199987e-03,  8.2787815117548282e-04,  4.1177126799925561e-04,  1.0358408654768070e-04,
         4.2031648967957693e-03,  1.0564867662633123e-03,  1.1307120648434643e-04,  2.8403158912345328e-05,
         4.1719433931249866e-03,  1.0484195264428195e-03, -1.3109695003317754e-04, -3.3060781994324240e-05,
         3.4251487633342710e-03,  8.6035825184318170e-04, -3.0006513050585907e-04, -7.5516445519920679e-05,
         2.2601283795130184e-03,  5.6740989895801171e-04, -3.7255970170474201e-04, -9.3617364876989796e-05,
         1.0490629032275287e-03,  2.6332332046408624e-04, -3.2664927703494792e-04, -8.1947103073430426e-05,
         2.4164487957367882e-04,  6.0693387203906175e-05, -1.3951373627015806e-04, -3.5041343440210460e-05,
                                                         
         2.2307627582824571e-03,  2.1865653271765166e-04,  1.2879314789924958e-03,  1.2624140802455308e-04,
         5.6952244942974333e-03,  5.5840533357266959e-04,  7.1227643689290886e-04,  6.9912653605944084e-05,
         7.2676304559715441e-03,  7.1278324924915857e-04,  1.9555256835501161e-04,  1.9217477566810153e-05,
         7.2133278602046344e-03,  7.0752461433476173e-04, -2.2690418663873251e-04, -2.2253551850207104e-05,
         5.9213079724602779e-03,  5.8079993948588189e-04, -5.1904384334882151e-04, -5.0910973286761460e-05,
         3.9067216146502770e-03,  3.8325114859430123e-04, -6.4407813263853044e-04, -6.3143874312577742e-05,
         1.8135960357381258e-03,  1.7808008939529483e-04, -5.6438848379409157e-04, -5.5311691945889258e-05,
         4.1802325329999405e-04,  4.1138714346210383e-05, -2.4134583782027904e-04, -2.3751447801897194e-05,
                                                         
         2.7426849265788276e-03,  7.6901868952394768e-05,  1.5834898806625327e-03,  4.4399314740898159e-05,
         7.0025973717428884e-03,  1.9640674915167951e-04,  8.7597171627717518e-04,  2.4596860011632871e-05,
         8.9364896785851626e-03,  2.5073307216640398e-04,  2.4056152766195870e-04,  6.7684572050008802e-06,
         8.8699248778119841e-03,  2.4891211971957313e-04, -2.7899273330690764e-04, -7.8197845903601858e-06,
         7.2812441275886316e-03,  2.0435956575817714e-04, -6.3823252549092191e-04, -1.7902644432337120e-05,
         4.8041572544601915e-03,  1.3488355963031703e-04, -7.9191424751587076e-04, -2.2209346407803159e-05,
         2.2306540994138199e-03,  6.2708495935571320e-05, -6.9389849181050321e-04, -1.9460946045136666e-05,
         5.1439332815130597e-04,  1.4500574311018977e-05, -2.9698512647750402e-04, -8.3719104818709799e-06,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    }
  } if (poly_order == 2) {
    if (bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.lo_type[1] == GKYL_POISSON_PERIODIC) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
	-9.5101835527468453e-02, -5.2237861001368492e-02, 4.3177375058357795e-02, 2.3387400381385207e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, -1.1937067064604241e-03, 5.2456792020110785e-03, 
	9.5101835527468453e-02, 5.2237861001368610e-02, 4.3177375058357802e-02, 2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905626e-03, -1.1937067064604172e-03, -5.2456792020110915e-03, 
	9.5101835527468384e-02, 5.2237861001368582e-02, -4.3177375058357843e-02, -2.3387400381385235e-02, -2.0675606649251488e-03, -9.0857828980905643e-03, 1.1937067064604181e-03, -5.2456792020110924e-03, 
	-9.5101835527468620e-02, -5.2237861001368589e-02, -4.3177375058357788e-02, -2.3387400381385214e-02, 2.0675606649251388e-03, 9.0857828980905712e-03, 1.1937067064604163e-03, 5.2456792020110898e-03, 
	-9.5101835527468509e-02, -5.2237861001368562e-02, 4.3177375058357823e-02, 2.3387400381385225e-02, 2.0675606649251596e-03, 9.0857828980905470e-03, -1.1937067064604137e-03, 5.2456792020110837e-03, 
	9.5101835527468453e-02, 5.2237861001368548e-02, 4.3177375058357816e-02, 2.3387400381385211e-02, -2.0675606649251488e-03, -9.0857828980905591e-03, -1.1937067064604259e-03, -5.2456792020110759e-03, 
	9.5101835527468537e-02, 5.2237861001368562e-02, -4.3177375058357795e-02, -2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905782e-03, 1.1937067064604241e-03, -5.2456792020110915e-03, 
	-9.5101835527468398e-02, -5.2237861001368506e-02, -4.3177375058357809e-02, -2.3387400381385193e-02, 2.0675606649251596e-03, 9.0857828980905504e-03, 1.1937067064604207e-03, 5.2456792020110776e-03, 
	-2.2959614113698959e-01, -2.1637630496127446e-02, 1.0423940445355719e-01, 9.6873784266194515e-03, 4.9915329982914147e-03, 2.1935020297347764e-02, -2.8818629202324655e-03, 2.1728314693314739e-03, 
	2.2959614113698990e-01, 2.1637630496127457e-02, 1.0423940445355717e-01, 9.6873784266194498e-03, -4.9915329982914251e-03, -2.1935020297347764e-02, -2.8818629202324387e-03, -2.1728314693314635e-03, 
	2.2959614113698973e-01, 2.1637630496127453e-02, -1.0423940445355724e-01, -9.6873784266194377e-03, -4.9915329982914043e-03, -2.1935020297347750e-02, 2.8818629202324517e-03, -2.1728314693314518e-03, 
	-2.2959614113698987e-01, -2.1637630496127373e-02, -1.0423940445355719e-01, -9.6873784266194654e-03, 4.9915329982914147e-03, 2.1935020297347754e-02, 2.8818629202324413e-03, 2.1728314693314496e-03, 
	-2.2959614113698981e-01, -2.1637630496127425e-02, 1.0423940445355721e-01, 9.6873784266194429e-03, 4.9915329982914147e-03, 2.1935020297347736e-02, -2.8818629202324378e-03, 2.1728314693314679e-03, 
	2.2959614113698976e-01, 2.1637630496127457e-02, 1.0423940445355720e-01, 9.6873784266194567e-03, -4.9915329982914113e-03, -2.1935020297347757e-02, -2.8818629202324551e-03, -2.1728314693314783e-03, 
	2.2959614113698987e-01, 2.1637630496127443e-02, -1.0423940445355716e-01, -9.6873784266194429e-03, -4.9915329982914113e-03, -2.1935020297347771e-02, 2.8818629202324551e-03, -2.1728314693314600e-03, 
	-2.2959614113698959e-01, -2.1637630496127443e-02, -1.0423940445355719e-01, -9.6873784266194619e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324560e-03, 2.1728314693314809e-03, 
	-2.2959614113698965e-01, 2.1637630496127432e-02, 1.0423940445355719e-01, -9.6873784266194515e-03, 4.9915329982914147e-03, 2.1935020297347771e-02, -2.8818629202324655e-03, -2.1728314693314657e-03, 
	2.2959614113698987e-01, -2.1637630496127495e-02, 1.0423940445355716e-01, -9.6873784266194689e-03, -4.9915329982914321e-03, -2.1935020297347757e-02, -2.8818629202324499e-03, 2.1728314693314679e-03, 
	2.2959614113698973e-01, -2.1637630496127470e-02, -1.0423940445355720e-01, 9.6873784266194654e-03, -4.9915329982914217e-03, -2.1935020297347750e-02, 2.8818629202324586e-03, 2.1728314693314540e-03, 
	-2.2959614113698976e-01, 2.1637630496127436e-02, -1.0423940445355717e-01, 9.6873784266194793e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324551e-03, -2.1728314693314583e-03, 
	-2.2959614113698976e-01, 2.1637630496127460e-02, 1.0423940445355719e-01, -9.6873784266194758e-03, 4.9915329982914286e-03, 2.1935020297347740e-02, -2.8818629202324586e-03, -2.1728314693314653e-03, 
	2.2959614113698981e-01, -2.1637630496127446e-02, 1.0423940445355721e-01, -9.6873784266194377e-03, -4.9915329982914251e-03, -2.1935020297347761e-02, -2.8818629202324577e-03, 2.1728314693314722e-03, 
	2.2959614113698987e-01, -2.1637630496127453e-02, -1.0423940445355717e-01, 9.6873784266194238e-03, -4.9915329982914321e-03, -2.1935020297347747e-02, 2.8818629202324586e-03, 2.1728314693314683e-03, 
	-2.2959614113698959e-01, 2.1637630496127436e-02, -1.0423940445355719e-01, 9.6873784266194619e-03, 4.9915329982914147e-03, 2.1935020297347750e-02, 2.8818629202324603e-03, -2.1728314693314774e-03, 
	-9.5101835527468509e-02, 5.2237861001368499e-02, 4.3177375058357795e-02, -2.3387400381385211e-02, 2.0675606649251735e-03, 9.0857828980905678e-03, -1.1937067064604276e-03, -5.2456792020110863e-03, 
	9.5101835527468342e-02, -5.2237861001368568e-02, 4.3177375058357781e-02, -2.3387400381385200e-02, -2.0675606649251449e-03, -9.0857828980905608e-03, -1.1937067064604207e-03, 5.2456792020110811e-03, 
	9.5101835527468342e-02, -5.2237861001368568e-02, -4.3177375058357788e-02, 2.3387400381385211e-02, -2.0675606649251509e-03, -9.0857828980905660e-03, 1.1937067064604163e-03, 5.2456792020110898e-03, 
	-9.5101835527468509e-02, 5.2237861001368534e-02, -4.3177375058357767e-02, 2.3387400381385211e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, 1.1937067064604215e-03, -5.2456792020110932e-03, 
	-9.5101835527468426e-02, 5.2237861001368534e-02, 4.3177375058357809e-02, -2.3387400381385193e-02, 2.0675606649251457e-03, 9.0857828980905556e-03, -1.1937067064604207e-03, -5.2456792020110776e-03, 
	9.5101835527468509e-02, -5.2237861001368527e-02, 4.3177375058357823e-02, -2.3387400381385228e-02, -2.0675606649251423e-03, -9.0857828980905574e-03, -1.1937067064604128e-03, 5.2456792020110846e-03, 
	9.5101835527468495e-02, -5.2237861001368555e-02, -4.3177375058357816e-02, 2.3387400381385221e-02, -2.0675606649251466e-03, -9.0857828980905504e-03, 1.1937067064604102e-03, 5.2456792020110854e-03, 
	-9.5101835527468453e-02, 5.2237861001368499e-02, -4.3177375058357823e-02, 2.3387400381385200e-02, 2.0675606649251527e-03, 9.0857828980905539e-03, 1.1937067064604337e-03, -5.2456792020110776e-03, 
	9.5101835527468231e-02, 5.2237861001368499e-02, -4.3177375058357809e-02, -2.3387400381385207e-02, -2.0675606649251488e-03, -9.0857828980905678e-03, 1.1937067064604233e-03, -5.2456792020110932e-03, 
	-9.5101835527468676e-02, -5.2237861001368610e-02, -4.3177375058357767e-02, -2.3387400381385211e-02, 2.0675606649251457e-03, 9.0857828980905591e-03, 1.1937067064604137e-03, 5.2456792020110846e-03, 
	-9.5101835527468565e-02, -5.2237861001368562e-02, 4.3177375058357836e-02, 2.3387400381385239e-02, 2.0675606649251319e-03, 9.0857828980905678e-03, -1.1937067064604172e-03, 5.2456792020110837e-03, 
	9.5101835527468426e-02, 5.2237861001368582e-02, 4.3177375058357823e-02, 2.3387400381385221e-02, -2.0675606649251509e-03, -9.0857828980905574e-03, -1.1937067064604154e-03, -5.2456792020110776e-03, 
	9.5101835527468509e-02, 5.2237861001368589e-02, -4.3177375058357753e-02, -2.3387400381385214e-02, -2.0675606649251423e-03, -9.0857828980905435e-03, 1.1937067064604198e-03, -5.2456792020110820e-03, 
	-9.5101835527468231e-02, -5.2237861001368478e-02, -4.3177375058357809e-02, -2.3387400381385228e-02, 2.0675606649251596e-03, 9.0857828980905383e-03, 1.1937067064604172e-03, 5.2456792020110707e-03, 
	-9.5101835527468370e-02, -5.2237861001368541e-02, 4.3177375058357739e-02, 2.3387400381385197e-02, 2.0675606649251527e-03, 9.0857828980905539e-03, -1.1937067064604207e-03, 5.2456792020110724e-03, 
	9.5101835527468315e-02, 5.2237861001368499e-02, 4.3177375058357753e-02, 2.3387400381385200e-02, -2.0675606649251509e-03, -9.0857828980905574e-03, -1.1937067064604215e-03, -5.2456792020110898e-03, 
	2.2959614113698942e-01, 2.1637630496127418e-02, -1.0423940445355720e-01, -9.6873784266194619e-03, -4.9915329982914113e-03, -2.1935020297347771e-02, 2.8818629202324586e-03, -2.1728314693314574e-03, 
	-2.2959614113699003e-01, -2.1637630496127470e-02, -1.0423940445355714e-01, -9.6873784266194619e-03, 4.9915329982914008e-03, 2.1935020297347740e-02, 2.8818629202324387e-03, 2.1728314693314609e-03, 
	-2.2959614113698998e-01, -2.1637630496127495e-02, 1.0423940445355724e-01, 9.6873784266194567e-03, 4.9915329982914008e-03, 2.1935020297347771e-02, -2.8818629202324447e-03, 2.1728314693314791e-03, 
	2.2959614113698981e-01, 2.1637630496127450e-02, 1.0423940445355726e-01, 9.6873784266194706e-03, -4.9915329982914286e-03, -2.1935020297347750e-02, -2.8818629202324638e-03, -2.1728314693314713e-03, 
	2.2959614113699001e-01, 2.1637630496127470e-02, -1.0423940445355717e-01, -9.6873784266194706e-03, -4.9915329982914355e-03, -2.1935020297347764e-02, 2.8818629202324621e-03, -2.1728314693314817e-03, 
	-2.2959614113698937e-01, -2.1637630496127439e-02, -1.0423940445355716e-01, -9.6873784266194151e-03, 4.9915329982914147e-03, 2.1935020297347722e-02, 2.8818629202324586e-03, 2.1728314693314783e-03, 
	-2.2959614113698959e-01, -2.1637630496127418e-02, 1.0423940445355706e-01, 9.6873784266194290e-03, 4.9915329982914286e-03, 2.1935020297347764e-02, -2.8818629202324586e-03, 2.1728314693314869e-03, 
	2.2959614113698953e-01, 2.1637630496127457e-02, 1.0423940445355712e-01, 9.6873784266194498e-03, -4.9915329982914182e-03, -2.1935020297347757e-02, -2.8818629202324560e-03, -2.1728314693314627e-03, 
	2.2959614113698945e-01, -2.1637630496127394e-02, -1.0423940445355723e-01, 9.6873784266194342e-03, -4.9915329982914147e-03, -2.1935020297347750e-02, 2.8818629202324517e-03, 2.1728314693314613e-03, 
	-2.2959614113699015e-01, 2.1637630496127450e-02, -1.0423940445355720e-01, 9.6873784266194619e-03, 4.9915329982914286e-03, 2.1935020297347750e-02, 2.8818629202324577e-03, -2.1728314693314653e-03, 
	-2.2959614113699009e-01, 2.1637630496127460e-02, 1.0423940445355727e-01, -9.6873784266194654e-03, 4.9915329982914425e-03, 2.1935020297347778e-02, -2.8818629202324517e-03, -2.1728314693314800e-03, 
	2.2959614113698970e-01, -2.1637630496127484e-02, 1.0423940445355724e-01, -9.6873784266194619e-03, -4.9915329982914008e-03, -2.1935020297347764e-02, -2.8818629202324517e-03, 2.1728314693314661e-03, 
	2.2959614113698992e-01, -2.1637630496127498e-02, -1.0423940445355714e-01, 9.6873784266194654e-03, -4.9915329982914182e-03, -2.1935020297347764e-02, 2.8818629202324447e-03, 2.1728314693314761e-03, 
	-2.2959614113698948e-01, 2.1637630496127401e-02, -1.0423940445355717e-01, 9.6873784266194342e-03, 4.9915329982914564e-03, 2.1935020297347750e-02, 2.8818629202324829e-03, -2.1728314693314635e-03, 
	-2.2959614113698953e-01, 2.1637630496127484e-02, 1.0423940445355713e-01, -9.6873784266194203e-03, 4.9915329982914425e-03, 2.1935020297347750e-02, -2.8818629202324794e-03, -2.1728314693314913e-03, 
	2.2959614113698965e-01, -2.1637630496127401e-02, 1.0423940445355709e-01, -9.6873784266194446e-03, -4.9915329982914251e-03, -2.1935020297347743e-02, -2.8818629202324517e-03, 2.1728314693314670e-03, 
	9.5101835527468259e-02, -5.2237861001368513e-02, -4.3177375058357836e-02, 2.3387400381385239e-02, -2.0675606649251405e-03, -9.0857828980905626e-03, 1.1937067064604181e-03, 5.2456792020110898e-03, 
	-9.5101835527468676e-02, 5.2237861001368631e-02, -4.3177375058357767e-02, 2.3387400381385214e-02, 2.0675606649251457e-03, 9.0857828980905608e-03, 1.1937067064604068e-03, -5.2456792020110776e-03, 
	-9.5101835527468620e-02, 5.2237861001368603e-02, 4.3177375058357795e-02, -2.3387400381385239e-02, 2.0675606649251319e-03, 9.0857828980905608e-03, -1.1937067064604102e-03, -5.2456792020110828e-03, 
	9.5101835527468356e-02, -5.2237861001368562e-02, 4.3177375058357823e-02, -2.3387400381385228e-02, -2.0675606649251492e-03, -9.0857828980905608e-03, -1.1937067064604163e-03, 5.2456792020110880e-03, 
	9.5101835527468453e-02, -5.2237861001368582e-02, -4.3177375058357760e-02, 2.3387400381385211e-02, -2.0675606649251466e-03, -9.0857828980905591e-03, 1.1937067064604172e-03, 5.2456792020110811e-03, 
	-9.5101835527468342e-02, 5.2237861001368485e-02, -4.3177375058357788e-02, 2.3387400381385207e-02, 2.0675606649251596e-03, 9.0857828980905574e-03, 1.1937067064604267e-03, -5.2456792020110863e-03, 
	-9.5101835527468287e-02, 5.2237861001368478e-02, 4.3177375058357788e-02, -2.3387400381385197e-02, 2.0675606649251735e-03, 9.0857828980905383e-03, -1.1937067064604241e-03, -5.2456792020110768e-03, 
	9.5101835527468398e-02, -5.2237861001368541e-02, 4.3177375058357746e-02, -2.3387400381385211e-02, -2.0675606649251379e-03, -9.0857828980905539e-03, -1.1937067064604189e-03, 5.2456792020110811e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) { {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         2.6174232732491770e-01,  1.4239797833436429e-01,  1.4239797833437121e-01,  7.7179579805169304e-02,
        -6.7537275729811939e-03, -6.7537275729829148e-03, -3.8992664322967511e-03, -3.8992664322970148e-03,
         6.6919163252610048e-01,  3.6607316345625468e-01,  9.5063998582039566e-02,  5.3241649763932503e-02,
        -1.5712542948565247e-02, -5.0333399450138384e-03, -1.2731080364168593e-03, -2.9060001721712362e-03,
         9.3261285108184055e-01,  5.1325052577226338e-01,  5.7325647937302174e-02,  3.1906370575560429e-02,
        -1.9514998592151618e-02, -4.7983776055374530e-03, -9.2224075298933392e-04, -2.7703446022362605e-03,
         1.0655299621632952e+00,  5.8724826482558512e-01,  1.9276481795971818e-02,  1.0736800514483670e-02,
        -2.1638935540250832e-02, -4.9049622075441517e-03, -3.0401481573749179e-04, -2.8318812508971261e-03,
         1.0655299621633028e+00,  5.8724826482558323e-01, -1.9276481795972602e-02, -1.0736800514480182e-02,
        -2.1638935540251626e-02, -4.9049622075482769e-03,  3.0401481573705659e-04, -2.8318812508936436e-03,
         9.3261285108186132e-01,  5.1325052577225383e-01, -5.7325647937298004e-02, -3.1906370575565557e-02,
        -1.9514998592151112e-02, -4.7983776055450077e-03,  9.2224075299059171e-04, -2.7703446022306617e-03,
         6.6919163252609015e-01,  3.6607316345626151e-01, -9.5063998582045506e-02, -5.3241649763929401e-02,
        -1.5712542948564903e-02, -5.0333399450088745e-03,  1.2731080364155342e-03, -2.9060001721745977e-03,
         2.6174232732491648e-01,  1.4239797833437010e-01, -1.4239797833436521e-01, -7.7179579805168486e-02,
        -6.7537275729838489e-03, -6.7537275729819659e-03,  3.8992664322963539e-03, -3.8992664322969480e-03,
                                                         
         6.6919163252608815e-01,  9.5063998582044021e-02,  3.6607316345626351e-01,  5.3241649763928853e-02,
        -5.0333399450085128e-03, -1.5712542948565542e-02, -2.9060001721746753e-03, -1.2731080364157741e-03,
         1.7423637960024951e+00,  2.5663155896154988e-01,  2.5663155896154877e-01,  3.9269506879060931e-02,
        -1.3304720861538551e-02, -1.3304720861538329e-02, -1.8694838265539103e-03, -1.8694838265541740e-03,
         2.4546088336295506e+00,  3.6679544630615912e-01,  1.5482377206844999e-01,  2.4122215546182844e-02,
        -1.8491615532474651e-02, -1.3118456058247588e-02, -1.1251712079693174e-03, -2.0332549321147595e-03,
         2.8136363432486524e+00,  4.2275290286990563e-01,  5.2063204566974094e-02,  8.1141391202813138e-03,
        -2.1072335114916858e-02, -1.3426476507161618e-02, -3.6480793765665813e-04, -2.0880173238903453e-03,
         2.8136363432486524e+00,  4.2275290286990524e-01, -5.2063204566972256e-02, -8.1141391202823199e-03,
        -2.1072335114916522e-02, -1.3426476507160893e-02,  3.6480793765687319e-04, -2.0880173238910131e-03,
         2.4546088336295533e+00,  3.6679544630615807e-01, -1.5482377206845077e-01, -2.4122215546182161e-02,
        -1.8491615532474186e-02, -1.3118456058246705e-02,  1.1251712079691587e-03, -2.0332549321155132e-03,
         1.7423637960024967e+00,  2.5663155896154877e-01, -2.5663155896155010e-01, -3.9269506879060175e-02,
        -1.3304720861539127e-02, -1.3304720861538711e-02,  1.8694838265534816e-03, -1.8694838265539018e-03,
         6.6919163252608682e-01,  9.5063998582043396e-02, -3.6607316345626262e-01, -5.3241649763929991e-02,
        -5.0333399450070591e-03, -1.5712542948565018e-02,  2.9060001721762860e-03, -1.2731080364161030e-03,
                                                         
         9.3261285108185710e-01,  5.7325647937300842e-02,  5.1325052577225627e-01,  3.1906370575564183e-02,
        -4.7983776055426555e-03, -1.9514998592151105e-02, -2.7703446022324602e-03, -9.2224075298994704e-04,
         2.4546088336295533e+00,  1.5482377206845188e-01,  3.6679544630615785e-01,  2.4122215546181956e-02,
        -1.3118456058246891e-02, -1.8491615532474508e-02, -2.0332549321152001e-03, -1.1251712079690640e-03,
         3.4763948230693575e+00,  2.2273928406436580e-01,  2.2273928406436491e-01,  1.4934123050414328e-02,
        -1.8796563154148561e-02, -1.8796563154148505e-02, -1.2450017281912052e-03, -1.2450017281912329e-03,
         3.9932287149419077e+00,  2.5750472735140589e-01,  7.4989196262071997e-02,  5.0536521254172630e-03,
        -2.1676786288648296e-02, -1.9312217015521736e-02, -4.1789587383843661e-04, -1.3101165429916409e-03,
         3.9932287149419086e+00,  2.5750472735140550e-01, -7.4989196262071983e-02, -5.0536521254171546e-03,
        -2.1676786288648310e-02, -1.9312217015521806e-02,  4.1789587383839465e-04, -1.3101165429915658e-03,
         3.4763948230693589e+00,  2.2273928406436544e-01, -2.2273928406436508e-01, -1.4934123050414191e-02,
        -1.8796563154148464e-02, -1.8796563154148630e-02,  1.2450017281912470e-03, -1.2450017281910922e-03,
         2.4546088336295528e+00,  1.5482377206845233e-01, -3.6679544630615840e-01, -2.4122215546181762e-02,
        -1.3118456058247220e-02, -1.8491615532474616e-02,  2.0332549321149230e-03, -1.1251712079691938e-03,
         9.3261285108185132e-01,  5.7325647937300522e-02, -5.1325052577225905e-01, -3.1906370575564613e-02,
        -4.7983776055399753e-03, -1.9514998592151157e-02,  2.7703446022344595e-03, -9.2224075299001112e-04,
                                                         
         1.0655299621633085e+00,  1.9276481795971488e-02,  5.8724826482558201e-01,  1.0736800514482708e-02,
        -4.9049622075492960e-03, -2.1638935540251314e-02, -2.8318812508928769e-03, -3.0401481573748951e-04,
         2.8136363432486569e+00,  5.2063204566973623e-02,  4.2275290286990502e-01,  8.1141391202819799e-03,
        -1.3426476507160820e-02, -2.1072335114916560e-02, -2.0880173238911450e-03, -3.6480793765673593e-04,
         3.9932287149419117e+00,  7.4989196262073329e-02,  2.5750472735140501e-01,  5.0536521254175961e-03,
        -1.9312217015521948e-02, -2.1676786288648508e-02, -1.3101165429914947e-03, -4.1789587383842935e-04,
         4.5910338085415718e+00,  8.6771187773204306e-02,  8.6771187773202557e-02,  1.7165201446689979e-03,
        -2.2348374449956559e-02, -2.2348374449956337e-02, -4.4280976908144193e-04, -4.4280976908130315e-04,
         4.5910338085415718e+00,  8.6771187773204195e-02, -8.6771187773203362e-02, -1.7165201446691534e-03,
        -2.2348374449956399e-02, -2.2348374449956399e-02,  4.4280976908152368e-04, -4.4280976908143911e-04,
         3.9932287149419095e+00,  7.4989196262072677e-02, -2.5750472735140451e-01, -5.0536521254171502e-03,
        -1.9312217015521778e-02, -2.1676786288648116e-02,  1.3101165429914377e-03, -4.1789587383835627e-04,
         2.8136363432486573e+00,  5.2063204566974289e-02, -4.2275290286990530e-01, -8.1141391202814977e-03,
        -1.3426476507160697e-02, -2.1072335114916685e-02,  2.0880173238911949e-03, -3.6480793765670535e-04,
         1.0655299621633088e+00,  1.9276481795971901e-02, -5.8724826482558179e-01, -1.0736800514483127e-02,
        -4.9049622075491390e-03, -2.1638935540251369e-02,  2.8318812508928517e-03, -3.0401481573744788e-04,
                                                         
         1.0655299621633059e+00, -1.9276481795971124e-02,  5.8724826482558301e-01, -1.0736800514483365e-02,
        -4.9049622075478649e-03, -2.1638935540251245e-02, -2.8318812508938908e-03,  3.0401481573757262e-04,
         2.8136363432486569e+00, -5.2063204566973290e-02,  4.2275290286990536e-01, -8.1141391202814873e-03,
        -1.3426476507160872e-02, -2.1072335114916612e-02, -2.0880173238910075e-03,  3.6480793765668437e-04,
         3.9932287149419130e+00, -7.4989196262072108e-02,  2.5750472735140517e-01, -5.0536521254173029e-03,
        -1.9312217015521708e-02, -2.1676786288648268e-02, -1.3101165429914407e-03,  4.1789587383857528e-04,
         4.5910338085415718e+00, -8.6771187773202793e-02,  8.6771187773202904e-02, -1.7165201446689685e-03,
        -2.2348374449956222e-02, -2.2348374449956222e-02, -4.4280976908148416e-04,  4.4280976908131210e-04,
         4.5910338085415727e+00, -8.6771187773202474e-02, -8.6771187773202779e-02,  1.7165201446694099e-03,
        -2.2348374449956278e-02, -2.2348374449956500e-02,  4.4280976908150861e-04,  4.4280976908139759e-04,
         3.9932287149419117e+00, -7.4989196262071275e-02, -2.5750472735140512e-01,  5.0536521254170332e-03,
        -1.9312217015521788e-02, -2.1676786288648237e-02,  1.3101165429913700e-03,  4.1789587383833242e-04,
         2.8136363432486591e+00, -5.2063204566973408e-02, -4.2275290286990502e-01,  8.1141391202817024e-03,
        -1.3426476507160754e-02, -2.1072335114916716e-02,  2.0880173238911632e-03,  3.6480793765672633e-04,
         1.0655299621633081e+00, -1.9276481795971284e-02, -5.8724826482558301e-01,  1.0736800514483185e-02,
        -4.9049622075484347e-03, -2.1638935540251134e-02,  2.8318812508933353e-03,  3.0401481573757126e-04,
                                                         
         9.3261285108185266e-01, -5.7325647937301265e-02,  5.1325052577225927e-01, -3.1906370575562747e-02,
        -4.7983776055406909e-03, -1.9514998592151191e-02, -2.7703446022337964e-03,  9.2224075298978842e-04,
         2.4546088336295555e+00, -1.5482377206845083e-01,  3.6679544630615885e-01, -2.4122215546182133e-02,
        -1.3118456058247074e-02, -1.8491615532474526e-02, -2.0332549321151190e-03,  1.1251712079692335e-03,
         3.4763948230693624e+00, -2.2273928406436483e-01,  2.2273928406436511e-01, -1.4934123050414449e-02,
        -1.8796563154148675e-02, -1.8796563154148564e-02, -1.2450017281912206e-03,  1.2450017281911187e-03,
         3.9932287149419121e+00, -2.5750472735140495e-01,  7.4989196262071955e-02, -5.0536521254169655e-03,
        -2.1676786288648216e-02, -1.9312217015521767e-02, -4.1789587383831805e-04,  1.3101165429915314e-03,
         3.9932287149419126e+00, -2.5750472735140467e-01, -7.4989196262072025e-02,  5.0536521254170704e-03,
        -2.1676786288648092e-02, -1.9312217015522087e-02,  4.1789587383837931e-04,  1.3101165429914169e-03,
         3.4763948230693611e+00, -2.2273928406436524e-01, -2.2273928406436566e-01,  1.4934123050414329e-02,
        -1.8796563154148464e-02, -1.8796563154148713e-02,  1.2450017281911525e-03,  1.2450017281911802e-03,
         2.4546088336295551e+00, -1.5482377206845158e-01, -3.6679544630615818e-01,  2.4122215546182102e-02,
        -1.3118456058246773e-02, -1.8491615532474696e-02,  2.0332549321152413e-03,  1.1251712079691884e-03,
         9.3261285108185532e-01, -5.7325647937301619e-02, -5.1325052577225805e-01,  3.1906370575563073e-02,
        -4.7983776055419217e-03, -1.9514998592151160e-02,  2.7703446022327573e-03,  9.2224075298971849e-04,
                                                         
         6.6919163252609370e-01, -9.5063998582043147e-02,  3.6607316345626068e-01, -5.3241649763929991e-02,
        -5.0333399450107766e-03, -1.5712542948565014e-02, -2.9060001721731752e-03,  1.2731080364163324e-03,
         1.7423637960024991e+00, -2.5663155896154932e-01,  2.5663155896154954e-01, -3.9269506879060966e-02,
        -1.3304720861538491e-02, -1.3304720861538602e-02, -1.8694838265540792e-03,  1.8694838265539704e-03,
         2.4546088336295564e+00, -3.6679544630615801e-01,  1.5482377206845152e-01, -2.4122215546181613e-02,
        -1.8491615532474564e-02, -1.3118456058247055e-02, -1.1251712079691414e-03,  2.0332549321152114e-03,
         2.8136363432486586e+00, -4.2275290286990402e-01,  5.2063204566973983e-02, -8.1141391202815896e-03,
        -2.1072335114916560e-02, -1.3426476507160459e-02, -3.6480793765670513e-04,  2.0880173238911949e-03,
         2.8136363432486595e+00, -4.2275290286990463e-01, -5.2063204566974060e-02,  8.1141391202815445e-03,
        -2.1072335114916383e-02, -1.3426476507160671e-02,  3.6480793765677414e-04,  2.0880173238913636e-03,
         2.4546088336295546e+00, -3.6679544630615818e-01, -1.5482377206845166e-01,  2.4122215546182019e-02,
        -1.8491615532474439e-02, -1.3118456058246849e-02,  1.1251712079690349e-03,  2.0332549321153233e-03,
         1.7423637960024978e+00, -2.5663155896154966e-01, -2.5663155896154949e-01,  3.9269506879060584e-02,
        -1.3304720861538579e-02, -1.3304720861538607e-02,  1.8694838265540319e-03,  1.8694838265540597e-03,
         6.6919163252609004e-01, -9.5063998582043216e-02, -3.6607316345626228e-01,  5.3241649763930428e-02,
        -5.0333399450094695e-03, -1.5712542948565080e-02,  2.9060001721740326e-03,  1.2731080364163240e-03,
                                                         
         2.6174232732492186e-01, -1.4239797833436652e-01,  1.4239797833436699e-01, -7.7179579805169304e-02,
        -6.7537275729839755e-03, -6.7537275729838663e-03, -3.8992664322955498e-03,  3.8992664322956201e-03,
         6.6919163252609259e-01, -3.6607316345626134e-01,  9.5063998582043535e-02, -5.3241649763929630e-02,
        -1.5712542948564886e-02, -5.0333399450095380e-03, -1.2731080364162786e-03,  2.9060001721740491e-03,
         9.3261285108185321e-01, -5.1325052577226005e-01,  5.7325647937300730e-02, -3.1906370575564294e-02,
        -1.9514998592151306e-02, -4.7983776055403804e-03, -9.2224075298998217e-04,  2.7703446022338254e-03,
         1.0655299621633147e+00, -5.8724826482558012e-01,  1.9276481795972494e-02, -1.0736800514482444e-02,
        -2.1638935540251224e-02, -4.9049622075506205e-03, -3.0401481573732217e-04,  2.8318812508918382e-03,
         1.0655299621633119e+00, -5.8724826482558135e-01, -1.9276481795972668e-02,  1.0736800514482532e-02,
        -2.1638935540251137e-02, -4.9049622075495979e-03,  3.0401481573736000e-04,  2.8318812508924367e-03,
         9.3261285108185554e-01, -5.1325052577225794e-01, -5.7325647937301508e-02,  3.1906370575563600e-02,
        -1.9514998592151313e-02, -4.7983776055423831e-03,  9.2224075298983992e-04,  2.7703446022324342e-03,
         6.6919163252608904e-01, -3.6607316345626240e-01, -9.5063998582043160e-02,  5.3241649763930317e-02,
        -1.5712542948565038e-02, -5.0333399450093533e-03,  1.2731080364162960e-03,  2.9060001721740885e-03,
         2.6174232732492109e-01, -1.4239797833436624e-01, -1.4239797833436613e-01,  7.7179579805169332e-02,
        -6.7537275729838767e-03, -6.7537275729839148e-03,  3.8992664322956569e-03,  3.8992664322956318e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         3.4754908174043031e-01,  1.9147096450810874e-01, -2.8048110764785784e-02, -1.3769366894035592e-02,
        -7.1159030934106020e-03,  2.2827046419852215e-03,  1.8777907284159622e-03,  1.3179201395296886e-03,
         2.6018094650847862e-01,  1.5014541606507584e-01, -2.4623528318145276e-02, -1.1377222883164622e-02,
        -5.4317382757281126e-05,  5.5564633406491573e-04,  2.1992176825353575e-03,  3.2080256054490733e-04,
         2.2061875008365231e-01,  1.3141643256470523e-01,  9.6631517600851008e-03,  5.1140605928803942e-03,
         3.1310264886016711e-03,  6.6601405906152981e-03, -3.6015854094473347e-04,  3.8452339628319973e-03,
         3.1728477488866691e-01,  1.7738358224001630e-01,  3.1338192874284236e-02,  1.2875214739491091e-02,
        -4.4933329941175878e-03, -4.8107149369187568e-03, -4.0417674588017032e-03, -2.7774675638256754e-03,
         2.9856868259017194e-01,  1.6904725884004743e-01, -4.4358254629947795e-02, -1.8966632025180653e-02,
        -2.5805304029333908e-03, -6.5259204756552055e-03,  5.1461245497285120e-03, -3.7677419433303117e-03,
         1.6472836762602577e-01,  1.0615968567242218e-01, -6.4622802488493898e-03, -2.0693364430828479e-03,
         8.5621733627070755e-03,  1.3963870956120502e-02,  1.2871184688642739e-03,  8.0620446554442865e-03,
         3.2567896066607988e-01,  1.8138329229246208e-01,  7.9694807939026877e-02,  3.4130327869920304e-02,
        -5.1491649123717653e-03, -1.2897504262264271e-03, -9.2033633129310552e-03, -7.4463775577016827e-04,
         4.5758071423995317e-01,  2.4387817433520742e-01, -1.7203978611667874e-02, -5.9370449568280142e-03,
        -1.5729094905946293e-02, -1.1872819639637478e-02,  3.0950378831334125e-03, -6.8547756149879950e-03,
                                                         
         9.2555564511589195e-01,  1.4283166346900319e-01, -5.7087203270645710e-02, -4.2107167636395056e-03,
        -6.6585943733794726e-03,  4.3810045656668519e-03,  9.3715478409904616e-04, -1.0646611375156228e-04,
         7.5643742119286139e-01,  1.3194920980432359e-01, -4.3151229134044317e-02, -9.9766889501246210e-04,
        -3.4774021071348250e-03,  2.3686176801918867e-03,  8.9950742716126829e-04,  7.2591693417456003e-04,
         7.0097359003902937e-01,  1.4059912442543343e-01,  2.5969251836348434e-02,  5.4595811792827933e-03,
        -9.8784417722106539e-04,  1.3863694644921702e-02,  5.3783951383761853e-04,  3.1373990954381765e-04,
         8.9395264296783827e-01,  1.6114696623034061e-01,  5.3336355793392402e-02,  4.9642255681518815e-03,
        -1.6240840546830664e-04, -1.1009301949378427e-02, -6.1273948817407056e-05, -8.0128831641327464e-04,
         8.5858667148480283e-01,  1.5815426669652310e-01, -7.7423345159838022e-02, -6.2531566299462403e-03,
         4.2102309020289271e-04, -1.3850860419405031e-02,  3.9811827989687796e-04, -4.6131410499152816e-04,
         6.0862841609240770e-01,  1.3886189996819126e-01, -1.1804099479237351e-02, -3.6257448687969356e-03,
        -1.6305922506315459e-04,  2.8818588277168135e-02, -7.3533836184461665e-04,  5.1433038859834352e-04,
         8.8573073934651869e-01,  1.4220407629839985e-01,  1.3180886462596331e-01,  5.2114498588328763e-03,
        -4.9622366030997230e-03, -2.1500466509340923e-03, -2.0354679892535497e-03,  2.4794549885235803e-04,
         1.1259773284051051e+00,  1.5150225754918686e-01, -2.1648595211938813e-02, -5.4796944887231718e-04,
        -8.3847823258310743e-03, -2.4261809779577561e-02,  5.9460294920793655e-05, -2.9801117729407584e-04,
                                                         
         1.3294757533414434e+00,  8.9788893137442261e-02, -6.2911497442604103e-02, -3.9614170768240779e-05,
        -7.1100499232079391e-03,  4.6679346543523435e-03,  2.4956613702374012e-04,  2.7212527769270804e-04,
         1.1558500674443988e+00,  9.5187518915131650e-02, -3.6801896226935131e-02,  3.8875289024224544e-03,
        -6.1608141969413492e-03,  5.0780393003260101e-03,  2.9847536506069259e-04,  8.3836836755803482e-04,
         1.1608685714961731e+00,  1.2037284313135879e-01,  5.3537968971594352e-02,  1.0607159984611206e-02,
        -4.5111155521692318e-03,  1.5797402263497122e-02,  6.5397859158025575e-04,  8.0268670457471416e-04,
         1.4515046707409818e+00,  1.5802986276729328e-01,  7.3764923340376418e-02,  7.7325221486257015e-03,
        -2.2739467835626078e-03, -1.5570616292196532e-02,  6.3765139919740651e-04, -1.8321877472712517e-03,
         1.4115464704823018e+00,  1.5776040893168652e-01, -9.9620570722380883e-02, -7.8173884689785901e-03,
        -2.1636612802878562e-03, -1.7728464227987777e-02, -5.7397803419402017e-04, -1.7774214977044778e-03,
         1.0684562036485115e+00,  1.2101078886542831e-01, -3.5201132471140130e-02, -9.9394198301684059e-03,
        -4.5077782114690034e-03,  3.1273725348626311e-02, -7.7939850703535667e-04,  9.0314366050543921e-04,
         1.3060178631336954e+00,  9.8310629669391275e-02,  1.3131146816260594e-01, -3.4378558002396307e-03,
        -6.6184622280532198e-03, -5.2033263638921626e-04, -4.3920547811377179e-04,  6.9297032614724626e-04,
         1.5457214754324096e+00,  9.2030059599753342e-02, -2.4079263611516431e-02, -9.9293276550431963e-04,
        -7.4607497917453056e-03, -2.4657338466129757e-02, -4.7089473518925407e-05,  6.9652583640923910e-05,
                                                         
         1.5410607675732244e+00,  3.1259713850908775e-02, -6.1609756908925520e-02,  4.0278881947996382e-04,
        -7.9698971022144766e-03,  5.4767589025950247e-03, -5.1276054205655116e-05,  1.9484961975728570e-04,
         1.3895291108288024e+00,  3.6049179432957129e-02, -2.3292866819651249e-02,  2.8177193660516302e-03,
        -9.0097499591924313e-03,  7.4784863973146746e-03, -5.4908327268816275e-04,  5.4753041006383430e-04,
         1.4814229261146523e+00,  5.3722083116000623e-02,  9.0436045201508081e-02,  7.5759200922297047e-03,
        -1.3014031031221703e-02,  1.8391205195993725e-02, -1.7627894821588457e-03,  6.9484611672705331e-04,
         1.9036696745231794e+00,  8.2014576705607065e-02,  1.0223309369326436e-01,  6.0293706608141534e-03,
        -1.8550660598603910e-02, -2.1202408094434324e-02, -1.4337850889725194e-03, -1.4193287657706831e-03,
         1.8635410509569001e+00,  8.2038435164303550e-02, -1.2804640133826681e-01, -6.0050462696767748e-03,
        -1.8554396836116659e-02, -2.3251244190818931e-02,  1.4316279712387684e-03, -1.4111570005111908e-03,
         1.3908819207316627e+00,  5.4053805794842834e-02, -7.0885965719392860e-02, -7.4371696780804368e-03,
        -1.3096490503921039e-02,  3.4035995955905070e-02,  1.7194957188661293e-03,  6.9165401818151314e-04,
         1.5458079707862993e+00,  3.6785251742612025e-02,  1.1844442680216077e-01, -2.7477367760387942e-03,
        -9.2110386667161315e-03,  1.5895856233272834e-03,  5.2377094526747942e-04,  5.2519154906815752e-04,
         1.7614161838416025e+00,  3.1685800667273945e-02, -2.7278574910696574e-02, -6.3584621477922470e-04,
        -8.0924625745515209e-03, -2.4216573332765313e-02,  1.2203926265282081e-04,  1.8482328475642772e-04,
                                                         
         1.5410607675732240e+00, -3.1259713850908713e-02, -6.1609756908925596e-02, -4.0278881947999954e-04,
        -7.9698971022143933e-03,  5.4767589025950534e-03, -5.1276054205712274e-05, -1.9484961975725637e-04,
         1.3895291108288028e+00, -3.6049179432956852e-02, -2.3292866819651072e-02, -2.8177193660514754e-03,
        -9.0097499591924035e-03,  7.4784863973146321e-03, -5.4908327268812123e-04, -5.4753041006383484e-04,
         1.4814229261146523e+00, -5.3722083116000553e-02,  9.0436045201507997e-02, -7.5759200922298131e-03,
        -1.3014031031221738e-02,  1.8391205195993847e-02, -1.7627894821589177e-03, -6.9484611672699509e-04,
         1.9036696745231787e+00, -8.2014576705607023e-02,  1.0223309369326418e-01, -6.0293706608139782e-03,
        -1.8550660598603833e-02, -2.1202408094434275e-02, -1.4337850889724136e-03,  1.4193287657707004e-03,
         1.8635410509568999e+00, -8.2038435164303550e-02, -1.2804640133826664e-01,  6.0050462696767818e-03,
        -1.8554396836116590e-02, -2.3251244190819000e-02,  1.4316279712386628e-03,  1.4111570005111837e-03,
         1.3908819207316627e+00, -5.4053805794842855e-02, -7.0885965719392902e-02,  7.4371696780805912e-03,
        -1.3096490503921063e-02,  3.4035995955905098e-02,  1.7194957188662130e-03, -6.9165401818148885e-04,
         1.5458079707862988e+00, -3.6785251742611901e-02,  1.1844442680216073e-01,  2.7477367760388098e-03,
        -9.2110386667160222e-03,  1.5895856233272539e-03,  5.2377094526745579e-04, -5.2519154906813107e-04,
         1.7614161838416023e+00, -3.1685800667273778e-02, -2.7278574910696542e-02,  6.3584621477932846e-04,
        -8.0924625745514220e-03, -2.4216573332765355e-02,  1.2203926265283721e-04, -1.8482328475641070e-04,
                                                         
         1.3294757533414430e+00, -8.9788893137442191e-02, -6.2911497442604464e-02,  3.9614170768292408e-05,
        -7.1100499232079313e-03,  4.6679346543524623e-03,  2.4956613702379482e-04, -2.7212527769268630e-04,
         1.1558500674443986e+00, -9.5187518915131830e-02, -3.6801896226934798e-02, -3.8875289024225706e-03,
        -6.1608141969413084e-03,  5.0780393003259120e-03,  2.9847536506066310e-04, -8.3836836755803027e-04,
         1.1608685714961731e+00, -1.2037284313135867e-01,  5.3537968971594171e-02, -1.0607159984610976e-02,
        -4.5111155521692517e-03,  1.5797402263497157e-02,  6.5397859158026453e-04, -8.0268670457479710e-04,
         1.4515046707409809e+00, -1.5802986276729319e-01,  7.3764923340376501e-02, -7.7325221486257874e-03,
        -2.2739467835625805e-03, -1.5570616292196462e-02,  6.3765139919740119e-04,  1.8321877472712586e-03,
         1.4115464704823009e+00, -1.5776040893168661e-01, -9.9620570722381063e-02,  7.8173884689786872e-03,
        -2.1636612802878735e-03, -1.7728464227987822e-02, -5.7397803419405682e-04,  1.7774214977044585e-03,
         1.0684562036485103e+00, -1.2101078886542833e-01, -3.5201132471139970e-02,  9.9394198301683799e-03,
        -4.5077782114689643e-03,  3.1273725348626304e-02, -7.7939850703529119e-04, -9.0314366050547369e-04,
         1.3060178631336949e+00, -9.8310629669391220e-02,  1.3131146816260625e-01,  3.4378558002397842e-03,
        -6.6184622280531843e-03, -5.2033263638922255e-04, -4.3920547811384877e-04, -6.9297032614728052e-04,
         1.5457214754324098e+00, -9.2030059599753106e-02, -2.4079263611516445e-02,  9.9293276550434868e-04,
        -7.4607497917453776e-03, -2.4657338466129785e-02, -4.7089473518925374e-05, -6.9652583640937327e-05,
                                                         
         9.2555564511589228e-01, -1.4283166346900242e-01, -5.7087203270645939e-02,  4.2107167636395108e-03,
        -6.6585943733794293e-03,  4.3810045656671034e-03,  9.3715478409907424e-04,  1.0646611375162348e-04,
         7.5643742119286161e-01, -1.3194920980432301e-01, -4.3151229134044268e-02,  9.9766889501250069e-04,
        -3.4774021071347677e-03,  2.3686176801920480e-03,  8.9950742716123533e-04, -7.2591693417441388e-04,
         7.0097359003902948e-01, -1.4059912442543290e-01,  2.5969251836348382e-02, -5.4595811792828132e-03,
        -9.8784417722099275e-04,  1.3863694644921844e-02,  5.3783951383765062e-04, -3.1373990954366011e-04,
         8.9395264296783794e-01, -1.6114696623034006e-01,  5.3336355793392118e-02, -4.9642255681519613e-03,
        -1.6240840546821416e-04, -1.1009301949378315e-02, -6.1273948817443689e-05,  8.0128831641332842e-04,
         8.5858667148480194e-01, -1.5815426669652288e-01, -7.7423345159838203e-02,  6.2531566299462507e-03,
         4.2102309020295776e-04, -1.3850860419404911e-02,  3.9811827989688620e-04,  4.6131410499162032e-04,
         6.0862841609240648e-01, -1.3886189996819115e-01, -1.1804099479237122e-02,  3.6257448687969642e-03,
        -1.6305922506312725e-04,  2.8818588277168312e-02, -7.3533836184465481e-04, -5.1433038859821081e-04,
         8.8573073934651858e-01, -1.4220407629839935e-01,  1.3180886462596375e-01, -5.2114498588326360e-03,
        -4.9622366030995808e-03, -2.1500466509339223e-03, -2.0354679892534604e-03, -2.4794549885219248e-04,
         1.1259773284051056e+00, -1.5150225754918639e-01, -2.1648595211938699e-02,  5.4796944887228096e-04,
        -8.3847823258309216e-03, -2.4261809779577131e-02,  5.9460294920692505e-05,  2.9801117729436375e-04,
                                                         
         3.4754908174043120e-01, -1.9147096450810924e-01, -2.8048110764785136e-02,  1.3769366894036281e-02,
        -7.1159030934106497e-03,  2.2827046419849114e-03,  1.8777907284160782e-03, -1.3179201395300791e-03,
         2.6018094650848222e-01, -1.5014541606507417e-01, -2.4623528318145321e-02,  1.1377222883164596e-02,
        -5.4317382757140756e-05,  5.5564633406385939e-04,  2.1992176825353384e-03, -3.2080256054575182e-04,
         2.2061875008365400e-01, -1.3141643256470487e-01,  9.6631517600846636e-03, -5.1140605928807351e-03,
         3.1310264886016273e-03,  6.6601405906146901e-03, -3.6015854094482351e-04, -3.8452339628325693e-03,
         3.1728477488866830e-01, -1.7738358224001571e-01,  3.1338192874284528e-02, -1.2875214739490477e-02,
        -4.4933329941175626e-03, -4.8107149369193128e-03, -4.0417674588015843e-03,  2.7774675638252209e-03,
         2.9856868259017150e-01, -1.6904725884004720e-01, -4.4358254629948746e-02,  1.8966632025180070e-02,
        -2.5805304029333561e-03, -6.5259204756554336e-03,  5.1461245497283941e-03,  3.7677419433300198e-03,
         1.6472836762602511e-01, -1.0615968567242212e-01, -6.4622802488490489e-03,  2.0693364430829689e-03,
         8.5621733627069246e-03,  1.3963870956119893e-02,  1.2871184688642873e-03, -8.0620446554448624e-03,
         3.2567896066608065e-01, -1.8138329229246203e-01,  7.9694807939027446e-02, -3.4130327869920409e-02,
        -5.1491649123716690e-03, -1.2897504262269553e-03, -9.2033633129309355e-03,  7.4463775576960654e-04,
         4.5758071423995794e-01, -2.4387817433520531e-01, -1.7203978611668343e-02,  5.9370449568276907e-03,
        -1.5729094905946251e-02, -1.1872819639639669e-02,  3.0950378831332307e-03,  6.8547756149861657e-03,
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_PERIODIC && bcs.up_type[0] == GKYL_POISSON_PERIODIC) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         3.4754908174042964e-01, -2.8048110764784386e-02,  1.9147096450811057e-01, -1.3769366894036857e-02,  2.2827046419850081e-03, -7.1159030934107581e-03,  1.3179201395300958e-03,  1.8777907284162560e-03,
         9.2555564511589383e-01, -5.7087203270645925e-02,  1.4283166346900303e-01, -4.2107167636395723e-03,  4.3810045656670115e-03, -6.6585943733795568e-03, -1.0646611375170161e-04,  9.3715478409909473e-04,
         1.3294757533414452e+00, -6.2911497442604436e-02,  8.9788893137442344e-02, -3.9614170768272533e-05,  4.6679346543522706e-03, -7.1100499232081170e-03,  2.7212527769271232e-04,  2.4956613702375823e-04,
         1.5410607675732264e+00, -6.1609756908925783e-02,  3.1259713850908692e-02,  4.0278881947991574e-04,  5.4767589025949666e-03, -7.9698971022144141e-03,  1.9484961975729643e-04, -5.1276054205724837e-05,
         1.5410607675732255e+00, -6.1609756908925964e-02, -3.1259713850908942e-02, -4.0278881947993422e-04,  5.4767589025949874e-03, -7.9698971022144575e-03, -1.9484961975728526e-04, -5.1276054205658064e-05,
         1.3294757533414445e+00, -6.2911497442604505e-02, -8.9788893137442580e-02,  3.9614170768346550e-05,  4.6679346543523348e-03, -7.1100499232081196e-03, -2.7212527769269334e-04,  2.4956613702381553e-04,
         9.2555564511589206e-01, -5.7087203270645800e-02, -1.4283166346900306e-01,  4.2107167636395082e-03,  4.3810045656669195e-03, -6.6585943733794492e-03,  1.0646611375155487e-04,  9.3715478409906871e-04,
         3.4754908174042604e-01, -2.8048110764785483e-02, -1.9147096450811166e-01,  1.3769366894035920e-02,  2.2827046419861504e-03, -7.1159030934108405e-03, -1.3179201395291872e-03,  1.8777907284159925e-03,
                                                         
         2.6018094650848411e-01, -2.4623528318145346e-02,  1.5014541606507359e-01, -1.1377222883164428e-02,  5.5564633406300114e-04, -5.4317382757119831e-05,  3.2080256054629734e-04,  2.1992176825352508e-03,
         7.5643742119286284e-01, -4.3151229134044358e-02,  1.3194920980432326e-01, -9.9766889501257138e-04,  2.3686176801920706e-03, -3.4774021071348232e-03,  7.2591693417436444e-04,  8.9950742716123674e-04,
         1.1558500674444006e+00, -3.6801896226935180e-02,  9.5187518915131983e-02,  3.8875289024226109e-03,  5.0780393003259051e-03, -6.1608141969413475e-03,  8.3836836755805694e-04,  2.9847536506072300e-04,
         1.3895291108288046e+00, -2.3292866819651030e-02,  3.6049179432956914e-02,  2.8177193660515751e-03,  7.4784863973144733e-03, -9.0097499591925388e-03,  5.4753041006382476e-04, -5.4908327268816091e-04,
         1.3895291108288044e+00, -2.3292866819650922e-02, -3.6049179432957101e-02, -2.8177193660515378e-03,  7.4784863973144950e-03, -9.0097499591925163e-03, -5.4753041006382487e-04, -5.4908327268817989e-04,
         1.1558500674444003e+00, -3.6801896226934715e-02, -9.5187518915131886e-02, -3.8875289024224817e-03,  5.0780393003259268e-03, -6.1608141969413258e-03, -8.3836836755802658e-04,  2.9847536506067546e-04,
         7.5643742119286161e-01, -4.3151229134043845e-02, -1.3194920980432370e-01,  9.9766889501253452e-04,  2.3686176801920966e-03, -3.4774021071348306e-03, -7.2591693417442061e-04,  8.9950742716124173e-04,
         2.6018094650848284e-01, -2.4623528318144457e-02, -1.5014541606507342e-01,  1.1377222883164968e-02,  5.5564633406298997e-04, -5.4317382757213696e-05, -3.2080256054627907e-04,  2.1992176825354850e-03,
                                                         
         2.2061875008365345e-01,  9.6631517600846272e-03,  1.3141643256470570e-01,  5.1140605928805252e-03,  6.6601405906149182e-03,  3.1310264886015996e-03,  3.8452339628324505e-03, -3.6015854094478215e-04,
         7.0097359003903059e-01,  2.5969251836348413e-02,  1.4059912442543326e-01,  5.4595811792830344e-03,  1.3863694644921848e-02, -9.8784417722103503e-04,  3.1373990954367497e-04,  5.3783951383765680e-04,
         1.1608685714961748e+00,  5.3537968971594505e-02,  1.2037284313135892e-01,  1.0607159984610891e-02,  1.5797402263497132e-02, -4.5111155521692725e-03,  8.0268670457474148e-04,  6.5397859158021802e-04,
         1.4814229261146545e+00,  9.0436045201507859e-02,  5.3722083116000588e-02,  7.5759200922298322e-03,  1.8391205195993732e-02, -1.3014031031221736e-02,  6.9484611672702631e-04, -1.7627894821588743e-03,
         1.4814229261146548e+00,  9.0436045201508108e-02, -5.3722083116000581e-02, -7.5759200922297211e-03,  1.8391205195993788e-02, -1.3014031031221748e-02, -6.9484611672701677e-04, -1.7627894821588838e-03,
         1.1608685714961748e+00,  5.3537968971594407e-02, -1.2037284313135886e-01, -1.0607159984611132e-02,  1.5797402263497066e-02, -4.5111155521693393e-03, -8.0268670457477877e-04,  6.5397859158024664e-04,
         7.0097359003903081e-01,  2.5969251836348378e-02, -1.4059912442543340e-01, -5.4595811792828314e-03,  1.3863694644921699e-02, -9.8784417722108534e-04, -3.1373990954372105e-04,  5.3783951383762341e-04,
         2.2061875008365178e-01,  9.6631517600844936e-03, -1.3141643256470645e-01, -5.1140605928809068e-03,  6.6601405906156555e-03,  3.1310264886016746e-03, -3.8452339628318906e-03, -3.6015854094490277e-04,
                                                         
         3.1728477488866852e-01,  3.1338192874284548e-02,  1.7738358224001591e-01,  1.2875214739490581e-02, -4.8107149369194429e-03, -4.4933329941175887e-03, -2.7774675638251212e-03, -4.0417674588016208e-03,
         8.9395264296783927e-01,  5.3336355793392055e-02,  1.6114696623034069e-01,  4.9642255681518581e-03, -1.1009301949378283e-02, -1.6240840546819605e-04, -8.0128831641333893e-04, -6.1273948817444122e-05,
         1.4515046707409838e+00,  7.3764923340376473e-02,  1.5802986276729358e-01,  7.7325221486260034e-03, -1.5570616292196622e-02, -2.2739467835626980e-03, -1.8321877472712773e-03,  6.3765139919743480e-04,
         1.9036696745231823e+00,  1.0223309369326505e-01,  8.2014576705607287e-02,  6.0293706608141795e-03, -2.1202408094434320e-02, -1.8550660598603982e-02, -1.4193287657706390e-03, -1.4337850889724951e-03,
         1.9036696745231827e+00,  1.0223309369326498e-01, -8.2014576705607231e-02, -6.0293706608142350e-03, -2.1202408094434320e-02, -1.8550660598604118e-02,  1.4193287657706867e-03, -1.4337850889725031e-03,
         1.4515046707409842e+00,  7.3764923340376751e-02, -1.5802986276729358e-01, -7.7325221486257995e-03, -1.5570616292196521e-02, -2.2739467835627301e-03,  1.8321877472712252e-03,  6.3765139919739501e-04,
         8.9395264296783961e-01,  5.3336355793392465e-02, -1.6114696623034064e-01, -4.9642255681518858e-03, -1.1009301949378327e-02, -1.6240840546830605e-04,  8.0128831641336322e-04, -6.1273948817420066e-05,
         3.1728477488866941e-01,  3.1338192874284715e-02, -1.7738358224001569e-01, -1.2875214739490671e-02, -4.8107149369198323e-03, -4.4933329941176468e-03,  2.7774675638248801e-03, -4.0417674588015757e-03,
                                                         
         2.9856868259017222e-01, -4.4358254629948406e-02,  1.6904725884004743e-01, -1.8966632025180216e-02, -6.5259204756554969e-03, -2.5805304029333379e-03, -3.7677419433299738e-03,  5.1461245497284617e-03,
         8.5858667148480361e-01, -7.7423345159837897e-02,  1.5815426669652335e-01, -6.2531566299461995e-03, -1.3850860419404974e-02,  4.2102309020287818e-04, -4.6131410499162384e-04,  3.9811827989688772e-04,
         1.4115464704823038e+00, -9.9620570722380813e-02,  1.5776040893168708e-01, -7.8173884689787462e-03, -1.7728464227987965e-02, -2.1636612802880054e-03, -1.7774214977044575e-03, -5.7397803419405899e-04,
         1.8635410509569048e+00, -1.2804640133826678e-01,  8.2038435164304063e-02, -6.0050462696768147e-03, -2.3251244190819049e-02, -1.8554396836116767e-02, -1.4111570005111832e-03,  1.4316279712387229e-03,
         1.8635410509569053e+00, -1.2804640133826684e-01, -8.2038435164303827e-02,  6.0050462696768702e-03, -2.3251244190819070e-02, -1.8554396836116722e-02,  1.4111570005111642e-03,  1.4316279712387932e-03,
         1.4115464704823049e+00, -9.9620570722380869e-02, -1.5776040893168702e-01,  7.8173884689786716e-03, -1.7728464227987951e-02, -2.1636612802879928e-03,  1.7774214977044958e-03, -5.7397803419400651e-04,
         8.5858667148480450e-01, -7.7423345159837939e-02, -1.5815426669652330e-01,  6.2531566299462082e-03, -1.3850860419404939e-02,  4.2102309020288029e-04,  4.6131410499163864e-04,  3.9811827989687552e-04,
         2.9856868259017277e-01, -4.4358254629948260e-02, -1.6904725884004765e-01,  1.8966632025180313e-02, -6.5259204756555290e-03, -2.5805304029333370e-03,  3.7677419433299339e-03,  5.1461245497284478e-03,
                                                         
         1.6472836762602650e-01, -6.4622802488492528e-03,  1.0615968567242214e-01, -2.0693364430828132e-03,  1.3963870956119633e-02,  8.5621733627069419e-03,  8.0620446554450047e-03,  1.2871184688642244e-03,
         6.0862841609240859e-01, -1.1804099479237192e-02,  1.3886189996819157e-01, -3.6257448687969846e-03,  2.8818588277168294e-02, -1.6305922506323095e-04,  5.1433038859822794e-04, -7.3533836184467519e-04,
         1.0684562036485135e+00, -3.5201132471140248e-02,  1.2101078886542867e-01, -9.9394198301683955e-03,  3.1273725348626248e-02, -4.5077782114690086e-03,  9.0314366050544658e-04, -7.7939850703530583e-04,
         1.3908819207316665e+00, -7.0885965719393249e-02,  5.4053805794843257e-02, -7.4371696780805782e-03,  3.4035995955904931e-02, -1.3096490503921115e-02,  6.9165401818148148e-04,  1.7194957188661963e-03,
         1.3908819207316667e+00, -7.0885965719393138e-02, -5.4053805794842896e-02,  7.4371696780806146e-03,  3.4035995955904966e-02, -1.3096490503921015e-02, -6.9165401818150035e-04,  1.7194957188661772e-03,
         1.0684562036485146e+00, -3.5201132471140234e-02, -1.2101078886542882e-01,  9.9394198301683487e-03,  3.1273725348626262e-02, -4.5077782114689956e-03, -9.0314366050543834e-04, -7.7939850703529628e-04,
         6.0862841609240925e-01, -1.1804099479237469e-02, -1.3886189996819157e-01,  3.6257448687969478e-03,  2.8818588277168287e-02, -1.6305922506316907e-04, -5.1433038859822100e-04, -7.3533836184460841e-04,
         1.6472836762602747e-01, -6.4622802488492458e-03, -1.0615968567242191e-01,  2.0693364430829832e-03,  1.3963870956119570e-02,  8.5621733627069974e-03, -8.0620446554450394e-03,  1.2871184688642598e-03,
                                                         
         3.2567896066608137e-01,  7.9694807939027321e-02,  1.8138329229246233e-01,  3.4130327869920436e-02, -1.2897504262272103e-03, -5.1491649123716135e-03, -7.4463775576943090e-04, -9.2033633129308540e-03,
         8.8573073934652091e-01,  1.3180886462596381e-01,  1.4220407629839993e-01,  5.2114498588327497e-03, -2.1500466509339136e-03, -4.9622366030995175e-03,  2.4794549885218521e-04, -2.0354679892533663e-03,
         1.3060178631336981e+00,  1.3131146816260650e-01,  9.8310629669391358e-02, -3.4378558002398189e-03, -5.2033263638920748e-04, -6.6184622280531470e-03,  6.9297032614725450e-04, -4.3920547811381597e-04,
         1.5458079707863022e+00,  1.1844442680216062e-01,  3.6785251742611789e-02, -2.7477367760390123e-03,  1.5895856233272446e-03, -9.2110386667161575e-03,  5.2519154906816750e-04,  5.2377094526740299e-04,
         1.5458079707863026e+00,  1.1844442680216057e-01, -3.6785251742611727e-02,  2.7477367760389013e-03,  1.5895856233273001e-03, -9.2110386667162338e-03, -5.2519154906811025e-04,  5.2377094526734563e-04,
         1.3060178631336978e+00,  1.3131146816260572e-01, -9.8310629669391830e-02,  3.4378558002395136e-03, -5.2033263638935265e-04, -6.6184622280532918e-03, -6.9297032614738189e-04, -4.3920547811387062e-04,
         8.8573073934651991e-01,  1.3180886462596297e-01, -1.4220407629839971e-01, -5.2114498588324946e-03, -2.1500466509341257e-03, -4.9622366030995643e-03, -2.4794549885213566e-04, -2.0354679892534569e-03,
         3.2567896066608099e-01,  7.9694807939026988e-02, -1.8138329229246228e-01, -3.4130327869920402e-02, -1.2897504262269091e-03, -5.1491649123716274e-03,  7.4463775576965577e-04, -9.2033633129309320e-03,
                                                         
         4.5758071423995789e-01, -1.7203978611669109e-02,  2.4387817433520595e-01, -5.9370449568272301e-03, -1.1872819639640254e-02, -1.5729094905946386e-02, -6.8547756149857598e-03,  3.0950378831330659e-03,
         1.1259773284051073e+00, -2.1648595211938900e-02,  1.5150225754918717e-01, -5.4796944887231913e-04, -2.4261809779577124e-02, -8.3847823258308921e-03, -2.9801117729441520e-04,  5.9460294920611305e-05,
         1.5457214754324125e+00, -2.4079263611516789e-02,  9.2030059599753189e-02, -9.9293276550425479e-04, -2.4657338466129764e-02, -7.4607497917453663e-03,  6.9652583640970205e-05, -4.7089473518950642e-05,
         1.7614161838416049e+00, -2.7278574910696709e-02,  3.1685800667273785e-02, -6.3584621477909753e-04, -2.4216573332765473e-02, -8.0924625745515469e-03,  1.8482328475636194e-04,  1.2203926265291815e-04,
         1.7614161838416049e+00, -2.7278574910696782e-02, -3.1685800667273903e-02,  6.3584621477904201e-04, -2.4216573332765518e-02, -8.0924625745517915e-03, -1.8482328475635242e-04,  1.2203926265290857e-04,
         1.5457214754324113e+00, -2.4079263611516525e-02, -9.2030059599753952e-02,  9.9293276550453235e-04, -2.4657338466129847e-02, -7.4607497917455797e-03, -6.9652583640970558e-05, -4.7089473518941840e-05,
         1.1259773284051051e+00, -2.1648595211938775e-02, -1.5150225754918667e-01,  5.4796944887201382e-04, -2.4261809779577238e-02, -8.3847823258309407e-03,  2.9801117729432342e-04,  5.9460294920663889e-05,
         4.5758071423995744e-01, -1.7203978611668745e-02, -2.4387817433520526e-01,  5.9370449568277878e-03, -1.1872819639639924e-02, -1.5729094905946272e-02,  6.8547756149861354e-03,  3.0950378831332199e-03,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_NEUMANN && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         2.9454103313902792e-03,  1.5802988321940406e-03, -8.1220407017440705e-05, -4.3533533969094911e-05, -9.3133332699638785e-05, -2.0155361167414094e-05,  2.6019398373932300e-06, -1.1636703195470535e-05,
         2.3977491683841680e-03,  1.2864694816254953e-03, -2.3080793010569082e-04, -1.2370484410927240e-04, -7.5811051152946431e-05, -1.6929953988298829e-05,  7.3990840765681019e-06, -9.7745134931449230e-06,
         1.3914765185268068e-03,  7.4658858891277018e-04, -3.4183579882356559e-04, -1.8318732905541317e-04, -4.3982182637681603e-05, -1.0479069282942790e-05,  1.0977321728721260e-05, -6.0500934710058721e-06,
         1.0469460860909087e-04,  5.6204603120590490e-05, -3.8859636591986353e-04, -2.0819395462031662e-04, -3.2849539730672778e-06, -8.0302664144464693e-07,  1.2519234196065596e-05, -4.6362764796970192e-07,
        -1.1954567969507462e-03, -6.4132674268215992e-04, -3.4538854730508360e-04, -1.8490854102460485e-04,  3.7854973652348051e-05,  1.2100089075750871e-05,  1.1232914092909681e-05,  6.9859896853320533e-06,
        -2.1527009307657449e-03, -1.1550160733844750e-03, -1.8646512118587568e-04, -9.9654613071152538e-05,  6.8045519919441081e-05,  2.8220599793438136e-05,  6.1976059213784106e-06,  1.6293170887256732e-05,
        -2.3222472150781618e-03, -1.2453787755312817e-03,  1.1362032473105673e-04,  6.1942147924073555e-05,  7.3874275224166040e-05,  4.7618666498975717e-05, -2.8323724771550466e-06,  2.7492649921608969e-05,
        -1.1659466124930729e-03, -6.2864063967241375e-04,  5.8282714795642951e-04,  3.1079238247111386e-04,  3.4484231093932982e-05,  6.9971209183505883e-05, -1.9909480104826275e-05,  4.0397896457623655e-05,
                                                         
         7.1908026223422268e-03,  8.7835451492739849e-04, -1.9820559523215311e-04, -2.4219787885346473e-05, -8.7265789944371484e-05, -5.0753567962180759e-05,  2.4378090406375824e-06, -6.0291797342080654e-06,
         5.8537844506442640e-03,  7.1503719826896035e-04, -5.6323778393837392e-04, -6.8826478025572259e-05, -7.1036608196370764e-05, -4.2632369469824507e-05,  6.9321134102981558e-06, -5.0647830039377193e-06,
         3.3971475313617555e-03,  4.1495831104339996e-04, -8.3413321262645803e-04, -1.0193410736794554e-04, -4.1215167803821793e-05, -2.6389859593325174e-05,  1.0285303227962514e-05, -3.1360055977131501e-06,
         2.5569422953826664e-04,  3.1227668027690434e-05, -9.4813271870638452e-04, -1.1587769002028600e-04, -3.0892913522580445e-06, -2.0264875026180285e-06,  1.1726681804438307e-05, -2.4273780957173269e-07,
        -2.9183838531053816e-03, -3.5648465667466106e-04, -8.4248196217494016e-04, -1.0298879473282828e-04,  3.5469990746443350e-05,  3.0460216528484907e-05,  1.0535530095005890e-05,  3.6142348418599042e-06,
        -5.2554375925067537e-03, -6.4196595003970119e-04, -4.5440357770413315e-04, -5.5604247794256093e-05,  6.3695248274986838e-05,  7.1059202625021415e-05,  5.7603299370450903e-06,  8.4397079892653313e-06,
        -5.6685849155469203e-03, -6.9247675367087595e-04,  2.7882913999687962e-04,  3.3872716670468030e-05,  6.9345551616678056e-05,  1.1982492778101243e-04, -2.4981257817164883e-06,  1.4195654466751501e-05,
        -2.8535365706562942e-03, -3.4823945630736106e-04,  1.4196182851664520e-03,  1.7380120617798016e-04,  3.2509335419496956e-05,  1.7650875436955982e-04, -1.8769273555653460e-05,  2.1111583934346482e-05,
                                                         
         9.1097511724811519e-03,  2.4472689915191017e-04, -2.5110926166449927e-04, -6.7480177895643811e-06, -7.5510547998584025e-05, -6.4110446364066200e-05,  2.1094887914742097e-06, -1.6824176066535231e-06,
         7.4159316667946106e-03,  1.9922299958803437e-04, -7.1357493292740430e-04, -1.9176063195632704e-05, -6.1467177351129740e-05, -5.3852499098817551e-05,  5.9984550321636152e-06, -1.4131618577037828e-06,
         4.3037083851823667e-03,  1.1561386512741955e-04, -1.0567820436852341e-03, -2.8400502557509977e-05, -3.5662066326326530e-05, -3.3336552782899689e-05,  8.9001327644746506e-06, -8.7466958526552817e-07,
         3.2391372804327112e-04,  8.6985512122039804e-06, -1.2012246270641767e-03, -3.2284200893221456e-05, -2.6712621232576688e-06, -2.5627210997466354e-06,  1.0147116922949568e-05, -6.6856802078990333e-08,
        -3.6972036585589349e-03, -9.9332502127231711e-05, -1.0673962190730847e-03, -2.8697440647027754e-05,  3.0694343198842332e-05,  3.8469492100302152e-05,  9.1165242914399565e-06,  1.0099225655426928e-06,
        -6.6579805374888416e-03, -1.7885544221809676e-04, -5.7579473695344371e-04, -1.5476621325387864e-05,  5.5126005058478898e-05,  8.9758924133812162e-05,  4.9891022599712436e-06,  2.3565812576061996e-06,
        -7.1813177459807169e-03, -1.9299534581858363e-04,  3.5309733185201524e-04,  9.4055531328710514e-06,  5.9976500775005430e-05,  1.5130177056208699e-04, -2.1886672523316694e-06,  3.9775091861413441e-06,
        -3.6146735134397743e-03, -9.6904873658027276e-05,  1.7988543280229861e-03,  4.8442709724383472e-05,  2.8092808946621387e-05,  2.2314456756612132e-04, -1.6219390807545056e-05,  5.8136153685654347e-06,
                                                         
         9.0177617763945398e-03, -2.7507821513658471e-04, -2.4857185935817300e-04,  7.5771827100076464e-06, -5.7881667695644288e-05, -6.3499984706895610e-05,  1.6169963684591646e-06,  2.0348678087510085e-06,
         7.3410449573774903e-03, -2.2393276551209854e-04, -7.0636416570195044e-04,  2.1531261468110171e-05, -4.7116925160056272e-05, -5.3339691789303345e-05,  4.5980306322194561e-06,  1.7092312958943321e-06,
         4.2602459445749237e-03, -1.2995859324763968e-04, -1.0461021797703569e-03,  3.1883966086193099e-05, -2.7336450616165367e-05, -3.3019082421789131e-05,  6.8222316703943620e-06,  1.0579611837125161e-06,
         3.2063450776089993e-04, -9.7868044816153818e-06, -1.1890828329960913e-03,  3.6236050292951021e-05, -2.0477075881192482e-06, -2.5381529188076893e-06,  7.7782309243154608e-06,  8.1041247957582063e-08,
        -3.6598843814004529e-03,  1.1162655413810952e-04, -1.0566030705150260e-03,  3.2180531622019254e-05,  2.3527586818445619e-05,  3.8103190436276860e-05,  6.9876721859188675e-06, -1.2214069298723191e-06,
        -6.5907643841046108e-03,  2.0104688140232032e-04, -5.6995693342905946e-04,  1.7344048752296198e-05,  4.2255442305614336e-05,  8.8904613048026854e-05,  3.8248602209425178e-06, -2.8498179929560688e-06,
        -7.1088551280792788e-03,  2.1676003123261188e-04,  3.4954380907202846e-04, -1.0794657828553013e-05,  4.5978279051241234e-05,  1.4987056547212117e-04, -1.6754794237056668e-06, -4.8038158300984427e-06,
        -3.5782390946798317e-03,  1.0947824330289588e-04,  1.7806367723506291e-03, -5.4075052491789290e-05,  2.1538131781093012e-05,  2.2096185809863191e-04, -1.2435046181748628e-05, -7.0738032671831658e-06,
                                                         
         7.3879864695543191e-03, -6.3552705278391004e-04, -2.0364876374977053e-04,  1.7511442224397156e-05, -3.4375800391546140e-05, -5.2008363013748118e-05,  9.6034257047601726e-07,  4.5998230692127934e-06,
         6.0142995354000800e-03, -5.1736168850744933e-04, -5.7870685045701179e-04,  4.9761142622007696e-05, -2.7982549695155280e-05, -4.3686943026266310e-05,  2.7308024400821805e-06,  3.8637858008648032e-06,
         3.4902869982107905e-03, -3.0024528106584891e-04, -8.5704661018547732e-04,  7.3690637200564862e-05, -1.6234829963312711e-05, -2.7044069648477544e-05,  4.0517467094614445e-06,  2.3917140493704468e-06,
         2.6267253097989252e-04, -2.2603492793735746e-05, -9.7418873251777860e-04,  8.3755249766217447e-05, -1.2157399348283212e-06, -2.0797959256945487e-06,  4.6195289614673165e-06,  1.8359128540123916e-07,
        -2.9984592237183578e-03,  2.5791313968629754e-04, -8.6565427151446292e-04,  7.4400352650767955e-05,  1.3973479886212411e-05,  3.1206398641377672e-05,  4.1499711909910582e-06, -2.7604576694575809e-06,
        -5.3996523945257051e-03,  4.6448890406098546e-04, -4.6696175890602277e-04,  4.0114845103294090e-05,  2.5095769483814482e-05,  7.2812691506560828e-05,  2.2714857021896471e-06, -6.4408572407877029e-06,
        -5.8241006381297899e-03,  5.0088327035022803e-04,  2.8636857138031270e-04, -2.4804063110138882e-05,  2.7302715375744872e-05,  1.2274483508726791e-04, -9.9730489739684769e-07, -1.0857231909561812e-05,
        -2.9315364229559847e-03,  2.5259891595726960e-04,  1.4588572995614490e-03, -1.2518220086764418e-04,  1.2787666311384714e-05,  1.8099693891043193e-04, -7.3829625872122989e-06, -1.5999953584265686e-05,
                                                         
         4.8513131773355034e-03, -7.9109069876015211e-04, -1.3372435092154613e-04,  2.1799769371536328e-05, -4.9943753964592190e-06, -3.4182706218737053e-05,  1.3952806239970073e-07,  5.6918246798688488e-06,
         3.9492814007634778e-03, -6.4400026897535223e-04, -3.8000312818291016e-04,  6.1947290139933206e-05, -4.0654689816510104e-06, -2.8713426374882169e-05,  3.9677630624171602e-07,  4.7811780685270075e-06,
         2.2918874561740788e-03, -3.7373700992882628e-04, -5.6277208845487690e-04,  9.1737969362744377e-05, -2.3586996557912928e-06, -1.7774863906246598e-05,  5.8862742348817751e-07,  2.9598643810806110e-06,
         1.7247110300047785e-04, -2.8132502262290692e-05, -6.3969099838744720e-04,  1.0427002022993412e-04, -1.7636778997424071e-07, -1.3669581959686877e-06,  6.7134246670242169e-07,  2.2796576974449874e-07,
        -1.9689587216544166e-03,  3.2104967547579499e-04, -5.6842050306369447e-04,  9.2627531330544835e-05,  2.0295270284150788e-06,  2.0510428102409549e-05,  6.0223150049873560e-07, -3.4148638004599835e-06,
        -3.5457053086019384e-03,  5.7819596092266944e-04, -3.0661118647039884e-04,  4.9958235643478029e-05,  3.6456508060224863e-06,  4.7856306249641441e-05,  3.3083799754669100e-07, -7.9677185052948154e-06,
        -3.8244625194617798e-03,  6.2348409059883561e-04,  1.8805221813274057e-04, -3.0857811037148491e-05,  3.9685754224339600e-06,  8.0684938806754262e-05, -1.4439738333350025e-07, -1.3426060530080551e-05,
        -1.9250439924812233e-03,  3.1439121673181155e-04,  9.5791921944402036e-04, -1.5588896344641382e-04,  1.8592359090243374e-06,  1.1890480858937924e-04, -1.0734303525460107e-06, -1.9848954571151450e-05,
                                                         
         2.1963273708129776e-03, -6.9624521968044599e-04, -6.0548832950326751e-05,  1.9176551580586097e-05,  3.0265924672526220e-05, -1.5330468825855125e-05, -8.4543670479713512e-07,  5.1925196537383204e-06,
         1.7879490740173164e-03, -5.6679004859827586e-04, -1.7206235631402029e-04,  5.4491471699017110e-05,  2.4637584790419636e-05, -1.2877994824182332e-05, -2.4040868412278171e-06,  4.3614126000035925e-06,
         1.0375939880119739e-03, -3.2892932134363578e-04, -2.5482318341185135e-04,  8.0691546713547118e-05,  1.4295392261313842e-05, -7.9731161290659173e-06, -3.5669807994622816e-06,  2.6991773366034914e-06,
         7.8062845570639860e-05, -2.4760850554311742e-05, -2.8966365570265039e-04,  9.1701427886025828e-05,  1.0732735433543565e-06, -6.1543594761109216e-07, -4.0668130016087937e-06,  2.0592580264676197e-07,
        -8.9143225927158321e-04,  2.8256048136170826e-04, -2.5741176802358226e-04,  8.1437949390533948e-05, -1.2300539453502850e-05,  9.1931546662751908e-06, -3.6545615322183870e-06, -3.1191670643847962e-06,
        -1.6053194730748213e-03,  5.0887146559029759e-04, -1.3891255610688619e-04,  4.3862503654998561e-05, -2.2082467539394576e-05,  2.1461720921825692e-05, -1.9930372813647170e-06, -7.2712024388681924e-06,
        -1.7314991145412693e-03,  5.4868866089157224e-04,  8.5020583194852107e-05, -2.7351028557261420e-05, -2.4072005061806998e-05,  3.6146429344417667e-05,  8.4437725723699960e-07, -1.2288259897304353e-05,
        -8.7122572656270051e-04,  2.7703640955630514e-04,  4.3399802853094253e-04, -1.3678522109887144e-04, -1.1304750375757150e-05,  5.3450562310730834e-05,  6.5268006725663523e-06, -1.7941072137430300e-05,
                                                         
         3.6940088467200675e-04, -3.0543905991253149e-04, -1.0152349145921873e-05,  8.4363937575269777e-06,  7.1390980649458711e-05, -3.1683804831397292e-06, -1.9945339406709136e-06,  1.8292653248364870e-06,
         3.0072239784789743e-04, -2.4864490788882663e-04, -2.8845334482960852e-05,  2.3976219401297706e-05,  5.8112372685613749e-05, -2.6619033040990169e-06, -5.6718739417184526e-06,  1.5368505891824747e-06,
         1.7453220377939475e-04, -1.4429055949825024e-04, -4.2702957969787832e-05,  3.5517819617943386e-05,  3.3713812407597968e-05, -1.6490019217011348e-06, -8.4146414026332655e-06,  9.5205170340452069e-07,
         1.3161502319387232e-05, -1.0847172268078632e-05, -4.8502637779981196e-05,  4.0393097358739686e-05,  2.5161807350256802e-06, -1.2938099737714092e-07, -9.5973196416051322e-06,  7.4698153711759745e-08,
        -1.4987992089873087e-04,  1.2399405082200117e-04, -4.3015520628986321e-05,  3.5948411410556158e-05, -2.9017040454096802e-05,  1.8952994167270875e-06, -8.6083941003507384e-06, -1.0942516284487451e-06,
        -2.6999480213896779e-04,  2.2323425366627318e-04, -2.3055629594970954e-05,  1.9455849036748828e-05, -5.2171163119635596e-05,  4.4338144327928987e-06, -4.7596448534476917e-06, -2.5598639563106976e-06,
        -2.9090801970764989e-04,  2.4103482194758066e-04,  1.4851075567464599e-05, -1.1412857194338685e-05, -5.6606748798579335e-05,  7.4312694328378820e-06,  2.1987582676958820e-06, -4.2904454074693965e-06,
        -1.4930099791391177e-04,  1.2028018408952197e-04,  7.1755510167657109e-05, -6.1104860468761771e-05, -2.6399193882683955e-05,  1.1187856913221623e-05,  1.5241581694556699e-05, -6.4593122005035123e-06,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-9) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_NEUMANN)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
        -1.1659466124925870e-03, -6.2864063967237656e-04, -5.8282714795653001e-04, -3.1079238247113566e-04,  3.4484231093643913e-05,  6.9971209183402152e-05,  1.9909480105078427e-05,  4.0397896457730226e-05,
        -2.3222472150785981e-03, -1.2453787755308061e-03, -1.1362032473107862e-04, -6.1942147924103668e-05,  7.3874275224223611e-05,  4.7618666499189353e-05,  2.8323724771029990e-06,  2.7492649921479451e-05,
        -2.1527009307649976e-03, -1.1550160733849353e-03,  1.8646512118581746e-04,  9.9654613071163177e-05,  6.8045519919442558e-05,  2.8220599793060780e-05, -6.1976059213587680e-06,  1.6293170887530521e-05,
        -1.1954567969514023e-03, -6.4132674268189385e-04,  3.4538854730495214e-04,  1.8490854102458541e-04,  3.7854973652330331e-05,  1.2100089075854109e-05, -1.1232914092940379e-05,  6.9859896852741129e-06,
         1.0469460860846501e-04,  5.6204603120574614e-05,  3.8859636591997499e-04,  2.0819395462014304e-04, -3.2849539730600425e-06, -8.0302664137044272e-07, -1.2519234196020526e-05, -4.6362764805102783e-07,
         1.3914765185255213e-03,  7.4658858891309176e-04,  3.4183579882334495e-04,  1.8318732905560228e-04, -4.3982182637659242e-05, -1.0479069282658081e-05, -1.0977321728757601e-05, -6.0500934712261464e-06,
         2.3977491683835508e-03,  1.2864694816254441e-03,  2.3080793010568892e-04,  1.2370484410934195e-04, -7.5811051153008407e-05, -1.6929953988485048e-05, -7.3990840765804194e-06, -9.7745134929980255e-06,
         2.9454103313896985e-03,  1.5802988321940428e-03,  8.1220407017759609e-05,  4.3533533968818222e-05, -9.3133332699603331e-05, -2.0155361167371271e-05, -2.6019398373247753e-06, -1.1636703195507935e-05,
                                                         
        -2.8535365706552347e-03, -3.4823945630734193e-04, -1.4196182851668080e-03, -1.7380120617789351e-04,  3.2509335418994639e-05,  1.7650875436956266e-04,  1.8769273556070163e-05,  2.1111583934301552e-05,
        -5.6685849155465925e-03, -6.9247675367087530e-04, -2.7882913999694348e-04, -3.3872716670495820e-05,  6.9345551616762041e-05,  1.1982492778101702e-04,  2.4981257816383076e-06,  1.4195654466760426e-05,
        -5.2554375925067303e-03, -6.4196595003967864e-04,  4.5440357770406545e-04,  5.5604247794242500e-05,  6.3695248274972893e-05,  7.1059202625060988e-05, -5.7603299370234622e-06,  8.4397079892324037e-06,
        -2.9183838531056006e-03, -3.5648465667464593e-04,  8.4248196217481829e-04,  1.0298879473287871e-04,  3.5469990746448053e-05,  3.0460216528486120e-05, -1.0535530095016834e-05,  3.6142348418588784e-06,
         2.5569422953777539e-04,  3.1227668027760893e-05,  9.4813271870631828e-04,  1.1587769002030146e-04, -3.0892913522687362e-06, -2.0264875026408742e-06, -1.1726681804436142e-05, -2.4273780954644288e-07,
         3.3971475313611765e-03,  4.1495831104344777e-04,  8.3413321262645413e-04,  1.0193410736792830e-04, -4.1215167803828800e-05, -2.6389859593362826e-05, -1.0285303227962569e-05, -3.1360055976790155e-06,
         5.8537844506435857e-03,  7.1503719826905554e-04,  5.6323778393841013e-04,  6.8826478025545804e-05, -7.1036608196370981e-05, -4.2632369469792659e-05, -6.9321134102941680e-06, -5.0647830039586368e-06,
         7.1908026223417220e-03,  8.7835451492739892e-04,  1.9820559523216205e-04,  2.4219787885353195e-05, -8.7265789944368827e-05, -5.0753567962192333e-05, -2.4378090406396369e-06, -6.0291797342018210e-06,
                                                         
        -3.6146735134388723e-03, -9.6904873657810910e-05, -1.7988543280232621e-03, -4.8442709724702708e-05,  2.8092808946372095e-05,  2.2314456756613335e-04,  1.6219390807745786e-05,  5.8136153686156205e-06,
        -7.1813177459804402e-03, -1.9299534581867776e-04, -3.5309733185210376e-04, -9.4055531328033261e-06,  5.9976500775041378e-05,  1.5130177056209683e-04,  2.1886672522954562e-06,  3.9775091861355224e-06,
        -6.6579805374887948e-03, -1.7885544221810332e-04,  5.7579473695337888e-04,  1.5476621325386200e-05,  5.5126005058467405e-05,  8.9758924133806091e-05, -4.9891022599626954e-06,  2.3565812576124270e-06,
        -3.6972036585590983e-03, -9.9332502127220734e-05,  1.0673962190730370e-03,  2.8697440647032161e-05,  3.0694343198842251e-05,  3.8469492100302328e-05, -9.1165242914415998e-06,  1.0099225655432040e-06,
         3.2391372804294538e-04,  8.6985512122400573e-06,  1.2012246270641422e-03,  3.2284200893222235e-05, -2.6712621232597889e-06, -2.5627210997372160e-06, -1.0147116922949034e-05, -6.6856802085648290e-08,
         4.3037083851819409e-03,  1.1561386512746682e-04,  1.0567820436852111e-03,  2.8400502557515543e-05, -3.5662066326329139e-05, -3.3336552782890148e-05, -8.9001327644753519e-06, -8.7466958527227997e-07,
         7.4159316667941648e-03,  1.9922299958807048e-04,  7.1357493292740094e-04,  1.9176063195631475e-05, -6.1467177351132559e-05, -5.3852499098818825e-05, -5.9984550321632874e-06, -1.4131618577021828e-06,
         9.1097511724807043e-03,  2.4472689915193754e-04,  2.5110926166450220e-04,  6.7480177895583096e-06, -7.5510547998584513e-05, -6.4110446364066688e-05, -2.1094887914729400e-06, -1.6824176066536736e-06,
                                                         
        -3.5782390946800420e-03,  1.0947824330271821e-04, -1.7806367723504401e-03,  5.4075052491834508e-05,  2.1538131781371344e-05,  2.2096185809865541e-04,  1.2435046181529038e-05, -7.0738032672267516e-06,
        -7.1088551280791635e-03,  2.1676003123251349e-04, -3.4954380907207237e-04,  1.0794657828602887e-05,  4.5978279051199770e-05,  1.4987056547211176e-04,  1.6754794237407190e-06, -4.8038158301038112e-06,
        -6.5907643841046420e-03,  2.0104688140230799e-04,  5.6995693342903322e-04, -1.7344048752291499e-05,  4.2255442305622854e-05,  8.8904613048028657e-05, -3.8248602209487630e-06, -2.8498179929576718e-06,
        -3.6598843814005847e-03,  1.1162655413811544e-04,  1.0566030705149933e-03, -3.2180531622011902e-05,  2.3527586818444599e-05,  3.8103190436277388e-05, -6.9876721859180035e-06, -1.2214069298727395e-06,
         3.2063450776066851e-04, -9.7868044815948819e-06,  1.1890828329960655e-03, -3.6236050292948168e-05, -2.0477075881198305e-06, -2.5381529188081099e-06, -7.7782309243160418e-06,  8.1041247958576030e-08,
         4.2602459445746253e-03, -1.2995859324761212e-04,  1.0461021797703442e-03, -3.1883966086191866e-05, -2.7336450616166916e-05, -3.3019082421789388e-05, -6.8222316703943747e-06,  1.0579611837134544e-06,
         7.3410449573771573e-03, -2.2393276551206748e-04,  7.0636416570194393e-04, -2.1531261468110605e-05, -4.7116925160056367e-05, -5.3339691789301366e-05, -4.5980306322190114e-06,  1.7092312958943696e-06,
         9.0177617763941911e-03, -2.7507821513655364e-04,  2.4857185935816975e-04, -7.5771827100067790e-06, -5.7881667695643983e-05, -6.3499984706895298e-05, -1.6169963684593899e-06,  2.0348678087514562e-06,
                                                         
        -2.9315364229562566e-03,  2.5259891595724883e-04, -1.4588572995613476e-03,  1.2518220086766616e-04,  1.2787666311537051e-05,  1.8099693891041217e-04,  7.3829625870840360e-06, -1.5999953584247227e-05,
        -5.8241006381298767e-03,  5.0088327035023312e-04, -2.8636857138029183e-04,  2.4804063110115609e-05,  2.7302715375720833e-05,  1.2274483508725964e-04,  9.9730489742333801e-07, -1.0857231909555901e-05,
        -5.3996523945257779e-03,  4.6448890406097527e-04,  4.6696175890601811e-04, -4.0114845103287584e-05,  2.5095769483824077e-05,  7.2812691506559567e-05, -2.2714857021969849e-06, -6.4408572407879129e-06,
        -2.9984592237184680e-03,  2.5791313968630372e-04,  8.6565427151444742e-04, -7.4400352650764811e-05,  1.3973479886211188e-05,  3.1206398641377354e-05, -4.1499711909899155e-06, -2.7604576694576410e-06,
         2.6267253097972637e-04, -2.2603492793718284e-05,  9.7418873251776332e-04, -8.3755249766214452e-05, -1.2157399348287045e-06, -2.0797959256938334e-06, -4.6195289614679814e-06,  1.8359128540090355e-07,
         3.4902869982105802e-03, -3.0024528106582506e-04,  8.5704661018546778e-04, -7.3690637200563778e-05, -1.6234829963313615e-05, -2.7044069648476379e-05, -4.0517467094610379e-06,  2.3917140493703914e-06,
         6.0142995353998458e-03, -5.1736168850742168e-04,  5.7870685045700724e-04, -4.9761142622006971e-05, -2.7982549695155148e-05, -4.3686943026265144e-05, -2.7308024400822708e-06,  3.8637858008644152e-06,
         7.3879864695540745e-03, -6.3552705278388185e-04,  2.0364876374976860e-04, -1.7511442224397298e-05, -3.4375800391547028e-05, -5.2008363013748226e-05, -9.6034257047616062e-07,  4.5998230692126876e-06,
                                                         
        -1.9250439924812317e-03,  3.1439121673178927e-04, -9.5791921944407717e-04,  1.5588896344645274e-04,  1.8592359090251569e-06,  1.1890480858938231e-04,  1.0734303525356414e-06, -1.9848954571156651e-05,
        -3.8244625194618687e-03,  6.2348409059885263e-04, -1.8805221813273238e-04,  3.0857811037138591e-05,  3.9685754224285652e-06,  8.0684938806755983e-05,  1.4439738334014451e-07, -1.3426060530080747e-05,
        -3.5457053086020095e-03,  5.7819596092267204e-04,  3.0661118647039786e-04, -4.9958235643475495e-05,  3.6456508060250766e-06,  4.7856306249640411e-05, -3.3083799754857883e-07, -7.9677185052943275e-06,
        -1.9689587216545025e-03,  3.2104967547580437e-04,  5.6842050306368851e-04, -9.2627531330543574e-05,  2.0295270284150072e-06,  2.0510428102409799e-05, -6.0223150049840589e-07, -3.4148638004594960e-06,
         1.7247110300037030e-04, -2.8132502262273717e-05,  6.3969099838743961e-04, -1.0427002022993214e-04, -1.7636778997413404e-07, -1.3669581959692599e-06, -6.7134246670271666e-07,  2.2796576974410103e-07,
         2.2918874561739478e-03, -3.7373700992880389e-04,  5.6277208845487137e-04, -9.1737969362743727e-05, -2.3586996557916888e-06, -1.7774863906246866e-05, -5.8862742348828677e-07,  2.9598643810798444e-06,
         3.9492814007633312e-03, -6.4400026897532968e-04,  3.8000312818290761e-04, -6.1947290139932623e-05, -4.0654689816514974e-06, -2.8713426374882010e-05, -3.9677630624155604e-07,  4.7811780685271430e-06,
         4.8513131773353473e-03, -7.9109069876012869e-04,  1.3372435092154383e-04, -2.1799769371536402e-05, -4.9943753964593757e-06, -3.4182706218736694e-05, -1.3952806239973906e-07,  5.6918246798690284e-06,
                                                         
        -8.7122572656279647e-04,  2.7703640955630882e-04, -4.3399802853090686e-04,  1.3678522109887521e-04, -1.1304750375731483e-05,  5.3450562310724884e-05, -6.5268006725849845e-06, -1.7941072137430320e-05,
        -1.7314991145413149e-03,  5.4868866089158492e-04, -8.5020583194851321e-05,  2.7351028557261854e-05, -2.4072005061809010e-05,  3.6146429344417531e-05, -8.4437725723438778e-07, -1.2288259897305080e-05,
        -1.6053194730748712e-03,  5.0887146559030431e-04,  1.3891255610688205e-04, -4.3862503655001990e-05, -2.2082467539394003e-05,  2.1461720921824973e-05,  1.9930372813635900e-06, -7.2712024388685252e-06,
        -8.9143225927163254e-04,  2.8256048136171899e-04,  2.5741176802358351e-04, -8.1437949390531183e-05, -1.2300539453503744e-05,  9.1931546662719044e-06,  3.6545615322186488e-06, -3.1191670643874198e-06,
         7.8062845570585460e-05, -2.4760850554297461e-05,  2.8966365570265337e-04, -9.1701427886020542e-05,  1.0732735433549475e-06, -6.1543594760877584e-07,  4.0668130016094222e-06,  2.0592580264880242e-07,
         1.0375939880119141e-03, -3.2892932134361545e-04,  2.5482318341184620e-04, -8.0691546713548108e-05,  1.4295392261314587e-05, -7.9731161290627782e-06,  3.5669807994617382e-06,  2.6991773366061849e-06,
         1.7879490740172451e-03, -5.6679004859825504e-04,  1.7206235631401552e-04, -5.4491471699019095e-05,  2.4637584790419202e-05, -1.2877994824181733e-05,  2.4040868412278150e-06,  4.3614126000037746e-06,
         2.1963273708128952e-03, -6.9624521968042691e-04,  6.0548832950325430e-05, -1.9176551580585481e-05,  3.0265924672525962e-05, -1.5330468825854478e-05,  8.4543670479707858e-07,  5.1925196537382645e-06,
                                                         
        -1.4930099791391114e-04,  1.2028018408953866e-04, -7.1755510167651973e-05,  6.1104860468769807e-05, -2.6399193882685622e-05,  1.1187856913230159e-05, -1.5241581694552498e-05, -6.4593122004950920e-06,
        -2.9090801970762918e-04,  2.4103482194761058e-04, -1.4851075567476117e-05,  1.1412857194324491e-05, -5.6606748798578000e-05,  7.4312694328325567e-06, -2.1987582676983663e-06, -4.2904454074716641e-06,
        -2.6999480213898167e-04,  2.2323425366628440e-04,  2.3055629594969809e-05, -1.9455849036741353e-05, -5.2171163119637358e-05,  4.4338144327932062e-06,  4.7596448534484058e-06, -2.5598639563097997e-06,
        -1.4987992089880115e-04,  1.2399405082198019e-04,  4.3015520628989580e-05, -3.5948411410557310e-05, -2.9017040454096246e-05,  1.8952994167542490e-06,  8.6083941003513550e-06, -1.0942516284285418e-06,
         1.3161502319421437e-05, -1.0847172268042416e-05,  4.8502637779978140e-05, -4.0393097358750420e-05,  2.5161807350257340e-06, -1.2938099739655025e-07,  9.5973196416042225e-06,  7.4698153697184452e-08,
         1.7453220377942224e-04, -1.4429055949822115e-04,  4.2702957969789695e-05, -3.5517819617936603e-05,  3.3713812407597874e-05, -1.6490019217184632e-06,  8.4146414026340905e-06,  9.5205170339006808e-07,
         3.0072239784788155e-04, -2.4864490788881427e-04,  2.8845334482958003e-05, -2.3976219401295368e-05,  5.8112372685614102e-05, -2.6619033040977036e-06,  5.6718739417178851e-06,  1.5368505891825668e-06,
         3.6940088467198788e-04, -3.0543905991251404e-04,  1.0152349145921298e-05, -8.4363937575268150e-06,  7.1390980649458413e-05, -3.1683804831397380e-06,  1.9945339406711283e-06,  1.8292653248362083e-06,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-9) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_NEUMANN && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         2.9454103313900420e-03, -8.1220407017360528e-05,  1.5802988321942330e-03, -4.3533533969087918e-05, -2.0155361167402259e-05, -9.3133332699637552e-05, -1.1636703195453850e-05,  2.6019398374247413e-06,
         7.1908026223423361e-03, -1.9820559523210294e-04,  8.7835451492740120e-04, -2.4219787885392877e-05, -5.0753567962165160e-05, -8.7265789944373733e-05, -6.0291797342224209e-06,  2.4378090406519917e-06,
         9.1097511724812664e-03, -2.5110926166449499e-04,  2.4472689915191136e-04, -6.7480177895632240e-06, -6.4110446364068518e-05, -7.5510547998585814e-05, -1.6824176066497919e-06,  2.1094887914739878e-06,
         9.0177617763946370e-03, -2.4857185935816872e-04, -2.7507821513659192e-04,  7.5771827100067790e-06, -6.3499984706893780e-05, -5.7881667695644532e-05,  2.0348678087499607e-06,  1.6169963684593903e-06,
         7.3879864695543833e-03, -2.0364876374976792e-04, -6.3552705278392219e-04,  1.7511442224397298e-05, -5.2008363013748253e-05, -3.4375800391546283e-05,  4.5998230692126664e-06,  9.6034257047631012e-07,
         4.8513131773355242e-03, -1.3372435092154242e-04, -7.9109069876016631e-04,  2.1799769371536473e-05, -3.4182706218736457e-05, -4.9943753964593952e-06,  5.6918246798694121e-06,  1.3952806239966510e-07,
         2.1963273708129382e-03, -6.0548832950323140e-05, -6.9624521968046508e-04,  1.9176551580588300e-05, -1.5330468825856626e-05,  3.0265924672526528e-05,  5.1925196537367144e-06, -8.4543670479524952e-07,
         3.6940088467188185e-04, -1.0152349145876114e-05, -3.0543905991257454e-04,  8.4363937575512824e-06, -3.1683804831222062e-06,  7.1390980649449387e-05,  1.8292653248490300e-06, -1.9945339406673447e-06,
                                                         
         2.3977491683849369e-03, -2.3080793010577910e-04,  1.2864694816250503e-03, -1.2370484410924019e-04, -1.6929953988743159e-05, -7.5811051152948952e-05, -9.7745134928128132e-06,  7.3990840765343790e-06,
         5.8537844506443412e-03, -5.6323778393839029e-04,  7.1503719826901152e-04, -6.8826478025534881e-05, -4.2632369469768569e-05, -7.1036608196368880e-05, -5.0647830039808468e-06,  6.9321134102861542e-06,
         7.4159316667947684e-03, -7.1357493292738706e-04,  1.9922299958802651e-04, -1.9176063195635161e-05, -5.3852499098824103e-05, -6.1467177351131122e-05, -1.4131618576971389e-06,  5.9984550321636160e-06,
         7.3410449573776170e-03, -7.0636416570193830e-04, -2.2393276551210901e-04,  2.1531261468109158e-05, -5.3339691789299821e-05, -4.7116925160056631e-05,  1.7092312958932896e-06,  4.5980306322191393e-06,
         6.0142995354001616e-03, -5.7870685045700615e-04, -5.1736168850746462e-04,  4.9761142622005094e-05, -4.3686943026266249e-05, -2.7982549695155219e-05,  3.8637858008641383e-06,  2.7308024400818845e-06,
         3.9492814007635029e-03, -3.8000312818290875e-04, -6.4400026897536893e-04,  6.1947290139933856e-05, -2.8713426374879743e-05, -4.0654689816513000e-06,  4.7811780685292284e-06,  3.9677630624177600e-07,
         1.7879490740172811e-03, -1.7206235631403549e-04, -5.6679004859829050e-04,  5.4491471699005909e-05, -1.2877994824194333e-05,  2.4637584790422119e-05,  4.3614125999929656e-06, -2.4040868412284838e-06,
         3.0072239784795717e-04, -2.8845334482907181e-05, -2.4864490788874721e-04,  2.3976219401359397e-05, -2.6619033040869899e-06,  5.8112372685623846e-05,  1.5368505892068682e-06, -5.6718739417108361e-06,
                                                         
         1.3914765185258561e-03, -3.4183579882341906e-04,  7.4658858891353021e-04, -1.8318732905546641e-04, -1.0479069282573469e-05, -4.3982182637698754e-05, -6.0500934712290136e-06,  1.0977321728746634e-05,
         3.3971475313620747e-03, -8.3413321262638485e-04,  4.1495831104337492e-04, -1.0193410736796484e-04, -2.6389859593333692e-05, -4.1215167803837290e-05, -3.1360055977081860e-06,  1.0285303227964386e-05,
         4.3037083851825975e-03, -1.0567820436851994e-03,  1.1561386512741261e-04, -2.8400502557514676e-05, -3.3336552782898530e-05, -3.5662066326327438e-05, -8.7466958526477833e-07,  8.9001327644751503e-06,
         4.2602459445751058e-03, -1.0461021797703379e-03, -1.2995859324766039e-04,  3.1883966086188762e-05, -3.3019082421787586e-05, -2.7336450616166147e-05,  1.0579611837119236e-06,  6.8222316703947881e-06,
         3.4902869982108981e-03, -8.5704661018546821e-04, -3.0024528106587151e-04,  7.3690637200562979e-05, -2.7044069648478385e-05, -1.6234829963312907e-05,  2.3917140493696006e-06,  4.0517467094615504e-06,
         2.2918874561741035e-03, -5.6277208845487245e-04, -3.7373700992885094e-04,  9.1737969362743320e-05, -1.7774863906239372e-05, -2.3586996557913038e-06,  2.9598643810861230e-06,  5.8862742348810996e-07,
         1.0375939880119213e-03, -2.5482318341186810e-04, -3.2892932134365936e-04,  8.0691546713532766e-05, -7.9731161290950010e-06,  1.4295392261310940e-05,  2.6991773365770677e-06, -3.5669807994648062e-06,
         1.7453220377912048e-04, -4.2702957969845511e-05, -1.4429055949834223e-04,  3.5517819617925015e-05, -1.6490019215428654e-06,  3.3713812407604711e-05,  9.5205170353912450e-07, -8.4146414026428166e-06,
                                                         
         1.0469460861058887e-04, -3.8859636591975815e-04,  5.6204603120045618e-05, -2.0819395462020145e-04, -8.0302664197519809e-07, -3.2849539730679389e-06, -4.6362764756131176e-07,  1.2519234196049730e-05,
         2.5569422953876234e-04, -9.4813271870625344e-04,  3.1227668027665605e-05, -1.1587769002035187e-04, -2.0264875025471974e-06, -3.0892913522517172e-06, -2.4273780963290975e-07,  1.1726681804449046e-05,
         3.2391372804365488e-04, -1.2012246270641442e-03,  8.6985512121576477e-06, -3.2284200893225467e-05, -2.5627210997615750e-06, -2.6712621232564135e-06, -6.6856802067339964e-08,  1.0147116922950308e-05,
         3.2063450776114599e-04, -1.1890828329960733e-03, -9.7868044816486329e-06,  3.6236050292946501e-05, -2.5381529188056594e-06, -2.0477075881180750e-06,  8.1041247955739184e-08,  7.7782309243161333e-06,
         2.6267253098003472e-04, -9.7418873251776874e-04, -2.2603492793763675e-05,  8.3755249766216688e-05, -2.0797959256965807e-06, -1.2157399348280904e-06,  1.8359128540073999e-07,  4.6195289614675409e-06,
         1.7247110300052504e-04, -6.3969099838744059e-04, -2.8132502262318143e-05,  1.0427002022993293e-04, -1.3669581959628747e-06, -1.7636778997442399e-07,  2.2796576974953898e-07,  6.7134246670252301e-07,
         7.8062845570577586e-05, -2.8966365570262740e-04, -2.4760850554352966e-05,  9.1701427886038458e-05, -6.1543594763093105e-07,  1.0732735433499284e-06,  2.0592580262690320e-07, -4.0668130016071132e-06,
         1.3161502319080330e-05, -4.8502637779986671e-05, -1.0847172268193755e-05,  4.0393097358703731e-05, -1.2938099725318061e-07,  2.5161807350095484e-06,  7.4698153814639228e-08, -9.5973196416087863e-06,
                                                         
        -1.1954567969495564e-03, -3.4538854730492959e-04, -6.4132674268225500e-04, -1.8490854102475103e-04,  1.2100089075559087e-05,  3.7854973652365025e-05,  6.9859896855152445e-06,  1.1232914092935815e-05,
        -2.9183838531045923e-03, -8.4248196217493842e-04, -3.5648465667479832e-04, -1.0298879473281048e-04,  3.0460216528527357e-05,  3.5469990746459681e-05,  3.6142348418117918e-06,  1.0535530095000664e-05,
        -3.6972036585584686e-03, -1.0673962190730737e-03, -9.9332502127297834e-05, -2.8697440647033212e-05,  3.8469492100284032e-05,  3.0694343198845416e-05,  1.0099225655558758e-06,  9.1165242914400633e-06,
        -3.6598843814001701e-03, -1.0566030705150247e-03,  1.1162655413806796e-04,  3.2180531622018800e-05,  3.8103190436277869e-05,  2.3527586818447662e-05, -1.2214069298743412e-06,  6.9876721859187693e-06,
        -2.9984592237181948e-03, -8.6565427151446151e-04,  2.5791313968626838e-04,  7.4400352650768673e-05,  3.1206398641375273e-05,  1.3973479886213371e-05, -2.7604576694575534e-06,  4.1499711909912555e-06,
        -1.9689587216543394e-03, -5.6842050306368872e-04,  3.2104967547577504e-04,  9.2627531330546922e-05,  2.0510428102411361e-05,  2.0295270284158576e-06, -3.4148638004576147e-06,  6.0223150049910437e-07,
        -8.9143225927156673e-04, -2.5741176802354974e-04,  2.8256048136169465e-04,  8.1437949390549534e-05,  9.1931546662631765e-06, -1.2300539453500705e-05, -3.1191670643951643e-06, -3.6545615322162961e-06,
        -1.4987992089868541e-04, -4.3015520628897802e-05,  1.2399405082203611e-04,  3.5948411410591198e-05,  1.8952994167578119e-06, -2.9017040454091117e-05, -1.0942516284136823e-06, -8.6083941003344906e-06,
                                                         
        -2.1527009307644061e-03, -1.8646512118600898e-04, -1.1550160733847874e-03, -9.9654613071125067e-05,  2.8220599793195932e-05,  6.8045519919463673e-05,  1.6293170887445170e-05,  6.1976059213554925e-06,
        -5.2554375925060745e-03, -4.5440357770421236e-04, -6.4196595003978554e-04, -5.5604247794220532e-05,  7.1059202625052870e-05,  6.3695248274996583e-05,  8.4397079892349228e-06,  5.7603299370468250e-06,
        -6.6579805374883992e-03, -5.7579473695345835e-04, -1.7885544221815801e-04, -1.5476621325388223e-05,  8.9758924133802350e-05,  5.5126005058482273e-05,  2.3565812576127653e-06,  4.9891022599713690e-06,
        -6.5907643841043436e-03, -5.6995693342907160e-04,  2.0104688140227861e-04,  1.7344048752298075e-05,  8.8904613048026841e-05,  4.2255442305616389e-05, -2.8498179929568506e-06,  3.8248602209426042e-06,
        -5.3996523945255568e-03, -4.6696175890603220e-04,  4.6448890406095863e-04,  4.0114845103294232e-05,  7.2812691506559039e-05,  2.5095769483816823e-05, -6.4408572407880002e-06,  2.2714857021901650e-06,
        -3.5457053086018543e-03, -3.0661118647040393e-04,  5.7819596092265860e-04,  4.9958235643479764e-05,  4.7856306249640905e-05,  3.6456508060240237e-06, -7.9677185052935991e-06,  3.3083799754672944e-07,
        -1.6053194730747523e-03, -1.3891255610688300e-04,  5.0887146559030246e-04,  4.3862503654999753e-05,  2.1461720921817611e-05, -2.2082467539391215e-05, -7.2712024388738133e-06, -1.9930372813661214e-06,
        -2.6999480213882576e-04, -2.3055629595024273e-05,  2.2323425366631928e-04,  1.9455849036699280e-05,  4.4338144328077108e-06, -5.2171163119625316e-05, -2.5598639562918557e-06, -4.7596448534612875e-06,
                                                         
        -2.3222472150779250e-03,  1.1362032473092164e-04, -1.2453787755311691e-03,  6.1942147924021161e-05,  4.7618666499018679e-05,  7.3874275224167368e-05,  2.7492649921588163e-05, -2.8323724771444113e-06,
        -5.6685849155464563e-03,  2.7882913999679830e-04, -6.9247675367089200e-04,  3.3872716670508288e-05,  1.1982492778101539e-04,  6.9345551616651534e-05,  1.4195654466749233e-05, -2.4981257817392477e-06,
        -7.1813177459803482e-03,  3.5309733185200396e-04, -1.9299534581860656e-04,  9.4055531328878938e-06,  1.5130177056208970e-04,  5.9976500774992101e-05,  3.9775091861434362e-06, -2.1886672523415170e-06,
        -7.1088551280790550e-03,  3.4954380907201540e-04,  2.1676003123256263e-04, -1.0794657828562771e-05,  1.4987056547212130e-04,  4.5978279051237229e-05, -4.8038158301019511e-06, -1.6754794237092915e-06,
        -5.8241006381297075e-03,  2.8636857138028245e-04,  5.0088327035020840e-04, -2.4804063110132956e-05,  1.2274483508726428e-04,  2.7302715375750029e-05, -1.0857231909560589e-05, -9.9730489739567751e-07,
        -3.8244625194617308e-03,  1.8805221813272528e-04,  6.2348409059883214e-04, -3.0857811037146757e-05,  8.0684938806753625e-05,  3.9685754224360386e-06, -1.3426060530079875e-05, -1.4439738333316115e-07,
        -1.7314991145412275e-03,  8.5020583194843203e-05,  5.4868866089156985e-04, -2.7351028557260157e-05,  3.6146429344417355e-05, -2.4072005061806466e-05, -1.2288259897304829e-05,  8.4437725723681421e-07,
        -2.9090801970761254e-04,  1.4851075567444313e-05,  2.4103482194757101e-04, -1.1412857194341387e-05,  7.4312694328424628e-06, -5.6606748798586050e-05, -4.2904454074661160e-06,  2.1987582676996568e-06,
                                                         
        -1.1659466124931000e-03,  5.8282714795635871e-04, -6.2864063967225589e-04,  3.1079238247121442e-04,  6.9971209183506967e-05,  3.4484231093877050e-05,  4.0397896457619813e-05, -1.9909480104869962e-05,
        -2.8535365706562447e-03,  1.4196182851663058e-03, -3.4823945630721974e-04,  1.7380120617804703e-04,  1.7650875436957223e-04,  3.2509335419638296e-05,  2.1111583934356826e-05, -1.8769273555533880e-05,
        -3.6146735134391776e-03,  1.7988543280231606e-03, -9.6904873658068883e-05,  4.8442709724350444e-05,  2.2314456756614801e-04,  2.8092808946594930e-05,  5.8136153685633849e-06, -1.6219390807542593e-05,
        -3.5782390946797293e-03,  1.7806367723505630e-03,  1.0947824330273572e-04, -5.4075052491861714e-05,  2.2096185809862445e-04,  2.1538131781132227e-05, -7.0738032672007858e-06, -1.2435046181720210e-05,
        -2.9315364229560697e-03,  1.4588572995613673e-03,  2.5259891595728174e-04, -1.2518220086761651e-04,  1.8099693891041654e-04,  1.2787666311392866e-05, -1.5999953584252699e-05, -7.3829625872116297e-06,
        -1.9250439924812228e-03,  9.5791921944401396e-04,  3.1439121673182614e-04, -1.5588896344640823e-04,  1.1890480858938340e-04,  1.8592359090147261e-06, -1.9848954571153225e-05, -1.0734303525530883e-06,
        -8.7122572656269011e-04,  4.3399802853093461e-04,  2.7703640955630682e-04, -1.3678522109886976e-04,  5.3450562310731660e-05, -1.1304750375758454e-05, -1.7941072137430394e-05,  6.5268006725655095e-06,
        -1.4930099791390629e-04,  7.1755510167653518e-05,  1.2028018408951893e-04, -6.1104860468759738e-05,  1.1187856913221945e-05, -2.6399193882684073e-05, -6.4593122005037097e-06,  1.5241581694556747e-05,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 5e-8) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_NEUMANN) &&
               (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] = {
        -1.1659466124932934e-03, -5.8282714795684855e-04, -6.2864063967216005e-04, -3.1079238247093031e-04,  6.9971209183639443e-05,  3.4484231093707806e-05,  4.0397896457465917e-05,  1.9909480104938588e-05,
        -2.8535365706562591e-03, -1.4196182851663028e-03, -3.4823945630696554e-04, -1.7380120617813520e-04,  1.7650875436955090e-04,  3.2509335419660020e-05,  2.1111583934421966e-05,  1.8769273555533724e-05,
        -3.6146735134380812e-03, -1.7988543280237509e-03, -9.6904873658122144e-05, -4.8442709724320046e-05,  2.2314456756617946e-04,  2.8092808946275449e-05,  5.8136153685287285e-06,  1.6219390807763102e-05,
        -3.5782390946801266e-03, -1.7806367723502590e-03,  1.0947824330270692e-04,  5.4075052491799996e-05,  2.2096185809862188e-04,  2.1538131781417463e-05, -7.0738032671858941e-06,  1.2435046181516502e-05,
        -2.9315364229558867e-03, -1.4588572995615559e-03,  2.5259891595718768e-04,  1.2518220086770341e-04,  1.8099693891041608e-04,  1.2787666311323092e-05, -1.5999953584266299e-05,  7.3829625872477057e-06,
        -1.9250439924815151e-03, -9.5791921944390847e-04,  3.1439121673182630e-04,  1.5588896344638343e-04,  1.1890480858937516e-04,  1.8592359090849746e-06, -1.9848954571144047e-05,  1.0734303525058341e-06,
        -8.7122572656279994e-04, -4.3399802853098210e-04,  2.7703640955631728e-04,  1.3678522109886429e-04,  5.3450562310724146e-05, -1.1304750375761637e-05, -1.7941072137439224e-05, -6.5268006725678338e-06,
        -1.4930099791394933e-04, -7.1755510167622631e-05,  1.2028018408956847e-04,  6.1104860468814029e-05,  1.1187856913243012e-05, -2.6399193882670721e-05, -6.4593122004783715e-06, -1.5241581694556352e-05,
                                                         
        -2.3222472150787334e-03, -1.1362032473084332e-04, -1.2453787755309766e-03, -6.1942147924150018e-05,  4.7618666499106953e-05,  7.3874275224159575e-05,  2.7492649921511310e-05,  2.8323724771691171e-06,
        -5.6685849155467226e-03, -2.7882913999692635e-04, -6.9247675367075322e-04, -3.3872716670541539e-05,  1.1982492778100981e-04,  6.9345551616658378e-05,  1.4195654466771952e-05,  2.4981257817306710e-06,
        -7.1813177459804046e-03, -3.5309733185209666e-04, -1.9299534581860401e-04, -9.4055531328531265e-06,  1.5130177056210770e-04,  5.9976500775014226e-05,  3.9775091861342900e-06,  2.1886672523183790e-06,
        -7.1088551280790602e-03, -3.4954380907209324e-04,  2.1676003123252243e-04,  1.0794657828594574e-05,  1.4987056547211805e-04,  4.5978279051207603e-05, -4.8038158301052012e-06,  1.6754794237310099e-06,
        -5.8241006381298463e-03, -2.8636857138028348e-04,  5.0088327035020753e-04,  2.4804063110122039e-05,  1.2274483508726146e-04,  2.7302715375748765e-05, -1.0857231909557101e-05,  9.9730489739923484e-07,
        -3.8244625194618865e-03, -1.8805221813273491e-04,  6.2348409059881881e-04,  3.0857811037153045e-05,  8.0684938806758653e-05,  3.9685754224312266e-06, -1.3426060530079016e-05,  1.4439738333697229e-07,
        -1.7314991145413944e-03, -8.5020583194847431e-05,  5.4868866089157560e-04,  2.7351028557253326e-05,  3.6146429344395346e-05, -2.4072005061812192e-05, -1.2288259897321156e-05, -8.4437725723589867e-07,
        -2.9090801970794327e-04, -1.4851075567502307e-05,  2.4103482194748590e-04,  1.1412857194310117e-05,  7.4312694329712092e-06, -5.6606748798579951e-05, -4.2904454073627322e-06, -2.1987582677042401e-06,
                                                         
        -2.1527009307651859e-03,  1.8646512118598223e-04, -1.1550160733848074e-03,  9.9654613071068146e-05,  2.8220599793311905e-05,  6.8045519919470124e-05,  1.6293170887319514e-05, -6.1976059213721782e-06,
        -5.2554375925066713e-03,  4.5440357770411104e-04, -6.4196595003967701e-04,  5.5604247794256595e-05,  7.1059202625017119e-05,  6.3695248274989264e-05,  8.4397079892732832e-06, -5.7603299370463194e-06,
        -6.6579805374887124e-03,  5.7579473695340045e-04, -1.7885544221810652e-04,  1.5476621325382080e-05,  8.9758924133817827e-05,  5.5126005058471634e-05,  2.3565812576037330e-06, -4.9891022599671389e-06,
        -6.5907643841045579e-03,  5.6995693342903311e-04,  2.0104688140229823e-04, -1.7344048752290777e-05,  8.8904613048027532e-05,  4.2255442305617853e-05, -2.8498179929565503e-06, -3.8248602209461516e-06,
        -5.3996523945257346e-03,  4.6696175890601562e-04,  4.6448890406096058e-04, -4.0114845103287008e-05,  7.2812691506559445e-05,  2.5095769483818270e-05, -6.4408572407884703e-06, -2.2714857021922919e-06,
        -3.5457053086020277e-03,  3.0661118647039792e-04,  5.7819596092265827e-04, -4.9958235643478964e-05,  4.7856306249641414e-05,  3.6456508060242685e-06, -7.9677185052931044e-06, -3.3083799754739176e-07,
        -1.6053194730749349e-03,  1.3891255610690820e-04,  5.0887146559029466e-04, -4.3862503654979696e-05,  2.1461720921818943e-05, -2.2082467539392662e-05, -7.2712024388737879e-06,  1.9930372813676503e-06,
        -2.6999480213906927e-04,  2.3055629594976582e-05,  2.2323425366628367e-04, -1.9455849036767998e-05,  4.4338144328155839e-06, -5.2171163119633333e-05, -2.5598639562880742e-06,  4.7596448534577155e-06,
                                                         
        -1.1954567969511302e-03,  3.4538854730502923e-04, -6.4132674268183737e-04,  1.8490854102465884e-04,  1.2100089076086387e-05,  3.7854973652334918e-05,  6.9859896850783830e-06, -1.1232914092940140e-05,
        -2.9183838531052914e-03,  8.4248196217491620e-04, -3.5648465667469651e-04,  1.0298879473282870e-04,  3.0460216528442132e-05,  3.5469990746440788e-05,  3.6142348418951585e-06, -1.0535530095007782e-05,
        -3.6972036585589322e-03,  1.0673962190730561e-03, -9.9332502127242391e-05,  2.8697440647027754e-05,  3.8469492100309260e-05,  3.0694343198843362e-05,  1.0099225655362619e-06, -9.1165242914393501e-06,
        -3.6598843814004858e-03,  1.0566030705150013e-03,  1.1162655413809776e-04, -3.2180531622016381e-05,  3.8103190436275193e-05,  2.3527586818444731e-05, -1.2214069298709838e-06, -6.9876721859177960e-06,
        -2.9984592237184285e-03,  8.6565427151444893e-04,  2.5791313968628822e-04, -7.4400352650765244e-05,  3.1206398641377970e-05,  1.3973479886212323e-05, -2.7604576694577499e-06, -4.1499711909905830e-06,
        -1.9689587216545142e-03,  5.6842050306368818e-04,  3.2104967547578957e-04, -9.2627531330544117e-05,  2.0510428102407723e-05,  2.0295270284153879e-06, -3.4148638004610372e-06, -6.0223150049890797e-07,
        -8.9143225927170561e-04,  2.5741176802357402e-04,  2.8256048136170148e-04, -8.1437949390539152e-05,  9.1931546662830869e-06, -1.2300539453500771e-05, -3.1191670643781148e-06,  3.6545615322155939e-06,
        -1.4987992089871799e-04,  4.3015520628923335e-05,  1.2399405082208473e-04, -3.5948411410592865e-05,  1.8952994166541755e-06, -2.9017040454095823e-05, -1.0942516285020779e-06,  8.6083941003399760e-06,
                                                         
         1.0469460860968653e-04,  3.8859636591991498e-04,  5.6204603120191877e-05,  2.0819395462023029e-04, -8.0302664154954582e-07, -3.2849539730121296e-06, -4.6362764792594457e-07, -1.2519234195995736e-05,
         2.5569422953822263e-04,  9.4813271870636913e-04,  3.1227668027657047e-05,  1.1587769002025282e-04, -2.0264875026315505e-06, -3.0892913522514445e-06, -2.4273780956273466e-07, -1.1726681804431036e-05,
         3.2391372804316346e-04,  1.2012246270641329e-03,  8.6985512121944022e-06,  3.2284200893228828e-05, -2.5627210997466735e-06, -2.6712621232558642e-06, -6.6856802080207916e-08, -1.0147116922949568e-05,
         3.2063450776077617e-04,  1.1890828329960644e-03, -9.7868044816164677e-06, -3.6236050292948690e-05, -2.5381529188086757e-06, -2.0477075881185540e-06,  8.1041247958263276e-08, -7.7782309243155539e-06,
         2.6267253097976892e-04,  9.7418873251776191e-04, -2.2603492793735177e-05, -8.3755249766214072e-05, -2.0797959256942010e-06, -1.2157399348280701e-06,  1.8359128540133437e-07, -4.6195289614675511e-06,
         1.7247110300035482e-04,  6.3969099838744070e-04, -2.8132502262291485e-05, -1.0427002022993157e-04, -1.3669581959691951e-06, -1.7636778997442494e-07,  2.2796576974392860e-07, -6.7134246670248564e-07,
         7.8062845570517656e-05,  2.8966365570265147e-04, -2.4760850554313606e-05, -9.1701427886023686e-05, -6.1543594760887685e-07,  1.0732735433517474e-06,  2.0592580264888424e-07,  4.0668130016089004e-06,
         1.3161502319048458e-05,  4.8502637779995114e-05, -1.0847172268212027e-05, -4.0393097358729251e-05, -1.2938099733080721e-07,  2.5161807350151427e-06,  7.4698153735111159e-08,  9.5973196416092420e-06,
                                                         
         1.3914765185263169e-03,  3.4183579882302766e-04,  7.4658858891292273e-04,  1.8318732905578478e-04, -1.0479069282939131e-05, -4.3982182637717829e-05, -6.0500934709877439e-06, -1.0977321728843863e-05,
         3.3971475313614463e-03,  8.3413321262635395e-04,  4.1495831104340495e-04,  1.0193410736797270e-04, -2.6389859593312221e-05, -4.1215167803816081e-05, -3.1360055977260237e-06, -1.0285303227970259e-05,
         4.3037083851820979e-03,  1.0567820436851799e-03,  1.1561386512743094e-04,  2.8400502557520168e-05, -3.3336552782903911e-05, -3.5662066326326869e-05, -8.7466958526250129e-07, -8.9001327644758008e-06,
         4.2602459445746964e-03,  1.0461021797703277e-03, -1.2995859324762573e-04, -3.1883966086187854e-05, -3.3019082421788250e-05, -2.7336450616165262e-05,  1.0579611837122543e-06, -6.8222316703946466e-06,
         3.4902869982106101e-03,  8.5704661018546214e-04, -3.0024528106583763e-04, -7.3690637200561745e-05, -2.7044069648476582e-05, -1.6234829963313171e-05,  2.3917140493707772e-06, -4.0517467094616275e-06,
         2.2918874561739383e-03,  5.6277208845487224e-04, -3.7373700992881511e-04, -9.1737969362741450e-05, -1.7774863906247947e-05, -2.3586996557919924e-06,  2.9598643810789287e-06, -5.8862742348835252e-07,
         1.0375939880118588e-03,  2.5482318341186436e-04, -3.2892932134363068e-04, -8.0691546713537333e-05, -7.9731161290558715e-06,  1.4295392261314382e-05,  2.6991773366117605e-06,  3.5669807994639816e-06,
         1.7453220377943444e-04,  4.2702957969832663e-05, -1.4429055949816410e-04, -3.5517819617933093e-05, -1.6490019217785784e-06,  3.3713812407599805e-05,  9.5205170334574603e-07,  8.4146414026363046e-06,
                                                         
         2.3977491683836250e-03,  2.3080793010601619e-04,  1.2864694816253582e-03,  1.2370484410882999e-04, -1.6929953988436089e-05, -7.5811051152950524e-05, -9.7745134930519137e-06, -7.3990840764270243e-06,
         5.8537844506436065e-03,  5.6323778393827656e-04,  7.1503719826903863e-04,  6.8826478025595203e-05, -4.2632369469812561e-05, -7.1036608196369422e-05, -5.0647830039444473e-06, -6.9321134102929034e-06,
         7.4159316667941665e-03,  7.1357493292736321e-04,  1.9922299958807395e-04,  1.9176063195635090e-05, -5.3852499098816223e-05, -6.1467177351132030e-05, -1.4131618577035615e-06, -5.9984550321638634e-06,
         7.3410449573771677e-03,  7.0636416570192323e-04, -2.2393276551206575e-04, -2.1531261468104317e-05, -5.3339691789301522e-05, -4.7116925160056001e-05,  1.7092312958943399e-06, -4.5980306322195806e-06,
         6.0142995353998537e-03,  5.7870685045700052e-04, -5.1736168850742613e-04, -4.9761142622004085e-05, -4.3686943026265049e-05, -2.7982549695155568e-05,  3.8637858008649633e-06, -2.7308024400819751e-06,
         3.9492814007633268e-03,  3.8000312818290967e-04, -6.4400026897533217e-04, -6.1947290139930319e-05, -2.8713426374882667e-05, -4.0654689816512509e-06,  4.7811780685258361e-06, -3.9677630624124989e-07,
         1.7879490740172357e-03,  1.7206235631401760e-04, -5.6679004859825255e-04, -5.4491471699023364e-05, -1.2877994824179879e-05,  2.4637584790421055e-05,  4.3614126000064665e-06,  2.4040868412266563e-06,
         3.0072239784787770e-04,  2.8845334482970725e-05, -2.4864490788881172e-04, -2.3976219401285248e-05, -2.6619033041074775e-06,  5.8112372685617497e-05,  1.5368505891731717e-06,  5.6718739417165061e-06,
                                                         
         2.9454103313885519e-03,  8.1220407017042382e-05,  1.5802988321944147e-03,  4.3533533969364281e-05, -2.0155361167078015e-05, -9.3133332699643664e-05, -1.1636703195740359e-05, -2.6019398375347192e-06,
         7.1908026223413629e-03,  1.9820559523205325e-04,  8.7835451492750843e-04,  2.4219787885412755e-05, -5.0753567962228477e-05, -8.7265789944387137e-05, -6.0291797341598904e-06, -2.4378090406525533e-06,
         9.1097511724805760e-03,  2.5110926166448328e-04,  2.4472689915198118e-04,  6.7480177895652484e-06, -6.4110446364049639e-05, -7.5510547998588145e-05, -1.6824176066647988e-06, -2.1094887914745764e-06,
         9.0177617763941547e-03,  2.4857185935816363e-04, -2.7507821513654036e-04, -7.5771827100049003e-06, -6.3499984706895637e-05, -5.7881667695645867e-05,  2.0348678087525679e-06, -1.6169963684599148e-06,
         7.3879864695540632e-03,  2.0364876374976535e-04, -6.3552705278387925e-04, -1.7511442224397156e-05, -5.2008363013747210e-05, -3.4375800391547042e-05,  4.5998230692123530e-06, -9.6034257047638593e-07,
         4.8513131773353421e-03,  1.3372435092154266e-04, -7.9109069876012804e-04, -2.1799769371535460e-05, -3.4182706218736620e-05, -4.9943753964595587e-06,  5.6918246798688327e-06, -1.3952806240011265e-07,
         2.1963273708128931e-03,  6.0548832950324576e-05, -6.9624521968042626e-04, -1.9176551580585735e-05, -1.5330468825854658e-05,  3.0265924672525651e-05,  5.1925196537383627e-06,  8.4543670479709690e-07,
         3.6940088467198571e-04,  1.0152349145921915e-05, -3.0543905991251355e-04, -8.4363937575262187e-06, -3.1683804831399350e-06,  7.1390980649458819e-05,  1.8292653248360236e-06,  1.9945339406707832e-06,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-7) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(rho_cellavg);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(perbuff);
}

void evalFunc2x_dirichletvarx_dirichletvary(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = -6.0;
}

void evalFunc2x_dirichletvarx_dirichletvary_bc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 1.0 + x*x + 2.0*y*y;
}

void
test_2x_varBC(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  // This test is meant to replicate a test for the Fenics Poisson solver here:
  //   https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1003.html
  double epsilon_0 = 1.0;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc2x_dirichletvarx_dirichletvary, NULL);

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
    phi_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);

  gkyl_array_clear(epsilon, 0.);
  gkyl_array_shiftc(epsilon, epsilon_0*pow(sqrt(2.),dim), 0);

  // Project the right-side source on the basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);
//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, rho_ho, "ctest_fem_poisson_2x_rho_1.gkyl");

  // Project the BC:
  struct gkyl_array *phibc = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phibc_ho = use_gpu? mkarr(false, basis.num_basis, localRange_ext.volume)
                                       : gkyl_array_acquire(phibc);
  gkyl_eval_on_nodes *projob_bc = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc2x_dirichletvarx_dirichletvary_bc, NULL);
  gkyl_eval_on_nodes_advance(projob_bc, 0.0, &localRange, phibc_ho);
  gkyl_array_copy(phibc, phibc_ho);

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bcs, NULL, epsilon, NULL, true, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho, phibc);
  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);
  gkyl_array_copy(phi_ho, phi);
//  gkyl_grid_sub_array_write(&grid, &localRange, NULL, phi_ho, "ctest_fem_poisson_2x_phi_1.gkyl");

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    cudaDeviceSynchronize();
  }
#endif

  if (poly_order == 1) {
    // Solution with N=8x8 (checked visually against published Fenics result):
    const double sol[256] = {
      2.0468750000000817e+00,  9.0210979560393176e-03, 1.8042195912143118e-02,  3.9857006584043132e-14, 
      2.1718749999999902e+00,  9.0210979561132863e-03, 5.4126587736507296e-02,  2.8865798640254078e-15, 
      2.4218750000000453e+00,  9.0210979560152154e-03, 9.0210979560931215e-02, -5.9581968988216755e-14, 
      2.7968749999997073e+00,  9.0210979561473875e-03, 1.2629537138498323e-01,  1.3589129821411919e-13, 
      3.2968749999995328e+00,  9.0210979562382800e-03, 1.6237976320972894e-01, -8.3414756583503458e-14, 
      3.9218749999999645e+00,  9.0210979560434185e-03, 1.9846415503403644e-01, -2.9124850679333284e-14, 
      4.6718750000000808e+00,  9.0210979560430352e-03, 2.3454854685825011e-01,  2.8717768903637391e-14, 
      5.5468749999999662e+00,  9.0210979560648285e-03, 2.7063293868260607e-01, -1.6135241291218946e-14, 
      
      2.1093750000000862e+00,  2.7063293868315248e-02, 1.8042195912137220e-02, -4.3150668223764432e-14, 
      2.2343749999999996e+00,  2.7063293868243968e-02, 5.4126587736515887e-02,  1.9613940101711104e-15, 
      2.4843749999999245e+00,  2.7063293868266659e-02, 9.0210979560847379e-02,  1.1176245114559913e-14, 
      2.8593749999998526e+00,  2.7063293868288197e-02, 1.2629537138522079e-01,  1.2952601953960163e-15, 
      3.3593749999998450e+00,  2.7063293868293582e-02, 1.6237976320958766e-01,  1.7023419710919072e-15, 
      3.9843749999999138e+00,  2.7063293868278712e-02, 1.9846415503396825e-01, -1.0288066694859786e-14, 
      4.7343749999999991e+00,  2.7063293868261659e-02, 2.3454854685830037e-01,  3.7007434154171896e-16, 
      5.6093749999999956e+00,  2.7063293868304222e-02, 2.7063293868262017e-01,  2.4054832200211730e-14, 
      
      2.2343749999999809e+00,  4.5105489780326961e-02, 1.8042195912166322e-02,  5.9915035895604302e-14, 
      2.3593749999999547e+00,  4.5105489780433362e-02, 5.4126587736521785e-02,  1.5173048003210476e-15, 
      2.6093749999999329e+00,  4.5105489780441314e-02, 9.0210979560872109e-02,  3.1086244689504391e-15, 
      2.9843749999999147e+00,  4.5105489780450675e-02, 1.2629537138522706e-01,  2.2944609175586574e-15, 
      3.4843749999999201e+00,  4.5105489780452854e-02, 1.6237976320958894e-01, -9.2518585385429738e-16, 
      4.1093749999999556e+00,  4.5105489780448746e-02, 1.9846415503394774e-01, -1.4802973661668758e-15, 
      4.8593750000000107e+00,  4.5105489780448496e-02, 2.3454854685830370e-01,  1.4062824978585319e-15, 
      5.7343750000001048e+00,  4.5105489780462082e-02, 2.7063293868267274e-01,  6.3652786745175660e-15, 
      
      2.4218749999998392e+00,  6.3147685692646088e-02, 1.8042195912224907e-02, -2.6090241078691185e-14, 
      2.5468749999999325e+00,  6.3147685692608410e-02, 5.4126587736532165e-02,  4.3298697960381121e-15, 
      2.7968749999999414e+00,  6.3147685692618152e-02, 9.0210979560879423e-02,  1.2212453270876726e-15, 
      3.1718749999999440e+00,  6.3147685692621358e-02, 1.2629537138523206e-01,  6.2912638062092217e-16, 
      3.6718749999999547e+00,  6.3147685692621733e-02, 1.6237976320958675e-01, -3.7007434154171896e-16, 
      4.2968749999999787e+00,  6.3147685692619429e-02, 1.9846415503394324e-01, -1.1842378929335006e-15, 
      5.0468750000000222e+00,  6.3147685692612768e-02, 2.3454854685830140e-01, -2.7385501274087204e-15, 
      5.9218750000001190e+00,  6.3147685692601485e-02, 2.7063293868267735e-01, -3.4786988104921582e-15, 
      
      2.6718749999999214e+00,  8.1189881604807806e-02, 1.8042195912175037e-02, -2.7755575615628921e-15, 
      2.7968749999999338e+00,  8.1189881604799077e-02, 5.4126587736535371e-02, -2.3684757858670013e-15, 
      3.0468749999999525e+00,  8.1189881604794983e-02, 9.0210979560881852e-02,  1.4802973661668758e-16, 
      3.4218749999999614e+00,  8.1189881604794983e-02, 1.2629537138523297e-01, -1.1102230246251568e-16, 
      3.9218749999999711e+00,  8.1189881604794470e-02, 1.6237976320958561e-01, -4.0708177569589086e-16, 
      4.5468749999999858e+00,  8.1189881604791653e-02, 1.9846415503393927e-01, -1.1102230246251569e-15, 
      5.2968750000000071e+00,  8.1189881604784978e-02, 2.3454854685829216e-01, -2.4424906541753452e-15, 
      6.1718750000000160e+00,  8.1189881604744732e-02, 2.7063293868263533e-01, -2.0872192862952949e-14, 
      
      2.9843749999999503e+00,  9.9232077516967068e-02, 1.8042195912180037e-02,  5.6991448597424721e-15, 
      3.1093749999999627e+00,  9.9232077516975783e-02, 5.4126587736530368e-02, -5.9211894646675032e-16, 
      3.3593749999999707e+00,  9.9232077516973480e-02, 9.0210979560880575e-02, -7.7715611723760978e-16, 
      3.7343749999999751e+00,  9.9232077516970912e-02, 1.2629537138523167e-01, -6.6613381477509412e-16, 
      4.2343749999999805e+00,  9.9232077516968734e-02, 1.6237976320958405e-01, -4.8109664400423467e-16, 
      4.8593749999999876e+00,  9.9232077516966943e-02, 1.9846415503393619e-01, -6.6613381477509412e-16, 
      5.6093749999999911e+00,  9.9232077516963613e-02, 2.3454854685828549e-01, -1.3322676295501882e-15, 
      6.4843749999999627e+00,  9.9232077516982584e-02, 2.7063293868262045e-01,  1.2286468139185069e-14, 
      
      3.3593749999999796e+00,  1.1727427342915954e-01, 1.8042195912185421e-02, -2.5535129566378608e-15, 
      3.4843749999999991e+00,  1.1727427342915479e-01, 5.4126587736528960e-02, -2.5905203907920326e-16, 
      3.7343749999999947e+00,  1.1727427342915031e-01, 9.0210979560875065e-02, -2.3684757858670013e-15, 
      4.1093749999999867e+00,  1.1727427342914544e-01, 1.2629537138522975e-01, -3.7007434154171896e-16, 
      4.6093749999999858e+00,  1.1727427342914454e-01, 1.6237976320958303e-01, -2.2204460492503136e-16, 
      5.2343749999999893e+00,  1.1727427342914377e-01, 1.9846415503393491e-01, -1.4802973661668758e-16, 
      5.9843749999999893e+00,  1.1727427342914505e-01, 2.3454854685828447e-01,  9.6219328800846933e-16, 
      6.8593749999999929e+00,  1.1727427342914480e-01, 2.7063293868264021e-01, -8.8817841970012543e-16, 
      
      3.7968750000000613e+00,  1.3531646934134881e-01, 1.8042195912204909e-02,  1.3692750637043602e-14, 
      3.9218750000000648e+00,  1.3531646934134434e-01, 5.4126587736500247e-02, -1.6320278461989804e-14, 
      4.1718750000000089e+00,  1.3531646934131908e-01, 9.0210979560873913e-02,  1.7393494052460791e-15, 
      4.5468749999999911e+00,  1.3531646934131880e-01, 1.2629537138522579e-01, -1.9243865760169387e-15, 
      5.0468749999999893e+00,  1.3531646934131855e-01, 1.6237976320958561e-01,  1.7763568394002509e-15, 
      5.6718749999999956e+00,  1.3531646934132113e-01, 1.9846415503393439e-01, -2.9605947323337516e-16, 
      6.4218749999999964e+00,  1.3531646934132061e-01, 2.3454854685828549e-01, -7.4014868308343790e-17, 
      7.2968749999999982e+00,  1.3531646934131958e-01, 2.7063293868263816e-01, -5.1810407815840652e-16, 
    };
    for (int j=0; j<cells[0]; j++) {
      for (int k=0; k<cells[1]; k++) {
        int idx0[] = {j+1,k+1};
        long linidx = gkyl_range_idx(&localRange, idx0);
        const double *phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    }
  } if (poly_order == 2) {
    // Solution with N=8x8 (checked visually against published Fenics result):
    const double sol[512] =  {
      2.0312500000001257e+00,  9.0210979560507252e-03,  1.8042195912104955e-02, -7.0647191800314145e-14, 
      2.3292374766316335e-03,  4.6584749531214647e-03, -8.1947194606234097e-14,  6.7645666878016016e-15, 
      2.1562499999998685e+00,  9.0210979561283280e-03,  5.4126587736559463e-02,  2.7940612786399781e-14, 
      2.3292374765457422e-03,  4.6584749532062484e-03,  3.2332739093037856e-14, -6.1077809810623154e-14, 
      2.4062500000003104e+00,  9.0210979559354754e-03,  9.0210979560834750e-02,  6.1469348130079521e-14, 
      2.3292374765589924e-03,  4.6584749529991823e-03, -2.4694159303746268e-14,  9.4367377650456223e-14, 
      2.7812500000000364e+00,  9.0210979561349929e-03,  1.2629537138530342e-01, -7.3533771664339552e-14, 
      2.3292374765769685e-03,  4.6584749531438747e-03,  3.5039438981605374e-14, -4.2046585965335703e-15, 
      3.2812500000002536e+00,  9.0210979559721821e-03,  1.6237976320941239e-01,  1.4236759919109928e-13, 
      2.3292374765516732e-03,  4.6584749529714631e-03, -4.9658348687980136e-14,  1.2192473212423358e-13, 
      3.9062499999998663e+00,  9.0210979562064860e-03,  1.9846415503412188e-01, -1.8303876932653420e-13, 
      2.3292374765487150e-03,  4.6584749531586927e-03,  4.7994499732987689e-14, -1.4333361791466009e-14, 
      4.6562500000001066e+00,  9.0210979560747407e-03,  2.3454854685823862e-01,  6.7131485555667820e-14, 
      2.3292374765842734e-03,  4.6584749531606165e-03, -2.7452598419239407e-14, -4.5254386700025249e-14, 
      5.5312500000000648e+00,  9.0210979561203795e-03,  2.7063293868261540e-01,  1.4802973661668758e-14, 
      2.3292374765635925e-03,  4.6584749531257815e-03,  1.5434742833184985e-14, -2.1410275698615774e-15, 
     
      2.0937499999995013e+00,  2.7063293868232811e-02,  1.8042195912383999e-02,  2.0465111087257056e-14, 
      2.3292374768579061e-03,  4.6584749531138068e-03, -2.4572480377165287e-13, -1.1238466779958611e-14, 
      2.2187499999999361e+00,  2.7063293868219863e-02,  5.4126587736566770e-02, -5.4030853865090965e-15, 
      2.3292374765130379e-03,  4.6584749531113261e-03,  4.6550895465749560e-14,  6.1908930800159793e-15, 
      2.4687499999999840e+00,  2.7063293868249391e-02,  9.0210979560905138e-02, -2.4794980883295169e-15, 
      2.3292374765755295e-03,  4.6584749531401936e-03, -1.0476388364290738e-14, -1.3035260273703731e-14, 
      2.8437500000000799e+00,  2.7063293868223197e-02,  1.2629537138524213e-01, -3.2196467714129548e-15, 
      2.3292374765624463e-03,  4.6584749531265136e-03,  2.9375630344784424e-15, -5.7340928323701605e-15, 
      3.3437500000000315e+00,  2.7063293868268494e-02,  1.6237976320957126e-01,  1.1361282285330771e-14, 
      2.3292374765646966e-03,  4.6584749531484067e-03, -1.5770822519325964e-15, -1.9685886830975228e-14, 
      3.9687500000000635e+00,  2.7063293868265288e-02,  1.9846415503393694e-01,  8.2156503822261600e-15, 
      2.3292374765540333e-03,  4.6584749531285623e-03, -4.7091552447998268e-15, -3.0592471352677708e-15, 
      4.7187500000001839e+00,  2.7063293868301995e-02,  2.3454854685834495e-01,  4.7147471112414996e-14, 
      2.3292374765693504e-03,  4.6584749531227406e-03,  1.3515389699622772e-14,  2.3466733138357839e-14, 
      5.5937500000001084e+00,  2.7063293868275416e-02,  2.7063293868255051e-01, -8.2896652505345043e-14, 
      2.3292374765787626e-03,  4.6584749531356972e-03, -8.1958484395797536e-15,  7.7194477785636871e-15, 
     
      2.2187499999999125e+00,  4.5105489780491989e-02,  1.8042195912075981e-02, -5.7842619582970668e-14, 
      2.3292374766908032e-03,  4.6584749531403532e-03, -1.3687086767597140e-13,  2.6637397095815581e-14, 
      2.3437499999998317e+00,  4.5105489780433487e-02,  5.4126587736560317e-02, -8.1786429480719886e-15, 
      2.3292374765212622e-03,  4.6584749531249662e-03,  3.8944952592307181e-14,  1.7178866998293063e-15, 
      2.5937499999999329e+00,  4.5105489780421053e-02,  9.0210979560899962e-02,  1.5173048003210476e-15, 
      2.3292374765734123e-03,  4.6584749531210961e-03, -8.9085791416444597e-15,  2.0618774948951493e-15, 
      2.9687499999999907e+00,  4.5105489780426777e-02,  1.2629537138524269e-01,  2.1094237467877978e-15, 
      2.3292374765610684e-03,  4.6584749531206832e-03,  1.8265392494382348e-15,  2.3682559923738673e-15, 
      3.4687500000000244e+00,  4.5105489780429213e-02,  1.6237976320958916e-01,  4.8109664400423467e-16, 
      2.3292374765635556e-03,  4.6584749531200448e-03, -3.5280701924203503e-16,  3.2852955422680583e-15, 
      4.0937500000000577e+00,  4.5105489780437803e-02,  1.9846415503395051e-01,  8.1416355139178173e-16, 
      2.3292374765561388e-03,  4.6584749531240442e-03, -3.9064366250724966e-15,  5.3228861362296299e-16, 
      4.8437500000000027e+00,  4.5105489780347209e-02,  2.3454854685824650e-01, -7.7715611723760983e-14, 
      2.3292374766086779e-03,  4.6584749531312641e-03,  3.4164249725952891e-14, -1.8577763747153151e-14, 
      5.7187500000011280e+00,  4.5105489780355758e-02,  2.7063293868331445e-01,  8.6597395920762235e-14, 
      2.3292374760671493e-03,  4.6584749531223867e-03, -3.4698372074465863e-13, -1.5290901529256266e-14, 
     
      2.4062499999994902e+00,  6.3147685692434133e-02,  1.8042195912316907e-02,  1.3748261788274859e-13, 
      2.3292374767801831e-03,  4.6584749531267998e-03, -1.8292727109120022e-13, -3.4440423165051116e-14, 
      2.5312499999997988e+00,  6.3147685692606731e-02,  5.4126587736557409e-02,  1.6653345369377352e-15, 
      2.3292374765244337e-03,  4.6584749531210865e-03,  3.5197377912472244e-14, -3.9016739951241898e-15, 
      2.7812499999998876e+00,  6.3147685692605829e-02,  9.0210979560903945e-02,  2.4054832200211730e-15, 
      2.3292374765721342e-03,  4.6584749531240477e-03, -7.6468443945834170e-15, -3.4678682376854966e-16, 
      3.1562499999999560e+00,  6.3147685692608743e-02,  1.2629537138524649e-01, -2.9605947323337516e-16, 
      2.3292374765614549e-03,  4.6584749531248014e-03,  1.5610309715696430e-15, -4.2619435446704727e-17, 
      3.6562500000000049e+00,  6.3147685692611685e-02,  1.6237976320959430e-01,  1.4062824978585319e-15, 
      2.3292374765619770e-03,  4.6584749531249203e-03, -1.2714864788379646e-15, -6.1321025639719978e-16, 
      4.2812500000000160e+00,  6.3147685692615488e-02,  1.9846415503392814e-01,  1.5543122344752196e-15, 
      2.3292374765739154e-03,  4.6584749531251197e-03,  7.9402782950453438e-15, -8.0345355108934640e-17, 
      5.0312500000001226e+00,  6.3147685692635055e-02,  2.3454854685832488e-01,  1.2730557349035132e-14, 
      2.3292374764988227e-03,  4.6584749531028034e-03, -5.1233949082405346e-14,  2.2886757802088086e-15, 
      5.9062499999982467e+00,  6.3147685692658814e-02,  2.7063293868152366e-01,  8.2156503822261600e-15, 
      2.3292374773235388e-03,  4.6584749531093901e-03,  5.2751130719945500e-13,  7.9465038372885637e-15, 
     
      2.6562499999996114e+00,  8.1189881604835173e-02,  1.8042195912257340e-02,  3.2936616397212987e-15, 
      2.3292374766196522e-03,  4.6584749531023021e-03, -4.7241819098139020e-14,  2.0252451952469991e-14, 
      2.7812499999997815e+00,  8.1189881604829109e-02,  5.4126587736569504e-02, -2.7311486405778857e-14, 
      2.3292374765548894e-03,  4.6584749531218602e-03,  9.8975008642476021e-15,  4.3918744393441411e-15, 
      3.0312499999998823e+00,  8.1189881604790695e-02,  9.0210979560897020e-02, -8.5117098554595358e-16, 
      2.3292374765662223e-03,  4.6584749531230242e-03, -3.2908957769632486e-15, -2.7012598421088404e-16, 
      3.4062499999999409e+00,  8.1189881604788641e-02,  1.2629537138524627e-01, -1.1102230246251569e-15, 
      2.3292374765615693e-03,  4.6584749531232436e-03,  4.5035036956445722e-16, -8.4542500933764381e-16, 
      3.9062499999999920e+00,  8.1189881604789502e-02,  1.6237976320959771e-01,  3.1456319031046112e-15, 
      2.3292374765636966e-03,  4.6584749531243452e-03,  7.5762521350236672e-16,  3.7951855339199665e-16, 
      4.5312500000000648e+00,  8.1189881604792485e-02,  1.9846415503396136e-01, -1.5543122344752196e-15, 
      2.3292374765527969e-03,  4.6584749531248318e-03, -7.0422838590550659e-15,  7.2720095134983539e-17, 
      5.2812500000000648e+00,  8.1189881604893516e-02,  2.3454854685827042e-01,  8.0232117246244667e-14, 
      2.3292374766189631e-03,  4.6584749531346364e-03,  4.5240515544967410e-14,  1.5894383950387577e-14, 
      6.1562500000014539e+00,  8.1189881604849495e-02,  2.7063293868345234e-01, -1.1886787850320013e-13, 
      2.3292374759688031e-03,  4.6584749531335930e-03, -4.2066895533795610e-13,  5.6532268422286223e-15, 
     
      2.9687499999999787e+00,  9.9232077517012532e-02,  1.8042195912102432e-02, -3.1826393372587830e-15, 
      2.3292374765248977e-03,  4.6584749531297012e-03,  2.2091223854892527e-14, -4.3996766967249189e-15, 
      3.0937499999998770e+00,  9.9232077516994241e-02,  5.4126587736536516e-02, -4.8109664400423461e-15, 
      2.3292374765628123e-03,  4.6584749531253236e-03, -2.3055680720774022e-16, -2.4129403865107740e-15, 
      3.3437499999999045e+00,  9.9232077516974146e-02,  9.0210979560885779e-02, -1.7763568394002509e-15, 
      2.3292374765619714e-03,  4.6584749531253695e-03, -2.7289534500209458e-16,  1.5253982230489157e-15, 
      3.7187499999999174e+00,  9.9232077516956160e-02,  1.2629537138524199e-01, -1.9984014443252822e-15, 
      2.3292374765613894e-03,  4.6584749531331940e-03, -4.4121826562303549e-17,  6.6464750601897660e-15, 
      4.2187499999999760e+00,  9.9232077516958631e-02,  1.6237976320959840e-01, -2.7755575615628921e-15, 
      2.3292374765628851e-03,  4.6584749531283593e-03,  7.9550314411133817e-16,  1.7544780801918967e-15, 
      4.8437500000000622e+00,  9.9232077516965222e-02,  1.9846415503396325e-01,  3.7007434154171895e-15, 
      2.3292374765532791e-03,  4.6584749531246574e-03, -6.2794631070053001e-15, -3.8670504428678561e-16, 
      5.5937500000002194e+00,  9.9232077516901204e-02,  2.3454854685836093e-01, -5.9211894646675032e-14, 
      2.3292374765787041e-03,  4.6584749531364223e-03,  2.0936574346453288e-14, -1.4756579258474011e-14, 
      6.4687500000001918e+00,  9.9232077516931375e-02,  2.7063293868253890e-01,  8.6671410789070575e-14, 
      2.3292374765506380e-03,  4.6584749531308678e-03, -3.7161118642530208e-14, -7.0362605745019720e-15, 
     
      3.3437499999999329e+00,  1.1727427342916795e-01,  1.8042195912195375e-02, -1.3655743202889429e-14, 
      2.3292374765998083e-03,  4.6584749531311514e-03, -3.2488553096641561e-14,  5.2323965287468765e-15, 
      3.4687499999999565e+00,  1.1727427342915786e-01,  5.4126587736523950e-02,  1.0029014655780583e-14, 
      2.3292374765598945e-03,  4.6584749531333215e-03,  9.4006294759552971e-15,  6.9863743267468417e-15, 
      3.7187499999998987e+00,  1.1727427342913505e-01,  9.0210979560836374e-02, -3.5120055012309129e-14, 
      2.3292374765645054e-03,  4.6584749531238872e-03, -6.7315349345457895e-15, -2.2957870379001291e-15, 
      4.0937499999999032e+00,  1.1727427342913478e-01,  1.2629537138523184e-01, -3.2936616397212987e-15, 
      2.3292374765535831e-03,  4.6584749530893133e-03,  3.7348972603798000e-16, -3.1998544279100648e-14, 
      4.5937499999999192e+00,  1.1727427342911641e-01,  1.6237976320963080e-01,  2.6275278249462047e-14, 
      2.3292374765617979e-03,  4.6584749531210596e-03,  4.4292513154912786e-15, -5.8905524963625637e-15, 
      5.2187500000000204e+00,  1.1727427342913538e-01,  1.9846415503395351e-01, -4.3668772301922834e-15, 
      2.3292374765652053e-03,  4.6584749531281998e-03, -2.3804519379379431e-15,  2.5914821573507062e-15, 
      5.9687500000000897e+00,  1.1727427342911410e-01,  2.3454854685829205e-01, -6.3652786745175660e-15, 
      2.3292374765628131e-03,  4.6584749531173291e-03,  1.0566464008386871e-15,  3.7395274599521350e-15, 
      6.8437500000000435e+00,  1.1727427342911957e-01,  2.7063293868260813e-01,  5.9952043329758469e-15, 
      2.3292374765711311e-03,  4.6584749531210371e-03,  3.6568088083192496e-15,  1.0642782935626107e-15, 
     
      3.7812499999999885e+00,  1.3531646934129221e-01,  1.8042195912098033e-02, -4.4927025063164683e-14, 
      2.3292374765743404e-03,  4.6584749531190968e-03, -3.4362057676687731e-14, -1.2196816090811993e-14, 
      3.9062500000000369e+00,  1.3531646934134264e-01,  5.4126587736587753e-02,  4.4075854077618727e-14, 
      2.3292374765543945e-03,  4.6584749530839669e-03,  2.2777841481676262e-14, -3.5549115766876916e-14, 
      4.1562500000001332e+00,  1.3531646934146771e-01,  9.0210979560846463e-02,  3.0716170347962673e-14, 
      2.3292374765685924e-03,  4.6584749530620860e-03, -1.4681564189098153e-14, -3.3529181380867156e-14, 
      4.5312499999993356e+00,  1.3531646934100616e-01,  1.2629537138525038e-01,  2.5905203907920327e-14, 
      2.3292374765597128e-03,  4.6584749534094626e-03,  9.6263135748099880e-15,  2.1685856667999877e-13, 
      5.0312500000000480e+00,  1.3531646934142214e-01,  1.6237976320956182e-01, -8.1194310534253132e-14, 
      2.3292374765639104e-03,  4.6584749530902015e-03, -7.2725414149003426e-15, -1.1931983865681503e-14, 
      5.6562499999999956e+00,  1.3531646934130417e-01,  1.9846415503395343e-01,  1.3248661427193539e-14, 
      2.3292374765591182e-03,  4.6584749531124311e-03,  4.6521938114153061e-15, -1.1698836831131545e-14, 
      6.4062500000000115e+00,  1.3531646934130512e-01,  2.3454854685828999e-01,  2.2204460492503139e-15, 
      2.3292374765648740e-03,  4.6584749531238256e-03, -1.2334173761407129e-15, -8.4690795679407940e-17, 
      7.2812500000000107e+00,  1.3531646934131161e-01,  2.7063293868263177e-01,  2.8865798640254078e-15, 
      2.3292374765635556e-03,  4.6584749531242212e-03,  6.3984964642330207e-17,  9.1210764523807938e-16, 
    };

    for (int j=0; j<cells[0]; j++) {
      for (int k=0; k<cells[1]; k++) {
        int idx0[] = {j+1,k+1};
        long linidx = gkyl_range_idx(&localRange, idx0);
        const double *phi_p  = gkyl_array_cfetch(phi, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-9) );
          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
          TEST_MSG("Produced: %.13e", phi_p[m]);
        }
      }
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_eval_on_nodes_release(projob_bc);
  gkyl_array_release(phibc);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phibc_ho);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
}

void
test_2x_bias(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  if ( ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1]==GKYL_POISSON_NEUMANN && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_NEUMANN))
      || ((bcs.lo_type[0]==GKYL_POISSON_NEUMANN && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
      || ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_NEUMANN) &&
          (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET))
     ) {
    lower[0] = 0.; lower[1] = 0.;
    upper[0] = 1.; upper[1] = 1.;
  }
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob;
  if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
      (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_dirichlety, NULL);
  } else if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_ho, *phi_ho;
  if (use_gpu) {
    rho_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
    phi_ho = mkarr(false, basis.num_basis, localRange_ext.volume);
  }
  else {
    rho_ho = gkyl_array_acquire(rho);
    phi_ho = gkyl_array_acquire(phi);
  }

  struct gkyl_array *epsilon = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  gkyl_array_clear(epsilon, 0.);
  gkyl_array_shiftc(epsilon, epsilon_0*pow(sqrt(2.),dim), 0);

  // Only used in the fully periodic case.
  struct gkyl_array *rho_cellavg = mkarr(use_gpu, 1, localRange_ext.volume);

  // Project the right-side source on the basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);
  gkyl_array_copy(rho, rho_ho);

  struct gkyl_array *perbuff = mkarr(use_gpu, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);

  gkyl_grid_sub_array_write(&grid, &localRange, NULL, rho_ho, "ctest_fem_poisson_2x_rho_1.gkyl");

//  // Specify a bias at a fixed x:
//  struct gkyl_poisson_bias_plane bias = {
//    .dir = 0,
//    .loc = -M_PI+2*2*M_PI/8, // Location of the plane in the 'dir' dimension.
//    .val = 0.2, // Biasing value.
//  };
//  struct gkyl_poisson_bias_plane_list bpl = {
//    .num_bias_plane = 1, // Number of bias planes.
//    .bp = &bias,
//  };
//  // Specify a bias at a fixed y:
//  struct gkyl_poisson_bias_plane bias = {
//    .dir = 1,
//    .loc = -1.0, // Location of the plane in the 'dir' dimension.
//    .val = 0., // Biasing value.
//  };
//  struct gkyl_poisson_bias_plane_list bpl = {
//    .num_bias_plane = 1, // Number of bias planes.
//    .bp = &bias,
//  };
//  // Specify a bias at a fixed x and another at a fixed y:
  struct gkyl_poisson_bias_plane bias[2];
  bias[0].dir = 0;
  bias[0].loc = -M_PI+2*2*M_PI/8; // Location of the plane in the 'dir' dimension.
  bias[0].val = 0.3; // Biasing value.
  bias[1].dir = 1;
  bias[1].loc = 0.8; // Location of the plane in the 'dir' dimension.
  bias[1].val = 0.3; // Biasing value.
  struct gkyl_poisson_bias_plane_list bpl = {
    .num_bias_plane = 2, // Number of bias planes.
    .bp = bias,
  };
//  struct gkyl_poisson_bias_plane_list bpl = {
//    .num_bias_plane = 0, // Number of bias planes.
//  };

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bcs, &bpl, epsilon, NULL, true, use_gpu);

  // Set the RHS source.
  gkyl_fem_poisson_set_rhs(poisson, rho, NULL);
  // Solve the problem.
  gkyl_fem_poisson_solve(poisson, phi);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    cudaDeviceSynchronize();
  }
#endif

  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);

  gkyl_grid_sub_array_write(&grid, &localRange, NULL, phi_ho, "ctest_fem_poisson_2x_phi_1.gkyl");

  gkyl_array_copy(phi_ho, phi);

  if (poly_order == 1) {
    if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8:
      const double sol[256] = {
         1.2895156953678269e-01,  7.4450223384485620e-02,  7.4450223384485828e-02,  4.2983856512260797e-02,
         2.8830517927384303e-01,  1.6645307286251634e-01,  1.7552626093544542e-02,  1.0134013400092890e-02,
         3.3225367768176517e-01,  1.9182675024880982e-01,  7.8210512927491659e-03,  4.5154860692143634e-03,
         3.5214664660491579e-01,  2.0331196121157141e-01,  3.6641596700126913e-03,  2.1155035718356517e-03,
         3.6972744538778313e-01,  2.1346224012142875e-01,  6.4861192398443946e-03,  3.7447626891202847e-03,
         4.9048086672757174e-01,  1.0997417968063671e-01,  6.3230901076250887e-02, -6.3493622242523848e-02,
         4.4653236831964949e-01,  8.4600502294343086e-02, -8.8604578462544659e-02,  4.8844122773216501e-02,
         1.4653236831964953e-01,  8.4600502294343058e-02, -8.4600502294343058e-02, -4.8844122773216522e-02,
                                
         4.2895156953678315e-01,  9.8754857372402388e-02,  7.4450223384485453e-02, -4.2983856512261012e-02,
         5.8830517927384274e-01,  6.7520078943712021e-03,  1.7552626093544479e-02, -1.0134013400092921e-02,
         6.3225367768176433e-01, -1.8621669491922586e-02,  7.8210512927488936e-03, -4.5154860692145091e-03,
         6.5214664660491373e-01, -3.0106880454684778e-02,  3.6641596700122750e-03, -2.1155035718358941e-03,
         6.6972744538778173e-01, -4.0257159364541752e-02,  6.4861192398452438e-03, -3.7447626891198025e-03,
         6.4048086672757265e-01, -2.3371639302192299e-02, -2.3371639302192552e-02,  1.3493622242524079e-02,
         5.9653236831965029e-01,  2.0020380841012666e-03, -2.0020380841012666e-03,  1.1558772267832309e-03,
         4.4653236831964949e-01,  8.8604578462544686e-02, -8.4600502294343058e-02,  4.8844122773216522e-02,
                                
         6.1722373660761676e-01,  9.9441289668576079e-03,  1.8314920972374435e-01,  1.0574124553587159e-01,
         1.0938631994762913e+00,  2.8513205116048679e-01,  9.2038712469883516e-02,  5.3138575420354009e-02,
         1.3261247268759631e+00,  4.1922830652708126e-01,  4.2057542896710598e-02,  2.4281933712870175e-02,
         1.3866182591485865e+00,  4.5415426366891287e-01, -7.1315857548800068e-03, -4.1174229553283363e-03,
         1.2598753892685446e+00,  3.8097923362580427e-01, -6.6043444288228045e-02, -3.8130200338019393e-02,
         8.7274239412724519e-01,  1.5746789466878686e-01, -1.5746789466878724e-01, -9.0914131375748317e-02,
         6.4048086672757232e-01,  2.3371639302192361e-02,  2.3371639302192489e-02,  1.3493622242524190e-02,
         4.9048086672757152e-01, -6.3230901076251081e-02, -1.0997417968063683e-01, -6.3493622242523959e-02,
                                
         7.0276798098732818e-01,  3.9444863553391919e-02,  4.0574328300088197e-01,  2.2773502590699067e-02,
         1.8104132743081691e+00,  1.2856832743154439e-01,  2.3375602526452760e-01,  2.8681953270464161e-02,
         2.4005612181889311e+00,  2.0109789095961778e-01,  1.0696604899673452e-01,  1.3193009756675185e-02,
         2.5504406343836705e+00,  2.1777889796197247e-01, -2.0433127710717744e-02, -3.5622258735116272e-03,
         2.2174000987841271e+00,  1.7184791516873588e-01, -1.7184791516873568e-01, -2.2956039406307645e-02,
         1.2598753892685450e+00,  6.6043444288228809e-02, -3.8097923362580405e-01, -3.8130200338019032e-02,
         6.6972744538778239e-01, -6.4861192398445074e-03,  4.0257159364541495e-02, -3.7447626891201633e-03,
         3.6972744538778263e-01, -6.4861192398443304e-03, -2.1346224012142911e-01,  3.7447626891202652e-03,
                                
         7.7717226524603089e-01,  3.5124699922315329e-03,  4.4870061654650595e-01,  2.0279254955354741e-03,
         2.0554433053316998e+00,  1.2899826939433492e-02,  2.8930951254425802e-01,  3.3918675649106675e-03,
         2.7899149178753726e+00,  2.3695572364325838e-02,  1.3473787066997003e-01,  2.8410589622536184e-03,
         2.9754662530247851e+00,  2.7609757372962000e-02, -2.7609757372962257e-02, -5.8120319385953278e-04,
         2.5504406343836705e+00,  2.0433127710717487e-02, -2.1777889796197247e-01, -3.5622258735116272e-03,
         1.3866182591485874e+00,  7.1315857548793658e-03, -4.5415426366891248e-01, -4.1174229553286512e-03,
         6.5214664660491461e-01, -3.6641596700128197e-03,  3.0106880454684556e-02, -2.1155035718355242e-03,
         3.5214664660491596e-01, -3.6641596700124029e-03, -2.0331196121157155e-01,  2.1155035718357644e-03,
                                
         7.2160737736820613e-01, -3.5592872965984977e-02,  4.1662021357275375e-01, -2.0549554788142962e-02,
         1.9135017578046214e+00, -9.4849817613384166e-02,  2.7152032781781293e-01, -1.3662458155387695e-02,
         2.6073728069988200e+00, -1.2908630921734573e-01,  1.2908630921734573e-01, -6.1039894882679970e-03,
         2.7899149178753717e+00, -1.3473787066997003e-01, -2.3695572364325967e-02,  2.8410589622536553e-03,
         2.4005612181889306e+00, -1.0696604899673465e-01, -2.0109789095961778e-01,  1.3193009756675074e-02,
         1.3261247268759631e+00, -4.2057542896710404e-02, -4.1922830652708110e-01,  2.4281933712870359e-02,
         6.3225367768176466e-01, -7.8210512927487010e-03,  1.8621669491922586e-02, -4.5154860692147311e-03,
         3.3225367768176545e-01, -7.8210512927493741e-03, -1.9182675024880988e-01,  4.5154860692143339e-03,
                                
         5.3333521029737307e-01, -7.3106113373274759e-02,  3.0792122723349458e-01, -4.2207834235467784e-02,
         1.4079437376021731e+00, -1.9703424144147386e-01,  1.9703424144147386e-01, -2.9342103864873393e-02,
         1.9135017578046212e+00, -2.7152032781781293e-01,  9.4849817613383916e-02, -1.3662458155387658e-02,
         2.0554433053316989e+00, -2.8930951254425813e-01, -1.2899826939433621e-02,  3.3918675649105383e-03,
         1.8104132743081678e+00, -2.3375602526452785e-01, -1.2856832743154434e-01,  2.8681953270464161e-02,
         1.0938631994762908e+00, -9.2038712469883932e-02, -2.8513205116048629e-01,  5.3138575420353953e-02,
         5.8830517927384285e-01, -1.7552626093544799e-02, -6.7520078943714267e-03, -1.0134013400092848e-02,
         2.8830517927384264e-01, -1.7552626093544701e-02, -1.6645307286251645e-01,  1.0134013400092904e-02,
                                
         2.0335585379548557e-01, -1.1740755693010950e-01,  1.1740755693010979e-01, -6.7785284598495046e-02,
         5.3333521029737341e-01, -3.0792122723349435e-01,  7.3106113373274675e-02, -4.2207834235467825e-02,
         7.2160737736820590e-01, -4.1662021357275386e-01,  3.5592872965984651e-02, -2.0549554788143062e-02,
         7.7717226524602934e-01, -4.4870061654650611e-01, -3.5124699922319496e-03,  2.0279254955354879e-03,
         7.0276798098732574e-01, -4.0574328300088236e-01, -3.9444863553392016e-02,  2.2773502590698942e-02,
         6.1722373660761443e-01, -1.8314920972374513e-01, -9.9441289668574153e-03,  1.0574124553587151e-01,
         4.2895156953678204e-01, -7.4450223384485856e-02, -9.8754857372401889e-02, -4.2983856512260686e-02,
         1.2895156953678202e-01, -7.4450223384485856e-02, -7.4450223384485856e-02,  4.2983856512260686e-02,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[256] = {
         2.0695556873727669e-01,  1.1948585332075820e-01, -1.9046453463662513e-02, -1.0996475034353147e-02, 
         1.4558864481965927e-01,  8.4055643277583011e-02, -1.6383756579513183e-02, -9.4591662715190448e-03, 
         1.2818040369987443e-01,  7.4004990580957775e-02,  6.3331038828874937e-03,  3.6564192315910174e-03, 
         1.8756046745906696e-01,  1.0828808637682398e-01,  2.7949991912979449e-02,  1.6136935354806071e-02, 
         2.0345037480984871e-01,  1.1746212866319607e-01, -1.8775949626607395e-02, -1.0840299571212131e-02, 
         3.8546473804810522e-01,  4.9343089518294872e-02,  1.2386199123859315e-01, -2.8488246016035241e-02, 
         4.7796329610410232e-01,  1.0274715691157675e-01, -7.0457923845310980e-02,  5.9321098701367460e-02, 
         2.9793579302427070e-01,  1.7201331030378694e-01, -3.3481003519366004e-02, -1.9330266394644992e-02, 
                                
         5.0695556873727465e-01,  5.3719227436128343e-02, -1.9046453463662225e-02,  1.0996475034353324e-02, 
         4.4558864481965854e-01,  8.9149437479304330e-02, -1.6383756579512680e-02,  9.4591662715193241e-03, 
         4.2818040369987509e-01,  9.9200090175930317e-02,  6.3331038828877582e-03, -3.6564192315908643e-03, 
         4.8756046745906639e-01,  6.4916994380063445e-02,  2.7949991912978502e-02, -1.6136935354806602e-02, 
         5.0345037480984711e-01,  5.5742952093690824e-02, -1.8775949626607009e-02,  1.0840299571212351e-02, 
         5.3546473804810524e-01,  3.7259450860149000e-02,  3.7259450860149736e-02, -2.1511753983964516e-02, 
         6.2796329610410262e-01, -1.6144616533132653e-02,  1.6144616533132653e-02, -9.3210987013675681e-03, 
         5.9793579302426936e-01,  1.1917704531000553e-03, -3.3481003519366739e-02,  1.9330266394644548e-02, 
                                
         7.8902649394019042e-01,  1.0913449716034193e-01, -1.8183182083654607e-02, -1.0498065070722228e-02, 
         7.4401342921382563e-01,  8.3146192123524101e-02, -7.8051229531621259e-03, -4.5062898380666509e-03, 
         7.8907637658376129e-01,  1.0916329691803331e-01,  3.3822227747672386e-02,  1.9527272294711247e-02, 
         9.1813274798934152e-01,  1.8367402768968197e-01,  4.0688503023975697e-02,  2.3491518173815857e-02, 
         8.7470690321784439e-01,  1.5860210452107232e-01, -6.5760426192586471e-02, -3.7966799764314260e-02, 
         6.8040325196245055e-01,  4.6420839164242632e-02, -4.6420839164242340e-02, -2.6801083987484216e-02, 
         7.2274050955105884e-01,  7.0864266229775458e-02,  7.0864266229775874e-02,  4.0913503183686570e-02, 
         8.3300085412723712e-01,  1.3452310585177427e-01, -7.2054266077784040e-03, -4.1600549916263156e-03, 
                                
         1.0664987158368842e+00,  5.1064164844357766e-02, -4.0682412354198595e-02, -2.4918715825357398e-03, 
         9.7165477291101243e-01,  4.8282598938734773e-02, -1.4075763625058893e-02,  8.8593375815563657e-04, 
         1.0940060220721883e+00,  6.6887915988607230e-02,  8.4715290263950777e-02,  9.8558510489465245e-03, 
         1.4003576031526708e+00,  9.4738622248794480e-02,  9.2156877539533236e-02,  6.2237617074937445e-03, 
         1.2974375859789402e+00,  8.5461568965760626e-02, -1.5157777715978749e-01, -1.1579870917739880e-02, 
         8.1744858731991787e-01,  3.2702322095573610e-02, -1.2554400042405800e-01, -1.8880694465004517e-02, 
         8.8097664729592318e-01,  2.0493410492748654e-02,  1.6222194295229936e-01,  1.1831876064601199e-02, 
         1.1494580078003347e+00,  4.8183517008371347e-02, -7.2141571926803982e-03,  4.1550143860829940e-03, 
                                
         1.1491297624369667e+00, -3.3571078463795392e-03, -4.8356641584218818e-02, -1.9388467292392932e-03, 
         1.0325050395128668e+00, -1.3150681135893643e-02, -1.8976673590175835e-02, -3.7154754457900223e-03, 
         1.1511743022204932e+00, -3.3881794055875360e-02,  8.7490404358958399e-02, -8.2536381791619359e-03, 
         1.4622888540779364e+00, -5.8982597855842096e-02,  9.2131665904746168e-02, -6.2383176516247864e-03, 
         1.3496078336183466e+00, -5.5341062447459850e-02, -1.5718808339965204e-01,  8.3407590999178289e-03, 
         8.3867504341288779e-01, -2.0447221956355628e-02, -1.3779910056327618e-01,  1.1805209100681080e-02, 
         9.1819497786477200e-01,  9.9460267996469374e-04,  1.8370995612501295e-01,  5.7423412501523663e-04, 
         1.2346379391368876e+00,  9.9513927834130609e-04, -1.0115272513946342e-03, -5.7392431979802494e-04, 
                                
         1.0495695703772063e+00, -5.4123995839893982e-02, -5.6960458456916822e-02, -3.0285692582711827e-03, 
         8.9218479705791143e-01, -6.7863248618228386e-02, -3.3905682785326366e-02, -4.9037920317645178e-03, 
         9.2801433171934899e-01, -9.4959674985309217e-02,  5.4591874267046404e-02, -1.0740337025346467e-02, 
         1.1455379463444140e+00, -1.2389362399016544e-01,  7.0995443191835192e-02, -5.9646862213262337e-03, 
         1.0657450080545758e+00, -1.0854721630471093e-01, -1.1706391759290431e-01,  1.4824938829617130e-02, 
         7.3149217752530427e-01, -4.1434834516356350e-02, -7.5917044090565114e-02,  2.3922412861846317e-02, 
         8.9385184002337748e-01, -1.5049119865613990e-02,  1.6965543893936455e-01, -8.6886134054795627e-03, 
         1.1679658292468793e+00, -3.9488299870777899e-02, -1.1395653472533528e-02, -5.4213537492755308e-03, 
                                
         7.6749864517429045e-01, -1.0872972875657638e-01, -5.7823729836924409e-02,  2.5301592946400673e-03, 
         5.9376001266374423e-01, -1.0443238098460009e-01, -4.2484316411676902e-02, -4.9084401688122760e-05, 
         5.6711835883546291e-01, -1.1340371210865438e-01,  2.7102750402261823e-02, -5.1305160377738985e-03, 
         7.1496566581413878e-01, -1.2469739807957991e-01,  5.8256932080838000e-02, -1.3898965976830624e-03, 
         6.9448847964657867e-01, -1.0579784031005225e-01, -7.0079441026924913e-02,  1.2301561363484743e-02, 
         5.8655366361095895e-01, -4.2245455508035185e-02,  7.7632459338270975e-03,  2.4390425109602477e-02, 
         7.9907462657642148e-01, -3.9670529831028617e-02,  1.1493578924272123e-01, -2.2903791076839465e-02, 
         9.3290076814391143e-01, -9.6226576434096420e-02, -3.7671230384121929e-02, -9.7488576537427662e-03, 
                                
         2.8958661533735852e-01, -1.6719291031873673e-01, -2.6720682693682486e-02,  1.5427193346128274e-02, 
         2.0643891142151358e-01, -1.1918756108042412e-01, -2.1284666544629894e-02,  1.2288707959153511e-02, 
         1.8534868384817985e-01, -1.0701111251368946e-01,  9.1082179778951937e-03, -5.2586321013756172e-03, 
         2.4949171838433232e-01, -1.4404411076977661e-01,  2.7924780278191910e-02, -1.6122379410675240e-02, 
         2.5562062244925482e-01, -1.4758263518149697e-01, -2.4386255866471529e-02,  1.4079411389034433e-02, 
         4.0669119414107574e-01, -6.1598189657512729e-02,  1.1160689109937504e-01,  3.5563731380358600e-02, 
         5.1518162667295164e-01, -1.2423517008429019e-01, -4.8969910672597566e-02, -7.1727208890983898e-02, 
         3.8311572436082308e-01, -2.2119196659049989e-01, -2.7278373578080651e-02,  1.5749176328359930e-02, 
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    }
  }
  if (poly_order == 2) {
    if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
        1.5533311401374375e-01,  8.0448181189757023e-02,  7.1383887553887471e-02,  3.5882581089513192e-02, 
       -7.1521872221691971e-03, -1.4173358881541273e-02, -4.1293172180145955e-03, -8.1829925655791530e-03, 
        3.1000404198245896e-01,  1.5667545777456296e-01,  3.1082184954682655e-02,  1.5729096122371413e-02, 
       -1.7277734490433301e-02, -3.9744222223506787e-03, -1.7166702230099677e-03, -2.2946337399476009e-03, 
        3.7857548166792060e-01,  1.8900099866514894e-01,  1.2087608801617626e-02,  5.0010090967487408e-03, 
       -2.2904558275500107e-02, -1.2013238572225022e-03, -1.5319780039809404e-03, -6.9358465235151343e-04, 
        4.0425568356501707e-01,  2.0122821599342283e-01,  3.3091265420904133e-03,  2.3876201154691768e-03, 
       -2.4917931882494315e-02, -7.5959936635899025e-04,  3.6955621000353653e-04, -4.3855489864395652e-04, 
        4.1727748169260132e-01,  1.9919888758010218e-01,  3.7258386804603460e-03, -3.8382137320551062e-03, 
       -3.2313368072759904e-02, -1.1338638972530678e-03, -4.6393132852280720e-03, -6.5463662630302550e-04, 
        4.3716603455009728e-01,  1.6861780147080760e-01,  5.1616434224149459e-02, -5.5162092692047124e-02, 
       -2.0174447197502901e-02,  3.2839645130219128e-02,  1.1647722520229270e-02, -3.2679799993408046e-02, 
        3.6859459486463331e-01,  1.3629226058022006e-01, -9.4786227980451360e-02,  3.4431987472926334e-02, 
       -1.4547623412436273e-02,  3.0066546765090796e-02, -8.3990742932384555e-03, -3.4280849081003914e-02, 
        1.6835491214132511e-01,  7.8418852776436732e-02, -7.8418852776436718e-02, -3.4431987472926508e-02, 
       -1.4547623412435470e-02, -1.4547623412435403e-02,  8.3990742932389239e-03, -8.3990742932387886e-03, 
                               
        3.6946581339115087e-01,  9.2756899567128556e-02,  1.3727843226277400e-01, -3.5882581089510680e-02, 
        3.1248837025134639e-02, -1.4173358881541398e-02, -3.3598253481106524e-02,  8.1829925655790611e-03, 
        6.3102752485654046e-01,  1.6529622982323057e-02,  2.6901062871493353e-02, -1.5729096122372332e-02, 
       -2.6679721856484408e-02, -3.9744222223490177e-03,  1.5318441703666890e-04,  2.2946337399485412e-03, 
        6.9398195984660038e-01, -1.5795917908265063e-02,  1.3025751711878395e-02, -5.0010090967484450e-03, 
       -2.9794544775782170e-02, -1.2013238572200467e-03, -1.9515282679712243e-03,  6.9358465235292582e-04, 
        7.2277418406744665e-01, -2.8023135236537666e-02,  4.1677105581693734e-03, -2.3876201154694843e-03, 
       -3.3199657075455628e-02, -7.5959936635750370e-04, -1.4414234866165383e-05,  4.3855489864478870e-04, 
        7.2190296554094235e-01, -2.5993806823209158e-02, -5.1538822407193605e-03,  3.8382137320559979e-03, 
       -3.4381947335499423e-02, -1.1338638972569219e-03, -6.6818136503040088e-04,  6.5463662630079952e-04, 
        6.8178871257800111e-01, -2.4280234173402714e-02, -3.1881507866188909e-02,  5.1620926920462346e-03, 
       -1.7769635704202077e-02, -1.1881714419775227e-02,  1.0259303957221907e-02,  6.8599110186927685e-03, 
        6.1883427758793386e-01,  8.0453067171828759e-03, -8.0453067171869994e-03,  1.5568012527073011e-02, 
       -1.4654812784905439e-02, -1.4654812784904112e-02, -8.4609601062880010e-03,  8.4609601062883271e-03, 
        3.6859459486463103e-01,  9.4786227980450402e-02, -1.3629226058021895e-01,  3.4431987472927521e-02, 
        3.0066546765090994e-02, -1.4547623412436293e-02,  3.4280849081004490e-02,  8.3990742932382664e-03, 
                               
        5.6355951897661249e-01,  1.9416433643301461e-02,  2.4441834236771329e-01,  9.7805239292666477e-02, 
        3.1336583368864267e-02, -1.7984485310210251e-02, -3.3547593105935421e-02, -1.0383347435086744e-02, 
        1.1384058619348891e+00,  2.7612252255954328e-01,  9.9223874625614897e-02,  5.7190569078997951e-02, 
       -2.6898835076978919e-02, -8.8794685107920785e-03, -7.4641410142641562e-05, -5.1265635349686976e-03, 
        1.3850274387252464e+00,  4.1402380087178647e-01,  4.6278677379411756e-02,  2.4225516380595231e-02, 
       -3.0373486296460757e-02, -6.4661984458473437e-03, -1.9314494067652782e-03, -3.7332614133465157e-03, 
        1.4544697477545541e+00,  4.4973913257369985e-01, -6.9271508916094642e-03, -4.0330819406447426e-03, 
       -3.3764054020119380e-02, -7.0401817996435913e-03, -2.6095781194674182e-05, -4.0646508571701365e-03, 
        1.3280762633280476e+00,  3.7870626470624230e-01, -6.8873390005292481e-02, -3.8610064659977676e-02, 
       -3.2260992258763090e-02, -9.2301436869407748e-03,  8.9388889372232631e-04, -5.3290259423093688e-03, 
        9.2841028936836756e-01,  1.6978278617843032e-01, -1.6978278617842624e-01, -8.6578178151638860e-02, 
       -1.5356365639260707e-02, -1.5356365639256775e-02,  8.8660018355975270e-03, -8.8660018356005350e-03, 
        6.8178871257800355e-01,  3.1881507866190303e-02,  2.4280234173395865e-02,  5.1620926920476311e-03, 
       -1.1881714419779777e-02, -1.7769635704201387e-02, -6.8599110186901499e-03, -1.0259303957222634e-02, 
        4.3716603455009401e-01, -5.1616434224148716e-02, -1.6861780147080513e-01, -5.5162092692047818e-02, 
        3.2839645130219648e-02, -2.0174447197503050e-02,  3.2679799993408386e-02, -1.1647722520229139e-02, 
                               
        7.4609938259350195e-01,  3.8300006873649907e-02,  3.9479280796571392e-01,  2.8156146829750164e-02, 
       -5.5907627764121559e-03, -2.7860593399204731e-02, -3.2278283939377078e-03,  4.6813737713599368e-03, 
        1.8293289250498705e+00,  1.3644217999604571e-01,  2.4285963049722542e-01,  2.2005598418681240e-02, 
       -1.6317777900311633e-02, -1.8372116865789129e-02, -2.9654166754468736e-03, -3.5401954811126055e-04, 
        2.4355022228369760e+00,  2.0016168776061879e-01,  1.0985282738580757e-01,  1.2766252509240556e-02, 
       -2.4414057689995498e-02, -1.6251161789052802e-02, -1.7089726405616346e-03, -1.9160898068638810e-03, 
        2.5933034580746872e+00,  2.1680376207396349e-01, -2.0816003865757875e-02, -3.4972802786778327e-03, 
       -2.6764066188286684e-02, -1.7854392635571970e-02,  3.5219460147529547e-04, -2.1789366800265826e-03, 
        2.2644815571737853e+00,  1.7477791651985289e-01, -1.7477791651985303e-01, -2.0820652819017429e-02, 
       -2.2307167309127181e-02, -2.2307167309126706e-02,  2.2209971661584305e-03, -2.2209971661588312e-03, 
        1.3280762633280376e+00,  6.8873390005293370e-02, -3.7870626470624885e-01, -3.8610064659975844e-02, 
       -9.2301436869363998e-03, -3.2260992258763520e-02,  5.3290259423125087e-03, -8.9388889372177130e-04, 
        7.2190296554093059e-01,  5.1538822407212193e-03,  2.5993806823215299e-02,  3.8382137320544804e-03, 
       -1.1338638972525981e-03, -3.4381947335499416e-02, -6.5463662630404866e-04,  6.6818136503074500e-04, 
        4.1727748169260243e-01, -3.7258386804597905e-03, -1.9919888758010085e-01, -3.8382137320548135e-03, 
       -1.1338638972534752e-03, -3.2313368072759835e-02,  6.5463662630354819e-04,  4.6393132852280521e-03, 
                               
        8.1840185220427497e-01,  3.9530659924141790e-03,  4.4611017446764928e-01,  1.7659483239659112e-03, 
       -5.1963206673014588e-03, -2.0444979630721963e-02, -3.0000971360621140e-03, -3.9996716592006934e-04, 
        2.0860799139066417e+00,  1.3887607915190749e-02,  2.8704596433183399e-01,  4.1266743013064410e-03, 
       -1.4695357771295256e-02, -1.9467540393028303e-02, -2.4841744929710204e-03, -2.7842352021709333e-04, 
        2.8177754775937482e+00,  2.4982469364969215e-02,  1.3467003300532823e-01,  2.2940277548075452e-03, 
       -2.0975940204581451e-02, -2.0031937337692138e-02, -1.1419214655543684e-03, -2.6674197388848914e-04, 
        3.0074496267139450e+00,  2.8078196400262975e-02, -2.8078196400262719e-02, -6.5628816075654272e-04, 
       -2.2291113936514215e-02, -2.2291113936514163e-02,  3.8260555739177558e-04, -3.8260555739177569e-04, 
        2.5933034580746881e+00,  2.0816003865758558e-02, -2.1680376207396365e-01, -3.4972802786774255e-03, 
       -1.7854392635572227e-02, -2.6764066188286983e-02,  2.1789366800265375e-03, -3.5219460147525676e-04, 
        1.4544697477545572e+00,  6.9271508916091649e-03, -4.4973913257369802e-01, -4.0330819406462224e-03, 
       -7.0401817996444040e-03, -3.3764054020119214e-02,  4.0646508571697921e-03,  2.6095781194406601e-05, 
        7.2277418406744964e-01, -4.1677105581695443e-03,  2.8023135236536021e-02, -2.3876201154679022e-03, 
       -7.5959936635850323e-04, -3.3199657075455476e-02, -4.3855489864451282e-04,  1.4414234865926383e-05, 
        4.0425568356501562e-01, -3.3091265420903547e-03, -2.0122821599342386e-01,  2.3876201154680596e-03, 
       -7.5959936635783113e-04, -2.4917931882494276e-02,  4.3855489864490373e-04, -3.6955621000352699e-04, 
                               
        7.6084873642469175e-01, -3.6904814421138479e-02,  4.1414515579073768e-01, -2.0061289360797218e-02, 
       -4.9820994058618140e-03, -1.9466440790086081e-02, -2.8764164331037327e-03,  9.6492682897369076e-04, 
        1.9431899809806050e+00, -9.4838695630359612e-02,  2.7013626166562887e-01, -1.3315992530960111e-02, 
       -1.3497480022339876e-02, -1.8183413090288651e-02, -2.0399408580720002e-03,  1.0198147641276472e-03, 
        2.6342354598592506e+00, -1.2809162129789270e-01,  1.2809162129789259e-01, -5.9085147528866377e-03, 
       -1.8762354610967265e-02, -1.8762354610967088e-02, -9.9973590292168867e-04,  9.9973590292161234e-04, 
        2.8177754775937496e+00, -1.3467003300532823e-01, -2.4982469364968788e-02,  2.2940277548078041e-03, 
       -2.0031937337692263e-02, -2.0975940204581451e-02,  2.6674197388852714e-04,  1.1419214655543684e-03, 
        2.4355022228369778e+00, -1.0985282738580707e-01, -2.0016168776061879e-01,  1.2766252509240687e-02, 
       -1.6251161789052726e-02, -2.4414057689995577e-02,  1.9160898068641100e-03,  1.7089726405617488e-03, 
        1.3850274387252517e+00, -4.6278677379411499e-02, -4.1402380087178442e-01,  2.4225516380595231e-02, 
       -6.4661984458494132e-03, -3.0373486296460778e-02,  3.7332614133449389e-03,  1.9314494067654499e-03, 
        6.9398195984660527e-01, -1.3025751711878397e-02,  1.5795917908262885e-02, -5.0010090967485561e-03, 
       -1.2013238572220258e-03, -2.9794544775782104e-02, -6.9358465235129171e-04,  1.9515282679713963e-03, 
        3.7857548166792104e-01, -1.2087608801618022e-02, -1.8900099866514941e-01,  5.0010090967483340e-03, 
       -1.2013238572224517e-03, -2.2904558275500139e-02,  6.9358465235103898e-04,  1.5319780039808881e-03, 
                               
        5.6675503083922718e-01, -7.5268518789291167e-02,  3.0700524568579995e-01, -4.1861368842359296e-02, 
       -5.0698457495897964e-03, -1.5655314361417098e-02, -2.9270768082761073e-03,  1.2354280405340475e-03, 
        1.4358116439022561e+00, -1.9781344991150729e-01,  1.9781344991150732e-01, -2.8145480425665254e-02, 
       -1.3278366801845595e-02, -1.3278366801845456e-02, -1.8121150308925303e-03,  1.8121150308926656e-03, 
        1.9431899809806050e+00, -2.7013626166562865e-01,  9.4838695630359529e-02, -1.3315992530959981e-02, 
       -1.8183413090288560e-02, -1.3497480022339674e-02, -1.0198147641276088e-03,  2.0399408580721147e-03, 
        2.0860799139066426e+00, -2.8704596433183377e-01, -1.3887607915189981e-02,  4.1266743013062745e-03, 
       -1.9467540393028414e-02, -1.4695357771295062e-02,  2.7842352021697933e-04,  2.4841744929711210e-03, 
        1.8293289250498728e+00, -2.4285963049722564e-01, -1.3644217999604569e-01,  2.2005598418681185e-02, 
       -1.8372116865788844e-02, -1.6317777900311383e-02,  3.5401954811151046e-04,  2.9654166754469412e-03, 
        1.1384058619348920e+00, -9.9223874625615535e-02, -2.7612252255954334e-01,  5.7190569078997333e-02, 
       -8.8794685107933587e-03, -2.6898835076978989e-02,  5.1265635349675900e-03,  7.4641410142385839e-05, 
        6.3102752485654290e-01, -2.6901062871494016e-02, -1.6529622982323668e-02, -1.5729096122371986e-02, 
       -3.9744222223503838e-03, -2.6679721856484789e-02, -2.2946337399474236e-03, -1.5318441703709033e-04, 
        3.1000404198245629e-01, -3.1082184954682957e-02, -1.5667545777456476e-01,  1.5729096122371902e-02, 
       -3.9744222223497307e-03, -1.7277734490433259e-02,  2.2946337399477939e-03,  1.7166702230100577e-03, 
                               
        2.2763558362451619e-01, -1.2270125405582197e-01,  1.2270125405582216e-01, -6.5804676243228263e-02, 
       -6.7577451130583438e-03, -6.7577451130582744e-03, -3.9015859601390134e-03,  3.9015859601392545e-03, 
        5.6675503083922774e-01, -3.0700524568579945e-01,  7.5268518789290390e-02, -4.1861368842360212e-02, 
       -1.5655314361417160e-02, -5.0698457495904582e-03, -1.2354280405342797e-03,  2.9270768082754663e-03, 
        7.6084873642469208e-01, -4.1414515579073780e-01,  3.6904814421139520e-02, -2.0061289360796382e-02, 
       -1.9466440790086130e-02, -4.9820994058621332e-03, -9.6492682897348031e-04,  2.8764164331032990e-03, 
        8.1840185220427752e-01, -4.4611017446764822e-01, -3.9530659924142328e-03,  1.7659483239655942e-03, 
       -2.0444979630721800e-02, -5.1963206673020955e-03,  3.9996716591999290e-04,  3.0000971360615619e-03, 
        7.4609938259350450e-01, -3.9479280796571403e-01, -3.8300006873650226e-02,  2.8156146829749602e-02, 
       -2.7860593399204932e-02, -5.5907627764129946e-03, -4.6813737713600227e-03,  3.2278283939370035e-03, 
        5.6355951897661616e-01, -2.4441834236771068e-01, -1.9416433643299941e-02,  9.7805239292669613e-02, 
       -1.7984485310210362e-02,  3.1336583368863830e-02,  1.0383347435086867e-02,  3.3547593105935511e-02, 
        3.6946581339115164e-01, -1.3727843226277259e-01, -9.2756899567130208e-02, -3.5882581089513220e-02, 
       -1.4173358881541453e-02,  3.1248837025135395e-02, -8.1829925655791651e-03,  3.3598253481107600e-02, 
        1.5533311401374345e-01, -7.1383887553887249e-02, -8.0448181189757509e-02,  3.5882581089513220e-02, 
       -1.4173358881541490e-02, -7.1521872221696117e-03,  8.1829925655791547e-03,  4.1293172180142459e-03, 
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
	  }
        }
      }
    } else if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
               (bcs.lo_type[1] == GKYL_POISSON_PERIODIC && bcs.up_type[1] == GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8 (checked visually against g2):
      const double sol[512] =  {
         2.1586267052960614e-01,  1.1494113900184900e-01, -2.2160891440115113e-02, -1.0033363538509546e-02, 
        -7.5036975921132943e-03,  3.4735460101775035e-03,  2.1388419622000235e-03,  2.0054527240191246e-03, 
         1.5189545122344580e-01,  8.8570694533605054e-02, -1.9733655611820895e-02, -8.0570412565080119e-03, 
         6.7685409640906536e-04, -3.7082045559042455e-04,  2.5842017572880105e-03, -2.1409328985628725e-04, 
         1.1533392259191108e-01,  7.1384648051902394e-02,  7.1644520061524379e-03,  3.0650158192709769e-03, 
         3.7154124105361753e-03,  6.2439289727240386e-03, -8.2988929667835050e-04,  3.6049340732027977e-03, 
         2.0164435976184503e-01,  1.0950275559554255e-01,  3.0200132477082454e-02,  1.1744791018291728e-02, 
        -5.3576293843106164e-03, -3.4127892749095893e-03, -4.4084338259451662e-03, -1.9703748065567693e-03, 
         2.2325502725940449e-01,  1.1425357838446824e-01, -1.8481830780721915e-02, -9.4398933959146669e-03, 
        -1.1342242188590531e-02, -4.0004161244795323e-03,  9.5321601239847461e-04, -2.3096413263394870e-03, 
         2.6140474630845567e-01,  8.6931415779306687e-02,  1.2058007860883976e-01, -2.6771353354680897e-02, 
        -4.8456118122639997e-03,  5.8023473555332633e-02,  2.7976152841991510e-03, -1.8139909872943488e-02, 
         3.9923702256698385e-01,  1.4670447600225761e-01, -5.9367977013402518e-02,  5.0678077012300753e-02, 
        -2.0186066294717149e-02,  4.3797695726909752e-02, -1.1654430809134692e-02, -2.6353166531615985e-02, 
         3.4662570879766441e-01,  1.6682943277557366e-01, -3.8200308246014211e-02, -1.1186232304250341e-02, 
        -2.5790206564719292e-02, -8.2498782638835838e-03,  8.4188789156725494e-03, -4.7630694364341362e-03, 
                                                          
         5.1689791587943612e-01,  5.8263941755036547e-02, -2.2859613871489739e-02,  1.0033363538509105e-02, 
        -7.9666733872370693e-03,  3.4735460101798705e-03,  2.4513201329911991e-03, -2.0054527240177637e-03, 
         4.5184750558026437e-01,  8.4634386223283203e-02, -1.9660313791848436e-02,  8.0570412565081368e-03, 
         6.9829603988514569e-04, -3.7082045559060377e-04,  2.5514022982777058e-03,  2.1409328985618688e-04, 
         4.1553825343001166e-01,  1.0182043270498713e-01,  7.2367620805570774e-03, -3.0650158192711053e-03, 
         3.6240328817590487e-03,  6.2439289727228946e-03, -8.6222734504381527e-04, -3.6049340732034773e-03, 
         5.0343323992976197e-01,  6.3702325161346079e-02,  3.1042662384792374e-02, -1.1744791018291617e-02, 
        -6.1576409161226298e-03, -3.4127892749100685e-03, -4.7852246552882657e-03,  1.9703748065564948e-03, 
         5.1062387227790884e-01,  5.8951502372419015e-02, -2.7649771863668387e-02,  9.4398933959141430e-03, 
        -5.6934179539989971e-03, -4.0004161244794638e-03,  5.0532439074343358e-03,  2.3096413263395108e-03, 
         4.9714949893947713e-01,  5.7406151518099226e-02,  4.2207809136240247e-02, -2.3228646645318689e-02, 
         1.5295286183578047e-03,  1.3302114005336851e-02, -8.8307375954184756e-04, -7.6799791017726259e-03, 
         6.3472949871672413e-01, -2.3669087048523519e-03,  1.8858640564821116e-02, -6.7807701230037298e-04, 
        -1.3698104391842032e-02, -9.2366382308450953e-04, -7.9086042580178602e-03,  5.3327755690076484e-04, 
         6.3324091899808888e-01,  6.3756479813146857e-03, -2.9176174639404251e-02,  1.1186232304250386e-02, 
        -1.9804346593439444e-02, -8.2498782638828674e-03,  4.3831636791885481e-03,  4.7630694364345456e-03, 
                                                          
         7.9669310875087918e-01,  1.0375414693349232e-01, -2.0266296955233722e-02, -9.4401510177539484e-03, 
        -7.5962157398736865e-03,  2.6659256900898546e-03,  1.7510539933307684e-03,  1.5391729148107644e-03, 
         7.4496312882187410e-01,  8.1917436794197662e-02, -1.0696170271041952e-02, -3.8001416042442735e-03, 
        -1.3765083617684950e-03,  1.8168520243091850e-03,  1.8398957356989957e-03,  1.0489600053126933e-03, 
         7.7101263805180820e-01,  9.9674385664237267e-02,  3.3254980887842074e-02,  1.8393286453245197e-02, 
         7.2827040512924985e-04,  7.6411369064037678e-03, -6.2470114804598874e-04,  4.4116124498280000e-03, 
         9.4196670141620098e-01,  1.8793144588299548e-01,  4.3756617849704256e-02,  2.0039930495653904e-02, 
        -7.3610860548943263e-03, -9.1589165610489118e-03, -4.0456909817194889e-03, -5.2879029420066122e-03, 
         8.9878642057493607e-01,  1.6656461618143703e-01, -7.1248302604095690e-02, -3.3854980021914928e-02, 
        -4.6009525242543502e-03, -1.1143075353917424e-02,  5.6392548186338003e-03, -6.4334575551852789e-03, 
         6.4345461989165043e-01,  2.8423512966575141e-02, -4.3435234660647325e-02, -2.7002793988992049e-02, 
         2.5832616692228463e-03,  1.4211291293031876e-02, -1.4914468201137644e-03,  8.2048928535644657e-03, 
         7.3875638680953348e-01,  6.7135559765420835e-02,  7.9130952232432605e-02,  3.8194881138876494e-02, 
        -1.0050761661718852e-02, -7.5917061052388143e-04, -5.8028099509534522e-03, -4.3830735634903646e-04, 
         8.8222345096854260e-01,  1.4310523382099843e-01, -1.0496546478960242e-02, -2.5300314548703986e-03, 
        -1.5365326773277715e-02, -1.1755610116889350e-02,  2.7344443531691315e-03, -6.7871046654745032e-03, 
                                                          
         1.0631958114771909e+00,  5.0444306239435549e-02, -3.8385918908420656e-02, -2.6237190677024800e-03, 
        -7.3382394066251553e-03,  4.9759462587561998e-03,  5.0976168975893060e-04, -2.0548191765771472e-04, 
         9.7626883686450805e-01,  4.4522171672997997e-02, -1.1454235674422020e-02,  6.7072891068743312e-04, 
        -6.8798661362913802e-03,  5.2448533868768932e-03, -2.4511975874243929e-04,  9.3019750414814512e-04, 
         1.0772723792237053e+00,  6.1975158927731977e-02,  8.4465886736531381e-02,  9.2088756428255498e-03, 
        -1.1022050775899156e-02,  1.6629277778719535e-02, -2.1463716579682504e-03,  7.7769310231773847e-04, 
         1.4592593018788609e+00,  9.6773808133432143e-02,  8.6635145383198373e-02,  7.3824276306553877e-03, 
        -1.8169628356648381e-02, -2.1666177752134325e-02, -1.9802841823310384e-03, -1.9331676734915578e-03, 
         1.3618658290516619e+00,  8.5682134779673375e-02, -1.4800034003590376e-01, -1.3793119643751742e-02, 
        -1.6306866669662572e-02, -2.5643767976196202e-02,  3.0557501437484154e-03, -1.9385212337234257e-03, 
         7.6245667412648599e-01,  2.9837773953860069e-02, -1.2937856967575645e-01, -1.6586431892270211e-02, 
        -5.5070760827273975e-03,  2.7563375792996905e-02,  3.1795118588102090e-03, -4.9606327326703458e-04, 
         8.8286367912837838e-01,  2.4272311064120894e-02,  1.5908799283094821e-01,  1.2706932497486551e-02, 
        -3.6932697103741656e-03, -3.2714536161833844e-03, -2.1323102614742352e-03, -1.0121599132489781e-03, 
         1.2013618239433934e+00,  5.0911208534270433e-02, -2.9699606561751119e-03,  3.0343059220695290e-03, 
        -7.8038559868938914e-03, -2.4779646387858571e-02, -2.4093783180159647e-04, -7.3232618150505940e-04, 
                                                          
         1.1414577541168380e+00, -4.4842413364417641e-03, -4.5401647444754153e-02, -2.2108034195797068e-03, 
        -6.7375309733318919e-03,  4.3032531205827850e-03, -9.7514142477539561e-05, -1.8289764674865404e-04, 
         1.0301082929209153e+00, -1.4574429930042739e-02, -1.6422458503697671e-02, -3.8590341960145040e-03, 
        -7.7601818536419875e-03,  6.2114994724910590e-03, -4.9291361855647415e-04, -3.7210412640766759e-04, 
         1.1285500148542398e+00, -3.3305613098992458e-02,  8.7284792291645957e-02, -7.1455058775012183e-03, 
        -1.1746771867359510e-02,  1.7076750786305099e-02, -1.8087451990120582e-03, -5.1934444093305811e-04, 
         1.5225609118376628e+00, -5.8805683445042833e-02,  8.9014061382215987e-02, -5.0842881575580900e-03, 
        -1.7068986026288629e-02, -2.2569774272950257e-02, -1.2640365783305148e-03,  1.4114759789596901e-03, 
         1.4123706053764706e+00, -5.5357476387942459e-02, -1.5753137114721999e-01,  7.2183584727072815e-03, 
        -1.5403909893680940e-02, -2.6364451060557272e-02,  2.2253687317128354e-03,  1.5224346609671899e-03, 
         7.8065106809459350e-01, -1.9678774248293812e-02, -1.3872732829273812e-01,  1.1388420362828854e-02, 
        -5.7747290923905788e-03,  2.6668112734835641e-02,  3.3340413959877086e-03, -2.0817094357912604e-05, 
         9.2370042204318070e-01,  1.0468664832981849e-03,  1.8019769892090598e-01,  4.8654833046057899e-04, 
        -2.3438657416457358e-03, -5.1826922176694433e-03, -1.3532315168823664e-03, -9.1294207804667053e-05, 
         1.2863851852374695e+00,  9.8583567633453483e-04,  1.5862527936420128e-03, -7.9369551534319593e-04, 
        -5.6281815036968632e-03, -2.6571005679194678e-02, -5.4296907244158628e-04, -3.0191558756314797e-04, 
                                                          
         1.0376461833777941e+00, -5.5129221073935433e-02, -5.4384076509311631e-02, -2.5299825176015598e-03, 
        -6.4879745936563564e-03,  3.9856936295335019e-03,  2.4735325927903245e-04, -4.4541089240756464e-07, 
         8.8500028578280476e-01, -6.6288092700731591e-02, -3.3167131347759046e-02, -4.7564009177981665e-03, 
        -5.5017489418357551e-03,  4.4341403371513126e-03,  3.2204438628128079e-04, -6.5405464882739204e-04, 
         8.9345440152655042e-01, -9.2173353509972092e-02,  5.1830438538806929e-02, -1.0139628277075470e-02, 
        -3.8044069399830838e-03,  1.5109880425834396e-02,  6.5791647539523612e-04, -6.1622869114581823e-04, 
         1.1686958437949939e+00, -1.2537800413393285e-01,  6.6368196919629655e-02, -5.4781225309281662e-03, 
        -1.4837875523777699e-03, -1.6425565266637097e-02,  6.8189375272534401e-04,  2.1358847447925373e-03, 
         1.0899279369791366e+00, -1.1163032235989877e-01, -1.1528576482924126e-01,  1.4113886282535629e-02, 
        -5.5135753115332156e-04, -1.9090871619691134e-02, -1.4355502897091374e-04,  2.6769683871889932e-03, 
         6.6184749305462187e-01, -4.1973757738133589e-02, -7.1464026093316230e-02,  2.3439979527187366e-02, 
        -4.0000106750650601e-04,  2.7696756840852817e-02,  2.3094072400100290e-04,  6.1470504586723948e-04, 
         8.9488445319278820e-01, -1.8082954491535259e-02,  1.6513342331838485e-01, -9.4143789520873197e-03, 
        -2.6530697323322171e-03, -3.9645368865213531e-03, -1.5317505241408177e-03,  7.9459651615785079e-04, 
         1.2139960141437611e+00, -4.3403409255577091e-02, -9.0310599971932699e-03, -5.2353526142323225e-03, 
        -6.1112712354048308e-03, -2.5759194735691038e-02, -4.6484304457016499e-04,  7.7061485432604099e-04, 
                                                          
         7.5785099050635130e-01, -1.0688886761459329e-01, -5.6977393425567437e-02,  1.9367699968464440e-03, 
        -6.8584322410197998e-03,  4.7933139496233734e-03,  9.4761939893941148e-04,  4.6672522009931088e-04, 
         5.9188466254119543e-01, -1.0026373031674923e-01, -4.2131274868565646e-02,  4.9950126553432094e-04, 
        -3.4269445401823024e-03,  2.2464678572515074e-03,  1.0335509488600284e-03, -6.0899864634152221e-04, 
         5.3798001690475372e-01, -1.0932146485925229e-01,  2.5812219731521820e-02, -5.1886423568986073e-03, 
        -9.0864446335326289e-04,  1.3712672492153478e-02,  4.2039027839748170e-04, -1.9044968547872484e-04, 
         7.3016238230855479e-01, -1.2625576691040866e-01,  5.3654241454717942e-02, -2.8170169464341132e-03, 
        -2.8034241360597916e-04, -1.0679437980498247e-02, -5.7639920843499212e-05,  1.1816433906575866e-03, 
         7.0176538868210925e-01, -1.1388579619395725e-01, -7.1687234088814028e-02,  1.0301200343465104e-02, 
        -1.6438229608980407e-03, -1.1948212390253175e-02, -7.2956594017040035e-04,  1.4468478416567725e-03, 
         5.1554237210244991e-01, -4.3855906746540885e-02,  1.4179017703572055e-02,  2.6791461107123344e-02, 
        -1.4537341183722822e-03,  2.6787579553157689e-02,  8.3931378457259087e-04, -1.1396187976591013e-03, 
         7.9085756509997984e-01, -4.6685696569033333e-02,  1.0486111165077257e-01, -2.8102425174488771e-02, 
        -6.3004124624559156e-03, -4.1290300990819368e-03, -3.6375448312048176e-03, -8.8956671670958789e-04, 
         9.6501348217330729e-01, -1.0607747254673591e-01, -2.7710688157637269e-02, -3.4208482351476936e-03, 
        -1.0550291055566489e-02, -2.2253462882684757e-02,  1.1838762814492037e-03,  1.2534203747138002e-03, 
                                                          
         2.9412461316925403e-01, -1.6090120390484253e-01, -2.9176619976449127e-02,  1.4867886025791300e-02, 
        -6.9029891588198930e-03,  2.8008528720044577e-03,  1.5315661299634606e-03, -1.6170731596124341e-03, 
         2.0573490727985275e-01, -1.1851843627656046e-01, -2.4701878441096427e-02,  1.1245346541835231e-02, 
        -2.0346162094144078e-04,  5.9582563002411070e-04,  2.3364078974739879e-03, -3.4400008788394194e-04, 
         1.6661155822244592e-01, -1.0005419388064173e-01,  9.9833575612669565e-03, -5.1283855845954841e-03, 
         2.9906913190757825e-03,  6.6914019803097670e-03, -4.9226283772221097e-04, -3.8632827345874148e-03, 
         2.6494596972064716e-01, -1.4747088028393204e-01,  3.2579048476099888e-02, -1.4042930491389011e-02, 
        -4.2569870539508094e-03, -4.3163857957252720e-03, -3.6921862219445266e-03,  2.4920665010887367e-03, 
         2.7375980358421392e-01, -1.4457823677619877e-01, -2.8012861892038388e-02,  1.6014654566958938e-02, 
        -1.0439285412608815e-02, -4.7210992088410593e-03,  1.2283460036275438e-04,  2.7257278990954230e-03, 
         2.7959914027655991e-01, -9.7090415484873821e-02,  1.1123131999185713e-01,  3.1969364884122366e-02, 
        -5.1132648219266094e-03,  5.7128210497171636e-02,  2.9521448213771064e-03,  1.8656790240568655e-02, 
         4.4007376548178490e-01, -1.7202365354967644e-01, -3.8258270923443416e-02, -6.3871557840247742e-02, 
        -1.8836662325988185e-02,  4.1886457125423347e-02, -1.0875352064543299e-02,  2.7456620652669513e-02, 
         4.3164907009174019e-01, -2.1872647698617886e-01, -3.3644094796196622e-02,  8.9456218975244150e-03, 
        -2.3614532081522272e-02, -1.0041237555219056e-02,  8.1168476750327311e-03,  5.7973112055029127e-03, 
      };

      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p  = gkyl_array_cfetch(phi_ho, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    }
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(rho_cellavg);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(perbuff);
}


// ......... CPU tests ............ //
void test_poisson_1x_p1_periodicx_ho() {
  int cells[] = {16}; // MF 2022/05/31: N=8 is broken.
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_1x(1, &cells[0], bc_tv, false);
}
void test_poisson_1x_p1_dirichletx_ho() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv, false);
}
void test_poisson_1x_p1_neumannx_dirichletx_ho() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv, false);
}
void test_poisson_1x_p1_dirichletx_neumannx_ho() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv, false);
}
void test_poisson_1x_p1_dirichletx_bias_ho() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(1, &cells[0], bc_tv, false);
}
void test_poisson_1x_p1_neumannx_dirichletx_bias_ho() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(1, &cells[0], bc_tv, false);
}

void test_poisson_1x_p2_periodicx_ho() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_1x(2, &cells[0], bc_tv, false);
}
void test_poisson_1x_p2_dirichletx_ho() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv, false);
}
void test_poisson_1x_p2_neumannx_dirichletx_ho() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv, false);
}
void test_poisson_1x_p2_dirichletx_neumannx_ho() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv, false);
}
void test_poisson_1x_p2_dirichletx_bias_ho() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(2, &cells[0], bc_tv, false);
}
void test_poisson_1x_p2_neumannx_dirichletx_bias_ho() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(2, &cells[0], bc_tv, false);
}

void test_poisson_2x_p1_periodicx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_dirichlety_ho() {
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
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_periodicx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_neumanny_dirichlety_ho() {
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
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_dirichlety_neumanny_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_neumannx_dirichletx_dirichlety_ho() {
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
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_neumannx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletvarx_dirichletvary_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  test_2x_varBC(1, cells, bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_dirichlety_bias_ho() {
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
  test_2x_bias(1, &cells[0], bc_tv, false);
}
void test_poisson_2x_p1_dirichletx_periodicy_bias_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x_bias(1, &cells[0], bc_tv, false);
}

void test_poisson_2x_p2_periodicx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_dirichlety_ho() {
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
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_periodicy_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_periodicx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_neumanny_dirichlety_ho() {
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
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_dirichlety_neumanny_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_neumannx_dirichletx_dirichlety_ho() {
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
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_neumannx_dirichlety_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletvarx_dirichletvary_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  test_2x_varBC(2, cells, bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_dirichlety_bias_ho() {
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
  test_2x_bias(2, &cells[0], bc_tv, false);
}
void test_poisson_2x_p2_dirichletx_periodicy_bias_ho() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x_bias(2, &cells[0], bc_tv, false);
}

#ifdef GKYL_HAVE_CUDA
// ......... GPU tests ............ //
void test_poisson_1x_p1_periodicx_dev() {
  int cells[] = {16}; // MF 2022/05/31: N=8 is broken.
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_1x(1, &cells[0], bc_tv, true);
}
void test_poisson_1x_p1_dirichletx_dev() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv, true);
}
void test_poisson_1x_p1_neumannx_dirichletx_dev() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv, true);
}
void test_poisson_1x_p1_dirichletx_neumannx_dev() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(1, &cells[0], bc_tv, true);
}
void test_poisson_1x_p1_dirichletx_bias_dev() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(1, &cells[0], bc_tv, true);
}
void test_poisson_1x_p1_neumannx_dirichletx_bias_dev() {
  int cells[] = {16};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(1, &cells[0], bc_tv, true);
}

void test_poisson_1x_p2_periodicx_dev() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  test_1x(2, &cells[0], bc_tv, true);
}
void test_poisson_1x_p2_dirichletx_dev() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv, true);
}
void test_poisson_1x_p2_neumannx_dirichletx_dev() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv, true);
}
void test_poisson_1x_p2_dirichletx_neumannx_dev() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x(2, &cells[0], bc_tv, true);
}
void test_poisson_1x_p2_dirichletx_bias_dev() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(2, &cells[0], bc_tv, true);
}
void test_poisson_1x_p2_neumannx_dirichletx_bias_dev() {
  int cells[] = {8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_1x_bias(2, &cells[0], bc_tv, true);
}

void test_poisson_2x_p1_periodicx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_dirichlety_dev() {
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
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_periodicx_dirichlety_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_neumanny_dirichlety_dev() {
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
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_dirichlety_neumanny_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_neumannx_dirichletx_dirichlety_dev() {
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
  test_2x(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_neumannx_dirichlety_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(1, &cells[0], bc_tv, true);
}
void gpu_test_2x_p1_dirichletvarx_dirichletvary() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  test_2x_varBC(1, cells, bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_dirichlety_bias_dev() {
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
  test_2x_bias(1, &cells[0], bc_tv, true);
}
void test_poisson_2x_p1_dirichletx_periodicy_bias_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x_bias(1, &cells[0], bc_tv, true);
}

void test_poisson_2x_p2_periodicx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  test_2x(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_dirichlety_dev() {
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
  test_2x(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_periodicy_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, true);
}
void gpu_test_2x_p2_periodicx_dirichlety() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_neumanny_dirichlety_dev() {
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
  test_2x(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_dirichlety_neumanny_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_neumannx_dirichletx_dirichlety_dev() {
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
  test_2x(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_neumannx_dirichlety_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_NEUMANN;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_2x(2, &cells[0], bc_tv, true);
}
void gpu_test_2x_p2_dirichletvarx_dirichletvary() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET_VARYING;
  test_2x_varBC(2, cells, bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_dirichlety_bias_dev() {
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
  test_2x_bias(2, &cells[0], bc_tv, true);
}
void test_poisson_2x_p2_dirichletx_periodicy_bias_dev() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  test_2x_bias(2, &cells[0], bc_tv, true);
}
#endif


TEST_LIST = {
  // 1x tests
  { "test_poisson_1x_p1_periodicx_ho", test_poisson_1x_p1_periodicx_ho },
  { "test_poisson_1x_p1_dirichletx_ho", test_poisson_1x_p1_dirichletx_ho },
  { "test_poisson_1x_p1_neumannx_dirichletx_ho", test_poisson_1x_p1_neumannx_dirichletx_ho },
  { "test_poisson_1x_p1_dirichletx_neumannx_ho", test_poisson_1x_p1_dirichletx_neumannx_ho },
  { "test_poisson_1x_p1_dirichletx_bias_ho", test_poisson_1x_p1_dirichletx_bias_ho },
  { "test_poisson_1x_p1_neumannx_dirichletx_bias_ho", test_poisson_1x_p1_neumannx_dirichletx_bias_ho },
  { "test_poisson_1x_p2_periodicx_ho", test_poisson_1x_p2_periodicx_ho },
  { "test_poisson_1x_p2_dirichletx_ho", test_poisson_1x_p2_dirichletx_ho },
  { "test_poisson_1x_p2_neumannx_dirichletx_ho", test_poisson_1x_p2_neumannx_dirichletx_ho },
  { "test_poisson_1x_p2_dirichletx_neumannx_ho", test_poisson_1x_p2_dirichletx_neumannx_ho },
  { "test_poisson_1x_p2_dirichletx_bias_ho", test_poisson_1x_p2_dirichletx_bias_ho },
  { "test_poisson_1x_p2_neumannx_dirichletx_bias_ho", test_poisson_1x_p2_neumannx_dirichletx_bias_ho },
  // 2x tests
  { "test_poisson_2x_p1_periodicx_periodicy_ho", test_poisson_2x_p1_periodicx_periodicy_ho },
  { "test_poisson_2x_p1_dirichletx_dirichlety_ho", test_poisson_2x_p1_dirichletx_dirichlety_ho },
  { "test_poisson_2x_p1_dirichletx_periodicy_ho", test_poisson_2x_p1_dirichletx_periodicy_ho },
  { "test_poisson_2x_p1_periodicx_dirichlety_ho", test_poisson_2x_p1_periodicx_dirichlety_ho },
  { "test_poisson_2x_p1_dirichletx_neumanny_dirichlety_ho", test_poisson_2x_p1_dirichletx_neumanny_dirichlety_ho },
  { "test_poisson_2x_p1_dirichletx_dirichlety_neumanny_ho", test_poisson_2x_p1_dirichletx_dirichlety_neumanny_ho },
  { "test_poisson_2x_p1_neumannx_dirichletx_dirichlety_ho", test_poisson_2x_p1_neumannx_dirichletx_dirichlety_ho },
  { "test_poisson_2x_p1_dirichletx_neumannx_dirichlety_ho", test_poisson_2x_p1_dirichletx_neumannx_dirichlety_ho },
  { "test_poisson_2x_p1_dirichletvarx_dirichletvary_ho", test_poisson_2x_p1_dirichletvarx_dirichletvary_ho },
  { "test_poisson_2x_p1_dirichletx_dirichlety_bias_ho", test_poisson_2x_p1_dirichletx_dirichlety_bias_ho },
  { "test_poisson_2x_p1_dirichletx_periodicy_bias_ho", test_poisson_2x_p1_dirichletx_periodicy_bias_ho },
  { "test_poisson_2x_p2_periodicx_periodicy_ho", test_poisson_2x_p2_periodicx_periodicy_ho },
  { "test_poisson_2x_p2_dirichletx_dirichlety_ho", test_poisson_2x_p2_dirichletx_dirichlety_ho },
  { "test_poisson_2x_p2_dirichletx_periodicy_ho", test_poisson_2x_p2_dirichletx_periodicy_ho },
  { "test_poisson_2x_p2_periodicx_dirichlety_ho", test_poisson_2x_p2_periodicx_dirichlety_ho },
  { "test_poisson_2x_p2_dirichletx_neumanny_dirichlety_ho", test_poisson_2x_p2_dirichletx_neumanny_dirichlety_ho },
  { "test_poisson_2x_p2_dirichletx_dirichlety_neumanny_ho", test_poisson_2x_p2_dirichletx_dirichlety_neumanny_ho },
  { "test_poisson_2x_p2_neumannx_dirichletx_dirichlety_ho", test_poisson_2x_p2_neumannx_dirichletx_dirichlety_ho },
  { "test_poisson_2x_p2_dirichletx_neumannx_dirichlety_ho", test_poisson_2x_p2_dirichletx_neumannx_dirichlety_ho },
  { "test_poisson_2x_p2_dirichletvarx_dirichletvary_ho", test_poisson_2x_p2_dirichletvarx_dirichletvary_ho },
  { "test_poisson_2x_p2_dirichletx_dirichlety_bias_ho", test_poisson_2x_p2_dirichletx_dirichlety_bias_ho },
  { "test_poisson_2x_p2_dirichletx_periodicy_bias_ho", test_poisson_2x_p2_dirichletx_periodicy_bias_ho },
#ifdef GKYL_HAVE_CUDA
  // 1x tests
  { "test_poisson_1x_p1_periodicx_dev", test_poisson_1x_p1_periodicx_dev },
  { "test_poisson_1x_p1_dirichletx_dev", test_poisson_1x_p1_dirichletx_dev },
  { "test_poisson_1x_p1_neumannx_dirichletx_dev", test_poisson_1x_p1_neumannx_dirichletx_dev },
  { "test_poisson_1x_p1_dirichletx_neumannx_dev", test_poisson_1x_p1_dirichletx_neumannx_dev },
  { "test_poisson_1x_p1_dirichletx_bias_dev", test_poisson_1x_p1_dirichletx_bias_dev },
  { "test_poisson_1x_p1_neumannx_dirichletx_bias_dev", test_poisson_1x_p1_neumannx_dirichletx_bias_dev },
  { "test_poisson_1x_p2_periodicx_dev", test_poisson_1x_p2_periodicx_dev },
  { "test_poisson_1x_p2_dirichletx_dev", test_poisson_1x_p2_dirichletx_dev },
  { "test_poisson_1x_p2_neumannx_dirichletx_dev", test_poisson_1x_p2_neumannx_dirichletx_dev },
  { "test_poisson_1x_p2_dirichletx_neumannx_dev", test_poisson_1x_p2_dirichletx_neumannx_dev },
  { "test_poisson_1x_p2_dirichletx_bias_dev", test_poisson_1x_p2_dirichletx_bias_dev },
  { "test_poisson_1x_p2_neumannx_dirichletx_bias_dev", test_poisson_1x_p2_neumannx_dirichletx_bias_dev },
  // 2x tests
  { "test_poisson_2x_p1_periodicx_periodicy_dev", test_poisson_2x_p1_periodicx_periodicy_dev },
  { "test_poisson_2x_p1_dirichletx_dirichlety_dev", test_poisson_2x_p1_dirichletx_dirichlety_dev },
  { "test_poisson_2x_p1_dirichletx_periodicy_dev", test_poisson_2x_p1_dirichletx_periodicy_dev },
  { "test_poisson_2x_p1_periodicx_dirichlety_dev", test_poisson_2x_p1_periodicx_dirichlety_dev },
  { "test_poisson_2x_p1_dirichletx_neumanny_dirichlety_dev", test_poisson_2x_p1_dirichletx_neumanny_dirichlety_dev },
  { "test_poisson_2x_p1_dirichletx_dirichlety_neumanny_dev", test_poisson_2x_p1_dirichletx_dirichlety_neumanny_dev },
  { "test_poisson_2x_p1_neumannx_dirichletx_dirichlety_dev", test_poisson_2x_p1_neumannx_dirichletx_dirichlety_dev },
  { "test_poisson_2x_p1_dirichletx_neumannx_dirichlety_dev", test_poisson_2x_p1_dirichletx_neumannx_dirichlety_dev },
  { "test_poisson_2x_p1_dirichletx_dirichlety_bias_dev", test_poisson_2x_p1_dirichletx_dirichlety_bias_dev },
  { "test_poisson_2x_p1_dirichletx_periodicy_bias_dev", test_poisson_2x_p1_dirichletx_periodicy_bias_dev },
  { "test_poisson_2x_p2_periodicx_periodicy_dev", test_poisson_2x_p2_periodicx_periodicy_dev },
  { "test_poisson_2x_p2_dirichletx_dirichlety_dev", test_poisson_2x_p2_dirichletx_dirichlety_dev },
  { "test_poisson_2x_p2_dirichletx_periodicy_dev", test_poisson_2x_p2_dirichletx_periodicy_dev },
  { "test_poisson_2x_p2_dirichletx_neumanny_dirichlety_dev", test_poisson_2x_p2_dirichletx_neumanny_dirichlety_dev },
  { "test_poisson_2x_p2_dirichletx_dirichlety_neumanny_dev", test_poisson_2x_p2_dirichletx_dirichlety_neumanny_dev },
  { "test_poisson_2x_p2_neumannx_dirichletx_dirichlety_dev", test_poisson_2x_p2_neumannx_dirichletx_dirichlety_dev },
  { "test_poisson_2x_p2_dirichletx_neumannx_dirichlety_dev", test_poisson_2x_p2_dirichletx_neumannx_dirichlety_dev },
  { "test_poisson_2x_p2_dirichletx_dirichlety_bias_dev", test_poisson_2x_p2_dirichletx_dirichlety_bias_dev },
  { "test_poisson_2x_p2_dirichletx_periodicy_bias_dev", test_poisson_2x_p2_dirichletx_periodicy_bias_dev },
#endif
  { NULL, NULL },
};