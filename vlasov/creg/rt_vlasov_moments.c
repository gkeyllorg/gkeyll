#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>

void
evalDistFunc1x1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  double n = 1.0*sin(2*M_PI*x);
  double ux = 0.1*cos(2*M_PI*x);
  double Txx = 0.75 + 0.25*cos(2*M_PI*x);
   
  double u2 = (vx-ux)*(vx-ux)/(2*Txx);
  fout[0] = n/sqrt(2*M_PI*Txx)*exp(-u2);
}

void
evalDistFunc1x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1], vy = xn[2];
  double n = 1.0*sin(2*M_PI*x);
  double ux = 0.1*cos(2*M_PI*x);
  double uy = 0.2*sin(2*M_PI*x);
  double Txx = 0.75 + 0.25*cos(2*M_PI*x);
  double Tyy = 0.75 + 0.25*sin(2*M_PI*x);
  double Txy = 0.1 + 0.01*sin(2*M_PI*x)*cos(2*M_PI*x);

  double detT = Txx*Tyy-Txy*Txy;
  double cx = vx-ux;
  double cy = vy-uy;

  double u2 = (cx*(cx*Tyy-cy*Txy)+cy*(cy*Txx-cx*Txy))/(2*detT);
  fout[0] = n/(2*M_PI*sqrt(detT))*exp(-u2);
}

void
evalDistFunc1x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  double n = 1.0*sin(2*M_PI*x);
  double ux = 0.1*cos(2*M_PI*x);
  double uy = 0.2*sin(2*M_PI*x);
  double uz = 0.1*cos(2*M_PI*x);
  double Txx = 0.75 + 0.25*cos(2*M_PI*x);
  double Tyy = 0.75 + 0.25*sin(2*M_PI*x);
  double Tzz = 0.75 + 0.1*sin(2*M_PI*x);
  double Txy = 0.5 + 0.1*sin(2*M_PI*x);
  double Txz = 0.25 + 0.1*sin(2*M_PI*x);
  double Tyz = 0.125 + 0.1*sin(2*M_PI*x);   

  double cx = vx-ux;
  double cy = vy-uy;
  double cz = vz-uz;

  double detT = Txx*(Tyy*Tzz-Tyz*Tyz)-Txy*(Txy*Tzz-Txz*Tyz)+Txz*(Txy*Tyz-Txz*Tyy);
  double u2 = cx*(cx*(Tyy*Tzz-Tyz*Tyz)+cy*(Txz*Tyz-Txy*Tzz)+cz*(Txy*Tyz-Txz*Tyy))+cy*(cx*(Txz*Tyz-Txy*Tzz)+cy*(Txx*Tzz-Txz*Txz)+cz*(Txy*Txz-Txx*Tyz))+cz*(cx*(Txy*Tyz-Txz*Tyy)+cy*(Txy*Txz-Txx*Tyz)+cz*(Txx*Tyy-Txy*Txy));
  u2 = u2/(2*detT);
  fout[0] = n/sqrt((2*M_PI)*(2*M_PI)*(2*M_PI)*detT)*exp(-u2);
}

void
evalDistFunc2x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3];
  double n = 1.0*sin(2*M_PI*x)*sin(2*M_PI*y);
  double ux = 0.1*cos(2*M_PI*x)*cos(2*M_PI*y);
  double uy = 0.2*sin(2*M_PI*x)*sin(2*M_PI*y);
  double Txx = 0.75 + 0.25*cos(2*M_PI*x)*cos(2*M_PI*y);
  double Tyy = 0.75 + 0.25*sin(2*M_PI*x)*sin(2*M_PI*y);
  double Txy = 0.1 + 0.01*sin(2*M_PI*x)*cos(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*y);

  double detT = Txx*Tyy-Txy*Txy;
  double cx = vx-ux;
  double cy = vy-uy;

  double u2 = (cx*(cx*Tyy-cy*Txy)+cy*(cy*Txx-cx*Txy))/(2*detT);
  fout[0] = n/(2*M_PI*sqrt(detT))*exp(-u2);
}

void
evalDistFunc2x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3], vz = xn[4];
  double n = 1.0*sin(2*M_PI*x)*sin(2*M_PI*y);
  double ux = 0.1*cos(2*M_PI*x)*cos(2*M_PI*y);
  double uy = 0.2*sin(2*M_PI*x)*sin(2*M_PI*y);
  double uz = 0.1*cos(2*M_PI*x)*cos(2*M_PI*y);

  double Txx = 0.75 + 0.25*cos(2*M_PI*x)*cos(2*M_PI*y);
  double Tyy = 0.75 + 0.25*sin(2*M_PI*x)*sin(2*M_PI*y);
  double Tzz = 0.75 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y);
  double Txy = 0.5 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y);
  double Txz = 0.25 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y);
  double Tyz = 0.125 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y); 

  double cx = vx-ux;
  double cy = vy-uy;
  double cz = vz-uz;

  double detT = Txx*(Tyy*Tzz-Tyz*Tyz)-Txy*(Txy*Tzz-Txz*Tyz)+Txz*(Txy*Tyz-Txz*Tyy);
  double u2 = cx*(cx*(Tyy*Tzz-Tyz*Tyz)+cy*(Txz*Tyz-Txy*Tzz)+cz*(Txy*Tyz-Txz*Tyy))+cy*(cx*(Txz*Tyz-Txy*Tzz)+cy*(Txx*Tzz-Txz*Txz)+cz*(Txy*Txz-Txx*Tyz))+cz*(cx*(Txy*Tyz-Txz*Tyy)+cy*(Txy*Txz-Txx*Tyz)+cz*(Txx*Tyy-Txy*Txy));
  u2 = u2/(2*detT);
  fout[0] = n/sqrt((2*M_PI)*(2*M_PI)*(2*M_PI)*detT)*exp(-u2);
}

void
evalDistFunc3x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3], vy = xn[4], vz = xn[5];
  double n = 1.0*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
  double ux = 0.1*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
  double uy = 0.2*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
  double uz = 0.1*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);

  double Txx = 0.75 + 0.25*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
  double Tyy = 0.75 + 0.25*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
  double Tzz = 0.75 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
  double Txy = 0.5 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
  double Txz = 0.25 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
  double Tyz = 0.125 + 0.1*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);

  double cx = vx-ux;
  double cy = vy-uy;
  double cz = vz-uz;

  double detT = Txx*(Tyy*Tzz-Tyz*Tyz)-Txy*(Txy*Tzz-Txz*Tyz)+Txz*(Txy*Tyz-Txz*Tyy);
  double u2 = cx*(cx*(Tyy*Tzz-Tyz*Tyz)+cy*(Txz*Tyz-Txy*Tzz)+cz*(Txy*Tyz-Txz*Tyy))+cy*(cx*(Txz*Tyz-Txy*Tzz)+cy*(Txx*Tzz-Txz*Txz)+cz*(Txy*Txz-Txx*Tyz))+cz*(cx*(Txy*Tyz-Txz*Tyy)+cy*(Txy*Txz-Txx*Tyz)+cz*(Txx*Tyy-Txy*Txy));
  u2 = u2/(2*detT);
  fout[0] = n/sqrt((2*M_PI)*(2*M_PI)*(2*M_PI)*detT)*exp(-u2);
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct moment_inp {
  int cdim, vdim, poly_order;
  int ccells[3], vcells[3];
  int nloop;
  bool use_gpu;
  evalf_t eval; // function to project
};

struct moment_inp
get_inp(int argc, char **argv)
{
  int c, cdim = 2, vdim = 2, poly_order = 2, nloop = 10;
  bool use_gpu = false;
  evalf_t eval = evalDistFunc2x2v;
  while ((c = getopt(argc, argv, "+hgc:v:p:n:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: app_vlasov_moments -c CDIM -v VDIM -p POLYORDER -n NLOOP -g\n");
        exit(-1);
        break;
      
      case 'c':
        cdim = atoi(optarg);
        break;

      case 'v':
        vdim = atoi(optarg);
        break;

      case 'p':
        poly_order = atoi(optarg);
        break;

      case 'n':
        nloop = atoi(optarg);
        break;          

      case 'g':
        use_gpu = true;
        break;

      case '?':
        break;
    }
  }
  
  if (cdim == 1 && vdim == 1)
    eval = evalDistFunc1x1v;
  else if (cdim == 1 && vdim == 2)
    eval = evalDistFunc1x2v;
  else if (cdim == 1 && vdim == 3)
    eval = evalDistFunc1x3v;
  else if (cdim == 2 && vdim == 3)
    eval = evalDistFunc2x3v;
  else if (cdim == 3 && vdim == 3)
    eval = evalDistFunc3x3v;
  
  return (struct moment_inp) {
    .cdim = cdim,
    .vdim = vdim,
    .poly_order = poly_order,
    .ccells = { 8, 8, 8 },
    .vcells = { 16, 16, 16 },
    .nloop = nloop,
    .eval = eval,
    .use_gpu = use_gpu,
   };
}

int
main(int argc, char **argv)
{
  struct moment_inp inp = get_inp(argc, argv);

#ifdef GKYL_HAVE_CUDA
  if (inp.use_gpu)
    printf("Running moment calculation on GPU with:\n");
  else
    printf("Running moment calculation on CPU with:\n");
#else
  printf("Running moment calculation on CPU with:\n");
#endif

  printf("cdim = %d; vdim = %d; poly_order = %d\n", inp.cdim, inp.vdim, inp.poly_order);
  printf("cells = [");
  for (int d=0; d<inp.cdim; ++d)
    printf("%d ", inp.ccells[d]);
  for (int d=0; d<inp.vdim; ++d)
    printf("%d ", inp.vcells[d]);
  printf("]\n");
    
  printf("nloop = %d\n", inp.nloop);
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    .lower = { -6.0, -6.0, -6.0 },
    .upper = { 6.0, 6.0, 6.0 }, 
    .cells = { inp.vcells[0], inp.vcells[1], inp.vcells[2] },
    .num_init = 1,     
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = inp.eval,
      .ctx_func = 0,
    },
    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2IJ },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "vlasov-moment",

    .cdim = inp.cdim, .vdim = inp.vdim,
    .lower = { -1.0, -1.0, -1.0 },
    .upper = { 1.0, 1.0, 1.0 },
    .cells = { inp.ccells[0], inp.ccells[1], inp.ccells[2] },
    .poly_order = inp.poly_order,

    .num_species = 1,
    .species = { elc },
    .field = field,

    .parallelism = {
      .use_gpu = inp.use_gpu,
    },
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);
  gkyl_vlasov_app_apply_ic(app, 0.0);

  struct timespec tm_start = gkyl_wall_clock();
  // time with volume term
  for (int i=0; i<inp.nloop; ++i)
    gkyl_vlasov_app_calc_mom(app);
  
  double tm_tot = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));
    
  printf("Time for full update: %g [s]\n", tm_tot/inp.nloop);

  gkyl_vlasov_app_write_mom(app, 0.0, 0);

  gkyl_vlasov_app_release(app);
  
  return 0;
}
