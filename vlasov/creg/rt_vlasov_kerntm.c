#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct kerntm_inp {
  int cdim, vdim, poly_order;
  int ccells[3], vcells[3];
  int nloop;
  bool use_gpu;
};

struct kerntm_inp
get_inp(int argc, char **argv)
{
  int c, cdim = 2, vdim = 2, poly_order = 2, nloop = 10;
  bool use_gpu = false;
  while ((c = getopt(argc, argv, "+hgc:v:p:n:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: app_vlasov_kerntm -c CDIM -v VDIM -p POLYORDER -n NLOOP -g\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
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

      case '?':
        break;
    }
  }
  
  return (struct kerntm_inp) {
    .cdim = cdim,
    .vdim = vdim,
    .poly_order = poly_order,
    .ccells = { 8, 8, 8 },
    .vcells = { 16, 16, 16 },
    .nloop = nloop,
    .use_gpu = use_gpu,
  };
}

int
main(int argc, char **argv)
{
  struct kerntm_inp inp = get_inp(argc, argv);

#ifdef GKYL_HAVE_CUDA
  if (inp.use_gpu)
    printf("Running kernel timers on GPU with:\n");
  else
    printf("Running kernel timers on CPU with:\n");
#else
  printf("Running kernel timers on CPU with:\n");
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
      .func = evalDistFunc,
      .ctx_func = 0,
    },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "kern-timer",

    .cdim = inp.cdim, .vdim = inp.vdim,
    .lower = { -1.0, -1.0, -1.0 },
    .upper = { 1.0, 1.0, 1.0 },
    .cells = { inp.ccells[0], inp.ccells[1], inp.ccells[2] },
    .poly_order = inp.poly_order,

    .parallelism = {
      .use_gpu = inp.use_gpu,
    },

    .num_species = 1,
    .species = { elc },
    .field = field
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);
  gkyl_vlasov_app_apply_ic(app, 0.0);

  struct timespec tm_start = gkyl_wall_clock();
  // time with volume term
  for (int i=0; i<inp.nloop; ++i)
    gkyl_vlasov_app_species_ktm_rhs(app, 1);

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  
  double tm_tot = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));

  tm_start = gkyl_wall_clock();
  // time without volume term
  for (int i=0; i<inp.nloop; ++i)
    gkyl_vlasov_app_species_ktm_rhs(app, 0);

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif  
  
  double tm_surf_only = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));

  int cdim = inp.cdim, vdim = inp.vdim, pdim = cdim + vdim;

  // compute how many volume updates we did per loop
  long nvol_updates = 1;
  for (int i=0; i<cdim; ++i)
    nvol_updates *= inp.ccells[i];
  for (int i=0; i<vdim; ++i)
    nvol_updates *= inp.vcells[i];

  long nsurf_updates[pdim];
  
  for (int d=0; d<pdim; ++d) {
    nsurf_updates[d] = 1;
    
    for (int i=0; i<cdim; ++i)
      if (i == d)
        nsurf_updates[d] *= inp.ccells[i]+1;
      else
        nsurf_updates[d] *= inp.ccells[i];

    for (int i=cdim; i<pdim; ++i)
      if (i == d)
        nsurf_updates[d] *= inp.vcells[i-cdim]-1;
      else
        nsurf_updates[d] *= inp.vcells[i-cdim];
  }

  printf("Volume kernel calls: %ld\n", nvol_updates);
  long nsurf_updates_tot = 0;
  printf("Surface kernels calls: [");
  for (int d=0; d<pdim; ++d) {
    printf("%ld ", nsurf_updates[d]);
    nsurf_updates_tot += nsurf_updates[d];
  }
  printf("]\n");
  printf("Net surface kernel calls: %ld\n", nsurf_updates_tot);

  double tm_vol = tm_tot - tm_surf_only;
  printf("Volume updates took:  %g [s]\n", tm_vol/inp.nloop);
  printf("Surface updates took: %g [s]\n", tm_surf_only/inp.nloop);
    
  printf("Time for full update: %g [s]\n", tm_tot/inp.nloop);

  gkyl_vlasov_app_release(app);
  
  return 0;
}
