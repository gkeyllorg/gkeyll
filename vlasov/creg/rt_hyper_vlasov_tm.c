#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

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
  int nx, ny, nz = 8;
  int nvx, nvy, nvz = 16;
  bool use_gpu = false;
  while ((c = getopt(argc, argv, "+hgc:d:p:n:x:y:z:u:v:w:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: app_vlasov_kerntm -c CDIM -d VDIM -p POLYORDER -x NX -y NY -z NZ -u VX -v VY -w VZ -n NLOOP -g\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;        
      
      case 'c':
        cdim = atoi(optarg);
        break;

      case 'd':
        vdim = atoi(optarg);
        break;

      case 'p':
        poly_order = atoi(optarg);
        break;

      case 'n':
        nloop = atoi(optarg);
        break;          

      case 'x':
        nx = atoi(optarg);
        break;          

      case 'y':
        ny = atoi(optarg);
        break;          

      case 'z':
        nz = atoi(optarg);
        break;          

      case 'u':
        nvx = atoi(optarg);
        break;          

      case 'v':
        nvy = atoi(optarg);
        break;          

      case 'w':
        nvz = atoi(optarg);
        break;          

      case '?':
        break;
    }
  }
  
  return (struct kerntm_inp) {
    .cdim = cdim,
    .vdim = vdim,
    .poly_order = poly_order,
    .ccells = { nx, ny, nz },
    .vcells = { nvx, nvy, nvz},
    .nloop = nloop,
    .use_gpu = use_gpu,
  };
}

int
main(int argc, char **argv)
{
  struct kerntm_inp inp = get_inp(argc, argv);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  if (inp.use_gpu) {
    printf("Running kernel timers on GPU with:\n");
    use_gpu = true;
  } else {
    printf("Running kernel timers on CPU with:\n");
  }
#else
  printf("Running kernel timers on CPU with:\n");
#endif

  int cdim = inp.cdim;
  int vdim = inp.vdim;
  int poly_order = inp.poly_order;

  int cells[6];
  int ghost[6];
  double lower[6];
  double upper[6];
  int velcells[6];
  int velghost[6];
  double vellower[6];
  double velupper[6];
  int up_dirs[GKYL_MAX_DIM];
  int zero_flux_flags[GKYL_MAX_DIM];
  
  printf("cdim = %d; vdim = %d; poly_order = %d\n", inp.cdim, inp.vdim, inp.poly_order);
  printf("cells = [");
  for (int d=0; d<inp.cdim; ++d) {
    printf("%d ", inp.ccells[d]);
    cells[d] = inp.ccells[d];
    lower[d] = 0.;
    upper[d] = 1.;
    ghost[d] = 1;
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
  }
  for (int d=0; d<inp.vdim; ++d) {
    printf("%d ", inp.vcells[d]);
    cells[d+cdim] = inp.vcells[d];
    lower[d+cdim] = 0.;
    upper[d+cdim] = 1.;
    ghost[d+cdim] = 0;
    up_dirs[d+cdim] = d+cdim;
    zero_flux_flags[d+cdim] = 1;

    velcells[d] = inp.vcells[d];
    vellower[d] = 0.;
    velupper[d] = 1.;
    velghost[d] = 0;
  }
  printf("]\n");
    
  printf("nloop = %d\n", inp.nloop);
  
  // initialize grid and ranges
  int pdim = cdim+vdim;

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid velGrid;
  struct gkyl_range velRange, velRange_ext;
  gkyl_rect_grid_init(&velGrid, vdim, vellower, velupper, velcells);
  gkyl_create_grid_ranges(&velGrid, velghost, &velRange_ext, &velRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  struct gkyl_basis basis, confBasis, velBasis; // phase-space, conf-space basis
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;

  int nem = confRange_ext.volume*confBasis.num_basis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.num_basis, confRange_ext.volume);
    qmem_d = qmem_h->data;
  } else {
    qmem_d = qmem->data;
  }
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(qmem, qmem_h);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr1(use_gpu, velBasis.num_basis, velRange.volume);
  struct gkyl_array *gamma_inv = mkarr1(use_gpu, velBasis.num_basis, velRange.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velRange, &velBasis, inp_vmap, false, use_gpu);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confRange, &confRange_ext, &confBasis, inp_pmap, use_gpu);

  gkyl_dg_vlasov_calc_hamil(&velGrid, &velBasis, &velRange, 
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, use_gpu);

  // Select the number of nodes, with case for hybrid-tensor.
  bool use_lo = false;
  int highorder = use_lo ? 0 : 1;
  int num_surf_vel_nodes = vdim*pow(poly_order+1+highorder,pdim - 1);
  if ((basis.b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
    num_surf_vel_nodes = (int) vdim*(pow(poly_order+1+highorder,vdim - 1) + pow(poly_order,cdim));
  }
  
  // Sturcture pointers for input objects (but not used)
  int num_pt_indices[3] = { 1 , 6, 18 }; 
  struct gkyl_array *poisson_tensor_conf = mkarr1(use_gpu, confBasis.num_basis*num_pt_indices[vdim-1], confRange.volume );
  struct gkyl_array *pot_tot = mkarr1(use_gpu, confBasis.num_basis*4, confRange_ext.volume );
  struct gkyl_array *vel_flux_surf = mkarr1(use_gpu, num_surf_vel_nodes, phaseRange_ext.volume );
  struct gkyl_array *f_no_J = mkarr1(use_gpu, fin->ncomp, fin->size); ;
  struct gkyl_array *rad = mkarr1(use_gpu, vdim*velBasis.num_basis, velRange.volume);
  
  struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
    .phase_grid = &phaseGrid, 
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .hamil_range = &velRange,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .use_lo = false,
    .use_gpu = use_gpu,
  }; 
  struct gkyl_dg_vlasov_vel_flux_surf *calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux); 

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confRange,
    .hamil_range = &velRange,
    .phase_range = &phaseRange,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .poisson_tensor_conf = poisson_tensor_conf,
    .hamil = hamil,
    .qmem = qmem, 
    .pot_tot = pot_tot, 
    .vel_flux_surf = vel_flux_surf, 
    .f_no_J = f_no_J, 
    .rad = rad, 
    .use_lo = false,
    .use_gpu = use_gpu,
  };  
  // Construct Vlasov equation and Hyper DG object for updating equation. 
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  // run hyper_dg_advance
  int nrep = inp.nloop;
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  struct timespec tm_start = gkyl_wall_clock();
  for(int n=0; n<nrep; n++) { 
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_vlasov_vel_flux_surf_advance(calc_vel_flux, &confRange, &phaseRange,
     poisson_tensor_conf, hamil, qmem, pot_tot, rad, 
     f_no_J, cflrate, vel_flux_surf);  
    gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs); 
  }

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  
  double tm_tot = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));
  printf("Avg time for vlasov hyper dg: %g [s]\n", tm_tot/inp.nloop);
  
  return 0;
}
