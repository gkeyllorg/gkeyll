#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_dist_type_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_dist_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Density
void 
eval_density(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;  
}

// Drift velocity
void
eval_drift_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 0.5;
  // double x = xn[0];
  // fout[0] = 3 * (1/cosh(x/0.2));
}

// T/m
void
eval_temperature_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

struct lte_test_ctx{double n; double drift; double T_over_m;};

void 
eval_lte_1x1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct lte_test_ctx *lte =ctx;
  double v = xn[1];
  double dv = v-lte->drift;
  fout[0] = lte-> n/sqrt(2.0*GKYL_PI*lte->T_over_m)*exp(-dv*dv/(2.0*lte->T_over_m));
}

void
test_1x1v_lte(int poly_order)
{
  double lower[] = {0.1, -8.0}, upper[]={1.0, 8.0};
  int cells[] = {2, 32};
  int cdim = 1, vdim =1;
  int pdim = cdim + vdim; 

  double conf_lower[] = {lower[0]}, conf_upper[] = {upper[0]};
  double vel_lower[] = {lower[1]}, vel_upper[] = {upper[1]};
  int conf_cells[] = {cells[0]};
  int vel_cells[] = {cells[1]};

  // Grids
  struct gkyl_rect_grid phase_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, conf_lower, conf_upper, conf_cells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, vel_lower, vel_upper, vel_cells);

  // Basis
  struct gkyl_basis phase_basis, conf_basis, vel_basis;
  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&vel_basis, vdim, poly_order); 

  // Ranges
  int conf_ghost[] = {1};
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[] = {0};
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int phase_ghost[] = {1,0};
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  // User provided LTE moments
  struct gkyl_array *density = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *drift = mkarr(vdim * conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *temperature = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *moms = mkarr((vdim+2) * conf_basis.num_basis, conf_local_ext.volume);  
  
  // Project the user provided moments onto conf space
  gkyl_proj_on_basis *proj_density = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, 1, eval_density, NULL);
  gkyl_proj_on_basis *proj_drift = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, vdim, eval_drift_1v, NULL);
  gkyl_proj_on_basis *proj_temperature = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, 1, eval_temperature_1v, NULL);

  gkyl_proj_on_basis_advance(proj_density, 0.0, &conf_local, density);
  gkyl_proj_on_basis_advance(proj_drift, 0.0, &conf_local, drift);
  gkyl_proj_on_basis_advance(proj_temperature, 0.0, &conf_local, temperature);

  gkyl_array_set_offset_range(moms, 1.0, density, 0*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms, 1.0, drift, 1*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms, 1.0, temperature, (vdim+1)*conf_basis.num_basis, &conf_local);

  struct gkyl_dist_proj_type *dist_proj = gkyl_vlasov_lte_new(cdim, pdim, poly_order, 
    conf_basis.num_basis, phase_basis.num_basis, GKYL_MODEL_DEFAULT);

  TEST_CHECK(dist_proj != NULL);
  TEST_CHECK(dist_proj->dist_id == GKYL_DIST_TYPE_LTE);
  TEST_CHECK(gkyl_dist_proj_type_num_mom(dist_proj) == vdim+2);

  // Kernel test
  double xc_test[] = {0.5,0.5};
  double dx_test[] = {phase_grid.dx[0], phase_grid.dx[1]};
  int idx_test[] = {1,1};
  double moms_test[] = {
    1.0,  // n
    0.5,  // Drift
    1.0   // Temperature
  };
  double f_test[1];
  gkyl_dist_proj_type_calc(dist_proj, xc_test, dx_test, idx_test, moms_test, f_test);

  double f_expected = 1.0/sqrt(2.0*GKYL_PI);
  TEST_CHECK(gkyl_compare_double(f_test[0], f_expected,1e-14));

 
  xc_test[1] = 1.5;
  gkyl_dist_proj_type_calc(dist_proj, xc_test, dx_test, idx_test, moms_test, f_test);
  
  f_expected = 1.0/sqrt(2.0*GKYL_PI) * exp(-0.5);
  TEST_CHECK(gkyl_compare_double(f_test[0], f_expected, 1e-14));


  // projection updater to compute the distribution
  struct gkyl_vlasov_dist_proj_on_basis_inp inp = {
    .dist_proj = dist_proj,
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &conf_basis,
    .vel_basis = &vel_basis,
    .phase_basis = &phase_basis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_local_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .model_id = GKYL_MODEL_DEFAULT,
    .use_gpu = false,
    .quad_type = GKYL_GAUSS_QUAD,
  };

  gkyl_vlasov_dist_proj_on_basis *proj_dist = gkyl_vlasov_dist_proj_on_basis_inew(&inp);
  TEST_CHECK(proj_dist != NULL);
  
  struct gkyl_array *distf = mkarr(phase_basis.num_basis, phase_local_ext.volume);

  gkyl_vlasov_dist_proj_on_basis_moments_advance(proj_dist, &conf_local, moms);
  gkyl_vlasov_dist_proj_on_basis_advance(proj_dist, &phase_local, &conf_local, distf);

  struct lte_test_ctx lte_ctx = {.n = 1.0, .drift = 0.5, .T_over_m = 1.0,};
  struct gkyl_array *distf_ref = mkarr(phase_basis.num_basis, phase_local_ext.volume);

  gkyl_proj_on_basis *proj_ref = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    poly_order+1, 1, eval_lte_1x1v, &lte_ctx);

  gkyl_proj_on_basis_advance(proj_ref, 0.0, &phase_local, distf_ref);



  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &phase_local);

  while (gkyl_range_iter_next(&iter)) {

    long idx = gkyl_range_idx(&phase_local_ext, iter.idx);
    const double *f = gkyl_array_cfetch(distf, idx);
    const double *f_ref = gkyl_array_cfetch(distf_ref, idx);

    for (int k=0; k<phase_basis.num_basis; ++k) {
      TEST_CHECK(gkyl_compare_double(f[k], f_ref[k], 1e-12));

    }
  }

  char fname[1024];
  sprintf(fname, "ctest_dist_proj_on_basis_lte_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&phase_grid, &phase_local, 0, distf, fname);

  gkyl_array_release(density);
  gkyl_array_release(drift);
  gkyl_array_release(temperature);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ref);
  gkyl_proj_on_basis_release(proj_density);
  gkyl_proj_on_basis_release(proj_drift);
  gkyl_proj_on_basis_release(proj_temperature);
  gkyl_proj_on_basis_release(proj_ref);
  gkyl_vlasov_dist_proj_on_basis_release(proj_dist);
  gkyl_dist_proj_type_release(dist_proj);  
}

// Velocity drifts
void
eval_bimax_drift_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 0.5;
  fout[1] = -0.25;
}

// T_parallel/m
void 
eval_bimax_Tpar(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

// T_perpendicular/m
void
eval_bimax_Tperp(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 2.0;
}

struct bimax_test_ctx {double n; double drift_par; double drift_perp; double Tpar_over_m; double Tperp_over_m;};


void
eval_bimax_1x2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct bimax_test_ctx *bimax = ctx;

  // xn[0] = x
  // xn[1] = v_parallel
  // xn[2] = v_perpendicular

  double vpar = xn[1]-bimax->drift_par;
  double vperp = xn[2]-bimax->drift_perp;

  double amplitude = bimax->n / (sqrt(2.0*GKYL_PI*bimax->Tpar_over_m) *
   sqrt(2.0*GKYL_PI*bimax->Tperp_over_m));
  fout[0] = amplitude *exp(-vpar*vpar/ (2.0*bimax->Tpar_over_m) -vperp*vperp / 
   (2.0*bimax->Tperp_over_m));
}

void
test_1x2v_bimax(int poly_order)
{
  double lower[] = {0.1, -8.0, -8.0}, upper[] = {1.0, 8.0, 8.0};
  int cells[] = {2, 16, 16};
  int cdim = 1, vdim = 2;
  int pdim = cdim+vdim;


  double conf_lower[] = {lower[0]}, conf_upper[] = {upper[0]};
  double vel_lower[] = {lower[1], lower[2]};
  double vel_upper[] = {upper[1], upper[2]};
  int conf_cells[] = {cells[0]};
  int vel_cells[] = {cells[1], cells[2]};

  // Grids
  struct gkyl_rect_grid phase_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, conf_lower, conf_upper, conf_cells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, vel_lower, vel_upper, vel_cells);

  // Basis
  struct gkyl_basis phase_basis, conf_basis, vel_basis;
  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&vel_basis, vdim, poly_order);

  // Ranges
  int conf_ghost[] = {1};
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[] = {0, 0};
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int phase_ghost[] = {1, 0, 0};
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  // User provided Bi-Max moments
  struct gkyl_array *density = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *drift = mkarr(vdim*conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *Tpar = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *Tperp = mkarr(conf_basis.num_basis, conf_local_ext.volume);
  struct gkyl_array *moms = mkarr((vdim+3)*conf_basis.num_basis, conf_local_ext.volume);

  // Project the user provided moments onto conf space
  gkyl_proj_on_basis *proj_density = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, 1, eval_density, NULL);
  gkyl_proj_on_basis *proj_drift = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, vdim, eval_bimax_drift_2v, NULL);
  gkyl_proj_on_basis *proj_Tpar = gkyl_proj_on_basis_new(&conf_grid,&conf_basis,
    poly_order+1, 1, eval_bimax_Tpar, NULL);
  gkyl_proj_on_basis *proj_Tperp = gkyl_proj_on_basis_new(&conf_grid, &conf_basis,
    poly_order+1, 1, eval_bimax_Tperp, NULL);


  gkyl_proj_on_basis_advance(proj_density, 0.0, &conf_local, density);
  gkyl_proj_on_basis_advance(proj_drift, 0.0, &conf_local, drift);
  gkyl_proj_on_basis_advance(proj_Tpar, 0.0, &conf_local, Tpar);
  gkyl_proj_on_basis_advance(proj_Tperp, 0.0, &conf_local, Tperp);

  gkyl_array_set_offset_range(moms, 1.0, density, 0*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms, 1.0, drift, 1*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms, 1.0, Tpar, (vdim+1)*conf_basis.num_basis, &conf_local);
  gkyl_array_set_offset_range(moms, 1.0, Tperp, (vdim+2)*conf_basis.num_basis, &conf_local);

  struct gkyl_dist_proj_type *dist_proj = gkyl_vlasov_bimax_new(cdim, pdim, poly_order,
    conf_basis.num_basis, phase_basis.num_basis, GKYL_MODEL_DEFAULT);

  TEST_CHECK(dist_proj != NULL);
  TEST_CHECK(dist_proj->dist_id == GKYL_DIST_TYPE_BIMAX);
  TEST_CHECK(gkyl_dist_proj_type_num_mom(dist_proj) == vdim+3);

  // Kernel test
  double xc_test[] = {0.5, 0.5, -0.25};
  double dx_test[] = {phase_grid.dx[0], phase_grid.dx[1], phase_grid.dx[2]};
  int idx_test[] = {1, 1, 1};
  double moms_test[] = {
    1.0,     // n
    0.5,     // V_parallel
    -0.25,   // V_perpendicular
    1.0,     // T_parallel/m
    2.0      // T_perpendicular/m
  };
  double f_test[1];
  gkyl_dist_proj_type_calc(dist_proj, xc_test, dx_test, idx_test, moms_test, f_test);

  double f_expected = 1.0/(sqrt(2.0*GKYL_PI*1.0)*sqrt(2.0*GKYL_PI*2.0));
  TEST_CHECK(gkyl_compare_double(f_test[0], f_expected, 1e-14));

  xc_test[1] = 1.5;
  xc_test[2] = 0.75;
  gkyl_dist_proj_type_calc(dist_proj, xc_test, dx_test, idx_test, moms_test, f_test);

  f_expected = 1.0/(sqrt(2.0*GKYL_PI*1.0)*sqrt(2.0*GKYL_PI*2.0))*exp(-0.75);
  TEST_CHECK(gkyl_compare_double(f_test[0], f_expected, 1e-14));

  // projection updater to compute the distribution

  struct gkyl_vlasov_dist_proj_on_basis_inp inp = {
    .dist_proj = dist_proj,
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &conf_basis,
    .vel_basis = &vel_basis,
    .phase_basis = &phase_basis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_local_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .model_id = GKYL_MODEL_DEFAULT,
    .use_gpu = false,
    .quad_type = GKYL_GAUSS_QUAD,
  };

  gkyl_vlasov_dist_proj_on_basis *proj_dist = gkyl_vlasov_dist_proj_on_basis_inew(&inp);
  TEST_CHECK(proj_dist != NULL);
  struct gkyl_array *distf = mkarr(phase_basis.num_basis, phase_local_ext.volume);

  gkyl_vlasov_dist_proj_on_basis_moments_advance(proj_dist, &conf_local, moms);
  gkyl_vlasov_dist_proj_on_basis_advance(proj_dist, &phase_local, &conf_local, distf);

  struct bimax_test_ctx bimax_ctx = {
    .n = 1.0,
    .drift_par = 0.5,
    .drift_perp = -0.25,
    .Tpar_over_m = 1.0,
    .Tperp_over_m = 2.0,
  };

  struct gkyl_array *distf_ref = mkarr(phase_basis.num_basis, phase_local_ext.volume);

  gkyl_proj_on_basis *proj_ref = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    poly_order+1, 1, eval_bimax_1x2v, &bimax_ctx);

  gkyl_proj_on_basis_advance(proj_ref, 0.0, &phase_local, distf_ref);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &phase_local);

  while (gkyl_range_iter_next(&iter)) {

    long idx = gkyl_range_idx(&phase_local_ext, iter.idx);
    const double *f = gkyl_array_cfetch(distf, idx);
    const double *f_ref = gkyl_array_cfetch(distf_ref, idx);

    for (int k=0; k<phase_basis.num_basis; ++k) {
      TEST_CHECK(gkyl_compare_double(f[k], f_ref[k],1e-12));
    }
  }
  char fname[1024];
  sprintf(fname,"ctest_dist_proj_on_basis_bimax_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&phase_grid, &phase_local, 0, distf, fname);

  gkyl_array_release(density);
  gkyl_array_release(drift);
  gkyl_array_release(Tpar);
  gkyl_array_release(Tperp);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ref);
  gkyl_proj_on_basis_release(proj_density);
  gkyl_proj_on_basis_release(proj_drift);
  gkyl_proj_on_basis_release(proj_Tpar);
  gkyl_proj_on_basis_release(proj_Tperp);
  gkyl_proj_on_basis_release(proj_ref);
  gkyl_vlasov_dist_proj_on_basis_release(proj_dist);
  gkyl_dist_proj_type_release(dist_proj);  
}

void
test_1x1v_lte_p2(void) {test_1x1v_lte(2);}

void
test_1x2v_bimax_p2(void) {test_1x2v_bimax(2);}

TEST_LIST = {
  {"test_1x1v_lte_p2", test_1x1v_lte_p2},
  {"test_1x2v_bimax_p2", test_1x2v_bimax_p2},
  {NULL, NULL},
};