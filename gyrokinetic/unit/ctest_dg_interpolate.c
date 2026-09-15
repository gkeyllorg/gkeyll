#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_array_integrate.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_ctx {
  double n0; // Density.
  double ux0, uy0, uz0; // Mean flow speed along x,y,z.
  double upar0; // Parallel flow speed.
  double T0; // Temperature.
  double mass; // Particle mass.
  double B0; // Magnetic field.
  int cdim; // Number of position space dimensions.
  int vdim; // Number of velocity space dimensions.
  double lower[GKYL_MAX_DIM]; // Grid lower limit in each direction.
  double upper[GKYL_MAX_DIM]; // Grid upper limit in each direction.
};

static void calc_int_moms(int num_mom, struct gkyl_rect_grid *confGrid, struct gkyl_basis *confBasis,
  struct gkyl_range *confLocal, bool use_gpu, struct gkyl_array *moms, double *int_moms)
{
  // Compute the volume integral of the moments.
  double *integrated_moms = use_gpu? gkyl_cu_malloc(num_mom*sizeof(double)) : gkyl_malloc(num_mom*sizeof(double));

  struct gkyl_array_integrate* integ_op = gkyl_array_integrate_new(confGrid, confBasis,
    num_mom, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);

  gkyl_array_integrate_advance(integ_op, moms, 1.0, moms, confLocal, confLocal, integrated_moms); 

  if (use_gpu)
    gkyl_cu_memcpy(int_moms, integrated_moms, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(int_moms, integrated_moms, sizeof(double[num_mom]));

  if (use_gpu)
    gkyl_cu_free(integrated_moms);
  else
    gkyl_free(integrated_moms);
  gkyl_array_integrate_release(integ_op);
}

void eval_fdonor_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  int cdim = tctx->cdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  fout[0] = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x));
}

void
test_1x(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  double lower[] = {0.0}, upper[] = {1.0};
  double mass = 1.0;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar0 = 1.2, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .mass = mass, // Particle mass.
    .cdim = ndim, // Number of position space dimensions.
    .lower = {lower[0]}, // Lower extents of the grid.
    .upper = {upper[0]}, // Upper extents of the grid.
  };

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create donot field.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_fdonor_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  int num_mom = 1;
  // Calculate the integral of the field.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &grid, &basis, &local, use_gpu, distf, int_moms);

//  // Write donor field to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_1x_p%d_N%d-N%d_do.gkyl", poly_order, cells[0], cells_tar[0]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  printf("\n  Donor int_moms  = %g\n", int_moms[0]);

  // Target grid.
  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);

  // Target range.
  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(ndim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &grid_tar, &basis, &local_tar, use_gpu, distf_tar, int_moms_tar);

//  // Write target field to file.
//  gkyl_array_copy(distf_tar_ho, distf_tar);
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_1x_p%d_N%d-N%d_tar.gkyl", poly_order, cells[0], cells_tar[0]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  printf("\n  Target int_moms = %g\n", int_moms_tar[0]);

  // Check that the moments and integrated moments are the same.
  for (int i=0; i<num_mom; i++) {
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );
    TEST_MSG( "Got: %g | Expected: %g\n", int_moms_tar[i], int_moms[i] );
  }

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_distf);
}

void eval_fdonor_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  int cdim = tctx->cdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  fout[0] = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x))*exp(-pow(y,2)/(2.0*pow(M_PI/3.0,2)));
}

void
test_2x(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  double lower[] = {0.0, -M_PI}, upper[] = {1.0, M_PI};
  double mass = 1.0;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar0 = 1.2, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .mass = mass, // Particle mass.
    .cdim = ndim, // Number of position space dimensions.
    .lower = {lower[0], lower[1]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1]}, // Upper extents of the grid.
  };

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create donot field.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_fdonor_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  int num_mom = 1;
  // Calculate the integral of the field.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &grid, &basis, &local, use_gpu, distf, int_moms);

//  // Write donor field to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_2x_p%d_N%dx%d-N%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  printf("\n  Donor int_moms  = %g\n", int_moms[0]);

  // Target grid.
  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);

  // Target range.
  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(ndim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &grid_tar, &basis, &local_tar, use_gpu, distf_tar, int_moms_tar);

//  // Write target field to file.
//  gkyl_array_copy(distf_tar_ho, distf_tar);
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_2x_p%d_N%dx%d-N%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  printf("\n  Target int_moms = %g\n", int_moms_tar[0]);

  // Check that the moments and integrated moments are the same.
  for (int i=0; i<num_mom; i++) {
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );
    TEST_MSG( "Got: %g | Expected: %g\n", int_moms_tar[i], int_moms[i] );
  }

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_distf);
}

void eval_distf_1x1v_vlasov(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double ux0 = tctx->ux0;
  double T0 = tctx->T0;
  double mass = tctx->mass;
  int cdim = tctx->cdim;
  int vdim = tctx->vdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  double den = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x));
  double ux = ux0;
  double temp = T0;

  double vtsq = temp/mass;

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vx-ux,2))/(2.0*vtsq));
}

static void calc_moms_vlasov(struct gkyl_rect_grid *grid, struct gkyl_rect_grid *velGrid,
  struct gkyl_basis *confBasis, struct gkyl_basis *basis,
  struct gkyl_range *confLocal, struct gkyl_range *local, struct gkyl_range *velLocal,
  bool use_gpu, struct gkyl_array *distf, struct gkyl_array *moms)
{
  int cdim = confBasis->ndim, pdim = basis->ndim, vdim = pdim - cdim;
  int poly_order = confBasis->poly_order;
  bool is_hybrid = basis->b_type == GKYL_BASIS_MODAL_HYBRID;

  // The p=1 hybrid pairs a p=1 tensor configuration-space basis with a p=2
  // tensor velocity-space basis; kernel selection keys on the conf basis type.
  struct gkyl_basis vbasis, momConfBasis;
  if (is_hybrid) {
    gkyl_cart_modal_tensor(&vbasis, vdim, 2);
    gkyl_cart_modal_tensor(&momConfBasis, cdim, poly_order);
  }
  else {
    gkyl_cart_modal_serendip(&vbasis, vdim, poly_order);
    momConfBasis = *confBasis;
  }

  // Identity velocity map and velocity-space Hamiltonian H = v^2/2.
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(velGrid,
    velLocal, &vbasis, inp_vmap, false, use_gpu);
  struct gkyl_array *hamil = mkarr(use_gpu, vbasis.num_basis, velLocal->volume);
  struct gkyl_array *gamma_inv = mkarr(use_gpu, vbasis.num_basis, velLocal->volume);
  gkyl_dg_vlasov_calc_hamil(velGrid, &vbasis, velLocal,
    GKYL_MODEL_DEFAULT, vvm, hamil, gamma_inv, use_gpu);

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &momConfBasis,
    .phase_basis = basis,
    .vel_range = velLocal,
    .hamil_range = velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .hamil_id = gkyl_hamil_id_from_model_id(GKYL_MODEL_DEFAULT),
    .mom_type = GKYL_F_MOMENT_M0M1M2,
    .use_gpu = use_gpu,
    .vel_map = vvm,
  };
  struct gkyl_mom_type *mtype = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *mcalc = gkyl_mom_calc_new(grid, mtype, use_gpu);
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    gkyl_mom_calc_advance_cu(mcalc, local, confLocal, distf, moms);
  else
    gkyl_mom_calc_advance(mcalc, local, confLocal, distf, moms);
#else
  gkyl_mom_calc_advance(mcalc, local, confLocal, distf, moms);
#endif
  gkyl_mom_calc_release(mcalc);
  gkyl_mom_type_release(mtype);
  gkyl_vlasov_velocity_map_release(vvm);
  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
}

void
test_1x1v_vlasov(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  int cdim = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double vx_min = -6.0;
  double vx_max =  6.0;
  double lower[] = {x_min, vx_min}, upper[] = {x_max, vx_max};
  double mass = 1.0;

  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .ux0 = 1.2, // Flow speed along x.
    .T0 = 1.0, // Temperature.
    .mass = mass, // Particle mass.
    .cdim = cdim, // Number of position space dimensions.
    .vdim = vdim, // Number of velocity space dimensions.
    .lower = {lower[0], lower[1]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1]}, // Upper extents of the grid.
  };

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x1v_vlasov, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  int num_mom = 2+vdim;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_ho = use_gpu? mkarr(false, moms->ncomp, moms->size)
                                      : gkyl_array_acquire(moms);
  calc_moms_vlasov(&grid, &velGrid, &confBasis, &basis, &confLocal, &local, &velLocal, use_gpu, distf, moms);
  gkyl_array_copy(moms_ho, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

//  // Write donor distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  // Write target moments to file.
//  char fname0m[1024];
//  sprintf(fname0m, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms_ho, fname0m);
//
//  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  struct gkyl_array *moms_tar_ho = use_gpu? mkarr(false, moms_tar->ncomp, moms_tar->size)
                                          : gkyl_array_acquire(moms_tar);
  calc_moms_vlasov(&grid_tar, &velGrid_tar, &confBasis, &basis, &confLocal_tar, &local_tar, &velLocal_tar, use_gpu, distf_tar, moms_tar);
  gkyl_array_copy(moms_tar_ho, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

//  // Write target distribution function to file.
//  gkyl_array_copy(distf_tar_ho, distf_tar);
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  // Write target moments to file.
//  char fname1m[1024];
//  sprintf(fname1m, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar_ho, fname1m);
//
//  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check that the moments and integrated moments are the same.
  bool same_conf_grid = true;
  for (int d=0; d<cdim; d++) same_conf_grid = same_conf_grid && (cells[d] == cells_tar[d]);
  if (same_conf_grid) {
    double m2_tol = 1e-10;
    if (poly_order == 2 && basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      m2_tol = 1e-2; // Because this is not a tensor basis.
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &confLocal);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&confLocal, iter.idx);
      const double *moms_c = gkyl_array_cfetch(moms_ho, linidx);
      const double *moms_tar_c = gkyl_array_cfetch(moms_tar_ho, linidx);
      for (int m=0; m<2*confBasis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], 1e-10) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
      for (int m=2*confBasis.num_basis; m<moms->ncomp; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], m2_tol) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
    }
  }

  for (int i=0; i<num_mom; i++) {
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );
    TEST_MSG( "Got: %g | Expected: %g\n", int_moms_tar[i], int_moms[i] );
  }

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar_ho);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_distf);
}

void eval_distf_1x2v_vlasov(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1], vy = xn[2];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double ux0 = tctx->ux0;
  double uy0 = tctx->uy0;
  double T0 = tctx->T0;
  double mass = tctx->mass;
  int cdim = tctx->cdim;
  int vdim = tctx->vdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  double den = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x));
  double ux = ux0;
  double uy = uy0;
  double temp = T0;

  double vtsq = temp/mass;

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vx-ux,2)+pow(vy-uy,2))/(2.0*vtsq));
}

void
test_1x2v_vlasov(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  int cdim = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double vx_min = -6.0;
  double vx_max =  6.0;
  double vy_min =  -0.5;
  double vy_max =   0.5;
  double lower[] = {x_min, vx_min, vy_min}, upper[] = {x_max, vx_max, vy_max};
  double mass = 1.0;

  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .ux0 = 0.0, // Flow speed along x.
    .uy0 = 1.2, // Flow speed along y.
    .T0 = 1.0, // Temperature.
    .mass = mass, // Particle mass.
    .cdim = cdim, // Number of position space dimensions.
    .vdim = vdim, // Number of velocity space dimensions.
    .lower = {lower[0], lower[1], lower[2]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1], upper[2]}, // Upper extents of the grid.
  };


  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x2v_vlasov, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  int num_mom = 2+vdim;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_ho = use_gpu? mkarr(false, moms->ncomp, moms->size)
                                      : gkyl_array_acquire(moms);
  calc_moms_vlasov(&grid, &velGrid, &confBasis, &basis, &confLocal, &local, &velLocal, use_gpu, distf, moms);
  gkyl_array_copy(moms_ho, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

//  // Write donor distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  // Write target moments to file.
//  char fname0m[1024];
//  sprintf(fname0m, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms_ho, fname0m);
//
//  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  struct gkyl_array *moms_tar_ho = use_gpu? mkarr(false, moms_tar->ncomp, moms_tar->size)
                                          : gkyl_array_acquire(moms_tar);
  calc_moms_vlasov(&grid_tar, &velGrid_tar, &confBasis, &basis, &confLocal_tar, &local_tar, &velLocal_tar, use_gpu, distf_tar, moms_tar);
  gkyl_array_copy(moms_tar_ho, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

//  // Write target distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  // Write target moments to file.
//  char fname1m[1024];
//  sprintf(fname1m, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar_ho, fname1m);
//
//  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check that the moments and integrated moments are the same.
  bool same_conf_grid = true;
  for (int d=0; d<cdim; d++) same_conf_grid = same_conf_grid && (cells[d] == cells_tar[d]);
  if (same_conf_grid) {
    double m2_tol = 1e-10;
    if (poly_order == 2 && basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      m2_tol = 1e-2; // Because this is not a tensor basis.
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &confLocal);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&confLocal, iter.idx);
      const double *moms_c = gkyl_array_cfetch(moms_ho, linidx);
      const double *moms_tar_c = gkyl_array_cfetch(moms_tar_ho, linidx);
      for (int m=0; m<2*confBasis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], 1e-10) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
      for (int m=2*confBasis.num_basis; m<moms->ncomp; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], m2_tol) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
    }
  }

  for (int i=0; i<num_mom; i++)
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar_ho);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_distf);
}

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void eval_bfield_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
}

void eval_bfield_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
}

void eval_bfield_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
}

void eval_distf_1x1v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double upar0 = tctx->upar0;
  double T0 = tctx->T0;
  double mass = tctx->mass;
  int cdim = tctx->cdim;
  int vdim = tctx->vdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  double den = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x));
  double upar = upar0;
  double temp = T0;

  double vtsq = temp/mass;

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2))/(2.0*vtsq));
}

static struct gk_geometry* init_gk_geo(int poly_order, struct gkyl_rect_grid confGrid, struct gkyl_basis confBasis,
  struct gkyl_range confLocal, struct gkyl_range confLocal_ext, void *bfield_ctx, bool use_gpu)
{
  // Initialize GK geometry.
  int cdim = confBasis.ndim;
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0, 0.0, 0.0},  .mapc2p = mapc2p,  .c2p_ctx = 0,
    .bfield_func = cdim==1? eval_bfield_1x : (cdim==2? eval_bfield_2x : eval_bfield_3x),  .bfield_ctx = bfield_ctx,
    .basis = confBasis,  .grid = confGrid,
    .local = confLocal,  .local_ext = confLocal_ext,
    .global = confLocal, .global_ext = confLocal_ext,
  };
  int geo_ghost[3] = {1, 1, 1};
  if (cdim < 3) {
    geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
    gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
    gkyl_create_grid_ranges(&geometry_input.geo_grid, geo_ghost, &geometry_input.geo_global_ext, &geometry_input.geo_global);
    memcpy(&geometry_input.geo_local, &geometry_input.geo_global, sizeof(struct gkyl_range));
    memcpy(&geometry_input.geo_local_ext, &geometry_input.geo_global_ext, sizeof(struct gkyl_range));
  }
  else {
    geometry_input.geo_grid = confGrid;
    geometry_input.geo_basis = confBasis;
    geometry_input.geo_global = confLocal;
    geometry_input.geo_global_ext = confLocal_ext;
    geometry_input.geo_local = confLocal;
    geometry_input.geo_local_ext = confLocal_ext;
  }
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();
  geometry_input.position_map = pmap;
  
  // Deflate geometry.
  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  struct gk_geometry *gk_geom = cdim < 3? gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input)
	                                : gkyl_gk_geometry_acquire(gk_geom_3d);
  gkyl_gk_geometry_release(gk_geom_3d);
  // If we are on the gpu, copy from host.
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }
  gkyl_position_map_release(pmap);
  return gk_geom;
}

static void calc_moms_gk(struct gkyl_rect_grid *grid, struct gkyl_basis *confBasis, struct gkyl_basis *basis,
  struct gkyl_range *confLocal, struct gkyl_range *local, double mass, double charge, struct gkyl_velocity_map *gvm,
  struct gk_geometry *gk_geom, bool use_gpu, struct gkyl_array *distf, struct gkyl_array *moms)
{
  struct gkyl_dg_updater_moment* mom_op = gkyl_dg_updater_moment_gyrokinetic_new(grid, confBasis,
    basis, confLocal, mass, charge, gvm, gk_geom, 0, GKYL_F_MOMENT_M0M1M2, false, use_gpu);
  gkyl_dg_updater_moment_gyrokinetic_advance(mom_op, local, confLocal, distf, moms);
  gkyl_dg_updater_moment_gyrokinetic_release(mom_op);
}

void
test_1x1v_gk(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  int cdim = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double vpar_min = -6.0;
  double vpar_max =  6.0;
  double lower[] = {x_min, vpar_min}, upper[] = {x_max, vpar_max};
  double mass = 1.0;
  double charge = 1.0;

  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar0 = 1.2, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .mass = mass, // Particle mass.
    .cdim = cdim, // Number of position space dimensions.
    .vdim = vdim, // Number of velocity space dimensions.
    .lower = {lower[0], lower[1]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1]}, // Upper extents of the grid.
  };

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                                      : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bfield_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x1v_gk, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize geometry.
  struct gk_geometry *gk_geom = init_gk_geo(poly_order, confGrid, confBasis, confLocal, confLocal_ext, &proj_ctx, use_gpu);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  int num_mom = 3;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_ho = use_gpu? mkarr(false, moms->ncomp, moms->size)
                                      : gkyl_array_acquire(moms);
  calc_moms_gk(&grid, &confBasis, &basis, &confLocal, &local, mass, charge, gvm, gk_geom, use_gpu, distf, moms);
  gkyl_array_copy(moms_ho, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

//  // Write donor distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  // Write target moments to file.
//  char fname0m[1024];
//  sprintf(fname0m, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms_ho, fname0m);
//
//  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target geometry.
  struct gk_geometry *gk_geom_tar = init_gk_geo(poly_order, confGrid_tar, confBasis,
    confLocal_tar, confLocal_tar_ext, &proj_ctx, use_gpu);

  // Target velocity space mapping.
  struct gkyl_velocity_map *gvm_tar = gkyl_velocity_map_new(c2p_in, grid_tar, velGrid_tar,
    local_tar, local_tar_ext, velLocal_tar, velLocal_tar_ext, use_gpu);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  struct gkyl_array *moms_tar_ho = use_gpu? mkarr(false, moms_tar->ncomp, moms_tar->size)
                                          : gkyl_array_acquire(moms_tar);
  calc_moms_gk(&grid_tar, &confBasis, &basis, &confLocal_tar, &local_tar,
    mass, charge, gvm_tar, gk_geom_tar, use_gpu, distf_tar, moms_tar);
  gkyl_array_copy(moms_tar_ho, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

//  // Write target distribution function to file.
//  gkyl_array_copy(distf_tar_ho, distf_tar);
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  // Write target moments to file.
//  char fname1m[1024];
//  sprintf(fname1m, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
//  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar_ho, fname1m);
//
//  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check that the moments and integrated moments are the same.
  bool same_conf_grid = true;
  for (int d=0; d<cdim; d++) same_conf_grid = same_conf_grid && (cells[d] == cells_tar[d]);
  if (same_conf_grid) {
    double m2_tol = 1e-10;
    if (poly_order == 2 && basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      m2_tol = 1e-2; // Because this is not a tensor basis.
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &confLocal);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&confLocal, iter.idx);
      const double *moms_c = gkyl_array_cfetch(moms_ho, linidx);
      const double *moms_tar_c = gkyl_array_cfetch(moms_tar_ho, linidx);
      for (int m=0; m<2*confBasis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], 1e-10) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
      for (int m=2*confBasis.num_basis; m<moms->ncomp; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], m2_tol) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
    }
  }

  for (int i=0; i<num_mom; i++) {
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );
    TEST_MSG( "Got: %g | Expected: %g\n", int_moms_tar[i], int_moms[i] );
  }

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar_ho);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm_tar);
  gkyl_gk_geometry_release(gk_geom_tar);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void eval_distf_1x2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1], mu = xn[2];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double upar0 = tctx->upar0;
  double T0 = tctx->T0;
  double mass = tctx->mass;
  int cdim = tctx->cdim;
  int vdim = tctx->vdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  double den = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x));
  double upar = upar0;
  double temp = T0;

  double vtsq = temp/mass;

  double bmag[1] = {-1.0};
  eval_bfield_1x(t, xn, bmag, ctx);

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2)+2.0*mu*bmag[0]/mass)/(2.0*vtsq));
}

void
test_1x2v_gk(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  int cdim = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double vpar_min = -6.0;
  double vpar_max =  6.0;
  double mu_max =  0.5;
  double lower[] = {x_min, vpar_min, 0.0}, upper[] = {x_max, vpar_max, mu_max};
  double mass = 1.0;
  double charge = 1.0;

  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar0 = 1.2, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .mass = mass, // Particle mass.
    .cdim = cdim, // Number of position space dimensions.
    .vdim = vdim, // Number of velocity space dimensions.
    .lower = {lower[0], lower[1], lower[2]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1], upper[2]}, // Upper extents of the grid.
  };


  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                                      : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bfield_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x2v_gk, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize geometry.
  struct gk_geometry *gk_geom = init_gk_geo(poly_order, confGrid, confBasis, confLocal, confLocal_ext, &proj_ctx, use_gpu);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  int num_mom = 3;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_ho = use_gpu? mkarr(false, moms->ncomp, moms->size)
                                      : gkyl_array_acquire(moms);
  calc_moms_gk(&grid, &confBasis, &basis, &confLocal, &local, mass, charge, gvm, gk_geom, use_gpu, distf, moms);
  gkyl_array_copy(moms_ho, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

//  // Write donor distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  // Write target moments to file.
//  char fname0m[1024];
//  sprintf(fname0m, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms_ho, fname0m);
//
//  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target geometry.
  struct gk_geometry *gk_geom_tar = init_gk_geo(poly_order, confGrid_tar, confBasis,
    confLocal_tar, confLocal_tar_ext, &proj_ctx, use_gpu);

  // Target velocity space mapping.
  struct gkyl_velocity_map *gvm_tar = gkyl_velocity_map_new(c2p_in, grid_tar, velGrid_tar,
    local_tar, local_tar_ext, velLocal_tar, velLocal_tar_ext, use_gpu);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  struct gkyl_array *moms_tar_ho = use_gpu? mkarr(false, moms_tar->ncomp, moms_tar->size)
                                          : gkyl_array_acquire(moms_tar);
  calc_moms_gk(&grid_tar, &confBasis, &basis, &confLocal_tar, &local_tar,
    mass, charge, gvm_tar, gk_geom_tar, use_gpu, distf_tar, moms_tar);
  gkyl_array_copy(moms_tar_ho, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

//  // Write target distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  // Write target moments to file.
//  char fname1m[1024];
//  sprintf(fname1m, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
//  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar_ho, fname1m);
//
//  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check that the moments and integrated moments are the same.
  bool same_conf_grid = true;
  for (int d=0; d<cdim; d++) same_conf_grid = same_conf_grid && (cells[d] == cells_tar[d]);
  if (same_conf_grid) {
    double m2_tol = 1e-10;
    if (poly_order == 2 && basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      m2_tol = 1e-2; // Because this is not a tensor basis.
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &confLocal);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&confLocal, iter.idx);
      const double *moms_c = gkyl_array_cfetch(moms_ho, linidx);
      const double *moms_tar_c = gkyl_array_cfetch(moms_tar_ho, linidx);
      for (int m=0; m<2*confBasis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], 1e-10) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
      for (int m=2*confBasis.num_basis; m<moms->ncomp; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], m2_tol) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
    }
  }

  for (int i=0; i<num_mom; i++)
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar_ho);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm_tar);
  gkyl_gk_geometry_release(gk_geom_tar);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void eval_distf_2x2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], vpar = xn[2], mu = xn[3];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double upar0 = tctx->upar0;
  double T0 = tctx->T0;
  double mass = tctx->mass;
  int cdim = tctx->cdim;
  int vdim = tctx->vdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  double den = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x)*cos((2*2.0*M_PI/Lx[1])*y));
  double upar = upar0;
  double temp = T0;

  double vtsq = temp/mass;

  double bmag[1] = {-1.0};
  eval_bfield_2x(t, xn, bmag, ctx);

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2)+2.0*mu*bmag[0]/mass)/(2.0*vtsq));
}

void
test_2x2v_gk(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  int cdim = 2;
//  double x_min = 0.0;
//  double x_max = 1.0;
//  double y_min = 0.0;
//  double y_max = 1.0;
//  double vpar_min = -6.0;
//  double vpar_max =  6.0;
//  double mu_max =  0.5;
  double x_min = 9.460546e-01;
  double x_max = 1.053945e+00;
  double y_min = -1.078908e-01;
  double y_max = 1.078908e-01;
  double vpar_min = -1.744683e+05;
  double vpar_max =  1.744683e+05;
  double mu_max =  9.047585e-17;
  double lower[] = {x_min, y_min, vpar_min, 0.0}, upper[] = {x_max, y_max, vpar_max, mu_max};
  double mass = 1.67e-27*3.973;
  double charge = 1.602e-19;

  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar0 = 1.2, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .mass = mass, // Particle mass.
    .cdim = cdim, // Number of position space dimensions.
    .vdim = vdim, // Number of velocity space dimensions.
    .lower = {lower[0], lower[1], lower[2], lower[3]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1], upper[2], upper[3]}, // Upper extents of the grid.
  };


  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                                      : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bfield_2x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_2x2v_gk, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize geometry.
  struct gk_geometry *gk_geom = init_gk_geo(poly_order, confGrid, confBasis, confLocal, confLocal_ext, &proj_ctx, use_gpu);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  int num_mom = 3;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_ho = use_gpu? mkarr(false, moms->ncomp, moms->size)
                                      : gkyl_array_acquire(moms);
  calc_moms_gk(&grid, &confBasis, &basis, &confLocal, &local, mass, charge, gvm, gk_geom, use_gpu, distf, moms);
  gkyl_array_copy(moms_ho, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

//  // Write donor distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_2x2v_p%d_N%dx%dx%dx%d-N%dx%dx%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  // Write target moments to file.
//  char fname0m[1024];
//  sprintf(fname0m, "ctest_dg_interp_2x2v_p%d_N%dx%dx%dx%d-N%dx%dx%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3]);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms_ho, fname0m);
//
//  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target geometry.
  struct gk_geometry *gk_geom_tar = init_gk_geo(poly_order, confGrid_tar, confBasis,
    confLocal_tar, confLocal_tar_ext, &proj_ctx, use_gpu);

  // Target velocity space mapping.
  struct gkyl_velocity_map *gvm_tar = gkyl_velocity_map_new(c2p_in, grid_tar, velGrid_tar,
    local_tar, local_tar_ext, velLocal_tar, velLocal_tar_ext, use_gpu);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  struct gkyl_array *moms_tar_ho = use_gpu? mkarr(false, moms_tar->ncomp, moms_tar->size)
                                          : gkyl_array_acquire(moms_tar);
  calc_moms_gk(&grid_tar, &confBasis, &basis, &confLocal_tar, &local_tar,
    mass, charge, gvm_tar, gk_geom_tar, use_gpu, distf_tar, moms_tar);
  gkyl_array_copy(moms_tar_ho, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

//  // Write target distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_2x2v_p%d_N%dx%dx%dx%d-N%dx%dx%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  // Write target moments to file.
//  char fname1m[1024];
//  sprintf(fname1m, "ctest_dg_interp_2x2v_p%d_N%dx%dx%dx%d-N%dx%dx%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3]);
//  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar_ho, fname1m);
//
//  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check that the moments and integrated moments are the same.
  bool same_conf_grid = true;
  for (int d=0; d<cdim; d++) same_conf_grid = same_conf_grid && (cells[d] == cells_tar[d]);
  if (same_conf_grid) {
    double m2_tol = 1e-10;
    if (poly_order == 2 && basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      m2_tol = 1e-2; // Because this is not a tensor basis.
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &confLocal);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&confLocal, iter.idx);
      const double *moms_c = gkyl_array_cfetch(moms_ho, linidx);
      const double *moms_tar_c = gkyl_array_cfetch(moms_tar_ho, linidx);
      for (int m=0; m<2*confBasis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], 1e-10) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
      for (int m=2*confBasis.num_basis; m<moms->ncomp; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], m2_tol) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
    }
  }

  for (int i=0; i<num_mom; i++)
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar_ho);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm_tar);
  gkyl_gk_geometry_release(gk_geom_tar);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void eval_distf_3x2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vpar = xn[3], mu = xn[4];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double upar0 = tctx->upar0;
  double T0 = tctx->T0;
  double mass = tctx->mass;
  int cdim = tctx->cdim;
  int vdim = tctx->vdim;
  double *lower = tctx->lower;
  double *upper = tctx->upper;

  double Lx[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  double den = n0*(1.0+0.5*sin((2.0*M_PI/Lx[0])*x)*cos((2*2.0*M_PI/Lx[1])*y));
  double upar = upar0;
  double temp = T0;

  double vtsq = temp/mass;

  double bmag[1] = {-1.0};
  eval_bfield_3x(t, xn, bmag, ctx);

  fout[0] = (den/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2)+2.0*mu*bmag[0]/mass)/(2.0*vtsq));
}

void
test_3x2v_gk(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  int cdim = 3;
//  double x_min = 0.0;
//  double x_max = 1.0;
//  double y_min = 0.0;
//  double y_max = 1.0;
//  double z_min = 0.0;
//  double z_max = 1.0;
//  double vpar_min = -6.0;
//  double vpar_max =  6.0;
//  double mu_max =  0.5;
//  double x_min = 9.460546e-01;
//  double x_max = 1.053945e+00;
  double y_min = -1.078908e-01;
  double y_max = 1.078908e-01;
  double z_min = -2.000000e+00;
  double z_max = 2.000000e+00;
  double vpar_min = -1.06096e+07;
  double vpar_max =  1.06096e+07;
  double mu_max =  9.047585e-17;
  double x_min = 9.460545916e-01;
  double x_max = 1.053945408e+00;
//  double y_min = -1.078908169e-01;
//  double y_max = 1.078908169e-01;
//  double z_min = -2.000000e+00;
//  double z_max = 2.000000e+00;
//  double vpar_min = -1.060964135e+07;
//  double vpar_max =  1.060964135e+07;
//  double mu_max =  9.047584868e-17;
  double lower[] = {x_min, y_min, z_min, vpar_min, 0.0}, upper[] = {x_max, y_max, z_max, vpar_max, mu_max};
  double mass = 1.67e-27*3.973;
  double charge = 1.602e-19;

  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = ndim-cdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar0 = 1.2, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .mass = mass, // Particle mass.
    .cdim = cdim, // Number of position space dimensions.
    .vdim = vdim, // Number of velocity space dimensions.
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]}, // Lower extents of the grid.
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]}, // Upper extents of the grid.
  };


  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<cdim; d++) confGhost[d] = 1;
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                                      : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bfield_3x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_3x2v_gk, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize geometry.
  struct gk_geometry *gk_geom = init_gk_geo(poly_order, confGrid, confBasis, confLocal, confLocal_ext, &proj_ctx, use_gpu);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  int num_mom = 3;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_ho = use_gpu? mkarr(false, moms->ncomp, moms->size)
                                      : gkyl_array_acquire(moms);
  calc_moms_gk(&grid, &confBasis, &basis, &confLocal, &local, mass, charge, gvm, gk_geom, use_gpu, distf, moms);
  gkyl_array_copy(moms_ho, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

//  // Write donor distribution function to file.
//  char fname0[1024];
//  sprintf(fname0, "ctest_dg_interp_3x2v_p%d_N%dx%dx%dx%dx%d-N%dx%dx%dx%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells[4], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3], cells_tar[4]);
//  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);
//
//  // Write target moments to file.
//  char fname0m[1024];
//  sprintf(fname0m, "ctest_dg_interp_3x2v_p%d_N%dx%dx%dx%dx%d-N%dx%dx%dx%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells[4], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3], cells_tar[4]);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms_ho, fname0m);
//
//  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target geometry.
  struct gk_geometry *gk_geom_tar = init_gk_geo(poly_order, confGrid_tar, confBasis,
    confLocal_tar, confLocal_tar_ext, &proj_ctx, use_gpu);

  // Target velocity space mapping.
  struct gkyl_velocity_map *gvm_tar = gkyl_velocity_map_new(c2p_in, grid_tar, velGrid_tar,
    local_tar, local_tar_ext, velLocal_tar, velLocal_tar_ext, use_gpu);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, &local, &local_tar, ghost, use_gpu);

  gkyl_dg_interpolate_advance(interp, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  struct gkyl_array *moms_tar_ho = use_gpu? mkarr(false, moms_tar->ncomp, moms_tar->size)
                                          : gkyl_array_acquire(moms_tar);
  calc_moms_gk(&grid_tar, &confBasis, &basis, &confLocal_tar, &local_tar,
    mass, charge, gvm_tar, gk_geom_tar, use_gpu, distf_tar, moms_tar);
  gkyl_array_copy(moms_tar_ho, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

//  // Write target distribution function to file.
//  char fname1[1024];
//  sprintf(fname1, "ctest_dg_interp_3x2v_p%d_N%dx%dx%dx%dx%d-N%dx%dx%dx%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells[4], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3], cells_tar[4]);
//  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);
//
//  // Write target moments to file.
//  char fname1m[1024];
//  sprintf(fname1m, "ctest_dg_interp_3x2v_p%d_N%dx%dx%dx%dx%d-N%dx%dx%dx%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells[3], cells[4], cells_tar[0], cells_tar[1], cells_tar[2], cells_tar[3], cells_tar[4]);
//  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar_ho, fname1m);
//
//  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check that the moments and integrated moments are the same.
  bool same_conf_grid = true;
  for (int d=0; d<cdim; d++) same_conf_grid = same_conf_grid && (cells[d] == cells_tar[d]);
  if (same_conf_grid) {
    double m2_tol = 1e-10;
    if (poly_order == 2 && basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      m2_tol = 1e-2; // Because this is not a tensor basis.
    }

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &confLocal);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&confLocal, iter.idx);
      const double *moms_c = gkyl_array_cfetch(moms_ho, linidx);
      const double *moms_tar_c = gkyl_array_cfetch(moms_tar_ho, linidx);
      for (int m=0; m<2*confBasis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], 1e-10) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
      for (int m=2*confBasis.num_basis; m<moms->ncomp; m++) {
        TEST_CHECK( gkyl_compare(moms_c[m], moms_tar_c[m], m2_tol) );
        TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, moms_tar_c[m], moms_c[m]);
      }
    }
  }

  for (int i=0; i<num_mom; i++)
    TEST_CHECK( gkyl_compare(int_moms[i], int_moms_tar[i], 1e-10) );

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar_ho);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm_tar);
  gkyl_gk_geometry_release(gk_geom_tar);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void test_1x_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6};
  int cells_tar0[] = {12};
  test_1x(cells_do0, cells_tar0, 1, use_gpu);
  test_1x(cells_do0, cells_tar0, 2, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16};
  int cells_tar1[] = {8};
  test_1x(cells_do1, cells_tar1, 1, use_gpu);
  test_1x(cells_do1, cells_tar1, 2, use_gpu);
}

void test_2x_hodev(bool use_gpu)
{
  // Refine along x.
//  int cells_do0[] = {6, 8};
//  int cells_tar0[] = {12, 8};
//  test_2x(cells_do0, cells_tar0, 1, use_gpu);
//  test_2x(cells_do0, cells_tar0, 2, use_gpu);
//
//  // Coarsen along x.
//  int cells_do1[] = {16, 8};
//  int cells_tar1[] = {8, 8};
//  test_2x(cells_do1, cells_tar1, 1, use_gpu);
//  test_2x(cells_do1, cells_tar1, 2, use_gpu);

  // Refine along vpar.
  int cells_do2[] = {96, 96};
  int cells_tar2[] = {128, 128};
  test_2x(cells_do2, cells_tar2, 1, use_gpu);
  test_2x(cells_do2, cells_tar2, 2, use_gpu);

//  // Coarsen along vpar.
//  int cells_do3[] = {8, 12};
//  int cells_tar3[] = {8, 6};
//  test_2x(cells_do3, cells_tar3, 1, use_gpu);
//  test_2x(cells_do3, cells_tar3, 2, use_gpu);
//
//  // Refine along x and vpar.
//  int cells_do4[] = {8, 8};
//  int cells_tar4[] = {32, 16};
//  test_2x(cells_do4, cells_tar4, 1, use_gpu);
//  test_2x(cells_do4, cells_tar4, 2, use_gpu);
//
//  // Coarsen along x and vpar.
//  int cells_do5[] = {8, 12};
//  int cells_tar5[] = {4, 6};
//  test_2x(cells_do5, cells_tar5, 1, use_gpu);
//  test_2x(cells_do5, cells_tar5, 2, use_gpu);
}

void test_1x1v_vlasov_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6, 8};
  int cells_tar0[] = {12, 8};
  test_1x1v_vlasov(cells_do0, cells_tar0, 1, use_gpu);
  test_1x1v_vlasov(cells_do0, cells_tar0, 2, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16, 8};
  int cells_tar1[] = {8, 8};
  test_1x1v_vlasov(cells_do1, cells_tar1, 1, use_gpu);
  test_1x1v_vlasov(cells_do1, cells_tar1, 2, use_gpu);

  // Refine along vx.
  int cells_do2[] = {8, 8};
  int cells_tar2[] = {8, 16};
  test_1x1v_vlasov(cells_do2, cells_tar2, 1, use_gpu);
  test_1x1v_vlasov(cells_do2, cells_tar2, 2, use_gpu);

  // Coarsen along vx.
  int cells_do3[] = {8, 12};
  int cells_tar3[] = {8, 6};
  test_1x1v_vlasov(cells_do3, cells_tar3, 1, use_gpu);
  test_1x1v_vlasov(cells_do3, cells_tar3, 2, use_gpu);

  // Refine along x and vx.
  int cells_do4[] = {8, 8};
  int cells_tar4[] = {32, 16};
  test_1x1v_vlasov(cells_do4, cells_tar4, 1, use_gpu);
  test_1x1v_vlasov(cells_do4, cells_tar4, 2, use_gpu);

  // Coarsen along x and vx.
  int cells_do5[] = {8, 12};
  int cells_tar5[] = {4, 6};
  test_1x1v_vlasov(cells_do5, cells_tar5, 1, use_gpu);
  test_1x1v_vlasov(cells_do5, cells_tar5, 2, use_gpu);
}

void test_1x2v_vlasov_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6, 8, 4};
  int cells_tar0[] = {12, 8, 4};
  test_1x2v_vlasov(cells_do0, cells_tar0, 1, use_gpu);
  test_1x2v_vlasov(cells_do0, cells_tar0, 2, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16, 8, 4};
  int cells_tar1[] = {8, 8, 4};
  test_1x2v_vlasov(cells_do1, cells_tar1, 1, use_gpu);
  test_1x2v_vlasov(cells_do1, cells_tar1, 2, use_gpu);

  // Refine along vpar.
  int cells_do2[] = {8, 8, 4};
  int cells_tar2[] = {8, 16, 4};
  test_1x2v_vlasov(cells_do2, cells_tar2, 1, use_gpu);
  test_1x2v_vlasov(cells_do2, cells_tar2, 2, use_gpu);

  // Coarsen along vpar.
  int cells_do3[] = {8, 12, 4};
  int cells_tar3[] = {8, 6, 4};
  test_1x2v_vlasov(cells_do3, cells_tar3, 1, use_gpu);
  test_1x2v_vlasov(cells_do3, cells_tar3, 2, use_gpu);

  // Refine along mu.
  int cells_do4[] = {8, 6, 4};
  int cells_tar4[] = {8, 6, 8};
  test_1x2v_vlasov(cells_do4, cells_tar4, 1, use_gpu);
  test_1x2v_vlasov(cells_do4, cells_tar4, 2, use_gpu);

  // Coarsen along mu.
  int cells_do5[] = {8, 6, 12};
  int cells_tar5[] = {8, 6, 4};
  test_1x2v_vlasov(cells_do5, cells_tar5, 1, use_gpu);
  test_1x2v_vlasov(cells_do5, cells_tar5, 2, use_gpu);

  // Refine along x and vpar.
  int cells_do6[] = {6, 8, 4};
  int cells_tar6[] = {12, 16, 4};
  test_1x2v_vlasov(cells_do6, cells_tar6, 1, use_gpu);
  test_1x2v_vlasov(cells_do6, cells_tar6, 2, use_gpu);

  // Coarsen along x and vpar.
  int cells_do7[] = {16, 8, 4};
  int cells_tar7[] = {8, 4, 4};
  test_1x2v_vlasov(cells_do7, cells_tar7, 1, use_gpu);
  test_1x2v_vlasov(cells_do7, cells_tar7, 2, use_gpu);

  // Refine along x and mu.
  int cells_do8[] = {6, 8, 4};
  int cells_tar8[] = {12, 8, 8};
  test_1x2v_vlasov(cells_do8, cells_tar8, 1, use_gpu);
  test_1x2v_vlasov(cells_do8, cells_tar8, 2, use_gpu);

  // Coarsen along x and mu.
  int cells_do9[] = {16, 4, 12};
  int cells_tar9[] = {8, 4, 4};
  test_1x2v_vlasov(cells_do9, cells_tar9, 1, use_gpu);
  test_1x2v_vlasov(cells_do9, cells_tar9, 2, use_gpu);

  // Refine along vpar and mu.
  int cells_do10[] = {8, 6, 4};
  int cells_tar10[] = {8, 12, 8};
  test_1x2v_vlasov(cells_do10, cells_tar10, 1, use_gpu);
  test_1x2v_vlasov(cells_do10, cells_tar10, 2, use_gpu);

  // Coarsen along vpar and mu.
  int cells_do11[] = {8, 16, 12};
  int cells_tar11[] = {8, 4, 4};
  test_1x2v_vlasov(cells_do11, cells_tar11, 1, use_gpu);
  test_1x2v_vlasov(cells_do11, cells_tar11, 2, use_gpu);
}

void test_1x1v_gk_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6, 8};
  int cells_tar0[] = {12, 8};
  test_1x1v_gk(cells_do0, cells_tar0, 1, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16, 8};
  int cells_tar1[] = {8, 8};
  test_1x1v_gk(cells_do1, cells_tar1, 1, use_gpu);

  // Refine along vpar.
  int cells_do2[] = {8, 8};
  int cells_tar2[] = {8, 16};
  test_1x1v_gk(cells_do2, cells_tar2, 1, use_gpu);

  // Coarsen along vpar.
  int cells_do3[] = {8, 12};
  int cells_tar3[] = {8, 6};
  test_1x1v_gk(cells_do3, cells_tar3, 1, use_gpu);

  // Refine along x and vpar.
  int cells_do4[] = {8, 8};
  int cells_tar4[] = {32, 16};
  test_1x1v_gk(cells_do4, cells_tar4, 1, use_gpu);

  // Coarsen along x and vpar.
  int cells_do5[] = {8, 12};
  int cells_tar5[] = {4, 6};
  test_1x1v_gk(cells_do5, cells_tar5, 1, use_gpu);
}

void test_1x2v_gk_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6, 8, 4};
  int cells_tar0[] = {12, 8, 4};
  test_1x2v_gk(cells_do0, cells_tar0, 1, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16, 8, 4};
  int cells_tar1[] = {8, 8, 4};
  test_1x2v_gk(cells_do1, cells_tar1, 1, use_gpu);

  // Refine along vpar.
  int cells_do2[] = {8, 8, 4};
  int cells_tar2[] = {8, 16, 4};
  test_1x2v_gk(cells_do2, cells_tar2, 1, use_gpu);

  // Coarsen along vpar.
  int cells_do3[] = {8, 12, 4};
  int cells_tar3[] = {8, 6, 4};
  test_1x2v_gk(cells_do3, cells_tar3, 1, use_gpu);

  // Refine along mu.
  int cells_do4[] = {8, 6, 4};
  int cells_tar4[] = {8, 6, 8};
  test_1x2v_gk(cells_do4, cells_tar4, 1, use_gpu);

  // Coarsen along mu.
  int cells_do5[] = {8, 6, 12};
  int cells_tar5[] = {8, 6, 4};
  test_1x2v_gk(cells_do5, cells_tar5, 1, use_gpu);

  // Refine along x and vpar.
  int cells_do6[] = {6, 8, 4};
  int cells_tar6[] = {12, 16, 4};
  test_1x2v_gk(cells_do6, cells_tar6, 1, use_gpu);

  // Coarsen along x and vpar.
  int cells_do7[] = {16, 8, 4};
  int cells_tar7[] = {8, 4, 4};
  test_1x2v_gk(cells_do7, cells_tar7, 1, use_gpu);

  // Refine along x and mu.
  int cells_do8[] = {6, 8, 4};
  int cells_tar8[] = {12, 8, 8};
  test_1x2v_gk(cells_do8, cells_tar8, 1, use_gpu);

  // Coarsen along x and mu.
  int cells_do9[] = {16, 4, 12};
  int cells_tar9[] = {8, 4, 4};
  test_1x2v_gk(cells_do9, cells_tar9, 1, use_gpu);

  // Refine along vpar and mu.
  int cells_do10[] = {8, 6, 4};
  int cells_tar10[] = {8, 12, 8};
  test_1x2v_gk(cells_do10, cells_tar10, 1, use_gpu);

  // Coarsen along vpar and mu.
  int cells_do11[] = {8, 16, 12};
  int cells_tar11[] = {8, 4, 4};
  test_1x2v_gk(cells_do11, cells_tar11, 1, use_gpu);
}

void test_2x2v_gk_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6, 6, 8, 4};
  int cells_tar0[] = {12, 6, 8, 4};
  test_2x2v_gk(cells_do0, cells_tar0, 1, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16, 6, 8, 4};
  int cells_tar1[] = {8, 6, 8, 4};
  test_2x2v_gk(cells_do1, cells_tar1, 1, use_gpu);

  // Refine along y.
  int cells_do2[] = {6, 6, 8, 4};
  int cells_tar2[] = {6, 12, 8, 4};
  test_2x2v_gk(cells_do2, cells_tar2, 1, use_gpu);

  // Coarsen along y.
  int cells_do3[] = {6, 16, 8, 4};
  int cells_tar3[] = {6, 8, 8, 4};
  test_2x2v_gk(cells_do3, cells_tar3, 1, use_gpu);

  // Refine along x and y.
  int cells_do4[] = {96, 96, 8, 4};
  int cells_tar4[] = {128, 128, 8, 4};
  test_2x2v_gk(cells_do4, cells_tar4, 1, use_gpu);

  // Coarsen along x and y.
  int cells_do5[] = {16, 16, 8, 4};
  int cells_tar5[] = {8, 8, 8, 4};
  test_2x2v_gk(cells_do5, cells_tar5, 1, use_gpu);
}

void test_3x2v_gk_hodev(bool use_gpu)
{
  // Refine along x.
  int cells_do0[] = {6, 6, 8, 8, 4};
  int cells_tar0[] = {12, 6, 8, 8, 4};
  test_3x2v_gk(cells_do0, cells_tar0, 1, use_gpu);

  // Coarsen along x.
  int cells_do1[] = {16, 6, 8, 8, 4};
  int cells_tar1[] = {8, 6, 8, 8, 4};
  test_3x2v_gk(cells_do1, cells_tar1, 1, use_gpu);

  // Refine along y.
  int cells_do2[] = {6, 6, 8, 8, 4};
  int cells_tar2[] = {6, 12, 8, 8, 4};
  test_3x2v_gk(cells_do2, cells_tar2, 1, use_gpu);

  // Coarsen along y.
  int cells_do3[] = {6, 16, 8, 8, 4};
  int cells_tar3[] = {6, 8, 8, 8, 4};
  test_3x2v_gk(cells_do3, cells_tar3, 1, use_gpu);

  // Refine along x and y.
  int cells_do4[] = {96, 96, 8, 8, 4};
  int cells_tar4[] = {128, 128, 8, 8, 4};
  test_3x2v_gk(cells_do4, cells_tar4, 1, use_gpu);

  // Coarsen along x and y.
  int cells_do5[] = {16, 16, 8, 8, 4};
  int cells_tar5[] = {8, 8, 8, 8, 4};
  test_3x2v_gk(cells_do5, cells_tar5, 1, use_gpu);
}

void test_1x_ho()
{
  test_1x_hodev(false);
}

void test_2x_ho()
{
  test_2x_hodev(false);
}

void test_1x1v_vlasov_ho()
{
  test_1x1v_vlasov_hodev(false);
}

void test_1x2v_vlasov_ho()
{
  test_1x2v_vlasov_hodev(false);
}

void test_1x1v_gk_ho()
{
  test_1x1v_gk_hodev(false);
}

void test_1x2v_gk_ho()
{
  test_1x2v_gk_hodev(false);
}

void test_2x2v_gk_ho()
{
  test_2x2v_gk_hodev(false);
}

void test_3x2v_gk_ho()
{
  test_3x2v_gk_hodev(false);
}

#ifdef GKYL_HAVE_CUDA
void test_1x_dev()
{
  test_1x_hodev(true);
}

void test_2x_dev()
{
  test_2x_hodev(true);
}

void test_1x1v_vlasov_dev()
{
  test_1x1v_vlasov_hodev(true);
}

void test_1x2v_vlasov_dev()
{
  test_1x2v_vlasov_hodev(true);
}

void test_1x1v_gk_dev()
{
  test_1x1v_gk_hodev(true);
}

void test_1x2v_gk_dev()
{
  test_1x2v_gk_hodev(true);
}

void test_2x2v_gk_dev()
{
  test_2x2v_gk_hodev(true);
}

void test_3x2v_gk_dev()
{
  test_3x2v_gk_hodev(true);
}
#endif

TEST_LIST = {
  { "test_1x_ho", test_1x_ho },
  { "test_2x_ho", test_2x_ho },
  { "test_1x1v_vlasov_ho", test_1x1v_vlasov_ho },
  { "test_1x2v_vlasov_ho", test_1x2v_vlasov_ho },
  { "test_1x1v_gk_ho", test_1x1v_gk_ho },
  { "test_1x2v_gk_ho", test_1x2v_gk_ho },
  { "test_2x2v_gk_ho", test_2x2v_gk_ho },
  { "test_3x2v_gk_ho", test_3x2v_gk_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_dev", test_1x_dev },
  { "test_2x_dev", test_2x_dev },
  { "test_1x1v_vlasov_dev", test_1x1v_vlasov_dev },
  { "test_1x2v_vlasov_dev", test_1x2v_vlasov_dev },
  { "test_1x1v_gk_dev", test_1x1v_gk_dev },
  { "test_1x2v_gk_dev", test_1x2v_gk_dev },
  { "test_2x2v_gk_dev", test_2x2v_gk_dev },
  { "test_3x2v_gk_dev", test_3x2v_gk_dev },
#endif
  { NULL, NULL },
};
