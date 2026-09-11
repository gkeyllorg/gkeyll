// A test for the integrated moments
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>

#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_radiation_read.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <math.h>


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

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 1.0;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 3.0e19;
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_vthsq(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // D ion mass
  fout[0] =150*eV/mi;
}

void
test_2x_option(bool use_gpu)
{
  int poly_order=1;
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // D ion mass
  double qi = eV; // ion charge
  double Ti = 150*eV;
  double B0 = 1.0; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19; // Particle density in 1/m^3
  double vtIon = sqrt(Ti/mi);
  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = 12.*mi*vtIon*vtIon/(2.0*B0);

  // Phase space and Configuration space extents and resolution
  double lower[] = {0.0, 0.0, -vpar_max_ion, 0.0};
  double upper[] = {1.0, 1.0, vpar_max_ion, mu_max_ion};
  int cells[] = {10, 16, 16, 20};
  int vdim = 2;
  int ndim = sizeof(cells)/sizeof(cells[0]);
  int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim], vLower[vdim], vUpper[vdim];
  int confCells[cdim], vCells[vdim];
  for (int i = 0; i < cdim; i++){
    confLower[i] = lower[i]; 
    confUpper[i] = upper[i];
    confCells[i] = cells[i];
  }
  for (int i = 0; i < vdim; i++){
    vLower[i] = lower[cdim+i]; 
    vUpper[i] = upper[cdim+i];
    vCells[i] = cells[cdim+i];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis, surf_vpar_basis, surf_mu_basis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
    gkyl_cart_modal_serendip(&surf_mu_basis, ndim-1, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
    // constant mu surface
    gkyl_cart_modal_gkhybrid(&surf_mu_basis, cdim, poly_order);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], confGhost[1], 0 , 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vGhost[] = {0, 0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0},
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func, // magnetic field magnitude
    .bfield_ctx = 0 ,
    .position_map = pmap,
    .grid = confGrid,
    .local = confLocal,
    .local_ext = confLocal_ext,
    .global = confLocal,
    .global_ext = confLocal_ext,
    .basis = confBasis,
  };

  int geo_ghost[3] = {1};
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, geo_ghost, &geometry_input.geo_global_ext, &geometry_input.geo_global);
  memcpy(&geometry_input.geo_local, &geometry_input.geo_global, sizeof(struct gkyl_range));
  memcpy(&geometry_input.geo_local_ext, &geometry_input.geo_global_ext, sizeof(struct gkyl_range));

  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Project n, udrift, and vt^2 based on input functions
  struct gkyl_array *m0 = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift = mkarr(false, vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_density, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_upar, 0);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vthsq, 0);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);
  
  // Projection routine expects the primitive moments as a single array.
  struct gkyl_array *prim_moms = mkarr(false, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms, 1.0, m0, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, vtsq, 2*confBasis.num_basis);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // Create distribution function array
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  // Maxwellian (or bi-Maxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &vLocal, 
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .mass = mi,
    .bimaxwellian = false, 
    .use_gpu = use_gpu,
  };
  struct gkyl_gk_maxwellian_proj_on_basis *proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  // If on GPUs, need to copy primitive moments onto device
  struct gkyl_array *prim_moms_dev, *m0_dev;
  if (use_gpu) {
    prim_moms_dev = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
    m0_dev = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    
    gkyl_array_copy(prim_moms_dev, prim_moms);
    gkyl_array_copy(m0_dev, m0);
    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms_dev, false, f);  
  } 
  else {
    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, f);  
  }

  // Initialize integrated moment calculator
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, mi, qi, gvm, gk_geom, NULL, GKYL_F_MOMENT_M0M1M2PARM2PERP, true, use_gpu);    

  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(mcalc);

  struct gkyl_array *marr = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  struct gkyl_array *marr_host = use_gpu? marr_host = mkarr(false, marr->ncomp, marr->size)
                                        : gkyl_array_acquire(marr);
    

  double *red_integ_diag_global;

  if (use_gpu)
    red_integ_diag_global = gkyl_cu_malloc(4*sizeof(double));
  else
    red_integ_diag_global = gkyl_malloc(4*sizeof(double));

  // Now calculate the integrated moments
  double avals_global[2+vdim];
  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc, &local, &confLocal, f, marr);
  gkyl_array_reduce_range(red_integ_diag_global, marr, GKYL_SUM, &confLocal);

  if (use_gpu)
    gkyl_cu_memcpy(avals_global, red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(avals_global, red_integ_diag_global, sizeof(double[2+vdim]));

  // Check the integrated moments are correct. Values computed by Akash Shukla on 2/26/24
  // Check of intM1 really just checks the drift velocity is close to zero, this will not be perfect.
  TEST_CHECK( gkyl_compare( avals_global[0]/3e19, 1.0, 1e-2));
  TEST_CHECK( gkyl_compare( avals_global[1]/3e19, 0.0, 1e-2));
  TEST_CHECK( gkyl_compare( avals_global[2]/2.14e+29, 1.0, 1e-2));
  TEST_MSG( "Got: %.9e | Expected: %.9e\n",avals_global[2]/2.14e+29, 1.0);
  TEST_CHECK( gkyl_compare( avals_global[3]/4.21e+29, 1.0, 1e-2));
  TEST_MSG( "Got: %.9e | Expected: %.9e\n",avals_global[3]/4.21e+29, 1.0);

  gkyl_array_release(m0);
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_array_release(prim_moms);

  gkyl_array_release(marr);
  gkyl_array_release(marr_host);

  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_gk_maxwellian_proj_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

  gkyl_array_release(f);

  if (use_gpu) {
    gkyl_array_release(m0_dev);
    gkyl_array_release(prim_moms_dev);
    gkyl_cu_free(red_integ_diag_global);
  }
  else {
    gkyl_free(red_integ_diag_global);
  }

  gkyl_dg_updater_moment_gyrokinetic_release(mcalc);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);
}

void test_integrated_moms_2x_ho() { test_2x_option(false); }

#ifdef GKYL_HAVE_CUDA

void test_integrated_moms_2x_dev() { test_2x_option(true); }
#endif

TEST_LIST = {
  { "test_integrated_moms_2x_ho", test_integrated_moms_2x_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_integrated_moms_2x_dev", test_integrated_moms_2x_dev },
#endif
  { NULL, NULL },
};
