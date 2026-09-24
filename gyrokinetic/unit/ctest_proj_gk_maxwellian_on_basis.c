#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_maxwellian_moments.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_den(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_udrift_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void
mapc2p_3x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}


void
bfield_func_3x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 1.0;
}

void
test_1x2v_gk(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, -6.0, 0.0}, upper[] = {1.0, 6.0, 6.0};
  int cells[] = {2, 16, 16};
  int vdim = 2;
  int ndim = sizeof(cells)/sizeof(cells[0]);
  int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *conf_basis_on_dev;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    conf_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) 
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(conf_basis_on_dev, cdim, poly_order);
#endif
  }
  else { 
    basis_on_dev = &basis;
    conf_basis_on_dev = &confBasis;
  }

  int confGhost[] = { 1, 1, 1 }; // 3 elements because it's used by geo.
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create primitive moment arrays
  struct gkyl_array *den, *udrift, *vtsq;
  den = mkarr(confBasis.num_basis, confLocal_ext.volume);
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_den = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_den, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_2v_gk, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_den, 0.0, &confLocal, den);
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // Projection routine expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., den, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0, 0.0},
    .mapc2p = mapc2p_3x, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x, // magnetic field magnitude
    .bfield_ctx = 0 ,
    .position_map = pmap,
    .grid = confGrid,
    .local = confLocal,
    .local_ext = confLocal_ext,
    .global = confLocal,
    .global_ext = confLocal_ext,
    .basis = confBasis,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, confGhost, &geometry_input.geo_local_ext, &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry* gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry if necessary
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Maxwellian (or bi-Maxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal, 
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .mass = mass,
    .bimaxwellian = false, 
    .use_gpu = use_gpu,
  };
  struct gkyl_gk_maxwellian_proj_on_basis *proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  if (use_gpu) {
    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, distf_cu);  
    gkyl_array_copy(distf, distf_cu);
  } 
  else {
    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, distf);  
  }

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p1_vals[] = {  
     7.2307139183122714e-03, 0.0000000000000000e+00, 1.9198293226362615e-04, -7.7970439910196674e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.0701958137127286e-05, 0.0000000000000000e+00, -1.4953406100022537e-04, 0.0000000000000000e+00, 1.6124599381836546e-05, 0.0000000000000000e+00,
    -8.2719200283232917e-19, 0.0000000000000000e+00, -3.4806248503322844e-20, 0.0000000000000000e+00, };
  double p2_vals[] = { 
    7.2307468609012666e-03, 0.0000000000000000e+00, 1.9198380692343289e-04, -7.8092230706225602e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.0734294852987710e-05, 3.6591823321385775e-18, -1.4953474226616330e-04, 3.7739922227981074e-05, 0.0000000000000000e+00, 7.0473141211557788e-19,
    0.0000000000000000e+00, -4.8789097761847700e-19, 1.6149786206441256e-05, 0.0000000000000000e+00, 1.0020339643610290e-06, 5.4210108624275222e-20,
    0.0000000000000000e+00, 0.0000000000000000e+00 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-2) );
    }
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-2) );
  }

  // write distribution function to file
  char fname[1024];
  if (use_gpu) {
    sprintf(fname, "ctest_proj_gkmaxwellian_on_basis_prim_mom_1x2v_p%d_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname, "ctest_proj_gkmaxwellian_on_basis_prim_mom_1x2v_p%d_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  gkyl_array_release(den); 
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_den);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(prim_moms_ho);
  if (use_gpu) {
    gkyl_array_release(prim_moms);
  }
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(distf);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
  }
  gkyl_gk_maxwellian_proj_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);
  gkyl_position_map_release(pmap);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(conf_basis_on_dev);
  }
#endif  
}

void eval_den_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = 1.0*(1.0+0.2*cos(x))*(1+0.2*cos(y))*(1.0+0.2*cos(z));
}

void
test_3x2v_gk(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, 0.1, 0.1, -6.0, 0.0}, upper[] = {1.0, 1.0, 1.0, 6.0, 6.0};
  int cells[] = {2, 2, 2, 16, 16};
  int vdim = 2;
  int ndim = sizeof(cells)/sizeof(cells[0]);
  int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *conf_basis_on_dev;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    conf_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) 
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(conf_basis_on_dev, cdim, poly_order);
#endif
  }
  else { 
    basis_on_dev = &basis;
    conf_basis_on_dev = &confBasis;
  }

  int confGhost[] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], confGhost[1], confGhost[2], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&velGrid, vGhost, &velLocal_ext, &velLocal);

  // Create primitive moment arrays
  struct gkyl_array *den, *udrift, *vtsq;
  den = mkarr(confBasis.num_basis, confLocal_ext.volume);
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_den = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_den_3x, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_2v_gk, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_den, 0.0, &confLocal, den);
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // Projection routine expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., den, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .mapc2p = mapc2p_3x, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x, // magnetic field magnitude
    .bfield_ctx =0 ,
    .position_map = pmap,
    .grid = confGrid,
    .local = confLocal,
    .local_ext = confLocal_ext,
    .global = confLocal,
    .global_ext = confLocal_ext,
    .basis = confBasis,
    .geo_grid = confGrid,
    .geo_local = confLocal,
    .geo_local_ext = confLocal_ext,
    .geo_global = confLocal,
    .geo_global_ext = confLocal_ext,
    .geo_basis = confBasis,
  };
  struct gk_geometry* gk_geom;
  gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Maxwellian (or bi-Maxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal, 
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .mass = mass,
    .bimaxwellian = false, 
    .use_gpu = use_gpu,
  };
  struct gkyl_gk_maxwellian_proj_on_basis *proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  if (use_gpu) {
    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, distf_cu);  
    gkyl_array_copy(distf, distf_cu);
  } 
  else {
    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, distf);  
  }

  // Values to compare at index (1, 1, 1, 9, 9) [remember, lower-left index is (1,1,1)]
  // They come from the gpu run results using the previous parallelization. 
  double p1_vals[] = {  
    2.4243050727223950e-02, -1.6882464773457025e-04, -1.6882464773457199e-04,
   -1.6882464773457155e-04,  6.4367806805923898e-04, -2.6141835388891624e-03,
    1.1756672872340893e-06,  1.1756672872319208e-06,  1.1756672872319208e-06,
   -4.4824689894557621e-06, -4.4824689894555452e-06, -4.4824689894551116e-06,
    1.8204747415343751e-05,  1.8204747415343968e-05,  1.8204747415343643e-05,
   -6.9409276447821217e-05, -8.1871550658971042e-09,  3.1215182306927454e-08,
    3.1215182306986178e-08,  3.1215182307144295e-08, -1.2677488918737411e-07,
   -1.2677488918650675e-07, -1.2677488918661517e-07,  4.8335487054231844e-07,
    4.8335487054236820e-07,  4.8335487054221002e-07, -2.1737743278163100e-10,
    8.8283963297070298e-10, -3.3660044136686707e-09, -3.3660044141023515e-09,
   -3.3660044141023515e-09,  2.3440305414670191e-11, -5.0135600263416188e-04,
    3.4913613590419644e-06,  3.4913613590417975e-06,  3.4913613590416391e-06,
    5.4062363023384633e-05, -2.4313270560376134e-08, -2.4313270561026655e-08,
   -2.4313270560918235e-08, -3.7648147074447582e-07, -3.7648147074447582e-07,
   -3.7648147074469266e-07,  1.6931364740599317e-10,  2.6217555038242103e-09,
    2.6217555040410507e-09,  2.6217555040410507e-09, -1.8257477671200626e-11
  };
  double p2_vals[] = { 
    2.4243271148107627e-02, -1.6844589780429920e-04, -1.6844589780430397e-04,
   -1.6844589780430353e-04,  6.4368392046165463e-04, -2.6182788029943962e-03,
    1.1703874577738927e-06,  1.1703874577719411e-06,  1.1703874577718327e-06,
   -4.4724127871178736e-06, -4.4724127871176025e-06, -4.4724127871178194e-06,
    1.8192195309699262e-05,  1.8192195309699208e-05,  1.8192195309699167e-05,
   -6.9518009944984270e-05, -2.9116288832330320e-05, -2.9116288832326417e-05,
   -2.9116288832326417e-05, -5.0136056102635131e-04,  1.2653453167179150e-04,
   -8.1320282624845902e-09,  3.1074997374563754e-08,  3.1074997374780594e-08,
    3.1074997374726384e-08, -1.2640211187902229e-07, -1.2640211187891387e-07,
   -1.2640211187896808e-07,  4.8302159915686031e-07,  4.8302159915691452e-07,
    4.8302159915684676e-07,  2.0230435831401641e-07,  2.0230435831618481e-07,
    2.0230435831542587e-07,  2.0230435831553429e-07,  2.0230435831596797e-07,
    2.0230435831607639e-07, -7.7306757946756566e-07, -7.7306757946751145e-07,
   -7.7306757946767408e-07,  3.4835286587286006e-06,  3.4835286587286006e-06,
    3.4835286587284921e-06,  3.1445658222361780e-06,  3.1445658222366930e-06,
    3.1445658222365439e-06,  5.4147054725953100e-05, -8.7918097605268546e-07,
   -8.7918097605265835e-07, -8.7918097605272612e-07,  3.3596226731406998e-06,
   -2.1591376005545164e-10,  8.7826090376862898e-10, -3.3561067908167493e-09,
   -3.3561067907896443e-09, -3.3561067907625392e-09, -1.4056411382700062e-09,
   -1.4056411383784264e-09, -1.4056411384868466e-09,  5.3713899286200678e-09,
    5.3713899286742779e-09,  5.3713899289995386e-09,  5.3713899289453285e-09,
    5.3713899286742779e-09,  5.3713899286742779e-09, -2.4204081572204428e-08,
   -2.4204081572123112e-08, -2.4204081572177323e-08, -2.1848916752811268e-08,
   -2.1848916752838373e-08, -2.1848916752960345e-08, -2.1848916752960345e-08,
   -2.1848916752770610e-08, -2.1848916752797715e-08,  8.3491474571778701e-08,
    8.3491474571765148e-08,  8.3491474571860016e-08, -3.7622188817055589e-07,
   -3.7622188817055589e-07, -3.7622188817074563e-07,  6.1086817838469061e-09,
    6.1086817838062485e-09,  6.1086817838469061e-09, -2.3343164130186359e-08,
   -2.3343164130145701e-08, -2.3343164130172806e-08,  2.3318735367687375e-11,
   -3.7321226043770821e-11, -3.7321226097980930e-11, -3.7321225881140495e-11,
    1.6817360273756146e-10,  1.5180956285367659e-10,  1.5180956285367659e-10,
    1.5180956288078164e-10, -5.8011133512704236e-10, -5.8011133512704236e-10,
   -5.8011133522191005e-10, -5.8011133516769994e-10, -5.8011133516769994e-10,
   -5.8011133519480499e-10,  2.6140463198124999e-09,  2.6140463198396050e-09,
    2.6140463198124999e-09, -4.2444040868023796e-11,  1.6219181891657465e-10,
    1.6219181869973422e-10,  1.6219181865907664e-10,  4.0307008461021814e-12,
    4.0307008732072357e-12,  4.0307008189971271e-12, -1.8162787487545735e-11,
   -1.1269332393403825e-12
  };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[5]) { 1, 1, 1, 9, 9 }));
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-2) );
    }
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-2) );
  }

  // Write distribution function to file
  char fname[1024];
  if (use_gpu) {
    sprintf(fname, "ctest_proj_gkmaxwellian_on_basis_prim_mom_3x2v_p%d_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname, "ctest_proj_gkmaxwellian_on_basis_prim_mom_3x2v_p%d_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Calculate the moments and copy from device to host.
  struct gkyl_gk_maxwellian_moments_inp inp_calc = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true,
    .mass = mass,
    .use_gpu = use_gpu,
  };
  gkyl_gk_maxwellian_moments *calc_moms = gkyl_gk_maxwellian_moments_inew( &inp_calc );

  struct gkyl_array *moms;
  moms = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_cu;
  if (use_gpu) {
    moms_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_gk_maxwellian_moments_advance(calc_moms, &local, &confLocal, distf_cu, moms_cu);
    gkyl_array_copy(moms, moms_cu);
  }
  else {
    gkyl_gk_maxwellian_moments_advance(calc_moms, &local, &confLocal, distf, moms);
  }
  gkyl_gk_maxwellian_moments_release(calc_moms);

  // Write moments to file
  char fname_moms[1024];
  if (use_gpu) {
    sprintf(fname_moms, "ctest_proj_gkmaxwellian_on_basis_prim_mom_3x2v_p%d_moms_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname_moms, "ctest_proj_gkmaxwellian_on_basis_prim_mom_3x2v_p%d_moms_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms, fname_moms);

  gkyl_array_release(den); 
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_den);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(prim_moms_ho);
  if (use_gpu) {
    gkyl_array_release(prim_moms);
  }
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(distf);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
  }
  gkyl_gk_maxwellian_proj_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);
  gkyl_position_map_release(pmap);

  gkyl_array_release(moms);
  if (use_gpu) 
    gkyl_array_release(moms_cu);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(conf_basis_on_dev);
  }
#endif  
}

void test_1x2v_p1_gk() { test_1x2v_gk(1, false); }
void test_3x2v_p1_gk() { test_3x2v_gk(1, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x2v_p1_gk_gpu() { test_1x2v_gk(1, true); }
void test_3x2v_p1_gk_gpu() { test_3x2v_gk(1, true); }
#endif

TEST_LIST = {
  { "test_1x2v_p1_gk", test_1x2v_p1_gk },
  { "test_3x2v_p1_gk", test_3x2v_p1_gk },

#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_p1_gk_gpu", test_1x2v_p1_gk_gpu },
  { "test_3x2v_p1_gk_gpu", test_3x2v_p1_gk_gpu },
#endif
  { NULL, NULL },
};
