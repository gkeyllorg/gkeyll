#include <gkyl_array.h>
#include <gkyl_util.h>
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_array_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
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

void eval_prim_moms_1x2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den = 1.0;
  double upar = 0.5;
  double tpar = 1.3;
  double tperp = 0.6;
  double mass = 1.0;
  fout[0] = den;  // Density.
  fout[1] = upar;  // Parallel drift speed.
  fout[2] = tpar/mass;  // Parallel temperature divided by mass (vtpar^2).
  fout[3] = tperp/mass;  // Perpendicular temperature divided by mass (vtperp^2).
}

void
test_1x2v_gk(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, -6.0, 0.0}, upper[] = {1.0, 6.0, 6.0};
  int cells[] = {2, 16, 16};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

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

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0, 0.0},
    .mapc2p = mapc2p_3x, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x, // magnetic field magnitude
    .bfield_ctx = 0,
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

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // create moment arrays
  struct gkyl_array *prim_moms, *prim_moms_ho;
  prim_moms_ho = mkarr(4*confBasis.num_basis, confLocal_ext.volume);
  if (use_gpu) { // create device copies
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 4*confBasis.num_basis, confLocal_ext.volume);
  } else {
    prim_moms = prim_moms_ho;
  }

  gkyl_proj_on_basis *proj_prim_moms = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 4, eval_prim_moms_1x2v_gk, NULL);

  gkyl_proj_on_basis_advance(proj_prim_moms, 0.0, &confLocal, prim_moms_ho);
  gkyl_array_copy(prim_moms, prim_moms_ho);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // bi-Maxwellian projection updater.
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
    .bimaxwellian = true, 
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
    1.2845117649060e-03, 3.4098924955929e-20, 2.6352311336353e-05, -2.2927175349421e-04, -1.0358294364538e-20, 
    1.1737441012087e-20, -4.7036086346421e-06, 3.2926920104084e-21, -2.0499212907186e-05, -6.6910161820819e-21, 
    3.6588925200118e-06, -6.7667527142105e-21
  };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
    }
  }

  // release memory for moment data object
  gkyl_velocity_map_release(gvm);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(prim_moms);
  gkyl_proj_on_basis_release(proj_prim_moms);
  gkyl_array_release(distf);
  if (use_gpu) {
    gkyl_array_release(prim_moms_ho);
    gkyl_array_release(distf_cu);
  }
  gkyl_gk_maxwellian_proj_on_basis_release(proj_max);
}

void test_proj_bimaxwellian_1x2v_p1_gk_ho() { test_1x2v_gk(1, false); }

#ifdef GKYL_HAVE_CUDA
void test_proj_bimaxwellian_1x2v_p1_gk_dev() { test_1x2v_gk(1, true); }
#endif

TEST_LIST = {
  { "test_proj_bimaxwellian_1x2v_p1_gk_ho", test_proj_bimaxwellian_1x2v_p1_gk_ho },

#ifdef GKYL_HAVE_CUDA
  { "test_proj_bimaxwellian_1x2v_p1_gk_dev", test_proj_bimaxwellian_1x2v_p1_gk_dev },
#endif
  { NULL, NULL },
};
