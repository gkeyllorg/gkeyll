#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>

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
test_3x2v_p1(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 3, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0};
  double lower[] = {0., 0., 0., -1., 0.};
  double upper[] = {1., 1., 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }
  struct gkyl_rect_grid velGrid;
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, surf_basis, confBasis; // phase-space, conf-space basis

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
    gkyl_cart_modal_gkhybrid(&surf_basis, cdim-1, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bfield_func = bfield_func, // magnetic field magnitude
      .bfield_ctx =0 ,
      .position_map = pmap,
      .grid = confGrid,
      .local = confRange,
      .local_ext = confRange_ext,
      .global = confRange,
      .global_ext = confRange_ext,
      .basis = confBasis,
      .geo_grid = confGrid,
      .geo_local = confRange,
      .geo_local_ext = confRange_ext,
      .geo_global = confRange,
      .geo_global_ext = confRange_ext,
      .geo_basis = confBasis,
  };

  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  // Initialize gyrokinetic variables
  struct gkyl_array *flux_surf = mkarr1(use_gpu, 4*surf_basis.num_basis, phaseRange_ext.volume);
  struct gkyl_array *phi = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  struct gkyl_array *apar = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  struct gkyl_array *apardot = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  struct gkyl_dg_gyrokinetic_auxfields aux = { .flux_surf = flux_surf, 
    .phi = phi, .apar = apar, .apardot = apardot };

  const bool is_zero_flux[GKYL_MAX_DIM] = {false};

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, phaseGrid, velGrid,
    phaseRange, phaseRange_ext, velLocal, velLocal_ext, false);

  struct gkyl_dg_updater_gyrokinetic* up;
  up = gkyl_dg_updater_gyrokinetic_new(&phaseGrid, &confBasis, &basis, &confRange, &phaseRange,
    is_zero_flux, 1.0, 1.0, GKYL_GK_COLLISIONLESS_EM, false, false, gk_geom, gvm, &aux, use_gpu);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);

  struct timespec tm = gkyl_wall_clock();
  // run hyper_dg_advance
  int nrep = 1;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_updater_gyrokinetic_advance(up, &phaseRange, fin, cflrate, rhs);
  }
  double gk_tm = gkyl_time_diff_now_sec(tm);

  // printf("\ngyrokinetic update on (%d, %d, %d, %d, %d) took %g sec\n", cells[0], cells[1], cells[2], cells[3], cells[4], gk_tm); 

  // clean up
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_position_map_release(pmap);
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(flux_surf);
  gkyl_array_release(phi);
  gkyl_array_release(apar);
  gkyl_array_release(apardot);

  gkyl_dg_updater_gyrokinetic_release(up);
}

void
test_gyrokinetic_3x2v_p1()
{
  test_3x2v_p1(false);
}

void
test_gyrokinetic_3x2v_p1_cu()
{
  test_3x2v_p1(true);
}

TEST_LIST = {
  { "test_gyrokinetic_3x2v_p1", test_gyrokinetic_3x2v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "test_gyrokinetic_3x2v_p1_cu", test_gyrokinetic_3x2v_p1_cu },
#endif
  { NULL, NULL },
};

