#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_gyrokinetic_priv.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_velocity_map.h>
#include <gkyl_basis.h>
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

void
test_dg_gyrokinetic_ho()
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
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
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

  double charge = 1.;
  double mass = 1.;

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, phaseGrid, velGrid,
    phaseRange, phaseRange_ext, velLocal, velLocal_ext, false);

  struct gkyl_dg_eqn* eqn = gkyl_dg_gyrokinetic_new(&confBasis, &basis, &confRange, &phaseRange, 
    charge, mass, 0, gk_geom, gvm, false);

  TEST_CHECK( eqn->num_equations == 1 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  TEST_CHECK( gyrokinetic->cdim == 3 );
  TEST_CHECK( gyrokinetic->pdim == 5 );
  TEST_CHECK( gyrokinetic->conf_range.volume == 512 );

  gkyl_gk_geometry_release(gk_geom);  
  gkyl_position_map_release(pmap);
  gkyl_velocity_map_release(gvm);
  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

#endif

TEST_LIST = {
  { "dg_gyrokinetic_ho", test_dg_gyrokinetic_ho },
#ifdef GKYL_HAVE_CUDA
#endif  
  { NULL, NULL },
};
