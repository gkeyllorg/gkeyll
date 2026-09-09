// Test construction of the gyrokinetic LBO drag DG equation object
// (gkyl_dg_lbo_gyrokinetic_drag_new). Verifies the base num_equations, the
// recorded cdim/pdim, species mass, number of conf basis functions, conf-range
// volume and stored geometry/velocity-map pointers, plus that the surface and
// volume kernels are wired up. Reuses the mapc2p geometry + identity
// velocity-map boilerplate from ctest_dg_gyrokinetic.c.
#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>
#include <gkyl_dg_lbo_gyrokinetic_drag_priv.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_position_map.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_velocity_map.h>

static void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

static void
bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 1.0;
}

struct lbo_env {
  struct gkyl_basis basis, confBasis;
  struct gkyl_range confRange, confRange_ext;
  struct gkyl_range phaseRange, phaseRange_ext;
  struct gkyl_rect_grid confGrid, phaseGrid;
  struct gkyl_position_map *pmap;
  struct gk_geometry *gk_geom;
  struct gkyl_velocity_map *gvm;
};

static void
lbo_env_init(struct lbo_env *e)
{
  int cdim = 3, vdim = 2;
  int pdim = cdim + vdim;

  int cells[] = {8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0};
  double lower[] = {0., 0., 0., -1., 0.};
  double upper[] = {1., 1., 1., 1., 1.};

  gkyl_rect_grid_init(&e->confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&e->confGrid, ghost, &e->confRange_ext, &e->confRange);

  gkyl_rect_grid_init(&e->phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&e->phaseGrid, ghost, &e->phaseRange_ext, &e->phaseRange);

  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }
  struct gkyl_rect_grid velGrid;
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int poly_order = 1;
  gkyl_cart_modal_gkhybrid(&e->basis, cdim, vdim);
  gkyl_cart_modal_serendip(&e->confBasis, cdim, poly_order);

  e->pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0, 0.0},
    .mapc2p = mapc2p,
    .c2p_ctx = 0,
    .bfield_func = bfield_func,
    .bfield_ctx = 0,
    .position_map = e->pmap,
    .grid = e->confGrid,
    .local = e->confRange,
    .local_ext = e->confRange_ext,
    .global = e->confRange,
    .global_ext = e->confRange_ext,
    .basis = e->confBasis,
    .geo_grid = e->confGrid,
    .geo_local = e->confRange,
    .geo_local_ext = e->confRange_ext,
    .geo_global = e->confRange,
    .geo_global_ext = e->confRange_ext,
    .geo_basis = e->confBasis,
  };
  e->gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  struct gkyl_mapc2p_inp c2p_in = { };
  e->gvm = gkyl_velocity_map_new(c2p_in, e->phaseGrid, velGrid,
    e->phaseRange, e->phaseRange_ext, velLocal, velLocal_ext, false);
}

static void
lbo_env_release(struct lbo_env *e)
{
  gkyl_gk_geometry_release(e->gk_geom);
  gkyl_position_map_release(e->pmap);
  gkyl_velocity_map_release(e->gvm);
}

void
test_lbo_drag_ctor()
{
  struct lbo_env e;
  lbo_env_init(&e);

  double mass = 2.5;
  struct gkyl_dg_eqn *eqn = gkyl_dg_lbo_gyrokinetic_drag_new(&e.confBasis, &e.basis,
    &e.confRange, &e.phaseGrid, mass, e.gk_geom, e.gvm, false);

  TEST_CHECK( eqn != NULL );
  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );
  TEST_CHECK( eqn->boundary_surf_term != NULL );

  struct dg_lbo_gyrokinetic_drag *drag =
    container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);

  TEST_CHECK( drag->cdim == 3 );
  TEST_CHECK( drag->pdim == 5 );
  TEST_CHECK( drag->mass == mass );
  TEST_CHECK( drag->num_cbasis == e.confBasis.num_basis );
  TEST_CHECK( drag->conf_range.volume == 512 );
  TEST_CHECK( drag->gk_geom != NULL );
  TEST_CHECK( drag->vel_map != NULL );
  TEST_CHECK( drag->surf[0] != NULL );
  TEST_CHECK( drag->boundary_surf[0] != NULL );

  gkyl_dg_eqn_release(eqn);
  lbo_env_release(&e);
}

void
test_lbo_drag_acquire()
{
  struct lbo_env e;
  lbo_env_init(&e);

  struct gkyl_dg_eqn *eqn = gkyl_dg_lbo_gyrokinetic_drag_new(&e.confBasis, &e.basis,
    &e.confRange, &e.phaseGrid, 1.0, e.gk_geom, e.gvm, false);

  struct gkyl_dg_eqn *eqn2 = gkyl_dg_eqn_acquire(eqn);
  TEST_CHECK( eqn2 == eqn );
  gkyl_dg_eqn_release(eqn);
  TEST_CHECK( eqn2->num_equations == 1 );
  gkyl_dg_eqn_release(eqn2);

  lbo_env_release(&e);
}

TEST_LIST = {
  { "lbo_drag_ctor", test_lbo_drag_ctor },
  { "lbo_drag_acquire", test_lbo_drag_acquire },
  { NULL, NULL },
};
