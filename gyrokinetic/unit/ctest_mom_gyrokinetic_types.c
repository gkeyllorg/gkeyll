// Test construction of gyrokinetic moment-type objects for every supported
// moment name, verifying the reported num_mom and basic dimension fields.
// Covers both the standard (gkyl_mom_gyrokinetic_new) and integrated
// (gkyl_int_mom_gyrokinetic_new) constructors. Host only.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_position_map.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <math.h>

static void
mapc2p_3x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

static void
bfield_func_3x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 1.0; // uniform Bz
}

// Shared fixture holding everything the moment constructor needs.
struct fixture {
  struct gkyl_basis basis, confBasis;
  struct gkyl_range confLocal;
  struct gk_geometry *gk_geom;
  struct gkyl_velocity_map *gvm;
  struct gkyl_position_map *pmap;
};

static void
make_fixture(struct fixture *fx, int poly_order, int vdim)
{
  const int cdim = 1;
  const int ndim = cdim + vdim;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  lower[0] = -M_PI; upper[0] = M_PI; cells[0] = 4;
  for (int d=0; d<vdim; d++) {
    lower[cdim+d] = (d==0)? -2.0 : 0.0;
    upper[cdim+d] = (d==0)?  2.0 : 2.0;
    cells[cdim+d] = 2;
  }

  double confLower[1] = {lower[0]}, confUpper[1] = {upper[0]};
  int confCells[1] = {cells[0]};
  double velLower[3], velUpper[3];
  int velCells[3];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  struct gkyl_rect_grid grid, confGrid, velGrid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  gkyl_cart_modal_serendip(&fx->confBasis, cdim, poly_order);
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&fx->basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&fx->basis, ndim, poly_order);

  int confGhost[] = {1, 1, 1};
  struct gkyl_range confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &fx->confLocal);

  int velGhost[3] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);
  int ghost[GKYL_MAX_DIM] = {0};
  ghost[0] = confGhost[0];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  fx->pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0, 0.0},
    .mapc2p = mapc2p_3x,
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x,
    .bfield_ctx = 0,
    .position_map = fx->pmap,
    .grid = confGrid,
    .local = fx->confLocal,
    .local_ext = confLocal_ext,
    .global = fx->confLocal,
    .global_ext = confLocal_ext,
    .basis = fx->confBasis,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, confGhost,
    &geometry_input.geo_local_ext, &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry *gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  fx->gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  struct gkyl_mapc2p_inp c2p_in = { };
  fx->gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid, local, local_ext,
    velLocal, velLocal_ext, false);
}

static void
free_fixture(struct fixture *fx)
{
  gkyl_gk_geometry_release(fx->gk_geom);
  gkyl_velocity_map_release(fx->gvm);
  gkyl_position_map_release(fx->pmap);
}

static struct gkyl_mom_type*
mk_mom(struct fixture *fx, enum gkyl_distribution_moments mt)
{
  return gkyl_mom_gyrokinetic_new(&fx->confBasis, &fx->basis, &fx->confLocal,
    1.0, 1.0, fx->gvm, fx->gk_geom, NULL, mt, false);
}

void
test_mom_types_1x2v()
{
  const int vdim = 2;
  struct fixture fx;
  make_fixture(&fx, 1, vdim);

  struct { enum gkyl_distribution_moments mt; int nm; } cases[] = {
    { GKYL_F_MOMENT_M0,             1 },
    { GKYL_F_MOMENT_M1,             1 },
    { GKYL_F_MOMENT_M2,             1 },
    { GKYL_F_MOMENT_M2PAR,          1 },
    { GKYL_F_MOMENT_M2PERP,         1 },
    { GKYL_F_MOMENT_M3PAR,          1 },
    { GKYL_F_MOMENT_M3PERP,         1 },
    { GKYL_F_MOMENT_M0M1M2,         3 },
    { GKYL_F_MOMENT_M0M1M2PARM2PERP, vdim+2 },
    { GKYL_F_MOMENT_HAMILTONIAN,    3 },
  };
  int ncase = sizeof(cases)/sizeof(cases[0]);

  for (int i=0; i<ncase; i++) {
    struct gkyl_mom_type *m = mk_mom(&fx, cases[i].mt);
    TEST_CHECK( m != NULL );
    TEST_MSG("moment case index %d", i);
    TEST_CHECK( m->cdim == 1 );
    TEST_CHECK( m->pdim == 3 );
    TEST_CHECK( m->poly_order == 1 );
    TEST_CHECK( m->num_config == fx.confBasis.num_basis );
    TEST_CHECK( m->num_phase == fx.basis.num_basis );
    TEST_CHECK( m->num_mom == cases[i].nm );
    TEST_CHECK( gkyl_mom_type_num_mom(m) == cases[i].nm );
    gkyl_mom_type_release(m);
  }

  free_fixture(&fx);
}

void
test_int_mom_types_1x2v()
{
  const int vdim = 2;
  struct fixture fx;
  make_fixture(&fx, 1, vdim);

  // The integrated-moment constructor combines several moments into a single
  // object; verify the multi-moment names produce the expected counts.
  struct { enum gkyl_distribution_moments mt; int nm; } cases[] = {
    { GKYL_F_MOMENT_M0M1M2,          3 },
    { GKYL_F_MOMENT_M0M1M2PARM2PERP, vdim+2 },
    { GKYL_F_MOMENT_HAMILTONIAN,     3 },
  };
  int ncase = sizeof(cases)/sizeof(cases[0]);

  for (int i=0; i<ncase; i++) {
    struct gkyl_mom_type *m = gkyl_int_mom_gyrokinetic_new(&fx.confBasis, &fx.basis,
      &fx.confLocal, 1.0, 1.0, fx.gvm, fx.gk_geom, NULL, cases[i].mt, false);
    TEST_CHECK( m != NULL );
    TEST_MSG("int moment case index %d", i);
    TEST_CHECK( m->cdim == 1 );
    TEST_CHECK( m->pdim == 3 );
    TEST_CHECK( m->num_mom == cases[i].nm );
    gkyl_mom_type_release(m);
  }

  free_fixture(&fx);
}

void
test_mom_types_1x1v()
{
  const int vdim = 1;
  struct fixture fx;
  make_fixture(&fx, 1, vdim);

  // With a single velocity dimension, the par/perp split collapses; check a
  // representative subset that is valid for 1v.
  struct { enum gkyl_distribution_moments mt; int nm; } cases[] = {
    { GKYL_F_MOMENT_M0,     1 },
    { GKYL_F_MOMENT_M1,     1 },
    { GKYL_F_MOMENT_M2,     1 },
    { GKYL_F_MOMENT_M2PAR,  1 },
    { GKYL_F_MOMENT_M3PAR,  1 },
    { GKYL_F_MOMENT_M0M1M2, 3 },
  };
  int ncase = sizeof(cases)/sizeof(cases[0]);

  for (int i=0; i<ncase; i++) {
    struct gkyl_mom_type *m = mk_mom(&fx, cases[i].mt);
    TEST_CHECK( m != NULL );
    TEST_MSG("1v moment case index %d", i);
    TEST_CHECK( m->cdim == 1 );
    TEST_CHECK( m->pdim == 2 );
    TEST_CHECK( m->num_mom == cases[i].nm );
    gkyl_mom_type_release(m);
  }

  free_fixture(&fx);
}

TEST_LIST = {
  { "mom_types_1x2v", test_mom_types_1x2v },
  { "int_mom_types_1x2v", test_int_mom_types_1x2v },
  { "mom_types_1x1v", test_mom_types_1x1v },
  { NULL, NULL },
};
