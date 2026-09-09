// Additional coverage for the position-map object: the null (identity) map,
// the inew() input-struct constructor, reference-count acquire/release, and
// the X-point compression setter.  These paths are not exercised by the
// existing ctest_position_map.c.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_position_map.h>
#include <gkyl_position_map_priv.h>
#include <math.h>

void
test_null_map_is_identity()
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  TEST_CHECK( pmap != NULL );
  TEST_CHECK( pmap->id == GKYL_PMAP_USER_INPUT );
  TEST_CHECK( pmap->to_optimize == false );

  // Each of the three maps must act as the identity, and each derivative as 1.
  for (double z = -1.0; z <= 1.0; z += 0.25) {
    for (int i=0; i<3; i++) {
      double x[1] = {z}, y[1] = {0.0};
      pmap->maps[i](0.0, x, y, pmap->ctxs[i]);
      TEST_CHECK( gkyl_compare(y[0], z, 1e-15) );

      double dy[1] = {0.0};
      pmap->map_derivs[i](0.0, x, dy, pmap->ctxs[i]);
      TEST_CHECK( gkyl_compare(dy[0], 1.0, 1e-15) );
    }
  }

  gkyl_position_map_release(pmap);
}

static void
nonuniform_map_1d(double t, const double *xn, double *fout, void *ctx)
{
  // A simple smooth monotone map z -> z + 0.1*sin(z) (identity-like near 0).
  fout[0] = xn[0] + 0.1*sin(xn[0]);
}

void
test_inew_constructor()
{
  int cells[] = {16};
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int dim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  struct gkyl_position_map_inew_inp inp = {
    .pmap_info = {
      .maps = {nonuniform_map_1d, nonuniform_map_1d, nonuniform_map_1d},
      .ctxs = {NULL, NULL, NULL},
    },
    .grid = grid,
    .local = local, .local_ext = local_ext,
    .global = local, .global_ext = local_ext,
    .basis = basis,
  };

  struct gkyl_position_map *pmap = gkyl_position_map_inew(inp);

  TEST_CHECK( pmap != NULL );
  TEST_CHECK( pmap->grid.ndim == 1 );
  TEST_CHECK( pmap->basis.poly_order == 1 );
  TEST_CHECK( pmap->id == GKYL_PMAP_USER_INPUT );

  // The user map must be wired in: check a couple of values.
  double x[1] = {0.3}, y[1];
  pmap->maps[0](0.0, x, y, pmap->ctxs[0]);
  TEST_CHECK( gkyl_compare(y[0], 0.3 + 0.1*sin(0.3), 1e-14) );

  gkyl_position_map_release(pmap);
}

void
test_acquire_refcount()
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_position_map *pmap2 = gkyl_position_map_acquire(pmap);
  TEST_CHECK( pmap2 == pmap );

  // Drop the first reference; the object must survive for use through pmap2.
  gkyl_position_map_release(pmap);

  double x[1] = {0.5}, y[1];
  pmap2->maps[0](0.0, x, y, pmap2->ctxs[0]);
  TEST_CHECK( gkyl_compare(y[0], 0.5, 1e-15) );

  gkyl_position_map_release(pmap2);
}

void
test_set_compression()
{
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // With both compression factors zero the maps stay at the identity backups,
  // but the geometric parameters must be recorded.
  pmap->xpt_ctx->compression_factor = 0.0;
  pmap->xpt_ctx->radial_compression_factor = 0.0;

  double zcut = 1.25, zcenter = 0.1, w = 0.05, psisep = 0.7;
  gkyl_position_map_set_compression(pmap, zcut, zcenter, w, psisep);

  TEST_CHECK( gkyl_compare(pmap->xpt_ctx->zcut, zcut, 1e-15) );
  TEST_CHECK( gkyl_compare(pmap->xpt_ctx->zcenter, zcenter, 1e-15) );
  TEST_CHECK( gkyl_compare(pmap->xpt_ctx->w, w, 1e-15) );
  TEST_CHECK( gkyl_compare(pmap->xpt_ctx->psisep, psisep, 1e-15) );

  // Maps 0 and 1 fall back to identity backups when compression is disabled.
  double x[1] = {0.4}, y[1];
  pmap->maps[0](0.0, x, y, pmap->ctxs[0]);
  TEST_CHECK( gkyl_compare(y[0], 0.4, 1e-15) );
  pmap->maps[1](0.0, x, y, pmap->ctxs[1]);
  TEST_CHECK( gkyl_compare(y[0], 0.4, 1e-15) );

  gkyl_position_map_release(pmap);
}

TEST_LIST = {
  { "null_map_is_identity", test_null_map_is_identity },
  { "inew_constructor", test_inew_constructor },
  { "acquire_refcount", test_acquire_refcount },
  { "set_compression", test_set_compression },
  { NULL, NULL },
};
