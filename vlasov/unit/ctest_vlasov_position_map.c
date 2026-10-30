#include <math.h>

#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Monotone maps of [-1, 1], different in each direction.
static void
pmap_cubic(double t, const double *xn, double *fout, void *ctx)
{
  double s = xn[0];
  fout[0] = 0.25 * s + 0.75 * s * s * s;
}

static void
pmap_quad(double t, const double *xn, double *fout, void *ctx)
{
  double s = xn[0];
  fout[0] = s + 0.3 * (s * s - 1.0);
}

static void
pmap_exact(int d, double s, double *xp)
{
  if (d == 0) {
    pmap_cubic(0.0, &s, xp, 0);
  } else {
    pmap_quad(0.0, &s, xp, 0);
  }
}

// Evaluate the computational-to-physical map on every rank of a decomposed
// grid at points inside its local cells: each rank must reproduce the map of
// the undecomposed grid (the cell lookup has to use the global cell index, not
// the rank's local range offset), and at cell vertices the C^0 linear map must
// interpolate the user's function.
static void
test_pmap_decomposed_2x(void)
{
  int cdim = 2;
  int cells[] = {8, 9};
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM] = {
    {.eval_pmap = pmap_cubic}, {.eval_pmap = pmap_quad}
  };

  int nghost[] = {1, 1};
  struct gkyl_range global, global_ext;
  gkyl_create_grid_ranges(&grid, nghost, &global_ext, &global);
  struct gkyl_vlasov_position_map *pm_global =
    gkyl_vlasov_position_map_new(&grid, &global, &global_ext, &basis, inp, false);

  int cuts[] = {2, 3};
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &global);
  const double xlog[] = {-0.9, -0.3, 0.4, 0.95};
  for (int r = 0; r < decomp->ndecomp; ++r) {
    struct gkyl_range local, local_ext;
    gkyl_create_ranges(&decomp->ranges[r], nghost, &local_ext, &local);
    struct gkyl_vlasov_position_map *pm_local =
      gkyl_vlasov_position_map_new(&grid, &local, &local_ext, &basis, inp, false);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      double xcc[GKYL_MAX_CDIM];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xcc);

      // Interior points: identical to the undecomposed map.
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          double xc[] = {xcc[0] + 0.5 * grid.dx[0] * xlog[i], xcc[1] + 0.5 * grid.dx[1] * xlog[j]};
          double xp_local[GKYL_MAX_CDIM], xp_global[GKYL_MAX_CDIM];
          gkyl_vlasov_position_map_eval_mc2p(pm_local, xc, xp_local);
          gkyl_vlasov_position_map_eval_mc2p(pm_global, xc, xp_global);
          for (int d = 0; d < cdim; ++d) {
            TEST_CHECK(xp_local[d] == xp_global[d]);
            TEST_MSG(
              "rank %d cell (%d,%d) dir %d: %.16e vs undecomposed %.16e", r, iter.idx[0],
              iter.idx[1], d, xp_local[d], xp_global[d]
            );
          }
        }
      }

      // Lower-left vertex: the C^0 map interpolates the user's function.
      double xv[] = {xcc[0] - 0.5 * grid.dx[0], xcc[1] - 0.5 * grid.dx[1]};
      double xp[GKYL_MAX_CDIM];
      gkyl_vlasov_position_map_eval_mc2p(pm_local, xv, xp);
      for (int d = 0; d < cdim; ++d) {
        double xe;
        pmap_exact(d, xv[d], &xe);
        TEST_CHECK(gkyl_compare_double(xp[d], xe, 1e-14));
        TEST_MSG("rank %d vertex dir %d: %.16e vs exact %.16e", r, d, xp[d], xe);
      }
    }
    gkyl_vlasov_position_map_release(pm_local);
  }
  gkyl_rect_decomp_release(decomp);
  gkyl_vlasov_position_map_release(pm_global);
}

TEST_LIST = {{"pmap_decomposed_2x", test_pmap_decomposed_2x}, {NULL, NULL}};
