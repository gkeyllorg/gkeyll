#include <acutest.h>
#include <gkyl_gk_collisionless_flux_kernels.h>

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
mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0];
  xp[1] = xc[1];
  xp[2] = xc[2];
}

void
bfield_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
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
  int pdim = cdim + vdim;

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
  for (int d = 0; d < vdim; d++) {
    velLower[d] = lower[cdim + d];
    velUpper[d] = upper[cdim + d];
    velCells[d] = cells[cdim + d];
  }
  struct gkyl_rect_grid velGrid;
  int velGhost[] = {0, 0};
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
    .bfield_ctx = 0,
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
  struct gkyl_mapc2p_inp c2p_in = {};
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(
    c2p_in, phaseGrid, velGrid, phaseRange, phaseRange_ext, velLocal, velLocal_ext, false
  );

  struct gkyl_dg_eqn *eqn = gkyl_dg_gyrokinetic_new(
    &confBasis, &basis, &confRange, &phaseRange, charge, mass, 0, gk_geom, gvm, false
  );

  TEST_CHECK(eqn->num_equations == 1);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  TEST_CHECK(gyrokinetic->cdim == 3);
  TEST_CHECK(gyrokinetic->pdim == 5);
  TEST_CHECK(gyrokinetic->conf_range.volume == 512);

  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);
  gkyl_velocity_map_release(gvm);
  gkyl_dg_eqn_release(eqn);
}

// The reduced 1D operator must ignore curvature while retaining electric and
// mirror acceleration. Exercise both volume and velocity-surface kernels.
static void
test_dg_gyrokinetic_parallel_only_1x(void)
{
  // Cover both v_parallel alone and (v_parallel, mu) velocity spaces.
  for (int vdim = 1; vdim <= 2; ++vdim) {
    double w[] = {0.1, 0.3, 0.5}, dx[] = {0.4, 0.5, 0.6};
    // Modal coefficients of the linear velocity maps and v_parallel squared.
    double vm[] = {0.3 * sqrt(2.0), 0.5 / sqrt(6.0), 0.5 * sqrt(2.0), 0.6 / sqrt(6.0)};
    double vs[] = {
      sqrt(2.0) * (0.09 + 0.25 / 12.0), sqrt(2.0 / 3.0) * 0.3 * 0.5, sqrt(2.0 / 5.0) * 0.25 / 6.0
    };
    double b[] = {2.0 * sqrt(2.0), 0.2}, phi[] = {0.0, 0.3}, b3[] = {sqrt(2.0), 0.0};
    double curl0[6] = {0.0}, curl1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, bi[6] = {0.0};
    double f[12], a[12] = {0.0}, bcur[12] = {0.0}, electric[12] = {0.0}, mirror[12] = {0.0},
                  flat[12] = {0.0};
    for (int k = 0; k < 12; ++k) {
      f[k] = 0.1 * (k + 1);
    }

    double (*vol)(const double *, const double *, const double *, const double *, double, double, const double *, const double *, const double *, const double *, const double *, const double *, double *) =
      vdim == 1 ? dg_gyrokinetic_vol_1x1v_ser_p1 : dg_gyrokinetic_vol_1x2v_ser_p1;

    // Changing curl(bhat) must leave the reduced 1D volume update unchanged.
    vol(w, dx, vm, vs, 1.0, 1.0, b, phi, curl0, b3, bi, f, a);
    vol(w, dx, vm, vs, 1.0, 1.0, b, phi, curl1, b3, bi, f, bcur);

    int nb = vdim == 1 ? 6 : 12;
    for (int k = 0; k < nb; ++k) {
      TEST_CHECK(gkyl_compare(a[k], bcur[k], 1e-14));
    }

    // Isolate electric and mirror acceleration relative to uniform B and zero phi.
    double zero[2] = {0.0}, uniform[] = {b[0], 0.0};
    vol(w, dx, vm, vs, 1.0, 1.0, uniform, zero, curl1, b3, bi, f, flat);
    vol(w, dx, vm, vs, 1.0, 1.0, uniform, phi, curl1, b3, bi, f, electric);
    vol(w, dx, vm, vs, 1.0, 1.0, b, zero, curl1, b3, bi, f, mirror);

    bool has_electric = false, has_mirror = false;
    for (int k = 0; k < nb; ++k) {
      has_electric |= fabs(electric[k] - flat[k]) > 1e-10;
      has_mirror |= fabs(mirror[k] - flat[k]) > 1e-10;
    }
    TEST_CHECK(has_electric);
    if (vdim == 2) {
      // Mirror acceleration requires the magnetic-moment coordinate mu.
      TEST_CHECK(has_mirror);
    }

    // Repeat the checks for fluxes through a v_parallel cell face.
    struct gkyl_dg_vol_geom dg[2] = {{.Jc = 1.0}, {.Jc = 1.0}};
    struct gkyl_gk_dg_vol_geom gk[2] = {{.B3 = 1.0, .bmag = 2.0}, {.B3 = 1.0, .bmag = 2.0}};
    double prime[] = {1.0, 1.0};
    double (*surf)(const double *, const double *, const double *, const double *, const double *, const double *, double, double, const struct gkyl_dg_vol_geom *, const struct gkyl_gk_dg_vol_geom *, const double *, const double *, const double *, const double *, double *) =
      vdim == 1 ? gk_collisionless_flux_surfvpar_1x1v_ser_p1 :
                  gk_collisionless_flux_surfvpar_1x2v_ser_p1;

    double sa[12] = {0.0}, sb[12] = {0.0}, se[12] = {0.0}, sm[12] = {0.0}, sf[12] = {0.0};
    surf(w, dx, prime, prime, vm, vs, 1, 1, dg, gk, b, phi, f, f, sa);
    // Perturb the curvature data at both configuration-space quadrature nodes.
    for (int q = 0; q < 2; ++q) {
      for (int d = 0; d < 3; ++d) {
        gk[q].dualcurlbhat.x[d] = 2.0 + d + q;
      }
    }
    surf(w, dx, prime, prime, vm, vs, 1, 1, dg, gk, b, phi, f, f, sb);
    surf(w, dx, prime, prime, vm, vs, 1, 1, dg, gk, uniform, zero, f, f, sf);
    surf(w, dx, prime, prime, vm, vs, 1, 1, dg, gk, uniform, phi, f, f, se);
    surf(w, dx, prime, prime, vm, vs, 1, 1, dg, gk, b, zero, f, f, sm);

    has_electric = has_mirror = false;
    for (int k = 0; k < 12; ++k) {
      TEST_CHECK(gkyl_compare(sa[k], sb[k], 1e-14));
      has_electric |= fabs(se[k] - sf[k]) > 1e-10;
      has_mirror |= fabs(sm[k] - sf[k]) > 1e-10;
    }
    TEST_CHECK(has_electric);
    if (vdim == 2) {
      TEST_CHECK(has_mirror);
    }
  }
}

#ifdef GKYL_HAVE_CUDA

#endif

TEST_LIST = {
  {"dg_gyrokinetic_parallel_only_1x", test_dg_gyrokinetic_parallel_only_1x},
  {"dg_gyrokinetic_ho", test_dg_gyrokinetic_ho},
#ifdef GKYL_HAVE_CUDA
#endif
  {NULL, NULL}
};
