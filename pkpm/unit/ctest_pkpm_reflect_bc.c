// Reflecting wall of the PKPM distributions (gkyl_bc_pkpm_species_reflect) on the eight blocks
// [F, G, M0_perp (x, y, z), N_perp (x, y, z)].
// 1. Mirror image: with the ghost field direction b' = s R b, R = I - 2 n n (s = +1: PEC image,
//    s = -1: symmetry image), every ghost block is the skin block mirrored in the wall direction,
//    with w -> s w, and the wall-normal components of M0_perp and N_perp negated; for oblique b,
//    every wall and both edges, in 1x, 2x and 3x.
// 2. Wall fluxes: for a skewed, anisotropic distribution with a flow, the mass, energy and
//    tangential momentum fluxes through a reflecting x wall (the moments of the kinetic flux
//    with the mirrored fluid state) vanish for both images at oblique b; with a copied field
//    (b' = b) the mass flux vanishes at every angle, and all three vanish for b along and across
//    the normal, where that image is a mirror.
#include <math.h>
#include <stdbool.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_bc_basic.h>
#include <gkyl_bc_pkpm_species_reflect.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

enum { NUM_BLOCKS = 8 };

static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  gkyl_array_clear(a, 0.0);
  return a;
}

// Constant configuration-space field: coefficient 0 of every component, others zero.
static void
set_const(
  struct gkyl_array *arr, const struct gkyl_range *range, int ncomp, int nc, const double *val,
  int cdim
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    double *a = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
    for (int c = 0; c < ncomp; ++c) {
      for (int k = 0; k < nc; ++k) {
        a[c * nc + k] = k == 0 ? val[c] * pow(sqrt(2.0), cdim) : 0.0;
      }
    }
  }
}

// R b for the wall normal along dir, times s.
static void
mirror_b(const double *b, int dir, double s, double *out)
{
  for (int i = 0; i < 3; ++i) {
    out[i] = s * (i == dir ? -b[i] : b[i]);
  }
}

static void
test_image(int cdim)
{
  int pdim = cdim + 1;
  double lower[GKYL_MAX_DIM] = {-1.0, -1.0, -1.0}, upper[GKYL_MAX_DIM] = {1.0, 1.0, 1.0};
  int cells[GKYL_MAX_DIM] = {3, 3, 3};
  lower[cdim] = -4.0;
  upper[cdim] = 4.0;
  cells[cdim] = 4;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_hybrid(&basis, cdim, 1);
  int np = basis.num_basis, nc = 1 << cdim;

  int ghost[GKYL_MAX_DIM] = {1, 1, 1, 0};
  ghost[cdim] = 0;
  struct gkyl_range local, local_ext, clocal, clocal_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&cgrid, ghost, &clocal_ext, &clocal);

  struct gkyl_array *f = mkarr(NUM_BLOCKS * np, local_ext.volume);
  struct gkyl_array *bvar = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_range_iter iter;
  const double probes[3][GKYL_MAX_DIM] = {
    {0.3, -0.6, 0.2, 0.8}, {-0.9, 0.4, -0.5, -0.1}, {0.1, 0.7, 0.9, -0.7}
  };
  const double b0[3] = {0.48, -0.6, 0.64}; // oblique to every wall
  double max_err = 0.0;
  for (int dir = 0; dir < cdim; ++dir) {
    for (int e = 0; e < 2; ++e) {
      for (int image = 0; image < 2; ++image) {
        double s = image == 0 ? 1.0 : -1.0, bm[3];
        mirror_b(b0, dir, s, bm);
        // skin and interior distributions; b in the interior and the image b in the ghost cells
        gkyl_array_clear(f, 0.0);
        gkyl_range_iter_init(&iter, &local);
        while (gkyl_range_iter_next(&iter)) {
          long loc = gkyl_range_idx(&local, iter.idx);
          double *fc = gkyl_array_fetch(f, loc);
          for (int k = 0; k < NUM_BLOCKS * np; ++k) {
            fc[k] = sin(1.3 * k + 0.7 * loc) + 0.1 * k;
          }
        }
        set_const(bvar, &clocal_ext, 3, nc, bm, cdim);
        set_const(bvar, &clocal, 3, nc, b0, cdim);

        enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
        struct gkyl_range skin_r, ghost_r;
        gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);
        struct gkyl_bc_pkpm_species_reflect *bc =
          gkyl_bc_pkpm_species_reflect_inew(&(struct gkyl_bc_pkpm_species_reflect_inp){
            .dir = dir,
            .edge = edge,
            .cdim = cdim,
            .phase_basis = &basis,
            .phase_basis_on_dev = &basis,
            .conf_num_basis = nc,
            .num_blocks = NUM_BLOCKS,
            .skin_r = &skin_r,
            .ghost_r = &ghost_r,
            .phase_range_ext = &local_ext,
            .conf_range_ext = &clocal_ext,
            .use_gpu = false,
          });
        gkyl_bc_pkpm_species_reflect_advance(bc, bvar, f);

        bool wflip = s < 0.0;
        gkyl_range_iter_init(&iter, &ghost_r);
        while (gkyl_range_iter_next(&iter)) {
          int sidx[GKYL_MAX_DIM];
          gkyl_copy_int_arr(pdim, iter.idx, sidx);
          sidx[dir] = e == 0 ? local.lower[dir] : local.upper[dir];
          bool corner = false; // ghost cells in another direction hold no interior data
          for (int d = 0; d < cdim; ++d) {
            corner = corner ||
                     (d != dir && (iter.idx[d] < local.lower[d] || iter.idx[d] > local.upper[d]));
          }
          if (corner) {
            continue;
          }
          if (wflip) {
            sidx[cdim] = local.lower[cdim] + local.upper[cdim] - iter.idx[cdim];
          }
          const double *fg = gkyl_array_cfetch(f, gkyl_range_idx(&local_ext, iter.idx));
          const double *fs = gkyl_array_cfetch(f, gkyl_range_idx(&local_ext, sidx));
          for (int p = 0; p < 3; ++p) {
            double zg[GKYL_MAX_DIM], zs[GKYL_MAX_DIM];
            for (int d = 0; d < pdim; ++d) {
              zg[d] = probes[p][d];
              zs[d] = (d == dir || (wflip && d == cdim)) ? -probes[p][d] : probes[p][d];
            }
            for (int k = 0; k < NUM_BLOCKS; ++k) {
              double sign = (k >= 2 && (k - 2) % 3 == dir) ? -1.0 : 1.0;
              double vg = basis.eval_expand(zg, &fg[k * np]);
              double vs = basis.eval_expand(zs, &fs[k * np]);
              max_err = fmax(max_err, fabs(vg - sign * vs));
            }
          }
        }
        gkyl_bc_pkpm_species_reflect_release(bc);
      }
    }
  }
  TEST_CHECK(max_err < 1e-12);
  TEST_MSG(
    "cdim = %d: largest mismatch of a ghost value against the mirror image %g", cdim, max_err
  );
  gkyl_array_release(bvar);
  gkyl_array_release(f);
}

// --------------------------------------------------------------------------------------------
// Wall fluxes in 1x: a reflecting upper x wall, uniform fluid state, a skewed distribution.
struct wall_ctx {
  double angle;
};

static void
eval_dist(double t, const double *xn, double *restrict fout, void *ctx)
{
  double w = xn[1];
  // two drifting Maxwellians (skewed, heat-flux carrying) and a w-dependent anisotropy
  double F = 0.6 * exp(-0.5 * (w - 0.8) * (w - 0.8) / 0.6) / sqrt(2.0 * M_PI * 0.6) +
             0.4 * exp(-0.5 * (w + 0.5) * (w + 0.5) / 1.2) / sqrt(2.0 * M_PI * 1.2);
  fout[0] = F;
  fout[1] = (1.3 + 0.2 * w) * F;
}

static struct gkyl_vlasov_position_map *
identity_pos_map(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *clocal,
  const struct gkyl_range *clocal_ext, const struct gkyl_basis *cbasis
)
{
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM] = {0};
  return gkyl_vlasov_position_map_new(cgrid, clocal, clocal_ext, cbasis, inp, false);
}

// Fluid fluxes [F_rho, F_x, F_y, F_z, F_E] through the upper x wall for b = (cos a, sin a, 0)
// in the interior; image 0: ghost b' = R b (PEC), 1: b' = -R b (symmetry), 2: b' = b (copied).
// absorb: zero ghost distributions instead of the reflected ones, with the fluid state and the field
// continued into the ghost cell (the outgoing flux of the wall cell; the face flux averages them).
void
pkpm_reflect_wall_flux(double angle, int image, bool absorb, double *out)
{
  int cdim = 1, nv = 64;
  double wmax = 8.0;
  double lower[2] = {0.0, -wmax}, upper[2] = {1.0, wmax};
  int cells[2] = {2, nv}, ghost[2] = {1, 0};
  struct gkyl_rect_grid cgrid, pgrid;
  gkyl_rect_grid_init(&cgrid, 1, lower, upper, cells);
  gkyl_rect_grid_init(&pgrid, 2, lower, upper, cells);
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_hybrid(&pbasis, 1, 1);
  int nc = cbasis.num_basis, np = pbasis.num_basis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext;
  gkyl_create_grid_ranges(&cgrid, ghost, &clocal_ext, &clocal);
  gkyl_create_grid_ranges(&pgrid, ghost, &plocal_ext, &plocal);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&cgrid, &clocal, &clocal_ext, &cbasis);

  double rho = 1.2, u[3] = {0.3, -0.2, 0.15}, xi = 1.0, theta0 = 1.3;
  // cos(pi/2) is not zero in floating point: an angle of pi/2 means a field exactly along the wall
  double b[3] = {fabs(cos(angle)) < 1e-14 ? 0.0 : cos(angle), sin(angle), 0.0}, bg[3];
  if (image == 2) {
    for (int i = 0; i < 3; ++i) {
      bg[i] = b[i];
    }
  } else {
    mirror_b(b, 0, image == 0 ? 1.0 : -1.0, bg);
  }
  double ug[3] = {-u[0], u[1], u[2]};
  double E = 0.5 * rho * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) + 1.5 * rho * xi * xi;
  double q[5] = {rho, rho * u[0], rho * u[1], rho * u[2], E};
  if (absorb) {
    for (int i = 0; i < 3; ++i) {
      bg[i] = b[i];
      ug[i] = u[i];
    }
  }
  double qg[5] = {rho, rho * ug[0], rho * ug[1], rho * ug[2], E};
  double xiv[3] = {xi, xi * xi, rho * xi * xi};

  struct gkyl_array *fluid = mkarr(5 * nc, clocal_ext.volume),
                    *uarr = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_array *bvar = mkarr(3 * nc, clocal_ext.volume),
                    *xiarr = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_array *theta = mkarr(nc, clocal_ext.volume);
  struct gkyl_array *f = mkarr(NUM_BLOCKS * np, plocal_ext.volume);
  // ghost values first (the whole extended range), then the interior
  set_const(fluid, &clocal_ext, 5, nc, qg, cdim);
  set_const(uarr, &clocal_ext, 3, nc, ug, cdim);
  set_const(bvar, &clocal_ext, 3, nc, bg, cdim);
  set_const(fluid, &clocal, 5, nc, q, cdim);
  set_const(uarr, &clocal, 3, nc, u, cdim);
  set_const(bvar, &clocal, 3, nc, b, cdim);
  set_const(xiarr, &clocal_ext, 3, nc, xiv, cdim);
  set_const(theta, &clocal_ext, 1, nc, &theta0, cdim);
  struct gkyl_array *fg2 = mkarr(2 * np, plocal_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&pgrid, &pbasis, 8, 2, eval_dist, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &plocal, fg2);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_set_offset(f, 1.0, fg2, 0);
  gkyl_array_release(fg2);

  if (!absorb) {
    struct gkyl_range skin_r, ghost_r;
    gkyl_skin_ghost_ranges(&skin_r, &ghost_r, 0, GKYL_UPPER_EDGE, &plocal_ext, ghost);
    struct gkyl_bc_pkpm_species_reflect *bc =
      gkyl_bc_pkpm_species_reflect_inew(&(struct gkyl_bc_pkpm_species_reflect_inp){
        .dir = 0,
        .edge = GKYL_UPPER_EDGE,
        .cdim = cdim,
        .phase_basis = &pbasis,
        .phase_basis_on_dev = &pbasis,
        .conf_num_basis = nc,
        .num_blocks = NUM_BLOCKS,
        .skin_r = &skin_r,
        .ghost_r = &ghost_r,
        .phase_range_ext = &plocal_ext,
        .conf_range_ext = &clocal_ext,
        .use_gpu = false,
      });
    gkyl_bc_pkpm_species_reflect_advance(bc, bvar, f);
    gkyl_bc_pkpm_species_reflect_release(bc);
  }

  struct gkyl_dg_pkpm_remap_flux_surf *up =
    gkyl_dg_pkpm_remap_flux_surf_inew(&(struct gkyl_dg_pkpm_remap_flux_surf_inp){
      .phase_grid = &pgrid,
      .conf_basis = &cbasis,
      .phase_basis = &pbasis,
      .pos_map = pm,
      .num_face_nodes = 3,
      .use_gpu = false,
    });
  int nk = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(up);
  struct gkyl_array *flux = mkarr(5, clocal_ext.volume), *cfl = mkarr(1, clocal_ext.volume);
  // the kinetic flux, the kinetic rate and the blend weight are computed alongside (not checked)
  struct gkyl_array *kin_flux = mkarr(2 * NUM_BLOCKS * nk, plocal_ext.volume);
  struct gkyl_array *kin_cfl = mkarr(1, plocal_ext.volume), *blend = mkarr(1, clocal_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(cfl, 0.0);
  gkyl_array_clear(kin_flux, 0.0);
  gkyl_array_clear(kin_cfl, 0.0);
  gkyl_array_clear(blend, 0.0);
  gkyl_dg_pkpm_remap_flux_surf_advance(
    up, &clocal, &clocal_ext, &plocal, &plocal_ext, bvar, fluid, uarr, xiarr, theta, f, kin_flux,
    flux, kin_cfl, cfl, blend
  );
  int gidx[1] = {clocal.upper[0] + 1};
  const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&clocal_ext, gidx));
  for (int c = 0; c < 5; ++c) {
    out[c] = fl[c];
  }
  gkyl_dg_pkpm_remap_flux_surf_release(up);
  gkyl_vlasov_position_map_release(pm);
  struct gkyl_array *arrs[] = {fluid, uarr, bvar, xiarr, theta, f, flux, cfl, kin_flux, kin_cfl, blend};
  for (int i = 0; i < 11; ++i) {
    gkyl_array_release(arrs[i]);
  }
}

// Leak of a component through the wall relative to the outgoing one-sided flux.
static double
leak(const double *refl, const double *abs, int c)
{
  return fabs(refl[c]) / fmax(fabs(abs[c]), 1e-300);
}

static void
test_wall_flux_images()
{
  const double angles[3] = {0.3, 0.7, 1.2};
  double worst = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int image = 0; image < 2; ++image) {
      double refl[5], abs[5];
      pkpm_reflect_wall_flux(angles[a], image, false, refl);
      pkpm_reflect_wall_flux(angles[a], image, true, abs);
      int comps[4] = {0, 2, 3, 4}; // mass, tangential momenta, energy
      for (int i = 0; i < 4; ++i) {
        worst = fmax(worst, leak(refl, abs, comps[i]));
      }
    }
  }
  TEST_CHECK(worst < 1e-12);
  TEST_MSG(
    "largest relative leak of mass, tangential momentum or energy through the wall %g", worst
  );
}

static void
test_wall_flux_copied_field()
{
  // the mass flux vanishes at every angle (w is flipped for b_n != 0)
  const double oblique[3] = {0.3, 0.7, 1.2};
  double worst_mass = 0.0;
  for (int a = 0; a < 3; ++a) {
    double refl[5], abs[5];
    pkpm_reflect_wall_flux(oblique[a], 2, false, refl);
    pkpm_reflect_wall_flux(oblique[a], 2, true, abs);
    worst_mass = fmax(worst_mass, leak(refl, abs, 0));
  }
  TEST_CHECK(worst_mass < 1e-12);
  TEST_MSG("largest relative mass leak with a copied oblique field %g", worst_mass);
  // b' = b is a mirror image (s = -1) for b along the normal and (s = +1) across it
  const double angles[2] = {0.0, 0.5 * M_PI};
  double worst = 0.0;
  for (int a = 0; a < 2; ++a) {
    double refl[5], abs[5];
    pkpm_reflect_wall_flux(angles[a], 2, false, refl);
    pkpm_reflect_wall_flux(angles[a], 2, true, abs);
    int comps[4] = {0, 2, 3, 4};
    for (int i = 0; i < 4; ++i) {
      if (fabs(abs[comps[i]]) > 1e-12) {
        worst = fmax(worst, leak(refl, abs, comps[i]));
      } else {
        worst = fmax(worst, fabs(refl[comps[i]]));
      }
    }
  }
  TEST_CHECK(worst < 1e-12);
  TEST_MSG("largest relative leak with a copied field along or across the normal %g", worst);
}

void
test_image_1x()
{
  test_image(1);
}
void
test_image_2x()
{
  test_image(2);
}
void
test_image_3x()
{
  test_image(3);
}

TEST_LIST = {
  {"test_image_1x", test_image_1x},
  {"test_image_2x", test_image_2x},
  {"test_image_3x", test_image_3x},
  {"test_wall_flux_images", test_wall_flux_images},
  {"test_wall_flux_copied_field", test_wall_flux_copied_field},
  {NULL, NULL}
};
