// Tests of the absorbing wall and of the position map of the DG Euler equation object:
//   wall_flux_*      : the flux of the particles of the local Maxwellian that move toward the wall
//                      (euler_wall_flux, iso_euler_wall_flux): the effusion flux of a gas at rest,
//                      F(q) for a supersonic flow toward the wall at any gas_gamma, no flux for one
//                      away from it, the mirror symmetry of the two walls, and a non-negative
//                      outgoing mass flux.
//   boundary_surf_1x : the boundary surface term of a uniform state takes the wall flux out through
//                      the wall face and F(q) in through the interior face.
//   wall_integral_2x : the wall flux integrated over the wall face of each boundary cell, with
//                      the physical face length on a position map.
//   mapped_terms_1x  : on a position map the volume, surface and boundary surface terms are those
//                      of the cell with its physical size J dx.
#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>

// Euler flux along x of q = [rho, rho ux, rho uy, rho uz, E].
static void
euler_flux_x(double gas_gamma, const double *q, double *f)
{
  double ux = q[1] / q[0];
  double pr = (gas_gamma - 1.0) * (q[4] - 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) / q[0]);
  f[0] = q[1];
  f[1] = q[1] * ux + pr;
  f[2] = q[2] * ux;
  f[3] = q[3] * ux;
  f[4] = (q[4] + pr) * ux;
}

// State of density rho, flow (un, ut1, ut2) and pressure p.
static void
euler_state(double gas_gamma, double rho, double un, double ut1, double ut2, double p, double *q)
{
  q[0] = rho;
  q[1] = rho * un;
  q[2] = rho * ut1;
  q[3] = rho * ut2;
  q[4] = 0.5 * rho * (un * un + ut1 * ut1 + ut2 * ut2) + p / (gas_gamma - 1.0);
}

static void
test_wall_flux_rest()
{
  // A gas at rest with a tangential drift: mass flux rho vt/sqrt(2 pi), normal momentum flux p/2,
  // energy flux (2 T + m ut^2/2) per particle (gas_gamma = 5/3).
  double gamma = 5.0 / 3.0, rho = 1.3, p = 0.7, ut1 = 0.4, ut2 = -0.2;
  double q[5], f[5];
  euler_state(gamma, rho, 0.0, ut1, ut2, p, q);
  double vt = sqrt(p / rho), gam_flux = rho * vt / sqrt(2.0 * M_PI);
  for (int s = -1; s <= 1; s += 2) {
    euler_wall_flux(gamma, s, q, f);
    TEST_CHECK(gkyl_compare(f[0], s * gam_flux, 1e-14));
    TEST_CHECK(gkyl_compare(f[1], 0.5 * p, 1e-14));
    TEST_CHECK(gkyl_compare(f[2], ut1 * s * gam_flux, 1e-14));
    TEST_CHECK(gkyl_compare(f[3], ut2 * s * gam_flux, 1e-14));
    double e_per_mass = 2.0 * vt * vt + 0.5 * (ut1 * ut1 + ut2 * ut2);
    TEST_CHECK(gkyl_compare(f[4], s * gam_flux * e_per_mass, 1e-14));
  }

  // Isothermal: the same mass and momentum fluxes at the thermal speed.
  double vti = 0.6, fi[4];
  double qi[4] = {rho, 0.0, rho * ut1, rho * ut2};
  iso_euler_wall_flux(vti, 1.0, qi, fi);
  TEST_CHECK(gkyl_compare(fi[0], rho * vti / sqrt(2.0 * M_PI), 1e-14));
  TEST_CHECK(gkyl_compare(fi[1], 0.5 * rho * vti * vti, 1e-14));
  TEST_CHECK(gkyl_compare(fi[2], ut1 * fi[0], 1e-14));
}

static void
test_wall_flux_supersonic()
{
  // Supersonic toward the wall (10 thermal speeds): F(q), for gas_gamma = 5/3, 1.4 and 3; away from
  // the wall: no flux.
  double gammas[3] = {5.0 / 3.0, 1.4, 3.0};
  for (int g = 0; g < 3; ++g) {
    double gamma = gammas[g], rho = 0.8, p = 0.5;
    double vt = sqrt(p / rho), un = 10.0 * sqrt(2.0) * vt;
    double q[5], f[5], fe[5];
    euler_state(gamma, rho, un, 0.3, 0.1, p, q);
    euler_wall_flux(gamma, 1.0, q, f);
    euler_flux_x(gamma, q, fe);
    for (int k = 0; k < 5; ++k) {
      TEST_CHECK(gkyl_compare(f[k], fe[k], 1e-12));
    }
    // the lower wall with the flow reversed: the mirror image
    double qm[5], fm[5];
    euler_state(gamma, rho, -un, 0.3, 0.1, p, qm);
    euler_wall_flux(gamma, -1.0, qm, fm);
    TEST_CHECK(gkyl_compare(fm[0], -f[0], 1e-12));
    TEST_CHECK(gkyl_compare(fm[1], f[1], 1e-12));
    TEST_CHECK(gkyl_compare(fm[2], -f[2], 1e-12));
    TEST_CHECK(gkyl_compare(fm[4], -f[4], 1e-12));
    // away from the wall
    euler_wall_flux(gamma, 1.0, qm, f);
    for (int k = 0; k < 5; ++k) {
      TEST_CHECK(fabs(f[k]) < 1e-40);
    }
  }
}

static void
test_wall_flux_monotone()
{
  // The outgoing mass flux is non-negative and grows with the flow toward the wall.
  double gamma = 5.0 / 3.0, rho = 1.0, p = 1.0, prev = -1.0;
  for (int i = -60; i <= 60; ++i) {
    double q[5], f[5];
    euler_state(gamma, rho, 0.1 * i, 0.0, 0.0, p, q);
    euler_wall_flux(gamma, 1.0, q, f);
    TEST_CHECK(f[0] >= 0.0);
    TEST_CHECK(f[0] >= prev);
    TEST_CHECK(f[4] >= 0.0);
    prev = f[0];
  }
}

// Grid, ranges and geometry of a p = 1 test, and an equation object on a position map (identity
// where pmap is NULL).
struct wall_test {
  struct gkyl_basis basis;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  struct gkyl_wave_geom *geom;
  struct gkyl_wv_eqn *wv_eqn;
  struct gkyl_vlasov_position_map *pos_map;
  struct gkyl_array *u, *p;
  struct gkyl_dg_eqn *eqn;
};

static struct wall_test
wall_test_new(int cdim, const int *cells, evalf_t pmap[GKYL_MAX_CDIM])
{
  struct wall_test wt;
  double lower[3] = {0.0, 0.0, 0.0}, upper[3] = {1.0, 2.0, 1.0};
  int nghost[3] = {1, 1, 1};
  gkyl_cart_modal_serendip(&wt.basis, cdim, 1);
  gkyl_rect_grid_init(&wt.grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&wt.grid, nghost, &wt.local_ext, &wt.local);
  wt.geom = gkyl_wave_geom_new(&wt.grid, &wt.local_ext, NULL, NULL, false);
  wt.wv_eqn = gkyl_wv_euler_inew(&(struct gkyl_wv_euler_inp){
    .gas_gamma = 5.0 / 3.0,
    .rp_type = WV_EULER_RP_HLLC,
  });
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    inp_pmap[d].eval_pmap = pmap ? pmap[d] : 0;
  }
  wt.pos_map =
    gkyl_vlasov_position_map_new(&wt.grid, &wt.local, &wt.local_ext, &wt.basis, inp_pmap, false);
  wt.u = gkyl_array_new(GKYL_DOUBLE, 3 * wt.basis.num_basis, wt.local_ext.volume);
  wt.p = gkyl_array_new(GKYL_DOUBLE, wt.basis.num_basis, wt.local_ext.volume);
  wt.eqn = gkyl_dg_euler_inew(&(struct gkyl_dg_euler_inp){
    .cbasis = &wt.basis,
    .conf_range = &wt.local,
    .wv_eqn = wt.wv_eqn,
    .geom = wt.geom,
    .pos_map = wt.pos_map,
    .u = wt.u,
    .p = wt.p,
    .use_gpu = false,
  });
  return wt;
}

static void
wall_test_release(struct wall_test *wt)
{
  gkyl_dg_eqn_release(wt->eqn);
  gkyl_array_release(wt->u);
  gkyl_array_release(wt->p);
  gkyl_vlasov_position_map_release(wt->pos_map);
  gkyl_wv_eqn_release(wt->wv_eqn);
  gkyl_wave_geom_release(wt->geom);
}

// DG coefficients of a uniform state (orthonormal basis: the first function is 2^(-cdim/2)).
static void
uniform_coeffs(int cdim, int nb, int ncomp, const double *q, double *coeffs)
{
  double c0 = pow(sqrt(2.0), cdim);
  for (int k = 0; k < ncomp * nb; ++k) {
    coeffs[k] = 0.0;
  }
  for (int k = 0; k < ncomp; ++k) {
    coeffs[k * nb] = q[k] * c0;
  }
}

static void
test_boundary_surf_1x()
{
  struct wall_test wt = wall_test_new(1, (int[]){8}, 0);
  int nb = wt.basis.num_basis;
  double gamma = 5.0 / 3.0, q[5], fe[5], fw[5];
  euler_state(gamma, 1.1, 0.3, 0.2, -0.1, 0.9, q);
  euler_flux_x(gamma, q, fe);
  double fluid[10];
  uniform_coeffs(1, nb, 5, q, fluid);
  double dx = wt.grid.dx[0], xc[1] = {0.5};

  for (int edge = -1; edge <= 1; edge += 2) {
    int idx_skin[1] = {edge == -1 ? 1 : 8}, idx_edge[1] = {idx_skin[0] - edge};
    double out[10] = {0.0};
    wt.eqn->boundary_surf_term(
      wt.eqn, 0, xc, xc, wt.grid.dx, wt.grid.dx, idx_edge, idx_skin, edge, fluid, fluid, out
    );
    euler_wall_flux(gamma, edge, q, fw);
    // rate of the cell averages: (flux in - flux out)/dx through the lower and upper faces
    for (int k = 0; k < 5; ++k) {
      double fl = edge == -1 ? fw[k] : fe[k], fr = edge == -1 ? fe[k] : fw[k];
      TEST_CHECK(gkyl_compare(out[k * nb] / sqrt(2.0), (fl - fr) / dx, 1e-12));
    }
  }
  wall_test_release(&wt);
}

// Position maps: x_p = 3 x_c (J = 3) and y_p = y_c + y_c^2/4 on [0, 2] (J from 1 to 2).
static void
pmap_scale3(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = 3.0 * xc[0];
}

static void
pmap_quad(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = xc[0] + 0.25 * xc[0] * xc[0];
}

static void
test_wall_integral_2x()
{
  // Walls in x on a grid mapped in y: the integral over the wall face of a boundary cell is the
  // wall flux times its physical length J_y dy.
  evalf_t pmap[GKYL_MAX_CDIM] = {0, pmap_quad};
  struct wall_test wt = wall_test_new(2, (int[]){4, 5}, pmap);
  int nb = wt.basis.num_basis;
  double gamma = 5.0 / 3.0, q[5], fw[5];
  euler_state(gamma, 0.9, -0.4, 0.1, 0.3, 1.2, q);
  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5 * nb, wt.local_ext.volume);
  struct gkyl_array *flux_int = gkyl_array_new(GKYL_DOUBLE, 5, wt.local_ext.volume);
  double coeffs[20];
  uniform_coeffs(2, nb, 5, q, coeffs);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &wt.local_ext);
  while (gkyl_range_iter_next(&iter)) {
    double *fd = gkyl_array_fetch(fluid, gkyl_range_idx(&wt.local_ext, iter.idx));
    for (int k = 0; k < 5 * nb; ++k) {
      fd[k] = coeffs[k];
    }
  }

  struct gkyl_range skin_lo, skin_up, ghost;
  gkyl_skin_ghost_ranges(&skin_lo, &ghost, 0, GKYL_LOWER_EDGE, &wt.local_ext, (int[]){1, 1});
  gkyl_skin_ghost_ranges(&skin_up, &ghost, 0, GKYL_UPPER_EDGE, &wt.local_ext, (int[]){1, 1});
  double length = 0.0;
  for (int edge = -1; edge <= 1; edge += 2) {
    const struct gkyl_range *skin = edge == -1 ? &skin_lo : &skin_up;
    gkyl_dg_euler_wall_flux(wt.eqn, &wt.grid, 0, edge, skin, fluid, flux_int);
    euler_wall_flux(gamma, edge, q, fw);
    length = 0.0;
    gkyl_range_iter_init(&iter, skin);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(skin, iter.idx);
      const double *fi = gkyl_array_cfetch(flux_int, loc);
      double y_lo = (iter.idx[1] - 1) * wt.grid.dx[1], y_up = iter.idx[1] * wt.grid.dx[1];
      double len = (y_up + 0.25 * y_up * y_up) - (y_lo + 0.25 * y_lo * y_lo);
      length += len;
      for (int k = 0; k < 5; ++k) {
        TEST_CHECK(gkyl_compare(fi[k], fw[k] * len, 1e-12));
      }
    }
  }
  TEST_CHECK(gkyl_compare(length, 3.0, 1e-13)); // the mapped domain [0, 3] in y

  gkyl_array_release(fluid);
  gkyl_array_release(flux_int);
  wall_test_release(&wt);
}

static void
test_mapped_terms_1x()
{
  // x_p = 3 x_c: every term of a cell is that of the uniform cell divided by J = 3.
  evalf_t pmap[GKYL_MAX_CDIM] = {pmap_scale3};
  struct wall_test wm = wall_test_new(1, (int[]){8}, pmap);
  struct wall_test wu = wall_test_new(1, (int[]){8}, 0);
  int nb = wm.basis.num_basis;
  // one set of terms for both meshes: the kernels take the Jacobian of the cell
  TEST_CHECK(wm.eqn->vol_term == wu.eqn->vol_term);
  TEST_CHECK(wm.eqn->surf_term == wu.eqn->surf_term);
  TEST_CHECK(wm.eqn->boundary_surf_term == wu.eqn->boundary_surf_term);

  // a state with slopes, its u and p (any values: the terms are linear in the kernels' dx factor)
  double fl[10], fc[10], fr[10];
  for (int k = 0; k < 5; ++k) {
    double base[5] = {1.0, 0.2, -0.1, 0.05, 2.0};
    fl[k * nb] = base[k] * sqrt(2.0);
    fl[k * nb + 1] = 0.01 * (k + 1);
    fc[k * nb] = 1.1 * base[k] * sqrt(2.0);
    fc[k * nb + 1] = -0.02 * (k + 1);
    fr[k * nb] = 0.9 * base[k] * sqrt(2.0);
    fr[k * nb + 1] = 0.015 * (k + 1);
  }
  struct wall_test *wts[2] = {&wm, &wu};
  for (int i = 0; i < 2; ++i) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &wts[i]->local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&wts[i]->local_ext, iter.idx);
      double *ud = gkyl_array_fetch(wts[i]->u, loc), *pd = gkyl_array_fetch(wts[i]->p, loc);
      ud[0] = 0.25;
      ud[1] = 0.01;
      ud[2] = ud[3] = ud[4] = ud[5] = 0.0;
      pd[0] = 0.8;
      pd[1] = -0.02;
    }
  }
  int idxl[1] = {3}, idxc[1] = {4}, idxr[1] = {5};
  double xc[1] = {0.5}, out_m[10] = {0.0}, out_u[10] = {0.0};
  double cfl_m = wm.eqn->vol_term(wm.eqn, xc, wm.grid.dx, idxc, fc, out_m);
  double cfl_u = wu.eqn->vol_term(wu.eqn, xc, wu.grid.dx, idxc, fc, out_u);
  TEST_CHECK(gkyl_compare(cfl_m, cfl_u / 3.0, 1e-14));
  wm.eqn->surf_term(
    wm.eqn, 0, xc, xc, xc, wm.grid.dx, wm.grid.dx, wm.grid.dx, idxl, idxc, idxr, fl, fc, fr, out_m
  );
  wu.eqn->surf_term(
    wu.eqn, 0, xc, xc, xc, wu.grid.dx, wu.grid.dx, wu.grid.dx, idxl, idxc, idxr, fl, fc, fr, out_u
  );
  for (int k = 0; k < 10; ++k) {
    TEST_CHECK(gkyl_compare(out_m[k], out_u[k] / 3.0, 1e-13));
  }
  // the boundary cell of the upper wall (its interior neighbor 7)
  int idxe[1] = {7}, idxs[1] = {8};
  double outb_m[10] = {0.0}, outb_u[10] = {0.0};
  wm.eqn->boundary_surf_term(
    wm.eqn, 0, xc, xc, wm.grid.dx, wm.grid.dx, idxe, idxs, 1, fl, fc, outb_m
  );
  wu.eqn->boundary_surf_term(
    wu.eqn, 0, xc, xc, wu.grid.dx, wu.grid.dx, idxe, idxs, 1, fl, fc, outb_u
  );
  for (int k = 0; k < 10; ++k) {
    TEST_CHECK(outb_u[k] != 0.0);
    TEST_CHECK(gkyl_compare(outb_m[k], outb_u[k] / 3.0, 1e-13));
  }
  wall_test_release(&wm);
  wall_test_release(&wu);
}

TEST_LIST = {
  {"wall_flux_rest", test_wall_flux_rest},
  {"wall_flux_supersonic", test_wall_flux_supersonic},
  {"wall_flux_monotone", test_wall_flux_monotone},
  {"boundary_surf_1x", test_boundary_surf_1x},
  {"wall_integral_2x", test_wall_integral_2x},
  {"mapped_terms_1x", test_mapped_terms_1x},
  {NULL, NULL}
};
