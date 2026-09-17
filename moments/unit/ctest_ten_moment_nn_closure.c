// Unit tests for the ten-moment neural-network closure geometry helper.
//
// gkyl_ten_moment_nn_closure_geom_calc() turns a fluid/EM stencil into the
// network input feature vector and a cache of magnetic-field-aligned geometry.
// These tests prescribe the density, pressure tensor, and magnetic field on the
// stencil and check that the field-aligned projections (p_par, p_perp), the
// unit field direction b = B/|B|, and the symmetric gradients come out as
// expected. The construction half (gkyl_ten_moment_nn_closure_construct) is also
// unit-tested below from a prescribed heat flux, independent of the network.

#include <acutest.h>
#include <math.h>

#include <gkyl_rect_grid.h>
#include <gkyl_ten_moment_nn_closure.h>

// Ten-moment conserved variable layout.
#define TM_RHO 0
#define TM_MX 1
#define TM_MY 2
#define TM_MZ 3
#define TM_P11 4
#define TM_P12 5
#define TM_P13 6
#define TM_P22 7
#define TM_P23 8
#define TM_P33 9

// EM field layout (EX,EY,EZ,BX,BY,BZ,phi,psi).
#define EM_BX 3
#define EM_BY 4
#define EM_BZ 5

static const double tol = 1.0e-5;

// Fill a ten-moment cell with zero bulk velocity so the pressure tensor equals
// the conserved P-components directly (p_ij = P_ij - rho u_i u_j = P_ij).
static void
set_fluid(double f[10], double rho, double p11, double p12, double p13, double p22, double p23, double p33)
{
  f[TM_RHO] = rho;
  f[TM_MX] = 0.0; f[TM_MY] = 0.0; f[TM_MZ] = 0.0;
  f[TM_P11] = p11; f[TM_P12] = p12; f[TM_P13] = p13;
  f[TM_P22] = p22; f[TM_P23] = p23; f[TM_P33] = p33;
}

static void
set_em(double e[8], double bx, double by, double bz)
{
  for (int i = 0; i < 8; i++) {
    e[i] = 0.0;
  }
  e[EM_BX] = bx; e[EM_BY] = by; e[EM_BZ] = bz;
}

static struct gkyl_ten_moment_nn_closure*
mk_closure_1d(int poly_order, double dx)
{
  static struct gkyl_rect_grid grid;
  // cells = 1/dx so that the grid spacing is exactly dx.
  int cells = (int)(1.0 / dx + 0.5);
  gkyl_rect_grid_init(&grid, 1, (double[]) { 0.0 }, (double[]) { 1.0 }, (int[]) { cells });
  return gkyl_ten_moment_nn_closure_new( (struct gkyl_ten_moment_nn_closure_inp) {
      .grid = &grid, .poly_order = poly_order, .k0 = 1.0, .ann = 0
    });
}

// Input/output feature counts for each supported configuration.
static void
test_nn_closure_dims(void)
{
  struct gkyl_rect_grid grid1, grid2;
  gkyl_rect_grid_init(&grid1, 1, (double[]) { 0.0 }, (double[]) { 1.0 }, (int[]) { 10 });
  gkyl_rect_grid_init(&grid2, 2, (double[]) { 0.0, 0.0 }, (double[]) { 1.0, 1.0 }, (int[]) { 10, 10 });

  struct gkyl_ten_moment_nn_closure *nn_1d_p1 = gkyl_ten_moment_nn_closure_new(
    (struct gkyl_ten_moment_nn_closure_inp) { .grid = &grid1, .poly_order = 1, .k0 = 1.0, .ann = 0 });
  struct gkyl_ten_moment_nn_closure *nn_1d_p2 = gkyl_ten_moment_nn_closure_new(
    (struct gkyl_ten_moment_nn_closure_inp) { .grid = &grid1, .poly_order = 2, .k0 = 1.0, .ann = 0 });
  struct gkyl_ten_moment_nn_closure *nn_2d_p1 = gkyl_ten_moment_nn_closure_new(
    (struct gkyl_ten_moment_nn_closure_inp) { .grid = &grid2, .poly_order = 1, .k0 = 1.0, .ann = 0 });

  TEST_CHECK(gkyl_ten_moment_nn_closure_n_in(nn_1d_p1) == 6);
  TEST_CHECK(gkyl_ten_moment_nn_closure_n_out(nn_1d_p1) == 4);
  TEST_CHECK(gkyl_ten_moment_nn_closure_n_in(nn_1d_p2) == 9);
  TEST_CHECK(gkyl_ten_moment_nn_closure_n_out(nn_1d_p2) == 6);
  TEST_CHECK(gkyl_ten_moment_nn_closure_n_in(nn_2d_p1) == 12);
  TEST_CHECK(gkyl_ten_moment_nn_closure_n_out(nn_2d_p1) == 8);

  gkyl_ten_moment_nn_closure_release(nn_1d_p1);
  gkyl_ten_moment_nn_closure_release(nn_1d_p2);
  gkyl_ten_moment_nn_closure_release(nn_2d_p1);
}

// Uniform stencil, B along x: b = (1,0,0), p_par = p_xx, p_perp = (p_yy+p_zz)/2,
// all gradients zero.
static void
test_geom_1d_p1_uniform_bx(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  double fL[10], fU[10], eL[8], eU[8];
  set_fluid(fL, 2.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_fluid(fU, 2.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_em(eL, 1.0, 0.0, 0.0);
  set_em(eU, 1.0, 0.0, 0.0);

  const double *fluid_d[2] = { fL, fU };
  const double *em_d[2] = { eL, eU };

  float in[6];
  struct gkyl_ten_moment_nn_closure_geom geom;
  gkyl_ten_moment_nn_closure_geom_calc(nn, fluid_d, em_d, in, &geom);

  TEST_CHECK( fabs(in[0] - 2.0) < tol );  // rho_avg
  TEST_CHECK( fabs(in[1] - 0.0) < tol );  // drho_dx
  TEST_CHECK( fabs(in[2] - 3.0) < tol );  // p_par = p_xx
  TEST_CHECK( fabs(in[3] - 0.0) < tol );  // p_par_dx
  TEST_CHECK( fabs(in[4] - 0.75) < tol ); // p_perp = (p_yy + p_zz)/2
  TEST_CHECK( fabs(in[5] - 0.0) < tol );  // p_perp_dx

  TEST_CHECK( fabs(geom.local_mag[0] - 1.0) < tol );
  TEST_CHECK( fabs(geom.local_mag[1] - 0.0) < tol );
  TEST_CHECK( fabs(geom.local_mag[2] - 0.0) < tol );
  TEST_CHECK( fabs(geom.local_mag_dx[0]) < tol );
  TEST_CHECK( fabs(geom.local_mag_dx[1]) < tol );
  TEST_CHECK( fabs(geom.local_mag_dx[2]) < tol );
  TEST_CHECK( fabs(geom.rho_avg - 2.0) < tol );
  TEST_CHECK( fabs(geom.B_avg[0] - 1.0) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Uniform stencil, B along z: b = (0,0,1), p_par = p_zz, p_perp = (p_xx+p_yy)/2.
static void
test_geom_1d_p1_uniform_bz(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  double fL[10], fU[10], eL[8], eU[8];
  set_fluid(fL, 2.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_fluid(fU, 2.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_em(eL, 0.0, 0.0, 1.0);
  set_em(eU, 0.0, 0.0, 1.0);

  const double *fluid_d[2] = { fL, fU };
  const double *em_d[2] = { eL, eU };

  float in[6];
  struct gkyl_ten_moment_nn_closure_geom geom;
  gkyl_ten_moment_nn_closure_geom_calc(nn, fluid_d, em_d, in, &geom);

  TEST_CHECK( fabs(in[2] - 0.5) < tol ); // p_par = p_zz
  TEST_CHECK( fabs(in[4] - 2.0) < tol ); // p_perp = (p_xx + p_yy)/2 = (3+1)/2
  TEST_CHECK( fabs(geom.local_mag[2] - 1.0) < tol );
  TEST_CHECK( fabs(geom.local_mag[0]) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Uniform stencil, B at 45 deg in the x-y plane: b = (1,1,0)/sqrt(2).
// p_par = 0.5 p_xx + 0.5 p_yy + p_xy, p_perp = (tr(p) - p_par)/2.
static void
test_geom_1d_p1_diagonal_b(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  double pxx = 3.0, pyy = 1.0, pzz = 0.5, pxy = 0.2;
  double fL[10], fU[10], eL[8], eU[8];
  set_fluid(fL, 2.0, pxx, pxy, 0.0, pyy, 0.0, pzz);
  set_fluid(fU, 2.0, pxx, pxy, 0.0, pyy, 0.0, pzz);
  set_em(eL, 1.0, 1.0, 0.0);
  set_em(eU, 1.0, 1.0, 0.0);

  const double *fluid_d[2] = { fL, fU };
  const double *em_d[2] = { eL, eU };

  float in[6];
  struct gkyl_ten_moment_nn_closure_geom geom;
  gkyl_ten_moment_nn_closure_geom_calc(nn, fluid_d, em_d, in, &geom);

  double p_par_exp = 0.5 * pxx + 0.5 * pyy + pxy;
  double p_perp_exp = 0.5 * ((pxx + pyy + pzz) - p_par_exp);
  double inv_sqrt2 = 1.0 / sqrt(2.0);

  TEST_CHECK( fabs(in[2] - p_par_exp) < tol );
  TEST_CHECK( fabs(in[4] - p_perp_exp) < tol );
  TEST_CHECK( fabs(geom.local_mag[0] - inv_sqrt2) < tol );
  TEST_CHECK( fabs(geom.local_mag[1] - inv_sqrt2) < tol );
  TEST_CHECK( fabs(geom.local_mag[2] - 0.0) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Density gradient only: B and pressure uniform along x, density linear across
// the stencil. drho_dx = (rho_U - rho_L)/dx; field-aligned pressure gradients 0.
static void
test_geom_1d_p1_density_gradient(void)
{
  double dx = 0.1;
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, dx);

  double fL[10], fU[10], eL[8], eU[8];
  set_fluid(fL, 1.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_fluid(fU, 3.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_em(eL, 1.0, 0.0, 0.0);
  set_em(eU, 1.0, 0.0, 0.0);

  const double *fluid_d[2] = { fL, fU };
  const double *em_d[2] = { eL, eU };

  float in[6];
  struct gkyl_ten_moment_nn_closure_geom geom;
  gkyl_ten_moment_nn_closure_geom_calc(nn, fluid_d, em_d, in, &geom);

  TEST_CHECK( fabs(in[0] - 2.0) < tol );           // rho_avg = (1+3)/2
  TEST_CHECK( fabs(in[1] - (2.0 / dx)) < tol );    // drho_dx = (3-1)/dx = 20
  TEST_CHECK( fabs(in[3] - 0.0) < tol );           // p_par_dx
  TEST_CHECK( fabs(in[5] - 0.0) < tol );           // p_perp_dx

  gkyl_ten_moment_nn_closure_release(nn);
}

// Pressure gradient with B uniform along x. With b = (1,0,0) and uniform b,
// p_par_dx = d(p_xx)/dx and p_perp_dx = 0.5( d tr(p)/dx - p_par_dx ).
static void
test_geom_1d_p1_pressure_gradient(void)
{
  double dx = 0.1;
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, dx);

  double fL[10], fU[10], eL[8], eU[8];
  set_fluid(fL, 2.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_fluid(fU, 2.0, 4.0, 0.0, 0.0, 2.0, 0.0, 0.5);
  set_em(eL, 1.0, 0.0, 0.0);
  set_em(eU, 1.0, 0.0, 0.0);

  const double *fluid_d[2] = { fL, fU };
  const double *em_d[2] = { eL, eU };

  float in[6];
  struct gkyl_ten_moment_nn_closure_geom geom;
  gkyl_ten_moment_nn_closure_geom_calc(nn, fluid_d, em_d, in, &geom);

  // p_xx avg = 3, p_yy avg = 1.5, p_zz = 0.5 -> p_par = 3, p_perp = 0.5(5-3) = 1.
  TEST_CHECK( fabs(in[2] - 3.0) < tol );
  TEST_CHECK( fabs(in[4] - 1.0) < tol );
  // d p_xx/dx = (4-2)/dx = 20 -> p_par_dx = 20.
  TEST_CHECK( fabs(in[3] - (2.0 / dx)) < tol );
  // d tr/dx = (dPxx + dPyy)/dx = (2+1)/dx = 30 -> p_perp_dx = 0.5(30-20) = 5.
  TEST_CHECK( fabs(in[5] - 5.0) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Vanishing magnetic field: the closure falls back to b = (1,0,0), so p_par
// reduces to p_xx (avoids a divide-by-zero in the field direction).
static void
test_geom_1d_p1_zero_b(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  double fL[10], fU[10], eL[8], eU[8];
  set_fluid(fL, 2.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_fluid(fU, 2.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.5);
  set_em(eL, 0.0, 0.0, 0.0);
  set_em(eU, 0.0, 0.0, 0.0);

  const double *fluid_d[2] = { fL, fU };
  const double *em_d[2] = { eL, eU };

  float in[6];
  struct gkyl_ten_moment_nn_closure_geom geom;
  gkyl_ten_moment_nn_closure_geom_calc(nn, fluid_d, em_d, in, &geom);

  TEST_CHECK( fabs(geom.local_mag[0] - 1.0) < tol );
  TEST_CHECK( fabs(in[2] - 3.0) < tol );  // p_par = p_xx
  TEST_CHECK( fabs(in[4] - 0.75) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Helper: zero a geometry cache and set a uniform field direction b along axis
// `ax` (0=x,1=y,2=z), with |B| = 1 (so the divQ B-prefactor is unity).
static void
set_geom_uniform_b(struct gkyl_ten_moment_nn_closure_geom *g, int ax)
{
  for (int i = 0; i < 3; i++) {
    g->local_mag[i] = 0.0; g->local_mag_dx[i] = 0.0; g->local_mag_dy[i] = 0.0; g->B_avg[i] = 0.0;
  }
  g->rho_avg = 1.0;
  for (int i = 0; i < 6; i++) {
    g->p_avg[i] = 0.0;
  }
  g->local_mag[ax] = 1.0;
  g->B_avg[ax] = 1.0;
}

// Coupling/sign test: with b = x and a prescribed heat flux (q_par, q_par_dx,
// q_perp, q_perp_dx), the pressure-tensor source must be d(P_ij)/dt = -d(q)/dx,
// with q_par -> Pxx and q_perp -> Pyy, Pzz (and no off-diagonal/mass/momentum).
static void
test_consume_1d_p1_sign_and_mapping(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  struct gkyl_ten_moment_nn_closure_geom geom;
  set_geom_uniform_b(&geom, 0); // b = x, uniform

  // pred = [q_par, q_par_dx, q_perp, q_perp_dx]
  float pred[4] = { 0.3f, 0.5f, 0.2f, 0.1f };
  double rhs[10] = { 0.0 };
  gkyl_ten_moment_nn_closure_construct(nn, &geom, pred, rhs);

  TEST_CHECK( fabs(rhs[TM_P11] - (-0.5)) < tol ); // -q_par_dx
  TEST_MSG("rhs[P11]=%g expected %g", rhs[TM_P11], -0.5);
  TEST_CHECK( fabs(rhs[TM_P22] - (-0.1)) < tol ); // -q_perp_dx
  TEST_CHECK( fabs(rhs[TM_P33] - (-0.1)) < tol ); // -q_perp_dx
  TEST_CHECK( fabs(rhs[TM_P12]) < tol );
  TEST_CHECK( fabs(rhs[TM_P13]) < tol );
  TEST_CHECK( fabs(rhs[TM_P23]) < tol );
  TEST_CHECK( fabs(rhs[TM_RHO]) < tol );
  TEST_CHECK( fabs(rhs[TM_MX]) < tol );
  TEST_CHECK( fabs(rhs[TM_MY]) < tol );
  TEST_CHECK( fabs(rhs[TM_MZ]) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Sign robustness: a negative parallel-flux gradient must flip the Pxx source.
static void
test_consume_1d_p1_sign_flip(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  struct gkyl_ten_moment_nn_closure_geom geom;
  set_geom_uniform_b(&geom, 0); // b = x

  float pred[4] = { 0.0f, -0.7f, 0.0f, 0.25f }; // q_par_dx<0, q_perp_dx>0
  double rhs[10] = { 0.0 };
  gkyl_ten_moment_nn_closure_construct(nn, &geom, pred, rhs);

  TEST_CHECK( fabs(rhs[TM_P11] - 0.7) < tol );    // -(-0.7)
  TEST_CHECK( fabs(rhs[TM_P22] - (-0.25)) < tol );
  TEST_CHECK( fabs(rhs[TM_P33] - (-0.25)) < tol );

  gkyl_ten_moment_nn_closure_release(nn);
}

// Divergence of a constant heat flux is zero: a uniform q (zero predicted
// gradients) must produce no pressure source, even with q != 0. This is why the
// network's constant offset in static regions does not inject a spurious source
// provided it also predicts q' ~ 0 there.
static void
test_consume_1d_p1_uniform_q_zero_source(void)
{
  struct gkyl_ten_moment_nn_closure *nn = mk_closure_1d(1, 0.1);

  struct gkyl_ten_moment_nn_closure_geom geom;
  set_geom_uniform_b(&geom, 0); // b = x

  float pred[4] = { 0.4f, 0.0f, 0.15f, 0.0f }; // nonzero q, zero gradients
  double rhs[10] = { 0.0 };
  gkyl_ten_moment_nn_closure_construct(nn, &geom, pred, rhs);

  for (int n = 0; n < 10; n++) {
    TEST_CHECK( fabs(rhs[n]) < tol );
    TEST_MSG("rhs[%d] = %g (expected 0)", n, rhs[n]);
  }

  gkyl_ten_moment_nn_closure_release(nn);
}

TEST_LIST = {
  { "nn_closure_dims", test_nn_closure_dims },
  { "geom_1d_p1_uniform_bx", test_geom_1d_p1_uniform_bx },
  { "geom_1d_p1_uniform_bz", test_geom_1d_p1_uniform_bz },
  { "geom_1d_p1_diagonal_b", test_geom_1d_p1_diagonal_b },
  { "geom_1d_p1_density_gradient", test_geom_1d_p1_density_gradient },
  { "geom_1d_p1_pressure_gradient", test_geom_1d_p1_pressure_gradient },
  { "geom_1d_p1_zero_b", test_geom_1d_p1_zero_b },
  { "consume_1d_p1_sign_and_mapping", test_consume_1d_p1_sign_and_mapping },
  { "consume_1d_p1_sign_flip", test_consume_1d_p1_sign_flip },
  { "consume_1d_p1_uniform_q_zero_source", test_consume_1d_p1_uniform_q_zero_source },
  { NULL, NULL },
};
