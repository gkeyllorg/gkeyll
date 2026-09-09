// Tests for the PKPM variables updater (gkyl_dg_calc_pkpm_vars).
//
// We exercise four host-side compute methods with analytically known inputs
// that are *uniform* in configuration space (only the cell-average / 0th
// orthonormal basis component is non-zero in each block).  For a 1D p1 conf
// basis the orthonormal cell-average basis is psi_0 = 1/sqrt(2), so a field
// with physical value V is represented by setting component 0 to V*sqrt(2).
//
// Methods tested:
//   pressure        : p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
//                     With b = (1,0,0) (so bxbx = 1, rest 0):
//                       Pxx = p_par, Pyy = Pzz = p_perp, off-diag = 0.
//   u               : weak-divide rhou_i / rho.  With uniform rho and rhou_i,
//                     u_i = rhou_i / rho (constant).
//   integrated_vars : cell integral of (rho, rhoux, rhouy, rhouz,
//                     rho ux^2, rho uy^2, rho uz^2, p_par, p_perp).
//   source          : Lorentz force q/m rho(E + u x B) momentum source.
//
#include <acutest.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wave_geom.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Set the cell-average (physical) value of block `blk` (each block has
// `nb` basis components) to `val` for *every* cell.  For a 1D p1 orthonormal
// basis psi_0 = 1/sqrt(2), so component 0 must hold val*sqrt(2).
static void
set_const_block(struct gkyl_array *arr, int blk, int nb, double val)
{
  gkyl_array_shiftc(arr, val*sqrt(2.0), blk*nb);
}

void
test_pressure_1x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {6};
  int cdim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis; // 2

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);

  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  // vlasov_pkpm_moms = [rho, p_par, p_perp] (3 blocks).
  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  // bvar = [bx, by, bz, bxbx, bxby, bxbz, byby, bybz, bzbz] (9 blocks).
  struct gkyl_array *bvar = mkarr(9*nb, local_ext.volume);
  struct gkyl_array *p_ij = mkarr(6*nb, local_ext.volume);

  double rho = 2.0, p_par = 3.5, p_perp = 1.25;
  gkyl_array_clear(moms, 0.0);
  set_const_block(moms, 0, nb, rho);
  set_const_block(moms, 1, nb, p_par);
  set_const_block(moms, 2, nb, p_perp);

  gkyl_array_clear(bvar, 0.0);
  set_const_block(bvar, 0, nb, 1.0); // bx = 1
  set_const_block(bvar, 3, nb, 1.0); // bxbx = 1

  gkyl_array_clear(p_ij, 0.0);
  gkyl_dg_calc_pkpm_vars_pressure(up, &local, bvar, moms, p_ij);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  double s2 = sqrt(2.0);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *P = gkyl_array_cfetch(p_ij, loc);
    // Block layout: Pxx(0), Pxy(1), Pxz(2), Pyy(3), Pyz(4), Pzz(5).
    double Pxx = P[0*nb]/s2, Pxy = P[1*nb]/s2, Pxz = P[2*nb]/s2;
    double Pyy = P[3*nb]/s2, Pyz = P[4*nb]/s2, Pzz = P[5*nb]/s2;
    TEST_CHECK( gkyl_compare(Pxx, p_par,  1e-12) );
    TEST_CHECK( gkyl_compare(Pyy, p_perp, 1e-12) );
    TEST_CHECK( gkyl_compare(Pzz, p_perp, 1e-12) );
    TEST_CHECK( gkyl_compare(Pxy, 0.0,    1e-12) );
    TEST_CHECK( gkyl_compare(Pxz, 0.0,    1e-12) );
    TEST_CHECK( gkyl_compare(Pyz, 0.0,    1e-12) );
  }

  gkyl_array_release(p_ij);
  gkyl_array_release(bvar);
  gkyl_array_release(moms);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(eqn);
}

void
test_u_1x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {6};
  int cdim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis;

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);

  // Updater built over local (this is the mem_range used for the u solve).
  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(3*nb, local_ext.volume); // [rhoux, rhouy, rhouz]
  struct gkyl_array *pkpm_u = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *cell_avg_prim = gkyl_array_new(GKYL_INT, 1, local_ext.volume);
  // Int array: zero it via memset (gkyl_array_clear only supports GKYL_DOUBLE).
  memset(cell_avg_prim->data, 0, cell_avg_prim->size*cell_avg_prim->esznc);

  double rho = 4.0;
  double rhoux = 6.0, rhouy = -2.0, rhouz = 10.0;
  gkyl_array_clear(moms, 0.0);
  set_const_block(moms, 0, nb, rho);

  gkyl_array_clear(euler, 0.0);
  set_const_block(euler, 0, nb, rhoux);
  set_const_block(euler, 1, nb, rhouy);
  set_const_block(euler, 2, nb, rhouz);

  gkyl_array_clear(pkpm_u, 0.0);
  gkyl_dg_calc_pkpm_vars_u(up, moms, euler, cell_avg_prim, pkpm_u);

  double ux_exp = rhoux/rho, uy_exp = rhouy/rho, uz_exp = rhouz/rho;
  double s2 = sqrt(2.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *u = gkyl_array_cfetch(pkpm_u, loc);
    TEST_CHECK( gkyl_compare(u[0*nb]/s2, ux_exp, 1e-12) );
    TEST_CHECK( gkyl_compare(u[1*nb]/s2, uy_exp, 1e-12) );
    TEST_CHECK( gkyl_compare(u[2*nb]/s2, uz_exp, 1e-12) );
    // Slopes should vanish for a uniform field.
    TEST_CHECK( gkyl_compare(u[0*nb+1], 0.0, 1e-12) );
    TEST_CHECK( gkyl_compare(u[1*nb+1], 0.0, 1e-12) );
    TEST_CHECK( gkyl_compare(u[2*nb+1], 0.0, 1e-12) );
  }

  gkyl_array_release(cell_avg_prim);
  gkyl_array_release(pkpm_u);
  gkyl_array_release(euler);
  gkyl_array_release(moms);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(eqn);
}

void
test_integrated_vars_1x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {6};
  int cdim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis;

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);

  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *prim = mkarr(9*nb, local_ext.volume); // [ux,uy,uz,...]
  struct gkyl_array *int_vars = mkarr(9, local_ext.volume);

  double rho = 4.0, p_par = 3.0, p_perp = 1.5;
  double ux = 1.5, uy = -0.5, uz = 2.0;
  double rhoux = rho*ux, rhouy = rho*uy, rhouz = rho*uz;

  gkyl_array_clear(moms, 0.0);
  set_const_block(moms, 0, nb, rho);
  set_const_block(moms, 1, nb, p_par);
  set_const_block(moms, 2, nb, p_perp);

  gkyl_array_clear(euler, 0.0);
  set_const_block(euler, 0, nb, rhoux);
  set_const_block(euler, 1, nb, rhouy);
  set_const_block(euler, 2, nb, rhouz);

  gkyl_array_clear(prim, 0.0);
  set_const_block(prim, 0, nb, ux);
  set_const_block(prim, 1, nb, uy);
  set_const_block(prim, 2, nb, uz);

  gkyl_array_clear(int_vars, 0.0);
  gkyl_dg_calc_pkpm_integrated_vars(up, &local, moms, euler, prim, int_vars);

  // The kernel accumulates 0.7071067811865476*<block>[0] for the linear
  // entries and 0.5*(c[1]*c[1]+c[0]*c[0]) for the quadratic energy entries.
  // With our uniform fields component 1 == 0, and component 0 == V*sqrt(2),
  // so 0.70710678*(V*sqrt(2)) = V and 0.5*(0 + (V*sqrt2)*(W*sqrt2)) = V*W.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *I = gkyl_array_cfetch(int_vars, loc);
    TEST_CHECK( gkyl_compare(I[0], rho,        1e-12) );
    TEST_CHECK( gkyl_compare(I[1], rhoux,      1e-12) );
    TEST_CHECK( gkyl_compare(I[2], rhouy,      1e-12) );
    TEST_CHECK( gkyl_compare(I[3], rhouz,      1e-12) );
    TEST_CHECK( gkyl_compare(I[4], rhoux*ux,   1e-12) );
    TEST_CHECK( gkyl_compare(I[5], rhouy*uy,   1e-12) );
    TEST_CHECK( gkyl_compare(I[6], rhouz*uz,   1e-12) );
    TEST_CHECK( gkyl_compare(I[7], p_par,      1e-12) );
    TEST_CHECK( gkyl_compare(I[8], p_perp,     1e-12) );
  }

  gkyl_array_release(int_vars);
  gkyl_array_release(prim);
  gkyl_array_release(euler);
  gkyl_array_release(moms);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(eqn);
}

void
test_source_1x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {6};
  int cdim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis;

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);

  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  // qmem = q/m*[Ex,Ey,Ez,Bx,By,Bz,phi,psi] (8 blocks).
  struct gkyl_array *qmem = mkarr(8*nb, local_ext.volume);
  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *rhs = mkarr(3*nb, local_ext.volume);

  // Choose B along z only, E along x only, so source = rho*E + rhou x B.
  double rho = 3.0;
  double Ex = 2.0, Ey = 0.0, Ez = 0.0;
  double Bx = 0.0, By = 0.0, Bz = 5.0;
  double rhoux = 1.0, rhouy = 4.0, rhouz = -2.0;

  gkyl_array_clear(qmem, 0.0);
  set_const_block(qmem, 0, nb, Ex);
  set_const_block(qmem, 1, nb, Ey);
  set_const_block(qmem, 2, nb, Ez);
  set_const_block(qmem, 3, nb, Bx);
  set_const_block(qmem, 4, nb, By);
  set_const_block(qmem, 5, nb, Bz);

  gkyl_array_clear(moms, 0.0);
  set_const_block(moms, 0, nb, rho);

  gkyl_array_clear(euler, 0.0);
  set_const_block(euler, 0, nb, rhoux);
  set_const_block(euler, 1, nb, rhouy);
  set_const_block(euler, 2, nb, rhouz);

  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_calc_pkpm_vars_source(up, &local, qmem, moms, euler, rhs);

  // Expected Lorentz-force cell-average source:
  //   out_x = rho*Ex + (rhouy*Bz - rhouz*By)
  //   out_y = rho*Ey + (rhouz*Bx - rhoux*Bz)
  //   out_z = rho*Ez + (rhoux*By - rhouy*Bx)
  double sx = rho*Ex + (rhouy*Bz - rhouz*By);
  double sy = rho*Ey + (rhouz*Bx - rhoux*Bz);
  double sz = rho*Ez + (rhoux*By - rhouy*Bx);
  double s2 = sqrt(2.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *o = gkyl_array_cfetch(rhs, loc);
    TEST_CHECK( gkyl_compare(o[0*nb]/s2, sx, 1e-12) );
    TEST_CHECK( gkyl_compare(o[1*nb]/s2, sy, 1e-12) );
    TEST_CHECK( gkyl_compare(o[2*nb]/s2, sz, 1e-12) );
  }

  gkyl_array_release(rhs);
  gkyl_array_release(euler);
  gkyl_array_release(moms);
  gkyl_array_release(qmem);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "pressure_1x_p1",        test_pressure_1x_p1 },
  { "u_1x_p1",               test_u_1x_p1 },
  { "integrated_vars_1x_p1", test_integrated_vars_1x_p1 },
  { "source_1x_p1",          test_source_1x_p1 },
  { NULL, NULL },
};
