// Additional tests for the PKPM variables updater in 2 configuration-space
// dimensions (gkyl_dg_calc_pkpm_vars, 2x p1).  These complement the 1x tests
// by exercising the 2x kernels and a different magnetic-field orientation.
//
//   pressure (b = y_hat): byby = 1 -> Pyy = p_par, Pxx = Pzz = p_perp.
//   integrated_vars      : cell integral of moments (2x normalization 1/2^cdim).
//   u                    : weak divide rhou_i / rho with uniform fields.
//
// For a 2x p1 orthonormal basis the cell-average basis is 1/sqrt(2^2) = 1/2,
// so a field with physical value V has component-0 equal to V*sqrt(4).
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

// set cell-average physical value `val` of block `blk` for cdim=2 (factor 2).
static void
set_const_block_2x(struct gkyl_array *arr, int blk, int nb, double val)
{
  gkyl_array_shiftc(arr, val*2.0, blk*nb);
}

static void
mk_setup(struct gkyl_rect_grid *grid, struct gkyl_basis *cbasis,
  struct gkyl_range *local, struct gkyl_range *local_ext,
  struct gkyl_wv_eqn **eqn, struct gkyl_wave_geom **geom)
{
  int poly_order = 1, cdim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {4, 4};
  gkyl_rect_grid_init(grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(cbasis, cdim, poly_order);
  int ghost[] = {1, 1};
  gkyl_create_grid_ranges(grid, ghost, local_ext, local);
  *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  *geom = gkyl_wave_geom_new(grid, local_ext, 0, 0, false);
}

void
test_pressure_2x_p1_by()
{
  struct gkyl_rect_grid grid; struct gkyl_basis cbasis;
  struct gkyl_range local, local_ext;
  struct gkyl_wv_eqn *eqn; struct gkyl_wave_geom *geom;
  mk_setup(&grid, &cbasis, &local, &local_ext, &eqn, &geom);
  int nb = cbasis.num_basis; // 4

  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *bvar = mkarr(9*nb, local_ext.volume);
  struct gkyl_array *p_ij = mkarr(6*nb, local_ext.volume);

  double rho = 1.0, p_par = 5.0, p_perp = 2.0;
  gkyl_array_clear(moms, 0.0);
  set_const_block_2x(moms, 0, nb, rho);
  set_const_block_2x(moms, 1, nb, p_par);
  set_const_block_2x(moms, 2, nb, p_perp);

  // b = y_hat: by = 1 (block 1), byby = 1 (block 6).
  gkyl_array_clear(bvar, 0.0);
  set_const_block_2x(bvar, 1, nb, 1.0); // by
  set_const_block_2x(bvar, 6, nb, 1.0); // byby

  gkyl_array_clear(p_ij, 0.0);
  gkyl_dg_calc_pkpm_vars_pressure(up, &local, bvar, moms, p_ij);

  double s = 2.0; // sqrt(2^cdim)
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *P = gkyl_array_cfetch(p_ij, loc);
    double Pxx = P[0*nb]/s, Pxy = P[1*nb]/s, Pxz = P[2*nb]/s;
    double Pyy = P[3*nb]/s, Pyz = P[4*nb]/s, Pzz = P[5*nb]/s;
    TEST_CHECK( gkyl_compare(Pyy, p_par,  1e-12) );
    TEST_CHECK( gkyl_compare(Pxx, p_perp, 1e-12) );
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
test_integrated_vars_2x_p1()
{
  struct gkyl_rect_grid grid; struct gkyl_basis cbasis;
  struct gkyl_range local, local_ext;
  struct gkyl_wv_eqn *eqn; struct gkyl_wave_geom *geom;
  mk_setup(&grid, &cbasis, &local, &local_ext, &eqn, &geom);
  int nb = cbasis.num_basis;

  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *prim = mkarr(9*nb, local_ext.volume);
  struct gkyl_array *int_vars = mkarr(9, local_ext.volume);

  double rho = 3.0, p_par = 2.5, p_perp = 0.75;
  double ux = 2.0, uy = 1.0, uz = -1.5;
  double rhoux = rho*ux, rhouy = rho*uy, rhouz = rho*uz;

  gkyl_array_clear(moms, 0.0);
  set_const_block_2x(moms, 0, nb, rho);
  set_const_block_2x(moms, 1, nb, p_par);
  set_const_block_2x(moms, 2, nb, p_perp);

  gkyl_array_clear(euler, 0.0);
  set_const_block_2x(euler, 0, nb, rhoux);
  set_const_block_2x(euler, 1, nb, rhouy);
  set_const_block_2x(euler, 2, nb, rhouz);

  gkyl_array_clear(prim, 0.0);
  set_const_block_2x(prim, 0, nb, ux);
  set_const_block_2x(prim, 1, nb, uy);
  set_const_block_2x(prim, 2, nb, uz);

  gkyl_array_clear(int_vars, 0.0);
  gkyl_dg_calc_pkpm_integrated_vars(up, &local, moms, euler, prim, int_vars);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *I = gkyl_array_cfetch(int_vars, loc);
    TEST_CHECK( gkyl_compare(I[0], rho,       1e-12) );
    TEST_CHECK( gkyl_compare(I[1], rhoux,     1e-12) );
    TEST_CHECK( gkyl_compare(I[2], rhouy,     1e-12) );
    TEST_CHECK( gkyl_compare(I[3], rhouz,     1e-12) );
    TEST_CHECK( gkyl_compare(I[4], rhoux*ux,  1e-12) );
    TEST_CHECK( gkyl_compare(I[5], rhouy*uy,  1e-12) );
    TEST_CHECK( gkyl_compare(I[6], rhouz*uz,  1e-12) );
    TEST_CHECK( gkyl_compare(I[7], p_par,     1e-12) );
    TEST_CHECK( gkyl_compare(I[8], p_perp,    1e-12) );
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
test_u_2x_p1()
{
  struct gkyl_rect_grid grid; struct gkyl_basis cbasis;
  struct gkyl_range local, local_ext;
  struct gkyl_wv_eqn *eqn; struct gkyl_wave_geom *geom;
  mk_setup(&grid, &cbasis, &local, &local_ext, &eqn, &geom);
  int nb = cbasis.num_basis;

  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis,
    &local, eqn, geom, 0.0, false);

  struct gkyl_array *moms = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *pkpm_u = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *cell_avg_prim = gkyl_array_new(GKYL_INT, 1, local_ext.volume);
  memset(cell_avg_prim->data, 0, cell_avg_prim->size*cell_avg_prim->esznc);

  double rho = 5.0, rhoux = 10.0, rhouy = -5.0, rhouz = 15.0;
  gkyl_array_clear(moms, 0.0);
  set_const_block_2x(moms, 0, nb, rho);
  gkyl_array_clear(euler, 0.0);
  set_const_block_2x(euler, 0, nb, rhoux);
  set_const_block_2x(euler, 1, nb, rhouy);
  set_const_block_2x(euler, 2, nb, rhouz);

  gkyl_array_clear(pkpm_u, 0.0);
  gkyl_dg_calc_pkpm_vars_u(up, moms, euler, cell_avg_prim, pkpm_u);

  double s = 2.0;
  double ux_exp = rhoux/rho, uy_exp = rhouy/rho, uz_exp = rhouz/rho;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *u = gkyl_array_cfetch(pkpm_u, loc);
    TEST_CHECK( gkyl_compare(u[0*nb]/s, ux_exp, 1e-12) );
    TEST_CHECK( gkyl_compare(u[1*nb]/s, uy_exp, 1e-12) );
    TEST_CHECK( gkyl_compare(u[2*nb]/s, uz_exp, 1e-12) );
  }

  gkyl_array_release(cell_avg_prim);
  gkyl_array_release(pkpm_u);
  gkyl_array_release(euler);
  gkyl_array_release(moms);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(eqn);
}

TEST_LIST = {
  { "pressure_2x_p1_by",     test_pressure_2x_p1_by },
  { "integrated_vars_2x_p1", test_integrated_vars_2x_p1 },
  { "u_2x_p1",               test_u_2x_p1 },
  { NULL, NULL },
};
