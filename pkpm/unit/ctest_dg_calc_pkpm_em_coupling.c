// Tests for the PKPM fluid-EM coupling updater (gkyl_dg_calc_pkpm_em_coupling).
//
// The updater performs an implicit, time-centered solve of the coupled
// momentum + Ampere system.  We test a fully analytically tractable regime:
//
//   * single species
//   * static self-consistent field (pkpm_field_static = true), so Ampere's
//     law reduces to Ẽ = epsilon0 * E^n and the electric field is unchanged.
//   * B = 0 (em B, external B all zero) and no external forces / currents.
//
// In this regime the time-centered solve gives
//     u_bar = u^n + 0.5 dt (q/m) E^n
//     u^{n+1} = 2 u_bar - u^n = u^n + dt (q/m) E^n        (a pure E-field push)
//     E^{n+1} = E^n                                        (field is static)
// and the output momentum is rho^n u^{n+1}.
//
// All fields are spatially uniform so we can compare cell-average values.
#include <acutest.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_em_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Set the cell-average (physical) value of block `blk` (each block has `nb`
// basis components) to `val` for every cell.  1D p1 orthonormal: psi_0=1/sqrt2.
static void
set_const_block(struct gkyl_array *arr, int blk, int nb, double val)
{
  gkyl_array_shiftc(arr, val*sqrt(2.0), blk*nb);
}

void
test_em_coupling_static_Epush_1x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {5};
  int cdim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis; // 2

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int num_species = 1;
  double qbym[GKYL_MAX_SPECIES] = {0.0};
  qbym[0] = -1.7;           // charge/mass of the single species
  double epsilon0 = 1.4;
  bool field_static = true;

  struct gkyl_dg_calc_pkpm_em_coupling *up =
    gkyl_dg_calc_pkpm_em_coupling_new(&cbasis, &local, num_species, qbym,
      epsilon0, field_static, false);
  TEST_CHECK( up != NULL );

  // Input arrays.
  struct gkyl_array *app_accel0 = mkarr(3*nb, local_ext.volume);
  struct gkyl_array *ext_em     = mkarr(6*nb, local_ext.volume);
  struct gkyl_array *app_current= mkarr(3*nb, local_ext.volume);
  struct gkyl_array *moms0      = mkarr(3*nb, local_ext.volume); // [rho,p_par,p_perp]
  struct gkyl_array *u0         = mkarr(3*nb, local_ext.volume); // [ux,uy,uz]
  struct gkyl_array *euler0     = mkarr(3*nb, local_ext.volume); // output [rhoux,...]
  struct gkyl_array *em         = mkarr(6*nb, local_ext.volume); // [E(3),B(3)]

  gkyl_array_clear(app_accel0, 0.0);
  gkyl_array_clear(ext_em, 0.0);
  gkyl_array_clear(app_current, 0.0);
  gkyl_array_clear(euler0, 0.0);

  double rho = 2.5;
  double ux0 = 0.3, uy0 = -1.0, uz0 = 0.7;
  double Ex = 1.2, Ey = -0.4, Ez = 2.0; // B = 0

  gkyl_array_clear(moms0, 0.0);
  set_const_block(moms0, 0, nb, rho);

  gkyl_array_clear(u0, 0.0);
  set_const_block(u0, 0, nb, ux0);
  set_const_block(u0, 1, nb, uy0);
  set_const_block(u0, 2, nb, uz0);

  gkyl_array_clear(em, 0.0);
  set_const_block(em, 0, nb, Ex);
  set_const_block(em, 1, nb, Ey);
  set_const_block(em, 2, nb, Ez);
  // B (blocks 3,4,5) left at zero.

  double dt = 0.25;

  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES] = {app_accel0};
  const struct gkyl_array *moms[GKYL_MAX_SPECIES] = {moms0};
  const struct gkyl_array *u[GKYL_MAX_SPECIES] = {u0};
  struct gkyl_array *euler[GKYL_MAX_SPECIES] = {euler0};

  gkyl_dg_calc_pkpm_em_coupling_advance(up, dt, app_accel, ext_em, app_current,
    moms, u, euler, em);

  // Expected new velocities (pure E push): u^{n+1} = u^n + dt*(q/m)*E^n.
  double uxn = ux0 + dt*qbym[0]*Ex;
  double uyn = uy0 + dt*qbym[0]*Ey;
  double uzn = uz0 + dt*qbym[0]*Ez;
  double rhoux = rho*uxn, rhouy = rho*uyn, rhouz = rho*uzn;
  double s2 = sqrt(2.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *eul = gkyl_array_cfetch(euler0, loc);
    TEST_CHECK( gkyl_compare(eul[0*nb]/s2, rhoux, 1e-12) );
    TEST_CHECK( gkyl_compare(eul[1*nb]/s2, rhouy, 1e-12) );
    TEST_CHECK( gkyl_compare(eul[2*nb]/s2, rhouz, 1e-12) );
    // Static field: E must be unchanged.
    const double *emd = gkyl_array_cfetch(em, loc);
    TEST_CHECK( gkyl_compare(emd[0*nb]/s2, Ex, 1e-12) );
    TEST_CHECK( gkyl_compare(emd[1*nb]/s2, Ey, 1e-12) );
    TEST_CHECK( gkyl_compare(emd[2*nb]/s2, Ez, 1e-12) );
  }

  gkyl_array_release(em);
  gkyl_array_release(euler0);
  gkyl_array_release(u0);
  gkyl_array_release(moms0);
  gkyl_array_release(app_current);
  gkyl_array_release(ext_em);
  gkyl_array_release(app_accel0);
  gkyl_dg_calc_pkpm_em_coupling_release(up);
}

void
test_em_coupling_new_2x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {4, 4};
  int cdim = 2;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double qbym[GKYL_MAX_SPECIES] = {0.0};
  qbym[0] = 1.0; qbym[1] = -1.0;

  struct gkyl_dg_calc_pkpm_em_coupling *up =
    gkyl_dg_calc_pkpm_em_coupling_new(&cbasis, &local, 2, qbym, 1.0, false, false);
  TEST_CHECK( up != NULL );
  gkyl_dg_calc_pkpm_em_coupling_release(up);
}

TEST_LIST = {
  { "em_coupling_static_Epush_1x_p1", test_em_coupling_static_Epush_1x_p1 },
  { "em_coupling_new_2x_p1",          test_em_coupling_new_2x_p1 },
  { NULL, NULL },
};
